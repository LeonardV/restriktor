# ic_weight and nchoosk are taken from ic.infer package
# slightly adjusted by LV (2024-08-19)
# check for small eigenvalues, abs(values) <= .Machine$double.eps are replaced by zero
# using ridging to deal with exactly singular corr matrix. A small constant (1e-05) 
# value is added to the diagonal of corr. 
# Optimized: analytical shortcuts for dim 1-3

make_positive_definite <- function(mat, tolerance = 1e-15, 
                                   ridge_constant = 1e-05) {
  ev <- eigen(mat)
  ev$values[ev$values <= tolerance] <- 0  
  ev_matrix <- diag(ev$values, nrow = length(ev$values))
  corr <- ev$vectors %*% ev_matrix %*% t(ev$vectors)
  corr_inv <- try(solve(corr), silent = TRUE)
  if (inherits(corr_inv, "try-error")) {
    diag(corr) <- diag(corr) + ridge_constant
    corr_inv <- solve(corr)
  } 
  list(corr = corr, corr_inv = corr_inv)
}


# =============================================================================
# Analytical ic_weights (Shapiro 1985, via ic.infer)
# Uses exact orthant probabilities for dim 1-3 (Sheppard/Plackett),
# pmvnorm for dim >= 4.
# =============================================================================

# P(X > 0 | X ~ N(0, sigma)): exact for dim 1-3, pmvnorm for dim >= 4
.orthant_prob <- function(sigma, lower, upper, ...) {
  k <- NROW(sigma)
  
  if (k == 1L) return(0.5)
  
  if (k == 2L) {
    rho <- sigma[1L, 2L] / sqrt(sigma[1L, 1L] * sigma[2L, 2L])
    rho <- max(-1, min(1, rho))
    return(0.25 + asin(rho) / (2 * pi))
  }
  
  if (k == 3L) {
    sd_inv <- 1 / sqrt(c(sigma[1L, 1L], sigma[2L, 2L], sigma[3L, 3L]))
    r12 <- max(-1, min(1, sigma[1L, 2L] * sd_inv[1L] * sd_inv[2L]))
    r13 <- max(-1, min(1, sigma[1L, 3L] * sd_inv[1L] * sd_inv[3L]))
    r23 <- max(-1, min(1, sigma[2L, 3L] * sd_inv[2L] * sd_inv[3L]))
    return(0.125 + (asin(r12) + asin(r13) + asin(r23)) / (4 * pi))
  }
  
  pmvnorm(lower = lower, upper = upper, sigma = sigma, ...)
}

ic_weights <- function(corr, tolerance, ridge_constant, ...) {
  
  if (!(is.matrix(corr) && nrow(corr) == ncol(corr))) 
    stop("corr must be a square matrix.")
  
  g <- nrow(corr)
  
  # Make corr positive definite; returns both corr and its inverse.
  # This avoids redundant eigen() and solve() calls.
  pd <- make_positive_definite(corr, tolerance = tolerance, 
                               ridge_constant = ridge_constant)
  corr     <- pd$corr
  corr_inv <- pd$corr_inv
  
  if (g == 1) {
    return(c("0" = 0.5, "1" = 0.5))
  }
  
  liste <- 1:g
  weights <- rep(0, g + 1)
  names(weights) <- rev(0:g)
  
  lower_g <- rep(0, g)
  upper_g <- rep(Inf, g)
  
  weights[1] <- .orthant_prob(corr, lower_g, upper_g, ...)
  weights[g + 1] <- .orthant_prob(corr_inv, lower_g, upper_g, ...)
  
  if (g > 4) {
    test_matrix <- try(matrix(0, floor((g - 2)/2), choose(g, floor((g - 2)/2))), silent = TRUE)
    if (!is.numeric(test_matrix)) 
      stop(paste("ic_weights will not work, corr too large, \n", 
                 "interim matrix with ", floor((g - 2)/2) * choose(g, floor((g - 2)/2)), 
                 " elements cannot be created.", sep = ""))
  }
  
  if (g == 2) {
    weights[2] <- 1 - sum(weights)
  }
  
  if (g == 3) {
    weights[2] <- 0.5 - weights[4]
    weights[3] <- 0.5 - weights[1]
  }
  
  if (g > 3) {
    half_floor <- floor((g - 2) / 2)
    
    bounds_cache <- list()
    for (size in 1:g) {
      bounds_cache[[size]] <- list(lower = rep(0, size), upper = rep(Inf, size))
    }
    
    for (k in 1:half_floor) {
      jetzt <- nchoosek(g, k)
      n_combinations <- ncol(jetzt)
      wjetzt <- matrix(0, n_combinations, 2)
      
      lower_k <- bounds_cache[[k]]$lower
      upper_k <- bounds_cache[[k]]$upper
      gk <- g - k
      lower_gk <- bounds_cache[[gk]]$lower
      upper_gk <- bounds_cache[[gk]]$upper
      
      is_symmetric <- (k == (g - 2) / 2)
      
      for (j in 1:n_combinations) {
        diese <- jetzt[, j]
        andere <- setdiff(liste, diese)
        
        corr_diese_inv <- solve(corr[diese, diese, drop = FALSE])
        
        corr_andere_diese <- corr[andere, diese, drop = FALSE]
        corr_diese_andere <- corr[diese, andere, drop = FALSE]
        corr_andere <- corr[andere, andere, drop = FALSE]
        
        hilf1 <- corr_andere - corr_andere_diese %*% (corr_diese_inv %*% corr_diese_andere)
        
        wjetzt[j, 1] <- .orthant_prob(corr_diese_inv, lower_k, upper_k, ...) * 
          .orthant_prob(hilf1, lower_gk, upper_gk, ...)
        
        if (!is_symmetric) {
          corr_andere_inv <- solve(corr[andere, andere, drop = FALSE])
          
          hilf2 <- corr[diese, diese, drop = FALSE] - 
            corr_diese_andere %*% (corr_andere_inv %*% corr_andere_diese)
          
          wjetzt[j, 2] <- .orthant_prob(corr_andere_inv, lower_gk, upper_gk, ...) * 
            .orthant_prob(hilf2, lower_k, upper_k, ...)
        }
      }
      
      weights[k + 1] <- sum(wjetzt[, 1])
      weights[g + 1 - k] <- sum(wjetzt[, 2])
    }
    
    if (g %% 2 == 0) {
      even.sum <- sum(weights[1 + 2 * ((g/2):0)])
      odd.sum <- sum(weights[2 * ((g/2):1)])
      if ((g/4) == floor(g/4)) {
        weights[g/2 + 1] <- 0.5 - even.sum
        weights[g/2 + 2] <- 0.5 - odd.sum
      } else {
        weights[g/2 + 1] <- 0.5 - odd.sum
        weights[g/2 + 2] <- 0.5 - even.sum
      }
    } else {
      even.sum <- sum(weights[2 * (((g + 1)/2):1)])
      odd.sum <- sum(weights[2 * (((g + 1)/2):1) - 1])
      if (((g + 1)/4) == floor((g + 1)/4)) {
        weights[(g + 1)/2] <- 0.5 - even.sum
        weights[(g + 3)/2] <- 0.5 - odd.sum
      } else {
        weights[(g + 1)/2] <- 0.5 - odd.sum
        weights[(g + 3)/2] <- 0.5 - even.sum
      }
    }
  }
  
  weights
}


# taken from ic.infer package, the function is not exported
nchoosek <- function (n, k) {
  if (!is.numeric(n) || !is.numeric(k) || is.na(n) || is.na(k) || 
      length(n) != 1 || length(k) != 1) 
    stop("arguments must be non-NA numeric scalars.")
  if (k > n || k < 0) 
    stop("Arguments must satisfy 0 <= k <= n.")
  nck = choose(n, k)
  res = matrix(NA, nrow = k, ncol = nck)
  res[, 1] = 1:k
  j = 2
  repeat {
    if (j > nck) 
      break
    res[, j] = res[, j - 1]
    i = k
    repeat {
      res[i, j] = res[i, j] + 1
      if (res[i, j] <= n - (k - i)) 
        break
      i = i - 1
      stopifnot(i >= 1)
    }
    if (i < k) 
      res[(i + 1):k, j] = res[i, j] + 1:(k - i)
    j = j + 1
  }
  stopifnot(all(res[, nck] == (n - k + 1):n))
  stopifnot(all(res <= n) && all(res >= 1))
  return(res)
}
