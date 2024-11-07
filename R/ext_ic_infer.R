# ic_weight and nchoosk are taken from ic.infer package
# slightly adjusted by LV (2024-08-19)
# check for small eigenvalues, abs(values) <= .Machine$double.eps are replaced by zero
# using ridging to deal with exactly singular corr matrix. A small constant (1e-05) 
# value is added to the diagonal of corr. 

make_positive_semi_definite <- function(mat, tolerance = 1e-15, 
                                        ridge_constant = 1e-05) {
  ev <- eigen(mat)
  ev$values[ev$values <= tolerance ] <- 0  
  ev_matrix <- diag(ev$values, nrow = length(ev$values))
  corr <- ev$vectors %*% ev_matrix %*% t(ev$vectors)
  corr_solve <- try(solve(corr), silent = TRUE)
  if (inherits(corr_solve, "try-error")) {
    diag(corr) <- diag(corr) + ridge_constant  
  } 
  return(corr)
}

ic_weights <- function(corr, tolerance, ridge_constant, ...) {
  if (!(is.matrix(corr) & nrow(corr) == ncol(corr))) 
    stop("corr must be a square matrix.")
  corr <- make_positive_semi_definite(corr, tolerance = tolerance, 
                                      ridge_constant = ridge_constant)
  if (!all(abs(eigen(corr)$values) > 0)) 
    stop("corr must be positive definite.")
  g <- nrow(corr)
  liste <- 1:g
  weights <- rep(0, g + 1)
  names(weights) <- rev(0:g)
  weights[1] <- pmvnorm(lower = rep(0, g), upper = rep(Inf, 
                                                       g), sigma = corr, ...)
  weights[g + 1] <- pmvnorm(lower = rep(0, g), upper = rep(Inf, 
                                                           g), sigma = solve(corr), ...)
  if (g > 4) {
    if (!is.numeric(try(matrix(0, floor((g - 2)/2), choose(g, 
                                                           floor((g - 2)/2))), silent = TRUE))) 
      stop(paste("ic_weights will not work, corr too large, \n", 
                 "interim matrix with ", floor((g - 2)/2) * choose(g, 
                                                                   floor((g - 2)/2)), " elements cannot be created.", 
                 sep = ""))
  }
  if (g == 2) 
    weights[2] <- 1 - sum(weights)
  if (g == 3) {
    weights[2] <- 0.5 - weights[4]
    weights[3] <- 0.5 - weights[1]
  }
  if (g > 3) {
    for (k in 1:floor((g - 2)/2)) {
      jetzt <- nchoosek(g, k)
      wjetzt <- matrix(0, choose(g, k), 2)
      for (j in 1:(choose(g, k))) {
        diese <- jetzt[, j]
        andere <- setdiff(liste, diese)
        hilf <- corr[andere, andere] - corr[andere, diese] %*% 
          solve(corr[diese, diese], matrix(corr[diese, 
                                                andere], k, g - k))
        wjetzt[j, 1] <- pmvnorm(lower = rep(0, k), upper = rep(Inf, 
                                                               k), sigma = solve(corr[diese, diese]), ...) * 
          pmvnorm(lower = rep(0, g - k), upper = rep(Inf, 
                                                     g - k), sigma = hilf, ...)
        if (!k == (g - 2)/2) {
          hilf <- corr[diese, diese] - corr[diese, andere] %*% 
            solve(corr[andere, andere], matrix(corr[andere, 
                                                    diese], g - k, k))
          wjetzt[j, 2] <- pmvnorm(lower = rep(0, g - 
                                                k), upper = rep(Inf, g - k), sigma = solve(corr[andere, 
                                                                                                andere]), ...) * pmvnorm(lower = rep(0, k), 
                                                                                                                         upper = rep(Inf, k), sigma = hilf, ...)
        }
      }
      weights[k + 1] <- sum(wjetzt[, 1])
      weights[g + 1 - k] <- sum(wjetzt[, 2])
    }
    if (g/2 == floor(g/2)) {
      even.sum <- sum(weights[1 + 2 * ((g/2):0)])
      odd.sum <- sum(weights[2 * ((g/2):1)])
      if (g/4 == floor(g/4)) {
        weights[g/2 + 1] <- 0.5 - even.sum
        weights[g/2 + 2] <- 0.5 - odd.sum
      }
      else {
        weights[g/2 + 1] <- 0.5 - odd.sum
        weights[g/2 + 2] <- 0.5 - even.sum
      }
    }
    else {
      even.sum <- sum(weights[2 * (((g + 1)/2):1)])
      odd.sum <- sum(weights[2 * (((g + 1)/2):1) - 1])
      if ((g + 1)/4 == floor((g + 1)/4)) {
        weights[(g + 1)/2] <- 0.5 - even.sum
        weights[(g + 3)/2] <- 0.5 - odd.sum
      }
      else {
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