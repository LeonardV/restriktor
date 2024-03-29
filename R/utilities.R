## utility functions
coef.restriktor <- function(object, ...)  {
  
  b.def <- c()
  b.restr <- object$b.restr
  
  if (any(object$parTable$op == ":=")) {
    b.def <- object$CON$def.function(object$b.restr)
  }
  
  if (inherits(object, "conMLM")) {
    OUT <- rbind(b.restr, b.def)
  } else {
    OUT <- c(b.restr, b.def)
  }
  
  return(OUT)
}


coef.con_goric <- function(object, ...)  {
  return(object$ormle$b.restr)
}

coef.gorica_est <- function(object, ...)  {
  return(object$b.restr)
}


logLik.restriktor <- function(object, ...) {
  return(object$loglik)
}


model.matrix.restriktor <- function(object, ...) {
  return(model.matrix(object$model.org))
}


tukeyChi <- function(x, c = 4.685061, deriv = 0, ...) {
  u <- x / c
  out <- abs(x) > c
  if (deriv == 0) { # rho function
    r <- 1 - (1 - u^2)^3
    r[out] <- 1
  } else if (deriv == 1) { # rho' = psi function
    r <- 6 * x * (1 - u^2)^2 / c^2
    r[out] <- 0
  } else if (deriv == 2) { # rho'' 
    r <- 6 * (1 - u^2) * (1 - 5 * u^2) / c^2
    r[out] <- 0
  } else {
    stop("deriv must be in {0,1,2}")
  }
  return(r)
}


# code taken from robustbase package.
# addapted by LV (3-12-2017).
robWeights <- function(w, eps = 0.1/length(w), eps1 = 0.001, ...) {
  stopifnot(is.numeric(w))
  cat("Robustness weights:", "\n")
  cat0 <- function(...) cat("", ...)
  n <- length(w)
  if (n <= 10) 
    print(w, digits = 5, ...)
  else {
    n1 <- sum(w1 <- abs(w - 1) < eps1)
    n0 <- sum(w0 <- abs(w) < eps)
    if (any(w0 & w1)) 
      warning("weights should not be both close to 0 and close to 1!\n", 
              "You should use different 'eps' and/or 'eps1'")
    if (n0 > 0 || n1 > 0) {
      if (n0 > 0) {
        formE <- function(e) formatC(e, digits = max(2, 
                                                     5 - 3), width = 1)
        i0 <- which(w0)
        maxw <- max(w[w0])
        c3 <- paste0("with |weight| ", if (maxw == 0) 
          "= 0"
          else paste("<=", formE(maxw)), " ( < ", formE(eps), 
          ");")
        cat0(if (n0 > 1) {
          cc <- sprintf("%d observations c(%s)", n0, 
                        strwrap(paste(i0, collapse = ",")))
          c2 <- " are outliers"
          paste0(cc, if (nchar(cc) + nchar(c2) + nchar(c3) > 
                         getOption("width")) 
            "\n\t", c2)
        }
        else sprintf("observation %d is an outlier", 
                     i0), c3, "\n")
      }
      if (n1 > 0) 
        cat0(ngettext(n1, "one weight is", sprintf("%s%d weights are", 
                                                   if (n1 == n) 
                                                     "All "
                                                   else "", n1)), "~= 1.")
      n.rem <- n - n0 - n1
      if (n.rem <= 0) {
        if (n1 > 0) 
          cat("\n")
        return(invisible())
      }
    }
  }
}




# function taken from 'bain' package 
expand_compound_constraints <- function(hyp) {
  equality_operators <- gregexpr("[=<>]", hyp)[[1]]
  if(length(equality_operators) > 1){
    string_positions <- c(0, equality_operators, nchar(hyp)+1)
    return(sapply(1:(length(string_positions)-2), function(pos) {
      substring(hyp, (string_positions[pos]+1), (string_positions[pos+2]-1))
    }))
  } else {
    return(hyp)
  }
}

# function taken from 'bain' package 
expand_parentheses <- function(hyp) {
  parenth_locations <- gregexpr("[\\(\\)]", hyp)[[1]]
  if (!parenth_locations[1] == -1 & !grepl("abs\\(.*\\)", hyp) ) {
    if (length(parenth_locations) %% 2 > 0) stop("Not all opening parentheses are matched by a closing parenthesis, or vice versa.")
    expanded_contents <- strsplit(substring(hyp, (parenth_locations[1]+1), (parenth_locations[2]-1)), ",")[[1]]
    if (length(parenth_locations) == 2){
      return(paste0(substring(hyp, 1, (parenth_locations[1]-1)), expanded_contents, substring(hyp, (parenth_locations[2]+1), nchar(hyp))))
    } else {
      return(apply(
        expand.grid(expanded_contents, expand_parentheses(substring(hyp, (parenth_locations[2]+1), nchar(hyp)))),
        1, paste, collapse = ""))
    }
  } else {
    return(hyp)
  }
}


format_numeric <- function(x, digits = 3) {
  if (abs(x) <= 1e-8) {
    format(0, nsmall = digits)
  } else if (abs(x) >= 1e3 || abs(x) <= 1e-3) {
    format(x, scientific = TRUE, digits = digits)
  } else {
    format(round(x, digits), nsmall = digits) 
  }
}


calculate_model_comparison_metrics <- function(x) {
  modelnames <- as.character(x$model)
  ## Log-likelihood
  LL = -2 * x$loglik
  delta_LL = LL - min(LL)
  loglik_weights = exp(0.5 * -delta_LL) / sum(exp(0.5 * -delta_LL))
  loglik_rw  <- loglik_weights %*% t(1/loglik_weights)
  diag(loglik_rw) = 1
  
  ## penalty
  penalty_weights = exp(-x$penalty) / sum(exp(-x$penalty))
  penalty_rw = penalty_weights %*% t(1/penalty_weights)
  diag(penalty_rw) = 1
  
  ## goric
  delta_goric = x$goric - min(x$goric)
  goric_weights = exp(0.5 * -delta_goric) / sum(exp(0.5 * -delta_goric))
  goric_rw = goric_weights %*% t(1/goric_weights)
  diag(goric_rw) = 1
  
  rownames(goric_rw)   = modelnames
  rownames(penalty_rw) = modelnames
  rownames(loglik_rw)  = modelnames
  colnames(goric_rw)   = paste0("vs. ", modelnames)
  colnames(penalty_rw) = paste0("vs. ", modelnames)
  colnames(loglik_rw)  = paste0("vs. ", modelnames)
  
  out <- list(loglik_weights  = loglik_weights, 
              penalty_weights = penalty_weights,
              goric_weights   = goric_weights,
              loglik_rw       = loglik_rw,
              penalty_rw      = penalty_rw,
              goric_rw        = goric_rw)
  
  return(out)
}

# this function is called from the goric_benchmark_anova() function
parallel_function <- function(i, samplesize, var.e, nr.iter, means_pop, 
                              hypos, PrefHypo, object, n.coef, sample, 
                              control, ...) {  
  # Sample residuals
  epsilon <- rnorm(sum(samplesize), sd = sqrt(var.e))
  # Generate data
  sample$y <- as.matrix(sample[, 2:(1 + n.coef)]) %*% matrix(means_pop, 
                                                             nrow = n.coef) + epsilon
  df <- data.frame(y = sample$y, sample[, 2:(1 + n.coef)])
  
  # Obtain fit
  fit <- lm(y ~ 0 + ., data = df)
  # GORICA or GORICA depending on what is done in data
  results.goric <- goric(fit,
                         hypotheses = hypos,
                         comparison = object$comparison,
                         type = object$type,
                         control = control, 
                         ...)
  
  # Return the relevant results
  list(
    #test  = attr(results.goric$objectList[[results.goric$objectNames]]$wt.bar, "mvtnorm"),
    goric = results.goric$result[PrefHypo, 7],
    gw    = results.goric$ratio.gw[PrefHypo, ],
    lw    = results.goric$ratio.lw[PrefHypo, ],
    ld    = (results.goric$result$loglik[PrefHypo] - results.goric$result$loglik)
  )
}

# compute_weights_ratioWeights <- function(x) {
#   IC <- 2*x
#   minIC <- min(IC)
#   weights <- exp(-0.5 * (IC - minIC)) / sum(exp(-0.5 * (IC - minIC)))
#   ratio_weights <- weights %*% t(1/weights)
# 
#   out <- list(weights = weights, ratio_weights = ratio_weights)
# 
#   return(out)
# }

# 
# remove_linear_dependent_rows_matrix <- function(Amat, bvec) {
#   ## remove any linear dependent rows from the constraint matrix. Amat must be of full row rank.
#   # remove any zero vectors
#   allZero.idx <- rowSums(abs(Amat)) == 0
#   Amat <- Amat[!allZero.idx, , drop = FALSE]
#   bvec <- bvec[!allZero.idx]
#   # what is the rank of Amat
#   rank <- qr(Amat)$rank
#   # decompose Amat using svd
#   s <- svd(Amat)
#   # continue untill Amat is of full-row rank
#   while (rank != length(s$d)) {
#     # check which singular values are zero
#     zero.idx <- which(zapsmall(s$d) <= 1e-16)
#     # remove linear dependent rows and reconstruct the constraint matrix
#     Amat <- s$u[-zero.idx, ] %*% diag(s$d) %*% t(s$v)
#     # zapping small ones to zero
#     Amat <- zapsmall(Amat)
#     bvec <- bvec[-zero.idx]
#     s <- svd(Amat)
#   }
#  
#   OUT <- list(Amat, bvec)
#   
#   OUT
# }
# 


#rankifremoved <- sapply(1:ncol(Amat), function (x) qr(Amat[-x, ])$rank)
#which(rankifremoved == max(rankifremoved))
