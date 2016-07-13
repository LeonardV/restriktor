## utility functions

coef.conLM <- function(object, ...)  {
  object$b.restr
}

model.matrix.conLM <- function(object, ...) {
  model.matrix(object$model.org)
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
  r
}


# input list
goricWeights <- function(x, ...) {
  if (!is.list(x)) {
    stop("restriktor ERROR: x must be a list.")
  }
  ll    <- unlist(lapply(x, function(x) attr(x$goric, "loglik")))
  PT    <- unlist(lapply(x, function(x) attr(x$goric, "penalty")))
  goric <- unlist(lapply(x, function(x) x$goric[1]))
  df    <- data.frame(loglik = ll, penalty = PT, goric)
  
  delta <- df$goric - min(df$goric)
  goric_weights <- exp(-delta / 2) / sum(exp(-delta / 2))
  df$goric_weights <- goric_weights
  
  df
}
