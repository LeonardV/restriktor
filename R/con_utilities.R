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


goric <- function(object, ..., digits = max(3, getOption("digits") - 2)) {
  if (inherits(object, "conLM")) {
    objectlist <- list(object, ...)
  } else {
    objectlist <- object
  }
  isconLM <- sapply(objectlist, function(x) inherits(x, "conLM"))
  conlist <- objectlist[isconLM]  
  isSummary <- lapply(conlist, function(x) summary(x))
  
  ll    <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik")))
  PT    <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
  goric <- unlist(lapply(isSummary, function(x) x$goric[1]))
  df    <- data.frame(loglik = ll, penalty = PT, goric)
  
  delta <- df$goric - min(df$goric)
  goric_weights <- exp(-delta / 2) / sum(exp(-delta / 2))
  df$goric_weights <- goric_weights
  
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  
  invisible(df)
}
