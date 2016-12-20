## utility functions
coef.conLM <- function(object, ...)  {
  
  b_def <- c()
  b_restr <- object$b.restr
  if (any(object$parTable$op == ":=")) {
    b_def <- object$CON$def.function(object$b.restr)
  }
  
  OUT <- c(b_restr, b_def)
  
  OUT
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