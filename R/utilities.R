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
  
  OUT
}

model.matrix.restriktor <- function(object, ...) {
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
