# theses are originally a function from the lavaan package.
lav_constraints_linear_idx <- function(func = NULL, npar = NULL) {
  if (is.null(func) || is.null(body(func))) 
    return(integer(0L))
  A0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))
  A1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))
  A0minA1 <- A0 - A1
  linear <- apply(A0minA1, 1, function(x) all(x == 0))
  which(linear)
}

lav_constraints_nonlinear_idx <- function (func = NULL, npar = NULL) {
  if (is.null(func) || is.null(body(func))) 
    return(integer(0L))
  A0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))
  A1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))
  A0minA1 <- A0 - A1
  linear <- apply(A0minA1, 1, function(x) all(x == 0))
  which(!linear)
}