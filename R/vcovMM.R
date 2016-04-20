#Adjusted code from the .vcov.avar1() function in the Robustbase package
#Computes the robust variance-covariance matrix
vcovMM <- function(X, resid0, resid, scale, ...) { 
  x <- X#as.matrix(obj$x)
  ## set psi and chi constants
  psi <- chi <- "bisquare"
  c.chi <- 1.54764
  c.psi <- 4.685061
  
  r0 <- resid0#as.numeric(obj$resid0)
  r <- as.numeric(resid)#as.numeric(resid(obj))
  scale <- scale#obj$s
  bb <- 1/2 ## this is always 1/2 for S estimates by convention
  ## scaled residuals
  n <- length(r)
  p <- ncol(x)
  r.s   <- r / scale # final   scaled residuals
  r0.s <- r0 / scale # initial scaled residuals
  w  <- robustbase:::Mpsi(r.s, cc = c.psi, psi = psi, deriv = 1)
  w0 <- robustbase:::Mchi(r0.s, cc = c.chi, psi = chi, deriv = 1)
  x.wx <- crossprod(x, x * w)
  if(inherits(A <- tryCatch(solve(x.wx) * scale,
                            error=function(e)e), "error")) {
    warning("X'WX is almost singular. Consider rather using cov = \".vcov.w\"")
    A <- tryCatch(solve(x.wx, tol = 0) * scale, error=function(e)e)
    if(inherits(A, "error"))
      stop("X'WX is singular. Rather use cov = \".vcov.w\"")
  }
  a <- A %*% (crossprod(x, w * r.s) / mean(w0 * r0.s))
  w <- robustbase:::Mpsi( r.s, cc = c.psi, psi = psi)
  
  ## 3) now the standard part  (w, x, r0.s,  n, A,a, c.chi, bb)
  w0 <- robustbase:::Mchi(r0.s, cc = c.chi, psi = chi)
  Xww <- crossprod(x, w*w0)
  u1 <- A %*% crossprod(x, x * w^2) %*% (n * A)
  u2 <- a %*% crossprod(Xww, A)
  u3 <- A %*% tcrossprod(Xww, a)
  u4 <- mean(w0^2 - bb^2) * tcrossprod(a)
  
  ## list(cov = matrix((u1 - u2 - u3 + u4)/n, p, p),
  ret <- (u1 - u2 - u3 + u4)/n
  
  ret
}
