con_solver <- function(b.unconstr, X, y, w, Amat, bvec, meq,
                       maxit = 10000, tol = sqrt(.Machine$double.eps)) {
  val <- 0
  y <- as.matrix(y)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  if (is.null(w)) {
    w <- rep(1, n)
  }
  W <- diag(w)
  
  for (i in 1:maxit) {
    Sigma <- (t(y - X%*%matrix(b.unconstr, ncol = ncol(y))) %*% W %*% 
                (y - X%*%matrix(b.unconstr, ncol = ncol(y)))) / (n-p)           #ML divided by n
    yVx <- kronecker(solve(Sigma), t(X)) %*% as.vector(y)
    dvec <- 2*yVx
    Dmat <- 2*kronecker(solve(Sigma), t(X) %*% X)
    out <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat),
                           bvec = bvec, meq = meq)

#     if(out$status == -1) {
#       stop("constraints are inconsistent, no solution!")
#     } else if (out$status == -2) {
#       stop("matrix D in quadratic function is not positive definite!")
#     }

    if (abs(out$value - val) <= tol) {
      break
    } else {
      val <- out$value
    }
    if (i == maxit & abs(out$value - val) > tol)
      warning("Maximum number of iterations reached without convergence.")
  }

  OUT <- out

  OUT
}


