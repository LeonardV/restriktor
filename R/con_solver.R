con_solver <- function(b.unconstr, X, y, Amat, bvec, meq,
                       maxit = 10000, tol = sqrt(.Machine$double.eps)) {
  val <- 0
  y <- as.matrix(y)
  for (i in 1:maxit) {
    Sigma <- (t(y - X%*%matrix(b.unconstr, ncol = ncol(y))) %*%
                (y - X%*%matrix(b.unconstr, ncol = ncol(y)))) / (nrow(X)) #ML!
    yVx <- kronecker(solve(Sigma), t(X)) %*% as.vector(y)
    dvec <- 2*yVx
    Dmat <- 2*kronecker(solve(Sigma), t(X) %*% X)
    out <- con_my_solve_QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat),
                           bvec = bvec, meq = meq)

    if(out$status == -1) {
      stop("constraints are inconsistent, no solution!")
    } else if (out$status == -2) {
      stop("matrix D in quadratic function is not positive definite!")
    }

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


