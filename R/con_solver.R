#REF: Zheng, S., Guo, J. Shi, N.Z, Tian, G.L (2012). Likelihood-based approaches 
# formultivariate linear models under inequality constraints for incomplete data.
# Journal of Statistical Planning and Inference 142, 2926-2942.
# Slightly adjusted by LV to deal with weigths.
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
                (y - X%*%matrix(b.unconstr, ncol = ncol(y)))) / (n - p)           #ML divided by n
    
    yVx <- kronecker(solve(Sigma), t(X)) %*% W %*% as.vector(y)
    dvec <- 2*yVx
    Dmat <- 2*kronecker(solve(Sigma), t(X) %*% W %*% X)
    out <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat),
                           bvec = bvec, meq = meq)

    if (abs(out$value - val) <= tol) {
      break
    } else {
      val <- out$value
    }
    
    if (i == maxit & abs(out$value - val) > tol) {
      warning(gettextf("'quadprog' failed to converge in %d steps", maxit), 
              domain = NA)
    }  
  }

  OUT <- out

  OUT
}


