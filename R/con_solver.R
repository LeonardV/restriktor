#REF: Zheng, S., Guo, J. Shi, N.Z, Tian, G.L (2012). Likelihood-based approaches 
# for multivariate linear models under inequality constraints for incomplete data.
# Journal of Statistical Planning and Inference 142, 2926-2942.
con_solver_lm <- function(X, y, b.unrestr, w, Amat, bvec, meq,
                          maxit = 10000, absval = sqrt(.Machine$double.eps)) {
  val <- 0
  y <- as.matrix(y)
  n <- dim(X)[1]
  
  if (is.null(w)) {
    w <- rep(1, n)
  }
  W <- diag(w)
  
  for (i in 1:maxit) {
    # dividing by n or (n-p-rank(meq)). Probably only needed for mlm.
    Sigma <- (t(y - X %*% matrix(b.unrestr, ncol = ncol(y))) %*% W %*%
                (y - X %*% matrix(b.unrestr, ncol = ncol(y)))) / n  
    
    yVx <- kronecker(solve(Sigma), t(X)) %*% W %*% as.vector(y)
    dvec <- 2 * yVx
    Dmat <- 2 * kronecker(solve(Sigma), t(X) %*% W %*% X)
    out <- solve.QP(Dmat = Dmat, 
                    dvec = dvec, 
                    Amat = t(Amat),
                    bvec = bvec, 
                    meq  = meq)

    if (abs(out$value - val) <= absval) {
      break
    } else {
      val <- out$value
    }
    
    if (i == maxit & abs(out$value - val) > absval) {
      warning(gettextf("'quadprog' failed to converge in %d steps", maxit), 
              domain = NA)
    }  
  }

  out
}


con_solver_rlm <- function(X, y, Amat, bvec, meq,
                           maxit = 10000, absval = sqrt(.Machine$double.eps)) {
  
  b.unrestr <- lm.fit(x = X, y)
  tBeta <- as.vector(coefficients(b.unrestr))
  invW <- crossprod(X)
  
  val <- 0
  for (i in 1:maxit) {
    Dmat <- 2 * invW
    dvec <- 2 * tBeta %*% invW
    out <- solve.QP(Dmat = Dmat, 
                    dvec = dvec, 
                    Amat = t(Amat),
                    bvec = bvec, 
                    meq  = meq)
    
    if (abs(out$value - val) <= absval) {
      break
    } else {
      val <- out$value
    }
    
    if (i == maxit & abs(out$value - val) > absval) {
      warning(gettextf("'quadprog' failed to converge in %d steps", maxit), 
              domain = NA)
    }  
  }
  
  out
}



con_solver_glm <- function(X, y, b.unrestr, Amat, bvec, meq, 
                           maxit = 10000, epsilon){
  b.unrestr <- lm.fit(x = X, y)
  tBeta <- as.vector(coefficients(b.unrestr))
  invW <- crossprod(X)#t(X) %*% X
  
  con_solver <- function(tBeta, invW, Amat, bvec, meq) {
    Dmat <- 2 * invW
    dvec <- 2 * tBeta %*% invW
    Amat <- t(Amat)
    solve.QP(Dmat, dvec, Amat = Amat, bvec = bvec, meq = meq)
  }
  
  b.restr <- tBeta
  val <- 0
  for (i in 1:maxit) {
    out <- con_solver(b.restr, invW, Amat, bvec, meq)
    b.restr <- out$solution
    if (abs(out$value - val) <= epsilon) {
      break
    } else {
      val <- out$value
    }
    
    if (i == maxit & abs(out$value - val) > epsilon) {
      warning(gettextf("'quadprog' failed to converge in %d steps", maxit), 
              domain = NA)
    }
  }
  
  out
}
