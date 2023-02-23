#REF: Zheng, S., Guo, J. Shi, N.Z, Tian, G.L (2012). Likelihood-based approaches 
# for multivariate linear models under inequality constraints for incomplete data.
# Journal of Statistical Planning and Inference 142, 2926-2942.
con_solver_lm <- function(X, y, w = NULL, Amat, bvec, meq,
                          maxit = 10000, absval = sqrt(.Machine$double.eps)) {
  val <- 0
  X <- as.matrix(X)
  y <- as.matrix(y)
  b.restr <- c(coef(lm.fit(x = X, y)))
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  for (i in 1:maxit) {
    if (!is.null(w)) {
      # weights  
      W <- diag(w)
      # dividing by n or (n-p-rank(meq)). Probably only needed for mlm.
      s2 <- (t(y - X %*% matrix(b.restr, ncol = ncol(y))) %*% W %*%
                  (y - X %*% matrix(b.restr, ncol = ncol(y)))) / (n - p)
      yVx <- kronecker(solve(s2), t(X)) %*% W %*% as.vector(y)  
      Dmat <- 2 * kronecker(solve(s2), t(X) %*% W %*% X)
    } else {
      # no weights
      s2 <- (t(y - X %*% matrix(b.restr, ncol = ncol(y))) %*%
                  (y - X %*% matrix(b.restr, ncol = ncol(y)))) / (n - p)  
      yVx <- kronecker(solve(s2), t(X)) %*% as.vector(y)  
      Dmat <- 2 * kronecker(solve(s2), t(X) %*% X)
    }
    
    dvec <- 2 * yVx
    
    out.qp <- solve.QP(Dmat = Dmat, 
                       dvec = dvec, 
                       Amat = t(Amat),
                       bvec = bvec, 
                       meq  = meq)
    
    b.restr <- out.qp$solution

    if (abs(out.qp$value - val) <= absval) {
      break
    } else {
      val <- out.qp$value
    }
    
    if (i == maxit & abs(out.qp$value - val) > absval) {
      warning(gettextf("'quadprog' failed to converge in %d steps", maxit), 
              domain = NA)
    }  
  }

  out <- list(qp = out.qp, s2 = s2, Niter = i)
  
  out
}



con_solver_rlm <- function(X, y, Amat, bvec, meq,
                           maxit = 10000, absval = sqrt(.Machine$double.eps)) {
  
  b.unrestr <- lm.fit(x = X, y)
  b.restr <- as.vector(coef(b.unrestr))
  invW <- crossprod(X)
  Dmat <- 2 * invW
  iact <- which(Amat %*% b.restr - bvec < 0)
  
  val <- 0 
  for (i in 1:maxit) {
    dvec <- 2 * b.restr %*% invW
    out.qp <- solve.QP(Dmat = Dmat, 
                       dvec = dvec, 
                       Amat = t(Amat),
                       bvec = bvec, 
                       meq  = meq)
      
    if (abs(out.qp$value - val) <= absval) {
      break
    } else {
      val <- out.qp$value
    }
    
    b.restr <- out.qp$solution
    
    if (i == maxit & abs(out.qp$value - val) > absval) {
      warning(gettextf("'quadprog' failed to converge in %d steps", maxit), 
              domain = NA)
    }  
  }
  
  out <- out.qp
  out$iact <- iact
  
  out
}



con_solver_glm <- function(X, y, Amat, bvec, meq, 
                           maxit = 10000, epsilon){
  
  b.restr <- lm.fit(x = X, y)
  b.restr <- as.vector(coef(b.restr))
  invW <- crossprod(X)
  
  # con_solver <- function(b.restr, invW, Amat, bvec, meq) {
  #   Dmat <- 2 * invW
  #   dvec <- 2 * b.restr %*% invW
  #   Amat <- t(Amat)
  #   solve.QP(Dmat, dvec, Amat = Amat, bvec = bvec, meq = meq)
  # }
  
  val <- 0
  for (i in 1:maxit) {
    #out <- con_solver(b.restr, invW, Amat, bvec, meq)
    Dmat <- 2 * invW
    dvec <- 2 * b.restr %*% invW
    
    out.qp <- solve.QP(Dmat = Dmat, 
                       dvec = dvec, 
                       Amat = t(Amat),
                       bvec = bvec, 
                       meq  = meq)
    
    b.restr <- out.qp$solution
    
    if (abs(out.qp$value - val) <= epsilon) {
      break
    } else {
      val <- out.qp$value
    }
    
    if (i == maxit & abs(out.qp$value - val) > epsilon) {
      warning(gettextf("'quadprog' failed to converge in %d steps", maxit), 
              domain = NA)
    }
  }
  
  out <- out.qp
  
  out
}


# compute restricted estimates based on unrestricted estimates and vcov
con_solver_gorica <- function(est, VCOV, Amat, bvec, meq) {
    
    VCOV <- as.matrix(VCOV)
    invVCOV <- solve(VCOV)
    #invVCOV <- ginv(VCOV)
    #invVCOV <- chol2inv(chol(VCOV))
    
    Dmat <- 2 * invVCOV 
    dvec <- 2 * (est %*% invVCOV)
    
    out.qp <- solve.QP(Dmat = Dmat, 
                       dvec = dvec, 
                       Amat = t(Amat),
                       bvec = bvec, 
                       meq  = meq)
    out <- out.qp
    
    out
}


