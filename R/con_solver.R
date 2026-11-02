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

    diff <- abs(out.qp$value - val)
    val <- out.qp$value

    if (diff <= absval) {
      break
    }

    if (i == maxit && diff > absval) {
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
      
    b.restr <- out.qp$solution

    diff <- abs(out.qp$value - val)
    val <- out.qp$value

    if (diff <= absval) {
      break
    }

    if (i == maxit && diff > absval) {
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

    diff <- abs(out.qp$value - val)
    val <- out.qp$value

    if (diff <= epsilon) {
      break
    }

    if (i == maxit && diff > epsilon) {
      warning(gettextf("'quadprog' failed to converge in %d steps", maxit),
              domain = NA)
    }
  }

  out <- out.qp

  out
}


# compute restricted estimates based on unrestricted estimates and vcov
#
# The QP is solved on a rescaled problem, because solve.QP() can wrongly report
# "constraints are inconsistent, no solution!" when the parameters are on very
# different scales (e.g. a quadratic time effect with time in days, giving
# variances of order 1e-9 next to order 1). Rescaling does not change the
# problem:
#  - parameters are standardised by their standard errors, s = sqrt(diag(VCOV)):
#      theta = S theta_s, with S = diag(s), so VCOV becomes the correlation
#      matrix R = S^-1 VCOV S^-1 and Amat becomes Amat S;
#  - each constraint row (and its rhs) is divided by its Euclidean norm.
# The solution is transformed back (theta = S theta_s); the objective value,
# the active set (iact) and the feasibility region are unchanged. The Lagrange
# multipliers are transformed back to the original constraint scale.
con_solver_gorica <- function(est, VCOV, Amat, bvec, meq) {

    VCOV <- as.matrix(VCOV)
    Amat <- as.matrix(Amat)
    est  <- c(est)

    # 1. standardise the parameters
    s <- sqrt(diag(VCOV))
    if (any(!is.finite(s)) || any(s <= 0)) {
      s <- rep(1, length(est)) # cannot standardise; solve the original problem
    }
    R      <- VCOV / tcrossprod(s)        # = S^-1 VCOV S^-1
    est_s  <- est / s
    Amat_s <- sweep(Amat, 2, s, `*`)      # = Amat S

    # 2. normalise the constraint rows (and their rhs)
    rnorm_A <- sqrt(rowSums(Amat_s^2))
    rnorm_A[rnorm_A == 0] <- 1
    Amat_s <- Amat_s / rnorm_A
    bvec_s <- bvec / rnorm_A

    invR <- solve(R)
    #invR <- ginv(R)
    #invR <- chol2inv(chol(R))

    Dmat <- 2 * invR
    dvec <- 2 * (est_s %*% invR)

    out.qp <- solve.QP(Dmat = Dmat,
                       dvec = dvec,
                       Amat = t(Amat_s),
                       bvec = bvec_s,
                       meq  = meq)

    # 3. back-transform to the original parameter / constraint scale
    out <- out.qp
    out$solution               <- c(out.qp$solution * s)
    out$unconstrained.solution <- c(out.qp$unconstrained.solution * s)
    out$Lagrangian             <- out.qp$Lagrangian / rnorm_A

    out
}
