con_augmented_hessian <- function(b.unconstr, b.constr, X, s2, constraints, bvec, meq) {

  A <- constraints
  Sigma <- 1/s2 * (t(X)%*%X)

  if(!isSymmetric(Sigma)) {
    stop("Information matrix Sigma is not symmetric.")
  }

  # lagrangean coefs
  L <- as.numeric(solve(A%*%solve(t(X)%*%X)%*%t(A)) %*% (A%*%b.unconstr-bvec))

  inactive.idx <- A %*% b.unconstr - bvec >= 0 * bvec
  if(meq > 0L) {
    inactive.idx[1:meq] <- FALSE
  }
  active.idx <- !inactive.idx
  if(meq > 0L) {
    active.idx[1:meq] <- FALSE
  }
  #equality constraints
  if(meq > 0L) {
    G <- rbind(A[1:meq,])
  } else {
    G <- matrix(0L, 1, nrow(Sigma))
  }
  #active inequality constraints
  H.A <- rbind(A[active.idx,])
  #inactive inequality constraints
  H.IA <- rbind(A[inactive.idx,])
  # diagonal matrix with Lagrangean coefficients
  Gamma <- diag(L, nrow=nrow(A), ncol=nrow(A))
  # slack parameters
  slacks <- A %*% b.constr - bvec
    slacks[abs(slacks) < sqrt(.Machine$double.eps)] <- 0L
  # diagonal matrix with slack parameters for the inactive constraints
  Z <-
    if(sum(inactive.idx) == 1L) {
      as.matrix(slacks[inactive.idx])
  } else {
      diag(slacks[inactive.idx])
  }

  #A <- A[!inactive.idx,,drop=FALSE]
  #L <- L[!inactive.idx]
  #H0 <- matrix(0,nrow(A),nrow(A))
  H12 <- matrix(0L, nrow(Sigma), ncol(Gamma))
  H13 <- matrix(0L, nrow(Sigma), ncol(Z))
  H23 <- matrix(0L, nrow(Gamma), ncol(Z))
  H24 <- matrix(0L, nrow(Gamma), nrow(G))
  H25 <- matrix(0L, nrow(Gamma), nrow(H.A))
  H26 <- matrix(0L, nrow(Gamma), nrow(H.IA))
  H33 <- matrix(0L, nrow(Z), ncol(Z))
  H34 <- matrix(0L, nrow(Z), nrow(G))
  H35 <- matrix(0L, nrow(Z), nrow(H.A))
  H44 <- matrix(0L, nrow(G), nrow(G))
  H45 <- matrix(0L, nrow(G), nrow(H.A))
  H46 <- matrix(0L, nrow(G), nrow(H.IA))
  H55 <- matrix(0L, nrow(H.A), nrow(H.A))
  H56 <- matrix(0L, nrow(H.A), nrow(H.IA))
  H66 <- matrix(0L, nrow(H.IA), nrow(H.IA))
  DL <- 2*diag(L, nrow(A), nrow(A))

#  E3 <- rbind( cbind(  Sigma,  H10, t(A)),
#               cbind( t(H10),   DL,  H0),
#               cbind(      A,   H0,  H0)  )

  # block matrix - Ronald Schoenberg (1997)
  E3 <- rbind( cbind(  Sigma,    H12,    H13,   t(G), t(H.A), t(H.IA)),
               cbind( t(H12),     DL,    H23,    H24,    H25,     H26),
               cbind( t(H13), t(H23),    H33,    H34,    H35,     2*Z),
               cbind(      G, t(H24), t(H34),    H44,    H45,     H46),
               cbind(    H.A, t(H25), t(H35), t(H45),    H55,     H56),
               cbind(   H.IA, t(H26),    2*Z, t(H46), t(H56),     H66)
  )


  #augmented hessian matrix
  p <- ncol(Sigma)
  aug.H <- MASS::ginv(E3)[1:p,1:p, drop = FALSE]
    aug.H[abs(aug.H) < sqrt(.Machine$double.eps)] <- 0L

  OUT <- aug.H

  OUT
}

