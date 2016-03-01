con_augmented_information <- function(model, s2, type, b.constr, constraints, 
                                      bvec, meq) {

  X <- model.matrix(model)[,,drop = FALSE]
  b.unconstr <- coef(model)
  H <- constraints
  information <- 1/fit.con$s2 * crossprod(X) #1/fit.con$s2 * crossprod(X) vcov(model)
  npar <- ncol(information)
  
  #information <- 
  #solve(vcovHC(model, type = type, adjust = TRUE))
  
#  if (!isSymmetric(information)) {
#    stop("Information matrix information is not symmetric.")
#  }

  # lagrangean coefs
  lambda <- as.numeric(solve(H%*%solve(t(X)%*%X)%*%t(H)) %*% (H%*%b.unconstr-bvec))

  inactive.idx <- which(H %*% b.unconstr - bvec >= 0 * bvec)
#  active.idx <- !inactive.idx
  
#  if (meq > 0L) {
#    inactive.idx[1:meq] <- FALSE
#  }
#  if (meq > 0L) {
#    active.idx[1:meq] <- FALSE
#  }
  #equality constraints
#  if (meq > 0L) {
#    G <- rbind(H[1:meq,])
#  } else {
#    G <- matrix(0L, 1, nrow(information))
#  }
  #active inequality constraints
#  H.H <- rbind(H[active.idx,])
  #inactive inequality constraints
#  H.IA <- rbind(H[inactive.idx,])
  # diagonal matrix with Lagrangean coefficients
#  Gamma <- diag(L, nrow=nrow(H), ncol=nrow(H))
  # slack parameters
#  slacks <- H %*% b.constr - bvec
#    slacks[abs(slacks) < sqrt(.Machine$double.eps)] <- 0L
  # diagonal matrix with slack parameters for the inactive constraints
#  Z <-
#    if (sum(inactive.idx) == 1L) {
#      as.matrix(slacks[inactive.idx])
#  } else {
#    diag(slacks[inactive.idx])
#  }

  #H <- H[!inactive.idx,,drop=FALSE]
  #L <- L[!inactive.idx]
  #H0 <- matrix(0,nrow(H),nrow(H))
#  H12 <- matrix(0L, nrow(information), ncol(Gamma))
#  H13 <- matrix(0L, nrow(information), ncol(Z))
#  H23 <- matrix(0L, nrow(Gamma), ncol(Z))
#  H24 <- matrix(0L, nrow(Gamma), nrow(G))
#  H25 <- matrix(0L, nrow(Gamma), nrow(H.H))
#  H26 <- matrix(0L, nrow(Gamma), nrow(H.IH))
#  H33 <- matrix(0L, nrow(Z), ncol(Z))
#  H34 <- matrix(0L, nrow(Z), nrow(G))
#  H35 <- matrix(0L, nrow(Z), nrow(H.H))
#  H44 <- matrix(0L, nrow(G), nrow(G))
#  H45 <- matrix(0L, nrow(G), nrow(H.H))
#  H46 <- matrix(0L, nrow(G), nrow(H.IH))
#  H55 <- matrix(0L, nrow(H.H), nrow(H.H))
#  H56 <- matrix(0L, nrow(H.H), nrow(H.IH))
#  H66 <- matrix(0L, nrow(H.IH), nrow(H.IH))
#  DL <- 2*diag(L, nrow(H), nrow(H))

##############################################################  
  
  if(length(inactive.idx) > 0L) {
    H <- H[-inactive.idx,,drop=FALSE]
    lambda <- lambda[-inactive.idx]
  }
  if(nrow(H) > 0L) {
    H0 <- matrix(0,nrow(H),nrow(H))
    H10 <- matrix(0, ncol(information), nrow(H))
    DL <- 2*diag(lambda, nrow(H), nrow(H))
    E3 <- rbind( cbind(     information,  H10, t(H)),
                 cbind(          t(H10),   DL,  H0),
                 cbind(               H,   H0,  H0)  )
    information <- E3   
  }
########################################################################3
  # block matrix - Ronald Schoenberg (1997)
#  E3 <- rbind( cbind(  information,    H12,    H13,   t(G), t(H.H), t(H.IH)),
#               cbind( t(H12),     DL,    H23,    H24,    H25,     H26),
#               cbind( t(H13), t(H23),    H33,    H34,    H35,     2*Z),
#               cbind(      G, t(H24), t(H34),    H44,    H45,     H46),
#               cbind(    H.H, t(H25), t(H35), t(H45),    H55,     H56),
#               cbind(   H.IH, t(H26),    2*Z, t(H46),  t(H56),     H66)
#  )

  information <- try( MASS::ginv(information)[1:npar, 1:npar, drop = FALSE], silent = TRUE )
  information[abs(information) < sqrt(.Machine$double.eps)] <- 0L

  OUT <- information

  OUT
}

