# Ref: Ronald Schoenberg (1997). Constrained Maximum Likelihood. 
# Computational Economics, 10, 251-266
con_augmented_information <- function(information  = NULL, 
                                      is.augmented = TRUE,
                                      X            = NULL, 
                                      b.unrestr    = NULL, 
                                      b.restr      = NULL,  
                                      Amat         = NULL, 
                                      bvec         = NULL, 
                                      meq          = NULL) {
  
  H <- Amat
  npar <- NCOL(information)
  OUT <- list()
  
  if (!isSymmetric(information)) {
    stop("Information matrix is not symmetric.")
  }
  
  if (!all(c(H) == 0L)) { 
    # lagrangean coefs
    lambda <- as.numeric(solve(H%*% solve(t(X)%*%X)%*%t(H)) %*% 
                           (H%*%b.unrestr-bvec))
    
    #equality constraints
    if (meq > 0L) {
      G <- rbind(H[1:meq, ])
    } else {
      G <- matrix(0L, 1, NROW(information))
    }
    
    if (meq > 0L) {
      H <- H[-c(1:meq),,drop=FALSE]
      bvec <- bvec[-c(1:meq)]
    }
    
    #inactive inequality constraints
    inactive.idx <- H %*% b.restr - bvec > 0 * bvec
    #active inequality constraints
    H.active <- H[!inactive.idx,,drop=FALSE]
    #inactive inequality constraints
    H.inactive <- H[inactive.idx,,drop=FALSE]
    # diagonal matrix with Lagrangean coefficients
    Gamma <- diag(lambda, NROW(H), NROW(H))
    # slack parameters
    slacks <- H %*% b.restr - bvec
    slacks[abs(slacks) < sqrt(.Machine$double.eps)] <- 0L
    # diagonal matrix with slack parameters for the inactive constraints
    Z <- matrix(0L, nrow = 0, ncol = 0)
    if (sum(inactive.idx) == 1L) {
      Z <- diag(slacks[inactive.idx,,drop=FALSE])  
    }
    if (sum(inactive.idx) > 1L) {
      Z <- diag(slacks[inactive.idx,,drop=TRUE])  
    }  
    
    H12 <- matrix(0L, NROW(information), NCOL(Gamma))
    H13 <- matrix(0L, NROW(information), NCOL(Z))
    H23 <- matrix(0L, NROW(Gamma), NCOL(Z))
    H24 <- matrix(0L, NROW(Gamma), NROW(G))
    H25 <- matrix(0L, NROW(Gamma), NROW(H.active))
    H26 <- matrix(0L, NROW(Gamma), NROW(H.inactive))
    H33 <- matrix(0L, NROW(Z), NCOL(Z))
    H34 <- matrix(0L, NROW(Z), NROW(G))
    H35 <- matrix(0L, NROW(Z), NROW(H.active))
    H44 <- matrix(0L, NROW(G), NROW(G))
    H45 <- matrix(0L, NROW(G), NROW(H.active))
    H46 <- matrix(0L, NROW(G), NROW(H.inactive))
    H55 <- matrix(0L, NROW(H.active), NROW(H.active))
    H56 <- matrix(0L, NROW(H.active), NROW(H.inactive))
    H66 <- matrix(0L, NROW(H.inactive), NROW(H.inactive))
    DL <- 2*diag(lambda, NROW(H), NROW(H))
    
    # block matrix 
    E3 <- rbind( cbind(  information,    H12,    H13,   t(G), t(H.active), t(H.inactive)),
                 cbind(       t(H12),     DL,    H23,    H24,         H25,           H26),
                 cbind(       t(H13), t(H23),    H33,    H34,         H35,           2*Z),
                 cbind(            G, t(H24), t(H34),    H44,         H45,           H46),
                 cbind(     H.active, t(H25), t(H35), t(H45),         H55,           H56),
                 cbind(   H.inactive, t(H26),    2*Z, t(H46),      t(H56),           H66)
               )
  }
  
  if (is.augmented) {
    OUT$information <- try( MASS::ginv(E3, tol = .Machine$double.eps^(3/4))[1:npar, 
                                                                            1:npar, 
                                                                            drop = FALSE], 
                            silent = TRUE )
    OUT$information[abs(OUT$information) < sqrt(.Machine$double.eps)] <- 0L
    OUT$information.augmented <- E3
  } else { 
    OUT$information <- try( solve(information), silent = TRUE )
  }
  
  OUT
}

