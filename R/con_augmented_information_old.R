con_augmented_information_old <- function(X = NULL, beta = NULL, 
                                      information = NULL,
                                      Amat = NULL, bvec = NULL,
                                      inverted    = TRUE,
                                      check.pd    = FALSE) {
  
  H <- Amat
  npar <- nrow(information)
  is.augmented <- FALSE
  
  inactive.idx <- H %*% beta - bvec >= 0 * bvec
#  H.active <- H[!inactive.idx,,drop=FALSE]
#  H.inactive <- H[inactive.idx,,drop=FALSE]
  
  # lagrangean coefs
  lambda <- as.numeric(solve(H%*%solve(t(X)%*%X)%*%t(H)) %*% (H%*%beta-bvec))
  
  # handle constraints
  if(nrow(H) > 0L) {
    #lambda <- lavmodel@con.lambda # lagrangean coefs
    if(length(inactive.idx) > 0L) {
      H <- H[-inactive.idx,,drop=FALSE]
      lambda <- lambda[-inactive.idx]
    }
    if(nrow(H) > 0L) {
      is.augmented <- TRUE   
      H0 <- matrix(0,nrow(H),nrow(H))
      H10 <- matrix(0, ncol(information), nrow(H))
      DL <- 2*diag(lambda, nrow(H), nrow(H))
      E3 <- rbind( cbind(     information,  H10, t(H)),
                   cbind(          t(H10),   DL,  H0),
                   cbind(               H,   H0,  H0)  )
      information <- E3   
    }
  }
  
  if(check.pd) {
    eigvals <- eigen(information, symmetric = TRUE, 
                     only.values = TRUE)$values
    if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
      warning("lavaan WARNING: matrix based on first order outer product of the derivatives is not positive definite; the model may not be identified")
    }
  }
  
  if(inverted) {
    if(is.augmented) {
      # note: default tol in MASS::ginv is sqrt(.Machine$double.eps)
      #       which seems a bit too conservative
      #       from 0.5-20, we changed this to .Machine$double.eps^(3/4)
      information <- 
        try( MASS::ginv(information, 
                        tol = .Machine$double.eps^(3/4))[1:npar, 
                                                         1:npar, 
                                                         drop = FALSE],
             silent = TRUE )
    } else {
      information <- try( solve(information), silent = TRUE )
    }
  }
  information[abs(information) < sqrt(.Machine$double.eps)] <- 0L
  
  # augmented/inverted information
  information
}
