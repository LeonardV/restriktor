robustWaldScores <- function(x, y, beta0, betaA, scale, Amat = NULL, bvec = NULL, 
                             meq = NULL) { 
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  #tukey bisquare tuning constant
  cc = 4.685061
  
  #Calculate M, Q, V 
  res0 <- y - x%*%beta0
  resA <- y - x%*%betaA
  
  rstar0 <- res0/scale
  rstar2 <- resA/scale
  
  psi0   <- tukeyChi(rstar0, cc, deriv=1)  
  psi1   <- tukeyChi(rstar2, cc, deriv=1) 
  psideriv0 <- tukeyChi(rstar0, cc, deriv=2) 
  psideriv1 <- tukeyChi(rstar2, cc, deriv=2) 
  
  #compute M 
  weightsM <- psideriv1 / scale             
  WM <- weightsM %*% rep(1, p)
  xwM <- x * WM
  M <- t(x) %*% xwM / n
  
  #compute Q 
  weightsQ <- psi1^2
  WQ <- weightsQ %*% rep(1, p)
  xwQ <- x * WQ
  Q <- t(x) %*% xwQ / n
  
  if (!is.null(Amat)) {
    idx1 <- which(colSums(abs(Amat)) > 0L)
    idx0 <- which(colSums(abs(Amat)) == 0L)
  }
  
  #Calculate V 
  Minv <- solve(M)
  #information matrix 
  V <- Minv %*% Q %*% t(Minv)
  
  # V[(p0+1):pa,(p0+1):pa]
  V22 <- V[idx1,idx1]
  
  #inverse information matrix
  #V.inv <- solve(V)
  #submatrix of V
  V22.inv <- solve(V22)
  
  #Schur complement?
  M221 <- M[idx1,idx1] - M[idx1,idx0] %*% solve(M[idx0,idx0,drop=FALSE], 
                                                M[idx0,idx1,drop=FALSE])
  
  #  M221 <- M[(p0+1):p,(p0+1):p] - M[(p0+1):p,1:p0] %*% solve(M[1:p0,1:p0,drop=FALSE], M[1:p0,(p0+1):p,drop=FALSE])
  
  weightsZ <- psi0
  Z <- t(x) %*% weightsZ / n  
  
  #Wald-type test statistic, Silvapulle (1996, eq. 2.6)
  if (!is.null(Amat)) {
    Tn <- sqrt(n)*(Amat%*%betaA)
    Dn <- Tn
    Dmat <- solve(Amat[,idx1,drop=FALSE]%*%V22%*%t(Amat[,idx1,drop=FALSE]))
    dvec <- t(Dn)%*%Dmat
    out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(diag(length(dvec))), bvec=bvec, meq=meq) 
    b <- out$solution
    TsWald  <- as.numeric((t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b)) )
    #TS.WB  <- t(Dn-b)%*%Dmat%*%(Dn-b) 
  } else {
    #unconstrained
    TsWald <- as.numeric(n * betaA[idx1] %*% solve(V22, betaA[idx1]) )
  }
  
  if (!is.null(Amat)) {
    #Score-type test (Silvapulle, 1996, eq. 2.6)
    An <- sqrt(n) * Amat[,idx1,drop=FALSE]%*%(solve(M221) %*% Z[idx1])
    Dn <- An
    Dmat <- solve(Amat[,idx1,drop=FALSE]%*%V22%*%t(Amat[,idx1,drop=FALSE]))
    dvec <- t(Dn)%*%Dmat
    out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(diag(length(dvec))), bvec=bvec, meq=meq) 
    b <- out$solution
    TsScore <- as.numeric((t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b)))
    #TS.SB <- (t(Dn-b)%*%Dmat%*%(Dn-b))
  } else {
    #unconstrained
    result.C <- M221 %*% V22 %*% t(M221)
    TsScore <- as.numeric(n * t(Z[idx1]) %*% solve(result.C, Z[idx1]))
  }
  
  OUT <- list(RWald = TsWald,
              Rscore = TsScore,
              V = V,
              V22 = Dmat)
  
  OUT
  
}



## robust Fm test statistic ##
robustFm <- function(x, y, beta0, betaA, scale, cc = 4.685061) {
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  #compute residuals under null and alternative model
  resid0 <- y - x%*%beta0
  resid1 <- y - x%*%betaA
  
  #residuals / scale
  rstar0 <- c(resid0/scale)                                               
  rstar1 <- c(resid1/scale)
  
  L0 <- sum(tukeyChi(rstar0, cc, deriv = 0))
  L1 <- sum(tukeyChi(rstar1, cc, deriv = 0))
  
  #first derivative psi function
  psi.prime.h1 <- tukeyChi(rstar1, cc, deriv = 1)                     
  #second derivative psi function
  psi.prime2.h1 <- tukeyChi(rstar1, cc, deriv = 2)                    
  
  #asymptotic covariance matrix standardizing constant
  l.h1 <- ( 0.5 * (1/(n - p)) * sum(psi.prime.h1^2) ) / ( (1/n)*sum(psi.prime2.h1) )
  out <- 1/l.h1 * (L0 - L1) 
    
  out
}


## robust Wald statistic, Silvapulle (1992) ##
robustWaldXX <- function(x, beta0, beta1, beta2, tau) {
  TsWald <- c(( (t(beta2-beta0)%*%(t(x)%*%x)%*%(beta2-beta0)) - 
              (t(beta2-beta1)%*%(t(x)%*%x)%*%(beta2-beta1)) ) / tau^2) 
  
  out <- TsWald
  
  out
}
