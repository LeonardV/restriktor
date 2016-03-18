# library(MASS)
# library(restriktor)
# n=100
# x1 <- rnorm(n, 0)
# x2 <- rnorm(n, 0)
# x3 <- rnorm(n, 0)
# x4 <- rnorm(n, 0)
# y <- 1 + 1*x1 + 0.1*x2 + 0.3*x3 + 0.5*x4 + rnorm(n)
# fit.rlm <- rlm(y~x1+x2+x3+x4, method="MM")
# fit0.con <- restriktor(fit.rlm, constraints = "x2 == 0; x3 == 0; x4 == 0", se = "const")
# fitA.con <- restriktor(fit.rlm, constraints = "x2 > 0; x3 > 0; x4 < 0", se = "const")
# 
# #conTest(fit0.con, test="wald", "A")
# 
# 
# 
# x <- model.matrix(fit.rlm)[,,drop=FALSE]
# beta0 <- coef(fit0.con)
# betaA <- coef(fitA.con) #coef(fit.rlm)
# scale <- fit.rlm$s
# Amat <- fitA.con$Amat
# bvec <- fitA.con$bvec
# meq  <- fitA.con$meq 


robustWaldScores <- function(x, y, type = "A", beta0, betaA, scale, 
                             Amat = NULL, bvec = NULL, meq = NULL, meq.alt = 0) { 
  
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
#  psideriv0 <- tukeyChi(rstar0, cc, deriv=2) 
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

  
  ############ <FIXME> ################
#   if (type == "A" | type == "global") {
#     idx0 <- which(colSums(abs(Amat[,1:p])) == 0L)
#     idx1 <- which(colSums(abs(Amat[,1:p])) > 0L)
#   } else if (type == "B") {
#     idx0 <- 1:p#which(colSums(abs(Amat[,1:p])) == 0L)
#     idx1 <- 1:p#which(colSums(abs(Amat[,1:p])) > 0L)
#   } 
  
  #Calculate V 
  Minv <- solve(M)
  #information matrix 
  V <- Minv %*% Q %*% t(Minv)
  
  #V[(p0+1):pa,(p0+1):pa]
#  V22 <- V[idx1,idx1]
  
  TsWald <- as.numeric(n * (betaA-beta0) %*% solve(V, (betaA-beta0)) )
  #inverse information matrix
  #V.inv <- solve(V)
  #submatrix of V
#  V22.inv <- solve(V22)
  
  #Schur complement?
  # M221 <- M[(p0+1):p,(p0+1):p] - M[(p0+1):p,1:p0] %*% solve(M[1:p0,1:p0,drop=FALSE], M[1:p0,(p0+1):p,drop=FALSE])
  #M221 <- M[idx1,idx1] - M[idx1,idx0] %*% solve(M[idx0,idx0,drop=FALSE], 
  #                                              M[idx0,idx1,drop=FALSE])

  #<FIXME>
  weightsZ <- psi0
  Z <- t(x) %*% weightsZ / n  

  result.C <- M %*% V %*% t(M)
  TsScore  <- as.numeric(n * t(Z) %*% solve(result.C, Z))
  
  
#   
#   #Wald-type test statistic, Silvapulle (1996, eq. 2.6)
#   #Tn <- sqrt(n)*(Amat%*%betaA)
#   Tn <- sqrt(n)*betaA[idx1]
#   Dn <- Tn
#   #Dmat <- solve(Amat[,idx1,drop=FALSE]%*%V22%*%t(Amat[,idx1,drop=FALSE]))
#   Dmat <- solve(V22)
#   dvec <- t(Dn)%*%Dmat
#   #out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(diag(length(dvec))), bvec=bvec, meq=meq) 
#   out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(Amat[,idx1]), bvec=bvec, meq=meq) 
#   b <- out$solution
#   TsWald  <- as.numeric((t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b)) )
#   
  #  TsWald  <- t(Dn-b)%*%Dmat%*%(Dn-b) 
  #unconstrained
  #TsWald <- as.numeric(n * betaA[idx1] %*% solve(V22, betaA[idx1]) )
  #as.numeric(n * (betaA[idx1]-beta0[idx1]) %*% solve(V22, (betaA[idx1]-beta0[idx1])) )
  #as.numeric(n * (betaA-beta0) %*% solve(V, (betaA-beta0)) )
  
  #Score-type test (Silvapulle, 1996, eq. 2.6)
  #An <- sqrt(n) * Amat[,idx1,drop=FALSE]%*%(solve(M221) %*% Z[idx1])
#   An <- sqrt(n) * solve(M221) %*% Z[idx1]
#   Dn <- An
#   #Dmat <- solve(Amat[,idx1,drop=FALSE]%*%V22%*%t(Amat[,idx1,drop=FALSE]))
#   Dmat <- solve(V22)
#   dvec <- t(Dn)%*%Dmat
#   out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(diag(length(dvec))), bvec=bvec, meq=meq) 
#   b <- out$solution
#   TsScore <- as.numeric((t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b)))
#   #TsScore <- (t(Dn-b)%*%Dmat%*%(Dn-b))
#   #unconstrained
#   #  result.C <- M221 %*% V22 %*% t(M221)
#   #  TsScore <- as.numeric(n * t(Z[idx1]) %*% solve(result.C, Z[idx1]))
  
  
  OUT <- list(RWald = TsWald,
              Rscore = TsScore,
              V = V)#,
              #V22 = Dmat)
  
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
