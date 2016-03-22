# library(MASS)
# library(restriktor)
# n=100
# x1 <- rnorm(n, 0)
# x2 <- rnorm(n, 0)
# x3 <- rnorm(n, 0)
# x4 <- rnorm(n, 0)
# y <- 1 + 1*x1 + 0.1*x2 + 0.3*x3 + 0.5*x4 + rnorm(n)
# fit.rlm <- rlm(y~x1+x2+x3+x4, method="MM")
# fit0.con <- restriktor(fit.rlm, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "const")
# fitA.con <- restriktor(fit.rlm, constraints = "x4 == 0", se = "const")
# 
# #conTest(fit0.con, test="wald", "A")
# x <- model.matrix(fit.rlm)[,,drop=FALSE]
# beta0 <- coef(fit0.con)
# betaA <- coef(fitA.con) #coef(fit.rlm)
# scale <- fit.rlm$s
# Amat <- fitA.con$Amat
# bvec <- fitA.con$bvec
# meq  <- fitA.con$meq 


robustWaldScores <- function(x, y, beta0, betaA, scale, 
                             Amat = NULL, bvec = NULL, meq = NULL) { 
  
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  #tukey bisquare tuning constant
  cc = 4.685061
    
  #Calculate M, Q, V 
  res0 <- y - x%*%beta0
  resA <- y - x%*%betaA
  
  rstar0 <- res0/scale
  rstarA <- resA/scale
  
  psi0   <- tukeyChi(rstar0, cc, deriv=1)  
  psiA   <- tukeyChi(rstarA, cc, deriv=1) 
  psideriv0 <- tukeyChi(rstar0, cc, deriv=2) 
  psiderivA <- tukeyChi(rstarA, cc, deriv=2) 
  
  #compute M 
  weightsM <- psiderivA / scale
#  weightsM <- psideriv0 / scale             
  WM <- weightsM %*% rep(1, p)
  xwM <- x * WM
  M <- t(x) %*% xwM / n
  
  #compute Q 
  weightsQ <- psiA^2
#  weightsQ <- psi0^2
  WQ <- weightsQ %*% rep(1, p)
  xwQ <- x * WQ
  Q <- t(x) %*% xwQ / n

  
  ############ <FIXME> ################
#   if (type == "A" | type == "global") {
#     idx0 <- which(colSums(abs(Amat)) == 0L)
#     idx1 <- which(colSums(abs(Amat)) > 0L)
#   } #else if (type == "B") {
     #idx0 <- which(colSums(abs(Amat)) == 0L)
     #idx0 <- idx1 <- which(colSums(abs(Amat)) > 0L)
   #} 
  
  #Calculate V 
  Minv <- solve(M)
  #information matrix 
  V <- Minv %*% Q %*% t(Minv)
  #V[(p0+1):pa,(p0+1):pa]
  #V[idx1,idx1]
  V22 <- Amat %*% V %*% t(Amat)
  
#   #Wald-type test statistic, Silvapulle (1996, eq. 2.6) & Silvapulle and Sen, 2005, p 154
#     Tn <- sqrt(n)*(Amat%*%betaA)
# #    #Tn <- sqrt(n)*betaA[idx1]
#     Dn <- Tn
# # #  Dmat <- solve(Amat[,idx1,drop=FALSE]%*%V22%*%t(Amat[,idx1,drop=FALSE]))
#     Dmat <- solve(Amat%*%V%*%t(Amat))
# # #   Dmat <- solve(V22)
#     dvec <- t(Dn)%*%Dmat
#     if (type == "A" | type == "global" | (type == "B" && meq.alt != 0)) {
#       out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(diag(length(dvec))), bvec=bvec, meq=meq)
#       b <- out$solution
#     } else {
#       b <- Dn
#     }
     
#    out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(Amat[,idx1]), bvec=bvec, meq=meq) 
#    TsWald <- as.numeric((t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b)) )
   #TsWald <- as.numeric(n * as.numeric(Amat%*%(betaA-beta0)) %*% solve(V22, as.numeric(t(Amat%*%(betaA-beta0))) ))
    TsWald <- as.numeric(n * c(betaA-beta0) %*% solve(V, c(betaA-beta0)))
   
  #unconstrained
  #TsWald <- as.numeric(n * betaA[idx1] %*% solve(V22, betaA[idx1]) )
  #as.numeric(n * (betaA[idx1]-beta0[idx1]) %*% solve(V22, (betaA[idx1]-beta0[idx1])) )
  #as.numeric(n * (betaA-beta0) %*% solve(V, (betaA-beta0)) )
  #  TsWald <- as.numeric(n * (betaA-beta0) %*% solve(V, (betaA-beta0)) )
   
   #Schur complement?
   # M221 <- M[(p0+1):p,(p0+1):p] - M[(p0+1):p,1:p0] %*% solve(M[1:p0,1:p0,drop=FALSE], M[1:p0,(p0+1):p,drop=FALSE])
#   M221 <- M[idx1,idx1] - M[idx1,idx0] %*% solve(M[idx0,idx0,drop=FALSE], 
#                                                 M[idx0,idx1,drop=FALSE])
   
   
   #<FIXME>
#   weightsZ0 <- psi0
#   Z0 <- (t(x) %*% weightsZ0) / n  
   #Z0 <- Amat%*%Z0
   
#   weightsZA <- psiA
#   ZA <- (t(x) %*% weightsZA) / n  
   #ZA <- Amat%*%ZA
   
   #n * t(Z0-ZA) %*% solve(V22, (Z0-ZA)) #n*(t(Z0 - ZA) %*% solve(V22) %*% (Z0 - ZA)) 
#   as.numeric(n * as.numeric(Amat%*%(ZA-Z0)) %*% solve(V22, as.numeric(t(Amat%*%(ZA-Z0))) ))
#   
   #result.C <- M221 %*% V22 %*% t(M221)
   #TsScore <- as.numeric(n * t(Z) %*% solve(result.C, Z))

   weightsZ0 <- psi0
   Z0 <- (t(x) %*% weightsZ0) / n  
   
   weightsZA <- psiA
   ZA <- (t(x) %*% weightsZA) / n  
   
#   #Score-type test (Silvapulle, 1996, eq. 2.6)
   #An <- sqrt(n) * solve(M221) %*% Z[idx1]
#   An <- sqrt(n) * solve(M) %*% Z0
#   Dn <- An
   #Dmat <- solve(Amat%*%V%*%t(Amat))
#   Dmat <- solve(V)
#   dvec <- t(Dn)%*%Dmat
   #out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(diag(length(dvec))), bvec=bvec, meq=meq) 
#   if (type == "A" | type == "global" | (type == "B" && meq.alt != 0)) {
#    out <- quadprog:::solve.QP(Dmat=Dmat, dvec=dvec, Amat=t(Amat), bvec=bvec, meq=meq) 
#    b <- out$solution
#   } else {
#     b <- Dn
#   }  
#   TsScore <- as.numeric((t(Dn)%*%Dmat%*%Dn) - (t(Dn-b)%*%Dmat%*%(Dn-b))) 
   
   result.C <- M %*% V %*% t(M)
   TsScore <- as.numeric(n * t(ZA-Z0) %*% solve(result.C, (ZA-Z0)))
   
   #TsScore <- as.numeric((t(Dn-b)%*%Dmat%*%(Dn-b)))
   #unconstrained
    #result.C <- M221 %*% V22 %*% t(M221)
    #TsScore <- as.numeric(n * t(Z[idx1]) %*% solve(result.C, Z[idx1]))
  
  
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
