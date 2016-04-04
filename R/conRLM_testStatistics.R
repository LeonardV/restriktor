robustWaldScores <- function(x, y, beta0, betaA, scale) { 
  
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
#  psideriv0 <- tukeyChi(rstar0, cc, deriv=2) 
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

  #Calculate V 
  Minv <- solve(M)
  #information matrix 
  V <- Minv %*% Q %*% t(Minv)

  # Wald test-statistic
  TsWald <- as.numeric(n * c(betaA-beta0) %*% solve(V, c(betaA-beta0)))
   
  # Score test-statistic
  weightsZ0 <- psi0
  Z0 <- (t(x) %*% weightsZ0) / n  

  weightsZA <- psiA
  ZA <- (t(x) %*% weightsZA) / n  
    
  result.C <- M %*% V %*% t(M)
  TsScore <- as.numeric(n * t(ZA-Z0) %*% solve(result.C, (ZA-Z0)))
   

  OUT <- list(RWald = TsWald,
              Rscore = TsScore,
              V = V)
  
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
  rstar0 <- as.numeric(resid0/scale)                                               
  rstar1 <- as.numeric(resid1/scale)
  
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
# robustWaldXX <- function(x, beta0, beta1, beta2, tau) {
#   TsWald <- c(( (t(beta2-beta0)%*%(t(x)%*%x)%*%(beta2-beta0)) - 
#                 (t(beta2-beta1)%*%(t(x)%*%x)%*%(beta2-beta1)) ) / tau^2) 
#   
#   out <- TsWald
#   
#   out
# }
