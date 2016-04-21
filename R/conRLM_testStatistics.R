# Data <- restriktor::ZelazoKolb1972
# idx <- which(Data$Group == 3)
# Data <- Data[-idx, ]
# Data$Group <- factor(Data$Group)
# 
# fit.rlm <- rlm(Age ~ 1+Group, data=Data, method="MM")
# fit0.con <- restriktor(fit.rlm, constraints="Group2 == 0; Group2 == Group4")
# fitA.con <- restriktor(fit.rlm, constraints="Group2 > 0; Group2 < Group4")
# 
# x <- model.matrix(fit.rlm)[,,drop=FALSE]
# y <- as.matrix(fit.rlm$model[, attr(fit.rlm$terms, "response")])
# beta0 <- coef(fit0.con)
# betaA <- coef(fitA.con)
# scale <- fit.rlm$s
# Amat <- fit0.con$Amat

robustWaldScores <- function(x, y, beta0, betaA, scale, test = "wald") { 
  
  test <- tolower(test)
  X <- as.matrix(x)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #tukey bisquare tuning constant
  cc = 4.685061
    
  #Calculate M, Q, V 
  res0 <- y - X %*% beta0
  resA <- y - X %*% betaA
  
  rstar0 <- res0 / scale
  rstarA <- resA / scale
  
  psi0   <- tukeyChi(rstar0, cc, deriv = 1)  
  psiA   <- tukeyChi(rstarA, cc, deriv = 1) 
#  psideriv0 <- tukeyChi(rstar0, cc, deriv=2)  
  psiderivA <- tukeyChi(rstarA, cc, deriv = 2) 
  
  #compute M 
  weightsM <- psiderivA / scale
#  weightsM <- psideriv0 / scale             
  WM <- weightsM %*% rep(1, p)
  xwM <- X * WM
  M <- t(X) %*% xwM / n
  
  #compute Q 
  weightsQ <- psiA^2
#  weightsQ <- psi0^2
  WQ <- weightsQ %*% rep(1, p)
  xwQ <- X * WQ
  Q <- t(X) %*% xwQ / n

  #Calculate V 
  Minv <- solve(M)
  #information matrix 
  V <- Minv %*% Q %*% t(Minv)

  # idx1 <- which(colSums(Amat) != 0L)
  # idx0 <- which(colSums(Amat) == 0L)
  # result.V22 <- V[idx1,idx1]
  # result.W <- n * betaA[idx1] %*% solve(result.V22, betaA[idx1])
  
  # Wald test-statistic
  if (test == "wald") {
    Ts <- as.numeric(n * c(betaA-beta0) %*% solve(V, c(betaA-beta0)))
  } else if (test == "score") {
   # Score test-statistic
    weightsZ0 <- psi0
    Z0 <- (t(X) %*% weightsZ0) / n  
  
    weightsZA <- psiA
    ZA <- (t(X) %*% weightsZA) / n  
      
    result.C <- M %*% V %*% t(M)
    Ts <- as.numeric(n * t(ZA-Z0) %*% solve(result.C, (ZA-Z0)))
  } 

  OUT <- list(Ts = Ts,
              V = V)
  
  OUT
  
}



## robust Fm test statistic ##
robustFm <- function(x, y, beta0, betaA, scale, cc = 4.685061) {
  
  X <- x
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #compute residuals under null and alternative model
  resid0 <- y - X %*% beta0
  resid1 <- y - X %*% betaA
  
  #residuals / scale
  rstar0 <- as.numeric(resid0 / scale)                                               
  rstar1 <- as.numeric(resid1 / scale)
  
  L0 <- sum(tukeyChi(rstar0, cc, deriv = 0))
  L1 <- sum(tukeyChi(rstar1, cc, deriv = 0))
  
  #first derivative psi function
  psi.prime.h1 <- tukeyChi(rstar1, cc, deriv = 1)                     
  #second derivative psi function
  psi.prime2.h1 <- tukeyChi(rstar1, cc, deriv = 2)                    
  
  #asymptotic covariance matrix standardizing constant
  l.h1 <- ( 0.5 * (1 / (n - p)) * sum(psi.prime.h1^2) ) / ( (1/n) * sum(psi.prime2.h1) )  
  out <- 1 / l.h1 * (L0 - L1) 
    
  OUT <- out
  
  OUT
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
