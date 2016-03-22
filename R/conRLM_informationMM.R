informationMM <- function(x, betaA, scale) {

  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  #tukey bisquare tuning constant
  cc = 4.685061
  
  #Calculate M, Q, V 
#  res0 <- y - x%*%beta0
  resA <- y - x%*%betaA
  
#  rstar0 <- res0/scale
  rstar2 <- resA/scale
  
#  psi0   <- tukeyChi(rstar0, cc, deriv=1)  
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
  
  #Calculate V 
  Minv <- solve(M)
  #information matrix 
  V <- Minv %*% Q %*% t(Minv)
  
  V
}