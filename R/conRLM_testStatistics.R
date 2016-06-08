robustWaldScores <- function(x, y, b.eqrestr, b.restr, b.unrestr, 
                             scale, test = "wald", cc = 4.685061) { 
  
  test <- tolower(test)
  X <- as.matrix(x)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Calculate M, Q, V 
  resA <- y - X %*% b.unrestr
  rstarA <- resA / scale
  # rho function
  psiA <- tukeyChi(rstarA, cc, deriv = 1) 
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
  # result.W <- n * b.restr[idx1] %*% solve(result.V22, b.restr[idx1])
  
  # Wald test-statistic
  if (test == "wald") {
    Ts <- as.numeric(n * c(b.restr-b.eqrestr) %*% solve(V, c(b.restr-b.eqrestr)))
  } else if (test == "score") {
    # Score test-statistic
    res0 <- y - X %*% b.eqrestr
    resA <- y - X %*% b.restr
    
    rstar0 <- res0 / scale
    rstarA <- resA / scale
    
    # rho functions
    psi0 <- tukeyChi(rstar0, cc, deriv = 1)  
    psiA <- tukeyChi(rstarA, cc, deriv = 1) 
    
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
robustFm <- function(x, y, b.eqrestr, b.restr, scale, cc = 4.685061) {
  
  X <- x
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #compute residuals under null and alternative model
  resid0 <- y - X %*% b.eqrestr
  residA <- y - X %*% b.restr
  
  #residuals / scale
  rstar0 <- as.numeric(resid0 / scale)                                               
  rstarA <- as.numeric(residA / scale)
  
  L0 <- sum(tukeyChi(rstar0, cc, deriv = 0))
  LA <- sum(tukeyChi(rstarA, cc, deriv = 0))
  
  #first derivative rho function
  psi.prime.hA <- tukeyChi(rstarA, cc, deriv = 1) 
  #second derivative rho function
  psi.prime2.hA <- tukeyChi(rstarA, cc, deriv = 2) 
  
  #asymptotic covariance matrix standardizing constant
  l.hA <- ( 0.5 * (1 / (n - p)) * sum(psi.prime.hA^2) ) / 
                                  ( (1/n) * sum(psi.prime2.hA) )  
  OUT <- 1 / l.hA * (L0 - LA) 
    
  OUT
}


## robust Wald statistic, Silvapulle (1992) ##
# robustWaldXX <- function(x, b.eqrestr, beta1, beta2, tau) {
#   TsWald <- c(( (t(beta2-b.eqrestr)%*%(t(x)%*%x)%*%(beta2-b.eqrestr)) - 
#                 (t(beta2-beta1)%*%(t(x)%*%x)%*%(beta2-beta1)) ) / tau^2) 
#   
#   out <- TsWald
#   
#   out
# }
