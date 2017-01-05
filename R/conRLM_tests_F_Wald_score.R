robustWaldScores <- function(x, y, b_eqrestr, b_restr, b_unrestr, 
                             scale, test = "wald", cc = 4.685061) { 
  
  test <- tolower(test)
  X <- as.matrix(x)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Calculate M, Q, V 
  resA <- y - X %*% b_unrestr
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

  # Wald test-statistic
  if (test == "wald") {
    Ts <- as.numeric(n * c(b_restr - b_eqrestr) %*% 
                       solve(V, c(b_restr - b_eqrestr)))
    
    test <- "Wald"
  } else if (test == "score") {
    # Score test-statistic
    res0 <- y - X %*% b_eqrestr
    resA <- y - X %*% b_restr
    
    rstar0 <- res0 / scale
    rstarA <- resA / scale
    
    # rho functions
    psi0 <- tukeyChi(rstar0, cc, deriv = 1)  
    psiA <- tukeyChi(rstarA, cc, deriv = 1) 
    
    weightsZ0 <- psi0
    Z0 <- (t(X) %*% weightsZ0) / n  
  
    weightsZA <- psiA
    ZA <- (t(X) %*% weightsZA) / n  
      
    result_C <- M %*% V %*% t(M)
    Ts <- as.numeric(n * t(ZA - Z0) %*% solve(result_C, (ZA - Z0)))
    
    test <- "Score"
  } 

  OUT <- list(test = test,
              Ts = Ts,
              V  = V)
  
  OUT
  
}



################# robust Fm test statistic ######################
# REF: Mervyn J. Silvapulle (1992). Robust Tests of Inequality 
# Constraints and One-Sided Hypotheses in the Linear Model, 
# Biometrika, 79, 3, 621-630.
robustFm <- function(x, y, b_unrestr, b_eqrestr, b_restr, scale, 
                     cc = 4.685061) {
  X <- x
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #compute residuals under null, restrikted and unconstrained model
  resid0 <- y - X %*% b_eqrestr
  resid1 <- y - X %*% b_restr
  resid2 <- y - X %*% b_unrestr
  
  #residuals / scale
  rstar0 <- as.numeric(resid0 / scale)                                               
  rstar1 <- as.numeric(resid1 / scale)
  rstar2 <- as.numeric(resid2 / scale)
  
  L0 <- sum(tukeyChi(rstar0, cc, deriv = 0))
  L1 <- sum(tukeyChi(rstar1, cc, deriv = 0))
  
  #first derivative rho function
  psi_prime_h2 <- tukeyChi(rstar2, cc, deriv = 1) 
  #second derivative rho function
  psi_prime2_h2 <- tukeyChi(rstar2, cc, deriv = 2) 
  
  #asymptotic covariance matrix standardizing constant
  lambda <- ( 0.5 * (1 / (n - p)) * sum(psi_prime_h2^2) ) / 
                                  ( (1/n) * sum(psi_prime2_h2) )  
  Ts <- 1 / lambda * (L0 - L1) 
  
  OUT <- list(test = "F",
              Ts   = Ts)
  
  OUT
}


## robust Wald statistic, Silvapulle (1992) ##
robustWaldXX <- function(x, b_eqrestr, b_restr, b_unrestr, tau) {
  #x <- scale(x, center = TRUE, scale = FALSE)
  Ts <- as.numeric(( (t(b_unrestr-b_eqrestr)%*%(t(x)%*%x)%*%(b_unrestr-b_eqrestr)) -
                     (t(b_unrestr-b_restr)%*%(t(x)%*%x)%*%(b_unrestr-b_restr)) ) / tau^2)

  OUT <- list(test = "Wald",
              Ts   = Ts)
  
  OUT
}
