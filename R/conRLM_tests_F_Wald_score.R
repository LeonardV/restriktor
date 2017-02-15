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

  # Wald test-statistic
  # REF: Silvapulle and Sen (2005, p 154)
  if (test == "wald") {
     # Ts <- as.numeric(n * t(Amat%*%(b.restr - b.eqrestr)) %*% 
     #                  Amat%*%solve(V)%*%t(Amat) %*% Amat%*%(b.restr - b.eqrestr))
     # 
    Dn  <- sqrt(n) * b.unrestr
    D0n <- sqrt(n) * b.eqrestr
    D1n <- sqrt(n) * b.restr
    Ts <- (t(Dn - D0n) %*% solve(V) %*% (Dn - D0n)) - (t(Dn - D1n) %*% solve(V) %*% (Dn - D1n))
    
    # Ts <- (t(Amat%*%(Dn - D0n)) %*% (Amat%*%solve(V)%*%t(Amat)) %*% Amat%*%(Dn - D0n)) - 
    #   (t(Amat%*%(Dn - D1n)) %*% (Amat%*%solve(V)%*%t(Amat)) %*% Amat%*%(Dn - D1n))
    # 
    test <- "Wald"
  } else if (test == "score") {
    # Global score test-statistic
    # REF: Robertson et al. (1988)
    res0 <- y - X %*% b.eqrestr
    resA <- y - X %*% b.restr
    res2 <- y - X %*% b.unrestr
    
    rstar0 <- res0 / scale
    rstarA <- resA / scale
    rstar2 <- res2 / scale
    
    # rho functions
    psi0 <- tukeyChi(rstar0, cc, deriv = 1)  
    psiA <- tukeyChi(rstarA, cc, deriv = 1) 
    psi2 <- tukeyChi(rstar2, cc, deriv = 1) 
    
    weightsZ0 <- psi0
    Z0 <- (t(X) %*% weightsZ0) / n  
    weightsZA <- psiA
    ZA <- (t(X) %*% weightsZA) / n  
    weightsZ2 <- psi2
    Z2 <- (t(X) %*% weightsZ2) / n  
    
    result.C <- M %*% V %*% t(M)
    
    #Ts <- as.numeric(n * t(Z0 - ZA) %*% solve(result.C) %*% c(Z0 - ZA))
    Z0 <- sqrt(n) * Z0
    ZA <- sqrt(n) * ZA
    Z2 <- sqrt(n) * Z2
    Ts <- (t(Z2 - Z0) %*% solve(result.C) %*% (Z2 - Z0)) - 
      (t(Z2 - ZA) %*% solve(result.C) %*% (Z2 - ZA))
    
    # Ts <- (t(Amat%*%(Z2 - Z0)) %*% (Amat%*%solve(result.C)%*%t(Amat)) %*% Amat%*%(Z2 - Z0)) - 
    #   (t(Amat%*%(Z2 - ZA)) %*% (Amat%*%solve(result.C)%*%t(Amat)) %*% Amat%*%(Z2 - ZA))
    # 
    test <- "Score"
  } 

  OUT <- list(test = test,
              Ts   = Ts,
              V    = V)
  
  OUT
  
}

################# robust Fm test statistic ######################
# REF: Silvapulle, M.J. (1992). Robust Tests of Inequality 
# Constraints and One-Sided Hypotheses in the Linear Model, 
# Biometrika, 79, 3, 621-630.
robustFm <- function(x, y, b.unrestr, b.eqrestr, b.restr, scale, 
                     cc = 4.685061) {
  X <- x
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #compute residuals under null, restricted and unconstrained model
  resid0 <- y - X %*% b.eqrestr
  resid1 <- y - X %*% b.restr
  resid2 <- y - X %*% b.unrestr
  
  #residuals / scale
  rstar0 <- as.numeric(resid0 / scale)                                               
  rstar1 <- as.numeric(resid1 / scale)
  rstar2 <- as.numeric(resid2 / scale)
  
  L0 <- sum(tukeyChi(rstar0, cc, deriv = 0))
  L1 <- sum(tukeyChi(rstar1, cc, deriv = 0))
  
  #first derivative rho function
  psi.prime.h2 <- tukeyChi(rstar2, cc, deriv = 1) 
  #second derivative rho function
  psi.prime2.h2 <- tukeyChi(rstar2, cc, deriv = 2) 
  
  #asymptotic covariance matrix standardizing constant
  lambda <- ( 0.5 * (1 / (n - p)) * sum(psi.prime.h2^2) ) / 
                                  ( (1/n) * sum(psi.prime2.h2) )  
  Ts <- 1 / lambda * (L0 - L1) 
  
  OUT <- list(test = "F",
              Ts   = Ts)
  
  OUT
}


## robust Wald statistic, Silvapulle (1992) ##
robustWaldXX <- function(x, Amat, b.eqrestr, b.restr, b.unrestr, tau) {
  X <- x
  n <- dim(X)[1]
  Ts <- as.numeric( 1/tau^2 * ( (t(b.unrestr - b.eqrestr) %*% (t(X)%*%X) %*% (b.unrestr - b.eqrestr)) -
                     (t(b.unrestr - b.restr) %*% (t(X)%*%X) %*% (b.unrestr - b.restr)) ) )

  OUT <- list(test = "Wald",
              Ts   = Ts)
  
  OUT
}
