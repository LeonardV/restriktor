## REF: Silvapull, M.J. (1996). Robust bounded influence tests against one-sided
## hypotheses in general parametric models. Statistics and Probability Letters, 31, 45 - 50.
robustWaldScores <- function(x, y, b.eqrestr, b.restr, b.unrestr, Amat,
                             scale, test = "wald", cc = 4.685061) { 
  
  test <- tolower(test)
  X <- as.matrix(x)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  ## Calculate M, Q, V 
  res2   <- y - X %*% b.unrestr
  rstar2 <- res2 / scale
  
  # rho function
  psi2      <- tukeyChi(rstar2, cc, deriv = 1) 
  psideriv2 <- tukeyChi(rstar2, cc, deriv = 2) 
  
  #compute M 
  weightsM <- psideriv2 / scale
  WM  <- weightsM %*% rep(1, p)
  xwM <- X * WM
  M <- t(X) %*% xwM / n
  
  #compute Q 
  weightsQ <- psi2^2
  WQ  <- weightsQ %*% rep(1, p)
  xwQ <- X * WQ
  Q <- t(X) %*% xwQ / n

  #Calculate V 
  Minv <- solve(M)
  #information matrix 
  V <- Minv %*% Q %*% t(Minv)

  idx1 <- which(colSums(abs(Amat)) > 0L)
  idx0 <- which(colSums(abs(Amat)) == 0L)
  
  V22 <- V[idx1, idx1] 
  V22.inv <- solve(V22)
  
  if (length(idx0) == 0) {
    M221 <- M
  } else {
    M221 <- M[idx1,idx1] - M[idx1,idx0] %*% solve(M[idx0,idx0,drop=FALSE], 
                                                  M[idx0,idx1,drop=FALSE])
  }
  
  result.C <- M221 %*% V22 %*% t(M221)
  
  # Wald test-statistic
  # REF: Silvapulle and Sen (2005, p 154)
  if (test == "wald") {
   Dn0 <- sqrt(n) * b.eqrestr
   Dn1 <- sqrt(n) * b.restr
   Dn2 <- sqrt(n) * b.unrestr
   
   #Ts <- (t(Dn2 - Dn0) %*% solve(V) %*% (Dn2 - Dn0)) - (t(Dn2 - Dn1) %*% solve(V) %*% (Dn2 - Dn1))
   Ts <- (t(Dn2 - Dn0)[idx1] %*% V22.inv %*% (Dn2 - Dn0)[idx1]) - 
     (t(Dn2 - Dn1)[idx1] %*% V22.inv %*% (Dn2 - Dn1)[idx1]) 
   
   test <- "Wald" 
  } else if (test == "score") {
    res0 <- y - X %*% b.eqrestr
    res1 <- y - X %*% b.restr
    res2 <- y - X %*% b.unrestr
    
    rstar0 <- res0 / scale
    rstar1 <- res1 / scale
    rstar2 <- res2 / scale
    
    # rho functions
    psi0 <- tukeyChi(rstar0, cc, deriv = 1)  
    psi1 <- tukeyChi(rstar1, cc, deriv = 1) 
    
    weightsZ0 <- psi0
    Z0 <- (t(X) %*% weightsZ0) / n  
    weightsZ1 <- psi1
    Z1 <- (t(X) %*% weightsZ1) / n  
    weightsZ2 <- psi2
    Z2 <- (t(X) %*% weightsZ2) / n  
    
    Z0 <- sqrt(n) * Z0
    Z1 <- sqrt(n) * Z1
    Z2 <- sqrt(n) * Z2
    
    Ts <- (t(Z2 - Z0)[idx1] %*% solve(result.C) %*% (Z2 - Z0)[idx1]) - 
        (t(Z2 - Z1)[idx1] %*% solve(result.C) %*% (Z2 - Z1)[idx1])
    
    #Ts <- as.numeric(n * t(Z0 - ZA) %*% solve(result.C) %*% c(Z0 - ZA))
    
    # Ts <- (t(Z2 - Z0) %*% solve(result.C) %*% (Z2 - Z0)) - 
    #   (t(Z2 - ZA) %*% solve(result.C) %*% (Z2 - ZA))
    
    test <- "Score"
  } 

  OUT <- list(test = test,
              Ts   = as.numeric(Ts),
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
              Ts   = as.numeric(Ts))
  
  OUT
}


## robust Wald statistic, Silvapulle (1992) ##
robustWaldXX <- function(x, b.eqrestr, b.restr, b.unrestr, Amat, tau) {
  X <- x
  W <- crossprod(X)
  idx1 <- which(colSums(abs(Amat)) > 0L)
  
  Ts <- as.numeric( 1/tau^2 * ( (t(b.unrestr - b.eqrestr)[idx1] %*% W[idx1, idx1] %*% (b.unrestr - b.eqrestr)[idx1]) -
                                (t(b.unrestr - b.restr)[idx1] %*% W[idx1, idx1] %*% (b.unrestr - b.restr)[idx1]) ) )

  OUT <- list(test = "Wald",
              Ts   = Ts)
  
  OUT
}
