## REF: Silvapull, M.J. (1996). Robust bounded influence tests against one-sided
## hypotheses in general parametric models. Statistics and Probability Letters, 31, 45 - 50.
robustScores <- function(x, y, b.eqrestr, b.restr, b.unrestr, Amat,
                         scale, test = "score", cc = 4.685061) { 
  
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
  #covariance matrix 
  V <- Minv %*% Q %*% t(Minv)

  idx1 <- which(colSums(abs(Amat)) > 0L)
  idx0 <- which(colSums(abs(Amat)) == 0L)
  
  V22 <- V[idx1, idx1] 
  #V22.inv <- solve(V22)
  
  if (length(idx0) == 0) {
    M221 <- M
  } else {
    M221 <- M[idx1,idx1] - M[idx1,idx0] %*% solve(M[idx0,idx0,drop=FALSE], 
                                                  M[idx0,idx1,drop=FALSE])
  }
  
  result.C <- M221 %*% V22 %*% t(M221)
  
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
  
  d0 <- sqrt(n) * Z0[idx1]
  d1 <- sqrt(n) * Z1[idx1]
  d2 <- sqrt(n) * Z2[idx1]
  
  Ts <- (t(d2 - d0) %*% solve(result.C, c(d2 - d0))) - 
      (t(d2 - d1) %*% solve(result.C, c(d2 - d1)))
    
  OUT <- list(test = "Score",
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
  p <- length(b.unrestr)
  
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
robustWald <- function(x, y, b.eqrestr, b.restr, b.unrestr, 
                       Amat, scale, cc = 4.685061) {
  X <- x
  XX <- crossprod(X)
  
  res   <- y - X %*% b.unrestr
  rstar <- res / scale
  
  # Yohai (1987, p. 648). Eq. 4.2, 4.3 and 4.4.
  # (for homoscedastic regression)
  A <- mean(tukeyChi(rstar, cc, deriv = 1)^2)
  B <- mean(tukeyChi(rstar, cc, deriv = 2))
  tau2 <- scale^2 * A/B^2
  
  U <- 1/tau2 * XX
  Ts <- as.numeric((t(b.unrestr - b.eqrestr) %*% U %*% c(b.unrestr - b.eqrestr)) - 
           (t(b.unrestr - b.restr) %*% U %*% c(b.unrestr - b.restr)))
  
  
  OUT <- list(test   = "Wald",
              Ts     = Ts,
              stddev = tau2)
  
  OUT
}
