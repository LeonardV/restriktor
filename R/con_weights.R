# compute weights based on multivariate normal distribution function with additional
# Monte Carlo steps.
# REF: Groeping, U. Inference with linear Equality and Inequality constraints
# using R: The Package ic.infer. 
con_weights <- function(cov, meq) {
  if (meq == 0L) {
    wt.bar <- try(ic.weights(cov))
  } else if (meq > 0) {
    wt.bar <- try(ic.weights(solve(solve(cov)[-c(1:meq), -c(1:meq)])))
  }
  if (inherits(wt.bar, "try-error")) {
    stop("Restriktor ERROR: the covariance matrix is too large. Try to set mix.weights = \"boot\".", call. = FALSE)
  }  
  wt.bar
}


## compute weights based on simulation.
## REF: Silvapulle and Sen (2005, p. 79). Constrained Statistical Inference: Order, 
## Inequality, and Shape Constraints. Hoboken, {NJ}: Wiley

## Add parallel option
con_weights_boot <- function(VCOV, Amat, meq, 
                             R = 99999L, seed = NULL, verbose = FALSE, ...) {
  
  ## seed
  if (!is.null(seed))
    set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  RNGstate <- .Random.seed
  
  ## weights do not depend on bvec.
  bvec <- rep(0L, nrow(Amat)) 
  invW <- solve(VCOV) 
  Dmat <- 2*invW
  Z <- rmvnorm(n = R, mean = rep(0, ncol(VCOV)), sigma = VCOV)
  dvec <- 2*(Z %*% invW)  
  ## create vector
  iact <- vector("numeric", R)
  for (i in 1:R) {
    if (all(Amat %*% Z[i, ] - bvec >= 0 * bvec) & meq == 0) {
      iact[i] <- 0L
    } else {
      QP <- try(solve.QP(Dmat = Dmat, 
                         dvec = dvec[i, ], 
                         Amat = t(Amat),
                         bvec = bvec, 
                         meq  = meq))
      if (verbose) {
        cat(" ...number of active constraints =", length(QP$iact), "\n")
      }
      if (inherits(QP, "try-error")) {
        if (verbose) {
          cat("quadprog FAILED\n")
        }
        iact[i] <- as.numeric(NA)
      }    
      iact[i] <- ifelse(QP$iact[1] == 0L, 0L, length(QP$iact)) 
    }
  }
  dimL   <- ncol(VCOV) - iact

  wt.bar <- sapply(0:ncol(VCOV), function(x) sum(x == dimL)) / R
  
  out <- wt.bar #dimWt
  out
}

