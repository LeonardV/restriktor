# compute weights based on multivariate normal distribution function with additional
# Monte Carlo steps.
# REF: Groemping, U. Inference with linear Equality and Inequality constraints
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
con_weights_boot <- function(VCOV, Amat, meq, 
                             R = 99999L, parallel = c("no", "multicore", "snow"),
                             ncpus = 1L, cl = NULL, seed = NULL, 
                             verbose = FALSE, ...) {
  
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore")
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow")
      have_snow <- TRUE
    if (!have_mc && !have_snow)
      ncpus <- 1L  
    }
  
  if (!is.null(seed))
    set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)
  RNGstate <- .Random.seed
  
  
  bvec <- rep(0L, nrow(Amat)) # weights do not depend on bvec.
  invW <- solve(VCOV) 
  Dmat <- 2*invW
  Z <- rmvnorm(n = R, mean = rep(0, ncol(VCOV)), sigma = VCOV)
  dvec <- 2*(Z %*% invW)
  
  
  fn <- function(b) {
    QP <- try(solve.QP(Dmat = Dmat, 
                       dvec = dvec[b, ], 
                       Amat = t(Amat),
                       bvec = bvec, 
                       meq  = meq), silent = TRUE)
    
    if (inherits(QP, "try-error")) {
      if (verbose) cat("quadprog FAILED\n")
      return(NULL)
    } else {
      if (verbose) {
        cat(" ...active inequality constraints =", QP$iact, "\n")
      }
      
    }    
    
    if (QP$iact[1] == 0L) {
      return(0L) 
    } else { 
      return(length(QP$iact))
    }  
  }
  
  RR <- sum(R)
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
    }
    else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost",
                                             ncpus))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        
        res <- parallel::parLapply(cl, seq_len(RR), fn)
        parallel::stopCluster(cl)
        res
      } else { 
        parallel::parLapply(cl, seq_len(RR), fn)
      }
    }
  }
  else lapply(seq_len(RR), fn) 
  
  error.idx <- integer(0)
  iact <- vector("numeric", ncol(Amat))
  
  for (b in seq_len(R)) {
    if (!is.null(res[[b]])) {
      iact[b] <- res[[b]]
    }
    else {
      error.idx <- c(error.idx, b)
    }
  }
  
  # compute the number of positive components of VCOV.
  # ncol(VCOV) = maximum number of constraints
  # iact       = number of active inequality constraints
  if (length(error.idx) > 0) {
    dimL <- ncol(VCOV) - iact[-error.idx]
  } else {
    dimL <- ncol(VCOV) - iact
  }
  wt.bar <- sapply(0:ncol(VCOV), function(x) sum(x == dimL)) / length(dimL)
  attr(wt.bar, "error.idx") <- error.idx
  
  wt.bar
}


