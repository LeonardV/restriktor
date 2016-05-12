# compute weights based on multivariate normal distribution function.
# REF: Groeping, U. Inference with linear Equality and Inequality constraints
# using R: The Package ic.infer. 
con_wt <- function(cov, meq) {
  if (meq == 0L) {
    wt.bar <- ic.weights(cov)
  } else if (meq > 0) {
    wt.bar <- ic.weights(solve(solve(cov)[-(1:meq), -(1:meq)]))
  }
  wt.bar
}


# compute weights based on simulation.
#REF: Silvapulle and Sen (2005, p. 79). Constrained Statistical Inference: Order, 
# Inequality, and Shape Constraints. Hoboken, {NJ}: Wiley
mix.boot <- function(VCOV, Amat, meq, 
                     R = 9999, parallel = c("no", "multicore", "snow"),
                     ncpus = 1L, cl = NULL, seed = NULL, verbose = FALSE, ...) {
  
  
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
  
  bvec <- rep(0L, nrow(Amat))
  W <- VCOV 
  invW <- solve(W) 
  Dmat <- 2*invW
  
  iact <- vector("numeric", ncol(Amat))
  fn <- function(b) {
    if (verbose) {
      cat("iteration =", b)
    }
    if (!is.null(seed))
      set.seed(seed + b)
    if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1)
    RNGstate <- .Random.seed
    
    Z <- rmvnorm(n = 1, mean = rep(0, ncol(W)), sigma = W)
    dvec <- 2*(Z %*% invW)
    QP <- try(solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat),
                       bvec = bvec, meq = meq))
    if (verbose) {
      cat(" ...active inequality constraints =", QP$iact, "\n")
    }
    
    if (inherits(QP, "try-error")) {
      if (verbose) cat("quadprog FAILED\n")
      return(NULL)
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
      }
      else parallel::parLapply(cl, seq_len(RR), fn)
    }
  }
  else lapply(seq_len(RR), fn)
  error.idx <- integer(0)
  for (b in seq_len(R)) {
    if (!is.null(res[[b]])) {
      iact[b] <- res[[b]]
    }
    else {
      error.idx <- c(error.idx, b)
    }
  }
  # compute the number of positive components of W.
  dimL <- ncol(W) - iact
  wt <- sapply(1:(ncol(W) + 1), function(x) sum(x == (dimL + 1))) / R
  
  OUT <- wt
  
  OUT
}

