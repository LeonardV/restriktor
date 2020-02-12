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
con_weights_boot <- function(VCOV, Amat, meq, 
                             R = 99999L, parallel = c("no", "multicore", "snow"),
                             ncpus = 1L, cl = NULL, seed = NULL, 
                             verbose = FALSE, ...) {
  
  ## seed
  if (!is.null(seed))
    set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  RNGstate <- .Random.seed
  
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
  
  ## weights do not depend on bvec.
  bvec <- rep(0L, nrow(Amat)) 
  invW <- solve(VCOV) 
  Dmat <- 2*invW
  Z <- rmvnorm(n = R, mean = rep(0, ncol(VCOV)), sigma = VCOV)
  dvec <- 2*(Z %*% invW)
  res <- vector("numeric", R)
  
  
  fn <- function(b) {
    
    QP <- solve.QP(Dmat = Dmat, 
                   dvec = b, 
                   Amat = t(Amat),
                   bvec = bvec, 
                   meq  = meq)
    if (verbose) {
      cat(" ...active inequality constraints =", QP$iact, "\n")
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
      parallel::parApply(cl, dvec, 1, FUN = fn)
    }
    else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost",
                                             ncpus))
        parallel::clusterExport(cl, c("solve.QP"))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parRapply(cl, dvec, fn)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parRapply(cl, dvec, fn)
    }
  }
  else apply(dvec, 1, fn)
  
  dimL   <- ncol(VCOV) - res
  wt.bar <- sapply(0:ncol(VCOV), function(x) sum(x == dimL)) / R
  
  
  OUT <- wt.bar
  OUT
}


