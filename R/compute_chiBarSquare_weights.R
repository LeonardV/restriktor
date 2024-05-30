# compute weights based on multivariate normal distribution function with additional
# Monte Carlo steps.
# REF: Groemping, U. Inference with linear Equality and Inequality constraints
# using R: The Package ic.infer. 
con_weights <- function(cov, meq) {
  if (meq == 0L) {
    wt_bar <- try(ic.weights(cov))
  } else if (meq > 0) {
    wt_bar <- try(ic.weights(solve(solve(cov)[-c(1:meq), -c(1:meq)])))
  }
  if (inherits(wt_bar, "try-error")) {
    stop("Restriktor ERROR: the covariance matrix is too large. Try to set mix_weights = \"boot\".", call. = FALSE)
  }  
  wt_bar
}


## compute weights based on simulation.
## REF: Silvapulle and Sen (2005, p. 79). Constrained Statistical Inference: Order, 
## Inequality, and Shape Constraints. Hoboken, {NJ}: Wiley

# con_weights_boot <- function(VCOV, Amat, meq,
#                              R = 99999L, parallel = c("no", "multicore", "snow"),
#                              ncpus = 1L, cl = NULL, seed = NULL,
#                              verbose = FALSE, ...) {
# 
#   parallel <- match.arg(parallel)
#   have_mc <- have_snow <- FALSE
#   if (parallel != "no" && ncpus > 1L) {
#     if (parallel == "multicore")
#       have_mc <- .Platform$OS.type != "windows"
#     else if (parallel == "snow")
#       have_snow <- TRUE
#     if (!have_mc && !have_snow)
#       ncpus <- 1L
#     }
# 
#   if (!is.null(seed)) set.seed(seed)
#   if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
# 
#   bvec <- rep(0L, nrow(Amat)) # weights do not depend on bvec.
#   invW <- solve(VCOV)
#   Dmat <- 2*invW
#   Z <- rmvnorm(n = R, mean = rep(0, ncol(VCOV)), sigma = VCOV)
#   dvec <- 2*(Z %*% invW)
# 
# 
#   fn <- function(b) {
#     QP <- try(solve.QP(Dmat = Dmat,
#                        dvec = dvec[b, ],
#                        Amat = t(Amat),
#                        bvec = bvec,
#                        meq  = meq), silent = TRUE)
# 
#     if (inherits(QP, "try-error")) {
#       if (verbose) cat("quadprog FAILED\n")
#       return(NULL)
#     } else {
#       if (verbose) {
#         cat(" ...active inequality constraints =", QP$iact, "\n")
#       }
#     }
# 
#     if (QP$iact[1] == 0L) {
#       return(0L)
#     } else {
#       return(length(QP$iact))
#     }
#   }
# 
#   RR <- sum(R)
#   res <- if (ncpus > 1L && (have_mc || have_snow)) {
#     if (have_mc) {
#       parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
#     } else if (have_snow) {
#        if (is.null(cl)) {
#         cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
#         if (RNGkind()[1L] == "L'Ecuyer-CMRG")
#           parallel::clusterSetRNGStream(cl)
#           res <- parallel::parLapply(cl, seq_len(RR), fn)
#           parallel::stopCluster(cl)
#           res
#        } else {
#           parallel::parLapply(cl, seq_len(RR), fn)
#        }
#       }
#   } else lapply(seq_len(RR), fn)
# 
#   # error.idx <- integer(0)
#   # iact <- vector("numeric", ncol(Amat))
#   #
#   # for (b in seq_len(R)) {
#   #   if (!is.null(res[[b]])) {
#   #     iact[b] <- res[[b]]
#   #   }
#   #   else {
#   #     error.idx <- c(error.idx, b)
#   #   }
#   # }
# 
#   # compute the number of positive components of VCOV.
#   # ncol(VCOV) = maximum number of constraints
#   # iact       = number of active inequality constraints
#   iact <- sapply(res, function(x) ifelse(is.null(x), NA, x))
#   error.idx <- which(is.na(iact))
# 
#   dimL <- ncol(VCOV) - iact
#   if (length(error.idx) > 0) {
#     dimL <- ncol(VCOV) - iact[-error.idx]
#   } else {
#     dimL <- ncol(VCOV) - iact
#   }
# 
#   wt_bar <- sapply(0:ncol(VCOV), function(x) sum(x == dimL)) / length(dimL)
#   attr(wt_bar, "error.idx") <- error.idx
# 
#   wt_bar
# }


## to avoid unnecessary runs we add samples of 'chunk_size' until 
## convergence is reached or the maximum number of bootstrap draws (R) is reached. 
## By default the rtmvnorm function is equal to the rmvnorm function (lower = -Inf, upper = Inf)
con_weights_boot <- function(VCOV, Amat, meq, R = 1e5L, 
                             chunk_size = 5000L, convergence_crit = 1e-03, 
                             seed = NULL, verbose = FALSE, ...) {
  
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  
  fn <- function(b, verbose) {
    QP <- try(quadprog::solve.QP(Dmat = Dmat,
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
  
  ldots <- list(...)
  prev_wt_bar    <- NULL
  chunk_wt_bar   <- NULL
  prev_res       <- NULL
  prev_error_idx <- NULL
  has_converged  <- FALSE
  chunk_size_org <- chunk_size
  
  bvec <- rep(0L, nrow(Amat)) # weights do not depend on bvec.
  invW <- solve(VCOV)
  Dmat <- 2*invW
  #Z <- mvtnorm::rmvnorm(n = R, mean = rep(0, ncol(VCOV)), sigma = VCOV)
  # truncated, the default (lower = -Inf and upper = Inf) is not truncated.
  Z <- tmvtnorm::rtmvnorm(n = R, mean = rep(0, ncol(VCOV)), sigma = VCOV, ...)
  dvec <- 2 * (Z %*% invW)
  RR <- sum(R)
  
  wt_bar <- numeric(ncol(VCOV) + 1)
  total_chunks <- ceiling(RR / chunk_size_org)
  chunk_iter <- 1
  
  while (chunk_iter <= total_chunks) { 
    start <- 1 + chunk_size_org * (chunk_iter - 1)
    end <- min(chunk_size_org * chunk_iter, RR)
    res <- lapply(start:end, fn, verbose = verbose)
    
    # Calculate wt_bar
    prev_res <- append(prev_res, res)
    iact <- sapply(prev_res, function(x) ifelse(is.null(x), NA, x))
    error_idx <- which(is.na(iact))
    if (length(error_idx) > 0) {
      dimL <- ncol(VCOV) - iact[-error_idx]
    } else {
      dimL <- ncol(VCOV) - iact
    }
    wt_bar <- sapply(0:ncol(VCOV), function(x) sum(x == dimL)) / length(dimL)
  
    # Check for convergence
    if (!is.null(prev_wt_bar)) { 
      if (all(abs(wt_bar - prev_wt_bar) < convergence_crit)) {
        has_converged  <- TRUE
        chunk_wt_bar   <- rbind(chunk_wt_bar, wt_bar)
        prev_error_idx <- c(error_idx, prev_error_idx)
        
        break
      }
    }
    
    # aanpassen argumenten namen in handleidingen
    prev_wt_bar    <- wt_bar
    chunk_wt_bar   <- rbind(chunk_wt_bar, wt_bar)
    prev_error_idx <- c(error_idx, prev_error_idx)
    chunk_iter     <- chunk_iter + 1L
    chunk_size     <- chunk_size + chunk_size_org
  }

  rownames(chunk_wt_bar) <- paste0("chunk_iter_", seq_len(nrow(chunk_wt_bar)))
  
  attr(wt_bar, "total_bootstrap_draws") <- length(iact)
  attr(wt_bar, "converged"            ) <- has_converged
  attr(wt_bar, "convergence_crit"     ) <- convergence_crit
  attr(wt_bar, "wt_bar_chunk"         ) <- chunk_wt_bar
  attr(wt_bar, "chunk_size"           ) <- chunk_size_org
  attr(wt_bar, "total_chunks"         ) <- total_chunks
  attr(wt_bar, "chunk_iter"           ) <- chunk_iter 
  attr(wt_bar, "error.idx"            ) <- prev_error_idx
  attr(wt_bar, "mvtnorm"              ) <- ldots
  
  return(wt_bar)
}


