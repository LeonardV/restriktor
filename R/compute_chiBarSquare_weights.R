# compute weights based on multivariate normal distribution function with additional
# Monte Carlo steps.
# REF: Groemping, U. Inference with linear Equality and Inequality constraints
# using R: The Package ic.infer. 
con_weights <- function(cov, meq, tolerance, ridge_constant, ...) {
  if (meq == 0L) {
    wt_bar <- try(ic_weights(cov, tolerance, ridge_constant, ...))
  } else if (meq > 0) {
    wt_bar <- try(ic_weights(solve(solve(cov)[-c(1:meq), -c(1:meq)]), 
                             tolerance, ridge_constant, ...))
  }
  if (inherits(wt_bar, "try-error")) {
    stop("restriktor ERROR: the covariance matrix is too large. Try to set mix_weights = \"boot\".", call. = FALSE)
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
  
  ldots <- list(...)
  p <- ncol(VCOV)
  
  bvec   <- rep(0L, nrow(Amat)) # weights do not depend on bvec.
  invW   <- solve(VCOV)
  Dmat   <- 2 * invW
  tAmat  <- t(Amat)             # pre-compute transpose once
  mu_vec <- rep(0, p)
  
  # inner function: returns number of active constraints, or NA on error
  fn <- function(b) {
    QP <- try(quadprog::solve.QP(Dmat = Dmat,
                                 dvec = dvec[b, ],
                                 Amat = tAmat,
                                 bvec = bvec,
                                 meq  = meq), silent = TRUE)
    
    if (inherits(QP, "try-error")) {
      if (verbose) cat("quadprog FAILED\n")
      return(NA_integer_)
    } else if (verbose) {
      cat(" ...active inequality constraints =", QP$iact, "\n")
    }
    
    if (QP$iact[1] == 0L) 0L else length(QP$iact)
  }
  
  prev_wt_bar    <- NULL
  chunk_wt_bar   <- NULL
  has_converged  <- FALSE
  chunk_size_org <- chunk_size
  
  # incremental level counts instead of storing all results
  level_counts   <- integer(p + 1)  # counts for levels 0:p
  n_valid        <- 0L
  error_idx_all  <- integer(0)
  total_generated <- 0L
  
  RR <- as.integer(R)
  total_chunks <- ceiling(RR / chunk_size_org)
  chunk_iter   <- 1L
  
  while (chunk_iter <= total_chunks) { 
    current_n <- min(chunk_size_org, RR - total_generated)
    
    # lazy sample generation: only generate what this chunk needs
    Z    <- tmvtnorm::rtmvnorm(n = current_n, mean = mu_vec, sigma = VCOV, ...)
    dvec <- 2 * (Z %*% invW)
    
    # process chunk with vapply (faster than lapply + sapply post-processing)
    chunk_iact <- vapply(seq_len(current_n), fn, integer(1))
    
    # update running counts incrementally (no need to recompute from scratch)
    valid_mask   <- !is.na(chunk_iact)
    n_chunk_err  <- sum(!valid_mask)
    
    if (n_chunk_err > 0L) {
      error_idx_all <- c(error_idx_all, which(!valid_mask) + total_generated)
    }
    
    valid_iact <- chunk_iact[valid_mask]
    dimL_chunk <- p - valid_iact  # values in 0:p
    
    # tabulate is O(n) vs sapply+sum which is O(n*p)
    if (length(dimL_chunk) > 0L) {
      level_counts <- level_counts + tabulate(dimL_chunk + 1L, nbins = p + 1)
    }
    n_valid         <- n_valid + sum(valid_mask)
    total_generated <- total_generated + current_n
    
    wt_bar <- level_counts / n_valid
    
    # Check for convergence
    if (!is.null(prev_wt_bar)) { 
      if (all(abs(wt_bar - prev_wt_bar) < convergence_crit)) {
        has_converged <- TRUE
        chunk_wt_bar  <- rbind(chunk_wt_bar, wt_bar)
        break
      }
    }
    
    prev_wt_bar  <- wt_bar
    chunk_wt_bar <- rbind(chunk_wt_bar, wt_bar)
    chunk_iter   <- chunk_iter + 1L
  }
  
  rownames(chunk_wt_bar) <- paste0("chunk_iter_", seq_len(nrow(chunk_wt_bar)))
  
  attr(wt_bar, "total_bootstrap_draws") <- total_generated
  attr(wt_bar, "converged"            ) <- has_converged
  attr(wt_bar, "convergence_crit"     ) <- convergence_crit
  attr(wt_bar, "wt_bar_chunk"         ) <- chunk_wt_bar
  attr(wt_bar, "chunk_size"           ) <- chunk_size_org
  attr(wt_bar, "total_chunks"         ) <- total_chunks
  attr(wt_bar, "chunk_iter"           ) <- chunk_iter 
  attr(wt_bar, "error.idx"            ) <- error_idx_all
  attr(wt_bar, "mvtnorm"              ) <- ldots
  
  return(wt_bar)
}
