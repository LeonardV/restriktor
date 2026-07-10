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


## to avoid unnecessary runs we add samples of 'chunk_size' until the Monte
## Carlo standard error of the weights is small enough (see convergence_crit)
## or the maximum number of bootstrap draws (R) is reached.
##
## Each draw projects z ~ N(0, VCOV) onto the cone {theta: Amat theta >= 0}
## (first meq rows treated as equalities) and counts the active constraints.
## The dual of this projection is a nonnegative least-squares problem with a
## fixed design matrix D = chol(VCOV) %*% t(Amat):
##
##    min_{lambda >= 0} || z_tilde - D lambda ||^2,   z_tilde ~ N(0, I_p),
##
## solved by Lawson-Hanson (nnls). The number of positive multipliers equals
## the number of active inequality constraints (what solve.QP reports as
## length(iact)), so dimL = p - meq - #positive. This is much faster per draw
## than solve.QP + rtmvnorm and, unlike the primal form based on
## chol(Amat VCOV t(Amat)), also works when Amat has more rows than columns
## or contains redundant rows. Equality rows get an unconstrained multiplier:
## their columns are projected out of D and z_tilde (QR residualization)
## before the NNLS solve.
##
## Sign convention: the response z_tilde corresponds to projecting
## z = -t(U) %*% z_tilde (the linear term of the dual flips sign). Since
## z_tilde is symmetrically distributed and drawn in antithetic pairs, the
## distribution of the counts -- and hence the weights -- is unaffected.
## This only matters when pairing individual draws with solve.QP results;
## with that pairing both engines give identical counts draw for draw.
##
## Antithetic pairs (z, -z): the pair members are negatively dependent, which
## markedly reduces the Monte Carlo variance of the expected level (the goric
## penalty). For individual level probabilities the effect varies by level
## (mostly a reduction, occasionally a slight increase); the convergence
## check below therefore estimates the standard errors from the observed
## between-pair variance rather than assuming independent draws.
##
## Arguments passed through ... are meant for rtmvnorm (truncation bounds);
## truncated sampling is not expressible in the NNLS reformulation, so in
## that case the original quadprog engine is used.
con_weights_boot <- function(VCOV, Amat, meq, R = 1e5L,
                             chunk_size = 5000L, convergence_crit = 5e-03,
                             seed = NULL, verbose = FALSE, ...) {

  # --- input validation (both engines) ---
  if (!is.matrix(VCOV) || !is.numeric(VCOV) || nrow(VCOV) != ncol(VCOV)) {
    stop("restriktor ERROR: VCOV must be a square numeric matrix.", call. = FALSE)
  }
  if (max(abs(VCOV - t(VCOV))) > 1e-08 * max(abs(VCOV))) {
    stop("restriktor ERROR: VCOV must be symmetric.", call. = FALSE)
  }
  VCOV <- (VCOV + t(VCOV)) / 2
  Amat <- rbind(Amat)          # allow a vector in case of one constraint
  if (!is.numeric(Amat) || ncol(Amat) != ncol(VCOV)) {
    stop("restriktor ERROR: ncol(Amat) must equal ncol(VCOV).", call. = FALSE)
  }
  if (!is.numeric(meq) || length(meq) != 1L || is.na(meq) ||
      meq < 0 || meq > nrow(Amat) || meq != round(meq)) {
    stop("restriktor ERROR: meq must be a single integer between 0 and nrow(Amat).", call. = FALSE)
  }
  meq <- as.integer(meq)
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || R < 1 ||
      R > .Machine$integer.max || R != round(R)) {
    stop("restriktor ERROR: R must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || is.na(chunk_size) ||
      chunk_size < 1 || chunk_size > .Machine$integer.max ||
      chunk_size != round(chunk_size)) {
    stop("restriktor ERROR: chunk_size must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(convergence_crit) || length(convergence_crit) != 1L ||
      is.na(convergence_crit) || convergence_crit < 0) {
    stop("restriktor ERROR: convergence_crit must be a single nonnegative number.", call. = FALSE)
  }

  ldots <- list(...)
  if (length(ldots) > 0L) {
    return(con_weights_boot_quadprog(VCOV = VCOV, Amat = Amat, meq = meq,
                                     R = R, chunk_size = chunk_size,
                                     convergence_crit = convergence_crit,
                                     seed = seed, verbose = verbose, ...))
  }

  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  p <- ncol(VCOV)
  g <- nrow(Amat)
  tol_active <- 1e-10          # threshold for counting a multiplier as positive

  U <- try(chol(VCOV), silent = TRUE)  # VCOV = U'U; PD required, as for solve.QP
  if (inherits(U, "try-error")) {
    stop("restriktor ERROR: VCOV must be positive definite.", call. = FALSE)
  }
  D <- U %*% t(Amat)           # p x g; crossprod(D) = Amat VCOV t(Amat)

  qr_eq  <- NULL
  D_ineq <- D
  if (meq > 0L) {
    qr_eq <- qr(D[, seq_len(meq), drop = FALSE])
    if (qr_eq$rank < meq) {
      stop("restriktor ERROR: the first meq rows of Amat (equality constraints) are linearly dependent.", call. = FALSE)
    }
    D_ineq <- qr.resid(qr_eq, D[, -seq_len(meq), drop = FALSE])
  }
  n_ineq <- g - meq

  # inner function: returns number of active constraints, or NA on error
  fn <- function(b) {
    fit <- try(nnls::nnls(D_ineq, Z[, b]), silent = TRUE)

    if (inherits(fit, "try-error")) {
      if (verbose) cat("nnls FAILED\n")
      return(NA_integer_)
    }
    iact <- meq + sum(fit$x > tol_active)
    if (verbose) {
      cat(" ...number of active constraints =", iact, "\n")
    }
    iact
  }

  # between-pair bookkeeping for the Monte Carlo standard errors: per level,
  # running sums of the pair averages a = (1{d1 = level} + 1{d2 = level}) / 2
  # and their squares, over complete antithetic pairs
  pair_sum   <- numeric(p + 1)
  pair_sumsq <- numeric(p + 1)
  n_pairs    <- 0L
  se_wt      <- rep(NA_real_, p + 1)
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

    # lazy sample generation: only generate what this chunk needs.
    # columns (j, j + half) form antithetic pairs; with odd current_n the
    # last generated draw stays unpaired (it still counts toward the weights,
    # only not toward the between-pair standard errors)
    half <- ceiling(current_n / 2)
    Z <- matrix(rnorm(p * half), nrow = p)
    Z <- cbind(Z, -Z)[, seq_len(current_n), drop = FALSE]
    if (meq > 0L) {
      Z <- qr.resid(qr_eq, Z)
    }

    # process chunk with vapply (faster than lapply + sapply post-processing)
    chunk_iact <- if (n_ineq == 0L) {
      rep(as.integer(g), current_n)
    } else {
      vapply(seq_len(current_n), fn, integer(1))
    }

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

    if (n_valid == 0L) {
      stop("restriktor ERROR: all bootstrap draws failed; no chi-bar-square weights can be computed.", call. = FALSE)
    }

    # between-pair bookkeeping over the complete pairs (j, j + half)
    n_pair_chunk <- current_n - half
    if (n_pair_chunk > 0L) {
      d1 <- p - chunk_iact[seq_len(n_pair_chunk)]
      d2 <- p - chunk_iact[half + seq_len(n_pair_chunk)]
      pair_ok <- !is.na(d1) & !is.na(d2)
      d1 <- d1[pair_ok]
      d2 <- d2[pair_ok]
      if (length(d1) > 0L) {
        t1 <- tabulate(d1 + 1L, nbins = p + 1)
        t2 <- tabulate(d2 + 1L, nbins = p + 1)
        same <- d1 == d2
        t_same <- tabulate(d1[same] + 1L, nbins = p + 1)
        # pair average a is 1 if both members hit the level, 0.5 if one does;
        # so sum(a) = (t1 + t2)/2 and sum(a^2) = t_same + (t1 + t2 - 2*t_same)/4
        pair_sum   <- pair_sum + (t1 + t2) / 2
        pair_sumsq <- pair_sumsq + t_same + (t1 + t2 - 2 * t_same) / 4
        n_pairs    <- n_pairs + length(d1)
      }
    }

    wt_bar <- level_counts / n_valid

    # convergence: stop once the 95% confidence half-width of every weight
    # (1.96 times its Monte Carlo standard error, estimated from the
    # between-pair variance) is below convergence_crit
    if (n_pairs >= 2L) {
      var_pair <- (pair_sumsq - pair_sum^2 / n_pairs) / (n_pairs - 1)
      var_pair[var_pair < 0] <- 0
      se_wt <- sqrt(var_pair / n_pairs)
      if (1.96 * max(se_wt) < convergence_crit) {
        has_converged <- TRUE
        chunk_wt_bar  <- rbind(chunk_wt_bar, wt_bar)
        break
      }
    }

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
  attr(wt_bar, "chunk_iter"           ) <- nrow(chunk_wt_bar)
  attr(wt_bar, "error.idx"            ) <- error_idx_all
  attr(wt_bar, "mvtnorm"              ) <- ldots
  attr(wt_bar, "engine"               ) <- "nnls"
  attr(wt_bar, "mc_se"                ) <- se_wt

  return(wt_bar)
}


## original simulation engine (rtmvnorm + solve.QP in theta-space). Kept for
## the case that truncation arguments are passed through ... to rtmvnorm,
## which the NNLS reformulation cannot express.
con_weights_boot_quadprog <- function(VCOV, Amat, meq, R = 1e5L,
                                      chunk_size = 5000L, convergence_crit = 5e-03,
                                      seed = NULL, verbose = FALSE, ...) {

  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  ldots <- list(...)
  Amat <- rbind(Amat) # allow a vector in case of one constraint
  p <- ncol(VCOV)
  se_wt <- rep(NA_real_, p + 1)

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

    if (n_valid == 0L) {
      stop("restriktor ERROR: all bootstrap draws failed; no chi-bar-square weights can be computed.", call. = FALSE)
    }

    wt_bar <- level_counts / n_valid

    # convergence: stop once the 95% confidence half-width of every weight
    # (1.96 times its binomial Monte Carlo standard error; the draws are
    # independent here) is below convergence_crit
    if (n_valid >= 2L) {
      se_wt <- sqrt(wt_bar * (1 - wt_bar) / n_valid)
      if (1.96 * max(se_wt) < convergence_crit) {
        has_converged <- TRUE
        chunk_wt_bar  <- rbind(chunk_wt_bar, wt_bar)
        break
      }
    }

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
  attr(wt_bar, "chunk_iter"           ) <- nrow(chunk_wt_bar)
  attr(wt_bar, "error.idx"            ) <- error_idx_all
  attr(wt_bar, "mvtnorm"              ) <- ldots
  attr(wt_bar, "engine"               ) <- "quadprog"
  attr(wt_bar, "mc_se"                ) <- se_wt

  return(wt_bar)
}
