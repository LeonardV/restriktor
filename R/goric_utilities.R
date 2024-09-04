coef.con_goric <- function(object, ...)  {
  return(object$ormle$b.restr)
}

coef.gorica_est <- function(object, ...)  {
  return(object$b.restr)
}

calculate_model_comparison_metrics <- function(x) {
  modelnames <- as.character(x$model)
  ## Log-likelihood
  LL = -2 * x$loglik
  delta_LL = LL - min(LL)
  loglik_weights = exp(0.5 * -delta_LL) / sum(exp(0.5 * -delta_LL))
  loglik_rw = loglik_weights %*% t(1/loglik_weights)
  diag(loglik_rw) = 1
  
  ## penalty
  penalty_weights = exp(-x$penalty) / sum(exp(-x$penalty))
  penalty_rw = penalty_weights %*% t(1/penalty_weights)
  diag(penalty_rw) = 1
  
  ## goric
  delta_goric = x$goric - min(x$goric)
  goric_weights = exp(0.5 * -delta_goric) / sum(exp(0.5 * -delta_goric))
  goric_rw = goric_weights %*% t(1/goric_weights)
  diag(goric_rw) = 1
  
  # if user specified hypotheses is >= 2 and comparison = unconstrained
  # add extra column with goric weights excluding unconstrained model.
  mn_unc_idx <- grep("unconstrained", modelnames)
  if (length(modelnames) > 2 && length(mn_unc_idx) > 0 && 
      which.max(goric_weights) != mn_unc_idx) {
    delta_goric = x$goric[-mn_unc_idx] - min(x$goric[mn_unc_idx])
    goric_weights_without_unc = exp(0.5 * -delta_goric) / 
      sum(exp(0.5 * -delta_goric))
    goric_weights_without_unc <- c(goric_weights_without_unc, NA)
  } else { goric_weights_without_unc <- NULL }
  
  rownames(goric_rw) = modelnames
  rownames(penalty_rw) = modelnames
  rownames(loglik_rw) = modelnames
  colnames(goric_rw) = paste0("vs. ", modelnames)
  colnames(penalty_rw) = paste0("vs. ", modelnames)
  colnames(loglik_rw) = paste0("vs. ", modelnames)
  
  out <- list(loglik_weights = loglik_weights, 
              penalty_weights = penalty_weights,
              goric_weights = goric_weights,
              goric_weights_without_unc = goric_weights_without_unc,
              loglik_rw = loglik_rw,
              penalty_rw = penalty_rw,
              goric_rw = goric_rw)
  
  return(out)
}

# compute penalty term, where range restrictions are treated as ceq.
PT_Amat_meq <- function(Amat, meq) {
  PT_meq  <- meq
  # check for linear dependence
  RREF <- GaussianElimination(t(Amat)) # qr(Amat)$rank
  # remove linear dependent rows
  PT_Amat <- Amat[RREF$pivot, , drop = FALSE] 
  
  if (nrow(Amat) > 1) {
    # check for range restrictions, e.g., -1 < beta < 1
    idx_range_restrictions <- detect_range_restrictions(Amat)
    # range restrictions are treated as equalities for computing PT (goric)
    n_range_restrictions <- nrow(idx_range_restrictions)
    PT_meq <- meq + n_range_restrictions
    # reorder PT_Amat: ceq first, ciq second, needed for QP.solve()
    meq_order_idx <- RREF$pivot %in% c(idx_range_restrictions)
    PT_Amat <- rbind(PT_Amat[meq_order_idx, ], PT_Amat[!meq_order_idx, ])
  }
  
  return(list(PT_meq = PT_meq, PT_Amat = PT_Amat, RREF = RREF))
}


# Create a function to sort the elements in each string
sort_combination <- function(combination) {
  split_combination <- strsplit(combination, " vs. ")[[1]]
  #sorted_combination <- sort(split_combination)
  # Check if the combination includes "complement"
  if ("complement" %in% split_combination) {
    # Find the other element that is not "complement"
    other_element <- split_combination[split_combination != "complement"]
    sorted_combination <- c(other_element, "complement")
  } else if ("unconstrained" %in% split_combination) {
    # Find the other element that is not "unconstrained"
    other_element <- split_combination[split_combination != "unconstrained"]
    sorted_combination <- c(other_element, "unconstrained")
  } else {
    sorted_combination <- sort(split_combination)
  }
  paste(sorted_combination, collapse = " vs. ")
}


calculate_weight_bar <- function(Amat, meq, VCOV, mix_weights, seed, control,
                                 verbose, ...) {
  wt.bar <- NA
  if (nrow(Amat) == meq) {
    # equality constraints only
    wt.bar <- rep(0L, ncol(VCOV) + 1)
    wt.bar.idx <- ncol(VCOV) - qr(Amat)$rank + 1
    wt.bar[wt.bar.idx] <- 1
  } else if (all(c(Amat) == 0)) { 
    # unrestricted case
    wt.bar <- c(rep(0L, ncol(VCOV)), 1)
  } else if (mix_weights == "boot") { 
    # compute chi-square-bar weights based on Monte Carlo simulation
    wt.bar <- con_weights_boot(VCOV = VCOV,
                               Amat = Amat, 
                               meq  = meq, 
                               R    = ifelse(is.null(control$mix_weights_bootstrap_limit), 
                                             1e5L, control$mix_weights_bootstrap_limit),
                               seed = seed,
                               convergence_crit = ifelse(is.null(control$convergence_crit), 
                                                         1e-03, control$convergence_crit),
                               chunk_size = ifelse(is.null(control$chunk_size), 
                                                   5000L, control$chunk_size),
                               verbose = verbose, ...)
    attr(wt.bar, "mix_weights_bootstrap_limit") <- control$mix_weights_bootstrap_limit 
  } else if (mix_weights == "pmvnorm" && meq < nrow(Amat)) {
    # compute chi-square-bar weights based on pmvnorm
    wt.bar <- rev(con_weights(Amat %*% VCOV %*% t(Amat), meq = meq, 
                              tolerance = ifelse(is.null(control$tolerance), 1e-15, control$tolerance), 
                              ridge_constant = ifelse(is.null(control$ridge_constant), 1e-05, control$ridge_constant), 
                              ...))
    
    # Check if wt.bar contains NaN values
    if (any(is.nan(wt.bar))) {
      mix_weights <- "boot"
      wt.bar <- con_weights_boot(VCOV = VCOV,
                                 Amat = Amat, 
                                 meq  = meq, 
                                 R    = ifelse(is.null(control$mix_weights_bootstrap_limit), 
                                               1e5L, control$mix_weights_bootstrap_limit),
                                 seed = seed,
                                 convergence_crit = ifelse(is.null(control$convergence_crit), 
                                                           1e-03, control$convergence_crit),
                                 chunk_size = ifelse(is.null(control$chunk_size), 
                                                     5000L, control$chunk_size),
                                 verbose = verbose, ...)
      attr(wt.bar, "mix_weights_bootstrap_limit") <- control$mix_weights_bootstrap_limit   
    } 
  }
  return(wt.bar)
}


