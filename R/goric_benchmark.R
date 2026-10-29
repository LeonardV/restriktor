benchmark_means  <- function(object, ...) UseMethod("benchmark_means")
benchmark_asymp  <- function(object, ...) UseMethod("benchmark_asymp")


# benchmark <- function(object, ...) {
#   benchmark_asymp(object, ...)
# }

benchmark <- function(object, model_type = c("asymp", "means"), ...) {

  model_type <- match.arg(model_type, c("asymp", "means"))
  if (is.null(model_type)) {
    stop("\nrestriktor ERROR: Please specify if you want to benchmark means or asymptotic result ",
         "In case of model_type = means, no intercept is allowed.")
  }

  if (model_type == "means") {
    # An intercept is allowed when goric() was run on an lm object: the
    # group means are then derived from the parameters (see
    # setup_means_benchmark() and group_means()). Without an lm object, the
    # estimates must be the group means themselves (no intercept); this is
    # checked in setup_means_benchmark().
    benchmark_means(object, ...)
  } else if (model_type == "asymp") {
    benchmark_asymp(object, ...)
  }
}


benchmark_means <- function(object, pop_es = NULL, ratio_pop_means = NULL,
                            group_size = NULL, alt_group_size = NULL,
                            quant = NULL, iter = NULL,
                            control = list(),
                            ncpus = 1, seed = NULL,
                            iter_adequacy_band = c(0.495, 0.505),
                            iter_stability_tol = 1,
                            iter_min = 500, iter_step = 100, iter_max = 2000,
                            # Threshold q used for the 'hypothesis_rate' output field. 
                            # Still allowing a different threshold to be requested at
                            # print time without rerunning the bootstrap.
                            hypo_rate_threshold = 1,
                            # Same mechanism but for the 'rate_rlw' output field
                            # Deliberately not named/called a "hypothesis rate": 
                            # unlike rgw, rlw isn't used (directly) for hypothesis selection.
                            threshold_rlw = 1,
                            # Prior model probabilities used in the (re-)calculated
                            # and simulated goric() calls; NULL (default) means: the
                            # ones of 'object' (object$priorICweights), such that the
                            # benchmarks are based on the same priors as the sample value.
                            priorICweights = NULL, ...) {

  # iter = NULL (the default): start at 500 draws and grow by 100 at a time,
  # up to 2000, stopping as soon as the "Observed" population's benchmark
  # looks adequate -- see run_benchmark_simulation(). A user-supplied numeric
  # 'iter' is used as-is (single fixed-size run, as before).
  user_iter <- iter

  if (length(control) == 1) {
    control <- object$objectList[[1]]$control
  }

  mix_weights <- attr(object$objectList[[1]]$wt.bar, "method")
  penalty_factor <- object$penalty_factor
  # Use the same prior model probabilities as in the goric object (unless
  # specified otherwise), in the re-calculated and the simulated goric() calls.
  if (is.null(priorICweights)) {
    priorICweights <- object$priorICweights
  }

  # Check:
  if (!inherits(object, "con_goric")) {
    stop(paste("\nrestriktor ERROR:",
               "The object should be of class 'con_goric' (a goric(a) object from",
               "the goric() function). However, it belongs to the following class(es):",
               paste(class(object), collapse = ", ")
    ), call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  # If ncpus > 1 but no parallel plan is active yet, set one up for the
  # duration of this call (restored automatically on exit) so the
  # future_lapply() calls below actually run in parallel, instead of ncpus
  # being silently ignored and everything running sequentially regardless.
  # A plan the user already set up themselves (sequential or not) is left
  # untouched.
  if (ncpus > 1 && inherits(future::plan(), "sequential")) {
    current_plan <- future::plan()
    on.exit(future::plan(current_plan), add = TRUE)
    if (.Platform$OS.type == "windows") {
      future::plan(future::multisession, workers = ncpus)
    } else {
      future::plan(future::multicore, workers = ncpus)
    }
  }

  # Hypotheses
  hypos <- object$hypotheses_usr
  nr_hypos <- dim(object$result)[1]
  Heq <- object$Heq
  # Capture add_Hc *before* 'object' is reassigned below (object <- goric(...))
  # -- otherwise it's lost and the recomputed/simulated goric() calls below
  # silently collapse the add_Hc row structure back down to a plain
  # comparison = "none" result (see goric_benchmark.R history for the bug
  # this caused: a dimension mismatch in benchmark()).
  add_Hc <- object$add_Hc

  # Group means, group sizes, and the link between the group means and the
  # parameters: group means = L %*% est. The model may contain an intercept,
  # multiple factors and/or covariates when goric() was run on an lm object;
  # see setup_means_benchmark() (goric_benchmark_utilities.R).
  setup <- setup_means_benchmark(object, group_size = group_size,
                                 alt_group_size = alt_group_size)
  est <- setup$est            # parameter estimates (hypotheses are on these)
  VCOV <- setup$VCOV          # their covariance matrix (based on N, not N-k)
  L <- setup$L
  N <- setup$N                # group sizes (alt_group_size, if specified)
  s2 <- setup$s2              # residual variance per group
  group_means <- setup$group_means
  ngroups <- length(group_means)
  if (!is.null(ratio_pop_means) && length(ratio_pop_means) != ngroups) {
    stop("\nrestriktor ERROR: 'ratio_pop_means' should be of length ", ngroups,
         " (the number of groups: ", paste(names(group_means), collapse = ", "),
         ").", call. = FALSE)
  }
  # The sample gorica(c) value must also be adjusted when alt_group_size is
  # specified; therefore, below, a new goric-object is fitted with est,
  # VCOV, and (if goricac) N.

  # (re-)calculate gorica(c) if:
  # - goric(c) not yet gorica(c)
  # - if adjusted sample size
  # Will do it anyway
  #
  # If needed, adjust type from goric(c) to gorica(c)
  type <- switch(object$type,
                 "goric" = {
                   message("\nrestriktor Message: 'goric' has been converted to 'gorica': the benchmarks are based on ",
                           "the estimates and their covariance matrix. For the penalties of the GORIC(A)(C), ",
                           "see the argument 'type' in ?goric.")
                   "gorica"
                 },
                 "goricc" = {
                   message("\nrestriktor Message: 'goricc' has been converted to 'goricac': the benchmarks are based on ",
                           "the estimates and their covariance matrix. For the penalties of the GORIC(A)(C), ",
                           "see the argument 'type' in ?goric.")
                   "goricac"
                 },
                 object$type)
  #
  # (re-)calculate gorica(c)
  object <-
    goric(
      est,
      CI = FALSE, # no CIs needed for the re-calculated object
      VCOV = VCOV, # based on N not N-k
      sample_nobs = sum(N), # Needed for type = "goricac" - daar genoeg aan sum(N) of moet het juist gehele N hebben?
      # number of nuisance parameters of the original model (type = "goricac")
      nr_nuis_param = nr_nuis_param_object(object),
      hypotheses = hypos,
      comparison = object$comparison,
      type = type,
      control = control,
      mix_weights = mix_weights,
      penalty_factor = penalty_factor,
      priorICweights = priorICweights,
      #Heq = FALSE,
      Heq = Heq,
      add_Hc = add_Hc,
      # already checked once on the user's original add_Hc fit; see the
      # NOTE above parallel_benchmark_draw()'s internal goric() call in
      # goric_benchmark_utilities.R for the full rationale.
      check_add_Hc = FALSE,
      ...
    )

  # Pre-parsed form of 'hypos', reused for every draw inside the Monte Carlo
  # loop below instead of re-parsing the same constraint-syntax text from
  # scratch on each one -- see precompute_hypos_for_simulation()'s comment
  # (goric_benchmark_utilities.R) for why this is safe and how much it
  # saves. 'hypos' itself (text) is left untouched and still used for the
  # one-off calculate_error_probability() call below.
  hypos_for_sim <- precompute_hypos_for_simulation(object, hypos, Heq)

  ## Compute observed Cohens f (of the group means)
  cohens_f_observed <- compute_cohens_f_s2(group_means, N, s2)

  # effect size population
  if (is.null(pop_es)) {
    pop_es <- c(0, cohens_f_observed)
    names(pop_es) <- c("No effect", "Observed")
  } else {
    pop_es <- sort(pop_es)
  }

  # Assign row names to pop_es if they are null or empty strings
  rnames <- names(pop_es)
  if (is.null(rnames)) {
    rnames <- paste0("PES_", seq_len(length(pop_es)))
    names(pop_es) <- rnames
  } else {
    empty_names <- rnames == ""
    if (any(empty_names)) {
      rnames[empty_names] <- paste0("PE_", seq_len(sum(empty_names)))
      names(pop_es) <- rnames
    }
  }

  es <- pop_es
  nr_es <- length(es)

  # Population values of the parameters (used in the simulation) such that
  # the implied group means have effect size pop_es and the same pattern as
  # the observed group means; see generate_pop_parameters().
  est_pop_all <- t(sapply(pop_es, function(x) {
    generate_pop_parameters(est, L, N, s2, target_f = x,
                            ratio_pop_means = ratio_pop_means)
  }))
  colnames(est_pop_all) <- names(est)
  rownames(est_pop_all) <- paste0("pop_es = ", pop_es)
  # implied population group means
  means_pop_all <- est_pop_all %*% t(L)
  colnames(means_pop_all) <- names(group_means)

  # preferred hypothesis
  pref_hypo <- which.max(object$result[, 7])
  pref_hypo_name <- object$result$model[pref_hypo]

  if (is.null(quant)) {
    quant <- c(.05, .35, .50, .65, .95)
    names_quant <- c("Sample", "5%", "35%", "50%", "65%", "95%")
  } else {
    names_quant <- c("Sample", paste0(as.character(quant*100), "%"))
  }

  # The draws use the same type (gorica or goricac) as the sample value.
  sim <- run_benchmark_simulation(
    type = type, sample_nobs = sum(N),
    nr_nuis_param = nr_nuis_param_object(object),
    nr_es = nr_es, rnames = rnames, name_prefix = "pop_es = ",
    center_matrix = est_pop_all, colnames_vec = names(est),
    VCOV = VCOV, hypos = hypos_for_sim, pref_hypo = pref_hypo,
    comparison = object$comparison, control = control,
    mix_weights = mix_weights, penalty_factor = penalty_factor, Heq = Heq,
    priorICweights = priorICweights,
    object = object, iter = user_iter,
    es_labels = paste0(es, " (", names(es), ")"),
    band = iter_adequacy_band,
    stability_tol = iter_stability_tol,
    iter_min = iter_min, iter_step = iter_step, iter_max = iter_max,
    add_Hc = add_Hc,
    ...
  )
  parallel_function_results <- sim$parallel_function_results
  iter <- sim$iter # final number of draws actually used, per category

  # get benchmark results
  benchmark_results <- get_results_benchmark(parallel_function_results,
                                             object, pref_hypo,
                                             pref_hypo_name, quant,
                                             names_quant, nr_hypos,
                                             hypo_rate_threshold = hypo_rate_threshold,
                                             threshold_rlw = threshold_rlw)

  if (!is.null(user_iter)) {
    # Fixed 'iter': run_benchmark_simulation() does not auto-grow or message
    # in this case, so do the (non-growing) adequacy check here instead.
    bias_check <- check_iter_adequacy(benchmark_results, "pop_es = Observed", iter,
                        band = iter_adequacy_band,
                        iter_min = iter_min, iter_step = iter_step, iter_max = iter_max,
                        stability_tol = iter_stability_tol,
                        control = control)
  } else {
    bias_check <- list(median_bias_check_gw = sim$median_bias_check_gw,
                       median_bias_check_lw = sim$median_bias_check_lw)
  }

  # compute error probability
  error_prob <- calculate_error_probability(object, hypos, pref_hypo,
                                            est = est,
                                            VCOV, control, ...)

  # residual error variance
  # var_e <- as.vector(diag(VCOV)[1])
  # var_e_data <- as.vector(diag(VCOV_orig)[1])
  # This is used as reference value, because the diag. elements will differ if group sizes differ
  
  OUT <- list(
    type = object$type,
    comparison = object$comparison,
    ngroups = ngroups,
    group_size = group_size_matrix(N, setup$definitions),
    group_means_observed = group_means,
    #ratio_group_means.data = ratio_data,
    cohens_f_observed = cohens_f_observed,
    #res_var_observed = var_e_data,
    pop_es = pop_es, 
    pop_group_means = means_pop_all,
    pop_parameters = est_pop_all,
    group_means_definitions = setup$definitions,
    continuous_vars = setup$continuous_vars,
    ratio_pop_means = ratio_pop_means,
    #res_var_pop = var_e,
    pref_hypo_name = pref_hypo_name,
    error_prob_pref_hypo = error_prob,
    # Grouped by family (rather than one flat field per output_type x family)
    # so the object doesn't sprawl into ~40 top-level names -- each family is
    # one list, indexed by output_type ("goric_weights", "ratio_goric_weights",
    # etc.), and each output_type's value is, as before, a list indexed by
    # pop_es/pop_est category, giving a matrix with one row per alternative
    # hypothesis (self-comparison row dropped) for the "matrix" output_types
    # (ratio_goric_weights, ratio_ll_weights, ratio_ll_weights_ge1,
    # ratio_goric_weights_log, ratio_ll_weights_log, difLL, absdifLL) --
    # nesting the family this way does not change that row structure at all,
    # so this still works the same with 3+ hypotheses (multiple pairwise
    # comparisons -> multiple rows) as it did with the flat field names.
    benchmarks = list(
      goric_weights = benchmark_results$benchmarks_gw,
      ll_weights = benchmark_results$benchmarks_lw,
      ratio_goric_weights = benchmark_results$benchmarks_rgw,
      ratio_ll_weights = benchmark_results$benchmarks_rlw,
      ratio_ll_weights_ge1 = benchmark_results$benchmarks_rlw_ge1,
      ratio_goric_weights_log = benchmark_results$benchmarks_rgw_log,
      ratio_ll_weights_log = benchmark_results$benchmarks_rlw_log,
      difLL = benchmark_results$benchmarks_difLL,
      absdifLL = benchmark_results$benchmarks_absdifLL
    ),
    combined_values = benchmark_results$combined_values,
    #
    pctl_Sample = list(
      goric_weights = benchmark_results$pctl_Sample_gw,
      ll_weights = benchmark_results$pctl_Sample_lw,
      ratio_goric_weights = benchmark_results$pctl_Sample_rgw,
      ratio_ll_weights = benchmark_results$pctl_Sample_rlw,
      ratio_ll_weights_ge1 = benchmark_results$pctl_Sample_rlw_ge1,
      ratio_goric_weights_log = benchmark_results$pctl_Sample_rgw_log,
      ratio_ll_weights_log = benchmark_results$pctl_Sample_rlw_log,
      difLL = benchmark_results$pctl_Sample_difLL,
      absdifLL = benchmark_results$pctl_Sample_absdifLL
    ),
    #
    pctl_medianRefPop = list(
      goric_weights = benchmark_results$pctl_medianRefPop_gw,
      ll_weights = benchmark_results$pctl_medianRefPop_lw,
      ratio_goric_weights = benchmark_results$pctl_medianRefPop_rgw,
      ratio_ll_weights = benchmark_results$pctl_medianRefPop_rlw,
      ratio_ll_weights_ge1 = benchmark_results$pctl_medianRefPop_rlw_ge1,
      ratio_goric_weights_log = benchmark_results$pctl_medianRefPop_rgw_log,
      ratio_ll_weights_log = benchmark_results$pctl_medianRefPop_rlw_log,
      difLL = benchmark_results$pctl_medianRefPop_difLL,
      absdifLL = benchmark_results$pctl_medianRefPop_absdifLL
    ),
    #
    # No ratio_ll_weights_ge1/absdifLL here -- overlap (vs. the reference
    # population) was never computed for those two folded/"ge1"-"ge0"
    # variants, same as before this restructuring. 'reference' (which
    # population "with <name>" refers to) lives here too now, rather than as
    # its own stray top-level field, since it's overlap-family metadata.
    overlap = list(
      goric_weights = benchmark_results$overlap_gw,
      ll_weights = benchmark_results$overlap_lw,
      ratio_goric_weights = benchmark_results$overlap_rgw,
      ratio_ll_weights = benchmark_results$overlap_rlw,
      ratio_goric_weights_log = benchmark_results$overlap_rgw_log,
      ratio_ll_weights_log = benchmark_results$overlap_rlw_log,
      difLL = benchmark_results$overlap_ld,
      reference = benchmark_results$overlap_reference_pop
    ),
    #
    hypothesis_rate = benchmark_results$hypothesis_rate,
    hypo_rate_threshold = benchmark_results$hypo_rate_threshold,
    rate_rlw = benchmark_results$rate_rlw,
    threshold_rlw = benchmark_results$threshold_rlw,
    hypothesis_rate_band = benchmark_results$hypothesis_rate_band,
    rate_rlw_band = benchmark_results$rate_rlw_band,
    #
    median_bias_check_gw = bias_check$median_bias_check_gw,
    median_bias_check_lw = bias_check$median_bias_check_lw,
    #
    iter = iter
  )

  class(OUT) <- c("benchmark_means", "benchmark", "list")
  return(OUT)
}




## asymp
benchmark_asymp <- function(object, pop_est = NULL, sample_size = NULL,
                            alt_sample_size = NULL, quant = NULL, iter = NULL,
                            control = list(),
                            ncpus = 1, seed = NULL,
                            iter_adequacy_band = c(0.495, 0.505),
                            iter_stability_tol = 1,
                            iter_min = 500, iter_step = 100, iter_max = 2000,
                            # See benchmark_means()'s 'hypo_rate_threshold' for the full
                            # rationale -- same argument, same default, same
                            # x$hypo_rate_threshold output field.
                            hypo_rate_threshold = 1,
                            # See benchmark_means()'s 'threshold_rlw' for the full
                            # rationale -- same argument, same default, same
                            # x$threshold_rlw output field.
                            threshold_rlw = 1,
                            # Prior model probabilities used in the (re-)calculated
                            # and simulated goric() calls; NULL (default) means: the
                            # ones of 'object' (object$priorICweights), such that the
                            # benchmarks are based on the same priors as the sample value.
                            priorICweights = NULL, ...) {

  # iter = NULL (the default): start at 500 draws and grow by 100 at a time,
  # up to 2000, stopping as soon as the "Observed" population's benchmark
  # looks adequate -- see run_benchmark_simulation(). A user-supplied numeric
  # 'iter' is used as-is (single fixed-size run, as before).
  user_iter <- iter

  # TO DO als in means nog 'group_size' argument dan ms hier zeggen dat dat sample_size moet zijn?
  #       NB alleen nodig als alt_sample_size nodig is...
    
  if (length(control) == 1) {
    control <- object$objectList[[1]]$control
  }
  
  mix_weights <- attr(object$objectList[[1]]$wt.bar, "method")
  penalty_factor <- object$penalty_factor
  # Use the same prior model probabilities as in the goric object (unless
  # specified otherwise), in the re-calculated and the simulated goric() calls.
  if (is.null(priorICweights)) {
    priorICweights <- object$priorICweights
  }
  
  # Check if object is of class con_goric
  if (!inherits(object, "con_goric")) {
    stop(paste("\nrestriktor ERROR:", 
               "The object should be of class 'con_goric' (a goric(a) object from the goric() function).",
               "However, it belongs to the following class(es):", 
               paste(class(object), collapse = ", ")), call. = FALSE)
  }
 
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  
  # If ncpus > 1 but no parallel plan is active yet, set one up for the
  # duration of this call (restored automatically on exit) so the
  # future_lapply() calls below actually run in parallel, instead of ncpus
  # being silently ignored and everything running sequentially regardless.
  # A plan the user already set up themselves (sequential or not) is left
  # untouched.
  if (ncpus > 1 && inherits(future::plan(), "sequential")) {
    current_plan <- future::plan()
    on.exit(future::plan(current_plan), add = TRUE)
    if (.Platform$OS.type == "windows") {
      future::plan(future::multisession, workers = ncpus)
    } else {
      future::plan(future::multicore, workers = ncpus)
    }
  }

  VCOV <- object$VCOV
  # Note that -- assuming an lm object was used -- VCOV is the unbiased cov.mx estimate.
  # It is also mentioned in tutorials, so if user specified it, they could have made this asjustment....
  hypos <- object$hypotheses_usr
  
  if (is.null(pop_est)) {
    check_rhs_constants(object$rhs)
  }
  
  nr_hypos <- dim(object$result)[1]
  Heq <- object$Heq
  comparison <- object$comparison
  # Capture add_Hc *before* 'object' is reassigned below (object <- goric(...))
  # -- otherwise it's lost and the recomputed/simulated goric() calls below
  # silently collapse the add_Hc row structure back down to a plain
  # comparison = "none" result (see goric_benchmark.R history for the bug
  # this caused: a dimension mismatch in benchmark()).
  add_Hc <- object$add_Hc
  type <- object$type
  # (for an mlm, the coefficient matrix becomes a named vector, in the order
  # of VCOV; see mlm_coef_vector())
  est_sample <- mlm_coef_vector(object$b.unrestr)
  n_coef <- length(est_sample)
  
  if (is.null(pop_est)) {
    #NE <- theta_restricted(theta = est_sample, V = VCOV, R = object$constraints[[1]], rhs = object$rhs[[1]])
    # TO determine pop_est, we need the constraint matrix for the preferred hypothesis.
    # In that constraint matrix / hypothesis, the inequalities will be set to equalities.
    # First determine the preferred hypothesis:
    NE_hypo <- which.max(object$result[, 7])

    # object$constraints/$rhs only have one entry per *specified* hypothesis
    # (one per name in object$objectNames) -- not one for an implicit
    # safeguard row. If the preferred row is such a safeguard, use:
    #  - add_Hc: the add_Hc reference hypothesis itself,
    #    since the complement is defined relative to it.
    #  - comparison = "complement": the (only) specified hypothesis.
    #  - comparison = "unconstrained": the best-supported (highest
    #    GORIC(A) weight) *among the specified* hypotheses,
    #    rendering a message since this may not be interesting for the user.
    if (NE_hypo > length(object$constraints)) {
      if (!is.null(object$add_Hc)) {
        fallback_name <- object$add_Hc
        #safeguard_desc <- paste0("the complement of ", sQuote(fallback_name))
        #fallback_reason <- "the add_Hc reference hypothesis"
      } else if (identical(object$comparison, "complement")) {
        fallback_name <- object$objectNames[1]
        #safeguard_desc <- "the complement (of that one hypothesis)"
        #fallback_reason <- "the only specified hypothesis"
      } else {
        best_among_specified <- which.max(object$result[seq_along(object$constraints), 7])
        fallback_name <- object$objectNames[best_among_specified]
        #safeguard_desc <- "the unconstrained model"
        #fallback_reason <- "the best-supported (highest GORIC(A) weight) of the specified hypotheses"
        message(
          "\nrestriktor Message: The preferred hypothesis is the unconstrained.",
          "The 'No inequalities' population is based on the best hypothesis of the specified ones.",
          "Supply your own 'pop_est' argument if you want a different 'No inequalities' population estimate."
        )
      }
      # # Message for all three cases, but only for unc of interest, so done there/above
      # message(paste0(
      #   "\nrestriktor Message: The preferred hypothesis/model is the safeguard (",
      #   safeguard_desc, "), which has no constraint matrix of its own to derive ",
      #   "a 'No-effect' population from. Falling back to ", sQuote(fallback_name),
      #   " (", fallback_reason, ") for this purpose instead. Supply your own ",
      #   "'pop_est' argument if you want a different 'No-effect' population estimate."
      # ))
      NE_hypo <- which(object$objectNames == fallback_name)
    }

    NE <- theta_restricted(theta = est_sample, V = VCOV, R = object$constraints[[NE_hypo]], rhs = object$rhs[[NE_hypo]])
    # Note that VCOV is the unbiased cov.mx estimate.
    pop_est <- matrix(rbind(NE, est_sample), nrow = 2)
    # Named "No inequalities" rather than "No effect": this population is
    # constructed by replacing the inequality constraints of the best
    # hypothesis by equalities (see theta_restricted() above), which is not
    # the same as a true absence of effect -- the equality-constrained point
    # can still represent a real, nonzero effect.
    row.names(pop_est) <- c("No inequalities", "Observed")
  } else {
    if (is.vector(pop_est)) {
      pop_est <- matrix(pop_est, nrow = 1, ncol = length(pop_est))
    }
  }
  
  if (ncol(pop_est) != length(est_sample)) {
    stop(paste("\nrestriktor ERROR:The number of columns in pop_est (", ncol(pop_est), 
               ") does not match the length of est_sample (", length(est_sample), ").", sep = ""), .call = FALSE)
  }
  
  rnames <- row.names(pop_est)
  if (is.null(rnames)) {
    rnames <- paste0("PE_", seq_len(nrow(pop_est)))
    row.names(pop_est) <- rnames
  } else {
    empty_names <- rnames == ""
    if (any(empty_names)) {
      rnames[empty_names] <- paste0("PE_", seq_len(sum(empty_names)))
      row.names(pop_est) <- rnames
    }
  }  
  
  colnames(pop_est) <- names(est_sample)
  N <- object$sample_nobs #length(object$model.org$residuals)
  
  # modeltype
  type <- switch(object$type,
                 "goric" = {
                   message("\nrestriktor Message: 'goric' has been converted to 'gorica': the benchmarks are based on ",
                           "the estimates and their covariance matrix. For the penalties of the GORIC(A)(C), ",
                           "see the argument 'type' in ?goric.")
                   "gorica"
                 },
                 "goricc" = {
                   message("\nrestriktor Message: 'goricc' has been converted to 'goricac': the benchmarks are based on ",
                           "the estimates and their covariance matrix. For the penalties of the GORIC(A)(C), ",
                           "see the argument 'type' in ?goric.")
                   "goricac"
                 },
                 object$type)
  
  
  # Controleer op alternatieve steekproefgrootte
  if (!is.null(alt_sample_size)) {
    # Controleer of de originele steekproefgrootte beschikbaar is
    if (is.null(N) || N == 0) {
      if (is.null(sample_size)) {
        stop("\nrestriktor ERROR: Please provide the original sample size(s) using the argument `sample_size`.", .call = FALSE)
      }
      N <- sample_size
    }
    VCOV <- VCOV * N / alt_sample_size
    N <- alt_sample_size
  } 

  # Herbereken met goric-functie
  object <- goric(
    est_sample,
    CI = FALSE, # no CIs needed for the re-calculated object
    VCOV = VCOV,
    sample_nobs = N[1], # Needed for type = "goricac"
    # number of nuisance parameters of the original model (type = "goricac")
    nr_nuis_param = nr_nuis_param_object(object),
    # TO DO, dit is toch maar een getal en als niet dan ws de sum (iig bij anova wel)
    #sample_nobs = sum(N), # Needed for type = "goricac" - daar genoeg aan sum(N) of moet het juist gehele N hebben?
    hypotheses = hypos,
    comparison = comparison,
    type = type,
    control = control,
    mix_weights = mix_weights,
    penalty_factor = penalty_factor,
    priorICweights = priorICweights,
    Heq = Heq,
    add_Hc = add_Hc,
    # already checked once on the user's original add_Hc fit; see the NOTE
    # above parallel_benchmark_draw()'s internal goric() call in
    # goric_benchmark_utilities.R for the full rationale.
    check_add_Hc = FALSE,
    ...
  )

  # Pre-parsed form of 'hypos', reused for every draw inside the Monte Carlo
  # loop below instead of re-parsing the same constraint-syntax text from
  # scratch on each one -- see precompute_hypos_for_simulation()'s comment
  # (goric_benchmark_utilities.R) for why this is safe and how much it
  # saves. 'hypos' itself (text) is left untouched and still used for the
  # one-off calculate_error_probability() call below.
  hypos_for_sim <- precompute_hypos_for_simulation(object, hypos, Heq)

  if (is.null(quant)) {
    quant <- c(.05, .35, .50, .65, .95)
    names_quant <- c("Sample", "5%", "35%", "50%", "65%", "95%")
  } else {
    names_quant <- c("Sample", paste0(as.character(quant*100), "%"))
  }
  
  pref_hypo <- which.max(object$result[, 7])
  pref_hypo_name <- object$result$model[pref_hypo]
  
  nr_es  <- nrow(pop_est)

  # The draws use the same type (gorica or goricac) as the sample value.
  sim <- run_benchmark_simulation(
    type = type, sample_nobs = N[1],
    nr_nuis_param = nr_nuis_param_object(object),
    nr_es = nr_es, rnames = rnames, name_prefix = "pop_est = ",
    center_matrix = pop_est, colnames_vec = names(est_sample),
    VCOV = VCOV, hypos = hypos_for_sim, pref_hypo = pref_hypo,
    comparison = comparison, control = control,
    mix_weights = mix_weights, penalty_factor = penalty_factor, Heq = Heq,
    priorICweights = priorICweights,
    object = object, iter = user_iter,
    band = iter_adequacy_band,
    stability_tol = iter_stability_tol,
    iter_min = iter_min, iter_step = iter_step, iter_max = iter_max,
    add_Hc = add_Hc,
    ...
  )
  parallel_function_results <- sim$parallel_function_results
  iter <- sim$iter # final number of draws actually used, per category

  benchmark_results <- get_results_benchmark(parallel_function_results, object, pref_hypo,
                                             pref_hypo_name, quant, names_quant, nr_hypos,
                                             hypo_rate_threshold = hypo_rate_threshold,
                                             threshold_rlw = threshold_rlw)

  if (!is.null(user_iter)) {
    # Fixed 'iter': run_benchmark_simulation() does not auto-grow or message
    # in this case, so do the (non-growing) adequacy check here instead.
    bias_check <- check_iter_adequacy(benchmark_results, "pop_est = Observed", iter,
                        band = iter_adequacy_band,
                        iter_min = iter_min, iter_step = iter_step, iter_max = iter_max,
                        stability_tol = iter_stability_tol,
                        control = control)
  } else {
    bias_check <- list(median_bias_check_gw = sim$median_bias_check_gw,
                       median_bias_check_lw = sim$median_bias_check_lw)
  }

  error_prob <- calculate_error_probability(object, hypos, pref_hypo,
                                            est = est_sample, VCOV, control, ...)
  
  OUT <- list(
    type = object$type,
    comparison = object$comparison,
    n_coef = n_coef,
    sample_size = N,
    pop_est = pop_est, 
    pop_VCOV = VCOV,
    pref_hypo_name = pref_hypo_name, 
    error_prob_pref_hypo = error_prob,
    # Grouped by family (rather than one flat field per output_type x family)
    # so the object doesn't sprawl into ~40 top-level names -- each family is
    # one list, indexed by output_type ("goric_weights", "ratio_goric_weights",
    # etc.), and each output_type's value is, as before, a list indexed by
    # pop_es/pop_est category, giving a matrix with one row per alternative
    # hypothesis (self-comparison row dropped) for the "matrix" output_types
    # (ratio_goric_weights, ratio_ll_weights, ratio_ll_weights_ge1,
    # ratio_goric_weights_log, ratio_ll_weights_log, difLL, absdifLL) --
    # nesting the family this way does not change that row structure at all,
    # so this still works the same with 3+ hypotheses (multiple pairwise
    # comparisons -> multiple rows) as it did with the flat field names.
    benchmarks = list(
      goric_weights = benchmark_results$benchmarks_gw,
      ll_weights = benchmark_results$benchmarks_lw,
      ratio_goric_weights = benchmark_results$benchmarks_rgw,
      ratio_ll_weights = benchmark_results$benchmarks_rlw,
      ratio_ll_weights_ge1 = benchmark_results$benchmarks_rlw_ge1,
      ratio_goric_weights_log = benchmark_results$benchmarks_rgw_log,
      ratio_ll_weights_log = benchmark_results$benchmarks_rlw_log,
      difLL = benchmark_results$benchmarks_difLL,
      absdifLL = benchmark_results$benchmarks_absdifLL
    ),
    combined_values = benchmark_results$combined_values,
    #
    pctl_Sample = list(
      goric_weights = benchmark_results$pctl_Sample_gw,
      ll_weights = benchmark_results$pctl_Sample_lw,
      ratio_goric_weights = benchmark_results$pctl_Sample_rgw,
      ratio_ll_weights = benchmark_results$pctl_Sample_rlw,
      ratio_ll_weights_ge1 = benchmark_results$pctl_Sample_rlw_ge1,
      ratio_goric_weights_log = benchmark_results$pctl_Sample_rgw_log,
      ratio_ll_weights_log = benchmark_results$pctl_Sample_rlw_log,
      difLL = benchmark_results$pctl_Sample_difLL,
      absdifLL = benchmark_results$pctl_Sample_absdifLL
    ),
    #
    pctl_medianRefPop = list(
      goric_weights = benchmark_results$pctl_medianRefPop_gw,
      ll_weights = benchmark_results$pctl_medianRefPop_lw,
      ratio_goric_weights = benchmark_results$pctl_medianRefPop_rgw,
      ratio_ll_weights = benchmark_results$pctl_medianRefPop_rlw,
      ratio_ll_weights_ge1 = benchmark_results$pctl_medianRefPop_rlw_ge1,
      ratio_goric_weights_log = benchmark_results$pctl_medianRefPop_rgw_log,
      ratio_ll_weights_log = benchmark_results$pctl_medianRefPop_rlw_log,
      difLL = benchmark_results$pctl_medianRefPop_difLL,
      absdifLL = benchmark_results$pctl_medianRefPop_absdifLL
    ),
    #
    # No ratio_ll_weights_ge1/absdifLL here -- overlap (vs. the reference
    # population) was never computed for those two folded/"ge1"-"ge0"
    # variants, same as before this restructuring. 'reference' (which
    # population "with <name>" refers to) lives here too now, rather than as
    # its own stray top-level field, since it's overlap-family metadata.
    overlap = list(
      goric_weights = benchmark_results$overlap_gw,
      ll_weights = benchmark_results$overlap_lw,
      ratio_goric_weights = benchmark_results$overlap_rgw,
      ratio_ll_weights = benchmark_results$overlap_rlw,
      ratio_goric_weights_log = benchmark_results$overlap_rgw_log,
      ratio_ll_weights_log = benchmark_results$overlap_rlw_log,
      difLL = benchmark_results$overlap_ld,
      reference = benchmark_results$overlap_reference_pop
    ),
    #
    hypothesis_rate = benchmark_results$hypothesis_rate,
    hypo_rate_threshold = benchmark_results$hypo_rate_threshold,
    rate_rlw = benchmark_results$rate_rlw,
    threshold_rlw = benchmark_results$threshold_rlw,
    hypothesis_rate_band = benchmark_results$hypothesis_rate_band,
    rate_rlw_band = benchmark_results$rate_rlw_band,
    #
    median_bias_check_gw = bias_check$median_bias_check_gw,
    median_bias_check_lw = bias_check$median_bias_check_lw,
    #
    iter = iter
  )

  class(OUT) <- c("benchmark_asymp", "benchmark", "list")
  return(OUT)
}
