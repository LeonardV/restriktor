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
    # Check if the model has an intercept. Since the goric function accepts both
    # a model/fit or only estimates+VCOV, we cannot rely on the original model fit.
    # So we only check if the names of the vector with parameter estimates includes
    # the word \(Intercept\).
    has_intercept <- detect_intercept(object)
    if (has_intercept) {
      stop("\nrestriktor ERROR: A model with an intercept is not allowed for model_type = means. ",
           "Please refit the model without an intercept.")
    }
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
                            # Threshold q used for the 'hypothesis_rate' output field (rate
                            # at which each alternative hypothesis's ratio-GORIC(A)-weight
                            # bootstrap draws exceed q, i.e. how often that alternative is
                            # preferred over the preferred hypothesis). Stored on the
                            # returned object as x$hypo_rate_threshold, which
                            # print.benchmark()'s own 'hypo_rate_threshold' argument
                            # defaults to -- so a print(x) matches what was used here,
                            # while still allowing a different threshold to be requested at
                            # print time (see print.benchmark()) without rerunning the
                            # bootstrap.
                            hypo_rate_threshold = 1,
                            # Same mechanism as 'hypo_rate_threshold' above, but for the
                            # 'rate_rlw' output field (rate at which each
                            # alternative's ratio-log-likelihood-weight (rlw) bootstrap
                            # draws exceed this threshold) -- kept as its own separate
                            # argument/threshold rather than reusing hypo_rate_threshold,
                            # since the GORIC(A)-weight and log-likelihood-weight ratios
                            # need not be evaluated at the same cutoff. Deliberately not
                            # named/called a "hypothesis rate": unlike rgw, rlw isn't the
                            # quantity GORIC(A) hypothesis selection is actually based on,
                            # so that label would be misleading here. Stored on the
                            # returned object as x$threshold_rlw.
                            threshold_rlw = 1, ...) {

  # iter = NULL (the default): start at iter_min (500) draws and grow by
  # iter_step (100) at a time, up to iter_max (2000), stopping as soon as the
  # "Observed" population's percentile has been stable for two consecutive
  # rounds -- see run_benchmark_simulation(). A user-supplied numeric
  # 'iter' is used as-is (single fixed-size run, as before).
  user_iter <- iter

  if (length(control) == 1) {
    control <- object$objectList[[1]]$control
  }

  mix_weights <- attr(object$objectList[[1]]$wt.bar, "method")
  penalty_factor <- object$penalty_factor

  # Check:
  if (!inherits(object, "con_goric")) {
    stop(paste("\nrestriktor ERROR:",
               "The object should be of class 'con_goric' (a goric(a) object from",
               "the goric() function). However, it belongs to the following class(es):",
               paste(class(object), collapse = ", ")
    ), call. = FALSE)
  }

  validate_iter_args(iter, iter_min = iter_min, iter_step = iter_step,
                     iter_max = iter_max, iter_stability_tol = iter_stability_tol,
                     iter_adequacy_band = iter_adequacy_band)

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

  # Unrestricted (adjusted) group_means
  group_means <- object$b.unrestr
  n_coef <- length(group_means)

  # original model fit (if exists)
  form_model_org <- formula(object$model.org)

  # Which coefficients are group means? 
  # NOTE (ANCOVA, e.g. lm(y ~ -1 + group + x)): only the coefficients that
  # belong to factor (or character/logical) terms of the original model fit
  # are treated as group means. The effect size (Cohen's f), the group sizes
  # and the scaling of the population means all refer to these group means
  # only; the coefficients of continuous covariates are nuisance parameters
  # here: they are kept fixed at their observed estimates in every
  # population (also under 'No-effect') and their (co)variances are rescaled
  # with the overall sample size (sum(N)/sum(alt_group_size)) when
  # alt_group_size is specified. All coefficients (group means and
  # covariates) are still simulated jointly with the full VCOV, as before.
  # Without a fitted model (input is est + VCOV) the type of each estimate is
  # unknown, and -- as before -- all estimates are treated as group means.
  fitLM <- object$model.org
  group_idx <- seq_len(n_coef)
  N_lm <- NULL
  if (!is.null(fitLM) && inherits(fitLM, "lm")) {
    group_info <- tryCatch({
      tt <- terms(fitLM)
      mf <- model.frame(fitLM)
      is_group_var <- function(v) is.factor(v) || is.character(v) || is.logical(v)
      fac <- attr(tt, "factors")
      vars <- rownames(fac)
      if (attr(tt, "response") > 0) vars <- vars[-attr(tt, "response")]
      # factor variables in the model frame (response and "(weights)" etc. excluded)
      group_vars <- intersect(rownames(fac)[-attr(tt, "response")], names(mf))
      group_vars <- group_vars[vapply(mf[group_vars], is_group_var, logical(1))]
      # terms that consist of factor variables only
      group_terms <- which(apply(fac[, , drop = FALSE] > 0, 2, function(v) {
        all(rownames(fac)[v] %in% group_vars)
      }))
      assign <- attr(model.matrix(fitLM), "assign")
      idx <- which(assign %in% group_terms)
      # group sizes: (interaction-)cell counts of the factor variables only.
      # table() returns a 1-D array with a 'dim' attribute; c() strips it to
      # a plain named vector (otherwise VCOV * (N - 1) in compute_cohens_f()
      # fails with "non-conformable arrays").
      counts <- if (length(group_vars) > 0) c(do.call(table, mf[group_vars])) else NULL
      list(idx = idx, counts = counts)
    }, error = function(e) NULL)
    if (!is.null(group_info) && length(group_info$idx) > 0) {
      group_idx <- group_info$idx
      # only usable if there is one count per group mean (e.g., not for an
      # additive model with multiple factors).
      if (length(group_info$counts) == length(group_idx)) {
        N_lm <- group_info$counts
      }
    }
  }
  # number of groups (covariates not included)
  ngroups <- length(group_idx)
  covariate_idx <- setdiff(seq_len(n_coef), group_idx)

  # # Number of subjects per group
  # NOTE: This is needed to rescale vcov based on alt_group_size.
  #       and also for calculating Cohens f. 
  if (!is.null(group_size)) { # So, user specified it as input
    if (length(group_size) == 1) {
      N <- rep(group_size, ngroups) 
    } else { # if vector and not scalar
      if (length(group_size) == ngroups) {
        N <- group_size
      } else { # so, length incorrect
        stop("\nrestriktor ERROR: The argument 'group_size' should be of length 1 or of length ", 
             ngroups, ". It is currently of length ", length(group_size), 
             ". It should be a scalar if all groups have the same size (e.g., group_size = 100)",
             " or a vector with for each group its group size (e.g., group_size = c(75, 100, 120)).", 
             call. = FALSE)
      }
    }
    # Check whether same as obtained from lm object.
    if (!is.null(N_lm) && (length(N) != length(N_lm) || any(N != N_lm))) {
      message("\nrestriktor Message: The argument 'group_size' differs from the group sizes ",
              "retrieved from the lm object. The function proceeded with the user-specified ",
              "'group_size'.\n",
              "Notably, based on group_size, N = ", paste(N, collapse = ", "), 
              "; and, based on the lm object, N = ", paste(N_lm, collapse = ", "), ".")
    }
  } else if (!is.null(N_lm)) { # so, not user specified
    N <- N_lm
  } else {
    stop("\nrestriktor ERROR: The group sizes could not be retrieved from the goric object ",
         "(e.g., because only estimates and their covariance matrix were used as input). ",
         "Please specify them by the argument 'group_size'; e.g., group_size = 100 ",
         "or group_size = c(75, 100, 120).", call. = FALSE)
  }
  N <- as.vector(N)
  names(N) <- names(group_means)[group_idx]
  
  VCOV <- VCOV_orig <- object$VCOV # Is already based on N (so, not N-k)
  
  # If alt_group_size specified, adjust VCOV accordingly
  # Notably, VCOV is based on N not N-k (i.e., sum(N) - ngroups)
  if (!is.null(alt_group_size)) {
  
    if (length(alt_group_size) != 1 && length(alt_group_size) != ngroups) {
      stop("\nrestriktor ERROR: The argument 'alt_group_size' should be of length 1 or ",
           ngroups, " (or NULL) but not of length ", length(alt_group_size), ".", 
           call. = FALSE)
    }
    alt_N <- rep_len(alt_group_size, ngroups)
    # scaling factor per coefficient: N/alt_N for the group means, and the
    # ratio of the total sample sizes for covariates (if any). The VCOV is
    # scaled symmetrically (for a diagonal VCOV, as in an ANOVA model, this
    # equals scaling each variance by N/alt_N).
    scale_coef <- rep(sum(N) / sum(alt_N), n_coef)
    scale_coef[group_idx] <- N / alt_N
    VCOV <- VCOV_orig * sqrt(outer(scale_coef, scale_coef))
    N <- alt_N
    names(N) <- names(group_means)[group_idx]
    #
    # The sample gorica(c) value must also be adjusted, 
    # thus we need to fit a new goric-object
    # with est, VCOV = new VCOV, and (if goricac) N = new N (so alt_group_size).
  }
  
  
  # (re-)calculate gorica(c) if:
  # - goric(c) not yet gorica(c)
  # - if adjusted sample size
  # Will do it anyway
  #
  # If needed, adjust type from goric(c) to gorica(c)
  type <- switch(object$type,
                 "goric" = {
                   message("\nrestriktor Message: 'goric' has been converted to 'gorica'.")
                   "gorica"
                 },
                 "goricc" = {
                   message("\nrestriktor Message: 'goricc' has been converted to 'goricac'.")
                   "goricac"
                 },
                 object$type)
  #
  # (re-)calculate gorica(c)
  object <-
    goric(
      group_means,
      VCOV = VCOV, # based on N not N-k
      sample_nobs = sum(N), # Needed for type = "goricac" - daar genoeg aan sum(N) of moet het juist gehele N hebben?
      hypotheses = hypos,
      comparison = object$comparison,
      type = type,
      control = control,
      mix_weights = mix_weights,
      penalty_factor = penalty_factor,
      # same prior weights as in the goric object, so that the preferred
      # hypothesis (and the weights) match those of the goric object
      priorICweights = object$priorICweights,
      #Heq = FALSE,
      Heq = Heq,
      ...
    )
  
  

  ## Compute observed Cohens f
  # (based on the group means only; see the note on covariates above)
  cohens_f_observed <- compute_cohens_f(group_means[group_idx], N, 
                                        VCOV[group_idx, group_idx, drop = FALSE])

  # effect size population
  if (is.null(pop_es)) {
    pop_es <- c(0, cohens_f_observed)
    names(pop_es) <- c("No-effect", "Observed")
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

  #means_pop_all <- compute_population_means(pop_es, ratio_pop_means, var_e, ngroups)
  # Only the group means are scaled to the targeted effect size; covariate
  # coefficients (if any) are kept at their observed estimates.
  means_pop_all <- t(sapply(pop_es, function(x) {
    means_pop <- group_means
    means_pop[group_idx] <- generate_scaled_means(group_means[group_idx], target_f = x, N,
                                                  VCOV[group_idx, group_idx, drop = FALSE])
    means_pop
  }))
  colnames(means_pop_all) <- colnames(coef(object))
  rownames(means_pop_all) <- paste0("pop_es = ", pop_es)

  # preferred hypothesis
  pref_hypo <- which.max(object$result[, 7])
  pref_hypo_name <- object$result$model[pref_hypo]

  if (is.null(quant)) {
    quant <- c(.05, .35, .50, .65, .95)
    names_quant <- c("Sample", "5%", "35%", "50%", "65%", "95%")
  } else {
    names_quant <- c("Sample", paste0(as.character(quant*100), "%"))
  }

  # TO DO dus niet ook goricac (evt type = type?); als default en het dus wel
  # kan, dan ook 'sample_nobs' nodig toch...
  sim <- run_benchmark_simulation(
    nr_es = nr_es, rnames = rnames, name_prefix = "pop_es = ",
    center_matrix = means_pop_all, colnames_vec = names(group_means),
    VCOV = VCOV, hypos = hypos, pref_hypo = pref_hypo,
    comparison = object$comparison, control = control,
    mix_weights = mix_weights, penalty_factor = penalty_factor, Heq = Heq,
    object = object, iter = user_iter,
    es_labels = paste0(es, " (", names(es), ")"),
    band = iter_adequacy_band,
    stability_tol = iter_stability_tol,
    iter_min = iter_min, iter_step = iter_step, iter_max = iter_max,
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
                                            est = group_means,
                                            VCOV, control, ...)

  # residual error variance
  # var_e <- as.vector(diag(VCOV)[1])
  # var_e_data <- as.vector(diag(VCOV_orig)[1])
  # This is used as reference value, because the diag. elements will differ if group sizes differ
  
  OUT <- list(
    type = object$type,
    comparison = object$comparison,
    ngroups = ngroups,
    group_size = N,
    group_means_observed = group_means,
    #ratio_group_means.data = ratio_data,
    cohens_f_observed = cohens_f_observed,
    #res_var_observed = var_e_data,
    pop_es = pop_es, 
    pop_group_means = means_pop_all,
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
                            threshold_rlw = 1, ...) {

  # iter = NULL (the default): start at iter_min (500) draws and grow by
  # iter_step (100) at a time, up to iter_max (2000), stopping as soon as the
  # "Observed" population's percentile has been stable for two consecutive
  # rounds -- see run_benchmark_simulation(). A user-supplied numeric
  # 'iter' is used as-is (single fixed-size run, as before).
  user_iter <- iter

  # TO DO als in means nog 'group_size' argument dan ms hier zeggen dat dat sample_size moet zijn?
  #       NB alleen nodig als alt_sample_size nodig is...
    
  if (length(control) == 1) {
    control <- object$objectList[[1]]$control
  }
  
  mix_weights <- attr(object$objectList[[1]]$wt.bar, "method")
  penalty_factor <- object$penalty_factor
  
  # Check if object is of class con_goric
  if (!inherits(object, "con_goric")) {
    stop(paste("\nrestriktor ERROR:", 
               "The object should be of class 'con_goric' (a goric(a) object from the goric() function).",
               "However, it belongs to the following class(es):", 
               paste(class(object), collapse = ", ")), call. = FALSE)
  }
 
  validate_iter_args(iter, iter_min = iter_min, iter_step = iter_step,
                     iter_max = iter_max, iter_stability_tol = iter_stability_tol,
                     iter_adequacy_band = iter_adequacy_band)

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
  type <- object$type
  est_sample <- object$b.unrestr
  n_coef <- length(est_sample)
  
  if (is.null(pop_est)) {
    #NE <- theta_restricted(theta = est_sample, V = VCOV, R = object$constraints[[1]], rhs = object$rhs[[1]])
    # TO determine pop_est, we need the constraint matrix for the preferred hypothesis.
    # In that constraint matrix / hypothesis, the inequalities will be set to equalities.
    # First determine the preferred hypothesis:
    pref_hypo <- which.max(object$result[, 7])
    if (pref_hypo > length(object$constraints) || 
        is.null(object$constraints[[pref_hypo]])) {
      stop("\nrestriktor ERROR: The preferred hypothesis is the ", object$result$model[pref_hypo], 
           " hypothesis, for which no constraint matrix is available. Hence, the default 'No-effect' population ",
           "estimates cannot be determined. Please specify the population estimates via the ",
           "argument 'pop_est'.", call. = FALSE)
    }
    NE <- theta_restricted(theta = est_sample, V = VCOV, R = object$constraints[[pref_hypo]], rhs = object$rhs[[pref_hypo]])
    # Note that VCOV is the unbiased cov.mx estimate.
    pop_est <- matrix(rbind(NE, est_sample), nrow = 2)
    row.names(pop_est) <- c("No-effect", "Observed")
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
                   message("\nrestriktor Message: 'goric' has been converted to 'gorica'.")
                   "gorica"
                 },
                 "goricc" = {
                   message("\nrestriktor Message: 'goricc' has been converted to 'goricac'.")
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
    VCOV = VCOV,
    sample_nobs = N[1], # Needed for type = "goricac"
    # TO DO, dit is toch maar een getal en als niet dan ws de sum (iig bij anova wel)
    #sample_nobs = sum(N), # Needed for type = "goricac" - daar genoeg aan sum(N) of moet het juist gehele N hebben?
    hypotheses = hypos,
    comparison = comparison,
    type = type,
    control = control,
    mix_weights = mix_weights,
    penalty_factor = penalty_factor,
    # same prior weights as in the goric object, so that the preferred
    # hypothesis (and the weights) match those of the goric object
    priorICweights = object$priorICweights,
    Heq = Heq,
    ...
  )
  
  
  
  if (is.null(quant)) {
    quant <- c(.05, .35, .50, .65, .95)
    names_quant <- c("Sample", "5%", "35%", "50%", "65%", "95%")
  } else {
    names_quant <- c("Sample", paste0(as.character(quant*100), "%"))
  }
  
  pref_hypo <- which.max(object$result[, 7])
  pref_hypo_name <- object$result$model[pref_hypo]
  
  nr_es  <- nrow(pop_est)

  sim <- run_benchmark_simulation(
    nr_es = nr_es, rnames = rnames, name_prefix = "pop_est = ",
    center_matrix = pop_est, colnames_vec = names(est_sample),
    VCOV = VCOV, hypos = hypos, pref_hypo = pref_hypo,
    comparison = comparison, control = control,
    mix_weights = mix_weights, penalty_factor = penalty_factor, Heq = Heq,
    object = object, iter = user_iter,
    band = iter_adequacy_band,
    stability_tol = iter_stability_tol,
    iter_min = iter_min, iter_step = iter_step, iter_max = iter_max,
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
    #
    median_bias_check_gw = bias_check$median_bias_check_gw,
    median_bias_check_lw = bias_check$median_bias_check_lw,
    #
    iter = iter
  )

  class(OUT) <- c("benchmark_asymp", "benchmark", "list")
  return(OUT)
}
