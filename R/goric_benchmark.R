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

  # Check:
  if (!inherits(object, "con_goric")) {
    stop(paste("\nrestriktor ERROR:",
               "The object should be of class 'con_goric' (a goric(a) object from",
               "the goric() function). However, it belongs to the following class(es):",
               paste(class(object), collapse = ", ")
    ), call. = FALSE)
  }
  check_benchmark_weights(object)

  # no user-specified control options: use those of the goric object
  # (previously 'length(control) == 1', which silently discarded a control
  # list with exactly one option)
  if (length(control) == 0) {
    control <- object$objectList[[1]]$control
  }

  mix_weights <- attr(object$objectList[[1]]$wt.bar, "method")
  penalty_factor <- object$penalty_factor

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
  if (is.null(hypos)) {
    stop("\nrestriktor ERROR: benchmark() requires hypotheses specified as text ",
         "(e.g., 'x1 > x2'); the GORIC(A) object was fitted with hypotheses given ",
         "as constraint matrices, which are not (yet) supported by benchmark().",
         call. = FALSE)
  }
  nr_hypos <- dim(object$result)[1]
  Heq <- object$Heq

  # Unrestricted (adjusted) group_means
  group_means <- object$b.unrestr
  n_coef <- length(group_means)

  # original model fit (if exists)
  form_model_org <- formula(object$model.org)

  # Which coefficients are group means? 
  # NOTE (covariates, e.g. lm(y ~ -1 + group + x) or lm(y ~ -1 + group + sex)):
  # only the coefficients of the factor (or character/logical) term that the
  # hypotheses refer to are treated as group means (the first factor term if
  # the hypotheses refer to no coefficient at all): without an intercept, R
  # codes the first factor term as cell means (one indicator per group),
  # while any further factor term is coded as contrasts (differences with its
  # reference level), which are not means and have no group size; hypotheses
  # on such contrasts, on covariates or on coefficients of several terms give
  # an error (use model_type = 'asymp' then). The effect size
  # (Cohen's f), the group sizes and the scaling of the population means all
  # refer to these group means only; all other coefficients (continuous
  # covariates and the contrasts of further factors) are nuisance parameters
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
  group_term_label <- NULL
  group_size_problem <- NULL
  if (!is.null(fitLM) && inherits(fitLM, "lm")) {
    # A model with an intercept (e.g. y ~ group) has an intercept and
    # contrasts as coefficients, not group means (benchmark() already rejects
    # this; benchmark_means() can also be called directly).
    if (isTRUE(attr(terms(fitLM), "intercept") == 1)) {
      stop("\nrestriktor ERROR: The model fit underlying the GORIC(A) object has an ",
           "intercept (e.g., lm(y ~ group)), so its coefficients are an intercept and ",
           "contrasts (differences with the reference group), not group means. ",
           "model_type = 'means' requires the group means as coefficients: please refit ",
           "the model without an intercept (e.g., lm(y ~ -1 + group)), formulate the ",
           "hypotheses in terms of the group means, and re-run goric().", call. = FALSE)
    }
    group_info <- tryCatch({
      tt <- terms(fitLM)
      mf <- model.frame(fitLM)
      is_group_var <- function(v) is.factor(v) || is.character(v) || is.logical(v)
      fac <- attr(tt, "factors")
      vars <- rownames(fac)
      if (attr(tt, "response") > 0) vars <- vars[-attr(tt, "response")]
      # factor variables in the model frame (response and "(weights)" etc. excluded)
      group_vars <- intersect(vars, names(mf))
      group_vars <- group_vars[vapply(mf[group_vars], is_group_var, logical(1))]
      # terms that consist of factor variables only; the first of these is
      # the (cell-means coded) grouping term
      factor_terms <- which(apply(fac[, , drop = FALSE] > 0, 2, function(v) {
        any(v) && all(rownames(fac)[v] %in% group_vars)
      }))
      if (length(factor_terms) == 0) stop("no factor term")
      assign <- attr(model.matrix(fitLM), "assign")
      # The grouping term is the factor term the hypotheses are about: the
      # one whose coefficients contain ALL coefficients referred to in the
      # hypotheses (the columns of the constraint matrices with a non-zero
      # entry; ':=' definitions are already resolved there). Previously the
      # first factor term was taken, so that for lm(y ~ -1 + sex + group)
      # with hypotheses on 'group' the sex means were benchmarked while the
      # hypothesised group coefficients were kept fixed in every population.
      ref_idx <- integer(0)
      if (is.list(object$constraints) && length(object$constraints) > 0) {
        Amat <- do.call(rbind, lapply(object$constraints, function(A) {
          A <- as.matrix(A)
          if (ncol(A) == n_coef) A else NULL
        }))
        if (!is.null(Amat)) ref_idx <- which(colSums(abs(Amat)) > 0)
      }
      hyp_problem <- NULL
      if (length(ref_idx) > 0) {
        candidates <- factor_terms[vapply(factor_terms, function(tm) {
          all(ref_idx %in% which(assign == tm))
        }, logical(1))]
        group_term <- if (length(candidates) > 0) candidates[1] else NA
        if (is.na(group_term)) {
          hyp_problem <- paste0(
            "The hypotheses refer to the coefficient(s) ",
            paste(names(group_means)[ref_idx], collapse = ", "),
            ", which do not belong to a single factor (grouping) term of the fitted ",
            "model (they include covariates or coefficients of different terms).")
        } else if (length(which(assign == group_term)) !=
                   length(c(do.call(table, mf[rownames(fac)[fac[, group_term] > 0]])))) {
          hyp_problem <- paste0(
            "The hypotheses refer to the coefficient(s) ",
            paste(names(group_means)[ref_idx], collapse = ", "),
            " of the factor term '", colnames(fac)[group_term], "', but this term is ",
            "coded as contrasts (differences with its reference level), not as group ",
            "means, since it is not the first factor term of the model (e.g., ",
            "lm(y ~ -1 + sex + group)); refit the model with this term first (e.g., ",
            "lm(y ~ -1 + group + sex)) and formulate the hypotheses in terms of its ",
            "group means.")
        }
      } else {
        group_term <- factor_terms[1]
      }
      if (!is.null(hyp_problem)) {
        # (no return(): that would return from benchmark_means() itself)
        list(hyp_problem = hyp_problem)
      } else {
        term_vars <- rownames(fac)[fac[, group_term] > 0]
        idx <- which(assign == group_term)
        # group sizes: (interaction-)cell counts of the factor variable(s) of
        # the grouping term. table() returns a 1-D array with a 'dim'
        # attribute; c() strips it to a plain named vector (otherwise
        # arithmetic with N in compute_cohens_f() fails with "non-conformable
        # arrays").
        counts <- c(do.call(table, mf[term_vars]))
        list(idx = idx, counts = counts, label = colnames(fac)[group_term])
      }
    }, error = function(e) NULL)
    if (!is.null(group_info$hyp_problem)) {
      stop("\nrestriktor ERROR: benchmark_means() requires hypotheses on the group means ",
           "of a single factor term (e.g., lm(y ~ -1 + group) with hypotheses on ",
           "group1, group2, ...). ", group_info$hyp_problem, " Alternatively, use ",
           "model_type = 'asymp' (benchmark_asymp()), which benchmarks the estimates ",
           "themselves.", call. = FALSE)
    }
    if (!is.null(group_info) && length(group_info$idx) > 0) {
      group_idx <- group_info$idx
      group_term_label <- group_info$label
      # only usable if there is one count per group mean (i.e., the term is
      # coded as cell means).
      if (length(group_info$counts) == length(group_idx)) {
        N_lm <- group_info$counts
      } else {
        group_size_problem <- paste0(
          "The group sizes could not be derived from the fitted model: the ",
          length(group_idx), " coefficient(s) of the term '", group_term_label,
          "' (", paste(names(group_means)[group_idx], collapse = ", "),
          ") do not correspond to the ", length(group_info$counts),
          " cell(s) of that term, so they are apparently not cell means.")
      }
    } else {
      group_size_problem <- paste0(
        "The group sizes could not be derived from the fitted model, since it has ",
        "no factor (grouping) variable, so it is unclear which coefficients are ",
        "group means. All coefficients are treated as group means.")
    }
  } else {
    group_size_problem <- paste0(
      "The group sizes could not be retrieved from the goric object, since only ",
      "estimates and their covariance matrix were used as input (all estimates are ",
      "treated as group means).")
  }
  # number of groups (covariates not included)
  ngroups <- length(group_idx)
  covariate_idx <- setdiff(seq_len(n_coef), group_idx)
  if (length(covariate_idx) > 0 && !is.null(group_term_label)) {
    message("\nrestriktor Message: The coefficients of the term '", group_term_label,
            "' (", paste(names(group_means)[group_idx], collapse = ", "),
            ") are treated as the group means. The other coefficient(s) (",
            paste(names(group_means)[covariate_idx], collapse = ", "),
            ") are treated as covariates: they are not part of Cohen's f and are kept ",
            "at their observed estimates in every population (also under 'No-effect').")
  }

  # Pattern of the population means (see generate_scaled_means()): by
  # default the observed group means; otherwise the user-specified
  # ratio_pop_means (one value per group).
  if (!is.null(ratio_pop_means)) {
    if (!is.numeric(ratio_pop_means) || length(ratio_pop_means) != ngroups ||
        anyNA(ratio_pop_means) || any(!is.finite(ratio_pop_means))) {
      stop("\nrestriktor ERROR: The argument 'ratio_pop_means' should be a numeric vector ",
           "of length ", ngroups, " (one value per group), e.g., ratio_pop_means = c(",
           paste(seq_len(ngroups), collapse = ", "), "). It is currently of length ",
           length(ratio_pop_means), ".", call. = FALSE)
    }
    ratio_pop_means <- as.vector(ratio_pop_means)
    names(ratio_pop_means) <- names(group_means)[group_idx]
  }

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
    stop("\nrestriktor ERROR: ", group_size_problem,
         " Please specify the group sizes by the argument 'group_size'; e.g., ",
         "group_size = 100 or group_size = c(75, 100, 120).", call. = FALSE)
  }
  N <- as.vector(N)
  names(N) <- names(group_means)[group_idx]
  
  VCOV <- VCOV_orig <- object$VCOV # Is already based on N (so, not N-k)

  # Residual (within-group) error variance sigma2, used for Cohen's f (the
  # observed f and the scaling of the population means to 'pop_es', see
  # compute_cohens_f()/generate_scaled_means()). It is the error variance
  # the simulation is based on: the draws are generated from VCOV, which for
  # a fitted lm model is vcov(fit) * (N - k) / N, i.e. based on the ML
  # estimate sigma2 = RSS / N (see VCOV.unbiased()), so that Var(mean_g) =
  # sigma2 / N_g. With a fitted lm model, sigma2 is therefore RSS / N (for a
  # one-way design this equals N_g * VCOV[g, g] for every group; for an
  # ANCOVA it is the residual variance after adjusting for the covariates,
  # which is what Cohen's f refers to); the observed f is then the usual
  # plug-in value sqrt(SS_between / SS_within). (Previously sigma(fit)^2 =
  # RSS / (N - k) was used here, so that the population means for a given
  # pop_es had an effective f of pop_es * sqrt(N / (N - k)) relative to the
  # draws' own error variance.) Without a model (input is est + VCOV) it is
  # derived from VCOV as the average of N_g * VCOV[g, g], assuming
  # independent group means with a common error variance -- see
  # residual_variance_from_vcov(), which warns if that assumption seems
  # violated. Computed from the ORIGINAL group sizes and VCOV: the product
  # N_g * VCOV[g, g] is unaffected by the alt_group_size rescaling below
  # (which scales VCOV[g, g] by N_g / alt_N_g, i.e. keeps Var(mean_g) =
  # sigma2 / N_g for the new group sizes), so the same sigma2 applies to the
  # alternative group sizes -- and to the simulation, which draws from VCOV.
  # N is the sample size of the FITTED MODEL (nobs(fitLM)): the VCOV of the
  # object is vcov(fit) * (N - k) / N with that N (see VCOV.unbiased()), also
  # when a different 'sample_nobs' was given to goric() -- which goric()
  # accepts with a message, and stores as object$sample_nobs. Using
  # object$sample_nobs here (as was done briefly) gave sigma2 = RSS /
  # sample_nobs, i.e. an f and population means that do not belong to the
  # draws.
  sigma2 <- NULL
  if (!is.null(fitLM) && inherits(fitLM, "lm") && !inherits(fitLM, "glm") &&
      !isTRUE(object$objectList[[1]]$missing == "fiml")) {
    # A weighted fit (lm(..., weights = w)) is not supported: its VCOV gives
    # Var(mean_g) = sigma2_w / sum(w_g), not sigma2 / N_g, so neither the
    # group sizes nor a single residual variance describe the draws.
    w_lm <- tryCatch(stats::weights(fitLM), error = function(e) NULL)
    if (!is.null(w_lm) && any(w_lm != 1)) {
      stop("\nrestriktor ERROR: The model fit underlying the GORIC(A) object is a ",
           "weighted fit (lm(..., weights = )). The covariance matrix of weighted ",
           "estimates does not correspond to the group sizes and a single residual ",
           "error variance, which benchmark_means() needs for Cohen's f and the ",
           "population means. Please use model_type = 'asymp' (benchmark_asymp()) ",
           "instead, or refit the model without weights.", call. = FALSE)
    }
    N_model <- tryCatch(stats::nobs(fitLM), error = function(e) NULL)
    # a vector of group sizes is summed (as goric() does)
    if (is.numeric(object$sample_nobs) && length(object$sample_nobs) >= 1 &&
        is.numeric(N_model) && length(N_model) == 1 && sum(object$sample_nobs) != N_model) {
      message("\nrestriktor Message: The 'sample_nobs' of the GORIC(A) object (",
              sum(object$sample_nobs), ") differs from the sample size of the fitted model (",
              N_model, "). The benchmark uses the model-based sample size, on which ",
              "the covariance matrix of the estimates is based.")
    }
    sigma2 <- tryCatch(stats::deviance(fitLM) / N_model, # RSS / N
                       error = function(e) NULL)
    if (!is.numeric(sigma2) || length(sigma2) != 1 || !is.finite(sigma2)) {
      sigma2 <- NULL
    }
  }
  if (is.null(sigma2)) {
    sigma2 <- residual_variance_from_vcov(N, VCOV_orig[group_idx, group_idx, drop = FALSE])
  }

  ## Compute observed Cohens f
  # (based on the group means only -- see the note on covariates above -- and
  # on the group sizes of the data, also when alt_group_size is specified)
  cohens_f_observed <- compute_cohens_f(group_means[group_idx], N, sigma2)
  cohens_f_alt_group_size <- NULL

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
    # Cohen's f of the observed group means with the alternative group sizes
    # (as weights); differs from cohens_f_observed only for unequal
    # alternative group sizes. Printed separately by print.benchmark().
    cohens_f_alt_group_size <- compute_cohens_f(group_means[group_idx], N, sigma2)
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
  
  check_benchmark_weights(object, refit = TRUE)

  # effect size population
  default_pop_es <- is.null(pop_es)
  if (default_pop_es) {
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
  rnames <- unique_population_names(rnames, "pop_es")
  names(pop_es) <- rnames

  es <- pop_es
  nr_es <- length(es)

  # Population means per pop_es: the pattern of the means (the observed group
  # means by default, or 'ratio_pop_means' if specified) scaled such that
  # Cohen's f equals pop_es, see generate_scaled_means(). Only the group means
  # are scaled to the targeted effect size; covariate coefficients (if any)
  # are kept at their observed estimates. The default 'Observed' population
  # (pop_es = NULL, no ratio_pop_means) uses the observed estimates as-is.
  pattern_means <- if (is.null(ratio_pop_means)) group_means[group_idx] else ratio_pop_means
  use_observed <- is.null(ratio_pop_means) & names(pop_es) == "Observed" & default_pop_es
  means_pop_all <- t(sapply(seq_along(pop_es), function(i) {
    means_pop <- group_means
    if (use_observed[i]) {
      means_pop[group_idx] <- group_means[group_idx]
    } else {
      means_pop[group_idx] <- generate_scaled_means(pattern_means, target_f = pop_es[i], N,
                                                    sigma2)
    }
    means_pop
  }))
  colnames(means_pop_all) <- names(group_means)
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

  # The draws are evaluated with the same criterion (gorica/goricac) and
  # sample size as the (refitted) object above, so that the 'Sample' value
  # and the benchmark distribution are based on the same criterion.
  sim <- run_benchmark_simulation(
    nr_es = nr_es, rnames = rnames, name_prefix = "pop_es = ",
    center_matrix = means_pop_all, colnames_vec = names(group_means),
    VCOV = VCOV, hypos = hypos, pref_hypo = pref_hypo,
    comparison = object$comparison, control = control,
    mix_weights = mix_weights, penalty_factor = penalty_factor, Heq = Heq,
    object = object, iter = user_iter,
    type = type, sample_nobs = sum(N),
    es_labels = paste0(es, " (", names(es), ")"),
    # the 'Observed' population has the observed effect size, but -- with
    # ratio_pop_means -- not the observed means (see B24/print.benchmark())
    observed_label = if (is.null(ratio_pop_means)) {
      "the 'Observed' population"
    } else {
      "the 'Observed' population (observed effect size; means from ratio_pop_means)"
    },
    band = iter_adequacy_band,
    stability_tol = iter_stability_tol,
    iter_min = iter_min, iter_step = iter_step, iter_max = iter_max,
    ...
  )
  parallel_function_results <- sim$parallel_function_results
  iter <- sim$iter # final number of SUCCESSFUL draws, per population

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
    bias_check <- check_iter_adequacy(benchmark_results, "pop_es = Observed", user_iter,
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
    cohens_f_alt_group_size = cohens_f_alt_group_size,
    res_var = sigma2, # residual error variance (RSS / N) used for Cohen's f
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
    # number of successful draws per population (a single number when it is
    # the same for every population); see also iter_requested and the
    # failed/warned draw counts from run_benchmark_simulation().
    iter = iter,
    iter_requested = sim$iter_requested,
    n_failed_draws = sim$n_failed_draws,
    n_warned_draws = sim$n_warned_draws,
    draw_errors = sim$draw_errors,
    draw_warnings = sim$draw_warnings
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
    
  # Check if object is of class con_goric
  if (!inherits(object, "con_goric")) {
    stop(paste("\nrestriktor ERROR:", 
               "The object should be of class 'con_goric' (a goric(a) object from the goric() function).",
               "However, it belongs to the following class(es):", 
               paste(class(object), collapse = ", ")), call. = FALSE)
  }
  check_benchmark_weights(object)

  # no user-specified control options: use those of the goric object
  # (previously 'length(control) == 1', which silently discarded a control
  # list with exactly one option)
  if (length(control) == 0) {
    control <- object$objectList[[1]]$control
  }
  
  mix_weights <- attr(object$objectList[[1]]$wt.bar, "method")
  penalty_factor <- object$penalty_factor
 
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
  if (is.null(hypos)) {
    stop("\nrestriktor ERROR: benchmark() requires hypotheses specified as text ",
         "(e.g., 'x1 > x2'); the GORIC(A) object was fitted with hypotheses given ",
         "as constraint matrices, which are not (yet) supported by benchmark().",
         call. = FALSE)
  }
  
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
    if (is.data.frame(pop_est)) {
      pop_est <- as.matrix(pop_est)
    }
    if (is.vector(pop_est) && is.null(dim(pop_est))) {
      pop_est <- matrix(pop_est, nrow = 1, ncol = length(pop_est),
                        dimnames = list(NULL, names(pop_est)))
    }
    if (!is.matrix(pop_est) || !is.numeric(pop_est) || anyNA(pop_est)) {
      stop("\nrestriktor ERROR: The argument 'pop_est' should be a numeric vector (one ",
           "population) or a numeric matrix / data.frame (one row per population) with ",
           "one column per estimate, without missing values.", call. = FALSE)
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
  rnames <- unique_population_names(rnames, "pop_est")
  row.names(pop_est) <- rnames
  
  colnames(pop_est) <- names(est_sample)
  N <- object$sample_nobs #length(object$model.org$residuals)
  # For an object based on a fitted lm model, the sample size is that of the
  # model (nobs(fit)): its VCOV is based on that N (see VCOV.unbiased()),
  # also when a different 'sample_nobs' was given to goric() (accepted with
  # a message and stored as object$sample_nobs), so the alt_sample_size
  # rescaling and the goricac should use the model's N as well.
  # (Not for a FIML fit: its VCOV and N are the FIML ones, not those of the lm.)
  if (inherits(object$model.org, "lm") && !inherits(object$model.org, "mlm") &&
      !isTRUE(object$objectList[[1]]$missing == "fiml")) {
    N_model <- tryCatch(stats::nobs(object$model.org), error = function(e) NULL)
    if (is.numeric(N_model) && length(N_model) == 1 && is.finite(N_model) && N_model > 0) {
      if (is.numeric(N) && length(N) >= 1 && sum(N) != N_model) {
        message("\nrestriktor Message: The 'sample_nobs' of the GORIC(A) object (",
                sum(N), ") differs from the sample size of the fitted model (", N_model,
                "). The benchmark uses the model-based sample size, on which the ",
                "covariance matrix of the estimates is based.")
      }
      N <- N_model
    }
  }
  
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
  
  
  # Original sample size: from the goric object (a fitted model, or
  # sample_nobs given to goric()), else from 'sample_size'. A user-specified
  # 'sample_size' that differs from the object's is used, with a message
  # (previously it was silently ignored whenever the object had one). A
  # vector of group sizes is accepted and summed (as goric() does).
  if (!is.null(sample_size)) {
    if (!is.numeric(sample_size) || length(sample_size) < 1 || anyNA(sample_size) ||
        any(!is.finite(sample_size)) || any(sample_size <= 0)) {
      stop("\nrestriktor ERROR: The argument 'sample_size' should be a positive number ",
           "(the total sample size) or a vector of positive group sizes.", call. = FALSE)
    }
    if (!(is.null(N) || all(N == 0)) && sum(N) != sum(sample_size)) {
      message("\nrestriktor Message: The argument 'sample_size' (total: ", sum(sample_size),
              ") differs from the sample size of the goric object (", sum(N), "). The ",
              "function proceeded with the user-specified 'sample_size'.")
    }
    N <- sample_size
  }
  # Alternative sample size: the covariance matrix of the estimates scales
  # with the ratio of the TOTAL sample sizes (a single factor, so the result
  # is symmetric like VCOV itself; the design proportions are kept). A vector
  # of alternative group sizes is summed. (Previously 'VCOV * N / alt' was
  # computed element-wise, which for a vector 'sample_size' recycled N over
  # the matrix: wrong scaling and a non-symmetric result.)
  if (!is.null(alt_sample_size)) {
    # Controleer of de originele steekproefgrootte beschikbaar is
    if (is.null(N) || all(N == 0)) {
      stop("\nrestriktor ERROR: Please provide the original sample size(s) using the argument `sample_size`.", call. = FALSE)
    }
    if (!is.numeric(alt_sample_size) || length(alt_sample_size) < 1 || anyNA(alt_sample_size) ||
        any(!is.finite(alt_sample_size)) || any(alt_sample_size <= 0)) {
      stop("\nrestriktor ERROR: The argument 'alt_sample_size' should be a positive number ",
           "(the alternative total sample size) or a vector of positive group sizes.",
           call. = FALSE)
    }
    VCOV <- VCOV * (sum(N) / sum(alt_sample_size))
    if (!isSymmetric(unname(VCOV))) {
      stop("\nrestriktor ERROR: The covariance matrix rescaled to 'alt_sample_size' is not ",
           "symmetric, so the covariance matrix of the GORIC(A) object is apparently not ",
           "symmetric.", call. = FALSE)
    }
    N <- alt_sample_size
  }
  # The goricac requires the sample size (also for every benchmark draw)
  if (type == "goricac" && (is.null(N) || all(N == 0))) {
    stop("\nrestriktor ERROR: The GORIC(A) object is of type '", object$type, "', so the ",
         "benchmark is based on the goricac, which requires the sample size. However, no ",
         "sample size is available from the GORIC(A) object. Please specify it via the ",
         "argument 'sample_size'.", call. = FALSE)
  }
  # overall sample size (a vector of group sizes is summed, as in goric())
  sample_nobs <- if (is.null(N)) NULL else sum(N)

  # Herbereken met goric-functie
  object <- goric(
    est_sample,
    VCOV = VCOV,
    sample_nobs = sample_nobs, # Needed for type = "goricac"
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
  check_benchmark_weights(object, refit = TRUE)
  
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
    # same criterion (gorica/goricac) and sample size as the refitted object
    type = type, sample_nobs = sample_nobs,
    band = iter_adequacy_band,
    stability_tol = iter_stability_tol,
    iter_min = iter_min, iter_step = iter_step, iter_max = iter_max,
    ...
  )
  parallel_function_results <- sim$parallel_function_results
  iter <- sim$iter # final number of SUCCESSFUL draws, per population

  benchmark_results <- get_results_benchmark(parallel_function_results, object, pref_hypo,
                                             pref_hypo_name, quant, names_quant, nr_hypos,
                                             hypo_rate_threshold = hypo_rate_threshold,
                                             threshold_rlw = threshold_rlw)

  if (!is.null(user_iter)) {
    # Fixed 'iter': run_benchmark_simulation() does not auto-grow or message
    # in this case, so do the (non-growing) adequacy check here instead.
    bias_check <- check_iter_adequacy(benchmark_results, "pop_est = Observed", user_iter,
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
    # see benchmark_means()
    iter = iter,
    iter_requested = sim$iter_requested,
    n_failed_draws = sim$n_failed_draws,
    n_warned_draws = sim$n_warned_draws,
    draw_errors = sim$draw_errors,
    draw_warnings = sim$draw_warnings
  )

  class(OUT) <- c("benchmark_asymp", "benchmark", "list")
  return(OUT)
}
