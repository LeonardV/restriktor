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
                            quant = NULL, iter = 2000,
                            control = list(),
                            ncpus = 1, seed = NULL, ...) {
# TO DO in asymp functie heet het niet group_size maar sample_size (dito sample_size) 
  #     en dat kan vervelend voor gebruiker zijn...
#       als argument ms toch sample_size en dan direct hierna:
  # group_size <- sample_size
  # alt_sample_size <- alt_group_size
  # Ws 'andersom' iets doen met sample_size <- sum(group_size)
  # TO DO dan wel ook help file en tutorial / example R scripts aanpassen ws!
  # TO DO waar sample_nobs oid nodig is, daar direct de sum meegeven! anders steeds message nl.
  
  # NOTE: group_size is needed to rescale vcov based on alt_group_size.
  #       and also for calculating Cohens f. 

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

  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  # keep current plan
  current_plan <- plan()
  on.exit(plan(current_plan), add = TRUE)
  if (inherits(current_plan, "sequential")) {
    if (.Platform$OS.type == "windows") {
      plan(future::multisession, workers = ncpus)
    } else {
      plan(future::multicore, workers = ncpus)
    }
  }

  # Hypotheses
  hypos <- object$hypotheses_usr
  nr_hypos <- dim(object$result)[1]

  # Unrestricted (adjusted) group_means
  group_means <- object$b.unrestr
  # number of groups
  ngroups <- length(group_means)

  # original model fit (if exists)
  form_model_org <- formula(object$model.org)

  # ES and ratio in data
  #if (is.null(object$model.org)) {
  #  # Number of subjects per group
    if (is.null(group_size)) {
      stop("\nrestriktor ERROR: please specify the group-size, e.g. group_size = 100.", call. = FALSE)
    } else if (length(group_size) == 1) {
      N <- rep(group_size, ngroups)
    } else {
      N <- group_size
    }
  #} else {
  #  # Number of subjects per group
  #  # Note that the next assumes equal group size, which does not need to be the case.
  #  N <- rep(object$sample_nobs/ngroups, ngroups) #colSums(model.matrix(object$model.org))
  #}

  VCOV <- VCOV_orig <- object$VCOV # Is already based on N (so, not N-k)
  
  # If alt_group_size specified, adjust VCOV accordingly
  # Notably, VCOV is based on N not N-k (i.e., sum(N) - ngroups)
  if (!is.null(alt_group_size)) {
  
    if (length(alt_group_size) != 1 && length(alt_group_size) != ngroups) {
      return(paste0("The argument alt_group_size should be of length 1 or ",
                    ngroups, " (or NULL) but not of length ", length(alt_group_size), "."))
    }
    VCOV <- VCOV_orig * N / alt_group_size
    N <- alt_group_size
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
      Heq = FALSE, # TO DO Default keuze of kan deze dus nooit TRUE zijn, ook niet als in het origneel gebruikt?
      ...
    )
  

  ## Compute observed Cohens f
  cohens_f_observed <- compute_cohens_f(group_means, N, VCOV)

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
  means_pop_all <- t(sapply(pop_es, function(x) generate_scaled_means(group_means, target_f = x, N, VCOV)))
  colnames(means_pop_all) <- colnames(coef(object))
  rownames(means_pop_all) <- paste0("pop_es = ", pop_es)

  # preferred hypothesis
  pref_hypo <- which.max(object$result[, 7])
  pref_hypo_name <- object$result$model[pref_hypo]

  nr_iter <- iter

  if (is.null(quant)) {
    quant <- c(.05, .35, .50, .65, .95)
    names_quant <- c("Sample", "5%", "35%", "50%", "65%", "95%")
  } else {
    names_quant <- c("Sample", paste0(as.character(quant*100), "%"))
  }

  parallel_function_results <- list()

  progressr::handlers(progressr::handler_txtprogressbar(char = ">"))

  progressr::with_progress({
    p <- progressr::progressor(along = seq_len(nr_iter * nr_es))

    for (teller_es in seq_len(nr_es)) {
        cat("Calculating means benchmark for effect-size =", es[teller_es],
            paste0("(", names(es)[teller_es], ")\n"))

      est <- mvtnorm::rmvnorm(n = iter, means_pop_all[teller_es, ], sigma = VCOV)
      # Note that VCOV is the unbiased cov.mx estimate.
      colnames(est) <- names(group_means)

      # Wrapper function for future_lapply
      wrapper_function_asymp <- function(i) {
        p() # Update progress
        parallel_function_asymp(i,
                                est = est, VCOV = VCOV,
                                hypos = hypos, pref_hypo = pref_hypo,
                                comparison = object$comparison, 
                                type = "gorica", # TO DO dus niet ook goricac (evt type = type?); als default en het dus wel kan, dan ook 'sample_nobs' nodig toch... 
                                #type = type,
                                #sample_nobs = sample_nobs,
                                control = control, mix_weights = mix_weights,
                                penalty_factor = penalty_factor, Heq = FALSE, ...)
      }

      name <- paste0("pop_es = ", rnames[teller_es])
      parallel_function_results[[name]] <- future_lapply(
        seq_len(nr_iter),
        wrapper_function_asymp,
        future.seed = TRUE  # Ensures safe and reproducible random number generation
      )
    }
  })


  # get benchmark results
  benchmark_results <- get_results_benchmark(parallel_function_results,
                                             object, pref_hypo,
                                             pref_hypo_name, quant,
                                             names_quant, nr_hypos)

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
    benchmarks_goric_weights = benchmark_results$benchmarks_gw,
    benchmarks_ratio_goric_weights = benchmark_results$benchmarks_rgw,
    benchmarks_ratio_ll_weights = benchmark_results$benchmarks_rlw,
    benchmarks_ratio_ll_ge1 = benchmark_results$benchmarks_rlw_ge1,
    benchmarks_difLL = benchmark_results$benchmarks_difLL,
    benchmarks_absdifLL = benchmark_results$benchmarks_absdifLL,
    combined_values = benchmark_results$combined_values,
    iter = iter
  )

  class(OUT) <- c("benchmark_means", "benchmark", "list")
  return(OUT)
}




## asymp
benchmark_asymp <- function(object, pop_est = NULL, sample_size = NULL, 
                            alt_sample_size = NULL, quant = NULL, iter = 2000,
                            control = list(), 
                            ncpus = 1, seed = NULL, ...) {

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
 
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  
  # keep current plan 
  current_plan <- plan()  
  on.exit(plan(current_plan), add = TRUE)
  if (inherits(current_plan, "sequential")) {
    if (.Platform$OS.type == "windows") {
      plan(future::multisession, workers = ncpus)
    } else {
      plan(future::multicore, workers = ncpus)
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
    Heq = FALSE,
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
  
  nr_iter <- iter
  nr_es  <- nrow(pop_est)
  parallel_function_results <- list()
  
  progressr::handlers(progressr::handler_txtprogressbar(char = ">"))
  progressr::with_progress({
    p <- progressr::progressor(along = seq_len(nr_iter * nr_es))  
    
    for (teller_es in seq_len(nr_es)) {
      cat("Calculating asymptotic benchmark for population estimates =", row.names(pop_est)[teller_es], "\n")
      
      est <- mvtnorm::rmvnorm(n = iter, pop_est[teller_es, ], sigma = VCOV)
      colnames(est) <- names(est_sample)
      
      # Wrapper function for future_lapply
      wrapper_function_asymp <- function(i) {
        p() # Update progress
        parallel_function_asymp(i, 
                                est = est, VCOV = VCOV,
                                hypos = hypos, pref_hypo = pref_hypo, 
                                comparison = comparison, type = "gorica",
                                control = control, mix_weights = mix_weights, 
                                penalty_factor = penalty_factor, Heq = FALSE, ...)
      }
      
      name <- paste0("pop_est = ", rnames[teller_es])
      parallel_function_results[[name]] <- future_lapply(
        seq_len(nr_iter),
        wrapper_function_asymp,
        future.seed = TRUE  # Ensures safe and reproducible random number generation
      )
    }
  })
  
  benchmark_results <- get_results_benchmark(parallel_function_results, object, pref_hypo, 
                                             pref_hypo_name, quant, names_quant, nr_hypos)
  
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
    benchmarks_goric_weights = benchmark_results$benchmarks_gw,
    benchmarks_ratio_goric_weights = benchmark_results$benchmarks_rgw,
    benchmarks_ratio_ll_weights = benchmark_results$benchmarks_rlw,
    benchmarks_ratio_ll_ge1 = benchmark_results$benchmarks_rlw_ge1,
    benchmarks_difLL = benchmark_results$benchmarks_difLL,
    benchmarks_absdifLL = benchmark_results$benchmarks_absdifLL,
    combined_values = benchmark_results$combined_values,
    iter = iter
  )
  
  class(OUT) <- c("benchmark_asymp", "benchmark", "list")
  return(OUT)
}
