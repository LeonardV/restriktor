benchmark_means  <- function(object, ...) UseMethod("benchmark_means")
benchmark_asymp  <- function(object, ...) UseMethod("benchmark_asymp")


benchmark <- function(object, model_type = c("asymp", "means"), ...) {

  args <- list(...)
  
  # if (!is.null(args$pop_es)) {
  #   model_type <- "means"
  # }
  
  model_type <- match.arg(model_type, c("asymp", "means"))
  if (is.null(model_type)) {
    stop("Restriktor ERROR: Please specify if you want to benchmark means or asymptotic result ",
         "In case of model_type = means (default), no intercept is allowed.")
  }
  
  # Check if the model has an intercept. Since the goric function accepts both 
  # a model/fit or only estimates+VCOV, we cannot rely on the original model fit.
  # So we only check if the names of the vector with parameter estimates includes
  # the word \(Intercept\).
  has_intercept <- detect_intercept(object)
  if (model_type == "means") {
    if (has_intercept) { 
      stop("Restriktor ERROR: A model with an intercept is not allowed for model_type = means. ",
           "Please refit the model without an intercept.")
    }
    benchmark_means(object, ...)  
  } else if (model_type == "asymp") {
    benchmark_asymp(object, ...)
  }
}
  

benchmark_means <- function(object, pop_es = NULL, ratio_pop_means = NULL, 
                            group_size = NULL, alt_group_size = NULL, 
                            quant = NULL, iter = 1000, 
                            control = list(convergence_crit = 1e-03, 
                                           chunk_size = 1e4), 
                            ncpus = 1, cl = NULL, seed = NULL, ...) {

  
  # group_size is needed to rescale vcov based on alt_group_size.
  
  # Check:
  if (!inherits(object, "con_goric")) {
    stop(paste("Restriktor ERROR:", 
               "The object should be of class 'con_goric' (a goric(a) object from",
               "the goric() function). However, it belongs to the following class(es):", 
      paste(class(object), collapse = ", ")
    ), call. = FALSE)
  }
  
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  # Hypotheses
  hypos <- object$hypotheses_usr
  nr_hypos <- dim(object$result)[1]
  pref_hypo <- which.max(object$result[, 7]) 
  pref_hypo_name <- object$result$model[pref_hypo]
  
  # number of groups
  ngroups <- length(coef(object))
  
  # does original model fit exist
  form_model_org <- formula(object$model.org)
  
  # ES and ratio in data
  if (is.null(object$model.org)) {
    # Number of subjects per group
    if (is.null(group_size)) {
      stop("Restriktor Error: please specify the sample-size, e.g. group_size = 100.", 
           call. = FALSE)
    } else if (length(group_size) == 1) {
      N <- rep(group_size, ngroups)
    } else {
      N <- group_size
    }
    # Unrestricted (adjusted) group_means
    group_means <- object$b.unrestr
    # residual variance
    VCOV <- object$Sigma
    var_e_data <- diag(VCOV * N)[1] # are the elements on the diagonal always equal?
  } else {
    # Number of subjects per group
    N <- colSums(model.matrix(object$model.org)) #summary(object$model.org$model[, 2])
    # Unrestricted group_means
    group_means <- coef(object$model.org)
    VCOV <- vcov(object$model.org)
    # residual variance
    var_e_data <- sum(object$model.org$residuals^2) / (sum(N) - ngroups)
  }

  ## Compute observed Cohens f
  cohens_f_observed <- compute_cohens_f(group_means, N, VCOV)
  
  # Compute ratio data based on group_means
  ratio_data <- compute_ratio_data(group_means)
  
  # ratio population means
  if (is.null(ratio_pop_means)) {
    # Then same as in data
    ratio_pop_means <- ratio_data 
  } else {
    if (length(ratio_pop_means) != ngroups) { 
      stop(paste0("The argument ratio_pop_means should be of length ", 
                  ngroups, " (or NULL) but not of length ", length(ratio_pop_means)))
    }
  }
  
  # Error variance
  # var_e <- var(object$model.org$residuals)
  # var_e <- 1
  var_e <- as.vector(var_e_data)
  
  # When determining pop.means, value does not matter: works exactly the same
  # choose first or last, then pop. means comparable to sample estimates
  
  # Possibly adjust var_e based on other sample size
  if (!is.null(alt_group_size)) {
    result_adjust_variance <- adjust_variance(var_e, N, alt_group_size = alt_group_size, ngroups)
    N <- result_adjust_variance$N
    var_e <- result_adjust_variance$var_e
  }

  # effect size population
  if (is.null(pop_es)) {
    pop_es <- c(0, round(cohens_f_observed, 3))
  } else {
    pop_es <- sort(pop_es)  
  }
  
  es <- pop_es
  nr_es <- length(es)
  
  means_pop_all <- compute_population_means(pop_es, ratio_pop_means, var_e, ngroups)
  colnames(means_pop_all) <- colnames(coef(object))
  rownames(means_pop_all) <- paste0("pop_es = ", pop_es)
  
  # Create dummies
  sample <- data.frame(D = as.factor(rep(1:ngroups, times = N)))
  sample <- data.frame(sample$D, model.matrix(~ D - 1, data = sample))
  colnames(sample)[-1] <- names(coef(object))
  
  nr_iter <- iter
  
  if (is.null(quant)) {
    quant <- c(.05, .35, .50, .65, .95)
    names_quant <- c("Sample", "5%", "35%", "50%", "65%", "95%")
  } else {
    names_quant <- c("Sample", paste0(as.character(quant*100), "%"))
  }

  # parallel backend
  if (is.null(cl)) {
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
  }
  on.exit(parallel::stopCluster(cl))
  
  
  parallel::clusterEvalQ(cl, {
    library(restriktor) 
  })

  # Export required variables to cluster nodes
  parallel::clusterExport(cl, c("N", "var_e", 
                                "hypos", "pref_hypo", 
                                "object", "ngroups", "sample",
                                "control", "form_model_org"), 
                          envir = environment())
  
  # Use pbapply::pblapply with arguments
  pbapply::pboptions(type = "timer", style = 1, char = ">")
  
  parallel_function_results <- list()  
  for (teller_es in 1:nr_es) {
    cat("Calculating means benchmark for effect-size =", pop_es[teller_es], "\n")
    
    # Update the means_pop variable within the cluster nodes
    means_pop <- means_pop_all[teller_es, ]
    parallel::clusterExport(cl, c("means_pop"), envir = environment())
    
    # main function
    # Create a function that wraps parallel_function
    wrapper_function_means <- function(i) {
      parallel_function_means(i, N = N, var_e = var_e, 
                              means_pop = means_pop, 
                              hypos = hypos, pref_hypo = pref_hypo, 
                              object = object, ngroups = ngroups, 
                              sample = sample, control = control, 
                              form_model_org = form_model_org, 
                              ...)
    }
    
    name <- paste0("pop_es = ", pop_es[teller_es])
    parallel_function_results[[name]] <- pbapply::pblapply(seq_len(nr_iter), 
                                                           wrapper_function_means, 
                                                           cl = cl)
  }
  
  # get benchmark results
  benchmark_results <- get_results_benchmark(parallel_function_results, 
                                             object, pref_hypo, 
                                             pref_hypo_name, quant, 
                                             names_quant, nr_hypos)

  # compute error probability
  error_prob <- calculate_error_probability(object, hypos, pref_hypo, 
                                            est = group_means, 
                                            VCOV, control, ...)

  OUT <- list(
    type = object$type,
    comparison = object$comparison,
    ngroups = ngroups,
    group_size = N,
    group_means_observed = group_means, 
    ratio_group_means.data = ratio_data, 
    cohens_f_observed = cohens_f_observed,
    res_var_observed = var_e_data,
    pop_es = pop_es, pop_group_means = means_pop_all,
    ratio_pop_means = ratio_pop_means,
    res_var_pop = var_e,
    pref_hypo_name = pref_hypo_name, 
    error_prob_pref_hypo_name = error_prob,
    benchmarks_goric_weights = benchmark_results$benchmarks_gw,
    benchmarks_ratio_goric_weights = benchmark_results$benchmarks_rgw,
    benchmarks_ratio_ll_weights = benchmark_results$benchmarks_rlw,
    benchmarks_ratio_ll_ge1 = benchmark_results$benchmarks_rlw_ge1,
    benchmarks_difLL = benchmark_results$benchmarks_difLL,
    benchmarks_absdifLL = benchmark_results$benchmarks_absdifLL,
    combined_values = benchmark_results$combined_values
    )
  
  class(OUT) <- c("benchmark_means", "benchmark", "list")
  return(OUT)
} 




benchmark_asymp <- function(object, pop_est = NULL, sample_size = NULL, 
                            alt_sample_size = NULL, quant = NULL, iter = 1000, 
                            control = list(convergence_crit = 1e-03, 
                                           chunk_size = 1e4), 
                            ncpus = 1, cl = NULL, seed = NULL, ...) {
  
  # Check if object is of class con_goric
  if (!inherits(object, "con_goric")) {
    stop(paste("Restriktor ERROR:", 
               "The object should be of class 'con_goric' (a goric(a) object from the goric() function).",
               "However, it belongs to the following class(es):", 
               paste(class(object), collapse = ", ")), call. = FALSE)
  }
  
  if (!is.null(seed)) set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  
  hypos <- object$hypotheses_usr
  # print warning if a constants is found in the constraints list
  check_rhs_constants(object$rhs)
  
  nr_hypos <- dim(object$result)[1]
  pref_hypo <- which.max(object$result[, 7])
  pref_hypo_name <- object$result$model[pref_hypo]
  comparison <- object$comparison
  type <- object$type
  
  est_sample <- object$b.unrestr
  n_coef <- length(est_sample)
  
  if (is.null(pop_est)) {
    pop_est <- matrix(rbind(rep(0, n_coef), round(est_sample, 3)), nrow = 2)
    row.names(pop_est) <- c("No-effect", "Observed")
  } else {
    if (is.vector(pop_est)) {
      pop_est <- matrix(pop_est, nrow = 1, ncol = length(pop_est))
    }
  }
  
  # Check if pop_est and est_sample have the same length
  if (ncol(pop_est) != length(est_sample)) {
    stop(paste("Restriktor Error: The number of columns in pop_est (", ncol(pop_est), 
               ") does not match the length of est_sample (", length(est_sample), ").", sep=""))
  }
  
  # Assign row names to pop_est if they are null or empty strings
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
  
  VCOV <- object$Sigma
  N <- length(object$model.org$residuals) 
  if (N == 0) {
    N <- ""
  }
  
  if (!is.null(alt_sample_size)) {
    if (is.null(sample_size)) {
      stop("Restriktor Error: Please provide the original sample size(s) using the argument `sample_size`. This information is required to rescale the variance-covariance matrix (VCOV).", call. = FALSE)
    }
    N <- sample_size
    VCOV <- VCOV * N / alt_sample_size
    N <- alt_sample_size
  }
  
  if (is.null(quant)) {
    quant <- c(.05, .35, .50, .65, .95)
    names_quant <- c("Sample", "5%", "35%", "50%", "65%", "95%")
  } else {
    names_quant <- c("Sample", paste0(as.character(quant*100), "%"))
  }
  
  # parallel backend
  if (is.null(cl)) {
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
  }
  on.exit(parallel::stopCluster(cl))
  
  parallel::clusterEvalQ(cl, {
    library(restriktor) 
  })
  
  # Export required variables to cluster nodes
  parallel::clusterExport(cl, c("VCOV", "hypos", "pref_hypo", "type",
                                "comparison", "control"), 
                          envir = environment())
  
  # Use pbapply::pblapply with arguments
  pbapply::pboptions(type = "timer", style = 1, char = ">")
  
  nr_iter <- iter
  nr_es  <- nrow(pop_est) #length(pop_est) / n_coef
  parallel_function_results <- list()  
  
    for (teller_es in seq_len(nr_es)) {
    cat("Calculating asymptotic benchmark for population estimates =", row.names(pop_est)[teller_es], "\n")
    
    est <- mvtnorm::rmvnorm(n = iter, pop_est[teller_es, ], sigma = VCOV)
    
    # Update the est variable within the cluster nodes
    parallel::clusterExport(cl, c("est"), envir = environment())
    
    # main function
    # Create a function that wraps parallel_function
    wrapper_function_asymp <- function(i) {
      parallel_function_asymp(i, 
                              est = est, VCOV = VCOV,
                              hypos = hypos, pref_hypo = pref_hypo, 
                              comparison = comparison, type = "gorica",
                              control = control, ...)
    }
    
    name <- paste0("pop_est = ", rnames[teller_es])
    parallel_function_results[[name]] <- pbapply::pblapply(seq_len(nr_iter), 
                                                           wrapper_function_asymp, 
                                                           cl = cl)
  }
  
  # get benchmark results
  benchmark_results <- get_results_benchmark(parallel_function_results, 
                                             object, pref_hypo, 
                                             pref_hypo_name, quant, 
                                             names_quant, nr_hypos)
  
  
  # compute error probability
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
    error_prob_pref_hypo_name = error_prob,
    benchmarks_goric_weights = benchmark_results$benchmarks_gw,
    benchmarks_ratio_goric_weights = benchmark_results$benchmarks_rgw,
    benchmarks_ratio_ll_weights = benchmark_results$benchmarks_rlw,
    benchmarks_ratio_ll_ge1 = benchmark_results$benchmarks_rlw_ge1,
    benchmarks_difLL = benchmark_results$benchmarks_difLL,
    benchmarks_absdifLL = benchmark_results$benchmarks_absdifLL,
    combined_values = benchmark_results$combined_values
  )
  
  class(OUT) <- c("benchmark_asymp", "benchmark", "list")
  return(OUT)
}


