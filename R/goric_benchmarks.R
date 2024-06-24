benchmark_means  <- function(object, ...) UseMethod("benchmark_means")
benchmark_asymp  <- function(object, ...) UseMethod("benchmark_asymp")


benchmark <- function(object, model_type = NULL, ...) {

  if (is.null(model_type)) {
    stop("Restriktor ERROR: Please specify if you want to benchmark means or not. ",
         "In case of model_type = means, no intercept is allowed.")
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

  
  # is het nodig om zowel N als other_N te gebruiken, volstaat other_N niet gewoon?
  # Als object = model dan wordt var_e herschaald ahdv other_N
  # Als object = est+vcov dan wordt voor var_e other_N gebruikt. 
  # Ja, je hebt group_size nodig om te herschalen obv alt_group_size.
  
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
  #group_means <- coef(model0)
  # Totale gemiddelde berekenen
  total_mean <- sum(group_means * N) / sum(N) 
  # Tussen-groep variantie berekenen
  ss_between <- sum(N * (group_means - total_mean)^2) 
  # Covariantiematrix (variances on the diagonal)
  cov_matrix <- VCOV * N
  # Binnen-groep variantie berekenen 
  ss_within <- sum((N - 1) * diag(cov_matrix)) 
  # Cohen's f berekenen
  cohens_f_observed <- sqrt(ss_between/ss_within)
  
  ratio_data <- rep(NA, ngroups)
  ratio_data[order(group_means) == 1] <- 1
  ratio_data[order(group_means) == 2] <- 2
  # The choice of the smallest and the second smallest mean makes the scaling 
  # more robust against changes in the other group means. Since these values 
  # represent the lower bound of the data, the scale is less sensitive to the 
  # spread of higher values.
  
  # For example:
  # The value of 2.28 indicates that this particular group mean is 2.28 times 
  # the scale factor d above the smallest mean. This means the group mean is 
  # further from the smallest mean compared to the second smallest mean, and 
  # helps in understanding the relative differences between the group means in 
  # a normalized manner.
  d <- group_means[order(group_means) == 2] - group_means[order(group_means) == 1]
  
  for (i in seq_len(ngroups)) {
    if (order(group_means)[i] > 2) {
      ratio_data[i] <- 1 + (group_means[i] - group_means[order(group_means) == 1]) / d
    }
  }
  # ratio population means
  if (is.null(ratio_pop_means)) {
    # Then same as in data
    ratio_pop_means <- ratio_data 
  } else {
    if (length(ratio_pop_means) != ngroups) { 
      return(paste0("The argument ratio_pop_means should be of length ", 
                    ngroups, " (or NULL) but not of length ", length(ratio_pop_means)))
    }
  }
  
  # Hypotheses
  hypos <- object$hypotheses_usr
  nr_hypos <- dim(object$result)[1]
  pref_hypo <- which.max(object$result[, 7]) 
  pref_hypo_name <- object$result$model[pref_hypo]
  
  # Error variance
  # var_e <- var(object$model.org$residuals)
  # var_e <- 1
  var_e <- as.vector(var_e_data)
  
  # When determining pop.means, value does not matter: works exactly the same
  # choose first or last, then pop. means comparable to sample estimates
  
  # Possibly adjust var_e based on other sample size
  if (!is.null(alt_group_size)) {
    if (length(alt_group_size) == 1) {
      var_e <- var_e * (sum(N) - ngroups)
      N <- rep(alt_group_size, ngroups)
      var_e <- var_e / (sum(N) - ngroups)
    } else if (length(alt_group_size) == ngroups) {
      var_e <- var_e * (sum(N) - ngroups)
      N <- alt_group_size
      var_e <- var_e / (sum(N) - ngroups)
    } else {
      return(paste0("The argument alt_group_size should be of length 1 or ", 
                    ngroups, " (or NULL) but not of length ", length(alt_group_size)))
    }
  }
  
  # effect size population
  if (is.null(pop_es)) {
    pop_es <- c(0, round(cohens_f_observed, 3))
  }
  
  es <- pop_es
  nr_es <- length(es)
  
  means_pop_all <- matrix(NA, ncol = ngroups, nrow = nr_es)
  for (teller_es in seq_len(nr_es)) {
    #teller_es = 1
    
    # Determine mean values, with ratio of ratio.m
    # Solve for x here
  
    # If all equal, then set population means to all 0
    if (length(unique(ratio_pop_means)) == 1) {
      means_pop <- rep(0, ngroups)
    } else {
      fun <- function (d) {
        means_pop = ratio_pop_means * d
        (1/sqrt(var_e)) * sqrt((1/ngroups) * sum((means_pop - mean(means_pop))^2)) - es[teller_es] #  AANPASSSEN NAAR NIEUWE FORMULE
      }
      d <- uniroot(fun, lower = 0, upper = 100)$root
      # Construct means_pop
      means_pop <- ratio_pop_means*d
    }
    means_pop_all[teller_es, ] <- means_pop
  }
  colnames(means_pop_all) <- colnames(coef(object))
  rownames(means_pop_all) <- paste0("pop_es = ", pop_es)

  # Create dummies
  sample <- data.frame(D = as.factor(rep(1:ngroups, times = N)))
  sample <- data.frame(sample$D, model.matrix(~ D - 1, data = sample))
  colnames(sample)[-1] <- names(coef(object))
  
  nr_iter <- iter
  
  if (is.null(quant)) {
    quant <- c(.025, .05, .35, .50, .65, .95, .975)
    names_quant <- c("Sample", "2.5%", "5%", "35%", "50%", "65%", "95%", "97.5%")
  } else {
    names_quant <- c("Sample", paste0(as.character(quant*100), "%"))
  }

  # parallel backend
  if (is.null(cl)) {
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
  }
  on.exit(parallel::stopCluster(cl))
  
  
  # Export required variables to cluster nodes
  parallel::clusterExport(cl, c("N", "var_e", "means_pop", 
                                "hypos", "pref_hypo", 
                                "object", "ngroups", "sample",
                                "control", "form_model_org"), 
                          envir = environment())
  
  parallel::clusterEvalQ(cl, {
    library(restriktor) 
  })

  # Use pbapply::pblapply with arguments
  pbapply::pboptions(type = "timer", style = 1, char = ">")
  
  parallel_function_results <- list()  
  for (teller_es in 1:nr_es) {
    cat("Calculating means benchmark for effect-size =", pop_es[teller_es], "\n")
    #teller_es = 1
    means_pop <- means_pop_all[teller_es, ]
    
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

  # Error probability based on complement of preferred hypothesis in data
  if (nr_hypos == 2 && object$comparison == "complement") { 
    if (object$type == 'goric') {
      error_prob <- 1 - object$result$goric.weights[pref_hypo]
    } else {
      error_prob <- 1 - object$result$gorica.weights[pref_hypo]
    }
  } else {
    if (pref_hypo == nr_hypos && object$comparison == "unconstrained") {
      error_prob <- "The unconstrained (i.e., the failsafe) containing all possible orderings is preferred."
    } else {
      H_pref <- hypos[[pref_hypo]]
      if (is.null(object$model.org)) {
        results_goric_pref <- goric(group_means, VCOV = VCOV,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = object$type, 
                                    control = control, 
                                    ...)
      } else {
        fit_data <- object$model.org
        results_goric_pref <- goric(fit_data,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = object$type,
                                    control = control, 
                                    ...)
      }
      if (object$type == 'goric') {
        error_prob <- results_goric_pref$result$goric.weights[2]
      } else {
        error_prob <- results_goric_pref$result$gorica.weights[2]
      }
    }
  }

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
  }
  
  rnames <- row.names(pop_est)
  if (is.null(rnames)) {
    rnames <- paste0("ES_", rep(1:nrow(pop_est)))
    row.names(pop_est) <- rnames
  }
  colnames(pop_est) <- names(est_sample)
  
  VCOV <- object$Sigma
  N <- length(object$model.org$residuals) 
  
  if (is.null(object$model.org$residuals)) {
    if (is.null(sample_size)) {
      stop("Restriktor Error: please specify the sample-size(s), e.g. sample_size = 100.", call. = FALSE)
    } else {
      N <- sample_size
    }
  }
  
  if (!is.null(alt_sample_size)) {
    VCOV <- VCOV * N / alt_sample_size
    N <- alt_sample_size
  }
  
  if (is.null(quant)) {
    quant <- c(.025, .05, .35, .50, .65, .95, .975)
    names_quant <- c("Sample", "2.5%", "5%", "35%", "50%", "65%", "95%", "97.5%")
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
  
  # Use pbapply::pblapply with arguments
  pbapply::pboptions(type = "timer", style = 1, char = ">")
  
  nr_iter <- iter
  nr_es  <- nrow(pop_est) #length(pop_est) / n_coef
  parallel_function_results <- list()  
  
  for (teller_es in seq_len(nr_es)) {
    cat("Calculating asymptotic benchmark for population values =", pop_est[teller_es, ], "\n")
    
    est <- mvtnorm::rmvnorm(n = iter, pop_est[teller_es, ], sigma = VCOV)
    
    # Export required variables to cluster nodes
    parallel::clusterExport(cl, c("est", "VCOV", "hypos", "pref_hypo", "type",
                                  "comparison", "control"), 
                            envir = environment())
    
    # main function
    # Create a function that wraps parallel_function
    wrapper_function_asymp <- function(i) {
      parallel_function_asymp(i, 
                              est = est, VCOV = VCOV,
                              hypos = hypos, pref_hypo = pref_hypo, 
                              comparison = comparison, type = type,
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
  
  # Error probability based on complement of preferred hypothesis in data
  if (nr_hypos == 2 && object$comparison == "complement") { 
    if (object$type == 'goric') {
      error_prob <- 1 - object$result$goric.weights[pref_hypo]
    } else {
      error_prob <- 1 - object$result$gorica.weights[pref_hypo]
    }
  } else {
    if (pref_hypo == nr_hypos && object$comparison == "unconstrained") {
      error_prob <- "The unconstrained (i.e., the failsafe) containing all possible orderings is preferred."
    } else {
      H_pref <- hypos[[pref_hypo]]
      if (is.null(object$model.org)) {
        results_goric_pref <- goric(est_sample, VCOV = VCOV,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = object$type, 
                                    control = control, 
                                    ...)
      } else {
        fit_data <- object$model.org
        results_goric_pref <- goric(fit_data,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = object$type,
                                    control = control, 
                                    ...)
      }
      if (object$type == 'goric') {
        error_prob <- results_goric_pref$result$goric.weights[2]
      } else {
        error_prob <- results_goric_pref$result$gorica.weights[2]
      }
    }
  }
  
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
