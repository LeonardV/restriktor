goric_benchmark_asymp <- function(object, pop_est = NULL, sample_size = NULL, 
                                  alt_sample_size = NULL, quant = NULL, iter = 1000, 
                                  control = list(convergence_crit = 1e-03, 
                                                 chunk_size = 1e4), 
                                  ncpus = 1, cl = NULL, seed = NULL, ...) {
        
  # Check if object is of class con_goric
  if (!inherits(object, "con_goric")) {
    stop(paste("Restriktor ERROR:", 
               "The object should be of class 'con_goric' (a goric object from the goric() function).",
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
  
  if (is.null(object$model.org)) {
    if (is.null(pop_est)) { 
      pop_est <- matrix(rbind(rep(0, n_coef), round(est_sample, 3)), nrow = 2)
    } else if (is.data.frame(pop_est)) {
      pop_est <- as.matrix(pop_est)
    } else if (!is.matrix(pop_est)) {
      pop_est <- matrix(pop_est, nrow = 1)
    }
    colnames(pop_est) <- names(est_sample)
    VCOV <- object$Sigma
    if (is.null(sample_size)) {
      alt_sample_size <- NULL
      stop("Restriktor Error: please specify the sample-size, e.g. group_size = 100.", 
           call. = FALSE)
    } else {
      N <- sample_size
    }
  } else {
    if (is.null(pop_est)) {
      pop_est <- est_sample
      pop_est <- matrix(pop_est, nrow = 1)
    }
    colnames(pop_est) <- names(object$model.org$coefficients)
    VCOV <- vcov(object$model.org)
    N <- length(object$model.org$residuals)
  }
  
  nr_es  <- nrow(pop_est) #length(pop_est) / n_coef
  
  if (!is.null(alt_sample_size)) {
    VCOV <- VCOV * N / alt_sample_size
    N <- alt_sample_size
  }
  
  nr_iter <- iter
  
  if (is.null(quant)) {
    quant <- c(.025, .05, .35, .50, .65, .95, .975)
    names_quant <- c("Sample", "2.5%", "5%", "35%", "50%", "65%", "95%", "97.5%")
  } else {
    names_quant <- c("Sample", paste0(as.character(quant*100), "%"))
  }
  
  rnames <- row.names(pop_est)
  if (is.null(rnames)) {
    rnames <- as.character(rep(1:nrow(pop_est)))
    row.names(pop_est) <- rnames
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
