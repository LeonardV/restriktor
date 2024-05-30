goric_benchmark_anova <- function(object, pop_es = 0, ratio_pop_means = NULL, 
                                  N = NULL, other_N = NULL, iter = 1000, 
                                  control = list(convergence_crit = 1e-03, 
                                                 chunk_size = 1e4), 
                                  ncpus = 1, seed.value = NULL, ...) {
  
  # Check:
  if (!inherits(object, "con_goric")) {
    return(paste0("Restriktor ERROR: The object should be of class con_goric (a goric object from the goric() function); it belongs to ", class(object)))
  }

  # number of groups
  n.coef <- length(coef(object))
  
  # ES and ratio in data
  if (is.null(object$model.org)) {
    # Number of subjects per group
    if (is.null(N)) {
      other_N <- NULL
      stop("Restriktor Error: please specify the sample-size, e.g. N = 100.")
      # TO DO MAAK ERROR
    } else if (length(N) == 1) {
      samplesize <- rep(N, n.coef)
    } else {
      samplesize <- N
    }
    
    # Unrestricted means
    # est_text <- paste0("object$objectList$", object$objectNames, "$b.unrestr")
    # means <- eval(parse(text = est_text))
    means <- object$objectList[[object$objectNames]]$b.unrestr
    
    # residual variance
    # vcov_text <- paste0("object$objectList$", object$objectNames, "$Sigma")
    # VCOV <- eval(parse(text = vcov_text))
    VCOV <- object$objectList[[object$objectNames]]$Sigma
    
    var.e_data_mx <- VCOV * samplesize
    var.e_data <- var.e_data_mx[1,1] # TO DO check always same elements on diagonal?
  } else {
    # Number of subjects per group
    samplesize <- summary(object$model.org$model[,2])
    # Unrestricted means
    means <- coef(object$model.org)
    # residual variance
    var.e_data <- (sum(object$model.org$residuals^2) / (sum(samplesize) - n.coef))
  }
  
  ES_data <- (1/sqrt(var.e_data)) * sqrt((1/n.coef) * sum((means - mean(means))^2))
  
  ratio_data <- rep(NA, n.coef)
  ratio_data[order(means) == 1] <- 1
  ratio_data[order(means) == 2] <- 2
  d <- means[order(means) == 2] - means[order(means) == 1]
  
  for (i in 1:n.coef) {
    if (order(means)[i] > 2) {
      ratio_data[i] <- 1 + (means[i] - means[order(means) == 1])/d
    }
  }
  
  # effect size population
  es <- pop_es
  nr.es <- length(es)
  
  # ratio population means
  if (is.null(ratio_pop_means)) {
    # Then same as in data
    ratio_pop_means <- ratio_data #coef(object$model.org)
  } else {
    if (length(ratio_pop_means) != n.coef) { 
      return(paste0("The argument ratio_pop_means should be of length ", 
                    n.coef, " (or NULL) but not of length ", length(ratio_pop_means)))
    }
  }
  
  # Hypotheses
  hypos     <- object$hypotheses_usr
  nr.hypos  <- dim(object$result)[1]
  PrefHypo  <- which.max(object$result[,7]) #which.max(object$result$goric.weights)
  pref.hypo <- object$result$model[PrefHypo]
  
  # Error variance
  # var.e <- var(object$model.org$residuals)
  # var.e <- 1
  var.e <- var.e_data
  #
  # When determining pop. means, value does not matter: works exactly the same
  # choose first or last, then pop. means comparable to sample estimates
  
  # Possibly adjust var.e based on other sample size
  if (!is.null(other_N)) {
    if (length(other_N) == 1) {
      var.e <- var.e * (sum(samplesize) - n.coef)
      samplesize <- rep(other_N, n.coef)
      var.e <- var.e / (sum(samplesize) - n.coef)
    } else if (length(other_N) == n.coef) {
      var.e <- var.e * (sum(samplesize) - n.coef)
      samplesize <- other_N
      var.e <- var.e / (sum(samplesize) - n.coef)
    } else {
      return(paste0("The argument other_N should be of length 1 or ", n.coef, " (or NULL) but not of length ", length(other_N)))
    }
  }
  
  means_pop_all <- matrix(NA, ncol = n.coef, nrow = nr.es)
  for (teller.es in 1:nr.es) {
    #teller.es = 1
    
    # Determine mean values, with ratio of ratio.m
    # Solve for x here
    #
    #If all equal, then set population means to all 0
    if (length(unique(ratio_pop_means)) == 1) {
      means_pop <- rep(0, n.coef)
    } else {
      fun <- function (d) {
        means_pop = ratio_pop_means*d
        (1/sqrt(var.e)) * sqrt((1/n.coef) * sum((means_pop - mean(means_pop))^2)) - es[teller.es]
      }
      d <- uniroot(fun, lower = 0, upper = 100)$root
      # Construct means_pop
      means_pop <- ratio_pop_means*d
    }
    
    means_pop_all[teller.es, ] <- means_pop
  }
  colnames(means_pop_all) <- colnames(coef(object))
  rownames(means_pop_all) <- paste0("pop_es = ", pop_es)

  # Create dummies
  sample <- NULL
  D <- as.factor(rep(1:n.coef, times = samplesize))
  sample$D <- D
  sample <- data.frame(D, model.matrix(~ D - 1, data = sample))
  #sample <- dummy_cols(sample, select_columns = 'D')
  colnames(sample)[-1] <- names(coef(object))
  
  nr.iter <- iter
  set.seed(seed.value)
  
  quant <- c(.025, .05, .35, .50, .65, .95, .975)
  names_quant <- c("Sample", "2.5%", "5%", "35%", "50%", "65%", "95%", "97.5%")
  
  CI.benchmarks_all <- NULL
  CI.benchmarks_gw_all <- NULL
  CI.benchmarks_lw_all <- NULL
  CI.benchmarks_lw_ge1_all <- NULL
  CI.benchmarks_ld_all <- NULL
  CI.benchmarks_ld_ge0_all <- NULL
  
  for (teller.es in 1:nr.es) {
    #teller.es = 1
    means_pop <- means_pop_all[teller.es, ]
    
    # parallel backend
    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))

    # Export required variables to cluster nodes
    parallel::clusterExport(cl, c("samplesize", "var.e", "nr.iter", "means_pop", 
                                  "hypos", "PrefHypo", "object", "n.coef", "sample",
                                  "control"), 
                            envir = environment())

    parallel::clusterEvalQ(cl, {
      library(restriktor) 
    })
    
    # Create a function that wraps parallel_function
    wrapper_function <- function(i) {
      parallel_function(i, samplesize = samplesize, var.e = var.e, nr.iter = nr.iter, 
                        means_pop = means_pop, hypos = hypos, PrefHypo = PrefHypo, 
                        object = object, n.coef = n.coef, sample = sample,
                        control = control, ...)
    }
    
    
    # Call the parallel function using 'parLapply'
    #results <- parLapply(cl, 1:10, wrapper_function)
    
    # Use pbapply::pblapply with arguments
    pbapply::pboptions(type = "timer", style = 1, char = ">")
    results <- pbapply::pblapply(1:nr.iter, wrapper_function, cl = cl)
    
    # Stop the parallel backend
    parallel::stopCluster(cl)
    
    #sapply(results, function(result) result$test)
    
    # combine results from parallel process
    goric <- sapply(results, function(result) result$goric)
    gw    <- do.call(cbind, lapply(results, function(result) result$gw))
    lw    <- do.call(cbind, lapply(results, function(result) result$lw))
    ld    <- do.call(cbind, lapply(results, function(result) result$ld))
    
    CI.benchmarks_goric <- matrix(c(object$result[PrefHypo,7], quantile(goric, quant)), 
                                  nrow = 1) # sample weight with calculated quantiles/percentiles
    colnames(CI.benchmarks_goric) <- names_quant
    rownames(CI.benchmarks_goric) <- pref.hypo
    #
    #
    CI.benchmarks_gw <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    CI.benchmarks_lw <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    CI.benchmarks_lw_ge1 <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    CI.benchmarks_ld <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    CI.benchmarks_ld_ge0 <- matrix(NA, nrow = nr.hypos, ncol = 1+length(quant))
    #
    CI.benchmarks_gw[,1] <- object$ratio.gw[PrefHypo,] # so in sample
    CI.benchmarks_lw[,1] <- object$ratio.lw[PrefHypo,] # so in sample
    
    for (j in 1:nr.hypos) {
      if (object$ratio.lw[PrefHypo,j] >= 1) {
        CI.benchmarks_lw_ge1[j,1] <- object$ratio.lw[PrefHypo,j] # so in sample
      } else {
        CI.benchmarks_lw_ge1[j,1] <- 1/object$ratio.lw[PrefHypo,j] # so in sample
      }
    }
    CI.benchmarks_ld[,1] <- (object$result$loglik[PrefHypo] - object$result$loglik) # so in sample
    CI.benchmarks_ld_ge0[,1] <- abs(object$result$loglik[PrefHypo] - object$result$loglik) # so in sample
    #
    lw_ge1 <- lw
    lw_ge1[lw < 1] <- 1/lw[lw < 1]
    ld_ge0 <- abs(ld)
    
    for (j in 1:nr.hypos) {
      CI.benchmarks_gw[j,2:(1+length(quant))] <- quantile(gw[j,], quant)
      CI.benchmarks_lw[j,2:(1+length(quant))] <- quantile(lw[j,], quant)
      CI.benchmarks_lw_ge1[j,2:(1+length(quant))] <- quantile(lw_ge1[j,], quant)
      CI.benchmarks_ld[j,2:(1+length(quant))] <- quantile(ld[j,], quant)
      CI.benchmarks_ld_ge0[j,2:(1+length(quant))] <- quantile(ld_ge0[j,], quant)
    }
    #
    colnames(CI.benchmarks_gw) <- colnames(CI.benchmarks_lw) <- 
      colnames(CI.benchmarks_lw_ge1) <- colnames(CI.benchmarks_ld) <- 
      colnames(CI.benchmarks_ld_ge0) <- names_quant
    #
    rownames(CI.benchmarks_gw) <- rownames(CI.benchmarks_lw) <- 
      rownames(CI.benchmarks_lw_ge1) <- rownames(CI.benchmarks_ld) <- 
      rownames(CI.benchmarks_ld_ge0) <- paste(pref.hypo, names(object$ratio.gw[PrefHypo,]))

    name <- paste0("pop.es = ", pop_es[teller.es])
    CI.benchmarks_all[[name]] <- CI.benchmarks_goric
    CI.benchmarks_gw_all[[name]] <- CI.benchmarks_gw
    CI.benchmarks_lw_all[[name]] <- CI.benchmarks_lw
    CI.benchmarks_lw_ge1_all[[name]] <- CI.benchmarks_lw_ge1
    CI.benchmarks_ld_all[[name]] <- CI.benchmarks_ld
    CI.benchmarks_ld_ge0_all[[name]] <- CI.benchmarks_ld_ge0
  }
  
  # Error probability based on complement of preferred hypothesis in data
  if (nr.hypos == 2 && object$comparison == "complement") {
    if (object$type == 'goric') {
      error.prob <- 1 - object$result$goric.weights[PrefHypo]
    } else {
      error.prob <- 1 - object$result$gorica.weights[PrefHypo]
    }
  } else {
    if (PrefHypo == nr.hypos && object$comparison == "unconstrained") {
      error.prob <- "The unconstrained (i.e., the failsafe) containing all possible orderings is preferred..."
    } else {
      H_pref <- hypos[[PrefHypo]]
      if (is.null(object$model.org)) {
        results.goric_pref <- goric(means, VCOV = VCOV,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = object$type, 
                                    control = control, 
                                    ...)
      } else {
        fit_data <- object$model.org
        results.goric_pref <- goric(fit_data,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = object$type,
                                    control = control, 
                                    ...)
      }
      if (object$type == 'goric') {
        error.prob <- results.goric_pref$result$goric.weights[2]
      } else {
        error.prob <- results.goric_pref$result$gorica.weights[2]
      }
    }
  }
  #
  # TO DO bepaal ook quantiles voor error prob. Ws verwerken in bovenstaande!
  # Is het zinnig? Op zich zou error prob al zinnig moeten zijn immers....
  
  final <- list(
    n.coef = n.coef,
    group.size = samplesize,
    means.data = means, ratio.means.data = ratio_data, ES.data = ES_data,
    res.var.data = var.e_data,
    pop_es = pop_es, pop.means = means_pop_all,
    ratio_pop_means = ratio_pop_means,
    res.var.pop = var.e,
    pref.hypo = pref.hypo, error.prob.pref.hypo = error.prob,
    benchmarks.weight = CI.benchmarks_all,
    benchmarks.ratios = CI.benchmarks_gw_all,
    benchmarks.LLratios = CI.benchmarks_lw_all,
    benchmarks.LLratios_ge1 = CI.benchmarks_lw_ge1_all,
    benchmarks.difLL = CI.benchmarks_ld_all,
    benchmarks.absdifLL = CI.benchmarks_ld_ge0_all)
  
  class(final) <- c("benchmarks", "list")
  final
  
} 
