# used in print.benchmark()
capitalize_first_letter <- function(input_string) {
  paste0(toupper(substring(input_string, 1, 1)), substring(input_string, 2))
}

# used in get_results_benchmark_means()
remove_single_value_rows <- function(data, value) {
  rows_to_keep <- apply(data, 1, function(row) !all(row == value))
  data[rows_to_keep, , drop = FALSE]
}  

remove_single_value_col <- function(data, value) {
  cols_to_keep <- apply(data, 2, function(col) !all(col == value))
  data[, cols_to_keep, drop = FALSE]
}  

# Function to filter columns based on exact matching of hypothesis_comparison
filter_columns <- function(data, hypothesis_comparison) {
  colnames_to_keep <- sapply(hypothesis_comparison, function(hypo) {
    grep(paste0("\\b", hypo, "\\b"), names(data), value = TRUE, fixed = FALSE)
  })
  colnames_to_keep <- unique(unlist(colnames_to_keep))
  colnames_to_keep <- intersect(names(data), colnames_to_keep) # Preserve original order
  return(data[, colnames_to_keep, drop = FALSE])
}

filter_vector <- function(data, hypothesis_comparison) {
  names_to_keep <- sapply(hypothesis_comparison, function(hypo) {
    grep(paste0("\\b", hypo, "\\b"), names(data), value = TRUE, fixed = FALSE)
  })
  names_to_keep <- unique(unlist(names_to_keep))
  names_to_keep <- intersect(names(data), names_to_keep) # Preserve original order
  return(data[names_to_keep])
}


# Functie om de power te berekenen voor een enkele density en critical value
calculate_power <- function(density_h1, critical_value) {
  sum(density_h1$y[density_h1$x > critical_value]) * mean(diff(density_h1$x))
}

# Function to calculate density
# calculate_density <- function(data, var, sample_value) {
#   dens <- density(data[[var]], kernel = "gaussian", na.rm = TRUE, bw = "nrd0")
#   data.frame(x = dens$x, y = dens$y, sample_value = sample_value)
# }


# Extract the parts within the parentheses
extract_in_parentheses <- function(name) {
  regmatches(name, regexpr("\\(.*\\)", name))
}

remove_spaces_in_parentheses <- function(string) {
  string <- gsub("\\(\\s+", "(", string)
  string <- gsub("\\s+\\)", ")", string)
  return(string)
}

construct_colnames <- function(list_name, colnames, pref_hypo_name) {
  remove_spaces_in_parentheses(paste0(list_name, " (", pref_hypo_name, " ", colnames, ")", sep = ""))
}

combine_matrices_cbind <- function(lst) {
  df_list <- lapply(lst, as.data.frame)
  combined_df <- do.call(cbind, df_list)
  return(combined_df)
}


combine_matrices_cbind <- function(lst) {
  # Vind het maximum aantal rijen in de lijst
  max_rows <- max(sapply(lst, nrow))
  
  # Zet de matrices om naar data frames en vul met NA waar nodig
  df_list <- lapply(lst, function(x) {
    as.data.frame(rbind(x, matrix(NA, nrow = max_rows - nrow(x), ncol = ncol(x))))
  })
  
  # Combineer de data frames kolomgewijs
  combined_df <- do.call(cbind, df_list)
  
  return(combined_df)
}


# model_type = "means" ----------------------------------------------------
detect_intercept <- function(model) {
  coefficients <- model$b.unrestr
  intercept_names <- c("(Intercept)", "Intercept", "(const)", "const",
                       "(Int)", "Int", "(Cons)", "Cons", "b0", "beta0")

  names_lower <- tolower(names(coefficients))
  intercept_names_lower <- tolower(intercept_names)

  detected_intercepts <- intercept_names[names_lower %in% intercept_names_lower]

  if (length(detected_intercepts) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Compute Cohen's f based on group_means, N, and VCOV
compute_cohens_f <- function(group_means, N, VCOV) {
  total_mean <- sum(group_means * N) / sum(N)
  ss_between <- sum(N * (group_means - total_mean)^2)
  cov_matrix <- VCOV * (N - 1) # covmx based on N instead of N-1
  ss_within <- sum(N * diag(cov_matrix)) # equates: summing over i = 1 to N
  cohens_f <- sqrt(ss_between/ss_within)
  
  return(cohens_f)
}


# Compute ratio data based on group_means
# compute_ratio_data <- function(group_means) {
#   # ratio_data <- rep(NA, ngroups)
#   # ratio_data[order(group_means) == 1] <- 1
#   # ratio_data[order(group_means) == 2] <- 2
#   # The choice of the smallest and the second smallest mean makes the scaling 
#   # more robust against changes in the other group means. Since these values 
#   # represent the lower bound of the data, the scale is less sensitive to the 
#   # spread of higher values.
#     
#   # For example:
#   # The value of 2.28 indicates that this particular group mean is 2.28 times 
#   # the scale factor d above the smallest mean. This means the group mean is 
#   # further from the smallest mean compared to the second smallest mean, and 
#   # helps in understanding the relative differences between the group means in 
#   # a normalized manner.
#   
#   # Aantal groepen
#   ngroups <- length(group_means)
#   # Lege vector voor ratio_data
#   ratio_data <- rep(NA, ngroups)
#   # Sorteer de indices van de waarden
#   sorted_indices <- order(group_means)
#   # Wijs 1 en 2 toe aan de kleinste en tweede kleinste waarden
#   ratio_data[sorted_indices[1]] <- 1
#   ratio_data[sorted_indices[2]] <- 2
#   # Bereken d: verschil tussen de tweede kleinste en de kleinste waarde
#   d <- group_means[sorted_indices[2]] - group_means[sorted_indices[1]]
#   # Bereken de ratio's voor de overige waarden
#   for (i in seq_len(ngroups)) {
#     if (!(i %in% sorted_indices[1:2])) {
#       ratio_data[i] <- 1 + (group_means[i] - group_means[sorted_indices[1]]) / d
#     }
#   }
#   return(ratio_data)
# }


generate_scaled_means <- function(group_means, target_f, N, VCOV) {
  if (target_f == 0) {
    # If targeted effect size f is 0, then set all the means to 0.
    new_means <- rep(0, length(group_means))
  } else {
    ratio_vector <- group_means / min(group_means)  # such that ratios remain the same
    #
    objective <- function(d) {
      means_new <- ratio_vector * d
      computed_f <- compute_cohens_f(means_new, N, VCOV)
      return(abs(computed_f - target_f))  # Minimize difference between calculated and desired Cohen's f
    }
    
    opt_result <- optimize(objective, interval = c(0, 100))
    d_optimal <- opt_result$minimum
    
    new_means <- ratio_vector * d_optimal
  }
  
  # Debugging output
  #cat(sprintf("Gevonden d: %.5f voor target f: %.5f\n", d_optimal, target_f))
  #cat("Oude means:", group_means, "\n")
  #cat("Nieuwe means:", new_means_ordered, "\n")
  
  return(new_means)
}

# Compute the population means based on the input parameters
# compute_population_means <- function(pop_es, ratio_pop_means, var_e, ngroups) {
#   means_pop_all <- matrix(NA, ncol = ngroups, nrow = length(pop_es))
#   nr_es <- length(pop_es)
#   for (teller_es in seq_len(nr_es)) {
#     #teller_es = 1
#     
#     # Determine mean values, with ratio of ratio.m
#     # Solve for x here
#     
#     # If all equal, then set population means to all 0
#     if (length(unique(ratio_pop_means)) == 1) {
#       means_pop <- rep(0, ngroups)
#     } else {
#       fun <- function(d) {
#         means_pop = ratio_pop_means * d 
#         (1/sqrt(var_e)) * sqrt((1/ngroups) * sum((means_pop - mean(means_pop))^2)) - pop_es[teller_es] #  AANPASSSEN NAAR NIEUWE FORMULE
#       }
#       d <- uniroot(fun, lower = 0, upper = 100)$root
#       # Construct means_pop
#       means_pop <- ratio_pop_means*d
#     }
#     means_pop_all[teller_es, ] <- means_pop
#   }
#   return(means_pop_all)  
# }

# this function is called from the goric_benchmark_anova() function
parallel_function_means <- function(i, N, var_e, means_pop,
                                    hypos, pref_hypo, comparison, ngroups, sample,
                                    control, form_model_org, mix_weights,
                                    penalty_factor, ...) {

  # Sample residuals
  #epsilon <- rnorm(sum(N), sd = sqrt(var_e/sum(N)))
  # TO DO ws delete:
  #VCOV <- diag(ngroups)
  #diag(VCOV) <- var_e
  # TO DO bovenstaande neemt nu gelijke varianties, wat niet klopt als ongelijke groepsgroottes
  # TO DO deze functie wordt denk ik niet meer gebruikt...
  VCOV <- diag(var_e, ngroups) * N[1]/N
  est <- as.vector(mvtnorm::rmvnorm(n = 1, mean = means_pop, sigma = VCOV))
  names(est) <- names(means_pop)

  # original model formula
  # if (length(form_model_org) > 0) {
  #   model <- form_model_org
  #   lhs <- all.vars(model)[1]
  #   sample[[lhs]] <- as.matrix(sample[, 2:(1 + ngroups)]) %*% matrix(means_pop,
  #                                                                    nrow = ngroups) + epsilon
  #   df_boot <- data.frame(lhs = sample[[lhs]], sample[, 2:(1 + ngroups)])
  #   colnames(df_boot)[1] <- lhs
  #
  #   has_intercept <- attr(terms(model), "intercept") == 1
  #   rhs <- as.character(attr(terms(model), "term.labels"))
  #
  #   # Create the RHS with all other variables and optionally the intercept
  #   if (has_intercept) {
  #     new_rhs <- "."
  #   } else {
  #     new_rhs <- "-1 + ."
  #   }
  #
  #   # Create the new formula
  #   new_model <- as.formula(paste(lhs, "~", new_rhs))
  # } else {
  #   new_model <- y ~ 0 + .
  #   # Generate data
  #   sample$y <- as.matrix(sample[, 2:(1 + ngroups)]) %*% matrix(means_pop,
  #                                                               nrow = ngroups) + epsilon
  #   df_boot <- data.frame(y = sample$y, sample[, 2:(1 + ngroups)])
  # }


  # Obtain fit
  #fit_boot <- lm(new_model, data = df_boot)

  results_goric <- tryCatch(
    {
      # Voer de goric functie uit
      goric(est,
            VCOV = VCOV,
            hypotheses = hypos,
            comparison = comparison,
            type = "gorica", # TO DO goricac ergens ook (als origineel dan goricc ws)?
            control = control,
            mix_weights = mix_weights,
            ...)
    },
    error = function(e) {
      # error message
      message(paste("\nrestriktor ERROR: Error in iteration", i, ":", e$message))
      return(NULL)
    },
    warning = function(w) {
      # warning message
      message(paste("\nrestriktor WARNING: Warning in iteration", i, ":", w$message))
      return(NULL)
    }
  )

  if (is.null(results_goric)) {
    return(NULL)
  }

  # Return the relevant results
  ld_names <- names(results_goric$ratio.gw[pref_hypo, ])
  ld <- results_goric$result$loglik[pref_hypo] - results_goric$result$loglik
  names(ld) <- ld_names

  list(
    #test  = attr(results.goric$objectList[[results.goric$objectNames]]$wt.bar, "mvtnorm"),
    gw  = results_goric$result[pref_hypo, 7], # goric(a) weight
    lw  = results_goric$result$loglik.weights[pref_hypo], # (unpenalized) log-likelihood weight
    rgw = results_goric$ratio.gw[pref_hypo, ], # ratio goric(a) weights
    rlw = results_goric$ratio.lw[pref_hypo, ], # ratio log-likelihood weights
    ld  = ld # loglik difference
  )
}


# model_type = "asymp" ----------------------------------------------------

# this function is called from the benchmark_asymp() function
parallel_function_asymp <- function(i, est, VCOV, hypos, pref_hypo, comparison,
                                    type, control, mix_weights, penalty_factor, ...) {  
  results_goric <- tryCatch(
    {
      # Voer de goric functie uit
      goric(est[i, ], VCOV = VCOV,
            hypotheses = hypos,
            comparison = comparison,
            type = type,
            control = control, 
            mix_weights = mix_weights,
            penalty_factor = penalty_factor,
            ...)
    },
    error = function(e) {
      # error message 
      message(paste("\nrestriktor ERROR: Error in iteration", i, ":", e$message))
      return(NULL)  
    },
    warning = function(w) {
      # warning message
      message(paste("\nrestriktor WARNING: Warning in iteration", i, ":", w$message))
      return(NULL)  
    }
  )
  
  if (is.null(results_goric)) {
    return(NULL)
  }
  
  ld_names <- names(results_goric$ratio.gw[pref_hypo, ])
  ld <- results_goric$result$loglik[pref_hypo] - results_goric$result$loglik
  names(ld) <- ld_names
  
  out <- list(
    gw  = results_goric$result[pref_hypo, 7], # goric(a) weight
    lw  = results_goric$result$loglik.weights[pref_hypo], # (unpenalized) log-likelihood weight
    rgw = results_goric$ratio.gw[pref_hypo, ], # ratio goric(a) weights
    rlw = results_goric$ratio.lw[pref_hypo, ], # ratio log-likelihood weights
    ld  = ld
  )

  return(out)
}


# Define a function to extract and combine values from all elements in each pop_es list
extract_and_combine_values <- function(pop_es_list, value_name) {
  empty_lists_count <- sum(sapply(pop_es_list, is.null))
  #print(empty_lists_count)
  # remove empty lists (is.na)
  #pop_es_list <- pop_es_list[!sapply(pop_es_list, function(x) any(is.na(x)))]
  out <- do.call(rbind, lapply(pop_es_list, function(sub_list) sub_list[[value_name]]))
  attr(out, "empty_lists_count") <- empty_lists_count
  
  return(out)
}


## 
get_results_benchmark <- function(x, object, pref_hypo, pref_hypo_name, 
                                  quant, names_quant, nr.hypos) {
  results <- x
  
    # Use lapply to apply the extract_and_combine_values function to each element in the results list
  gw_combined  <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "gw"))
  lw_combined  <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "lw"))
  rgw_combined <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "rgw"))
  rlw_combined <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "rlw"))
  ld_combined  <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "ld"))
  
  # Calculate CI_benchmarks_gw for each pop_es category
  CI_benchmarks_gw <- lapply(gw_combined, function(gw_values) {
    CI_benchmarks_gw <- matrix(c(object$result[pref_hypo, 7], quantile(gw_values, 
                                                                       quant, na.rm = TRUE)), 
                               nrow = 1)
    colnames(CI_benchmarks_gw) <- names_quant
    rownames(CI_benchmarks_gw) <- pref_hypo_name
    CI_benchmarks_gw

  })

  # Percentile of the observed goric(a) weight within its benchmark
  # distribution, for each pop_es category (0-100 scale, matching the others)
  percentile_gw <- lapply(gw_combined, function(gw_values) {
    Fn <- ecdf(gw_values)
    percentile_gw <- matrix(Fn(object$result[pref_hypo, 7]) * 100, nrow = 1)
    colnames(percentile_gw) <- "percentile"
    rownames(percentile_gw) <- pref_hypo_name
    percentile_gw
  })

  # Same as CI_benchmarks_gw/percentile_gw above, but for the (unpenalized)
  # log-likelihood weight -- used together with percentile_gw by the
  # iter-adequacy check (see goric_percentile_test()/check_iter_adequacy()),
  # since a user may ultimately be interested in either output_type.
  CI_benchmarks_lw <- lapply(lw_combined, function(lw_values) {
    CI_benchmarks_lw <- matrix(c(object$result$loglik.weights[pref_hypo], quantile(lw_values,
                                                                       quant, na.rm = TRUE)),
                               nrow = 1)
    colnames(CI_benchmarks_lw) <- names_quant
    rownames(CI_benchmarks_lw) <- pref_hypo_name
    CI_benchmarks_lw
  })

  percentile_lw <- lapply(lw_combined, function(lw_values) {
    Fn <- ecdf(lw_values)
    percentile_lw <- matrix(Fn(object$result$loglik.weights[pref_hypo]) * 100, nrow = 1)
    colnames(percentile_lw) <- "percentile"
    rownames(percentile_lw) <- pref_hypo_name
    percentile_lw
  })


  # Initialize matrices to store CI benchmarks for current pop_es category
  CI_benchmarks_rgw <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_rlw <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_rlw_ge1 <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_ld <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_ld_ge0 <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  
  # Fill the first column with sample values
  CI_benchmarks_rgw[, 1] <- object$ratio.gw[pref_hypo,] 
  CI_benchmarks_rlw[, 1] <- object$ratio.lw[pref_hypo,] 
  for (j in seq_len(nr.hypos)) {
    if (object$ratio.lw[pref_hypo, j] >= 1) {
      CI_benchmarks_rlw_ge1[j, 1] <- object$ratio.lw[pref_hypo, j] 
    } else {
      CI_benchmarks_rlw_ge1[j, 1] <- 1 / object$ratio.lw[pref_hypo, j] 
    }
  }
  CI_benchmarks_ld[, 1] <- object$result$loglik[pref_hypo] - object$result$loglik 
  CI_benchmarks_ld_ge0[, 1] <- abs(object$result$loglik[pref_hypo] - object$result$loglik) 
  
  CI_benchmarks_rgw_all <- list()
  CI_benchmarks_rlw_all <- list()
  CI_benchmarks_rlw_ge1_all <- list()
  CI_benchmarks_ld_all <- list()
  CI_benchmarks_ld_ge0_all <- list()

  percentile_rgw_all <- list()
  percentile_rlw_all <- list()
  percentile_rlw_ge1_all <- list()
  percentile_ld_all <- list()
  percentile_ld_ge0_all <- list()
  
  # Loop through each pop_es category to fill in the CI benchmark lists
  for (name in names(results)) {
    rgw_combined_values <- rgw_combined[[name]]
    rlw_combined_values  <- rlw_combined[[name]]
    ld_combined_values  <- ld_combined[[name]]
    
    # Prepare rlw_ge1 and ld_ge0 matrices
    rlw_ge1 <- rlw_combined_values
    rlw_ge1[rlw_combined_values < 1] <- 1 / rlw_combined_values[rlw_combined_values < 1]
    ld_ge0 <- abs(ld_combined_values)
    
    # Loop through the hypotheses and calculate the quantiles
    for (j in seq_len(nr.hypos)) {
      CI_benchmarks_rgw[j, 2:(1 + length(quant))] <- quantile(rgw_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_rlw[j, 2:(1 + length(quant))] <- quantile(rlw_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_rlw_ge1[j, 2:(1 + length(quant))] <- quantile(rlw_ge1[, j], quant, na.rm = TRUE)
      CI_benchmarks_ld[j, 2:(1 + length(quant))] <- quantile(ld_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_ld_ge0[j, 2:(1 + length(quant))] <- quantile(ld_ge0[, j], quant, na.rm = TRUE)
    }
    # Loop through the hypotheses and calculate the percentile of the sample finding
    percentile_rgw <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_rlw <- matrix(NA, nrow = nr.hypos, ncol = 1 )
    percentile_rlw_ge1 <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_ld <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_ld_ge0 <- matrix(NA, nrow = nr.hypos, ncol = 1)
    for (j in seq_len(nr.hypos)) {
      Fn <- ecdf(rgw_combined_values[, j])
      percentile_rgw[j, 1] <- Fn(CI_benchmarks_rgw[j, 1]) * 100
      #
      Fn <- ecdf(rlw_combined_values[, j])
      percentile_rlw[j, 1] <- Fn(CI_benchmarks_rlw[j, 1]) * 100
      #
      Fn <- ecdf(rlw_ge1[, j])
      percentile_rlw_ge1[j, 1] <- Fn(CI_benchmarks_rlw_ge1[j, 1]) * 100
      #
      Fn <- ecdf(ld_combined_values[, j])
      percentile_ld[j, 1] <- Fn(CI_benchmarks_ld[j, 1]) * 100
      #
      Fn <- ecdf(ld_ge0[, j])
      percentile_ld_ge0[j, 1] <- Fn(CI_benchmarks_ld_ge0[j, 1]) * 100
    }

    # Label the percentiles (percentage of the benchmark distribution at or
    # below the sample value, i.e. on a 0-100 scale)
    percentile_names <- paste(pref_hypo_name, names(object$ratio.gw[pref_hypo, ]))
    rownames(percentile_rgw) <- rownames(percentile_rlw) <-
      rownames(percentile_rlw_ge1) <- rownames(percentile_ld) <-
      rownames(percentile_ld_ge0) <- percentile_names
    colnames(percentile_rgw) <- colnames(percentile_rlw) <-
      colnames(percentile_rlw_ge1) <- colnames(percentile_ld) <-
      colnames(percentile_ld_ge0) <- "percentile"

    # Store this pop_es category's percentiles so they survive past this
    # loop iteration (mirrors the CI_benchmarks_*_all pattern below)
    percentile_rgw_all[[name]] <- percentile_rgw
    percentile_rlw_all[[name]] <- percentile_rlw
    percentile_rlw_ge1_all[[name]] <- percentile_rlw_ge1
    percentile_ld_all[[name]] <- percentile_ld
    percentile_ld_ge0_all[[name]] <- percentile_ld_ge0
    
    # Set column names for the CI benchmarks
    colnames(CI_benchmarks_rgw) <- colnames(CI_benchmarks_rlw) <- 
      colnames(CI_benchmarks_rlw_ge1) <- colnames(CI_benchmarks_ld) <- 
      colnames(CI_benchmarks_ld_ge0) <- names_quant
    
    # Set row names for the CI benchmarks
    rownames(CI_benchmarks_rgw) <- rownames(CI_benchmarks_rlw) <- 
      rownames(CI_benchmarks_rlw_ge1) <- rownames(CI_benchmarks_ld) <- 
      rownames(CI_benchmarks_ld_ge0) <- paste(pref_hypo_name, names(object$ratio.gw[pref_hypo, ]))
    
    # Store CI benchmarks in lists
    CI_benchmarks_rgw_all[[name]] <- CI_benchmarks_rgw
    CI_benchmarks_rlw_all[[name]] <- CI_benchmarks_rlw
    CI_benchmarks_rlw_ge1_all[[name]] <- CI_benchmarks_rlw_ge1
    CI_benchmarks_ld_all[[name]] <- CI_benchmarks_ld
    CI_benchmarks_ld_ge0_all[[name]] <- CI_benchmarks_ld_ge0
  } 
  
  
  CI_benchmarks_rgw_all_cleaned <- lapply(CI_benchmarks_rgw_all, function(pop_es_list) {
    remove_single_value_rows(pop_es_list, 1)
  })
  
  CI_benchmarks_rlw_all_cleaned <- lapply(CI_benchmarks_rlw_all, function(pop_es_list) {
    remove_single_value_rows(pop_es_list, 1)
  })
  
  CI_benchmarks_rlw_ge1_all_cleaned <- lapply(CI_benchmarks_rlw_ge1_all, function(pop_es_list) {
    remove_single_value_rows(pop_es_list, 1)
  })
  
  CI_benchmarks_ld_all_cleaned <- lapply(CI_benchmarks_ld_all, function(pop_es_list) {
    remove_single_value_rows(pop_es_list, 0)
  })
  
  CI_benchmarks_ld_ge0_all_cleaned <- lapply(CI_benchmarks_ld_ge0_all, function(pop_es_list) {
    remove_single_value_rows(pop_es_list, 0)
  })

  # remove_single_value_rows() drops the preferred hypothesis' self-comparison
  # row from the benchmarks_* matrices above (ratio/difference vs. itself is
  # always 1 or 0). The percentile_*_all matrices were never subject to that
  # filter (a percentile is essentially never exactly 1 or 0), so without this
  # they'd have one row more than their benchmarks_* counterpart, and cbind()-ing
  # them together (e.g. in print/summary) would silently misalign rows instead
  # of erroring. Subset by rowname to keep both sets of matrices in lockstep.
  align_rows <- function(percentile_list, cleaned_list) {
    Map(function(perc, bench) perc[rownames(bench), , drop = FALSE],
        percentile_list, cleaned_list)
  }
  percentile_rgw_all_cleaned <- align_rows(percentile_rgw_all, CI_benchmarks_rgw_all_cleaned)
  percentile_rlw_all_cleaned <- align_rows(percentile_rlw_all, CI_benchmarks_rlw_all_cleaned)
  percentile_rlw_ge1_all_cleaned <- align_rows(percentile_rlw_ge1_all, CI_benchmarks_rlw_ge1_all_cleaned)
  percentile_ld_all_cleaned <- align_rows(percentile_ld_all, CI_benchmarks_ld_all_cleaned)
  percentile_ld_ge0_all_cleaned <- align_rows(percentile_ld_ge0_all, CI_benchmarks_ld_ge0_all_cleaned)


  rgw_combined <- lapply(rgw_combined, function(pop_es_list) {
    remove_single_value_col(pop_es_list, 1)
  })
  
  rlw_combined <- lapply(rlw_combined, function(pop_es_list) {
    remove_single_value_col(pop_es_list, 1)
  })
  
  ld_combined <- lapply(ld_combined, function(pop_es_list) {
    remove_single_value_col(pop_es_list, 0)
  })
  
  
  OUT <- list(
    benchmarks_gw = CI_benchmarks_gw,
    benchmarks_lw = CI_benchmarks_lw,
    benchmarks_rgw = CI_benchmarks_rgw_all_cleaned,
    benchmarks_rlw = CI_benchmarks_rlw_all_cleaned,
    benchmarks_rlw_ge1 = CI_benchmarks_rlw_ge1_all_cleaned,
    benchmarks_difLL = CI_benchmarks_ld_all_cleaned,
    benchmarks_absdifLL = CI_benchmarks_ld_ge0_all_cleaned,
    percentile_gw  = percentile_gw,
    percentile_lw  = percentile_lw,
    percentile_rgw = percentile_rgw_all_cleaned,
    percentile_rlw = percentile_rlw_all_cleaned,
    percentile_rlw_ge1 = percentile_rlw_ge1_all_cleaned,
    percentile_difLL = percentile_ld_all_cleaned,
    percentile_absdifLL = percentile_ld_ge0_all_cleaned,
    combined_values = list(gw_combined = gw_combined,
                           lw_combined = lw_combined,
                           rgw_combined = rgw_combined,
                           rlw_combined = rlw_combined,
                           ld_combined = ld_combined)
  )

  return(OUT)
}


# Calculate the error probability
calculate_error_probability <- function(object, hypos, pref_hypo, est, 
                                        VCOV, control, ...) {
  # Error probability based on complement of preferred hypothesis in data
  nr_hypos <- dim(object$result)[1]
  if (nr_hypos == 2 && object$comparison == "complement") { 
    if (object$type == 'goric') {
      # TO DO also here re-run with GORICA, as we do for sample value as well?
      #       is ws al opgelost als we goric en gorica resultaten gelijk maken!!!
      #       Dus dan laten staan + re-run met gorica niet nodig dan ook!
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
        results_goric_pref <- goric(est, VCOV = VCOV,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = "gorica", 
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
  return(error_prob)
}


# Diagnostic check: is 'iter' large enough for a stable benchmark?
#
# Under the "Observed" population (only present when pop_es/pop_est was left
# NULL, so it is auto-added), the benchmark draws are simulated centered
# exactly on the observed estimate. So the observed goric(a) weight is
# expected to sit close to the 50th percentile of its own benchmark
# distribution; systematic deviation from that is a sign that 'iter' bootstrap
# draws is not (yet) enough to have converged, rather than a real finding.
#
# This is checked using 'gw' rather than 'rgw'/'rlw' on purpose: gw is a
# single, bounded [0, 1] weight per draw, whereas rgw/rlw = gw_pref / gw_k is
# a ratio of two such weights (see calc_ICweights() in
# goric_calculate_IC_weights.R). When gw_k is small, that ratio amplifies
# small (essentially unavoidable, e.g. optimizer-precision-level)
# fluctuations in gw_k multiplicatively, so rgw/rlw are a noisier basis for
# this particular check than gw itself -- even though all of them should, in
# principle, lead to the same conclusion about whether 'iter' is adequate.
check_iter_adequacy <- function(benchmark_results, observed_name, iter,
                                band = c(0.495, 0.505), control = list(), ...) {
  if (!observed_name %in% names(benchmark_results$percentile_gw)) {
    # No "Observed" population in this run (user supplied a custom pop_es/
    # pop_est without an "Observed" category) -- nothing to check.
    return(invisible(NULL))
  }

  sample_gw <- benchmark_results$benchmarks_gw[[observed_name]][1, 1]
  sample_lw <- benchmark_results$benchmarks_lw[[observed_name]][1, 1]
  gw_draws  <- benchmark_results$combined_values$gw_combined[[observed_name]]
  lw_draws  <- benchmark_results$combined_values$lw_combined[[observed_name]]

  chk_gw <- goric_percentile_test(gw_draws, sample_gw, band = band, control = control, ...)
  chk_lw <- goric_percentile_test(lw_draws, sample_lw, band = band, control = control, ...)

  if (!(chk_gw$converged && chk_lw$converged)) {
    # Same style/format as the auto-'iter' messages in run_benchmark_simulation()
    # (percentile + support per output_type, band mentioned once, a closing
    # suggestion) -- this is the fixed-'iter' counterpart of those: since
    # 'iter' was user-specified here, there is no growing/stopping decision
    # to report, just whether the fixed 'iter' the user chose looks adequate.
    # "close/equal" (rather than just "close") for a converged output_type,
    # since printing e.g. percentile = 50.0 right next to "close to" reads
    # oddly when it is, in fact, (numerically) equal.
    describe_type <- function(chk) {
      if (chk$converged) {
        "-- close/equal to the 50th percentile."
      } else {
        "-- not close to the 50th percentile."
      }
    }
    # The closing suggestion is tailored to which output_type(s) are actually
    # inadequate at this 'iter' -- mirroring the "Notably, you could
    # user-specify a lower 'iter' ..." clause in the auto-'iter' converged
    # message, but for the opposite situation: here, if only one of 'gw'/'lw'
    # is the problem, raising 'iter' only matters to someone who cares about
    # that one; someone only interested in the other (already-adequate) one
    # doesn't need to.
    closing <- if (!chk_gw$converged && !chk_lw$converged) {
      "Consider increasing 'iter' for a more stable benchmark."
    } else if (!chk_gw$converged) {
      paste0(
        "Consider increasing 'iter' for a more stable benchmark when you are (also) interested ",
        "in GORIC(A) weights ('gw') and not (only) in log-likelihood weights ('lw')."
      )
    } else {
      paste0(
        "Consider increasing 'iter' for a more stable benchmark when you are (also) interested ",
        "in log-likelihood weights ('lw') and not (only) in GORIC(A) weights ('gw')."
      )
    }
    message(
      "\nrestriktor Message: For the user-specified 'iter' = ", iter, ", the value based on ",
      "your data (called the 'Sample' value in the output) is not close to the 50th ",
      "percentile of the benchmark distribution under the 'Observed' population.\n",
      "This is checked for the weight-type output, where the GORICA determined the 'support', ",
      "that is, the GORICA weight, for the percentile being between ", band[1], " and ", band[2],
      " versus outside of that range:\n",
      "For output_type = 'gw', the percentile is ", sprintf("%.1f", chk_gw$percentile),
      " (support = ", sprintf("%.3f", chk_gw$gw), ") ", describe_type(chk_gw), "\n",
      "For output_type = 'lw', the percentile is ", sprintf("%.1f", chk_lw$percentile),
      " (support = ", sprintf("%.3f", chk_lw$gw), ") ", describe_type(chk_lw), "\n",
      closing
    )
  }
  invisible(NULL)
}


# Core of the Monte Carlo adequacy check: is 'sample_value' consistent with
# being the median of 'draws'? Rather than a binomial test against p = 0.5,
# this expresses "close to the median" as a small, fixed interval around
# 0.5 (band -- a region of practical equivalence, deliberately NOT scaled to
# the percentile estimate's own standard error, since doing so would just
# reconstruct a Wald interval that goric(a) would then re-analyze with that
# same standard error a second time) and lets goric(a) itself -- type =
# "gorica", using the percentile estimate's own sampling variance as VCOV --
# weigh the evidence for "the true percentile lies in that band" against its
# complement. This keeps the whole diagnostic inside goric(a)'s own
# machinery rather than a separate hypothesis-testing framework. Returns the
# percentile (0-100) 'sample_value' falls at within 'draws', and the
# resulting gorica(a) weight for the "in-band" hypothesis (>= 0.5 is taken
# to mean converged, i.e. goric(a) favors "close to the median" over its
# complement).
goric_percentile_test <- function(draws, sample_value, band = c(0.495, 0.505),
                                  control = list(), ...) {
  draws <- draws[!is.na(draws)]
  n <- length(draws)
  k <- sum(draws <= sample_value)
  phat <- k / n

  # Continuity-corrected proportion, used for the variance only, so that
  # phat = 0 or 1 (possible at small iter) doesn't collapse VCOV to exactly
  # zero.
  phat_v <- (k + 0.5) / (n + 1)
  VCOV <- matrix(phat_v * (1 - phat_v) / n)
  rownames(VCOV) <- colnames(VCOV) <- "p"
  est <- c(p = phat)

  H1 <- paste(band[1], "< p <", band[2])
  fit <- goric(est, VCOV = VCOV, hypotheses = list(H1 = H1),
              comparison = "complement", type = "gorica",
              control = control, ...)
  gw_H1 <- fit$result$gorica.weights[fit$result$model == "H1"]

  list(percentile = 100 * phat, n = n, band = band,
      gw = gw_H1, converged = isTRUE(gw_H1 >= 0.5))
}


# Run the pop_es/pop_est simulation loop used by benchmark_means()/
# benchmark_asymp(), growing the number of draws adaptively when the user
# leaves 'iter' unspecified (iter = NULL): start at iter_min draws and, if
# the "Observed" population's sample value is not close to its own 50th
# percentile for output_type = 'gw' and/or 'lw' (see goric_percentile_test()
# above for how "close" is defined and tested), add iter_step more draws --
# WITHOUT discarding or redrawing the ones already computed -- repeating until
# either it looks adequate or iter_max is reached. If the user supplies a
# fixed numeric 'iter' instead, this runs exactly one round of that many
# draws (the original, non-adaptive behaviour); benchmark_means()/
# benchmark_asymp() then call check_iter_adequacy() themselves afterwards
# for that fixed-iter case (this function does not, to avoid messaging
# twice).
#
# Returns list(parallel_function_results = <as before, one element per
# pop_es/pop_est category>, iter = <final number of draws used>).

run_benchmark_simulation <- function(nr_es, rnames, name_prefix, center_matrix,
                                     colnames_vec, VCOV, hypos, pref_hypo,
                                     comparison, control, mix_weights,
                                     penalty_factor, Heq, object, iter,
                                     es_labels = rnames,
                                     iter_min = 500, iter_step = 100,
                                     iter_max = 2000, band = c(0.495, 0.505),
                                     stability_tol = 1, ...) {

  auto_iter <- is.null(iter)
  sample_gw <- object$result[pref_hypo, 7]
  sample_lw <- object$result$loglik.weights[pref_hypo]
  obs_pos   <- which(rnames == "Observed")

  parallel_function_results <- vector("list", nr_es)
  names(parallel_function_results) <- paste0(name_prefix, rnames)
  est_accum <- vector("list", nr_es)

  n_done <- 0L
  target <- if (auto_iter) min(iter_min, iter_max) else iter

  # First (smallest) number of draws at which each output_type individually
  # already showed adequate support, tracked separately since 'gw' and 'lw'
  # need not converge at the same iter -- reported in the adequacy message so
  # the user can see whether a smaller 'iter' would have sufficed had they
  # only cared about one of the two.
  min_iter_gw <- NA_integer_
  min_iter_lw <- NA_integer_

  # Previous round's percentile per output_type, used to detect when a
  # not-yet-converged percentile has stopped moving (see 'no_progress'
  # below) -- NA until a second round exists to compare against.
  prev_percentile_gw <- NA_real_
  prev_percentile_lw <- NA_real_

  progressr::handlers(progressr::handler_txtprogressbar(char = ">"))

  repeat {
    n_new <- target - n_done

    progressr::with_progress({
      p <- progressr::progressor(along = seq_len(n_new * nr_es))

      for (teller_es in seq_len(nr_es)) {
        cat("Calculating benchmark for", es_labels[teller_es],
            "-- draws", n_done + 1, "to", target, "\n")

        new_draws <- mvtnorm::rmvnorm(n = n_new, center_matrix[teller_es, ], sigma = VCOV)
        colnames(new_draws) <- colnames_vec
        est_accum[[teller_es]] <- rbind(est_accum[[teller_es]], new_draws)
        est_full <- est_accum[[teller_es]]

        # Wrapper function for future_lapply
        wrapper_function_asymp <- function(i) {
          p() # Update progress
          parallel_function_asymp(i,
                                  est = est_full, VCOV = VCOV,
                                  hypos = hypos, pref_hypo = pref_hypo,
                                  comparison = comparison, type = "gorica",
                                  control = control, mix_weights = mix_weights,
                                  penalty_factor = penalty_factor,
                                  Heq = Heq, ...)
        }

        new_results <- future_lapply(
          seq(n_done + 1, target),
          wrapper_function_asymp,
          future.seed = TRUE # Ensures safe and reproducible random number generation
        )

        key <- paste0(name_prefix, rnames[teller_es])
        parallel_function_results[[key]] <- c(parallel_function_results[[key]], new_results)
      }
    })

    n_done <- target

    if (!auto_iter) break # fixed iter: exactly one round, done

    if (length(obs_pos) == 1) {
      obs_key <- paste0(name_prefix, "Observed")
      gw_draws <- vapply(parallel_function_results[[obs_key]], function(x) {
        if (is.null(x)) NA_real_ else as.numeric(x$gw)
      }, numeric(1))
      lw_draws <- vapply(parallel_function_results[[obs_key]], function(x) {
        if (is.null(x)) NA_real_ else as.numeric(x$lw)
      }, numeric(1))
      gw_draws <- gw_draws[!is.na(gw_draws)]
      lw_draws <- lw_draws[!is.na(lw_draws)]
      chk_gw <- goric_percentile_test(gw_draws, sample_gw, band = band, control = control)
      chk_lw <- goric_percentile_test(lw_draws, sample_lw, band = band, control = control)
      converged <- chk_gw$converged && chk_lw$converged
      if (is.na(min_iter_gw) && chk_gw$converged) min_iter_gw <- n_done
      if (is.na(min_iter_lw) && chk_lw$converged) min_iter_lw <- n_done

      # Has a not-yet-converged percentile stopped moving with more draws?
      # If every currently-unconverged output_type is stable, more draws
      # are unlikely to bring it closer to the 50th percentile either --
      # that looks like a genuine feature of the benchmark distribution
      # rather than a Monte Carlo adequacy problem, so growth can stop early
      # instead of running all the way to iter_max. FALSE (not NA) on the
      # very first round, since there is no
      # earlier percentile yet to compare against. This is a first-cut
      # heuristic based on a single round-to-round comparison -- if it
      # turns out to call "stalled" too eagerly (a percentile can drift by
      # a percentage point or two by chance alone), requiring stability
      # across two consecutive rounds instead of one would be the natural
      # tightening.
      stable_gw <- isTRUE(!is.na(prev_percentile_gw) &&
                            abs(chk_gw$percentile - prev_percentile_gw) < stability_tol)
      stable_lw <- isTRUE(!is.na(prev_percentile_lw) &&
                            abs(chk_lw$percentile - prev_percentile_lw) < stability_tol)
      no_progress <- (chk_gw$converged || stable_gw) && (chk_lw$converged || stable_lw)
      prev_percentile_gw <- chk_gw$percentile
      prev_percentile_lw <- chk_lw$percentile
    } else {
      # No "Observed" category (custom pop_es/pop_est) -- nothing to check
      # against, so stop growing after the first (iter_min-sized) batch.
      converged <- TRUE
      no_progress <- FALSE
      chk_gw <- chk_lw <- NULL
    }

    if (converged || (no_progress && !converged) || n_done >= iter_max) {
      if (auto_iter && !is.null(chk_gw)) {
        if (converged) {
          message(
            "\nrestriktor Message: 'iter' was not specified, so it was set automatically.\n",
            "Using iter = ", n_done, " draws, the value based on your data (called the ",
            "'Sample' value in the output) is close to the 50th percentile of the benchmark ",
            "distribution under the 'Observed' population.\n",
            "This is checked for the weight-type output, where the GORICA determined the ",
            "'support', that is, the GORICA weight, for the percentile being between ",
            band[1], " and ", band[2], " versus outside of that range:\n",
            "For output_type = 'gw', the percentile is ", sprintf("%.1f", chk_gw$percentile),
            " (support = ", sprintf("%.3f", chk_gw$gw), "; minimum iter = ", min_iter_gw,
            "); for output_type = 'lw', the percentile is ", sprintf("%.1f", chk_lw$percentile),
            " (support = ", sprintf("%.3f", chk_lw$gw), "; minimum iter = ", min_iter_lw, ").\n",
            "So, no further draws were added.",
            if (min_iter_gw != min_iter_lw) paste0(
              " Notably, you could user-specify a lower 'iter' if you are not interested in ",
              "both GORIC(A) weights ('gw') and log-likelihood weights ('lw')."
            ) else ""
          )
        } else if (no_progress) {
          describe_type <- function(chk) {
            if (chk$converged) {
              "-- converged (close to the 50th percentile)."
            } else {
              "-- stabilized, but not close to the 50th percentile."
            }
          }
          message(
            "\nrestriktor Message: The default for 'iter' was used, since it was not ",
            "specified. Draws were added up to iter = ", n_done, " (below the maximum of ",
            iter_max, "). Over the last ", iter_step, " draws, the percentile(s) changed by ",
            "less than ", stability_tol, " percentage point(s). This may imply that increasing ",
            "'iter' to ", iter_max, " will not help.\n",
            "The percentile values are not near 50 yet. This is checked for the weight-type ",
            "output, where the GORICA determined the 'support', that is, the GORICA weight, ",
            "for the percentile being between ", band[1], " and ", band[2],
            " versus outside of that range:\n",
            "For output_type = 'gw', the percentile is ", sprintf("%.1f", chk_gw$percentile),
            " (support = ", sprintf("%.3f", chk_gw$gw), ") ", describe_type(chk_gw), "\n",
            "For output_type = 'lw', the percentile is ", sprintf("%.1f", chk_lw$percentile),
            " (support = ", sprintf("%.3f", chk_lw$gw), ") ", describe_type(chk_lw), "\n",
            "Since more draws are unlikely to change this, consider inspecting the 'Observed' ",
            "benchmark distribution directly rather than increasing 'iter' further."
          )
        } else {
          message(
            "\nrestriktor Message: 'iter' was not specified, so it was increased ",
            "automatically up to its maximum of iter = ", iter_max, " draws. The 'Observed' ",
            "population's sample value is still not close to the 50th percentile of its ",
            "benchmark distribution for output_type = 'gw' (percentile is ",
            sprintf("%.1f", chk_gw$percentile), ", support = ", sprintf("%.3f", chk_gw$gw),
            ") and/or output_type = 'lw' (percentile is ", sprintf("%.1f", chk_lw$percentile),
            ", support = ", sprintf("%.3f", chk_lw$gw), "). Consider re-running with a ",
            "manually specified, larger 'iter' for a more stable benchmark."
          )
        }
      }
      break
    }

    target <- min(n_done + iter_step, iter_max)
  }

  list(parallel_function_results = parallel_function_results, iter = n_done)
}


calculate_hypothesis_rate <- function(x, q = 1) {
  return(colMeans(x > q))
}


# called by the benchmark.print() function
print_section <- function(header, content_printer, nchar, text_color, reset) {
  cat("\n")
  cat(strrep("=", nchar), "\n")
  cat(paste0(text_color, header, reset), "\n")
  cat(strrep("-", nchar), "\n")
  content_printer()
  #cat(strrep("-", nchar), "\n")
  cat("\n")
}


format_value <- function(value) {
  if (is.na(value)) {
    return("")  
  }
  if (abs(value) >= 1000 || (abs(value) <= 0.001 && value != 0)) {
    return(sprintf("%.3e", value))  
  } else {
    return(sprintf("%.3f", value))  
  }
}


# called by the benchmark.print() function
print_rounded_es_value <- function(df, pop_es, model_type, text_color, reset) {
  if (model_type == "benchmark_asymp") {
    pop_es_value <- gsub("pop_est = ", "", pop_es)
    cat(sprintf("Population estimates = %s%s%s\n", text_color, pop_es_value, reset))
  } else {
    pop_es_value <- gsub("pop_es = ", "", pop_es)
    cat(sprintf("Population effect-size = %s%s%s\n", text_color, pop_es_value, reset))
  }
  
  #formatted_column <- sprintf("%.3f", df)
  formatted_values <- sapply(as.numeric(df), format_value)
  formatted_df <- `dim<-`(formatted_values, dim(df))
  rownames(formatted_df) <- rownames(df)
  colnames(formatted_df) <- colnames(df)
  print(formatted_df, row.names = TRUE, quote = FALSE)
  cat("\n")
}


print_formatted_matrix <- function(mat, text_color, reset) {
  row_width <- max(nchar(rownames(mat))) + 0  
  col_widths <- apply(mat, 2, function(col) max(nchar(as.character(col))))
  col_widths <- pmax(col_widths, nchar(colnames(mat))) + 2  
  
  cat("Population Estimates (PE):\n")
  #cat(sprintf("%-10s", ""))
  cat(sprintf(paste0("%-", row_width, "s"), "")) 
  for (k in 1:ncol(mat)) {
    cat(sprintf(paste0("%", col_widths[k], "s"), colnames(mat)[k]))
  }
  cat("\n")
  
  # Print rows
  for (i in 1:nrow(mat)) {
    cat(sprintf(paste0("%-", row_width, "s"), rownames(mat)[i]))
    for (j in 1:ncol(mat)) {
      cat(sprintf(paste0("%s%", col_widths[j], "s%s"), text_color, mat[i, j], reset))
    }
    cat("\n")
  }
}

# called by the benchmark_asymp() function
check_rhs_constants <- function(rhs_list) {
  constants_check <- lapply(rhs_list, function(element) {
    return(any(element != 0))
  })
  hypotheses_with_constants <- names(constants_check)[unlist(constants_check)]
  if (length(hypotheses_with_constants) > 0) {
    warning_message <- paste0("\nrestriktor WARNING: The following hypotheses contain constants",
                              " greater or less than 0: ", 
                              paste(hypotheses_with_constants, collapse = ", "),
                              ". The default population estimates are likely incorrect.",
                              " Consider providing custom population estimates via the",
                              " pop_est argument.")
    warning(warning_message, call. = FALSE)
  }
}

# restricted least squares
theta_restricted <- function(theta, V, R, rhs) {
  correction <- V %*% t(R) %*% solve(R %*% V %*% t(R)) %*% (R %*% theta - rhs)
  as.vector(theta - correction)
}
