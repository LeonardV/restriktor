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


# Adjust the variance based on alternative group sizes
adjust_variance <- function(var_e, N, alt_group_size, ngroups) {
  # Possibly adjust var_e based on other sample size
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
  return(list(var_e = var_e, N = N))
}

generate_scaled_means <- function(group_means, target_f, N, VCOV) {
  # TO DO evt delete:
  #original_order <- order(group_means)
  #sorted_means <- sort(group_means)
  #min_value <- sorted_means[1]
  #ratio_vector <- sorted_means / min_value  # Behoudt de verhoudingen
  ##
  #objective <- function(d) {
  #  means_new <- ratio_vector * d
  #  computed_f <- compute_cohens_f(means_new, N, VCOV)
  #  # TO DO nu staan de means in andere volgorde en hoort het dus niet meer bij die VCOV (als ongelijke n's)
  #  #       waarom volgorde eigenlijk aanpassen?
  #  return(abs(computed_f - target_f))  # Minimaliseer het verschil tussen berekende en gewenste Cohen's f
  #}
  ##
  #opt_result <- optimize(objective, interval = c(0, 100))
  #d_optimal <- opt_result$minimum
  ##
  #new_means <- ratio_vector * d_optimal
  #new_means_ordered <- new_means[original_order]  # Behoud de oorspronkelijke volgorde
  ##
  ## Debugging output
  ##cat(sprintf("Gevonden d: %.5f voor target f: %.5f\n", d_optimal, target_f))
  ##cat("Oude means:", group_means, "\n")
  ##cat("Nieuwe means:", new_means_ordered, "\n")
  ##
  ##return(new_means_ordered)
  ##
  ##
  # TO DO Suggestion RMK (works much better :-))
  if (target_f == 0) {
    # If targeted effect size f is 0, then set all the means to 0.
    new_means <- rep(0, length(group_means))
  } else {
    ratio_vector <- group_means / min(group_means)  # Behoudt de verhoudingen
    #
    objective <- function(d) {
      means_new <- ratio_vector * d
      computed_f <- compute_cohens_f(means_new, N, VCOV)
      return(abs(computed_f - target_f))  # Minimaliseer het verschil tussen berekende en gewenste Cohen's f
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
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)
    },
    warning = function(w) {
      # warning message
      message(paste("Warning in iteration", i, ":", w$message))
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
      message(paste("Error in iteration", i, ":", e$message))
      return(NULL)  
    },
    warning = function(w) {
      # warning message
      message(paste("Warning in iteration", i, ":", w$message))
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
    benchmarks_rgw = CI_benchmarks_rgw_all_cleaned,
    benchmarks_rlw = CI_benchmarks_rlw_all_cleaned,
    benchmarks_rlw_ge1 = CI_benchmarks_rlw_ge1_all_cleaned,
    benchmarks_difLL = CI_benchmarks_ld_all_cleaned,
    benchmarks_absdifLL = CI_benchmarks_ld_ge0_all_cleaned,
    combined_values = list(gw_combined = gw_combined, 
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
    warning_message <- paste0("Restriktor Warning: The following hypotheses contain constants",
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
