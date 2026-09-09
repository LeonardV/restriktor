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

# Overlapping coefficient (OVL) between two samples' distributions: the
# proportion of area shared by their kernel density estimates -- 1 = fully
# overlapping (identical) distributions, 0 = no overlap at all. This is the
# numeric counterpart of what overlaying, say, the 'Observed' and 'No-effect'
# population's benchmark draws in a density/histogram plot shows visually.
# Both densities are evaluated on the SAME x-grid (spanning the combined
# range of both samples, with a small margin so neither tail gets clipped)
# so their y-values are directly comparable point-by-point; Gaussian kernel
# with bw = "nrd0" matches the density estimate already used elsewhere in
# this file (see calculate_power() above). Like a plain density()-based plot,
# this does not respect the [0, 1] bounds that 'gw'/'lw' are constrained to
# -- a Gaussian kernel can put a little mass just outside that range, same
# as it would visually spill past the axis limits in a density plot.
compute_overlap <- function(draws1, draws2, n = 512) {
  draws1 <- draws1[!is.na(draws1)]
  draws2 <- draws2[!is.na(draws2)]
  if (length(draws1) < 2 || length(draws2) < 2 ||
      sd(draws1) == 0 || sd(draws2) == 0) {
    # Can't fit a (non-degenerate) density on too few draws, or on draws
    # that are all identical (zero variance) -- not enough information for
    # an overlap coefficient to be meaningful.
    return(NA_real_)
  }
  rng <- range(c(draws1, draws2))
  pad <- diff(rng) * 0.1
  from <- rng[1] - pad
  to   <- rng[2] + pad
  d1 <- tryCatch(
    density(draws1, kernel = "gaussian", bw = "nrd0", from = from, to = to, n = n),
    error = function(e) NULL
  )
  d2 <- tryCatch(
    density(draws2, kernel = "gaussian", bw = "nrd0", from = from, to = to, n = n),
    error = function(e) NULL
  )
  if (is.null(d1) || is.null(d2)) return(NA_real_)
  dx <- mean(diff(d1$x))
  min(sum(pmin(d1$y, d2$y)) * dx, 1) # cap at 1 for the (rare) numerical-integration overshoot
}

# Which population the "Overlap with ..." column is computed against.
# Prefers the "Observed" category when one is present (the usual case: the
# default 'No-effect' plus the user's actually-observed estimates). When
# there is no "Observed" category -- a custom pop_es/pop_est run, which
# never includes one -- falls back to the LAST population in
# 'pop_names' (i.e. the last pop_es/pop_est value supplied), so overlap is
# still shown rather than silently omitted. Returns NULL only when there are
# no populations at all.
determine_overlap_reference <- function(pop_names) {
  observed_name <- pop_names[grepl("= Observed$", pop_names)]
  if (length(observed_name) == 1) {
    return(observed_name)
  }
  if (length(pop_names) == 0) {
    return(NULL)
  }
  pop_names[length(pop_names)]
}

# Overlap of the reference population's (see determine_overlap_reference())
# benchmark distribution against EACH other population present -- by default
# just 'No-effect', but there can be more than one if the user supplied
# multiple pop_es/pop_est values. 'draws_combined' is a named list of draw
# vectors keyed by population (e.g. gw_combined/lw_combined from
# get_results_benchmark()), with names like "pop_es = No-effect"/
# "pop_es = Observed" (or "pop_est = ..." for benchmark_asymp). Returns a
# named numeric vector (named by the other population, plus the reference
# population's own entry fixed at 1 -- its overlap with itself), or NULL if
# there are no populations at all.
compute_overlap_vs_observed <- function(draws_combined) {
  reference_name <- determine_overlap_reference(names(draws_combined))
  if (is.null(reference_name)) {
    return(NULL)
  }
  reference_draws <- draws_combined[[reference_name]]
  other_names <- setdiff(names(draws_combined), reference_name)
  overlaps <- vapply(other_names, function(nm) {
    compute_overlap(reference_draws, draws_combined[[nm]])
  }, numeric(1))
  names(overlaps) <- other_names
  overlaps[reference_name] <- 1
  overlaps
}

# Same idea as compute_overlap_vs_observed(), but for the matrix-valued
# per-draw statistics ('rgw'/'rlw'/'ld'): each population's entry is a
# matrix with one column per alternative hypothesis being compared against
# the preferred one (e.g. column "H2"), rather than a single vector. Overlap
# is computed per shared column, so the result is a list (named by the other
# population, plus the reference population's own entry -- a same-length
# vector of 1s, one per column) of named numeric vectors (named by
# hypothesis column). Returns NULL under the same conditions as
# compute_overlap_vs_observed() (no populations at all), or when there are
# no columns to compare -- e.g. after remove_single_value_col() has dropped
# every column because there is only one alternative hypothesis and it is
# constant across draws).
compute_overlap_vs_observed_matrix <- function(draws_combined) {
  reference_name <- determine_overlap_reference(names(draws_combined))
  if (is.null(reference_name)) {
    return(NULL)
  }
  reference_mat <- draws_combined[[reference_name]]
  other_names <- setdiff(names(draws_combined), reference_name)
  cols <- colnames(reference_mat)
  if (is.null(cols) || length(cols) == 0) {
    return(NULL)
  }
  overlaps <- lapply(other_names, function(nm) {
    other_mat <- draws_combined[[nm]]
    shared_cols <- intersect(cols, colnames(other_mat))
    vapply(shared_cols, function(cn) {
      compute_overlap(reference_mat[, cn], other_mat[, cn])
    }, numeric(1))
  })
  names(overlaps) <- other_names
  self_overlap <- rep(1, length(cols))
  names(self_overlap) <- cols
  overlaps[[reference_name]] <- self_overlap
  overlaps
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

  # Which population the "median ref. pop." columns below (and the "Overlap
  # with ..." column -- see overlap_reference_pop further down, which now
  # just reuses this) are computed against: "Observed" when present, else
  # the last pop_es/pop_est population supplied. Determined once here, since
  # it doesn't depend on gw/lw/rgw/etc. -- see determine_overlap_reference().
  reference_pop_name <- determine_overlap_reference(names(gw_combined))

  # Calculate CI_benchmarks_gw for each pop_es category
  CI_benchmarks_gw <- lapply(gw_combined, function(gw_values) {
    CI_benchmarks_gw <- matrix(c(object$result[pref_hypo, 7], quantile(gw_values, 
                                                                       quant, na.rm = TRUE)), 
                               nrow = 1)
    colnames(CI_benchmarks_gw) <- names_quant
    rownames(CI_benchmarks_gw) <- pref_hypo_name
    CI_benchmarks_gw

  })

  # Percentile of the observed goric(a) weight ('Sample' value) within its
  # own benchmark distribution, for each pop_es category (0-100 scale,
  # matching the others). Column/print header is "Pctl. Sample".
  pctl_Sample_gw <- lapply(gw_combined, function(gw_values) {
    Fn <- ecdf(gw_values)
    pctl_Sample_gw <- matrix(Fn(object$result[pref_hypo, 7]) * 100, nrow = 1)
    colnames(pctl_Sample_gw) <- "Pctl. Sample"
    rownames(pctl_Sample_gw) <- pref_hypo_name
    pctl_Sample_gw
  })

  # Percentile of the reference population's own median goric(a) weight,
  # located within *each* population's own benchmark distribution (0-100
  # scale): where the reference/last population's typical draw falls
  # relative to this population's distribution. For the reference
  # population's own row this is hardcoded to exactly 50 -- ecdf() of a
  # finite/even-n sample evaluated at its own median need not land on
  # exactly 50. Column/print header is "Pctl. median ref. pop.".
  median_gw_ref <- median(gw_combined[[reference_pop_name]], na.rm = TRUE)
  pctl_medianRefPop_gw <- lapply(names(gw_combined), function(nm) {
    value <- if (identical(nm, reference_pop_name)) {
      50
    } else {
      Fn <- ecdf(gw_combined[[nm]])
      Fn(median_gw_ref) * 100
    }
    out <- matrix(value, nrow = 1)
    colnames(out) <- "Pctl. median ref. pop."
    rownames(out) <- pref_hypo_name
    out
  })
  names(pctl_medianRefPop_gw) <- names(gw_combined)

  # Same as CI_benchmarks_gw/pctl_Sample_gw above, but for the (unpenalized)
  # log-likelihood weight -- used together with pctl_Sample_gw by the
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

  pctl_Sample_lw <- lapply(lw_combined, function(lw_values) {
    Fn <- ecdf(lw_values)
    pctl_Sample_lw <- matrix(Fn(object$result$loglik.weights[pref_hypo]) * 100, nrow = 1)
    colnames(pctl_Sample_lw) <- "Pctl. Sample"
    rownames(pctl_Sample_lw) <- pref_hypo_name
    pctl_Sample_lw
  })

  median_lw_ref <- median(lw_combined[[reference_pop_name]], na.rm = TRUE)
  pctl_medianRefPop_lw <- lapply(names(lw_combined), function(nm) {
    value <- if (identical(nm, reference_pop_name)) {
      50
    } else {
      Fn <- ecdf(lw_combined[[nm]])
      Fn(median_lw_ref) * 100
    }
    out <- matrix(value, nrow = 1)
    colnames(out) <- "Pctl. median ref. pop."
    rownames(out) <- pref_hypo_name
    out
  })
  names(pctl_medianRefPop_lw) <- names(lw_combined)


  # Initialize matrices to store CI benchmarks for current pop_es category
  CI_benchmarks_rgw <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_rlw <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_rlw_ge1 <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_rgw_log <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_rlw_log <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_ld <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_ld_ge0 <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))

  # Fill the first column with sample values
  CI_benchmarks_rgw[, 1] <- object$ratio.gw[pref_hypo,]
  CI_benchmarks_rlw[, 1] <- object$ratio.lw[pref_hypo,]
  CI_benchmarks_rgw_log[, 1] <- log(object$ratio.gw[pref_hypo,])
  CI_benchmarks_rlw_log[, 1] <- log(object$ratio.lw[pref_hypo,])
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
  CI_benchmarks_rgw_log_all <- list()
  CI_benchmarks_rlw_log_all <- list()
  CI_benchmarks_ld_all <- list()
  CI_benchmarks_ld_ge0_all <- list()

  pctl_Sample_rgw_all <- list()
  pctl_Sample_rlw_all <- list()
  percentile_rlw_ge1_all <- list()
  pctl_Sample_rgw_log_all <- list()
  pctl_Sample_rlw_log_all <- list()
  pctl_Sample_ld_all <- list()
  percentile_ld_ge0_all <- list()

  medianRefPop_rgw_all <- list()
  medianRefPop_rlw_all <- list()
  medianRefPop_rgw_log_all <- list()
  medianRefPop_rlw_log_all <- list()
  medianRefPop_ld_all <- list()

  # Reference population's own per-hypothesis median value, for each
  # output_type that the "Pctl. median ref. pop." column is computed for.
  # Fixed across all iterations of the loop below (it doesn't depend on
  # 'name'), so computed once here rather than per-population. rgw_log/
  # rlw_log need their own log() here since the rgw_log_combined/
  # rlw_log_combined lists aren't built until after this loop (see below).
  ref_median_rgw <- apply(rgw_combined[[reference_pop_name]], 2, median, na.rm = TRUE)
  ref_median_rlw <- apply(rlw_combined[[reference_pop_name]], 2, median, na.rm = TRUE)
  ref_median_rgw_log <- apply(log(rgw_combined[[reference_pop_name]]), 2, median, na.rm = TRUE)
  ref_median_rlw_log <- apply(log(rlw_combined[[reference_pop_name]]), 2, median, na.rm = TRUE)
  ref_median_ld <- apply(ld_combined[[reference_pop_name]], 2, median, na.rm = TRUE)

  # Loop through each pop_es category to fill in the CI benchmark lists
  for (name in names(results)) {
    rgw_combined_values <- rgw_combined[[name]]
    rlw_combined_values  <- rlw_combined[[name]]
    ld_combined_values  <- ld_combined[[name]]
    
    # Prepare rlw_ge1 and ld_ge0 matrices
    rlw_ge1 <- rlw_combined_values
    rlw_ge1[rlw_combined_values < 1] <- 1 / rlw_combined_values[rlw_combined_values < 1]
    ld_ge0 <- abs(ld_combined_values)

    # log(rgw)/log(rlw) -- a scale-invariant, symmetric-around-0 alternative
    # to rgw/rlw themselves: log(rgw)/log(rlw) equals the log-odds (logit) of
    # the corresponding pairwise-rescaled weight (e.g. lw_pref/(lw_pref+lw_k)),
    # since the shared (lw_pref+lw_k) normalizer cancels out of the ratio --
    # so it carries no new information beyond rgw/rlw, but is far better
    # suited to things like the 'Observed'-vs-'No Effect' overlap comparison:
    # unbounded and symmetric around 0 (rather than living on the heavily
    # right-skewed (0, Inf) scale that rgw/rlw do, where '1' isn't a natural
    # visual center), and not distorted by how much weight other hypotheses
    # in the set happen to be carrying (that dilution cancels out of the
    # ratio -- and so out of its log -- exactly, per draw). Self-comparison
    # (pref vs pref) is log(1) = 0 here, rather than the 1 that rgw/rlw use,
    # so it gets cleaned via the 0-baseline (like ld) rather than the
    # 1-baseline used for rgw/rlw below.
    rgw_log_combined_values <- log(rgw_combined_values)
    rlw_log_combined_values <- log(rlw_combined_values)
    
    # Loop through the hypotheses and calculate the quantiles
    for (j in seq_len(nr.hypos)) {
      CI_benchmarks_rgw[j, 2:(1 + length(quant))] <- quantile(rgw_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_rlw[j, 2:(1 + length(quant))] <- quantile(rlw_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_rlw_ge1[j, 2:(1 + length(quant))] <- quantile(rlw_ge1[, j], quant, na.rm = TRUE)
      CI_benchmarks_rgw_log[j, 2:(1 + length(quant))] <- quantile(rgw_log_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_rlw_log[j, 2:(1 + length(quant))] <- quantile(rlw_log_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_ld[j, 2:(1 + length(quant))] <- quantile(ld_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_ld_ge0[j, 2:(1 + length(quant))] <- quantile(ld_ge0[, j], quant, na.rm = TRUE)
    }
    # Loop through the hypotheses and calculate the percentile of the sample
    # finding ('Sample' value; pctl_Sample_* below) as well as the percentile
    # of the reference population's own median value within this
    # population's own distribution ("Pctl. median ref. pop."; medianRefPop_*
    # below).
    percentile_rgw <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_rlw <- matrix(NA, nrow = nr.hypos, ncol = 1 )
    percentile_rlw_ge1 <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_rgw_log <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_rlw_log <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_ld <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_ld_ge0 <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_rgw <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_rlw <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_rgw_log <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_rlw_log <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_ld <- matrix(NA, nrow = nr.hypos, ncol = 1)
    # Hardcode this population's own row to exactly 50 when it IS the
    # reference population -- see the comment above ref_median_rgw etc.
    is_ref_pop <- identical(name, reference_pop_name)
    for (j in seq_len(nr.hypos)) {
      Fn <- ecdf(rgw_combined_values[, j])
      percentile_rgw[j, 1] <- Fn(CI_benchmarks_rgw[j, 1]) * 100
      medianRefPop_rgw[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_rgw[j]) * 100
      #
      Fn <- ecdf(rlw_combined_values[, j])
      percentile_rlw[j, 1] <- Fn(CI_benchmarks_rlw[j, 1]) * 100
      medianRefPop_rlw[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_rlw[j]) * 100
      #
      Fn <- ecdf(rlw_ge1[, j])
      percentile_rlw_ge1[j, 1] <- Fn(CI_benchmarks_rlw_ge1[j, 1]) * 100
      #
      Fn <- ecdf(rgw_log_combined_values[, j])
      percentile_rgw_log[j, 1] <- Fn(CI_benchmarks_rgw_log[j, 1]) * 100
      medianRefPop_rgw_log[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_rgw_log[j]) * 100
      #
      Fn <- ecdf(rlw_log_combined_values[, j])
      percentile_rlw_log[j, 1] <- Fn(CI_benchmarks_rlw_log[j, 1]) * 100
      medianRefPop_rlw_log[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_rlw_log[j]) * 100
      #
      Fn <- ecdf(ld_combined_values[, j])
      percentile_ld[j, 1] <- Fn(CI_benchmarks_ld[j, 1]) * 100
      medianRefPop_ld[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_ld[j]) * 100
      #
      Fn <- ecdf(ld_ge0[, j])
      percentile_ld_ge0[j, 1] <- Fn(CI_benchmarks_ld_ge0[j, 1]) * 100
    }

    # Label the percentiles (percentage of the benchmark distribution at or
    # below the sample value / reference-population median, i.e. on a 0-100
    # scale)
    percentile_names <- paste(pref_hypo_name, names(object$ratio.gw[pref_hypo, ]))
    rownames(percentile_rgw) <- rownames(percentile_rlw) <-
      rownames(percentile_rlw_ge1) <- rownames(percentile_rgw_log) <-
      rownames(percentile_rlw_log) <- rownames(percentile_ld) <-
      rownames(percentile_ld_ge0) <-
      rownames(medianRefPop_rgw) <- rownames(medianRefPop_rlw) <-
      rownames(medianRefPop_rgw_log) <- rownames(medianRefPop_rlw_log) <-
      rownames(medianRefPop_ld) <- percentile_names
    colnames(percentile_rgw) <- colnames(percentile_rlw) <-
      colnames(percentile_rgw_log) <- colnames(percentile_rlw_log) <-
      colnames(percentile_ld) <- "Pctl. Sample"
    colnames(percentile_rlw_ge1) <- colnames(percentile_ld_ge0) <- "percentile"
    colnames(medianRefPop_rgw) <- colnames(medianRefPop_rlw) <-
      colnames(medianRefPop_rgw_log) <- colnames(medianRefPop_rlw_log) <-
      colnames(medianRefPop_ld) <- "Pctl. median ref. pop."

    # Store this pop_es category's percentiles so they survive past this
    # loop iteration (mirrors the CI_benchmarks_*_all pattern below)
    pctl_Sample_rgw_all[[name]] <- percentile_rgw
    pctl_Sample_rlw_all[[name]] <- percentile_rlw
    percentile_rlw_ge1_all[[name]] <- percentile_rlw_ge1
    pctl_Sample_rgw_log_all[[name]] <- percentile_rgw_log
    pctl_Sample_rlw_log_all[[name]] <- percentile_rlw_log
    pctl_Sample_ld_all[[name]] <- percentile_ld
    percentile_ld_ge0_all[[name]] <- percentile_ld_ge0

    medianRefPop_rgw_all[[name]] <- medianRefPop_rgw
    medianRefPop_rlw_all[[name]] <- medianRefPop_rlw
    medianRefPop_rgw_log_all[[name]] <- medianRefPop_rgw_log
    medianRefPop_rlw_log_all[[name]] <- medianRefPop_rlw_log
    medianRefPop_ld_all[[name]] <- medianRefPop_ld

    # Set column names for the CI benchmarks
    colnames(CI_benchmarks_rgw) <- colnames(CI_benchmarks_rlw) <-
      colnames(CI_benchmarks_rlw_ge1) <- colnames(CI_benchmarks_rgw_log) <-
      colnames(CI_benchmarks_rlw_log) <- colnames(CI_benchmarks_ld) <-
      colnames(CI_benchmarks_ld_ge0) <- names_quant

    # Set row names for the CI benchmarks
    rownames(CI_benchmarks_rgw) <- rownames(CI_benchmarks_rlw) <-
      rownames(CI_benchmarks_rlw_ge1) <- rownames(CI_benchmarks_rgw_log) <-
      rownames(CI_benchmarks_rlw_log) <- rownames(CI_benchmarks_ld) <-
      rownames(CI_benchmarks_ld_ge0) <- paste(pref_hypo_name, names(object$ratio.gw[pref_hypo, ]))

    # Store CI benchmarks in lists
    CI_benchmarks_rgw_all[[name]] <- CI_benchmarks_rgw
    CI_benchmarks_rlw_all[[name]] <- CI_benchmarks_rlw
    CI_benchmarks_rlw_ge1_all[[name]] <- CI_benchmarks_rlw_ge1
    CI_benchmarks_rgw_log_all[[name]] <- CI_benchmarks_rgw_log
    CI_benchmarks_rlw_log_all[[name]] <- CI_benchmarks_rlw_log
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

  # rgw_log/rlw_log's self-comparison row is log(1) = 0, not 1 -- clean
  # against the 0 baseline (like ld/ld_ge0) rather than the 1 baseline used
  # for rgw/rlw/rlw_ge1 above.
  CI_benchmarks_rgw_log_all_cleaned <- lapply(CI_benchmarks_rgw_log_all, function(pop_es_list) {
    remove_single_value_rows(pop_es_list, 0)
  })

  CI_benchmarks_rlw_log_all_cleaned <- lapply(CI_benchmarks_rlw_log_all, function(pop_es_list) {
    remove_single_value_rows(pop_es_list, 0)
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
  pctl_Sample_rgw_all_cleaned <- align_rows(pctl_Sample_rgw_all, CI_benchmarks_rgw_all_cleaned)
  pctl_Sample_rlw_all_cleaned <- align_rows(pctl_Sample_rlw_all, CI_benchmarks_rlw_all_cleaned)
  percentile_rlw_ge1_all_cleaned <- align_rows(percentile_rlw_ge1_all, CI_benchmarks_rlw_ge1_all_cleaned)
  pctl_Sample_rgw_log_all_cleaned <- align_rows(pctl_Sample_rgw_log_all, CI_benchmarks_rgw_log_all_cleaned)
  pctl_Sample_rlw_log_all_cleaned <- align_rows(pctl_Sample_rlw_log_all, CI_benchmarks_rlw_log_all_cleaned)
  pctl_Sample_ld_all_cleaned <- align_rows(pctl_Sample_ld_all, CI_benchmarks_ld_all_cleaned)
  percentile_ld_ge0_all_cleaned <- align_rows(percentile_ld_ge0_all, CI_benchmarks_ld_ge0_all_cleaned)

  # Same alignment treatment for the new "Pctl. median ref. pop." matrices.
  pctl_medianRefPop_rgw_all_cleaned <- align_rows(medianRefPop_rgw_all, CI_benchmarks_rgw_all_cleaned)
  pctl_medianRefPop_rlw_all_cleaned <- align_rows(medianRefPop_rlw_all, CI_benchmarks_rlw_all_cleaned)
  pctl_medianRefPop_rgw_log_all_cleaned <- align_rows(medianRefPop_rgw_log_all, CI_benchmarks_rgw_log_all_cleaned)
  pctl_medianRefPop_rlw_log_all_cleaned <- align_rows(medianRefPop_rlw_log_all, CI_benchmarks_rlw_log_all_cleaned)
  pctl_medianRefPop_ld_all_cleaned <- align_rows(medianRefPop_ld_all, CI_benchmarks_ld_all_cleaned)

  # rgw_log/rlw_log combined draws, for combined_values/overlap below -- same
  # log() transform as CI_benchmarks_rgw_log/rlw_log above, just derived
  # directly from the (still self-column-including) rgw_combined/rlw_combined
  # matrices before those get trimmed just below.
  rgw_log_combined <- lapply(rgw_combined, function(m) log(m))
  rlw_log_combined <- lapply(rlw_combined, function(m) log(m))

  rgw_combined <- lapply(rgw_combined, function(pop_es_list) {
    remove_single_value_col(pop_es_list, 1)
  })

  rlw_combined <- lapply(rlw_combined, function(pop_es_list) {
    remove_single_value_col(pop_es_list, 1)
  })

  # rgw_log/rlw_log's self-column is log(1) = 0 -- trim against the 0
  # baseline (like ld) rather than the 1 baseline used for rgw/rlw above.
  rgw_log_combined <- lapply(rgw_log_combined, function(pop_es_list) {
    remove_single_value_col(pop_es_list, 0)
  })

  rlw_log_combined <- lapply(rlw_log_combined, function(pop_es_list) {
    remove_single_value_col(pop_es_list, 0)
  })

  ld_combined <- lapply(ld_combined, function(pop_es_list) {
    remove_single_value_col(pop_es_list, 0)
  })


  # Overlap (0-1 overlapping coefficient) between the reference population's
  # benchmark distribution and each other population's -- the numeric
  # counterpart of the overlap visible when plotting them together. The
  # reference population is "Observed" when present, else the last
  # pop_es/pop_est population supplied (see determine_overlap_reference());
  # it is the same population across gw/lw/rgw/rlw/rgw_log/rlw_log/ld since
  # they're all built from the same set of populations, so it only needs
  # computing once here and is exposed below as overlap_reference_pop for
  # print.R to build the "Overlap with <name>" column header.
  # Computed separately for gw, lw, rgw, rlw and ld (rather than assuming
  # gw and rgw -- or lw and rlw -- give the same overlap): rgw is only a
  # bijective transform of gw (and so shares its overlap) in the special
  # case of exactly 2 hypotheses; with more hypotheses, or for lw/rlw
  # (lw is not normalized to sum to 1 the way gw is), that need not hold.
  overlap_reference_pop <- reference_pop_name
  overlap_gw  <- compute_overlap_vs_observed(gw_combined)
  overlap_lw  <- compute_overlap_vs_observed(lw_combined)
  overlap_rgw <- compute_overlap_vs_observed_matrix(rgw_combined)
  overlap_rlw <- compute_overlap_vs_observed_matrix(rlw_combined)
  overlap_rgw_log <- compute_overlap_vs_observed_matrix(rgw_log_combined)
  overlap_rlw_log <- compute_overlap_vs_observed_matrix(rlw_log_combined)
  overlap_ld  <- compute_overlap_vs_observed_matrix(ld_combined)

  OUT <- list(
    benchmarks_gw = CI_benchmarks_gw,
    benchmarks_lw = CI_benchmarks_lw,
    benchmarks_rgw = CI_benchmarks_rgw_all_cleaned,
    benchmarks_rlw = CI_benchmarks_rlw_all_cleaned,
    benchmarks_rlw_ge1 = CI_benchmarks_rlw_ge1_all_cleaned,
    benchmarks_rgw_log = CI_benchmarks_rgw_log_all_cleaned,
    benchmarks_rlw_log = CI_benchmarks_rlw_log_all_cleaned,
    benchmarks_difLL = CI_benchmarks_ld_all_cleaned,
    benchmarks_absdifLL = CI_benchmarks_ld_ge0_all_cleaned,
    pctl_Sample_gw  = pctl_Sample_gw,
    pctl_Sample_lw  = pctl_Sample_lw,
    pctl_Sample_rgw = pctl_Sample_rgw_all_cleaned,
    pctl_Sample_rlw = pctl_Sample_rlw_all_cleaned,
    percentile_rlw_ge1 = percentile_rlw_ge1_all_cleaned,
    pctl_Sample_rgw_log = pctl_Sample_rgw_log_all_cleaned,
    pctl_Sample_rlw_log = pctl_Sample_rlw_log_all_cleaned,
    pctl_Sample_difLL = pctl_Sample_ld_all_cleaned,
    percentile_absdifLL = percentile_ld_ge0_all_cleaned,
    pctl_medianRefPop_gw  = pctl_medianRefPop_gw,
    pctl_medianRefPop_lw  = pctl_medianRefPop_lw,
    pctl_medianRefPop_rgw = pctl_medianRefPop_rgw_all_cleaned,
    pctl_medianRefPop_rlw = pctl_medianRefPop_rlw_all_cleaned,
    pctl_medianRefPop_rgw_log = pctl_medianRefPop_rgw_log_all_cleaned,
    pctl_medianRefPop_rlw_log = pctl_medianRefPop_rlw_log_all_cleaned,
    pctl_medianRefPop_difLL = pctl_medianRefPop_ld_all_cleaned,
    overlap_gw = overlap_gw,
    overlap_lw = overlap_lw,
    overlap_rgw = overlap_rgw,
    overlap_rlw = overlap_rlw,
    overlap_rgw_log = overlap_rgw_log,
    overlap_rlw_log = overlap_rlw_log,
    overlap_ld = overlap_ld,
    overlap_reference_pop = overlap_reference_pop,
    combined_values = list(gw_combined = gw_combined,
                           lw_combined = lw_combined,
                           rgw_combined = rgw_combined,
                           rlw_combined = rlw_combined,
                           rgw_log_combined = rgw_log_combined,
                           rlw_log_combined = rlw_log_combined,
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
# NOTE on what "stable" means here -- and what it deliberately does NOT mean.
# Under the "Observed" population (only present when pop_es/pop_est was left
# NULL, so it is auto-added), the benchmark draws are simulated centered
# exactly on the observed estimate. It is tempting to then expect the
# observed goric(a) weight to sit at the 50th percentile of its own benchmark
# distribution -- but that is NOT generally true, even with unlimited draws.
# gw/lw/rgw/rlw/ld are nonlinear functions of the (multivariate normal)
# bootstrap draws, computed through an order-restricted/inequality-constrained
# log-likelihood, which involves projecting each draw onto a constraint cone.
# That projection is nonlinear, and once more than one constraint is
# simultaneously relevant -- i.e. the observed estimates sit close to more
# than one boundary of the hypothesis at once -- the resulting distribution
# need not be symmetric around the plug-in value, so its median can genuinely
# differ from the observed value. This is the same mechanism behind
# chi-bar-square (mixture) distributions in order-restricted inference, and
# is analogous to why bootstrap bias-correction (the "BC" in BCa) exists at
# all: a nonlinear statistic's bootstrap median need not equal the original
# estimate. So a percentile far from 50 is not necessarily a sign that 'iter'
# is too low -- it can be a genuine feature of the benchmark distribution,
# especially when (some of) the hypotheses are close to being (an) equality
# constraint(s) -- and is not something to strive to eliminate by adjusting
# 'iter'. Because of this, "is the percentile near 50" is NOT what is
# checked (or acted on) here; see goric_percentile_test() for a diagnostic
# that DOES test that directly, kept purely informational for now (returned,
# not printed, and not used to decide anything).
#
# What IS checked here is whether the percentile has stabilized: whether
# growing the number of draws still meaningfully changes where the sample
# value falls in the benchmark distribution. If it does not, more draws are
# unlikely to change anything else about the benchmark either, regardless of
# whether that percentile happens to be near 50.
#
# This is checked using 'gw' rather than 'rgw'/'rlw' on purpose: gw is a
# single, bounded [0, 1] weight per draw, whereas rgw/rlw = gw_pref / gw_k is
# a ratio of two such weights (see calc_ICweights() in
# goric_calculate_IC_weights.R). When gw_k is small, that ratio amplifies
# small (essentially unavoidable, e.g. optimizer-precision-level)
# fluctuations in gw_k multiplicatively, so rgw/rlw are a noisier basis for
# this particular check than gw itself.
check_iter_adequacy <- function(benchmark_results, observed_name, iter,
                                band = c(0.495, 0.505),
                                iter_min = 500, iter_step = 100, iter_max = 2000,
                                stability_tol = 1,
                                control = list(), ...) {
  if (!observed_name %in% names(benchmark_results$pctl_Sample_gw)) {
    # No "Observed" population in this run (user supplied a custom pop_es/
    # pop_est without an "Observed" category) -- nothing to check.
    return(invisible(NULL))
  }

  sample_gw <- benchmark_results$benchmarks_gw[[observed_name]][1, 1]
  sample_lw <- benchmark_results$benchmarks_lw[[observed_name]][1, 1]
  gw_draws  <- benchmark_results$combined_values$gw_combined[[observed_name]]
  lw_draws  <- benchmark_results$combined_values$lw_combined[[observed_name]]
  gw_draws  <- gw_draws[!is.na(gw_draws)]
  lw_draws  <- lw_draws[!is.na(lw_draws)]

  # Informational only, for now: a formal GORICA-based test of whether the
  # sample value sits near the 50th percentile of its own ('Observed')
  # benchmark distribution -- see the note above this function for why that
  # need not hold even asymptotically. Returned on the object (see
  # median_bias_check_gw/median_bias_check_lw in benchmark_means()/
  # benchmark_asymp()) but not printed and not used below to decide anything.
  median_bias_check_gw <- goric_percentile_test(gw_draws, sample_gw, band = band, control = control, ...)
  median_bias_check_lw <- goric_percentile_test(lw_draws, sample_lw, band = band, control = control, ...)

  # Suggests running with 'iter' left at its default (iter = NULL) instead of
  # manually raising a too-low fixed 'iter' -- only sensible if the adaptive
  # procedure could actually try more draws than the user's fixed 'iter', i.e.
  # if that 'iter' is below iter_max; if they already used iter_max or more,
  # the automatic procedure wouldn't go any further either.
  suggest_default <- function() {
    if (iter >= iter_max) return("")
    paste0(
      " Instead of manually raising 'iter', you could also run this with 'iter' left at its ",
      "default (iter = NULL): draws are then added automatically, starting at iter_min = ",
      iter_min, " and increasing by iter_step = ", iter_step, " at a time, up to a maximum of ",
      "iter_max = ", iter_max, "."
    )
  }

  # Stability check: since 'iter' here is a single, fixed, non-growing run
  # (no rounds to compare across draw-by-draw, unlike the auto-'iter' case in
  # run_benchmark_simulation()), stability is instead assessed by comparing
  # the percentile computed from the first 80% of the draws against the
  # percentile computed from the full 'iter' draws -- i.e. did the last 20%
  # of draws still meaningfully move where the sample value falls in the
  # benchmark distribution?
  n <- length(gw_draws) # == length(lw_draws)
  n80 <- max(1, floor(0.8 * n))
  percentile_gw_80  <- 100 * mean(gw_draws[seq_len(n80)] <= sample_gw)
  percentile_lw_80  <- 100 * mean(lw_draws[seq_len(n80)] <= sample_lw)
  percentile_gw_full <- 100 * mean(gw_draws <= sample_gw)
  percentile_lw_full <- 100 * mean(lw_draws <= sample_lw)
  stable_gw <- abs(percentile_gw_full - percentile_gw_80) < stability_tol
  stable_lw <- abs(percentile_lw_full - percentile_lw_80) < stability_tol

  if (!(stable_gw && stable_lw)) {
    # The closing suggestion is tailored to which output_type(s) actually
    # have not (yet) stabilized at this 'iter': if only one of 'gw'/'lw' is
    # the problem, raising 'iter' only matters to someone who cares about
    # that one; someone only interested in the other (already-stable) one
    # doesn't need to.
    describe_type <- function(pct_80, pct_full) {
      paste0(
        "percentile went from ", sprintf("%.1f", pct_80), " (based on the first 80% of ",
        "'iter') to ", sprintf("%.1f", pct_full), " (based on the full 'iter')."
      )
    }
    closing <- if (!stable_gw && !stable_lw) {
      "Consider increasing 'iter' for a more stable benchmark."
    } else if (!stable_gw) {
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
      "\nrestriktor Message: For the user-specified 'iter' = ", iter, ", the percentile of the ",
      "value based on your data (called the 'Sample' value in the output) within the benchmark ",
      "distribution under the 'Observed' population has not (yet) stabilized: it changed by ",
      stability_tol, " percentage point(s) or more over the last 20% of the draws.\n",
      "For output_type = 'gw', ", describe_type(percentile_gw_80, percentile_gw_full), "\n",
      "For output_type = 'lw', ", describe_type(percentile_lw_80, percentile_lw_full), "\n",
      closing, suggest_default()
    )
  } else if (iter < iter_min) {
    # This stability check itself has little power at a very low 'iter'
    # (comparing an 80%/20% split of, say, 10 draws is barely more than
    # comparing single draws) -- so a "looks stable" verdict can slip
    # through undetected. Flag any 'iter' below iter_min (the adaptive
    # procedure's own starting point) outright, regardless of what the
    # stability check itself concluded, rather than silently trusting a
    # check that may not have had the power to catch real instability.
    message(
      "\nrestriktor Message: The user-specified 'iter' = ", iter, " is below the recommended ",
      "starting point of iter_min = ", iter_min, " draws. In that case, there may be too little ",
      "power to detect whether the percentile of the 'Sample value' under the 'Observed' ",
      "population has actually stabilized.\n",
      "You may want to consider increasing 'iter'; either manually or by running the code with ",
      "'iter' left at its default (iter = NULL): draws are then added automatically, starting ",
      "at iter_min = ", iter_min, " and increasing by iter_step = ", iter_step, " at a time, up ",
      "to a maximum of iter_max = ", iter_max, "."
    )
  }
  invisible(list(median_bias_check_gw = median_bias_check_gw,
                median_bias_check_lw = median_bias_check_lw))
}


# Is 'sample_value' consistent with being the median of 'draws'? Rather than
# a binomial test against p = 0.5, this expresses "close to the median" as a
# small, fixed interval around 0.5 (band -- a region of practical
# equivalence, deliberately NOT scaled to the percentile estimate's own
# standard error, since doing so would just reconstruct a Wald interval that
# goric(a) would then re-analyze with that same standard error a second
# time) and lets goric(a) itself -- type = "gorica", using the percentile
# estimate's own sampling variance as VCOV -- weigh the evidence for "the
# true percentile lies in that band" against its complement. This keeps the
# whole diagnostic inside goric(a)'s own machinery rather than a separate
# hypothesis-testing framework. Returns the percentile (0-100) 'sample_value'
# falls at within 'draws', and the resulting gorica(a) weight for the
# "in-band" hypothesis (>= 0.5 is taken to mean converged, i.e. goric(a)
# favors "close to the median" over its complement).
#
# NOTE: despite the name, this is NOT a check for whether 'iter' is large
# enough -- see the long note above check_iter_adequacy() for why the
# 'Observed' population's median need not be at the 50th percentile at all
# (a real feature of the benchmark distribution near constraint boundaries,
# not a sign of too few draws). This function is kept around as an
# informational, GORICA-based diagnostic of that median-vs-sample-value gap
# specifically -- both check_iter_adequacy() and run_benchmark_simulation()
# call it and return its result on the benchmark object (median_bias_check_gw
# / median_bias_check_lw) for possible future use, but neither currently acts
# on it or prints it; the adequacy checks themselves are based on whether the
# percentile has stabilized with more draws instead (see both functions).
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
# the "Observed" population's percentile (for output_type = 'gw' and/or
# 'lw') is still changing meaningfully round to round, add iter_step more
# draws -- WITHOUT discarding or redrawing the ones already computed --
# repeating until the percentile has stabilized or iter_max is reached.
# Growth is deliberately based on STABILITY of the percentile, not on
# whether it is close to 50 -- see the long note above check_iter_adequacy()
# for why the 'Observed' population's median need not be at the 50th
# percentile at all, so that is not something more draws can be expected to
# fix. If the user supplies a fixed numeric 'iter' instead, this runs
# exactly one round of that many draws (the original, non-adaptive
# behaviour); benchmark_means()/benchmark_asymp() then call
# check_iter_adequacy() themselves afterwards for that fixed-iter case (this
# function does not, to avoid messaging twice).
#
# Returns list(parallel_function_results = <as before, one element per
# pop_es/pop_est category>, iter = <final number of draws used>,
# median_bias_check_gw/median_bias_check_lw = <the last round's
# goric_percentile_test() result, informational only -- see the note above
# that function>).

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

  # Previous round's percentile per output_type, used to detect when it has
  # stopped moving (see 'stabilized' below) -- NA until a second round
  # exists to compare against.
  prev_percentile_gw <- NA_real_
  prev_percentile_lw <- NA_real_

  # NULL unless/until computed inside the repeat loop below -- stays NULL for
  # a fixed, user-specified 'iter' (auto_iter = FALSE), since that path
  # breaks out of the loop before reaching the percentile-stability check
  # (benchmark_means()/benchmark_asymp() do that check themselves afterwards,
  # via check_iter_adequacy(), for the fixed-iter case).
  chk_gw <- NULL
  chk_lw <- NULL

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
      # Informational only, for now -- see the note above goric_percentile_test().
      chk_gw <- goric_percentile_test(gw_draws, sample_gw, band = band, control = control)
      chk_lw <- goric_percentile_test(lw_draws, sample_lw, band = band, control = control)

      # Has the percentile stopped moving with more draws? This -- not
      # closeness to the 50th percentile -- is what growth is based on; see
      # the long note above check_iter_adequacy() for why the percentile
      # need not be near 50 at all, even once it has fully stabilized.
      # FALSE (not NA) on the very first round, since there is no earlier
      # percentile yet to compare against. This is a first-cut heuristic
      # based on a single round-to-round comparison -- if it turns out to
      # call "stabilized" too eagerly (a percentile can drift by a
      # percentage point or two by chance alone), requiring stability across
      # two consecutive rounds instead of one would be the natural
      # tightening.
      stable_gw <- isTRUE(!is.na(prev_percentile_gw) &&
                            abs(chk_gw$percentile - prev_percentile_gw) < stability_tol)
      stable_lw <- isTRUE(!is.na(prev_percentile_lw) &&
                            abs(chk_lw$percentile - prev_percentile_lw) < stability_tol)
      stabilized <- stable_gw && stable_lw
      prev_percentile_gw <- chk_gw$percentile
      prev_percentile_lw <- chk_lw$percentile
    } else {
      # No "Observed" category (custom pop_es/pop_est) -- nothing to check
      # against, so stop growing after the first (iter_min-sized) batch.
      stabilized <- TRUE
      chk_gw <- chk_lw <- NULL
    }

    if (stabilized || n_done >= iter_max) {
      if (auto_iter && !is.null(chk_gw)) {
        if (stabilized) {
          message(
            "\nrestriktor Message: 'iter' was not specified, so it was set automatically.\n",
            "Using iter = ", n_done, " draws, the percentile of the value based on your data ",
            "(called the 'Sample' value in the output) within the benchmark distribution under ",
            "the 'Observed' population has stabilized: over the last ", iter_step, " draws, it ",
            "changed by less than ", stability_tol, " percentage point(s).\n",
            "For output_type = 'gw', the percentile is ", sprintf("%.1f", chk_gw$percentile),
            "; for output_type = 'lw', the percentile is ", sprintf("%.1f", chk_lw$percentile),
            ".\n",
            "So, no further draws were added. Note that this percentile is not necessarily ",
            "expected to be near 50 -- that need not indicate a problem, particularly when ",
            "(some of) the hypotheses are close to being (an) equality constraint(s)."
          )
        } else {
          message(
            "\nrestriktor Message: 'iter' was not specified, so it was increased automatically ",
            "up to its maximum of iter = ", iter_max, " draws, since the percentile of the value ",
            "based on your data (called the 'Sample' value in the output) within the benchmark ",
            "distribution under the 'Observed' population had not (yet) stabilized: over the ",
            "last ", iter_step, " draws, it changed by ", stability_tol,
            " percentage point(s) or more.\n",
            "For output_type = 'gw', the percentile is ", sprintf("%.1f", chk_gw$percentile),
            "; for output_type = 'lw', the percentile is ", sprintf("%.1f", chk_lw$percentile),
            ".\n",
            "Consider re-running with a manually specified, larger 'iter' for a more stable ",
            "benchmark."
          )
        }
      }
      break
    }

    target <- min(n_done + iter_step, iter_max)
  }

  list(parallel_function_results = parallel_function_results, iter = n_done,
      median_bias_check_gw = chk_gw, median_bias_check_lw = chk_lw)
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

# Used only for the "Overlap with ..." column: an overlap of exactly 1 there
# is always the reference population's own (self-)entry, fixed at 1 by
# construction rather than actually computed (see compute_overlap_vs_observed()/
# compute_overlap_vs_observed_matrix() in this file) -- printing "1.000" would
# suggest a computed value carrying that much precision, which is misleading.
# Every other value in this column (and every value in every other column)
# still goes through the normal format_value().
format_overlap_value <- function(value) {
  if (!is.na(value) && value == 1) {
    return("1")
  }
  format_value(value)
}

# Used only for the "Pctl. median ref. pop." column: the reference
# population's own row there is hardcoded to exactly 50 by construction
# (see get_results_benchmark()'s pctl_medianRefPop_* computation, and its
# is_ref_pop/is_reference handling) rather than actually computed -- an
# ecdf() evaluated at a finite/even-n sample's own median need not land on
# exactly 50, so printing "50.000" there would suggest a computed value
# carrying that much precision, which is misleading. Every other value in
# this column (and every value in every other column) still goes through
# format_value().
format_median_ref_pop_value <- function(value) {
  if (!is.na(value) && value == 50) {
    return("50")
  }
  format_value(value)
}


# Builds the "Overlap with Observed" column added to each output_type's
# per-population benchmark table in print.benchmark(). 'overlap_source' is
# the relevant overlap_* field on the benchmark object: for gw/lw it is a
# plain vector named by population (a single overlap value per population);
# for rgw/rlw/ld it is a list named by population, each entry itself a
# vector named by alternative hypothesis (there can be more than one per
# population). 'n_rows' is the number of rows in the benchmark table this
# column is being attached to (nrow() of e.g. benchmarks_ratio_goric_weights
# [[pop_es]]). Deliberately aligns by POSITION rather than by matching
# names against the table's rownames: those rownames are formatted as
# "<preferred hypothesis> <alternative hypothesis>" (e.g. "H1 H2") for
# display, and parsing that back apart is fragile -- e.g. R drops the
# alternative hypothesis's name entirely when there is exactly one
# alternative (a 1x1 ratio.gw matrix, as in a single hypothesis vs. its
# complement), leaving no name to parse. Positional alignment instead
# matches the existing convention already used for the 'hypothesis_rate'
# column added the same way (see print.benchmark()): both the benchmark
# table's rows and the overlap data's entries are built from -- and
# filtered the same way as -- the same underlying combined draws, so their
# order lines up. Returns a column of NA when there's nothing to report --
# no "Observed" category in this run, the 'Observed' population's own row
# (no self-overlap to report), or (defensively) a length mismatch.
overlap_column <- function(overlap_source, pop_es, n_rows) {
  na_col <- rep(NA_real_, n_rows)
  if (is.null(overlap_source) || !(pop_es %in% names(overlap_source))) {
    # Note: for the gw/lw case, overlap_source is a plain named vector, and
    # `[[` on that with a non-matching name errors ("subscript out of
    # bounds") rather than returning NULL the way it would for a list (the
    # rgw/rlw/ld case) -- so membership is checked explicitly first.
    return(na_col)
  }
  vals <- overlap_source[[pop_es]]
  if (is.null(vals) || length(vals) == 0) {
    return(na_col)
  }
  vals <- unname(vals)
  if (length(vals) == n_rows) {
    return(vals)
  }
  # gw/lw case: a single overlap value for the whole population, repeated
  # across every row (normally just one: the preferred hypothesis).
  if (length(vals) == 1) {
    return(rep(vals, n_rows))
  }
  # Length mismatch that isn't the gw/lw broadcast case -- shouldn't
  # normally happen, but pad/truncate defensively rather than risk silently
  # mis-aligning a row with the wrong hypothesis's overlap value.
  out <- na_col
  out[seq_len(min(n_rows, length(vals)))] <- vals[seq_len(min(n_rows, length(vals)))]
  out
}


# Rebuilds the "Sample" + percentile% columns of an already-built benchmark
# table (e.g. x$benchmarks_goric_weights[[pop_es]]) at PRINT time, for a
# caller-supplied set of percentiles, instead of using the percentiles that
# happened to be requested via 'quant' back when the (often expensive,
# bootstrap-based) benchmark object was computed. Mirrors exactly what
# plot.benchmark()'s own 'percentiles' argument already does (see
# goric_benchmark_plot.R): both recompute quantile() fresh from the raw
# combined draws stored on the object (x$combined_values), so you can look at
# different percentiles without rerunning benchmark_means()/benchmark_asymp().
# 'existing_mat' is the current table (used only for its "Sample" column and
# its row names/count -- both stay unchanged); 'combined_data' is the
# matching raw-draws entry from x$combined_values for this population: a
# plain numeric vector for gw/lw, or a matrix (one column per row of
# 'existing_mat', in the same order -- gw/lw and rgw/rlw/rgw_log/rlw_log/ld
# align this way throughout the codebase, e.g. overlap_column() above relies
# on the same positional correspondence) for rgw/rlw/rgw_log/rlw_log/ld.
recompute_percentile_table <- function(existing_mat, combined_data, percentiles) {
  sample_col <- existing_mat[, 1, drop = FALSE]
  pct_names <- paste0(percentiles * 100, "%")
  if (is.null(dim(combined_data))) {
    # gw/lw case: a single row, one shared set of draws.
    q <- unname(quantile(combined_data, probs = percentiles, na.rm = TRUE))
    new_mat <- matrix(c(sample_col, q), nrow = nrow(existing_mat))
  } else {
    # rgw/rlw/rgw_log/rlw_log/ld case: one column of draws per row (hypothesis).
    q <- t(vapply(seq_len(nrow(existing_mat)), function(j) {
      quantile(combined_data[, j], probs = percentiles, na.rm = TRUE)
    }, numeric(length(percentiles))))
    new_mat <- cbind(sample_col, q)
  }
  colnames(new_mat) <- c("Sample", pct_names)
  rownames(new_mat) <- rownames(existing_mat)
  new_mat
}


# Column-header display text for print_grouped_header_table() below, keyed
# off the (internal, formatter-matching) column names already set on the
# benchmark tables by print.benchmark() -- "Sample", a percentile like "5%",
# "Pctl. Sample", "Pctl. median ref. pop.", "hypothesis_rate", or an
# "Overlap with <population>" header. Returns, per column: 'line1' (blank
# for a column that instead shares a merged group label -- see 'group'
# below), 'line2', and 'group' (NA for a column with its own two-line
# header; otherwise the shared label text for the columns sharing a merged
# top cell, e.g. "Percentiles distributions" for the quantile columns).
# Layout per Rebecca's spec:
#   Sample                  -> "Sample" / "value" (own two-line header)
#   5%, 35%, ... (any N%)   -> grouped "Percentiles distributions", showing
#                              the percentile itself on line 2
#   Pctl. Sample            -> grouped "Percentiles (pctl.) comparisons",
#                              "Sample value" on line 2
#   Pctl. median ref. pop.  -> same group, "Median ref. pop." on line 2
#   hypothesis_rate         -> "Hypothesis" / "rate" (own two-line header)
#   Overlap with <pop>      -> "Overlap" / "with <pop>" (own two-line header,
#                              second line already dynamic via overlap_header)
# Anything unrecognized (shouldn't normally occur) falls back to a blank
# line 1 and the column's own name on line 2, so nothing is ever dropped.
header_lines_for_columns <- function(colnames_vec) {
  n <- length(colnames_vec)
  line1 <- character(n)
  line2 <- character(n)
  group <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    cn <- colnames_vec[i]
    if (identical(cn, "Sample")) {
      line1[i] <- "Sample"
      line2[i] <- "value"
    } else if (grepl("^[0-9.]+%$", cn)) {
      group[i] <- "Percentiles distributions"
      line2[i] <- cn
    } else if (identical(cn, "Pctl. Sample")) {
      group[i] <- "Percentiles (pctl.) comparisons"
      line2[i] <- "Sample value"
    } else if (identical(cn, "Pctl. median ref. pop.")) {
      group[i] <- "Percentiles (pctl.) comparisons"
      line2[i] <- "Median ref. pop."
    } else if (identical(cn, "hypothesis_rate")) {
      line1[i] <- "Hypothesis"
      line2[i] <- "rate"
    } else if (grepl("^Overlap with ", cn)) {
      line1[i] <- "Overlap"
      line2[i] <- sub("^Overlap with ", "with ", cn)
    } else {
      line2[i] <- cn
    }
  }
  list(line1 = line1, line2 = line2, group = group)
}

# Prints 'formatted_df' (a character matrix, already formatted -- see
# print_rounded_es_value() below) with a 2-line column header, where
# columns sharing the same non-NA 'group' (from header_lines_for_columns())
# get a single label centered across their combined width on header line 1
# -- a text "merged cell" -- instead of each getting a separate line-1 cell.
# This deliberately does NOT use R's own print(): base R print() has no
# notion of a header spanning two lines or multiple columns, and computes
# its column widths/line-wrapping independently of any second header line,
# so a naive two-line header would drift out of alignment with the data
# the moment the console width changes (print() may then wrap the table
# into column blocks, and nothing keeps a second header line in sync with
# that). Here, every column width -- and the width each merged group label
# needs -- is computed and padded by hand up front, so the two header lines
# and the data stay aligned regardless of getOption("width").
print_grouped_header_table <- function(formatted_df, rn) {
  ncol_ <- ncol(formatted_df)
  colnames_vec <- colnames(formatted_df)
  hdr <- header_lines_for_columns(colnames_vec)

  col_w <- vapply(seq_len(ncol_), function(j) {
    max(nchar(hdr$line1[j]), nchar(hdr$line2[j]), nchar(formatted_df[, j]))
  }, integer(1))

  # If a group's shared label is wider than the columns it spans (plus the
  # single-space separator between them), widen that group's last column to
  # make room, rather than truncating or overflowing the label.
  group_ids <- unique(stats::na.omit(hdr$group))
  for (g in group_ids) {
    idx <- which(hdr$group == g)
    span_w <- sum(col_w[idx]) + (length(idx) - 1)
    if (nchar(g) > span_w) {
      col_w[idx[length(idx)]] <- col_w[idx[length(idx)]] + (nchar(g) - span_w)
    }
  }

  rn_w <- max(nchar(rn), 0)
  pad_left <- function(s, w) formatC(s, width = -w)
  pad_center <- function(s, w) {
    total <- w - nchar(s)
    left <- total %/% 2
    right <- total - left
    paste0(strrep(" ", left), s, strrep(" ", right))
  }

  # Header line 1, built in column order: a merged group's columns collapse
  # into one centered label spanning their combined width; an ungrouped
  # column shows its own (possibly blank) line1 text.
  line1_parts <- character(0)
  j <- 1
  while (j <= ncol_) {
    if (!is.na(hdr$group[j])) {
      idx <- which(hdr$group == hdr$group[j])
      span_w <- sum(col_w[idx]) + (length(idx) - 1)
      line1_parts <- c(line1_parts, pad_center(hdr$group[j], span_w))
      j <- max(idx) + 1
    } else {
      line1_parts <- c(line1_parts, pad_left(hdr$line1[j], col_w[j]))
      j <- j + 1
    }
  }
  line2_parts <- mapply(pad_left, hdr$line2, col_w)

  cat(pad_left("", rn_w), " ", paste(line1_parts, collapse = " "), "\n", sep = "")
  cat(pad_left("", rn_w), " ", paste(line2_parts, collapse = " "), "\n", sep = "")
  for (i in seq_len(nrow(formatted_df))) {
    cat(pad_left(rn[i], rn_w), " ",
        paste(mapply(pad_left, formatted_df[i, ], col_w), collapse = " "), "\n", sep = "")
  }
}


# called by the benchmark.print() function. is_reference: TRUE when 'pop_es'
# is the reference population that overlap/pctl_medianRefPop are computed
# against (see determine_overlap_reference()/x$overlap_reference) -- appends
# "(Reference population)" to the printed population label, e.g. "Population
# effect-size = Observed (Reference population)".
print_rounded_es_value <- function(df, pop_es, model_type, text_color, reset,
                                   is_reference = FALSE) {
  if (model_type == "benchmark_asymp") {
    pop_es_value <- gsub("pop_est = ", "", pop_es)
    label <- "Population estimates"
  } else {
    pop_es_value <- gsub("pop_es = ", "", pop_es)
    label <- "Population effect-size"
  }
  if (is_reference) {
    pop_es_value <- paste0(pop_es_value, " (Reference population)")
  }
  cat(sprintf("%s = %s%s%s\n", label, text_color, pop_es_value, reset))

  #formatted_column <- sprintf("%.3f", df)
  # The "Overlap with ..." column (if present -- see overlap_column() /
  # print.benchmark()) is formatted with format_overlap_value() instead of
  # format_value(), so its self-overlap entries print as "1" rather than
  # "1.000" (that value isn't actually computed, it's fixed by construction --
  # see format_overlap_value()'s own comment). Likewise the "Pctl. median
  # ref. pop." column (if present) uses format_median_ref_pop_value() so the
  # reference population's own hardcoded-50 entry prints as "50" rather than
  # "50.000". Every other column keeps using format_value() as before.
  is_overlap_col <- grepl("^Overlap with ", colnames(df))
  is_median_ref_col <- colnames(df) == "Pctl. median ref. pop."
  formatted_cols <- lapply(seq_len(ncol(df)), function(j) {
    col_formatter <- if (is_overlap_col[j]) {
      format_overlap_value
    } else if (is_median_ref_col[j]) {
      format_median_ref_pop_value
    } else {
      format_value
    }
    vapply(df[, j], col_formatter, character(1))
  })
  formatted_df <- do.call(cbind, formatted_cols)
  rownames(formatted_df) <- rownames(df)
  colnames(formatted_df) <- colnames(df)
  print_grouped_header_table(formatted_df, rownames(df))
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
