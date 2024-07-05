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

plot_all_groups <- function(plot_df, groups, title, xlabel, x_lim = NULL, 
                            alpha = 0.5, distr_grid = FALSE, 
                            percentiles = c(0.05, 0.95)) {
  plot_list <- list()
  for (group in groups) {
    plot <- create_density_plot(plot_df, group, title, xlabel, x_lim, alpha, 
                                distr_grid, percentiles)
    plot_list[[group]] <- plot
  }
  return(plot_list)
}

# benchmark plots
create_density_plot <- function(plot_df, group_comparison, title, xlabel,
                                x_lim = NULL, alpha = 0.5, distr_grid = FALSE,
                                percentiles = c(0.05, 0.95)) {
  
  df_subset <- subset(plot_df, Group_hypo_comparison == group_comparison)
  
  if (!is.null(df_subset$Group_hypo_comparison)) {
    title <- paste(title, "vs.", sub(".*vs\\. ", "", unique(df_subset$Group_hypo_comparison)))
  }
  
  # Function to calculate density
  calculate_density <- function(data, var, sample_value) {
    dens <- density(data[[var]], kernel = "gaussian", from = 0)
    data.frame(x = dens$x, y = dens$y, sample_value = sample_value)
  }
  
  # Calculate densities for each group
  # density_data <- df_subset %>%
  #   group_by(Group_pop_values) %>%
  #   do(calculate_density(., "Value", unique(.$sample_value))) %>%
  #   ungroup()
  
  group_levels <- unique(df_subset$Group_pop_values)
  density_data <- do.call(rbind, lapply(group_levels, function(group) {
    data <- subset(df_subset, Group_pop_values == group)
    dens <- calculate_density(data, var = "Value", sample_value = data$sample_value[1])
    dens$Group_pop_values <- group
    dens$sample_value <- unique(data$sample_value)
    dens
  }))

  percentile_label <- paste0("[", percentiles[1]*100, ", ", percentiles[2]*100, "]th Percentiles")
  percentile_first_group_lower <- df_subset$lb_first_group[1]
  percentile_first_group_upper <- df_subset$ub_first_group[1]
  formatted_sample_value <- sprintf("Sample Value = %.3f", unique(df_subset$sample_value)[1])
  linetype_labels <- unique(c(formatted_sample_value, percentile_label))
  
  p <- ggplot(density_data, aes(x = x, y = y, fill = Group_pop_values)) +
    geom_ribbon(aes(ymin = 0, ymax = y), alpha = alpha) +
    geom_segment(data = df_subset, aes(x = sample_value, xend = sample_value, 
                                       y = 0, yend = Inf, linetype = formatted_sample_value), 
                 color = "red", linewidth = 1) +
    geom_segment(data = df_subset, aes(x = percentile_first_group_lower, 
                                       xend = percentile_first_group_lower, y = 0, 
                                       yend = Inf, linetype = percentile_label), 
                 color = df_subset$first_group_color[1], linewidth = 1) + 
    geom_segment(data = df_subset, aes(x = percentile_first_group_upper, 
                                       xend = percentile_first_group_upper, y = 0, 
                                       yend = Inf, linetype = percentile_label), 
                 color = df_subset$first_group_color[1], linewidth = 1) + 
    #scale_x_continuous(expand = c(0, 0)) +  
    #scale_y_continuous(expand = c(0, 0)) +  
    ggtitle(title) +
    xlab(xlabel) + 
    ylab("Density") + 
    theme(axis.text = element_text(size = 11),
          axis.title.x = element_text(size = 12, margin = margin(t = 10)),
          axis.title.y = element_text(size = 12, margin = margin(r = 10)),
          plot.title = element_text(size = 12)) +
    scale_fill_brewer(palette = "Set2", name = "Distribution under:") +
    scale_linetype_manual(values = setNames(rep("solid", length(linetype_labels)), linetype_labels), name = "") +
    theme(legend.key = element_rect(fill = "white")) +
    labs(fill = "Distribution under:", linetype = "Legend") +
    guides(fill = guide_legend(order = 1),
           linetype = guide_legend(order = 2)) 
 
  if (distr_grid) {
    p <- p + facet_grid(. ~ Group_pop_values, scales = "free_x")
  }

  if (!is.null(x_lim) && length(x_lim) == 2) {
    p <- p + coord_cartesian(xlim = x_lim)
  } else {
    iqr <- IQR(df_subset$Value)
    #q1 <- quantile(df_subset$Value, 0.05)
    q3 <- quantile(df_subset$Value, 0.95)
    #lower_limit <- q1 - 1.5 * iqr
    upper_limit <- q3 + 1.5 * iqr
    p <- p + coord_cartesian(xlim = c(0, upper_limit))
  }
  
  return(p)
}

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
  cov_matrix <- VCOV * N
  ss_within <- sum((N - 1) * diag(cov_matrix)) 
  cohens_f <- sqrt(ss_between/ss_within)
  
  return(cohens_f)
}


# Compute ratio data based on group_means
compute_ratio_data <- function(group_means) {
  ngroups <- length(group_means)
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
  return(ratio_data)
}


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


# Compute the population means based on the input parameters
compute_population_means <- function(pop_es, ratio_pop_means, var_e, ngroups) {
  means_pop_all <- matrix(NA, ncol = ngroups, nrow = length(pop_es))
  nr_es <- length(pop_es)
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
        (1/sqrt(var_e)) * sqrt((1/ngroups) * sum((means_pop - mean(means_pop))^2)) - pop_es[teller_es] #  AANPASSSEN NAAR NIEUWE FORMULE
      }
      d <- uniroot(fun, lower = 0, upper = 100)$root
      # Construct means_pop
      means_pop <- ratio_pop_means*d
    }
    means_pop_all[teller_es, ] <- means_pop
  }
  return(means_pop_all)  
}


# Calculate the error probability
calculate_error_probability <- function(object, hypos, pref_hypo, est, 
                                        VCOV, control, ...) {
  # Error probability based on complement of preferred hypothesis in data
  nr_hypos <- dim(object$result)[1]
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



# this function is called from the goric_benchmark_anova() function
parallel_function_means <- function(i, N, var_e, means_pop, 
                                    hypos, pref_hypo, object, ngroups, sample, 
                                    control, form_model_org, ...) {  
  # Sample residuals
  epsilon <- rnorm(sum(N), sd = sqrt(var_e))
  
  # original model formula 
  if (length(form_model_org) > 0) {
    model <- form_model_org
    lhs <- all.vars(model)[1]
    sample[[lhs]] <- as.matrix(sample[, 2:(1 + ngroups)]) %*% matrix(means_pop, 
                                                                     nrow = ngroups) + epsilon
    df_boot <- data.frame(lhs = sample[[lhs]], sample[, 2:(1 + ngroups)])
    colnames(df_boot)[1] <- lhs
    
    has_intercept <- attr(terms(model), "intercept") == 1
    rhs <- as.character(attr(terms(model), "term.labels"))
    
    # Create the RHS with all other variables and optionally the intercept
    if (has_intercept) {
      new_rhs <- "."
    } else {
      new_rhs <- "-1 + ."
    }
    
    # Create the new formula
    new_model <- as.formula(paste(lhs, "~", new_rhs))
  } else {
    new_model <- y ~ 0 + .
    # Generate data
    sample$y <- as.matrix(sample[, 2:(1 + ngroups)]) %*% matrix(means_pop, 
                                                                nrow = ngroups) + epsilon
    df_boot <- data.frame(y = sample$y, sample[, 2:(1 + ngroups)])
  }
  
  
  # Obtain fit
  fit_boot <- lm(new_model, data = df_boot)  
  # GORICA or GORICA depending on what is done in data
  results_goric <- goric(fit_boot,
                         hypotheses = hypos,
                         comparison = object$comparison,
                         type = object$type,
                         control = control, 
                         ...)
  
  # Return the relevant results
  ld_names <- names(results_goric$ratio.gw[pref_hypo, ])
  ld <- results_goric$result$loglik[pref_hypo] - results_goric$result$loglik
  names(ld) <- ld_names
  
  list(
    #test  = attr(results.goric$objectList[[results.goric$objectNames]]$wt.bar, "mvtnorm"),
    gw  = results_goric$result[pref_hypo, 7], # goric(a) weight
    rgw = results_goric$ratio.gw[pref_hypo, ], # ratio goric(a) weights
    rlw = results_goric$ratio.lw[pref_hypo, ], # ratio likelihood weights
    ld  = ld # loglik difference
  )
}


# this function is called from the benchmark_asymp() function
parallel_function_asymp <- function(i, est, VCOV, hypos, pref_hypo, comparison,
                                    type, control, ...) {  
  
  results_goric <- goric(est[i, ], VCOV = VCOV,
                         hypotheses = hypos,
                         comparison = comparison,
                         type = type,
                         control = control, 
                         ...)
  
  # Return the relevant results
  ld_names <- names(results_goric$ratio.gw[pref_hypo, ])
  ld <- results_goric$result$loglik[pref_hypo] - results_goric$result$loglik
  names(ld) <- ld_names
  
  list(
    #test  = attr(results.goric$objectList[[results.goric$objectNames]]$wt.bar, "mvtnorm"),
    gw  = results_goric$result[pref_hypo, 7], # goric(a) weight
    rgw = results_goric$ratio.gw[pref_hypo, ], # ratio goric(a) weights
    rlw = results_goric$ratio.lw[pref_hypo, ], # ratio likelihood weights
    ld  = ld
  )
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
  if (abs(value) >= 1000 || abs(value) <= 0.001 && value != 0) {
    return(sprintf("%.3e", value))  
  } else {
    return(sprintf("%.3f", value))  
  }
}

# called by the benchmark.print() function
print_rounded_es_value <- function(df, pop_es, model_type, text_color, reset) {
  if (model_type == "benchmark_asymp") {
    pop_es_value <- gsub("pop_est = ", "", pop_es)
  } else {
    pop_es_value <- gsub("pop_es = ", "", pop_es)
  }
  cat(sprintf("Population estimates = %s%s%s\n", text_color, pop_es_value, reset))

  #formatted_column <- sprintf("%.3f", df)
  formatted_values <- sapply(as.numeric(df), format_value)
  formatted_df <- `dim<-`(formatted_values, dim(df))
  rownames(formatted_df) <- rownames(df)
  colnames(formatted_df) <- colnames(df)
  print(formatted_df, row.names = TRUE, quote = FALSE)
  cat("\n")
}


# print_rounded_es_value <- function(df, pop_es, model_type, text_color, reset) {
#   # Extract population effect estimate
#   if (model_type == "benchmark_asymp") {
#     pop_es_value <- gsub("pop_est = ", "", pop_es)
#   } else {
#     pop_es_value <- gsub("pop_es = ", "", pop_es)  
#   }
#   cat(sprintf("Population effect estimates = %s%s%s\n", text_color, pop_es_value, reset))
#   
#   # Function to format and right-align values
#   format_value <- function(value) {
#     formatted <- sprintf("%.3f", as.numeric(value))
#     return(formatted)
#   }
#   
#   # Apply formatting to the dataframe
#   formatted_values <- apply(df, c(1, 2), format_value)
#   formatted_df <- as.data.frame(formatted_values)
#   rownames(formatted_df) <- rownames(df)
#   colnames(formatted_df) <- colnames(df)
#   
#   # Function to print the dataframe with aligned columns and rows
#   print_aligned_df <- function(df) {
#     # Voeg de rijnamen als een aparte kolom toe
#     #df_with_row_names <- cbind(Row = rownames(df), df)
#     
#     # Vind de maximale breedte van elke kolom inclusief kolomnamen
#     col_widths <- sapply(df, function(col) max(nchar(as.character(col))))
#     col_widths <- pmax(col_widths, nchar(names(df)))
#     
#     # Formatteer de kolomnamen
#     colnames_formatted <- mapply(format, names(df), width = col_widths, SIMPLIFY = FALSE)
#     
#     # Print de kolomnamen
#     cat(paste(unlist(colnames_formatted), collapse = "  "), "\n")
#     
#     # Print de rijen met de rijnamen
#     for (i in 1:nrow(df)) {
#       row_values <- mapply(format, as.character(df[i, ]), width = col_widths, SIMPLIFY = FALSE)
#       cat(paste(unlist(row_values), collapse = "  "), "\n")
#     }
#   }
#   
#   # Gebruik de functie om de geformatteerde data frame af te drukken
#   print_aligned_df(formatted_df)
#   cat("\n")
# }



print_formatted_matrix <- function(mat, text_color, reset) {
  col_widths <- apply(mat, 2, function(col) max(nchar(as.character(col))))
  col_widths <- pmax(col_widths, nchar(colnames(mat))) + 2  
  
  cat("Population Estimates (PE):\n")
  cat(sprintf("%-10s", ""))
  for (k in 1:ncol(mat)) {
    cat(sprintf(paste0("%", col_widths[k], "s"), colnames(mat)[k]))
  }
  cat("\n")
  
  # Print de rijen
  for (i in 1:nrow(mat)) {
    cat(sprintf("%-10s", rownames(mat)[i]))
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
                              "greater or less than 0: ", 
                              paste(hypotheses_with_constants, collapse = ", "),
                              ". The default population estimates are likely incorrect.",
                              "Consider providing custom population estimates via the",
                              "pop_est argument.")
    warning(warning_message, call. = FALSE)
  }
}
