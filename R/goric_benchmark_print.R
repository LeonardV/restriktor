print.benchmark <- function(x, output_type = c("rgw", "gw", "rlw", "ld", "all"), 
                            hypo_rate_threshold = 1, color = TRUE, ...) {

  # Ensure the object is of class 'benchmark_means'
  if (!inherits(x, "benchmark")) {
    stop("Invalid object. The object should be of class 'benchmark'.", call. = FALSE)
  }
  
  output_type <- match.arg(output_type, c("rgw", "gw", "rlw", "ld", "all"))
  
  ldots <- list(...)
  
  # compute hypothesis rate
  x$hypothesis_rate <- lapply(x$combined_values$rgw_combined, 
                              calculate_hypothesis_rate, q = hypo_rate_threshold)
  
  # no effect names
  NE_names <- names(x$hypothesis_rate) %in% c("pop_est = No-effect", "pop_es = 0")
  x$hypothesis_rate[NE_names] <- as.numeric(NA)
  
  # number of failed bootstrap runs
  empty_lists_count <- sapply(x$combined_values$gw_combined, function(x) { 
    attr(x, "empty_lists_count") } )
  #max_empty_lists_count <- max(empty_lists_count)
  
  model_type <- class(x)[1]
  goric_type <- toupper(x$type)
  pref_hypo <- x$pref_hypo_name
  error_prob_pref_hypo <- x$error_prob_pref_hypo
  error_prob_pref_hypo <- 
    if (error_prob_pref_hypo < 0.001) {
      "<.001"
    } else {
      sprintf("%.3f", error_prob_pref_hypo)
    }
  
  # means specific
  if (inherits(x, "benchmark_means")) {
    formatted_values <- sprintf("%.3f", x$pop_es)
    ratio_pop_means <- x$ratio_pop_means
    group_size <- x$group_size
    ngroups <- x$ngroups
    cohens_f_observed <- x$cohens_f_observed
  } else if (inherits(x, "benchmark_asymp")) {
    pop_est <- x$pop_est
    formatted_values <- sprintf("%.3f", pop_est)
    dim(formatted_values) <- dim(pop_est)
    dimnames(formatted_values) <- dimnames(pop_est)
    group_size <- x$sample_size
    ngroups <- x$n_coef
  } 
  
  output_format <- ldots$output_format
  if (is.null(output_format)) {
    output_format <- "console"
  }
  
  # ANSI escape codes for console colors
  if (color) {
    # ANSI escape codes for console colors
    if (output_format == "console") {
      blue <- "\033[34m"
      green <- "\033[32m"
      orange <- "\033[38;5;214m"
      reset <- "\033[0m"
    } else if (output_format == "html") {
      blue <- "<span style='color:blue'>"
      green <- "<span style='color:green'>"
      orange <- "<span style='color:orange'>"
      reset <- "</span>"
    } else if (output_format == "latex") {
      blue <- "\\textcolor{blue}{"
      green <- "\\textcolor{green}{"
      orange <- "\\textcolor{orange}{"
      reset <- "}"
    }
  } else {
    # R Markdown cannot handle color codes (console, html and latex)
    blue <- green <- orange <- reset <- ""
  }
  
  text_gw  <- paste0("Benchmark: Percentiles of ", orange, goric_type, " Weights", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_rgw <- paste0("Benchmark: Percentiles of ", orange, "Ratio-of-", goric_type, "-weights", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_rlw <- paste0("Benchmark: Percentiles of ", orange, "Ratio-of-log-likelihood-weights", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_ld  <- paste0("Benchmark: Percentiles of ", orange, "Differences in Log-likelihood Values", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  
  cat("\n")
  #cat(strrep("=", 70), "\n")
  cat(paste0(blue, "Benchmark Results", reset), "\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("Preferred Hypothesis: %s%s%s\n", green, pref_hypo, reset))
  cat(sprintf("Error probability Preferred Hypothesis vs. Complement: %s%s%s\n", green, error_prob_pref_hypo, reset))
  if (inherits(x, "benchmark_means")) {
    cat(sprintf("Number of Groups: %s%s%s\n", green, ngroups, reset))
  } else {
    #if (group_size != "") {
    if (!is.null(group_size)) {
      cat(sprintf("Sample Size: %s%s%s\n", green, paste(group_size, collapse = ", "), reset))
    } 
    cat(sprintf("Number of Parameters: %s%s%s\n", green, ngroups, reset))
  }
  if (inherits(x, "benchmark_means")) {
    cat(sprintf("Group Sizes: %s%s%s\n", green, paste(group_size, collapse = ", "), reset))
    cat(sprintf("Ratio of Population Means: %s%s%s\n", green, paste(sprintf("%.3f", ratio_pop_means), collapse = ", "), reset))
    cat(sprintf("Population Effect-Sizes (Cohens f): %s%s%s\n", green,  paste(formatted_values, collapse = ", "), reset))
    cat(sprintf("Observed Effect-Size (Cohens f): %s%s%s\n", green, sprintf("%.3f", cohens_f_observed), reset))
  } else {
    print_formatted_matrix(formatted_values, green, reset)
  }
  #cat(strrep("-", 70), "\n")
  cat("\n")
  
# -------------------------------------------------------------------------
  R <- x$iter
  succesful_draws <- R - empty_lists_count
  if (any(succesful_draws < R)) { 
    cat("Number of requested bootstrap draws:", R, "\n")
    for (i in seq_along(succesful_draws)) {
      cat(sprintf("  Number of successful bootstrap draws for %*s: %s\n", 
                  10, names(succesful_draws)[i], succesful_draws[i]))
    }
    text_msg <- paste("Advise: If a substantial number of bootstrap draws fail to converge,", 
                      "it is advisable to increase the number of bootstrap iterations.")
    message("---\n", text_msg)
  }

  # Normalize output_type to lowercase
  output_type <- tolower(output_type)
  
  # Print benchmarks for Goric-weights percentiles
  if ("all" %in% output_type || "gw" %in% output_type) { 
    print_section(
      text_gw,
      function() {
        for (pop_es in names(x$benchmarks_goric_weights)) {
          print_rounded_es_value(x$benchmarks_goric_weights[[pop_es]], pop_es, 
                                 model_type, green, reset)
        }
      }, nchar(text_gw), text_color = blue, reset = reset
    )
  }
  
  # Print benchmarks for ratio Goric-weights percentiles
  if ("all" %in% output_type || "rgw" %in% output_type) {
    
    # Loop over alle sets in both benchmarks_ratio_goric_weights and hypothesis_rate
    for (pop_es_name in names(x$benchmarks_ratio_goric_weights)) {
      x$benchmarks_ratio_goric_weights[[pop_es_name]] <- cbind(
        x$benchmarks_ratio_goric_weights[[pop_es_name]],
        hypothesis_rate = x$hypothesis_rate[[pop_es_name]] # x$hypothesis_rate
      )
    }
    
    if (any(NE_names)) {
      x$benchmarks_ratio_goric_weights[NE_names][[1]] <- 
        x$benchmarks_ratio_goric_weights[NE_names][[1]][, -7, drop = FALSE]
    }
    
    print_section(
      text_rgw,
      function() {
        for (pop_es in names(x$benchmarks_ratio_goric_weights)) {
          print_rounded_es_value(x$benchmarks_ratio_goric_weights[[pop_es]], pop_es, 
                                 model_type, green, reset)
        }
      }, nchar(text_rgw), text_color = blue, reset = reset
    )
  }
  
  # Print benchmarks for ratio log-likelihood-weights percentiles
  if ("all" %in% output_type || "rlw" %in% output_type) {
    print_section(
      text_rlw,
      function() {
        for (pop_es in names(x$benchmarks_ratio_ll_weights)) {
          print_rounded_es_value(x$benchmarks_ratio_ll_weights[[pop_es]], pop_es, 
                                 model_type, green, reset)
        }
      }, nchar(text_rlw), text_color = blue, reset = reset
    )
  }
  
  # Print benchmarks for difference log-likelihood-values percentiles
  if ("all" %in% output_type || "ld" %in% output_type) {
    print_section(
      text_ld,
      function() {
        for (pop_es in names(x$benchmarks_difLL)) {
          print_rounded_es_value(x$benchmarks_difLL[[pop_es]], pop_es, model_type, 
                                 green, reset)
        }
      }, nchar(text_ld), text_color = blue, reset = reset
    )
  }
  return(invisible(x))
}
