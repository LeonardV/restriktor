print.benchmark <- function(x, output_type = c("rgw", "gw", "rlw", "ld", "overlap", "all"),
                            hypo_rate_threshold = 1, color = TRUE, ...) {

  # Ensure the object is of class 'benchmark_means'
  if (!inherits(x, "benchmark")) {
    stop("Invalid object. The object should be of class 'benchmark'.", call. = FALSE)
  }

  output_type <- match.arg(output_type, c("rgw", "gw", "rlw", "ld", "overlap", "all"))
  
  ldots <- list(...)
  
  # compute hypothesis rate
  x$hypothesis_rate <- lapply(x$combined_values$rgw_combined, 
                              calculate_hypothesis_rate, q = hypo_rate_threshold)
  
  # no effect names
  # Was c("pop_est = No-effect", "pop_es = 0") -- but benchmark_means() names
  # this category "pop_es = No-effect" (matching benchmark_asymp's
  # "pop_est = No-effect"), never "pop_es = 0", so this never matched for
  # benchmark_means objects and the hypothesis_rate column was never hidden
  # for their No-effect category.
  NE_names <- names(x$hypothesis_rate) %in% c("pop_est = No-effect", "pop_es = No-effect")
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
    if (!is.numeric(error_prob_pref_hypo)) {
      "The unconstrained is preferred and has no complement, it contains all possible orderings."
    } else if (error_prob_pref_hypo < 0.001) {
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
  text_overlap <- paste0("Benchmark: ", orange, "Overlap", blue, " of the 'Observed' Population's Benchmark Distribution with each 'Null' Population")

  cat("\n")
  #cat(strrep("=", 70), "\n")
  cat(paste0(blue, "Benchmark Results", reset), "\n")
  cat(strrep("-", 70), "\n")
  # TO DO Geef ook het aantal iteraties (iter) wat gebruikt is.
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
    for (pop_es in names(x$benchmarks_goric_weights)) {
      x$benchmarks_goric_weights[[pop_es]] <- cbind(
        x$benchmarks_goric_weights[[pop_es]],
        x$percentile_goric_weights[[pop_es]]
      )
    }
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
        x$percentile_ratio_goric_weights[[pop_es_name]],
        hypothesis_rate = x$hypothesis_rate[[pop_es_name]] # x$hypothesis_rate
      )
    }
    # TO DO nu bij No-effect ook hypothesis_rate maar die zouden we dacht ik niet meer laten zien omdat het verwarrend is wat het betekent
    # Is nu juist weer weg, als het goed is.
    #       daarnaast heeft anders echt beschrijving nodig, want het is steun info hypo onder NE en dan in gehele set?
    #       Ws overleggen of dit handig is - nu ineens denk ik dat het zo gek nog niet is :-).
    
    if (any(NE_names)) {
      # Was hardcoded as column -7 (assuming Sample + 5 default quantiles as
      # columns 1-6, hypothesis_rate as 7); now selected by name instead,
      # since adding the percentile column shifted hypothesis_rate to column 8
      # and a positional index would silently drop the wrong column.
      x$benchmarks_ratio_goric_weights[NE_names][[1]] <-
        x$benchmarks_ratio_goric_weights[NE_names][[1]][
          , colnames(x$benchmarks_ratio_goric_weights[NE_names][[1]]) != "hypothesis_rate",
          drop = FALSE]
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
    for (pop_es in names(x$benchmarks_ratio_ll_weights)) {
      x$benchmarks_ratio_ll_weights[[pop_es]] <- cbind(
        x$benchmarks_ratio_ll_weights[[pop_es]],
        x$percentile_ratio_ll_weights[[pop_es]]
      )
    }
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
    for (pop_es in names(x$benchmarks_difLL)) {
      x$benchmarks_difLL[[pop_es]] <- cbind(
        x$benchmarks_difLL[[pop_es]],
        x$percentile_difLL[[pop_es]]
      )
    }
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

  # Print overlap of the 'Observed' population's benchmark distribution
  # against each 'null' population (by default just 'No-effect'). One table
  # for gw/lw (a single overlap value per population), plus -- whenever
  # present -- a block per matrix-valued statistic (rgw/rlw/ld, one value
  # per alternative hypothesis), since gw and rgw (and lw and rlw) need not
  # have the same overlap: rgw is only a bijective transform of gw (and so
  # shares its overlap) in the special case of exactly 2 hypotheses; with
  # more hypotheses, or for lw/rlw (lw is not normalized to sum to 1 the
  # way gw is), that need not hold. Skipped (no error) when there's no
  # "Observed" category to compare against (a custom pop_es/pop_est without
  # one), same as the other overlap-dependent diagnostics.
  if ("all" %in% output_type || "overlap" %in% output_type) {
    if (!is.null(x$overlap_goric_weights) && length(x$overlap_goric_weights) > 0) {
      print_section(
        text_overlap,
        function() {
          other_names <- names(x$overlap_goric_weights)
          overlap_tab <- cbind(
            gw = x$overlap_goric_weights[other_names],
            lw = x$overlap_ll_weights[other_names]
          )
          formatted_values <- sapply(as.numeric(overlap_tab), format_value)
          formatted_tab <- `dim<-`(formatted_values, dim(overlap_tab))
          rownames(formatted_tab) <- gsub("^pop_es(t)? = ", "", other_names)
          colnames(formatted_tab) <- colnames(overlap_tab)
          print(formatted_tab, quote = FALSE)
          cat("(0 = no overlap, 1 = fully overlapping distributions, versus the 'Observed' population)\n")

          print_matrix_overlap <- function(overlap_list, label) {
            if (is.null(overlap_list) || length(overlap_list) == 0) return(invisible(NULL))
            cat("\n", label, ":\n", sep = "")
            for (pop_es in names(overlap_list)) {
              vals <- overlap_list[[pop_es]]
              if (length(vals) == 0) next
              formatted <- sapply(vals, format_value)
              names(formatted) <- names(vals)
              cat(sprintf("  %s: ", gsub("^pop_es(t)? = ", "", pop_es)))
              cat(paste(names(formatted), formatted, sep = " = ", collapse = ", "), "\n")
            }
          }
          print_matrix_overlap(x$overlap_ratio_goric_weights, "Ratio-of-GORIC(A)-weights (rgw)")
          print_matrix_overlap(x$overlap_ratio_ll_weights, "Ratio-of-log-likelihood-weights (rlw)")
          print_matrix_overlap(x$overlap_difLL, "Differences in log-likelihood values (ld)")
        }, nchar(text_overlap), text_color = blue, reset = reset
      )
    }
  }

  return(invisible(x))
}
