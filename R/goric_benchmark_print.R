print.benchmark <- function(x, output_type = c("rgw", "gw", "lw", "rlw", "ld",
                                               "rgw_log", "rlw_log", "all"),
                            hypo_rate_threshold = 1, color = TRUE,
                            percentiles = NULL, ...) {

  # Ensure the object is of class 'benchmark_means'
  if (!inherits(x, "benchmark")) {
    stop("Invalid object. The object should be of class 'benchmark'.", call. = FALSE)
  }

  output_type <- match.arg(output_type, c("rgw", "gw", "lw", "rlw", "ld",
                                          "rgw_log", "rlw_log", "all"))
  
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
  text_lw  <- paste0("Benchmark: Percentiles of ", orange, "Log-likelihood Weights", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_rgw <- paste0("Benchmark: Percentiles of ", orange, "Ratio-of-", goric_type, "-weights", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_rlw <- paste0("Benchmark: Percentiles of ", orange, "Ratio-of-log-likelihood-weights", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_ld  <- paste0("Benchmark: Percentiles of ", orange, "Differences in Log-likelihood Values", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_rgw_log <- paste0("Benchmark: Percentiles of ", orange, "Log(Ratio-of-", goric_type, "-weights)", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_rlw_log <- paste0("Benchmark: Percentiles of ", orange, "Log(Ratio-of-log-likelihood-weights)", blue, " for the Preferred Hypothesis '", pref_hypo, "'")

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

  # Header for the overlap column below: "Overlap with <reference population>".
  # x$overlap_reference is the population the overlap values were computed
  # against (see determine_overlap_reference() in goric_benchmark_utilities.R)
  # -- "Observed" when that category is present, else the last pop_es/pop_est
  # population supplied (e.g. "PES_2"). Strip the "pop_es = "/"pop_est = "
  # prefix the same way print_rounded_es_value() does for row labels, so the
  # header reads "Overlap with PES_2" rather than "Overlap with pop_es = PES_2".
  overlap_header <- if (!is.null(x$overlap_reference)) {
    paste0("Overlap with ", gsub("^(pop_es|pop_est) = ", "", x$overlap_reference))
  } else {
    "Overlap with Observed"
  }

  # 'percentiles', if supplied, overrides which percentiles are shown in
  # every section below -- the "Sample" + percentile% columns are recomputed
  # from the raw bootstrap draws stored on the object (x$combined_values), the
  # same way plot.benchmark()'s own 'percentiles' argument already works (see
  # recompute_percentile_table() in goric_benchmark_utilities.R). This lets
  # you look at different percentiles than the ones 'quant' was set to back
  # when the (often expensive) benchmark object was computed, without
  # rerunning benchmark_means()/benchmark_asymp(). Leave at the default
  # (NULL) to keep using the percentiles the object was built with.
  recompute_section <- function(benchmarks_list, combined_list) {
    if (is.null(percentiles)) {
      return(benchmarks_list)
    }
    for (pop_es in names(benchmarks_list)) {
      benchmarks_list[[pop_es]] <- recompute_percentile_table(
        benchmarks_list[[pop_es]], combined_list[[pop_es]], percentiles
      )
    }
    benchmarks_list
  }

  # Print benchmarks for Goric-weights percentiles
  if ("all" %in% output_type || "gw" %in% output_type) {
    x$benchmarks_goric_weights <- recompute_section(x$benchmarks_goric_weights, x$combined_values$gw_combined)
    for (pop_es in names(x$benchmarks_goric_weights)) {
      x$benchmarks_goric_weights[[pop_es]] <- cbind(
        x$benchmarks_goric_weights[[pop_es]],
        x$percentile_goric_weights[[pop_es]],
        overlap_column(
          x$overlap_goric_weights, pop_es, nrow(x$benchmarks_goric_weights[[pop_es]])
        )
      )
      colnames(x$benchmarks_goric_weights[[pop_es]])[
        ncol(x$benchmarks_goric_weights[[pop_es]])
      ] <- overlap_header
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

  # Print benchmarks for (unpenalized) log-likelihood-weight percentiles --
  # was computed all along (used internally by check_iter_adequacy()'s
  # too-low/not-converged diagnostics) but never exposed to the user or
  # printed here; only 'gw' (the GORIC(A) weight) had its own section.
  if ("all" %in% output_type || "lw" %in% output_type) {
    x$benchmarks_ll_weights <- recompute_section(x$benchmarks_ll_weights, x$combined_values$lw_combined)
    for (pop_es in names(x$benchmarks_ll_weights)) {
      x$benchmarks_ll_weights[[pop_es]] <- cbind(
        x$benchmarks_ll_weights[[pop_es]],
        x$percentile_ll_weights[[pop_es]],
        overlap_column(
          x$overlap_ll_weights, pop_es, nrow(x$benchmarks_ll_weights[[pop_es]])
        )
      )
      colnames(x$benchmarks_ll_weights[[pop_es]])[
        ncol(x$benchmarks_ll_weights[[pop_es]])
      ] <- overlap_header
    }
    print_section(
      text_lw,
      function() {
        for (pop_es in names(x$benchmarks_ll_weights)) {
          print_rounded_es_value(x$benchmarks_ll_weights[[pop_es]], pop_es,
                                 model_type, green, reset)
        }
      }, nchar(text_lw), text_color = blue, reset = reset
    )
  }

  # Print benchmarks for ratio Goric-weights percentiles
  if ("all" %in% output_type || "rgw" %in% output_type) {
    x$benchmarks_ratio_goric_weights <- recompute_section(x$benchmarks_ratio_goric_weights, x$combined_values$rgw_combined)

    # Loop over alle sets in both benchmarks_ratio_goric_weights and hypothesis_rate
    for (pop_es_name in names(x$benchmarks_ratio_goric_weights)) {
      x$benchmarks_ratio_goric_weights[[pop_es_name]] <- cbind(
        x$benchmarks_ratio_goric_weights[[pop_es_name]],
        x$percentile_ratio_goric_weights[[pop_es_name]],
        hypothesis_rate = x$hypothesis_rate[[pop_es_name]], # x$hypothesis_rate
        overlap_column(
          x$overlap_ratio_goric_weights, pop_es_name,
          nrow(x$benchmarks_ratio_goric_weights[[pop_es_name]])
        )
      )
      colnames(x$benchmarks_ratio_goric_weights[[pop_es_name]])[
        ncol(x$benchmarks_ratio_goric_weights[[pop_es_name]])
      ] <- overlap_header
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

  # Print benchmarks for log(ratio-GORIC(A)-weights) percentiles -- same
  # information as ratio-GORIC(A)-weights (rgw) above (it's just log(rgw)),
  # but on a scale-invariant, symmetric-around-0 scale that is better suited
  # to judging whether the preferred hypothesis and an alternative are
  # equally supported (rgw itself lives on a heavily right-skewed (0, Inf)
  # scale where '1' isn't a natural visual center; see the note above
  # get_results_benchmark()'s rgw_log/rlw_log computation in
  # goric_benchmark_utilities.R for the full rationale). No hypothesis_rate
  # column here, matching rlw/ld's (not rgw's) simpler layout.
  if ("all" %in% output_type || "rgw_log" %in% output_type) {
    x$benchmarks_ratio_goric_weights_log <- recompute_section(x$benchmarks_ratio_goric_weights_log, x$combined_values$rgw_log_combined)
    for (pop_es in names(x$benchmarks_ratio_goric_weights_log)) {
      x$benchmarks_ratio_goric_weights_log[[pop_es]] <- cbind(
        x$benchmarks_ratio_goric_weights_log[[pop_es]],
        x$percentile_ratio_goric_weights_log[[pop_es]],
        overlap_column(
          x$overlap_ratio_goric_weights_log, pop_es,
          nrow(x$benchmarks_ratio_goric_weights_log[[pop_es]])
        )
      )
      colnames(x$benchmarks_ratio_goric_weights_log[[pop_es]])[
        ncol(x$benchmarks_ratio_goric_weights_log[[pop_es]])
      ] <- overlap_header
    }
    print_section(
      text_rgw_log,
      function() {
        for (pop_es in names(x$benchmarks_ratio_goric_weights_log)) {
          print_rounded_es_value(x$benchmarks_ratio_goric_weights_log[[pop_es]], pop_es,
                                 model_type, green, reset)
        }
      }, nchar(text_rgw_log), text_color = blue, reset = reset
    )
  }

  # Print benchmarks for ratio log-likelihood-weights percentiles
  if ("all" %in% output_type || "rlw" %in% output_type) {
    x$benchmarks_ratio_ll_weights <- recompute_section(x$benchmarks_ratio_ll_weights, x$combined_values$rlw_combined)
    for (pop_es in names(x$benchmarks_ratio_ll_weights)) {
      x$benchmarks_ratio_ll_weights[[pop_es]] <- cbind(
        x$benchmarks_ratio_ll_weights[[pop_es]],
        x$percentile_ratio_ll_weights[[pop_es]],
        overlap_column(
          x$overlap_ratio_ll_weights, pop_es, nrow(x$benchmarks_ratio_ll_weights[[pop_es]])
        )
      )
      colnames(x$benchmarks_ratio_ll_weights[[pop_es]])[
        ncol(x$benchmarks_ratio_ll_weights[[pop_es]])
      ] <- overlap_header
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

  # Print benchmarks for log(ratio-log-likelihood-weights) percentiles --
  # same idea as rgw_log above, but for lw/rlw.
  if ("all" %in% output_type || "rlw_log" %in% output_type) {
    x$benchmarks_ratio_ll_weights_log <- recompute_section(x$benchmarks_ratio_ll_weights_log, x$combined_values$rlw_log_combined)
    for (pop_es in names(x$benchmarks_ratio_ll_weights_log)) {
      x$benchmarks_ratio_ll_weights_log[[pop_es]] <- cbind(
        x$benchmarks_ratio_ll_weights_log[[pop_es]],
        x$percentile_ratio_ll_weights_log[[pop_es]],
        overlap_column(
          x$overlap_ratio_ll_weights_log, pop_es,
          nrow(x$benchmarks_ratio_ll_weights_log[[pop_es]])
        )
      )
      colnames(x$benchmarks_ratio_ll_weights_log[[pop_es]])[
        ncol(x$benchmarks_ratio_ll_weights_log[[pop_es]])
      ] <- overlap_header
    }
    print_section(
      text_rlw_log,
      function() {
        for (pop_es in names(x$benchmarks_ratio_ll_weights_log)) {
          print_rounded_es_value(x$benchmarks_ratio_ll_weights_log[[pop_es]], pop_es,
                                 model_type, green, reset)
        }
      }, nchar(text_rlw_log), text_color = blue, reset = reset
    )
  }

  # Print benchmarks for difference log-likelihood-values percentiles
  if ("all" %in% output_type || "ld" %in% output_type) {
    x$benchmarks_difLL <- recompute_section(x$benchmarks_difLL, x$combined_values$ld_combined)
    for (pop_es in names(x$benchmarks_difLL)) {
      x$benchmarks_difLL[[pop_es]] <- cbind(
        x$benchmarks_difLL[[pop_es]],
        x$percentile_difLL[[pop_es]],
        overlap_column(
          x$overlap_difLL, pop_es, nrow(x$benchmarks_difLL[[pop_es]])
        )
      )
      colnames(x$benchmarks_difLL[[pop_es]])[
        ncol(x$benchmarks_difLL[[pop_es]])
      ] <- overlap_header
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

  return(invisible(x))
}
