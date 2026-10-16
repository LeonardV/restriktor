print.benchmark <- function(x, output_type = c("rgw", "gw", "lw", "rlw", "ld",
                                               "rgw_log", "rlw_log", "all"),
                            hypo_rate_threshold = NULL,
                            threshold_rlw = NULL, color = TRUE,
                            percentiles = NULL, ...) {

  # Ensure the object is of class 'benchmark_means'
  if (!inherits(x, "benchmark")) {
    stop("Invalid object. The object should be of class 'benchmark'.", call. = FALSE)
  }

  output_type <- match.arg(output_type, c("rgw", "gw", "lw", "rlw", "ld",
                                          "rgw_log", "rlw_log", "all"))

  # x$benchmarks/x$pctl_Sample/x$pctl_medianRefPop/x$overlap are nested list
  # fields (e.g. x$benchmarks$goric_weights) since the output was
  # restructured this way -- an older benchmark_means()/benchmark_asymp()
  # object (from before that change -- still held in an R session, or
  # reloaded from a saved .RData/.rds file) instead has flat fields like
  # x$benchmarks_goric_weights, so x$benchmarks$goric_weights silently
  # evaluates to NULL for it (R's '$' does not error on a missing list
  # field). Left unchecked that would print every section as empty rather
  # than fail loudly, so catch it here instead (see plot.benchmark()'s
  # matching check for the full rationale).
  if (is.null(x$benchmarks) || is.null(x$benchmarks$goric_weights)) {
    stop("\nrestriktor ERROR: 'x' does not have the expected x$benchmarks$goric_weights ",
         "field. This usually means 'x' is a benchmark object computed with an older ",
         "version of benchmark_means()/benchmark_asymp() (from before its output was ",
         "restructured into nested benchmarks/pctl_Sample/pctl_medianRefPop/overlap ",
         "lists) -- e.g. one still held in your R session, or reloaded from a saved ",
         ".RData/.rds file. Please rerun benchmark_means()/benchmark_asymp() to get a ",
         "fresh benchmark object, then call print() on that.", call. = FALSE)
  }

  ldots <- list(...)

  # Threshold q for the printed hypothesis rate. Defaults (NULL) to whatever
  # threshold the object itself was built with -- x$hypo_rate_threshold,
  # set via benchmark_means()/benchmark_asymp()'s own 'hypo_rate_threshold'
  # argument, falling back to 1 for objects built before that argument
  # existed (x$hypo_rate_threshold missing) -- so a plain print(x) matches
  # x$hypothesis_rate without having to repeat the threshold. Passing
  # hypo_rate_threshold explicitly here still overrides it for just this
  # print call, without rerunning the (often expensive) bootstrap, as shown
  # in the Guidelines vignette.
  if (is.null(hypo_rate_threshold)) {
    hypo_rate_threshold <- if (!is.null(x$hypo_rate_threshold)) x$hypo_rate_threshold else 1
  }
  # Same idea, but for the rlw-based 'rate_rlw' column (own
  # threshold, independent of hypo_rate_threshold.
  if (is.null(threshold_rlw)) {
    threshold_rlw <- if (!is.null(x$threshold_rlw)) x$threshold_rlw else 1
  }

  # compute hypothesis rate
  # x$hypothesis_rate already exists on the object itself (benchmark_means()/
  # benchmark_asymp() now store it via get_results_benchmark(), computed at
  # x$hypo_rate_threshold), but it's recomputed here regardless so that a
  # caller-supplied hypo_rate_threshold different from that is honored for
  # this print call -- this local recomputation is not saved back onto the
  # original object.
  x$hypothesis_rate <- lapply(x$combined_values$rgw_combined,
                              calculate_hypothesis_rate, q = hypo_rate_threshold)
  # Same recomputation, for the rlw-based rate_rlw column.
  x$rate_rlw <- lapply(x$combined_values$rlw_combined,
                                 calculate_hypothesis_rate, q = threshold_rlw)

  # no effect names
  # Was c("pop_est = No-effect", "pop_es = 0") -- but benchmark_means() names
  # this category "pop_es = No-effect" (matching benchmark_asymp's
  # "pop_est = No-effect"), never "pop_es = 0", so this never matched for
  # benchmark_means objects and the hypothesis_rate column was never hidden
  # for their No-effect category.
  #
  # This NA-ing (and the column-drop further below) only ever affects what
  # gets printed/returned by THIS print() call's local copy of x -- the
  # object's own x$hypothesis_rate field (set by benchmark_means()/
  # benchmark_asymp()) keeps the real computed value for "No-effect" too,
  # for a user who wants it without it being printed.
  NE_names <- names(x$hypothesis_rate) %in% c("pop_est = No-effect", "pop_es = No-effect")
  x$hypothesis_rate[NE_names] <- as.numeric(NA)
  # rate_rlw is hidden for "No-effect" for the same reason.
  x$rate_rlw[NE_names] <- as.numeric(NA)
  
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
  # "(rgw)" spells out the abbreviation used elsewhere for this same
  # quantity -- e.g. the "Hypothesis rate (rgw > <threshold>)" column header
  # further down this section -- so it's clear at a glance what 'rgw'
  # refers to. Kept inside the orange-colored span (same colors as before,
  # just extended slightly) since it's clarifying the highlighted weight
  # name, not a separate piece of the sentence.
  text_rgw <- paste0("Benchmark: Percentiles of ", orange, "Ratio-of-", goric_type, "-weights (rgw)", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_rlw <- paste0("Benchmark: Percentiles of ", orange, "Ratio-of-log-likelihood-weights", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_ld  <- paste0("Benchmark: Percentiles of ", orange, "Differences in Log-likelihood Values", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_rgw_log <- paste0("Benchmark: Percentiles of ", orange, "Log(Ratio-of-", goric_type, "-weights)", blue, " for the Preferred Hypothesis '", pref_hypo, "'")
  text_rlw_log <- paste0("Benchmark: Percentiles of ", orange, "Log(Ratio-of-log-likelihood-weights)", blue, " for the Preferred Hypothesis '", pref_hypo, "'")

  # Sub-headers for the two blocks every output_type section below is split
  # into (see combine_population_tables()/print_block1_table()/
  # print_block2_table()/print_block_header() in
  # goric_benchmark_utilities.R): block 1 is the raw "Sample" value plus the
  # percentile-distribution columns, as ONE shared table across every
  # population (with "Population" -- shown once per population -- and
  # "Hypotheses comparison" as leading columns); block 2 is the
  # percentile-rank comparisons (where the Sample value and the reference
  # population's median fall within those distributions), the overlap with
  # the reference population, and the hypothesis/rate column, in the same
  # combined layout.
  block1_header <- "Percentiles distributions"
  block2_header <- "Percentile & overlap comparisons"

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
  # x$overlap$reference is the population the overlap values were computed
  # against (see determine_overlap_reference() in goric_benchmark_utilities.R)
  # -- "Observed" when that category is present, else the last pop_es/pop_est
  # population supplied (e.g. "PES_2"). Strip the "pop_es = "/"pop_est = "
  # prefix the same way print_rounded_es_value() does for row labels, so the
  # header reads "Overlap with PES_2" rather than "Overlap with pop_es = PES_2".
  overlap_header <- if (!is.null(x$overlap$reference)) {
    paste0("Overlap with ", gsub("^(pop_es|pop_est) = ", "", x$overlap$reference))
  } else {
    "Overlap with Observed"
  }

  # Reference population used both for the "Overlap with ..." header above
  # and for combine_population_tables() below (its "(ref. pop.)" suffix on
  # the printed "Population" column).
  reference_pop_name <- x$overlap$reference

  # Raw (not percentile) value of the reference population's own median, for
  # the "Median ref. pop. (<value>)" column-header annotation in block 2 --
  # see print_block2_table()'s own comment for exactly when this is shown.
  # Recomputed fresh from x$combined_values (the raw bootstrap draws) rather
  # than rerunning the bootstrap. Same reference_pop_name used for the
  # "Overlap with ..." header above.
  get_median_ref_raw <- function(combined_list) {
    if (is.null(combined_list) || is.null(reference_pop_name) ||
        !(reference_pop_name %in% names(combined_list))) {
      return(NULL)
    }
    v <- combined_list[[reference_pop_name]]
    if (is.null(v)) return(NULL)
    if (is.null(dim(v))) {
      median(v, na.rm = TRUE)
    } else {
      apply(v, 2, median, na.rm = TRUE)
    }
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
    x$benchmarks$goric_weights <- recompute_section(x$benchmarks$goric_weights, x$combined_values$gw_combined)
    median_ref_raw_gw <- get_median_ref_raw(x$combined_values$gw_combined)
    for (pop_es in names(x$benchmarks$goric_weights)) {
      x$benchmarks$goric_weights[[pop_es]] <- cbind(
        x$benchmarks$goric_weights[[pop_es]],
        x$pctl_Sample$goric_weights[[pop_es]],
        x$pctl_medianRefPop$goric_weights[[pop_es]],
        overlap_column(
          x$overlap$goric_weights, pop_es, nrow(x$benchmarks$goric_weights[[pop_es]])
        )
      )
      colnames(x$benchmarks$goric_weights[[pop_es]])[
        ncol(x$benchmarks$goric_weights[[pop_es]])
      ] <- overlap_header
    }
    combined_gw <- combine_population_tables(x$benchmarks$goric_weights, reference_pop_name)
    print_section(
      text_gw,
      function() {
        cat("\n")
        print_block_header(block1_header, blue, reset)
        print_block1_table(combined_gw)
        print_block_header(block2_header, blue, reset)
        print_block2_table(combined_gw, median_ref_raw = median_ref_raw_gw)
      }, nchar(text_gw), text_color = blue, reset = reset
    )
  }

  # Print benchmarks for (unpenalized) log-likelihood-weight percentiles --
  # was computed all along (used internally by check_iter_adequacy()'s
  # too-low/not-converged diagnostics) but never exposed to the user or
  # printed here; only 'gw' (the GORIC(A) weight) had its own section.
  if ("all" %in% output_type || "lw" %in% output_type) {
    x$benchmarks$ll_weights <- recompute_section(x$benchmarks$ll_weights, x$combined_values$lw_combined)
    median_ref_raw_lw <- get_median_ref_raw(x$combined_values$lw_combined)
    for (pop_es in names(x$benchmarks$ll_weights)) {
      x$benchmarks$ll_weights[[pop_es]] <- cbind(
        x$benchmarks$ll_weights[[pop_es]],
        x$pctl_Sample$ll_weights[[pop_es]],
        x$pctl_medianRefPop$ll_weights[[pop_es]],
        overlap_column(
          x$overlap$ll_weights, pop_es, nrow(x$benchmarks$ll_weights[[pop_es]])
        )
      )
      colnames(x$benchmarks$ll_weights[[pop_es]])[
        ncol(x$benchmarks$ll_weights[[pop_es]])
      ] <- overlap_header
    }
    combined_lw <- combine_population_tables(x$benchmarks$ll_weights, reference_pop_name)
    print_section(
      text_lw,
      function() {
        cat("\n")
        print_block_header(block1_header, blue, reset)
        print_block1_table(combined_lw)
        print_block_header(block2_header, blue, reset)
        print_block2_table(combined_lw, median_ref_raw = median_ref_raw_lw)
      }, nchar(text_lw), text_color = blue, reset = reset
    )
  }

  # Print benchmarks for ratio Goric-weights percentiles
  if ("all" %in% output_type || "rgw" %in% output_type) {
    x$benchmarks$ratio_goric_weights <- recompute_section(x$benchmarks$ratio_goric_weights, x$combined_values$rgw_combined)
    median_ref_raw_rgw <- get_median_ref_raw(x$combined_values$rgw_combined)

    # Loop over alle sets in both benchmarks_ratio_goric_weights and hypothesis_rate
    for (pop_es_name in names(x$benchmarks$ratio_goric_weights)) {
      x$benchmarks$ratio_goric_weights[[pop_es_name]] <- cbind(
        x$benchmarks$ratio_goric_weights[[pop_es_name]],
        x$pctl_Sample$ratio_goric_weights[[pop_es_name]],
        x$pctl_medianRefPop$ratio_goric_weights[[pop_es_name]],
        hypothesis_rate = x$hypothesis_rate[[pop_es_name]], # x$hypothesis_rate
        overlap_column(
          x$overlap$ratio_goric_weights, pop_es_name,
          nrow(x$benchmarks$ratio_goric_weights[[pop_es_name]])
        )
      )
      colnames(x$benchmarks$ratio_goric_weights[[pop_es_name]])[
        ncol(x$benchmarks$ratio_goric_weights[[pop_es_name]])
      ] <- overlap_header
    }

    # Note: the "No-effect" population's hypothesis_rate column used to be
    # dropped outright here (it has no alternative hypothesis to compare
    # against, so the rate is meaningless there) -- it's no longer dropped:
    # x$hypothesis_rate[NE_names] was already set to NA above, and
    # combine_population_tables()/format_value() now render that NA as a
    # blank cell within the shared "Hypothesis rate" column instead, so
    # every population keeps the same column set for the combined table
    # below.
    combined_rgw <- combine_population_tables(x$benchmarks$ratio_goric_weights, reference_pop_name)
    print_section(
      text_rgw,
      function() {
        cat("\n")
        print_block_header(block1_header, blue, reset)
        print_block1_table(combined_rgw)
        print_block_header(block2_header, blue, reset)
        print_block2_table(combined_rgw, hypo_rate_threshold = hypo_rate_threshold,
                           median_ref_raw = median_ref_raw_rgw)
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
    x$benchmarks$ratio_goric_weights_log <- recompute_section(x$benchmarks$ratio_goric_weights_log, x$combined_values$rgw_log_combined)
    median_ref_raw_rgw_log <- get_median_ref_raw(x$combined_values$rgw_log_combined)
    for (pop_es in names(x$benchmarks$ratio_goric_weights_log)) {
      x$benchmarks$ratio_goric_weights_log[[pop_es]] <- cbind(
        x$benchmarks$ratio_goric_weights_log[[pop_es]],
        x$pctl_Sample$ratio_goric_weights_log[[pop_es]],
        x$pctl_medianRefPop$ratio_goric_weights_log[[pop_es]],
        overlap_column(
          x$overlap$ratio_goric_weights_log, pop_es,
          nrow(x$benchmarks$ratio_goric_weights_log[[pop_es]])
        )
      )
      colnames(x$benchmarks$ratio_goric_weights_log[[pop_es]])[
        ncol(x$benchmarks$ratio_goric_weights_log[[pop_es]])
      ] <- overlap_header
    }
    combined_rgw_log <- combine_population_tables(x$benchmarks$ratio_goric_weights_log, reference_pop_name)
    print_section(
      text_rgw_log,
      function() {
        cat("\n")
        print_block_header(block1_header, blue, reset)
        print_block1_table(combined_rgw_log)
        print_block_header(block2_header, blue, reset)
        print_block2_table(combined_rgw_log, median_ref_raw = median_ref_raw_rgw_log)
      }, nchar(text_rgw_log), text_color = blue, reset = reset
    )
  }

  # Print benchmarks for ratio log-likelihood-weights percentiles
  if ("all" %in% output_type || "rlw" %in% output_type) {
    x$benchmarks$ratio_ll_weights <- recompute_section(x$benchmarks$ratio_ll_weights, x$combined_values$rlw_combined)
    median_ref_raw_rlw <- get_median_ref_raw(x$combined_values$rlw_combined)
    for (pop_es in names(x$benchmarks$ratio_ll_weights)) {
      x$benchmarks$ratio_ll_weights[[pop_es]] <- cbind(
        x$benchmarks$ratio_ll_weights[[pop_es]],
        x$pctl_Sample$ratio_ll_weights[[pop_es]],
        x$pctl_medianRefPop$ratio_ll_weights[[pop_es]],
        rate_rlw = x$rate_rlw[[pop_es]],
        overlap_column(
          x$overlap$ratio_ll_weights, pop_es, nrow(x$benchmarks$ratio_ll_weights[[pop_es]])
        )
      )
      colnames(x$benchmarks$ratio_ll_weights[[pop_es]])[
        ncol(x$benchmarks$ratio_ll_weights[[pop_es]])
      ] <- overlap_header
    }
    # See the matching comment in the "rgw" section above -- rate_rlw is no
    # longer dropped for the "No-effect" population; its NA value now prints
    # blank within the shared "Rate" column instead.
    combined_rlw <- combine_population_tables(x$benchmarks$ratio_ll_weights, reference_pop_name)
    print_section(
      text_rlw,
      function() {
        cat("\n")
        print_block_header(block1_header, blue, reset)
        print_block1_table(combined_rlw)
        print_block_header(block2_header, blue, reset)
        print_block2_table(combined_rlw, threshold_rlw = threshold_rlw,
                           median_ref_raw = median_ref_raw_rlw)
      }, nchar(text_rlw), text_color = blue, reset = reset
    )
  }

  # Print benchmarks for log(ratio-log-likelihood-weights) percentiles --
  # same idea as rgw_log above, but for lw/rlw.
  if ("all" %in% output_type || "rlw_log" %in% output_type) {
    x$benchmarks$ratio_ll_weights_log <- recompute_section(x$benchmarks$ratio_ll_weights_log, x$combined_values$rlw_log_combined)
    median_ref_raw_rlw_log <- get_median_ref_raw(x$combined_values$rlw_log_combined)
    for (pop_es in names(x$benchmarks$ratio_ll_weights_log)) {
      x$benchmarks$ratio_ll_weights_log[[pop_es]] <- cbind(
        x$benchmarks$ratio_ll_weights_log[[pop_es]],
        x$pctl_Sample$ratio_ll_weights_log[[pop_es]],
        x$pctl_medianRefPop$ratio_ll_weights_log[[pop_es]],
        overlap_column(
          x$overlap$ratio_ll_weights_log, pop_es,
          nrow(x$benchmarks$ratio_ll_weights_log[[pop_es]])
        )
      )
      colnames(x$benchmarks$ratio_ll_weights_log[[pop_es]])[
        ncol(x$benchmarks$ratio_ll_weights_log[[pop_es]])
      ] <- overlap_header
    }
    combined_rlw_log <- combine_population_tables(x$benchmarks$ratio_ll_weights_log, reference_pop_name)
    print_section(
      text_rlw_log,
      function() {
        cat("\n")
        print_block_header(block1_header, blue, reset)
        print_block1_table(combined_rlw_log)
        print_block_header(block2_header, blue, reset)
        print_block2_table(combined_rlw_log, median_ref_raw = median_ref_raw_rlw_log)
      }, nchar(text_rlw_log), text_color = blue, reset = reset
    )
  }

  # Print benchmarks for difference log-likelihood-values percentiles
  if ("all" %in% output_type || "ld" %in% output_type) {
    x$benchmarks$difLL <- recompute_section(x$benchmarks$difLL, x$combined_values$ld_combined)
    median_ref_raw_ld <- get_median_ref_raw(x$combined_values$ld_combined)
    for (pop_es in names(x$benchmarks$difLL)) {
      x$benchmarks$difLL[[pop_es]] <- cbind(
        x$benchmarks$difLL[[pop_es]],
        x$pctl_Sample$difLL[[pop_es]],
        x$pctl_medianRefPop$difLL[[pop_es]],
        overlap_column(
          x$overlap$difLL, pop_es, nrow(x$benchmarks$difLL[[pop_es]])
        )
      )
      colnames(x$benchmarks$difLL[[pop_es]])[
        ncol(x$benchmarks$difLL[[pop_es]])
      ] <- overlap_header
    }
    combined_ld <- combine_population_tables(x$benchmarks$difLL, reference_pop_name)
    print_section(
      text_ld,
      function() {
        cat("\n")
        print_block_header(block1_header, blue, reset)
        print_block1_table(combined_ld)
        print_block_header(block2_header, blue, reset)
        print_block2_table(combined_ld, median_ref_raw = median_ref_raw_ld)
      }, nchar(text_ld), text_color = blue, reset = reset
    )
  }

  # Caution about the hypothesis rate (rgw) and the rlw rate for the
  # 'Observed' population (shown in the rgw and rlw sections): these are
  # proportions of draws generated around the sample estimates (i.e., a
  # parametric bootstrap) in which the preferred hypothesis is supported.
  # Preacher and Merkle (2012) showed that such bootstrap 'ranking
  # proportions' are unreliable as estimates of how often a replication
  # would prefer the same hypothesis: they depend strongly on the specific
  # sample, especially when the IC values of the hypotheses are close; only
  # when (nearly) all draws agree are they trustworthy. In contrast,
  # intervals/percentiles of the bootstrap distributions of (differences in)
  # IC values do perform well. There is no established cutoff for 'close', so
  # this is a standard note; the rate itself is the natural check (a rate
  # not close to 0 or 1 indicates that the IC values are close).
  # For other (user-specified) populations, the draws are generated around
  # those population values, so the rate is then a genuine simulation result
  # (the probability of supporting the hypothesis if that population were
  # true), and the caution does not apply.
  rate_shown <- any(c("all", "rgw", "rlw") %in% output_type)
  obs_present <- any(grepl("Observed", names(x$hypothesis_rate), fixed = TRUE))
  if (rate_shown && obs_present) {
    message("---\nNote: The hypothesis rate and the rlw rate for the 'Observed' population are",
            " the proportions of draws, generated around the sample estimates, in which the",
            " preferred hypothesis is supported (ratio of weights > threshold). When the",
            " hypotheses in a row are close, that is, when the ratio of weights is close to the",
            " threshold, these rates are unreliable: they depend strongly on the specific sample",
            " (Preacher & Merkle, 2012) and are not the probability that a replication would",
            " prefer the same hypothesis. A rate that is not close to 0 or 1 signals this for",
            " that row. The percentiles of the benchmark distributions are not affected, and",
            " neither are the rates for other (specified) populations.")
  }

  return(invisible(x))
}
