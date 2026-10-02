# [CHANGE 2026-10 | Rebecca] new output types lw/rgw_log/rlw_log; thresholds default NULL (taken from object); new args threshold_rlw, percentiles
print.benchmark <- function(x, output_type = c("rgw", "gw", "lw", "rlw", "ld",
                                               "rgw_log", "rlw_log", "all"),
                            hypo_rate_threshold = NULL,
                            threshold_rlw = NULL, color = TRUE,
                            percentiles = NULL, ...) {
# [/CHANGE 2026-10]

  # Ensure the object is of class 'benchmark_means'
  if (!inherits(x, "benchmark")) {
    stop("Invalid object. The object should be of class 'benchmark'.", call. = FALSE)
  }

  # [CHANGE 2026-10 | Rebecca] match.arg incl. new output types
  output_type <- match.arg(output_type, c("rgw", "gw", "lw", "rlw", "ld",
                                          "rgw_log", "rlw_log", "all"))

  # [CHANGE 2026-10 | Rebecca] informative error for benchmark objects with the old flat field names
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
  # [/CHANGE 2026-10]

  ldots <- list(...)

  # [CHANGE 2026-10 | Rebecca] hypo_rate_threshold / threshold_rlw default to the values stored in the object
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
  # threshold, independent of hypo_rate_threshold -- see
  # benchmark_means()'s 'threshold_rlw' for the full rationale).
  if (is.null(threshold_rlw)) {
    threshold_rlw <- if (!is.null(x$threshold_rlw)) x$threshold_rlw else 1
  }

  # [/CHANGE 2026-10]
  # compute hypothesis rate
  # [CHANGE 2026-10 | Rebecca] comment on recomputation; rate_rlw (rlw-based rate) recomputed with threshold_rlw
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
  # [/CHANGE 2026-10]

  # no effect names
  # [CHANGE 2026-10 | Rebecca] No-effect category is named "pop_es = No-effect" (was "pop_es = 0"); hide rate_rlw for No-effect too
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
  # [/CHANGE 2026-10]
  
  # [CHANGE 2026-10 | audit] A9: failed draws from x$n_failed_draws; requested vs. successful draws (x$iter counts successful draws)
  # number of failed bootstrap draws per population (goric() error in that
  # draw, see parallel_function_asymp()); x$n_failed_draws for a benchmark
  # object that stores it, else (older object) the "empty_lists_count"
  # attribute of the combined draws
  empty_lists_count <- if (!is.null(x$n_failed_draws)) {
    x$n_failed_draws
  } else {
    sapply(x$combined_values$gw_combined, function(x) {
      attr(x, "empty_lists_count") })
  }
  # number of requested draws (x$iter counts the successful draws)
  iter_requested <- if (!is.null(x$iter_requested)) x$iter_requested else max(x$iter)
  # [/CHANGE 2026-10]
  
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
    # [CHANGE 2026-10 | audit] B27: Cohen's f of the observed means under alt_group_size
    cohens_f_alt_group_size <- x$cohens_f_alt_group_size
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
  # [CHANGE 2026-10 | Rebecca] new section headers lw, rgw_log, rlw_log; "(rgw)" added to the rgw header
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
  # [/CHANGE 2026-10]

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
    # [CHANGE 2026-10 | audit] B27: print Cohen's f of the observed means w.r.t. alt_group_size separately
    # with (unequal) alternative group sizes, the observed means have another
    # Cohen's f w.r.t. those group sizes -- see benchmark_means()
    if (!is.null(cohens_f_alt_group_size) &&
        !isTRUE(all.equal(cohens_f_alt_group_size, cohens_f_observed))) {
      cat(sprintf("Effect-Size (Cohens f) of the Observed Means with alt_group_size: %s%s%s\n",
                  green, sprintf("%.3f", cohens_f_alt_group_size), reset))
    }
    # [/CHANGE 2026-10]
  } else {
    print_formatted_matrix(formatted_values, green, reset)
  }
  #cat(strrep("-", 70), "\n")
  cat("\n")
  
# -------------------------------------------------------------------------
  # [CHANGE 2026-10 | audit] A9: requested draws (x$iter = successful draws)
  R <- iter_requested
  succesful_draws <- R - empty_lists_count
  if (any(succesful_draws < R)) { 
    cat("Number of requested bootstrap draws:", R, "\n")
    for (i in seq_along(succesful_draws)) {
      cat(sprintf("  Number of successful bootstrap draws for %*s: %s\n", 
                  10, names(succesful_draws)[i], succesful_draws[i]))
      # [CHANGE 2026-10 | audit] A9: report the first error message of failed draws
      # why the draws failed (first error message per population), if stored
      if (!is.null(x$draw_errors) && !is.na(x$draw_errors[i]) && succesful_draws[i] < R) {
        cat(sprintf("    (first error: %s)\n", trimws(x$draw_errors[i])))
      }
      # [/CHANGE 2026-10]
    }
    text_msg <- paste("Advise: If a substantial number of bootstrap draws fail to converge,", 
                      "it is advisable to increase the number of bootstrap iterations.")
    message("---\n", text_msg)
  }
  # [CHANGE 2026-10 | audit] A9: report draws with a muffled goric() warning (draws kept)
  # draws in which goric() gave a (muffled) warning; these draws were kept
  if (!is.null(x$n_warned_draws) && any(x$n_warned_draws > 0)) {
    cat("Number of bootstrap draws with a (muffled) goric() warning:",
        paste0(names(x$n_warned_draws), ": ", x$n_warned_draws, collapse = "; "), "\n")
    first_w <- x$draw_warnings[!is.na(x$draw_warnings)][1]
    cat("  (these draws were kept; first warning:", trimws(first_w), ")\n")
  }
  # [/CHANGE 2026-10]

  # Normalize output_type to lowercase
  output_type <- tolower(output_type)

  # [CHANGE 2026-10 | Rebecca] overlap column header named after the reference population (x$overlap$reference)
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
  # [/CHANGE 2026-10]
  # [CHANGE 2026-10 | audit] B24: label of the reference population; notes that means follow ratio_pop_means
  # Label appended to the reference population (see print_rounded_es_value()).
  # With ratio_pop_means, the default 'Observed' population has the observed
  # effect size, but its means follow the pattern in ratio_pop_means (not
  # the observed means) -- say so.
  reference_label <- if (inherits(x, "benchmark_means") && !is.null(x$ratio_pop_means) &&
                         !is.null(x$overlap$reference) &&
                         grepl("= Observed$", x$overlap$reference)) {
    "(Reference population; observed effect size, means from ratio_pop_means)"
  } else {
    "(Reference population)"
  }
  # [/CHANGE 2026-10]

  # [CHANGE 2026-10 | Rebecca] 'percentiles' argument: recompute percentile tables from the stored draws
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
        # [CHANGE 2026-10 | audit] N7: pass pref_hypo_name so the recomputed table is aligned by hypothesis name
        benchmarks_list[[pop_es]], combined_list[[pop_es]], percentiles,
        pref_hypo_name = pref_hypo
      )
    }
    benchmarks_list
  }
  # [/CHANGE 2026-10]

  # Print benchmarks for Goric-weights percentiles
  if ("all" %in% output_type || "gw" %in% output_type) {
    # [CHANGE 2026-10 | Rebecca] gw section: nested fields; add pctl_Sample, pctl_medianRefPop and overlap columns
    x$benchmarks$goric_weights <- recompute_section(x$benchmarks$goric_weights, x$combined_values$gw_combined)
    for (pop_es in names(x$benchmarks$goric_weights)) {
      # [CHANGE 2026-10 | audit] E5/O10: overlap column via overlap_column() (NA-safe, with notes)
      ov_col <- overlap_column(
        x$overlap$goric_weights, pop_es, nrow(x$benchmarks$goric_weights[[pop_es]])
      )
      x$benchmarks$goric_weights[[pop_es]] <- cbind(
        x$benchmarks$goric_weights[[pop_es]],
        x$pctl_Sample$goric_weights[[pop_es]],
        x$pctl_medianRefPop$goric_weights[[pop_es]],
        # [CHANGE 2026-10 | audit] E5/O10: overlap column via overlap_column()
        ov_col
      )
      colnames(x$benchmarks$goric_weights[[pop_es]])[
        ncol(x$benchmarks$goric_weights[[pop_es]])
      ] <- overlap_header
      # [CHANGE 2026-10 | audit] E5/O10: keep the reasons for NA overlaps as attribute
      # why an overlap is NA (if any), per row -- see overlap_column()
      attr(x$benchmarks$goric_weights[[pop_es]], "overlap_notes") <- attr(ov_col, "notes")
    }
    # [/CHANGE 2026-10]
    print_section(
      text_gw,
      function() {
        # [CHANGE 2026-10 | Rebecca] nested field x$benchmarks$goric_weights
        for (pop_es in names(x$benchmarks$goric_weights)) {
          print_rounded_es_value(x$benchmarks$goric_weights[[pop_es]], pop_es,
                                 model_type, green, reset,
                                 # [CHANGE 2026-10 | audit] B24: mark the reference population in the printed table
                                 is_reference = identical(pop_es, x$overlap$reference),
                                 reference_label = reference_label)
        }
      }, nchar(text_gw), text_color = blue, reset = reset
    )
  }

  # [CHANGE 2026-10 | Rebecca] new 'lw' section (log-likelihood weights)
  # Print benchmarks for (unpenalized) log-likelihood-weight percentiles --
  # was computed all along (used internally by check_iter_adequacy()'s
  # too-low/not-converged diagnostics) but never exposed to the user or
  # printed here; only 'gw' (the GORIC(A) weight) had its own section.
  if ("all" %in% output_type || "lw" %in% output_type) {
    x$benchmarks$ll_weights <- recompute_section(x$benchmarks$ll_weights, x$combined_values$lw_combined)
    for (pop_es in names(x$benchmarks$ll_weights)) {
      # [CHANGE 2026-10 | audit] E5/O10: overlap column via overlap_column() (NA-safe, with notes)
      ov_col <- overlap_column(
        x$overlap$ll_weights, pop_es, nrow(x$benchmarks$ll_weights[[pop_es]])
      )
      x$benchmarks$ll_weights[[pop_es]] <- cbind(
        x$benchmarks$ll_weights[[pop_es]],
        x$pctl_Sample$ll_weights[[pop_es]],
        x$pctl_medianRefPop$ll_weights[[pop_es]],
        # [CHANGE 2026-10 | audit] E5/O10: overlap column via overlap_column()
        ov_col
      )
      colnames(x$benchmarks$ll_weights[[pop_es]])[
        ncol(x$benchmarks$ll_weights[[pop_es]])
      ] <- overlap_header
      # [CHANGE 2026-10 | audit] E5/O10: keep the reasons for NA overlaps as attribute
      # why an overlap is NA (if any), per row -- see overlap_column()
      attr(x$benchmarks$ll_weights[[pop_es]], "overlap_notes") <- attr(ov_col, "notes")
    }
    print_section(
      text_lw,
      function() {
        for (pop_es in names(x$benchmarks$ll_weights)) {
          print_rounded_es_value(x$benchmarks$ll_weights[[pop_es]], pop_es,
                                 model_type, green, reset,
                                 # [CHANGE 2026-10 | audit] B24: mark the reference population in the printed table
                                 is_reference = identical(pop_es, x$overlap$reference),
                                 reference_label = reference_label)
        }
      }, nchar(text_lw), text_color = blue, reset = reset
    )
  }
  # [/CHANGE 2026-10]

  # Print benchmarks for ratio Goric-weights percentiles
  if ("all" %in% output_type || "rgw" %in% output_type) {
    # [CHANGE 2026-10 | Rebecca] rgw section: recompute percentiles if requested (nested field)
    x$benchmarks$ratio_goric_weights <- recompute_section(x$benchmarks$ratio_goric_weights, x$combined_values$rgw_combined)

    # Loop over alle sets in both benchmarks_ratio_goric_weights and hypothesis_rate
    # [CHANGE 2026-10 | Rebecca] rgw section: nested fields; add pctl_Sample, pctl_medianRefPop, hypothesis_rate and overlap columns
    for (pop_es_name in names(x$benchmarks$ratio_goric_weights)) {
      # [CHANGE 2026-10 | audit] E5/O10 + N7: overlap column via overlap_column(), aligned by hypothesis name
      ov_col <- overlap_column(
        x$overlap$ratio_goric_weights, pop_es_name,
        nrow(x$benchmarks$ratio_goric_weights[[pop_es_name]]),
        row_names = rownames(x$benchmarks$ratio_goric_weights[[pop_es_name]]),
        pref_hypo_name = pref_hypo
      )
      # [/CHANGE 2026-10]
      x$benchmarks$ratio_goric_weights[[pop_es_name]] <- cbind(
        x$benchmarks$ratio_goric_weights[[pop_es_name]],
        x$pctl_Sample$ratio_goric_weights[[pop_es_name]],
        x$pctl_medianRefPop$ratio_goric_weights[[pop_es_name]],
        # [CHANGE 2026-10 | audit] N7: hypothesis_rate aligned with the table rows by hypothesis name
        # aligned with the table's rows by hypothesis name
        hypothesis_rate = align_by_hypothesis(
          x$hypothesis_rate[[pop_es_name]],
          rownames(x$benchmarks$ratio_goric_weights[[pop_es_name]]), pref_hypo
        ),
        # [/CHANGE 2026-10]
        # [CHANGE 2026-10 | audit] E5/O10: overlap column via overlap_column()
        ov_col
      )
      colnames(x$benchmarks$ratio_goric_weights[[pop_es_name]])[
        ncol(x$benchmarks$ratio_goric_weights[[pop_es_name]])
      ] <- overlap_header
      # [CHANGE 2026-10 | audit] E5/O10: keep the reasons for NA overlaps as attribute
      # why an overlap is NA (if any), per row -- see overlap_column()
      attr(x$benchmarks$ratio_goric_weights[[pop_es_name]], "overlap_notes") <- attr(ov_col, "notes")
    # [/CHANGE 2026-10]
    }
    # TO DO nu bij No-effect ook hypothesis_rate maar die zouden we dacht ik niet meer laten zien omdat het verwarrend is wat het betekent
    # [CHANGE 2026-10 | Rebecca] note on the TO DO above
    # Is nu juist weer weg, als het goed is.
    #       daarnaast heeft anders echt beschrijving nodig, want het is steun info hypo onder NE en dan in gehele set?
    #       Ws overleggen of dit handig is - nu ineens denk ik dat het zo gek nog niet is :-).
    
    if (any(NE_names)) {
      # [CHANGE 2026-10 | Rebecca] comment: hypothesis_rate column selected by name instead of position
      # Was hardcoded as column -7 (assuming Sample + 5 default quantiles as
      # columns 1-6, hypothesis_rate as 7); now selected by name instead,
      # since adding the percentile column shifted hypothesis_rate to column 8
      # and a positional index would silently drop the wrong column.
      # [/CHANGE 2026-10]
      # [CHANGE 2026-10 | audit] N7: drop the hypothesis_rate column by name for No-effect, keeping the overlap_notes attribute
      ne_tab <- x$benchmarks$ratio_goric_weights[NE_names][[1]]
      ne_notes <- attr(ne_tab, "overlap_notes") # dropped by `[`, see overlap_column()
      ne_tab <- ne_tab[, colnames(ne_tab) != "hypothesis_rate", drop = FALSE]
      attr(ne_tab, "overlap_notes") <- ne_notes
      x$benchmarks$ratio_goric_weights[NE_names][[1]] <- ne_tab
      # [/CHANGE 2026-10]
    }
    
    print_section(
      text_rgw,
      function() {
        # [CHANGE 2026-10 | Rebecca] nested field; pass is_reference and hypo_rate_threshold to print_rounded_es_value()
        for (pop_es in names(x$benchmarks$ratio_goric_weights)) {
          print_rounded_es_value(x$benchmarks$ratio_goric_weights[[pop_es]], pop_es,
                                 model_type, green, reset,
                                 is_reference = identical(pop_es, x$overlap$reference),
                                 # [CHANGE 2026-10 | audit] B24: label of the reference population
                                 reference_label = reference_label,
                                 hypo_rate_threshold = hypo_rate_threshold)
        # [/CHANGE 2026-10]
        }
      }, nchar(text_rgw), text_color = blue, reset = reset
    )
  }

  # [CHANGE 2026-10 | Rebecca] new 'rgw_log' section (log of ratio of GORIC(A) weights)
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
    for (pop_es in names(x$benchmarks$ratio_goric_weights_log)) {
      # [CHANGE 2026-10 | audit] E5/O10: overlap column via overlap_column()
      ov_col <- overlap_column(
        x$overlap$ratio_goric_weights_log, pop_es,
        # [CHANGE 2026-10 | audit] N7: overlap aligned by hypothesis name
        nrow(x$benchmarks$ratio_goric_weights_log[[pop_es]]),
        row_names = rownames(x$benchmarks$ratio_goric_weights_log[[pop_es]]),
        pref_hypo_name = pref_hypo
      )
      # [CHANGE 2026-10 | audit] E5/O10: cbind with NA-safe overlap column (ov_col)
      x$benchmarks$ratio_goric_weights_log[[pop_es]] <- cbind(
        x$benchmarks$ratio_goric_weights_log[[pop_es]],
        x$pctl_Sample$ratio_goric_weights_log[[pop_es]],
        x$pctl_medianRefPop$ratio_goric_weights_log[[pop_es]],
        ov_col
      # [/CHANGE 2026-10]
      )
      colnames(x$benchmarks$ratio_goric_weights_log[[pop_es]])[
        ncol(x$benchmarks$ratio_goric_weights_log[[pop_es]])
      ] <- overlap_header
      # [CHANGE 2026-10 | audit] E5/O10: keep the reasons for NA overlaps as attribute
      # why an overlap is NA (if any), per row -- see overlap_column()
      attr(x$benchmarks$ratio_goric_weights_log[[pop_es]], "overlap_notes") <- attr(ov_col, "notes")
    }
    print_section(
      text_rgw_log,
      function() {
        for (pop_es in names(x$benchmarks$ratio_goric_weights_log)) {
          print_rounded_es_value(x$benchmarks$ratio_goric_weights_log[[pop_es]], pop_es,
                                 model_type, green, reset,
                                 # [CHANGE 2026-10 | audit] B24: mark the reference population in the printed table
                                 is_reference = identical(pop_es, x$overlap$reference),
                                 reference_label = reference_label)
        }
      }, nchar(text_rgw_log), text_color = blue, reset = reset
    )
  }
  # [/CHANGE 2026-10]

  # Print benchmarks for ratio log-likelihood-weights percentiles
  if ("all" %in% output_type || "rlw" %in% output_type) {
    # [CHANGE 2026-10 | Rebecca] rlw section: nested fields; add pctl_Sample, pctl_medianRefPop, rate_rlw and overlap columns; threshold_rlw passed to print
    x$benchmarks$ratio_ll_weights <- recompute_section(x$benchmarks$ratio_ll_weights, x$combined_values$rlw_combined)
    for (pop_es in names(x$benchmarks$ratio_ll_weights)) {
      # [CHANGE 2026-10 | audit] E5/O10 + N7: overlap column via overlap_column(), aligned by hypothesis name
      ov_col <- overlap_column(
        x$overlap$ratio_ll_weights, pop_es,
        nrow(x$benchmarks$ratio_ll_weights[[pop_es]]),
        row_names = rownames(x$benchmarks$ratio_ll_weights[[pop_es]]),
        pref_hypo_name = pref_hypo
      )
      # [/CHANGE 2026-10]
      x$benchmarks$ratio_ll_weights[[pop_es]] <- cbind(
        x$benchmarks$ratio_ll_weights[[pop_es]],
        x$pctl_Sample$ratio_ll_weights[[pop_es]],
        x$pctl_medianRefPop$ratio_ll_weights[[pop_es]],
        # [CHANGE 2026-10 | audit] N7: rate_rlw aligned with the table rows by hypothesis name
        # aligned with the table's rows by hypothesis name
        rate_rlw = align_by_hypothesis(
          x$rate_rlw[[pop_es]], rownames(x$benchmarks$ratio_ll_weights[[pop_es]]), pref_hypo
        ),
        # [/CHANGE 2026-10]
        # [CHANGE 2026-10 | audit] E5/O10: overlap column via overlap_column()
        ov_col
      )
      colnames(x$benchmarks$ratio_ll_weights[[pop_es]])[
        ncol(x$benchmarks$ratio_ll_weights[[pop_es]])
      ] <- overlap_header
      # [CHANGE 2026-10 | audit] E5/O10: keep the reasons for NA overlaps as attribute
      # why an overlap is NA (if any), per row -- see overlap_column()
      attr(x$benchmarks$ratio_ll_weights[[pop_es]], "overlap_notes") <- attr(ov_col, "notes")
    }
    if (any(NE_names)) {
      # [CHANGE 2026-10 | audit] N7: drop the rate_rlw column by name for No-effect, keeping the overlap_notes attribute
      ne_tab <- x$benchmarks$ratio_ll_weights[NE_names][[1]]
      ne_notes <- attr(ne_tab, "overlap_notes") # dropped by `[`, see overlap_column()
      ne_tab <- ne_tab[, colnames(ne_tab) != "rate_rlw", drop = FALSE]
      attr(ne_tab, "overlap_notes") <- ne_notes
      x$benchmarks$ratio_ll_weights[NE_names][[1]] <- ne_tab
      # [/CHANGE 2026-10]
    }
    print_section(
      text_rlw,
      function() {
        for (pop_es in names(x$benchmarks$ratio_ll_weights)) {
          print_rounded_es_value(x$benchmarks$ratio_ll_weights[[pop_es]], pop_es,
                                 model_type, green, reset,
                                 is_reference = identical(pop_es, x$overlap$reference),
                                 # [CHANGE 2026-10 | audit] B24: label of the reference population
                                 reference_label = reference_label,
                                 threshold_rlw = threshold_rlw)
    # [/CHANGE 2026-10]
        }
      }, nchar(text_rlw), text_color = blue, reset = reset
    )
  }

  # [CHANGE 2026-10 | Rebecca] new 'rlw_log' section (log of ratio of log-likelihood weights)
  # Print benchmarks for log(ratio-log-likelihood-weights) percentiles --
  # same idea as rgw_log above, but for lw/rlw.
  if ("all" %in% output_type || "rlw_log" %in% output_type) {
    x$benchmarks$ratio_ll_weights_log <- recompute_section(x$benchmarks$ratio_ll_weights_log, x$combined_values$rlw_log_combined)
    for (pop_es in names(x$benchmarks$ratio_ll_weights_log)) {
      # [CHANGE 2026-10 | audit] E5/O10: overlap column via overlap_column()
      ov_col <- overlap_column(
        x$overlap$ratio_ll_weights_log, pop_es,
        # [CHANGE 2026-10 | audit] N7: overlap aligned by hypothesis name
        nrow(x$benchmarks$ratio_ll_weights_log[[pop_es]]),
        row_names = rownames(x$benchmarks$ratio_ll_weights_log[[pop_es]]),
        pref_hypo_name = pref_hypo
      )
      # [CHANGE 2026-10 | audit] E5/O10: cbind with NA-safe overlap column (ov_col)
      x$benchmarks$ratio_ll_weights_log[[pop_es]] <- cbind(
        x$benchmarks$ratio_ll_weights_log[[pop_es]],
        x$pctl_Sample$ratio_ll_weights_log[[pop_es]],
        x$pctl_medianRefPop$ratio_ll_weights_log[[pop_es]],
        ov_col
      # [/CHANGE 2026-10]
      )
      colnames(x$benchmarks$ratio_ll_weights_log[[pop_es]])[
        ncol(x$benchmarks$ratio_ll_weights_log[[pop_es]])
      ] <- overlap_header
      # [CHANGE 2026-10 | audit] E5/O10: keep the reasons for NA overlaps as attribute
      # why an overlap is NA (if any), per row -- see overlap_column()
      attr(x$benchmarks$ratio_ll_weights_log[[pop_es]], "overlap_notes") <- attr(ov_col, "notes")
    }
    print_section(
      text_rlw_log,
      function() {
        for (pop_es in names(x$benchmarks$ratio_ll_weights_log)) {
          print_rounded_es_value(x$benchmarks$ratio_ll_weights_log[[pop_es]], pop_es,
                                 model_type, green, reset,
                                 # [CHANGE 2026-10 | audit] B24: mark the reference population in the printed table
                                 is_reference = identical(pop_es, x$overlap$reference),
                                 reference_label = reference_label)
        }
      }, nchar(text_rlw_log), text_color = blue, reset = reset
    )
  }
  # [/CHANGE 2026-10]

  # Print benchmarks for difference log-likelihood-values percentiles
  if ("all" %in% output_type || "ld" %in% output_type) {
    # [CHANGE 2026-10 | Rebecca] ld section: nested fields; add pctl_Sample, pctl_medianRefPop and overlap columns
    x$benchmarks$difLL <- recompute_section(x$benchmarks$difLL, x$combined_values$ld_combined)
    for (pop_es in names(x$benchmarks$difLL)) {
      # [CHANGE 2026-10 | audit] E5/O10 + N7: overlap column via overlap_column(), aligned by hypothesis name
      ov_col <- overlap_column(
        x$overlap$difLL, pop_es,
        nrow(x$benchmarks$difLL[[pop_es]]),
        row_names = rownames(x$benchmarks$difLL[[pop_es]]),
        pref_hypo_name = pref_hypo
      # [/CHANGE 2026-10]
      )
      # [CHANGE 2026-10 | audit] E5/O10: cbind with NA-safe overlap column (ov_col)
      x$benchmarks$difLL[[pop_es]] <- cbind(
        x$benchmarks$difLL[[pop_es]],
        x$pctl_Sample$difLL[[pop_es]],
        x$pctl_medianRefPop$difLL[[pop_es]],
        ov_col
      # [/CHANGE 2026-10]
      )
      colnames(x$benchmarks$difLL[[pop_es]])[
        ncol(x$benchmarks$difLL[[pop_es]])
      ] <- overlap_header
      # [CHANGE 2026-10 | audit] E5/O10: keep the reasons for NA overlaps as attribute
      # why an overlap is NA (if any), per row -- see overlap_column()
      attr(x$benchmarks$difLL[[pop_es]], "overlap_notes") <- attr(ov_col, "notes")
    }
    print_section(
      text_ld,
      function() {
        for (pop_es in names(x$benchmarks$difLL)) {
          print_rounded_es_value(x$benchmarks$difLL[[pop_es]], pop_es, model_type,
                                 green, reset,
                                 # [CHANGE 2026-10 | audit] B24: mark the reference population in the printed table
                                 is_reference = identical(pop_es, x$overlap$reference),
                                 reference_label = reference_label)
    # [/CHANGE 2026-10]
        }
      }, nchar(text_ld), text_color = blue, reset = reset
    )
  }

  return(invisible(x))
}
