print.con_goric <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  type <- x$type
  comparison <- x$comparison
  
  # # TO DO HIER
  # # Determine best hypo vs its complement (and thus conditional error prob.),
  # # but only if it is not Hunc or when complement&Heq not Hc.
  # best_Hypo_ICw <- which.max(x$result$goric.weights)
  # best_Hypo_ICw_name <- x$result$model[best_Hypo_ICw]
  # messagePrW <- FALSE
  # if ( (comparison != "complement" && best_Hypo_ICw_name != "unconstrained") || (x$Heq && best_Hypo_ICw_name != "complement") ) {
  #   Hypo_best <- unlist(x$hypotheses_usr[best_Hypo_ICw], use.names = F)
  #   if (any(x$priorICweights != x$priorICweights[1])) {
  #     priorICweights_best <- c(x$priorICweights[best_Hypo_ICw], 1-x$priorICweights[best_Hypo_ICw]) 
  #     messagePrW <- TRUE
  #     # TO DO make note mbt hoe priorICweights gemaakt nu! als niet vooraf gelijk dan
  #   } else {
  #     priorICweights_best <- c(0.5, 0.5)
  #   }
  #   inputList <- list(x$model.org, #x$objectList$H1,
  #                       type = type, sample_nobs = x$sample_nobs,
  #                     hypotheses = list(Hbest = Hypo_best),
  #                     comparison = "complement",
  #                     priorICweights = priorICweights_best,
  #                     Heq = FALSE,
  #                     penalty_factor = x$penalty_factor) 
  #   goric_HbestVsHc <- do.call(goric.default, inputList)
  #   #
  #   ICratio_HbestVsHc <- goric_HbestVsHc$ratio.gw[1,2]
  #   errorProb_Hbest <- goric_HbestVsHc$result$goric.weights[2]
  # } else if (comparison == "complement" && !x$Heq) {
  #   # If already H1 vs Hc, then we can also report on cond. error prob.
  #   ICratio_HbestVsHc <- NULL # Null, such that Future info does not appear in output
  #   errorProb_Hbest <- x$result$goric.weights[2]
  # } else {
  #   ICratio_HbestVsHc <- NULL
  #   errorProb_Hbest <- NULL
  # }
  # # TO DO HIER Verwerk bovenstaande in output (print en summary)
  # # print: dan had dit ms al eerder gedaan moeten worden?
  # # summary: hieronder poging, maar niet mooi nog.
  # # Gedachte (maar ws niet): Ms ook onderdeel van 'goric_object$result' laten zijn, 
  # #                          al wordt die tabel dan wel nog weer langer... dus ws niet doen, 
  # #                          dan ook minder duidelijk dat het posthoc is...
  # #          Wel opvraagbaar laten zijn!
  
  
  dig <- paste0("%6.", digits, "f")
  IC_lab <- toupper(type) # e.g., GORIC, GORICA, GORICC, GORICAC (used in the text)
  #x2 <- lapply(x$result[-1], sprintf, fmt = dig)
  x2 <- as.data.frame(lapply(x$result[, -1], function(column) {
    sapply(column, function(val) {
      if (is.na(val)) {
        return("")  # Zet NA om naar een lege string
      } else {
        return(sprintf(dig, val))  # Anders sprintf toepassen
      }
    })
  }))
  
  df <- data.frame(model = x$result$model, x2)
  objectnames <- as.character(df$model)

  # Confidence intervals (if computed by goric(); see goric_CI()):
  # extra columns in the Values and Weights tables, and the CIs for the ratios
  # of weights in the Conclusion text.
  has_CI <- !is.null(x$IC.CI) && !is.null(x$weights.CI) && !is.null(x$ratio.gw.lower)
  if (has_CI) {
    CI_pct <- paste0(format(100 * x$CI$level), "%")
    fmt_ci <- function(v) vapply(v, function(val) if (is.na(val)) "" else sprintf(dig, val), character(1))
  }
  # text " (95%CI: [lower, upper])" for the ratio of weights of hypothesis
  # 'm' versus hypothesis 'k' (names as in x$result$model); "" if not available
  ci_txt <- function(m, k) {
    if (!has_CI) return("")
    mi <- match(m, x$result$model); ki <- match(k, x$result$model)
    if (is.na(mi) || is.na(ki)) return("")
    lo <- x$ratio.gw.lower[mi, ki]; up <- x$ratio.gw.upper[mi, ki]
    if (is.na(lo) || is.na(up)) return("")
    paste0(" (", CI_pct, "CI: [", format_value(lo), ", ", format_value(up), "])")
  }
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  
  if (type == "goric") {
    cat("generalized order-restricted information criterion: \n\n")
  } else if (type == "gorica") {
    cat("generalized order-restricted information criterion approximation:\n\n")
  } else if (type == "goricc") {
    cat("small sample generalized order-restricted information criterion:\n\n")
  } else if (type == "goricac") {
    cat("small sample generalized order-restricted information criterion approximation:\n\n")
  }
  
  wt_bar_attributes <- lapply(x$objectList, function(obj) {
    list(
      method = attr(obj$wt.bar, "method"),
      converged = attr(obj$wt.bar, "converged"),
      total_draws = attr(obj$wt.bar, "total_bootstrap_draws"),
      errors = attr(obj$wt.bar, "error.idx")
    )
  })
  
  # Compute indicators
  wt_bar <- vapply(wt_bar_attributes, function(attr) attr$method == "boot", logical(1))
  #ceq_only <- vapply(x$objectList, function(obj) nrow(obj$constraints) == obj$neq, logical(1))
  ceq_only <- vapply(x$objectList, function(obj) nrow(obj$PT_Amat) == obj$PT_meq, logical(1))
  wt_bar <- wt_bar & !ceq_only
  
  if (any(wt_bar)) {
    wt_bar_attributes <- wt_bar_attributes[wt_bar]
    wt_method_boot <- x$objectList[wt_bar]
    #wt_attributes <- wt_bar_attributes
    max_nchar <- max(nchar(names(wt_method_boot)))
    
    # Summarize bootstrap information
    bootstrap_summary <- vapply(wt_bar_attributes, function(attr) {
      successful_draws <- attr$total_draws - length(attr$errors)
      paste0(successful_draws, ifelse(attr$converged, " (Converged)", " (Not converged)"))
    }, character(1))
    
    converged <- vapply(wt_bar_attributes, function(attr) attr$converged, logical(1))
    total_bootstrap_draws <- vapply(wt_bar_attributes, function(attr) attr$total_draws, integer(1))
    wt_bootstrap_errors <- sapply(wt_bar_attributes, function(attr) attr$errors)
    
    if (length(wt_method_boot) > 0) {
      successful_draws <- total_bootstrap_draws - sapply(wt_bootstrap_errors, length)
      
      hypo_messages <- names(x$objectList)
      
      if (length(hypo_messages) > 0) {
        messages_info <- identify_messages(x)
        
        rank_messages <- sapply(messages_info, function(x) x == "mix_weights_rank")
        NaN_messages  <- sapply(messages_info, function(x) x == "mix_weights_NaN")
      
        if (any(rank_messages)) {
          text_msg_1 <- paste("Note: Since the constraint matrix for hypotheses", paste0(sQuote(names(rank_messages)[rank_messages]), collapse = " and "), 
                              "is not full row-rank, we used the 'boot' method for calculating", 
                              "the penalty term value. For additional details, see '?goric' or the Vignette.")
        }
        if (any(NaN_messages)) {
          text_msg_2 <- paste("Note: Some returned mixing weights for hypotheses", paste0(sQuote(names(NaN_messages)[NaN_messages]), collapse = " and "), 
                              "are NaN (not a number). Continued the analysis with mix_weights = 'boot' method.")
        }
      }
      
      has_errors <- vapply(wt_bootstrap_errors, function(errors) length(errors) > 0, logical(1))
      not_all_converged <- !all(converged)
      not_all_draws_successful <- !all(successful_draws == total_bootstrap_draws)
      
      # for testing purposes only
      # has_errors <- TRUE
      # not_all_draws_successful <- TRUE
      
      if (any(has_errors) || not_all_converged) { 
        if (not_all_draws_successful || not_all_converged) {
          cat("Bootstrap-based penalty term calculation:\n")
          cat("  Number of bootstrap draws:", sapply(wt_bar_attributes, `[[`, "total_draws"), "\n")
          for (i in seq_along(bootstrap_summary)) {
            cat(sprintf("  Number of successful bootstrap draws for %*s: %s\n", 
                        max_nchar, names(wt_method_boot)[i], bootstrap_summary[i]))
          }
          cat("\n")
        } 
        text_msg_3 <- paste("Advise: If a substantial number of bootstrap draws fail to converge,", 
                            "the resulting penalty term may become unreliable. In such cases, it is advisable", 
                            "to increase the maximum number of bootstrap draws, e.g., control = list(mix_weights_bootstrap_limit = 1e5)")
      }
    }
  }
  
  cat("Results:\n\n")
  # Split the Results table into two blocks: the first with the raw values
  # (model, loglik, penalty, and the IC value itself), the second with the
  # corresponding weights (loglik.weights, penalty.weights, IC weight, and,
  # if present, the IC weight excluding unconstrained/Heq) -- rather than one
  # wide table mixing both kinds of quantities.
  values_cols  <- 1:4
  weights_cols <- c(1, 5:ncol(df))

  # Columns that appear in BOTH blocks -- model, and each raw value paired
  # with its corresponding weight (loglik/loglik.weights,
  # penalty/penalty.weights, the IC value/its weight) -- are each given one
  # shared printed width, so the two tables line up as one consistent grid
  # when stacked, even though e.g. "loglik" and "loglik.weights" are
  # different lengths as header text. The optional weights_without_unc/_heq
  # column, present only in the weights block with no counterpart to align
  # with, keeps its own natural width.
  shared_pairs <- list(c(1, 1), c(2, 5), c(3, 6), c(4, 7))
  widths_by_name <- list()
  for (pair in shared_pairs) {
    w <- max(nchar(names(df)[pair]),
             vapply(df[pair], function(col) max(nchar(col)), integer(1)))
    for (idx in pair) widths_by_name[[names(df)[idx]]] <- w
  }

  pad_block <- function(sub_df) {
    out <- format(sub_df, digits = digits, scientific = FALSE)
    for (cn in names(out)) {
      w <- widths_by_name[[cn]]
      if (!is.null(w)) {
        out[[cn]] <- sprintf("%*s", w, out[[cn]])
        names(out)[names(out) == cn] <- sprintf("%*s", w, cn)
      }
    }
    out
  }

  values_df  <- pad_block(df[, values_cols])
  weights_df <- pad_block(df[, weights_cols])
  if (has_CI) {
    # CI columns right after the IC values and right after the IC weights
    g_lo  <- fmt_ci(x$IC.CI[, "lower"]);      g_up  <- fmt_ci(x$IC.CI[, "upper"])
    gw_lo <- fmt_ci(x$weights.CI[, "lower"]); gw_up <- fmt_ci(x$weights.CI[, "upper"])
    # same printed width for g.lower & gw.lower and for g.upper & gw.upper, as
    # for the other corresponding columns of the two blocks
    pad_pair <- function(a, b, na, nb) {
      w <- max(nchar(c(a, b, na, nb)))
      list(a = sprintf("%*s", w, a), b = sprintf("%*s", w, b),
           na = sprintf("%*s", w, na), nb = sprintf("%*s", w, nb))
    }
    lo <- pad_pair(g_lo, gw_lo, "g.lower", "gw.lower")
    up <- pad_pair(g_up, gw_up, "g.upper", "gw.upper")
    values_df[[lo$na]] <- lo$a
    values_df[[up$na]] <- up$a
    w_ci <- data.frame(lo$b, up$b, stringsAsFactors = FALSE)
    names(w_ci) <- c(lo$nb, up$nb)
    weights_df <- cbind(weights_df[, 1:4, drop = FALSE], w_ci,
                        weights_df[, -(1:4), drop = FALSE])
  }
  cat("Values:\n")
  print(values_df, print.gap = 2, quote = FALSE)
  cat("\n")
  cat("Weights:\n")
  print(weights_df, print.gap = 2, quote = FALSE)
  if (has_CI) {
    cat(sprintf("\n.lower, .upper: bounds of the %s confidence interval for the %s values (g.xxx) and %s weights (gw.xxx); the CI methods used are stored in $CI_methods.\n",
                CI_pct, toupper(type), toupper(type)))
  }
  #cat("---\n")
  
  if (exists("text_msg_1")) message("---\n", text_msg_1)
  if (exists("text_msg_2")) message("---\n", text_msg_2)
  if (exists("text_msg_3")) message("---\n", text_msg_3)
  
  # Calculate the absolute difference between the logliks
  loglik <- x$result$loglik
  #-2 * (as.numeric(A)-as.numeric(B)))
  loglik_diff <- as.matrix(dist(loglik, diag = TRUE))
  # we only need the lower or upper part
  loglik_diff[upper.tri(loglik_diff)] <- NA
  diag(loglik_diff) <- NA
  # Create a matrix with the logical vector and assign dimension names
  loglik_diff_mat <- matrix(loglik_diff, nrow(loglik_diff), ncol(loglik_diff), 
                            dimnames = list(rownames(x$ratio.gw), colnames(x$ratio.gw)))
  # which hypotheses overlap, i.e., have equal log-likelihood values
  loglik_overlap <- which(loglik_diff_mat == 0, arr.ind = TRUE)
  # get unique combination from row- and columnnames
  overlap_unique_combinations <- apply(loglik_overlap, 1, function(x) {
    get_names <- c(rownames(loglik_diff_mat)[x[1]], gsub("vs. ", "", colnames(loglik_diff_mat)[x[2]]))
    paste(get_names, collapse = " vs. ")
  })
  
  overlap_sorted_vector <- sapply(overlap_unique_combinations, sort_combination)
  overlap_unique_combinations <- unique(overlap_sorted_vector)
  
  # remove all combinations involving unconstrained. They overlap by definition
  overlap_unique_combinations <- overlap_unique_combinations[!grepl("unconstrained", overlap_unique_combinations)]
  
  overlap_hypo <- gsub("vs\\.", "", overlap_unique_combinations)
  overlap_hypo <- strsplit(overlap_hypo, " ")
  overlap_hypo <- unique(unlist(overlap_hypo))
  overlap_hypo <- overlap_hypo[overlap_hypo != ""]
  
  best_hypo <- which.max(x$result[, 7])
  best_hypo_name <- x$result$model[best_hypo]
  best_hypo_overlap <- best_hypo_name %in% overlap_hypo
  
  # check if the log-likelihood of models are equal, this means that the hypotheses overlap (i.e., subset)
  overlap_note_printed <- length(overlap_hypo) > 0 & !x$Heq
  if (overlap_note_printed) {
    message("---\nNote: Hypotheses ", paste0(sQuote(overlap_hypo), collapse = " and "),
            " have equal log-likelihood values. This indicates that they overlap and",
            " that the ratio of their ", IC_lab, " weights reached its maximum.")
  }
  if (best_hypo_overlap) {
    message("Since the best hypothesis overlaps with one or more other hypotheses,",
" we recommend evaluating the overlap (or best hypothesis) vs. its complement.",
          " Run vignette(\"Guidelines_GORIC_output\") in your console for more",
          " information and an example.")
  }

  # Does the best hypothesis (by GORIC(A) weight, 'best_hypo') also have the
  # highest log-likelihood (loglik) value -- among the actual, substantive
  # hypotheses being compared? With unequal priorICweights, or with penalty
  # terms that differ enough between hypotheses, the hypothesis with the
  # highest GORIC(A) weight need not be the one with the highest raw loglik
  # -- which can be a sign that there isn't enough data to distinguish the
  # hypotheses well.
  #
  # "unconstrained" and "Heq" are deliberately excluded from the comparison
  # here (though 'best_hypo' itself may still BE one of them): both are
  # reference/safeguard categories, not substantive alternative hypotheses,
  # and each can trivially attain the global maximum loglik by construction
  # (unconstrained always does; Heq is nested inside every order-restricted
  # hypothesis, so it can never exceed them anyway) -- comparing against
  # those would make this fire on essentially every analysis with an active
  # constraint, which isn't the (rarer, more informative) condition intended
  # here. A tie for the highest loglik among the remaining hypotheses
  # (best_hypo's own value equal to that maximum, up to floating-point
  # tolerance) does NOT trigger this -- "the highest" includes a tie for
  # highest -- and nor does best_hypo's own loglik simply exceeding it (e.g.
  # when best_hypo IS "unconstrained" or "Heq" itself).
  reference_categories <- c("unconstrained", "Heq")
  genuine_idx <- !(x$result$model %in% reference_categories)
  best_hypo_has_max_loglik <- if (!any(genuine_idx)) {
    TRUE
  } else {
    max_loglik_genuine <- max(x$result$loglik[genuine_idx], na.rm = TRUE)
    isTRUE(all.equal(x$result$loglik[best_hypo], max_loglik_genuine)) ||
      x$result$loglik[best_hypo] > max_loglik_genuine
  }

  # Is the loglik(.weight) of the best hypothesis "close" to that of another
  # hypothesis -- other than an exact tie (already reported above via
  # overlap_hypo) and other than the reference/safeguard categories
  # "unconstrained"/"Heq" (which, like in the mismatch check above, are
  # excluded rather than flagged the same way -- being close to a trivial
  # reference category isn't the same kind of finding as being close to
  # another substantive hypothesis)? Flagged via the RATIO of loglik.weights
  # ('ratio.lw', the same quantity reported elsewhere) rather than the raw
  # loglik difference, so the cutoff doesn't depend on the scale of the
  # data. "Close" here means that ratio falls within [0.9, 1.1] (i.e. within
  # +/-10%) -- a simple, symmetric-in-practice heuristic rather than a
  # formal likelihood-ratio/chi-bar-square test: the hypotheses compared
  # here are order restrictions, not equality constraints, so they are
  # typically not nested in the way a classical LRT assumes -- a
  # chi-bar-square-based cutoff would need its own per-pair mixing-weight
  # computation (as con_pbar()/wt.bar already do for a single hypothesis
  # vs. its complement) rather than one general-purpose threshold,
  # considerably more machinery for a print-time heads-up message.
  close_idx <- integer(0)
  if (!is.null(x$ratio.lw) && nrow(x$result) > 1) {
    ratio_lw_best <- x$ratio.lw[best_hypo, ]
    is_self <- names(ratio_lw_best) == paste0("vs. ", best_hypo_name)
    is_reference <- x$result$model %in% reference_categories
    is_exact_tie <- x$result$model %in% overlap_hypo
    close_idx <- which(!is_self & !is_reference & !is_exact_tie &
                       ratio_lw_best >= 0.9 & ratio_lw_best <= 1.1)
  }
  close_names <- if (length(close_idx) > 0) sQuote(x$result$model[close_idx]) else character(0)
  close_hypo_text <- if (length(close_names) == 1) {
    paste0("Hypothesis ", close_names)
  } else if (length(close_names) > 1) {
    paste0("Hypotheses ", paste(close_names[-length(close_names)], collapse = ", "),
          ", and ", close_names[length(close_names)])
  } else {
    NULL
  }

  # These two findings can co-occur (the best hypothesis need not have the
  # top loglik among the genuine hypotheses AND still be close in loglik to
  # one or more of them) -- printed as ONE combined message when they do,
  # rather than two separate ones repeating much of the same context.
  if (!best_hypo_has_max_loglik && !is.null(close_hypo_text)) {
    message("---\nNote: Hypothesis ", sQuote(best_hypo_name),
            " is the best hypothesis but note that it does not have the highest",
            " log-likelihood (loglik) value. This may indicate that the sample size",
            " might be too low to base a conclusion on. Additionally, its",
            " log-likelihood (loglik) value is close to the log-likelihood (loglik)",
            " value of ", close_hypo_text, ". You may want to use the benchmark()",
            " function to investigate whether this is the case. If so, then this",
            " indicates that there is support for their boundary. More details can",
            " be found in 'Decision tree - Interpretation.pdf' and/or 'Guidelines",
            " output GORIC.html', available from",
            " https://github.com/rebeccakuiper/Tutorials.")
  } else if (!best_hypo_has_max_loglik) {
    message("---\nNote: Hypothesis ", sQuote(best_hypo_name),
            " is the best hypothesis but note that it does not have the highest",
            " log-likelihood (loglik) value. This may indicate that the sample size",
            " might be too low to base a conclusion on. More details can be found in",
            " 'Decision tree - Interpretation.pdf' and/or 'Guidelines output GORIC.html',",
            " available from https://github.com/rebeccakuiper/Tutorials.")
  } else if (!is.null(close_hypo_text)) {
    message("---\nNote: The log-likelihood (loglik) value of the best hypothesis",
            " Hypothesis ", sQuote(best_hypo_name), " is close to the",
            " log-likelihood (loglik) value of ", close_hypo_text, ". You may want to",
            " use the benchmark() function to investigate whether this is the case.",
            " If so, then this indicates that there is support for their boundary.",
            " More details can be found in 'Decision tree - Interpretation.pdf'",
            " and/or 'Guidelines output GORIC.html', available from",
            " https://github.com/rebeccakuiper/Tutorials.")
  } else if (!overlap_note_printed && length(df$model) > 1) {
    # None of the loglik-related notes above fired (no exact tie, no
    # mismatch, no close ratio) -- so, unlike those cases, a general pointer
    # to check the loglik(.weights) values isn't redundant here and is shown
    # instead.
    message("---\nTip: It is advisable to compare the log-likelihood values",
            " (loglik) and/or weights (loglik.weights) of the best hypothesis",
            " with that of the other specified/informative hypotheses. More",
            " details can be found in 'Decision tree - Interpretation.pdf'",
            " and/or 'Guidelines output GORIC.html', available from",
            " https://github.com/rebeccakuiper/Tutorials.")
  }

  if (comparison == "unconstrained" && length(df$model) == 2 &&
      # If comparison = "complement" and ceq only, than comparison is set to unconstrained internally.
      nrow(x$objectList[[1]]$constraints) != x$objectList[[1]]$neq) { 
    message("---\nAdvise: Are you certain you wish to assess the order-restricted hypothesis", 
            " in comparison to the unconstrained one, rather than its complement?", 
            " Note that the order-restricted hypothesis overlaps (is contained in) the unconstrained,",
            " while it does not overlap with its complement (except for their boundary).")
  }

  
  if (length(df$model) > 1) {
    cat("\nConclusion:\n")
  }
  
  if (comparison == "complement" && length(overlap_unique_combinations) == 0 && !x$Heq) {
    formatted_numbers <- sprintf("%.3f %.3f", x$result[[7]][1], x$result[[7]][2])
    numbers <- strsplit(formatted_numbers, " ")[[1]]
    if (as.numeric(numbers[1]) / as.numeric(numbers[2]) > 1) {
      if (x$ratio.gw[1, 2] > 10000) {
        support_ratio <- sprintf("%.2e", x$ratio.gw[1, 2])
      } else {
        support_ratio <- sprintf("%.2f", x$ratio.gw[1, 2])
      }
      cat(paste0("The order-restricted hypothesis ", sQuote(objectnames[1]), 
                 " has ", support_ratio, " times more support than its complement",
                 ci_txt(objectnames[1], objectnames[2]), ".\n"))
    } else if (as.numeric(numbers[1]) / as.numeric(numbers[2]) < 1) {
      result <- paste(numbers[1], "/", numbers[2], "< 1", sep = " ")
      gw_c <- as.numeric(numbers[2]) / as.numeric(numbers[1])
      if (gw_c > 10000) {
        support_ratio_c <- sprintf("%.2e", gw_c)
      } else {
        support_ratio_c <- sprintf("%.2f", gw_c)
      }
      cat("The order-restricted hypothesis", sQuote(objectnames[1]), 
          "has less support (namely, ", result, 
          "times more support) than its complement.",
      "That is, the complement has", support_ratio_c, "times more support than", 
      paste0(sQuote(objectnames[1]), ci_txt(objectnames[2], objectnames[1]), ".\n"))      
    } else {
      result <- paste(numbers[1], "/", numbers[2], "= 1", sep = " ")
      cat("The order-restricted hypothesis", sQuote(objectnames[1]), "and the complement have equal support:", result)      
    }
    # cat(paste0("---\nThe order-restricted hypothesis ", objectname1, 
    #            " has ", support_ratio, " times more support than its complement.\n\n"))
    #
    # # TO DO HIER maak mooier en ms nog onderdeel van een opvraagbare lijst ofzo
    # if (!is.null(errorProb_Hbest)) {
    #   message(paste0("\nThe conditional error probability, that is, the probability ",
    #   "(given the data) that this conclusion is incorrect, is ", round(errorProb_Hbest, 2), "."))
    # }
  } else if (comparison == "complement" && x$Heq) { 
    modelnames <- x$result$model[!x$result$model == "Heq"]
    if (best_hypo_name != "Heq") {
      goric_weights_without_heq <- x$result[, 8][!is.na(x$result[, 8])]
      goric_rw_without_heq <- goric_weights_without_heq %*% t(1/goric_weights_without_heq)
      diag(goric_rw_without_heq) <- 1L
      colnames(goric_rw_without_heq) <- paste0("vs. ", modelnames)
      goric_rw_without_heq_best_hypo <- goric_rw_without_heq[best_hypo-1, ]
      goric_rw_without_heq_best_hypo <- goric_rw_without_heq_best_hypo[goric_rw_without_heq_best_hypo != 1]
      goric_rw_without_heq_best_hypo <- sapply(goric_rw_without_heq_best_hypo, format_value)
      best_hypos_rest <- paste(df$model[!df$model %in% c(best_hypo_name, "Heq")])
      
      if (best_hypo_name == "complement") {
        message <- paste0("- The complement is the best in the set, as it has the highest ", IC_lab, " weight.",
                          " Since the complement has a higher ", IC_lab, " weight than the equality-restricted", 
                          " hypothesis Heq, we can now inspect the relative support for the complement",
                          " against the order-restricted hypothesis ", modelnames[modelnames != "complement"], ":")
      } else if (best_hypo_name != "complement") {
        message <- paste0("- The order-restricted hypothesis ", modelnames[modelnames != "complement"], " is the best",
                          " in the set, as it has the highest ", IC_lab, " weight. Since it has a higher ", IC_lab, " weight",
                          " than the equality-restricted hypothesis (Heq), we can now inspect the relative support for", 
                          " the order-restricted hypothesis against its complement:")
      } 

      for (i in seq_along(best_hypos_rest)) {
        message <- paste0(message, "\n  * ", sQuote(best_hypo_name), " is ", 
                          goric_rw_without_heq_best_hypo[i], 
                          " times more supported than ", sQuote(best_hypos_rest[i]),
                          ci_txt(best_hypo_name, best_hypos_rest[i]))

        # Voeg punt toe, behalve als het de laatste hypothese is
        if (i < length(best_hypos_rest)) {
          message <- paste0(message, ".")
        } else {
          message <- paste0(message, ".")
        }
      }
      cat(paste0(message, "\n"))
    } else {
      message <- paste0("\n- The equality-constrained hypothesis (Heq) is the best in the set,", 
                        " as it has the highest ", IC_lab, " weight.")
      
      message <- paste0(message, "\n- Since the order-restricted hypotheses contain", 
                        " the equality-constrained hypothesis (Heq), inspecting", 
                         " the relative support is not meaningful.")
      
      cat(paste0(message, "\n"))
    }
  } else if (comparison == "none" && length(overlap_unique_combinations) == 0 && length(df$model) == 2) {
    class(x$ratio.gw) <- "numeric"
    if (x$ratio.gw[1, 2] > 10000) {
      support_ratio <- sprintf("%.2e", x$ratio.gw[1, 2])
    } else {
      support_ratio <- sprintf("%.2f", x$ratio.gw[1, 2])
    }
    objectname1 <- sQuote(objectnames[1])
    objectname2 <- sQuote(objectnames[2])
    cat(paste0("The order-restricted hypothesis ", objectname1, 
               " has ", support_ratio, " times more support than ", objectname2,
               ci_txt(objectnames[1], objectnames[2]), ".\n\n"))
  } else if (comparison == "unconstrained" && length(overlap_unique_combinations) == 0 && length(df$model) == 2) { 
    formatted_numbers <- sprintf("%.3f %.3f", x$result[[7]][1], x$result[[7]][2])
    numbers <- strsplit(formatted_numbers, " ")[[1]]
    if (as.numeric(numbers[1]) / as.numeric(numbers[2]) > 1) {
      result <- paste(numbers[1], "/", numbers[2], "> 1", sep = " ")
      cat("The order-restricted hypothesis", sQuote(objectnames[1]), "has", result, 
          paste0("times more support than the unconstrained", ci_txt(objectnames[1], objectnames[2]), ".\n\n"))
    } else if (as.numeric(numbers[1]) / as.numeric(numbers[2]) < 1) {
      result <- paste(numbers[1], "/", numbers[2], "< 1", sep = " ")
      cat("The order-restricted hypothesis", sQuote(objectnames[1]), 
          "has less support (namely, ", result, 
          paste0("times more support) than the unconstrained", ci_txt(objectnames[1], objectnames[2]), ".\n\n"))
    } else {
      result <- paste(numbers[1], "/", numbers[2], "= 1", sep = " ")
      cat("The order-restricted hypothesis", sQuote(objectnames[1]), "and the unconstrained have equal support:", result, "\n\n")      
    }
  } else if ( (comparison == "unconstrained" && length(df$model) > 2) )  {
    #best_hypo <- which.max(x$result[, 7])
    #best_hypo_name <- x$result$model[best_hypo]
    modelnames <- x$result$model[!x$result$model == "unconstrained"]
    if (best_hypo_name != "unconstrained") {
      goric_weights_without_unc <- x$result[, 8][!is.na(x$result[, 8])]
      goric_rw_without_unc <- goric_weights_without_unc %*% t(1/goric_weights_without_unc)
      diag(goric_rw_without_unc) <- 1L
      colnames(goric_rw_without_unc) <- paste0("vs. ", modelnames)
      goric_rw_without_unc_best_hypo <- goric_rw_without_unc[best_hypo, ]
      goric_rw_without_unc_best_hypo <- goric_rw_without_unc_best_hypo[goric_rw_without_unc_best_hypo != 1]
      goric_rw_without_unc_best_hypo <- sapply(goric_rw_without_unc_best_hypo, format_value)
      best_hypos_rest <- paste(df$model[!df$model %in% c(best_hypo_name, "unconstrained")])
      # Step 1: Check if the best hypothesis in the set is not weak
      message <- paste0("- The order-restricted hypothesis ", sQuote(best_hypo_name), 
                        " is the best in the set, as it has the highest ", IC_lab, " weight.")
      
      # Step 2: if not weak, compare it against all other hypotheses in the set
      message <- paste0(message, "\n- Since ", sQuote(best_hypo_name), " has a higher", 
      " ", IC_lab, " weight than the unconstrained hypothesis, it is not considered weak.", 
      " We can now inspect the relative support for ", sQuote(best_hypo_name), " against",
      " the other order-restricted hypotheses:")
      
      for (i in seq_along(best_hypos_rest)) {
        if (best_hypo_overlap & best_hypos_rest[i] %in% overlap_hypo) {
          # Als er overlap is, voeg toe dat de relatieve support zijn maximum heeft bereikt
          message <- paste0(message, "\n  * ", sQuote(best_hypo_name), " is ", 
                            goric_rw_without_unc_best_hypo[i], 
                            " times more supported than ", sQuote(best_hypos_rest[i]), 
                            ci_txt(best_hypo_name, best_hypos_rest[i]),
                            " (This relative support reached its maximum, see Note)")
        } else {
          # Als er geen overlap is, geef normale relatieve support
          message <- paste0(message, "\n  * ", sQuote(best_hypo_name), " is ", 
                            goric_rw_without_unc_best_hypo[i], 
                            " times more supported than ", sQuote(best_hypos_rest[i]),
                            ci_txt(best_hypo_name, best_hypos_rest[i]))
        }
        
        # Voeg punt toe, behalve als het de laatste hypothese is
        if (i < length(best_hypos_rest)) {
          message <- paste0(message, ".")
        } else {
          message <- paste0(message, ".")
        }
      }
      cat(paste0(message, "\n"))
    } else {
      message <- paste0("\n- The unconstrained hypothesis is the best in the set,", 
                        " as it has the highest ", IC_lab, " weight. As a result, the order-restricted", 
                        " hypotheses are considered weak.")

      message <- paste0(message, "\n- Since all the order-restricted hypotheses are weak,",
                        " inspecting their relative support is not meaningful.")
            
      #cat(paste0("---\n", message, "\n"))
      cat(paste0(message, "\n"))
    }
  } else if (comparison == "none" && length(df$model) > 1) {
    #best_hypo <- which.max(x$result[, 7])
    #best_hypo_name <- x$result$model[best_hypo]
    modelnames <- x$result$model
    
    goric_rw <- x$ratio.gw
    goric_rw_best_hypo <- goric_rw[best_hypo, ]
    #goric_rw_best_hypo <- goric_rw_best_hypo[goric_rw_best_hypo != 1]
    # Should not filter on the value 1, but on the name 'best_hypo_name':
    filterout_best <- names(goric_rw[best_hypo, ]) != paste0("vs. ", best_hypo_name)
    goric_rw_best_hypo <- goric_rw_best_hypo[filterout_best]
    # 
    goric_rw_best_hypo <- sapply(goric_rw_best_hypo, format_value)
    best_hypos_rest <- paste(df$model[!df$model %in% best_hypo_name])
    
    message <- ""
    for (i in seq_along(best_hypos_rest)) {
      if (best_hypo_overlap & best_hypos_rest[i] %in% overlap_hypo) {
        # Als er overlap is, voeg toe dat de relatieve support zijn maximum heeft bereikt
        message <- paste0(message, "  * ", sQuote(best_hypo_name), " is ", 
                          goric_rw_best_hypo[i], 
                          " times more supported than ", sQuote(best_hypos_rest[i]), 
                          ci_txt(best_hypo_name, best_hypos_rest[i]),
                          " (This relative support reached its maximum, see Note)")
      } else {
        # Als er geen overlap is, geef normale relatieve support
        message <- paste0(message, "  * ", sQuote(best_hypo_name), " is ", 
                          goric_rw_best_hypo[i], 
                          " times more supported than ", sQuote(best_hypos_rest[i]),
                          ci_txt(best_hypo_name, best_hypos_rest[i]))
      }
      # Voeg punt toe, behalve als het de laatste hypothese is
      if (i < length(best_hypos_rest)) {
        message <- paste0(message, ".\n")
      } else {
        message <- paste0(message, ".")
      }
    }
    cat(paste0(message, "\n"))
  } else if (length(overlap_unique_combinations) == 0 && length(df$model) > 2) {
    if (!is.null(x$ratio.gw)) {
      if (type == "goric") {
        cat("---\n\nRatio GORIC-weights:\n")
      } else if (type == "gorica") {
        cat("---\n\nRatio GORICA-weights:\n")
      } else if (type == "goricc") {
        cat("---\n\nRatio GORICC-weights:\n")
      } else if (type == "goricac") {
        cat("---\n\nRatio GORICAC-weights:\n") 
      }
      
      ratio.gw <- apply(x$ratio.gw, 2, sprintf, fmt = dig)
      rownames(ratio.gw) <- rownames(x$ratio.gw)
      class(ratio.gw) <- "numeric"
      
      if (max(ratio.gw, na.rm = TRUE) >= 1e4) {
        print(format(ratio.gw, digits = digits, scientific = TRUE, trim = TRUE), 
              print.gap = 2, quote = FALSE, right = TRUE) 
      } else {
        print(format(ratio.gw, digits = digits, scientific = FALSE, trim = TRUE), 
              print.gap = 2, quote = FALSE, right = TRUE)
      }
    }
  }

  # Tip pointing at the benchmarks function, for labeling/interpreting the
  # quantified output above (the ratio(s) of GORIC(A)/loglik/penalty
  # weights): a ratio on its own (e.g. 5.2) doesn't say how strong that level
  # of support actually is in this specific setup -- benchmark_means()/
  # benchmark_asymp() simulate the ratio's own distribution under, among
  # others, a "no effect" population, which is what such labeling requires.
  # Shown once, whenever there is more than one hypothesis to compare (i.e.
  # whenever the Conclusion section above -- which is where a ratio of
  # weights gets reported -- was printed at all).
  if (length(df$model) > 1) {
    message("---\nTip: The benchmark() function can be used to label the height",
            " of the output, for example, the height of the ratio(s) of weights",
            " (e.g., to say whether that value denotes medium, high, or very",
            " large evidence). More details can be found in",
            " 'Guidelines_GORIC-benchmarks.html', available from",
            " https://github.com/rebeccakuiper/Tutorials.")
  }

  # # TO DO HIER maak mooier en ms nog onderdeel van een opvraagbare lijst ofzo
  # if (!is.null(ICratio_HbestVsHc) && !is.null(errorProb_Hbest)) {
  #   message(paste0("\nrestriktor Future research information based on a post-hoc analysis: \n",
  #   "The best hypothesis ", best_hypo_name, " is ", round(ICratio_HbestVsHc, 2), 
  #   "times more likely than its complement. \n",
  #   "The conditional error probability, that is, the probability ",
  #   "(given the data) that this conclusion is incorrect, is ", round(errorProb_Hbest, 2), "."))
  #   #
  #   if (messagePrW) {
  #     message(paste0("\nThis result is based on rescaling the priorICweights to: ", priorICweights_best))
  #   }
  # }
  
  # Above the best hypothesis is determined based on ICweights.
  # In the case of unequal priorICweights, this can differ from conclusion based on ICvalues.
  # Check and, if, give message:
  best_hypo_IC <- which.min(x$result[, 4])
  best_hypo_IC_name <- x$result$model[best_hypo_IC]
  if (best_hypo_IC != best_hypo) {
    message(paste0("\nrestriktor Note: The best hypothesis ", best_hypo_name, " is based on the highest IC weight. \n",
                   "Based on the (smallest) IC values, the best hypothesis is ", best_hypo_IC_name, ". \n",
                   "This difference in conclusion is due to the specified priorICweights."))
  }
  
  if (x$penalty_factor != 2) {
    message(sprintf("\nrestriktor Message: Note that a penalty factor of %s (default is 2) is used in the calculation of the %s value, that is -2 x log-likelihood + penalty_factor x penalty",
                    x$penalty_factor, x$type))
  }
}

# TO DO print ook de warnings hier nog eens, anders mis je die ms
#       Kijk of messages en/of warnings en zag dan onderaan dat die er zijn en geef dan de code om die te zien!
# Ws zeg dat er messgaes zijn en zeg hoe die in te zien.
# vang ze dus op en geef ze niet direct (niet tijdens runnen en niet bij print).
#
# TO DO veel output (bvb warnings, conclusion) is doorlopende tekst zonder 'line breaks' en dus in de pdf zie je maar de helft; het is beter om line breaks toe te voegen, zodat de output nooit meer is dan 76 karakters.