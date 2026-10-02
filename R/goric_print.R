print.con_goric <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  type <- x$type
  comparison <- x$comparison
  
  dig <- paste0("%6.", digits, "f")
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
  
  cat("Results:\n")
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
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
  # [CHANGE 2026-10 | audit] B1: exact match on the model name 'unconstrained' (no substring grep); as.character for empty overlap
  # (exact match on the model name; a hypothesis named, e.g., 'unconstrained1' is kept)
  # (as.character: an empty list when there is no overlap)
  overlap_unique_combinations <- as.character(overlap_unique_combinations)
  overlap_unique_combinations <- overlap_unique_combinations[!vapply(
    strsplit(overlap_unique_combinations, " vs. ", fixed = TRUE), 
    function(s) "unconstrained" %in% s, logical(1))]
    # [/CHANGE 2026-10]
  
  overlap_hypo <- gsub("vs\\.", "", overlap_unique_combinations)
  overlap_hypo <- strsplit(overlap_hypo, " ")
  overlap_hypo <- unique(unlist(overlap_hypo))
  overlap_hypo <- overlap_hypo[overlap_hypo != ""]
  
  best_hypo <- which.max(x$result[, 7])
  # [CHANGE 2026-10 | audit] B31: no conclusion when all IC weights are NaN
  if (length(best_hypo) == 0L) {
    # all IC weights are NaN (e.g., non-finite penalties): no conclusion
    message("---\nNote: The IC weights are not available (NaN), e.g., because the ",
            "penalty term is not finite. No conclusion can be drawn.")
    return(invisible(x))
  }
  # [/CHANGE 2026-10]
  best_hypo_name <- x$result$model[best_hypo]
  best_hypo_overlap <- best_hypo_name %in% overlap_hypo
  
  # check if the log-likelihood of models are equal, this means that the hypotheses overlap (i.e., subset)
  if (length(overlap_hypo) > 0 & !x$Heq) {
    message("---\nNote: Hypotheses ", paste0(sQuote(overlap_hypo), collapse = " and "), 
            " have equal log-likelihood values. This indicates that they overlap and",
            " that the ratio of their GORIC(A) weights reached its maximum.")
  }
  if (best_hypo_overlap) {
    message("Since the best hypothesis overlaps with one or more other hypotheses,", 
" we recommend evaluating the overlap (or best hypothesis) vs. its complement.", 
          " Run vignette(\"Guidelines_GORIC_output\") in your console for more", 
          " information and an example.")
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
    # [CHANGE 2026-10 | audit] A2/B10: conclusion via support_sentence(), ratio taken by name from x$ratio.gw
    cat(support_sentence(objectnames[1], "complement", 
                         x$ratio.gw[objectnames[1], "vs. complement"],
                         a_label = paste("The order-restricted hypothesis", 
                                         sQuote(objectnames[1])),
                         b_label = "its complement"), ".", sep = "")
                         # [/CHANGE 2026-10]
  } else if (comparison == "complement" && x$Heq) { 
    modelnames <- x$result$model[!x$result$model == "Heq"]
    if (best_hypo_name != "Heq") {
      # [CHANGE 2026-10 | audit] A2: ratios of the best hypothesis by name (only the best itself excluded, ties kept)
      # ratios of the best hypothesis vs. the other models, taken by name 
      # from the ratio matrix (the ratios do not depend on the renormalisation
      # without Heq); only the best hypothesis itself is left out
      best_hypos_rest <- paste(df$model[!df$model %in% c(best_hypo_name, "Heq")])
      ratios_best <- x$ratio.gw[best_hypo_name, paste0("vs. ", best_hypos_rest)]
      # [/CHANGE 2026-10]
      
      if (best_hypo_name == "complement") {
        message <- paste0("- The complement is the best in the set, as it has the highest GORIC(A) weight.",
                          " Since the complement has a higher GORIC(A) weight than the equality-restricted", 
                          " hypothesis Heq, we can now inspect the relative support for the complement",
                          " against the order-restricted hypothesis ", modelnames[modelnames != "complement"], ":")
      } else if (best_hypo_name != "complement") {
        message <- paste0("- The order-restricted hypothesis ", modelnames[modelnames != "complement"], " is the best",
                          " in the set, as it has the highest GORIC(A) weight. Since it has a higher GORIC(A) weight",
                          " than the equality-restricted hypothesis (Heq), we can now inspect the relative support for", 
                          " the order-restricted hypothesis against its complement:")
      } 

      for (i in seq_along(best_hypos_rest)) {
        # [CHANGE 2026-10 | audit] A2/B10: sentence via support_sentence()
        message <- paste0(message, "\n  * ", 
                          support_sentence(best_hypo_name, best_hypos_rest[i], 
                                           ratios_best[i]), ".")
                                           # [/CHANGE 2026-10]
      }
      cat(paste0(message, "\n"))
    } else {
      message <- paste0("\n- The equality-constrained hypothesis (Heq) is the best in the set,", 
                        " as it has the highest GORIC(A) weight.")
      
      message <- paste0(message, "\n- Since the order-restricted hypotheses contain", 
                        " the equality-constrained hypothesis (Heq), inspecting", 
                         " the relative support is not meaningful.")
      
      cat(paste0(message, "\n"))
    }
  } else if (comparison == "none" && length(overlap_unique_combinations) == 0 && length(df$model) == 2) {
    # [CHANGE 2026-10 | audit] A2/B10: conclusion via support_sentence(), ratio taken by name from x$ratio.gw
    cat(paste0(support_sentence(objectnames[1], objectnames[2], 
                                x$ratio.gw[objectnames[1], paste0("vs. ", objectnames[2])],
                                a_label = paste("The order-restricted hypothesis", 
                                                sQuote(objectnames[1]))),
               ".\n\n"))
               # [/CHANGE 2026-10]
  } else if (comparison == "unconstrained" && length(overlap_unique_combinations) == 0 && length(df$model) == 2) { 
    # [CHANGE 2026-10 | audit] A2/B10: conclusion via support_sentence(), ratio taken by name from x$ratio.gw
    cat(paste0(support_sentence(objectnames[1], "unconstrained", 
                                x$ratio.gw[objectnames[1], "vs. unconstrained"],
                                a_label = paste("The order-restricted hypothesis", 
                                                sQuote(objectnames[1])),
                                b_label = "the unconstrained"),
               ".\n\n"))
               # [/CHANGE 2026-10]
  } else if ( (comparison == "unconstrained" && length(df$model) > 2) )  {
    #best_hypo <- which.max(x$result[, 7])
    #best_hypo_name <- x$result$model[best_hypo]
    modelnames <- x$result$model[!x$result$model == "unconstrained"]
    if (best_hypo_name != "unconstrained") {
      # [CHANGE 2026-10 | audit] A2: ratios of the best hypothesis by name (only the best itself excluded, ties kept)
      # ratios of the best hypothesis vs. the other hypotheses, taken by name
      # from the ratio matrix (the ratios do not depend on the renormalisation
      # without the unconstrained model); only the best hypothesis itself is
      # left out, so ties (ratio 1) are reported as such
      best_hypos_rest <- paste(df$model[!df$model %in% c(best_hypo_name, "unconstrained")])
      ratios_best <- x$ratio.gw[best_hypo_name, paste0("vs. ", best_hypos_rest)]
      # [/CHANGE 2026-10]
      # Step 1: Check if the best hypothesis in the set is not weak
      message <- paste0("- The order-restricted hypothesis ", sQuote(best_hypo_name), 
                        " is the best in the set, as it has the highest GORIC(A) weight.")
      
      # Step 2: if not weak, compare it against all other hypotheses in the set
      message <- paste0(message, "\n- Since ", sQuote(best_hypo_name), " has a higher", 
      " GORIC(A) weight than the unconstrained hypothesis, it is not considered weak.", 
      " We can now inspect the relative support for ", sQuote(best_hypo_name), " against",
      " the other order-restricted hypotheses:")
      
      for (i in seq_along(best_hypos_rest)) {
        # [CHANGE 2026-10 | audit] A2/B10: sentence via support_sentence() (max_note for overlap)
        # in case of overlap, add that the relative support reached its maximum
        max_note <- if (best_hypo_overlap & best_hypos_rest[i] %in% overlap_hypo) {
          " (This relative support reached its maximum, see Note)"
        } else { "" }
        message <- paste0(message, "\n  * ", 
                          support_sentence(best_hypo_name, best_hypos_rest[i], 
                                           ratios_best[i], max_note = max_note), ".")
                                           # [/CHANGE 2026-10]
      }
      cat(paste0(message, "\n"))
    } else {
      message <- paste0("\n- The unconstrained hypothesis is the best in the set,", 
                        " as it has the highest GORIC(A) weight. As a result, the order-restricted", 
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
    
    # [CHANGE 2026-10 | audit] A2: ratios of the best hypothesis by name
    # ratios of the best hypothesis vs. the other hypotheses, by name
    best_hypos_rest <- paste(df$model[!df$model %in% best_hypo_name])
    ratios_best <- x$ratio.gw[best_hypo_name, paste0("vs. ", best_hypos_rest)]
    # [/CHANGE 2026-10]
    
    message <- ""
    for (i in seq_along(best_hypos_rest)) {
      # [CHANGE 2026-10 | audit] A2/B10: sentence via support_sentence() (max_note for overlap)
      # in case of overlap, add that the relative support reached its maximum
      max_note <- if (best_hypo_overlap & best_hypos_rest[i] %in% overlap_hypo) {
        " (This relative support reached its maximum, see Note)"
      } else { "" }
      message <- paste0(message, "  * ", 
                        support_sentence(best_hypo_name, best_hypos_rest[i], 
                                         ratios_best[i], max_note = max_note), 
                        if (i < length(best_hypos_rest)) ".\n" else ".")
                        # [/CHANGE 2026-10]
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
      
      # [CHANGE 2026-10 | audit] fmt_ratio_matrix(); R5: '(best)' only in the printed rownames
      ratio.gw <- fmt_ratio_matrix(x$ratio.gw, dig)
      rownames(ratio.gw) <- mark_best_hypo(rownames(x$ratio.gw), x)
      
      if (max(ratio.gw, na.rm = TRUE) >= 1e4) {
        print(format(ratio.gw, digits = digits, scientific = TRUE, trim = TRUE), 
              print.gap = 2, quote = FALSE, right = TRUE) 
      } else {
        print(format(ratio.gw, digits = digits, scientific = FALSE, trim = TRUE), 
              print.gap = 2, quote = FALSE, right = TRUE)
      }
    }
  } 
  
  # [CHANGE 2026-10 | Rebecca] note when the best hypothesis by IC weight differs from the one by IC value (priorICweights)
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
  # [/CHANGE 2026-10]
  
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

# [CHANGE 2026-10 | audit] new function format_support_ratio(): formatting of ratios in the conclusion text
# Format a ratio of IC weights for the conclusion text: 3 decimals, scientific
# notation for very large (>= 1e4) or very small (< 1e-2) ratios.
format_support_ratio <- function(r) {
  if (!is.finite(r)) {
    return(as.character(r))
  }
  if (r >= 1e4 || (r > 0 && r < 1e-2)) {
    sprintf("%.2e", r)
  } else {
    sprintf("%.3f", r)
  }
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] A2/B10: new function support_sentence(): relative-support sentence (ties, r < 1, 0 and Inf)
# Sentence describing the relative support of model 'a' versus model 'b',
# given the ratio r = weight(a) / weight(b) (taken from x$ratio.gw by name).
# Handles ties (r == 1), r < 1 (phrased from the better supported model),
# a zero weight (r = Inf or 0) and two zero weights (r = NaN).
# 'a_label' is used for 'a' at the start of the sentence (e.g., "The 
# order-restricted hypothesis 'H1'"); 'b_label' is used for 'b' (e.g., 
# "its complement"). When the sentence starts with 'b' (r == 0), its label
# is capitalised.
support_sentence <- function(a, b, r, a_label = sQuote(a), b_label = sQuote(b), 
                             max_note = "") {
  if (is.nan(r) || is.na(r)) {
    return(paste0(a_label, " and ", b_label, " both have an IC weight of 0; ",
                  "their relative support is undefined"))
  }
  if (r == 1) {
    return(paste0(a_label, " and ", b_label, " have equal support", max_note))
  }
  if (r > 1) {
    if (is.infinite(r)) {
      return(paste0(a_label, " is infinitely more supported than ", b_label,
                    " (which has an IC weight of 0)"))
    }
    return(paste0(a_label, " is ", format_support_ratio(r), 
                  " times more supported than ", b_label, max_note))
  }
  # r < 1: phrase from the better supported model
  if (r == 0) {
    b_start <- paste0(toupper(substr(b_label, 1, 1)), substring(b_label, 2))
    return(paste0(b_start, " is infinitely more supported than ", sQuote(a),
                  " (which has an IC weight of 0)"))
  }
  paste0(a_label, " has less support than ", b_label, ": ", b_label, " is ",
         format_support_ratio(1 / r), " times more supported than ", sQuote(a),
         max_note)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] R5: new function mark_best_hypo(): '(best)' for printed ratio matrices only
# add " (best)" to the name of the best hypothesis (highest GORIC(A) weight);
# only used for the printed copies of the ratio matrices.
mark_best_hypo <- function(names_ratio, x) {
  best <- x$best_hypo
  if (is.null(best)) {
    best <- which.max(x$result[, 7])
  }
  if (length(best) == 1L && best >= 1L && best <= length(names_ratio)) {
    names_ratio[best] <- paste0(names_ratio[best], " (best)")
  }
  names_ratio
}
# [/CHANGE 2026-10]
