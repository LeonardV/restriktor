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
  overlap_unique_combinations <- overlap_unique_combinations[!grepl("unconstrained", overlap_unique_combinations)]
  
  overlap_hypo <- gsub("vs\\.", "", overlap_unique_combinations)
  overlap_hypo <- strsplit(overlap_hypo, " ")
  overlap_hypo <- unique(unlist(overlap_hypo))
  overlap_hypo <- overlap_hypo[overlap_hypo != ""]
  
  # check if the log-likelihood of models are equal, this means that the hypotheses overlap (i.e., subset)
  if (length(overlap_hypo) > 0 & !x$Heq) {
    message("---\nNote: Hypotheses ", paste0(sQuote(overlap_hypo), collapse = " and "), 
            " have equal log-likelihood values. This indicates that they overlap and",
            " that the ratio of their GORIC(A) weights reached its maximum. If of", 
            " interest, we recommend evaluating the overlap versus its complement.", 
            " Run vignette(\"Guidelines_GORIC_output\") in your console for more", 
            " information and an example.")
  }

  if (comparison == "unconstrained" && length(df$model) == 2 && 
      # If comparison = "complement" and ceq only, than comparison is set to unconstrained internally.
      nrow(x$objectList[[1]]$constraints) != x$objectList[[1]]$neq) { 
    message("---\nAdvise: Are you certain you wish to assess the order-restricted hypothesis", 
            " in comparison to the unconstrained one, rather than its complement?")
  }

  best_hypo <- which.max(x$result[, 7])
  best_hypo_name <- x$result$model[best_hypo]
  best_hypo_overlap <- best_hypo_name %in% overlap_hypo
  
  if (length(df$model) > 1) {
    cat("\nConclusion:\n")
  }
  
  if (comparison == "complement" && length(overlap_unique_combinations) == 0 && !x$Heq) {
    formatted_numbers <- sprintf("%.3f %.3f", x$result[[7]][1], x$result[[7]][2])
    numbers <- strsplit(formatted_numbers, " ")[[1]]
    if (as.numeric(numbers[1]) / as.numeric(numbers[2]) > 1) {
      support_ratio <- sprintf("%.2f", x$ratio.gw[1, 2])
      cat(paste0("The order-restricted hypothesis ", sQuote(objectnames[1]), 
                 " has ", support_ratio, " times more support than its complement."))
    } else if (as.numeric(numbers[1]) / as.numeric(numbers[2]) < 1) {
      result <- paste(numbers[1], "/", numbers[2], "< 1", sep = " ")
      cat("The order-restricted hypothesis", sQuote(objectnames[1]), "has", result, "times more support than its complement.",
      "That is, the complement has", sprintf("%.2f", as.numeric(numbers[2]) / as.numeric(numbers[1])), "times more support than", sQuote(objectnames[1]))      
    } else {
      result <- paste(numbers[1], "/", numbers[2], "= 1", sep = " ")
      cat("The order-restricted hypothesis", sQuote(objectnames[1]), "and the complement have equal support:", result)      
    }
    # cat(paste0("---\nThe order-restricted hypothesis ", objectname1, 
    #            " has ", support_ratio, " times more support than its complement.\n\n"))
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
        message <- paste0(message, "\n  * ", sQuote(best_hypo_name), " is ", 
                          goric_rw_without_heq_best_hypo[i], 
                          " times more supported than ", sQuote(best_hypos_rest[i]))

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
                        " as it has the highest GORIC(A) weight.")
      
      message <- paste0(message, "\n- Since the order-restricted hypotheses contain", 
                        " the equality-constrained hypothesis (Heq), inspecting", 
                         " the relative support is not meaningful.")
      
      cat(paste0(message, "\n"))
    }
  } else if (comparison == "none" && length(overlap_unique_combinations) == 0 && length(df$model) == 2) {
    class(x$ratio.gw) <- "numeric"
    support_ratio <- sprintf("%.2f", x$ratio.gw[1, 2])
    objectname1 <- sQuote(objectnames[1])
    objectname2 <- sQuote(objectnames[2])
    cat(paste0("The order-restricted hypothesis ", objectname1, 
               " has ", support_ratio, " times more support than ", objectname2, ".\n\n"))
  } else if (comparison == "unconstrained" && length(overlap_unique_combinations) == 0 && length(df$model) == 2) { 
    formatted_numbers <- sprintf("%.3f %.3f", x$result[[7]][1], x$result[[7]][2])
    numbers <- strsplit(formatted_numbers, " ")[[1]]
    if (as.numeric(numbers[1]) / as.numeric(numbers[2]) > 1) {
      result <- paste(numbers[1], "/", numbers[2], "> 1", sep = " ")
      cat("The order-restricted hypothesis", sQuote(objectnames[1]), "has", result, "times more support than the unconstrained.\n\n")
    } else if (as.numeric(numbers[1]) / as.numeric(numbers[2]) < 1) {
      result <- paste(numbers[1], "/", numbers[2], "< 1", sep = " ")
      cat("The order-restricted hypothesis", sQuote(objectnames[1]), "has", result, "times more support than the unconstrained.\n\n")      
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
                        " is the best in the set, as it has the highest GORIC(A) weight.")
      
      # Step 2: if not weak, compare it against all other hypotheses in the set
      message <- paste0(message, "\n- Since ", sQuote(best_hypo_name), " has a higher", 
      " GORIC(A) weight than the unconstrained hypothesis, it is not considered weak.", 
      " We can now inspect the relative support for ", sQuote(best_hypo_name), " against",
      " the other order-restricted hypotheses:")
      
      for (i in seq_along(best_hypos_rest)) {
        if (best_hypos_rest[i] %in% overlap_hypo) {
          # Als er overlap is, voeg toe dat de relatieve support zijn maximum heeft bereikt
          message <- paste0(message, "\n  * ", sQuote(best_hypo_name), " is ", 
                            goric_rw_without_unc_best_hypo[i], 
                            " times more supported than ", sQuote(best_hypos_rest[i]), 
                            " (This relative support reached its maximum, see Note)")
        } else {
          # Als er geen overlap is, geef normale relatieve support
          message <- paste0(message, "\n  * ", sQuote(best_hypo_name), " is ", 
                            goric_rw_without_unc_best_hypo[i], 
                            " times more supported than ", sQuote(best_hypos_rest[i]))
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
    
    goric_rw <- x$ratio.gw
    goric_rw_best_hypo <- goric_rw[best_hypo, ]
    goric_rw_best_hypo <- goric_rw_best_hypo[goric_rw_best_hypo != 1]
    goric_rw_best_hypo <- sapply(goric_rw_best_hypo, format_value)
    best_hypos_rest <- paste(df$model[!df$model %in% best_hypo_name])
    
    message <- ""
    for (i in seq_along(best_hypos_rest)) {
      if (best_hypos_rest[i] %in% overlap_hypo) {
        # Als er overlap is, voeg toe dat de relatieve support zijn maximum heeft bereikt
        message <- paste0(message, "  * ", sQuote(best_hypo_name), " is ", 
                          goric_rw_best_hypo[i], 
                          " times more supported than ", sQuote(best_hypos_rest[i]), 
                          " (This relative support reached its maximum, see Note)")
      } else {
        # Als er geen overlap is, geef normale relatieve support
        message <- paste0(message, "  * ", sQuote(best_hypo_name), " is ", 
                          goric_rw_best_hypo[i], 
                          " times more supported than ", sQuote(best_hypos_rest[i]))
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
  if (x$penalty_factor != 2) {
    message(sprintf("\nrestriktor Message: Note that a penalty factor of %s (default is 2) is used in the calculation of the %s value, that is -2 x log-likelihood + penalty_factor x penalty",
                    x$penalty_factor, x$type))
  }
}

# TO DO print ook de warnings hier nog eens, anders mis je die ms
#       Kijk of messages en/of warnings en zag dan onderaan dat die er zijn en geef dan de code om die te zien!
