print.con_goric <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  type <- x$type
  comparison <- x$comparison
  dig <- paste0("%6.", digits, "f")
  x2 <- lapply(x$result[-1], sprintf, fmt = dig)
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
                            "to increase the number of bootstrap draws, i.e., control = list(mix_weights_bootstrap_limit = 1e5)")
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
  # which hypotheses overlap, i.e., have equal likelihood values
  loglik_overlap <- which(loglik_diff_mat == 0, arr.ind = TRUE)
  # get unique combination from row- and columnnames
  overlap_unique_combinations <- apply(loglik_overlap, 1, function(x) {
    get_names <- c(rownames(loglik_diff_mat)[x[1]], gsub("vs. ", "", colnames(loglik_diff_mat)[x[2]]))
    paste(get_names, collapse = " vs. ")
  })
  
  # Create a function to sort the elements in each string
  sort_combination <- function(combination) {
    split_combination <- strsplit(combination, " vs. ")[[1]]
    #sorted_combination <- sort(split_combination)
    # Check if the combination includes "complement"
    if ("complement" %in% split_combination) {
      # Find the other element that is not "complement"
      other_element <- split_combination[split_combination != "complement"]
      sorted_combination <- c(other_element, "complement")
    } else if ("unconstrained" %in% split_combination) {
      # Find the other element that is not "unconstrained"
      other_element <- split_combination[split_combination != "unconstrained"]
      sorted_combination <- c(other_element, "unconstrained")
    } else {
      sorted_combination <- sort(split_combination)
    }
    paste(sorted_combination, collapse = " vs. ")
  }

  overlap_sorted_vector <- sapply(overlap_unique_combinations, sort_combination)
  overlap_unique_combinations <- unique(overlap_sorted_vector)
  
  combined_string <- gsub("vs\\.", "", overlap_unique_combinations)
  combined_string <- strsplit(combined_string, " ")
  combined_string <- unique(unlist(combined_string))
  combined_string <- combined_string[combined_string != ""]
  
  # check if the likelihood of models are equal, this means that the hypotheses overlap (i.e., subset)
  if (length(combined_string) > 0) {
    message("---\nNote: Hypotheses ", paste0(sQuote(combined_string), collapse = " and "), 
            " have equal likelihood values (i.e., the hypotheses overlap). The GORIC(A) weights,", 
            " are influenced only by changes in penalty values and have a maximum limit.",
            " If any overlapping hypothesis has the lowest GORIC(A) value in the entire set,",
            " it is recommended to further evaluate it against its complement.")
  }

  if (comparison == "unconstrained" && length(df$model) == 2) {
    message("---\nAdvise: Are you certain you wish to assess the order-restricted hypothesis", 
            " in comparison to the unconstrained one, rather than its complement?")
  }
  
  if (comparison == "complement" && length(overlap_unique_combinations) == 0) { 
    class(x$ratio.gw) <- "numeric"
    objectname1 <- sQuote(objectnames[1])
    support_ratio <- sprintf("%.2f", x$ratio.gw[1, 2])
    cat(glue("---\nThe order-restricted hypothesis {objectname1} has {support_ratio} times more support than its complement.\n\n"))
  } else if (comparison == "none" && length(overlap_unique_combinations) == 0 && length(df$model) == 2) {
    class(x$ratio.gw) <- "numeric"
    support_ratio <- sprintf("%.2f", x$ratio.gw[1, 2])
    objectname1 <- sQuote(objectnames[1])
    objectname2 <- sQuote(objectnames[2])
    cat(glue("---\nThe order-restricted hypothesis {objectname1} has {support_ratio} times more support than {objectname2}.\n\n"))
  } else if (comparison == "unconstrained" && length(overlap_unique_combinations) == 0 && length(df$model) == 2) { 
    formatted_numbers <- sprintf("%.3f %.3f", x$result[[7]][1], x$result[[7]][2])
    numbers <- strsplit(formatted_numbers, " ")[[1]]
    if (as.numeric(numbers[1]) / as.numeric(numbers[2]) > 1) {
      result <- paste(numbers[1], "/", numbers[2], "> 1", sep = " ")
      cat("---\nThe order-restricted hypothesis", sQuote(objectnames[1]), "has", result, "times more support than the unconstrained.\n\n")
    } else if (as.numeric(numbers[1]) / as.numeric(numbers[2]) < 1) {
      result <- paste(numbers[1], "/", numbers[2], "< 1", sep = " ")
      cat("---\nThe order-restricted hypothesis", sQuote(objectnames[1]), "has", result, "times more support than the unconstrained.\n\n")      
    } else {
      result <- paste(numbers[1], "/", numbers[2], "= 1", sep = " ")
      cat("---\nThe order-restricted hypothesis", sQuote(objectnames[1]), "and the unconstrained have equal support:", result, "\n\n")      
    }
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
}
