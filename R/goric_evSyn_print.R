print.evSyn <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  # Determine label: default to "GORICA" unless type is explicitly given
  # Use [[ ]] to avoid partial matching (e.g., $type matching $type_ev)
  if (!is.null(x[["type"]])) {
    type_label <- switch(
      x$type,
      goric   = "GORIC",
      goricc  = "GORICC",
      gorica  = "GORICA",
      goricac = "GORICAC"
    )
  } else {
    type_label <- "GORICA"
  }
  
  cat(sprintf("restriktor (%s): %s Evidence Synthesis results:\n", 
              packageDescription("restriktor", fields = "Version"), x[["type_ev"]]))
  
  if (!is.null(x[["messageAdded"]])) {
    cat(x[["messageAdded"]])
    # TO DO
  }
  
  # Display input type information
  if (inherits(x, "evSyn_gorica")) {
    cat(paste("\nInput type 'gorica' detected: "))
    cat("Parameter estimates and covariance matrix inherited from gorica object\n")
  } else if (inherits(x, "evSyn_escalc")) {
    cat(paste("\nInput type 'escalc' detected: "))
    cat("Parameter estimates and covariance matrix inherited from escalc object\n")
  } else if (inherits(x, "evSyn_est")) {
    cat(paste("\nInput type 'est_vcov' detected: "))
    cat("Parameter estimates and covariance matrix\n")
  } else if (inherits(x, "evSyn_LL")) {
    cat(paste("\nInput type 'll_pt' detected: "))
    cat("Log-likelihood and penalty values\n")
  } else if (inherits(x, "evSyn_ICvalues")) {
    cat(paste("\nInput type 'icvalues' detected: "))
    cat("Information criteria values\n")
  } else if (inherits(x, "evSyn_ICweights")) {
    cat(paste("\nInput type 'icweights' detected: "))
    cat("Information criteria weights (summing to 1)\n")
  } else if (inherits(x, "evSyn_ICratios")) {
    cat(paste("\nInput type 'icratios' detected: "))
    #cat("Ratio of information criteria weights (each vector ends with 1)\n") # Note: not necessary!
    cat(paste0("Ratio of information criteria weights (each ", 
               .evSyn_unit_label(x[["type_ev"]], capital = FALSE), " has at least one 1)\n"))
  } 
  
  cat("\n")
  cat("---\n")
  
  # Prior IC weights (priorICweights): shown when they are unequal (i.e., specified by the user)
  prior_w <- x[["priorICweights"]]
  show_prior <- !is.null(prior_w) && !is.null(x[["Cumulative_GORICA_weights"]]) &&
    length(prior_w) == ncol(x[["Cumulative_GORICA_weights"]]) &&
    length(unique(round(prior_w, 12))) > 1
  if (show_prior) {
    pr_fmt <- sapply(prior_w, format_numeric, digits = digits)
    pr_name <- paste0("prior ", type_label, " weights")   # e.g., prior GORICA weights
    prior_note <- paste0("(", pr_name, ": the priorICweights, which are used in the ", type_label,
                         if (isTRUE(x[["priorICweights_study"]])) " weights (study-specific, cumulative, and final) and their ratios, but not in the " else
                           " weights (cumulative and final) and their ratios, but not in the study-specific ones (nor in the measures of variability based on these), the ",
                         type_label, " values, the log-likelihood (LL) and penalty values, and the LL weights.)")
  }
  
  if (!is.null(x[["Cumulative_LL_weights"]]) && !is.null(x[["Cumulative_GORICA_weights"]])) {
    cat(paste0("\nFinal log-likelihood (LL) & ", type_label,  " weights:\n"))
    clw <- sapply(x[["Cumulative_LL_weights"]]["Final", , drop = FALSE], 
                  FUN = function(x) format_numeric(x, digits = digits))
    cgw <- sapply(x[["Cumulative_GORICA_weights"]]["Final", , drop = FALSE], 
                  FUN = function(x) format_numeric(x, digits = digits))
    cw <- t(rbind(clw, cgw))
    rownames(cw) <- colnames(x[["Cumulative_LL_weights"]])
    #colnames(cw) <- c("log-likelihood", type_label)
    colnames(cw) <- c("LL weights", paste0(type_label, " weights"))
    # Prior IC weights (priorICweights; shown when unequal): column right before the IC weights
    if (show_prior) {
      pr_m <- matrix(pr_fmt, ncol = 1, dimnames = list(NULL, pr_name))
      cw <- cbind(cw[, 1, drop = FALSE], pr_m, cw[, -1, drop = FALSE])
    }
    # Confidence intervals for the final IC weights (if evSyn(..., CI = TRUE); see evSyn_CI())
    ci <- x[["CI"]]
    if (!is.null(ci)) {
      idx <- match(rownames(cw), ci$weights$model)
      if (!anyNA(idx)) {
        cw <- cbind(cw, 
                    lower = sapply(ci$weights$lower[idx], format_numeric, digits = digits),
                    upper = sapply(ci$weights$upper[idx], format_numeric, digits = digits))
      }
    }
    print(cw, print.gap = 2, quote = FALSE, right = TRUE)
    if (show_prior) cat(paste0("\n", prior_note, "\n"))
    if (!is.null(ci)) {
      cat(paste0("\n", .evSyn_CI_note(ci), "\n"))
    }
    cat("\n")
    cat("---\n")
    #cat("\n")
  } else if (!is.null(x[["Cumulative_LL_weights"]])) {
    cat(paste0("\nFinal ", "log-likelihood  (LL)", " weights:\n"))
    clw <- sapply(x[["Cumulative_LL_weights"]]["Final", , drop = FALSE], 
                  FUN = function(x) format_numeric(x, digits = digits))
    # names(clw) <- colnames(x[["Cumulative_LL_weights"]])
    # print(clw, print.gap = 2, quote = FALSE, right = TRUE)
    cw <- matrix(clw, ncol = 1)
    cw <- as.data.frame(cw)
    rownames(cw) <- colnames(x[["Cumulative_LL_weights"]])
    colnames(cw) <- "LL weights" #NULL #paste0("Final ", "log-likelihood", " weights")
    print(cw, print.gap = 2, quote = FALSE, right = TRUE, col.names = FALSE)
    #cat("---\n")
    cat("\n")
  } else if (!is.null(x[["Cumulative_GORICA_weights"]])) {
    cat(paste0("\nFinal ", type_label, " weights:\n"))
    cgw <- sapply(x[["Cumulative_GORICA_weights"]]["Final", , drop = FALSE], 
                  FUN = function(x) format_numeric(x, digits = digits))
    # names(cgw) <- colnames(x[["Cumulative_GORICA_weights"]])
    # print(cgw, print.gap = 2, quote = FALSE, right = TRUE)
    cw <- matrix(cgw, ncol = 1)
    cw <- as.data.frame(cw)
    rownames(cw) <- colnames(x[["Cumulative_GORICA_weights"]])
    colnames(cw) <- paste0(type_label, " weights") #NULL #paste0("Final ", type_label, " weights")
    if (show_prior) {
      pr_df <- data.frame(pr_fmt, stringsAsFactors = FALSE); names(pr_df) <- pr_name
      cw <- data.frame(pr_df, cw, check.names = FALSE, row.names = rownames(cw))
    }
    print(cw, print.gap = 2, quote = FALSE, right = TRUE, col.names = FALSE)
    if (show_prior) cat(paste0("\n", prior_note, "\n"))
    #write.table(format(cw, justify="right"), col.names = FALSE, quote = F, sep = "              ")
    cat("\n")
    cat("---\n")
  }
  
  cat(paste0("\nRatio final ", type_label, " weights:\n"))
  formatted_weights <- apply(x[["Final_ratio_GORICA_weights"]], c(1, 2), 
                             function(val) format_numeric(val, digits = digits))
  print(formatted_weights, print.gap = 2, quote = FALSE, right = TRUE)
  cat("\n")
  
  # Measures of variability over the studies (or analyses), incl. conclusion 
  # (as in the summary; only when requested, i.e., evSyn(..., variability = TRUE), 
  # and with at least two studies/analyses; see also evSyn_variability())
  if (isTRUE(x[["show_variability"]]) && !is.null(x[["n_studies"]]) && x[["n_studies"]] > 1) {
    cat("---\n")
    v <- suppressMessages(evSyn_variability(x))
    .evSyn_print_variability(v, type_ev = x[["type_ev"]], type_label = type_label, digits = digits)
    cat("\n")
  }
  
  message(x[["messages"]]$mix_weights)
  
  if (!is.null(x[["penalty_factor"]]) && x[["penalty_factor"]] != 2) {
    message(sprintf("\nrestriktor Message: Note that a penalty factor of %s (default is 2) is used in the calculation of the %s values, that is -2 x log-likelihood + penalty_factor x penalty",
                    x[["penalty_factor"]], type_label))
  }
}
