print.evSyn <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  cat(sprintf("restriktor (%s): %s Evidence Synthesis results:\n", 
              packageDescription("restriktor", fields = "Version"), x$type_ev))
  
  if (!is.null(x$messageAdded)) {
    cat(x$messageAdded)
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
    cat(paste("\nInput type 'icvalues' detected: "))
    cat("Information criteria icweights (summing to 1)\n")
  } else if (inherits(x, "evSyn_ICratios")) {
    cat(paste("\nInput type 'icratios' detected: "))
    cat("Ratio of information criteria weights (each vector ends with 1)\n")
  } 
  
  if (!is.null(x$Cumulative_GORICA_weights)) {
    if (x$type == "gorica") {
      cat("\nFinal GORICA weights:\n")
    } else if (x$type == "goricac") {
      cat("\nFinal GORICAC weights:\n") 
    }
    cgw <- sapply(x$Cumulative_GORICA_weights["Final", , drop = FALSE], 
                  FUN = function(x) format_numeric(x, digits = digits))
    names(cgw) <- colnames(x$Cumulative_GORICA_weights)
    print(cgw, print.gap = 2, quote = FALSE, right = TRUE)
    cat("---\n")
  }
  
  if (x$type == "gorica") {
    cat("\nRatio final GORICA weights:\n")  
  } else if (x$type == "goricac") {
    cat("\nRatio final GORICAC weights:\n")  
  }
  formatted_weights <- apply(x$Final_ratio_GORICA_weights, c(1, 2), 
                             function(val) format_numeric(val, digits = digits))
  print(formatted_weights, print.gap = 2, quote = FALSE, right = TRUE)
  cat("\n")
  
  message(x$messages$mix_weights)
}
