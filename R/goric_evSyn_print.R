print.evSyn <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  cat(sprintf("restriktor (%s): %s Evidence Synthesis results:\n", 
              packageDescription("restriktor", fields = "Version"), x$type))
  
  # Display input type information
  cat("\nInput type detected: ")
  if (inherits(x, "evSyn_gorica")) {
    cat("Parameter estimates and covariance matrix inherited from gorica object\n")
  } else if (inherits(x, "evSyn_escalc")) {
    cat("Parameter estimates and covariance matrix inherited from escalc object\n")
  } else if (inherits(x, "evSyn_est")) {
    cat("Parameter estimates and covariance matrix\n")
  } else if (inherits(x, "evSyn_LL")) {
    cat("Log-likelihood and penalty values\n")
  } else if (inherits(x, "evSyn_ICvalues")) {
    cat("Information criteria values\n")
  } else if (inherits(x, "evSyn_ICweights")) {
    cat("Information criteria weights (summing to 1)\n")
  } else if (inherits(x, "evSyn_ICratios")) {
    cat("Ratio of information criteria weights (each vector ends with 1)\n")
  } 
  
  if (!is.null(x$Cumulative_GORICA_weights)) {
    cat("\nFinal GORICA weights:\n") 
    cgw <- sapply(x$Cumulative_GORICA_weights["Final", , drop = FALSE], 
                  FUN = function(x) format_numeric(x, digits = digits))
    names(cgw) <- colnames(x$Cumulative_GORICA_weights)
    print(cgw, print.gap = 2, quote = FALSE, right = TRUE)
    cat("---\n")
  }
  
  cat("\nRatio final GORICA weights:\n")  
  formatted_weights <- apply(x$Final_ratio_GORICA_weights, c(1, 2), 
                             function(val) format_numeric(val, digits = digits))
  print(formatted_weights, print.gap = 2, quote = FALSE, right = TRUE)
  cat("\n")
  
  message(x$messages$mix_weights)
}
