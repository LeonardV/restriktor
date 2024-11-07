print.evSyn <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  # added or equal approach
  #type <- x$type
  
  cat(sprintf("restriktor (%s): %s Evidence Synthesis results:\n", 
              packageDescription("restriktor", fields = "Version"), x$type))
  
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
