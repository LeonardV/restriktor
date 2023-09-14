print.con_goric <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  type <- x$type
  comparison <- x$comparison
  dig <- paste0("%6.", digits, "f")
  x2 <- lapply(x$result[-1], sprintf, fmt = dig)
  df <- data.frame(model = x$result$model, x2)
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  
  if (type == "goric") {
    cat("generalized order-restricted information criterion: \n")
  } else if (type == "gorica") {
    cat("generalized order-restricted information criterion approximation:\n")
  } else if (type == "goricc") {
    cat("small sample generalized order-restricted information criterion:\n")
  } else if (type == "goricac") {
    cat("small sample generalized order-restricted information criterion approximation:\n")
  }
  
  wt_bar <- sapply(x$objectList, function(x) attr(x$wt.bar, "method") == "boot")
  # we need a check if the hypothesis is equalities only
  ceq_only <- sapply(x$objectList, function(x) nrow(x$constraints) == x$neq)
  wt_bar <- as.logical(wt_bar * !ceq_only)
  
  if (sum(wt_bar) > 0) {
    wt_method_boot <- x$objectList[wt_bar]
    wt_bootstrap_draws  <- sapply(wt_method_boot, function(x) attr(x$wt.bar, "mix.bootstrap"))
    wt_bootstrap_errors <- lapply(wt_method_boot, function(x) attr(x$wt.bar, "error.idx"))
    max_nchar <- max(nchar(names(wt_method_boot)))
    
    len <- length(wt_method_boot)
    if (len > 0) { 
      cat("\n")
      cat("Level probabilities:\n")
      cat("  Number of requested bootstrap draws", wt_bootstrap_draws[1], "\n")
      for (i in 1:len) {
        #cat("Number of successful bootstrap draws for", names(wt_method_boot)[i], ":", (wt_bootstrap_draws[1] - length(wt_bootstrap_errors[[i]])), "\n")
        cat(paste0("  Number of successful bootstrap draws for ", sprintf(paste0("%", max_nchar, "s"), names(wt_method_boot)[i]), ": ", (wt_bootstrap_draws[1] - length(wt_bootstrap_errors[[i]])), "\n"))
      }
    }
  }
  
  cat("\nResults:\n")
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  
  if (comparison == "complement") {
    objectnames <- as.character(df$model)
    class(x$ratio.gw) <- "numeric"
    cat("---", "\nThe order-restricted hypothesis", sQuote(objectnames[1]), "has", sprintf("%.3f", x$ratio.gw[1,2]), "times more support than its complement.\n\n")
  } else {
    cat("---\n")  
  }
  
  message(x$messages$mix_weights)
  
  invisible(x)
}
