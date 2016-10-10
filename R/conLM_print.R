print.conLM <- function(x, digits = max(3, getOption("digits") - 2), 
                        ...) {

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  
  if (class(x)[1] == "conLM") {
    cat("\nRestriktor: restrikted linear model:\n\n")
  } else if (class(x)[1] == "conRLM") {
    cat("\nRestriktor: restrikted robust linear model:\n\n")
  }
  
  if (length(coef(x)) > 0L) {
    cat("Coefficients:\n")
    print(coef(x), digits = digits, scientific = FALSE, print.gap = 2L,
          quote = FALSE)
  } else {
    cat("No coefficients\n")
  }  
  invisible(x)
}