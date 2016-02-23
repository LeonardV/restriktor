print.conLM <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nRestriktor: constrained linear model:\n\n")
  if (length(coef(x)) > 0L) {
    cat("Coefficients:\n")
    print(coef(x), digits = digits, scientific = FALSE, print.gap = 2L,
          quote = FALSE)
  } else {
    cat("No coefficients\n")
  }  
  invisible(x)
}
