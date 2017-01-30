print.restriktor <- function(x, digits = max(3, getOption("digits") - 2), ...) {

  #cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nCall:\n", gsub("\\n", "\n", paste(deparse(x$call), sep = "\n", collapse="\n"), 
                        fixed = TRUE), "\n\n", sep = "")
  
  if (inherits(x, "conRLM")) {
    cat("Restriktor: restricted robust linear model:\n\n")
  } else if (inherits(x,"conGLM")) {
    cat("Restriktor: restricted generalized linear model:\n\n")
  } else if (inherits(x,"conLM")) {
    cat("Restriktor: restricted linear model:\n\n")
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