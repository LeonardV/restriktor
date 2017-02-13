print.restriktor <- function(x, digits = max(3, getOption("digits") - 2), ...) {

  #cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nCall:\n", gsub("\\n", "\n", paste(deparse(x$call), sep = "\n", collapse="\n"), 
                        fixed = TRUE), "\n\n", sep = "")
  
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  
  if (inherits(x, "conRLM")) {
    cat("restricted robust linear model:\n")
    if (x$converged) {
      cat("Converged in", length(x$iter), "iterations\n\n")
    } else {
      cat("Ran", length(x$iter), "iterations without convergence\n\n")
    }
  } else if (inherits(x,"conGLM")) {
    cat("restricted generalized linear model:\n")
    if (x$converged) {
      cat("Converged in", x$iter, "iterations\n\n")
    } else {
      cat("Ran", x$iter, "iterations without convergence\n\n")
    }
  } else if (inherits(x,"conLM")) {
    cat("restricted linear model:\n\n")
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
