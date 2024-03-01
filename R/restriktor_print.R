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
  } else if (inherits(x,"conMLM")) {
    cat("restricted multivariate linear model:\n\n")
  }
  
  if (attr(x$wt.bar, "method") == "boot") {
    total_bootstrap_draws <- attr(x$wt.bar, "total_bootstrap_draws")
    wt_boostrap_errors <- attr(x$wt.bar, "error.idx")
    converged <- attr(x$wt.bar, "converged")
    successful_draws <- total_bootstrap_draws - length(wt_boostrap_errors)
    cat("Bootstrap-based penalty term calculation:", paste0(ifelse(converged, "(Converged)", "(Not converged)")), "\n")
    cat("  Number of bootstrap draws:", total_bootstrap_draws, "\n")
    cat("  Number of successful bootstrap draws:", paste0(successful_draws), "\n\n")
  }
  
  coef <- coef(x)
  ny <- ncol(coef(x$model.org))
  if (!is.null(ny) && ny > 1L) {
    ynames <- colnames(coef(x))
    if (is.null(ynames)) {
      lhs <- x$model.org$terms[[2L]]
      if (mode(lhs) == "call" && lhs[[1L]] == "cbind") {
        ynames <- as.character(lhs)[-1L]
      } else {
        ynames <- paste0("Y", seq_len(ny))
      }
    }
    ind <- ynames == ""
    if (any(ind)) 
      ynames[ind] <- paste0("Y", seq_len(ny))[ind]
    colnames(coef) <- ynames
  }
  
  if (length(coef) > 0L) { 
    cat("Coefficients:\n")
    print(coef, digits = digits, scientific = FALSE, print.gap = 2L,
          quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  
  if (inherits(x, "conRLM")) {
    robWeights(x$wgt)
    cat("\n")
  }
  
  message(x$messages$mix_weights)
  
  invisible(x)
}
