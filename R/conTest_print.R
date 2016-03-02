print.conTest <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  if (!("conTest" %in% class(x))) {
    stop("x must be of class \"conTest\"")
  }
  if (nrow(x$Amat) != x$meq) {
    cat("\nConstrained hypothesis test type", x$type, "\n\n")
    out.test <- c(x$Ts, if (x$pvalue < 1e-04) { "<0.0001" } else { x$pvalue })
      names(out.test) <- c(" Test statistic", "p-value")
  } else if (x$test == "F") {
    cat("\nConstrained hypothesis test type\n\n")
    out.test <- c(x$Ts, x$df, x$df.residual, if (x$pvalue < 1e-04) { "<0.0001" } else { x$pvalue})
      names(out.test) <- c(" Test statistic", "df", "df.residual", "p-value")
  } else if (x$test == "Wald") {
    cat("\nConstrained hypothesis test type\n\n")
    out.test <- c(x$Ts, x$df, if (x$pvalue < 1e-04) { "<0.0001" } else { x$pvalue})
    names(out.test) <- c(" Test statistic", "df", "p-value")
  } 
  print(out.test, digits = digits)
}


# add overall test
#summary.iclm <- function() {
#
#a
#}
