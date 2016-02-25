print.conTest <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  if (!("conTest" %in% class(x))) {
    stop("x must be of class \"conTest\"")
  }
  object <- x
  cat("\nConstrained hypothesis test type", object$type, "\n\n")
  out.test <- as.numeric(c(sprintf("%.4f", object$Ts), if (object$pvalue < 1e-04) { "<0.0001" } else { sprintf("%.4f", object$pvalue) }))
    names(out.test) <- c(" Test statistic", "p-value")
  print(out.test)
}


# add overall test
#summary.iclm <- function() {
#
#a
#}
