# this function is largely based on the summary.orlm() function 
# from the ic.infer package.

summary.conTest.conLM <- function(object, test = "F", ...) {
  
  x <- object
  if (!("conLM" %in% class(x))) {
    stop("x must be of class \"conLM\" or \"conRLM\"")
  }
  ldots <- list(...)
  ldots$type <- NULL
  CALL <- c(list(object = x, test = test), ldots)
  CALL$type <- "global"
  out0 <- do.call("conTest", CALL)
  CALL$type <- "A"
  out1 <- do.call("conTest", CALL)
  CALL$type <- "B"
  out2 <- do.call("conTest", CALL)
  
  cat("\nRestriktor: hypothesis tests (", x$df.residual, "error degrees of freedom ):\n", 
      "\n")
  cat("Global test under the inequality restriktions:", "\n")
  ###
  cat("       Test statistic: ", out0$Ts, ",   p-value: ", 
      if (out0$pvalue < 1e-04) {
        "<0.0001"
      } else {
        format(round(out0$pvalue, 4), nsmall = 4)}, "\n\n", sep = "")
  ###
  cat("Type A test: H0: all restriktions active (=)", "\n", 
      "        vs. HA: at least one restriktion strictly true (>)", "\n")
  cat("       Test statistic: ", out1$Ts, ",   p-value: ", 
      if (out1$pvalue < 1e-04) {
        "<0.0001"
      } else {
        format(round(out1$pvalue, 4), nsmall = 4)}, "\n\n", sep = "")
  ###
  cat("Type B test: H0: all restriktions true", "\n", 
      "        vs. HA: at least one restriktion false", "\n")
  cat("       Test statistic: ", out2$Ts, ",   p-value: ", 
      if (out2$pvalue < 1e-04) {
        "<0.0001"
      } else {
        format(round(out2$pvalue, 4), nsmall = 4)}, "\n\n", sep = "")
  ###
  if (x$neq == 0) {
    CALL$type <- "C"
    out3 <- do.call("conTest", CALL)
    cat("Type C test: H0: at least one restriktion false or active (=)", 
        "\n", "        vs. HA: all restriktions strictly true (>)", "\n")
    cat("       Test statistic: ", out3$Ts, ",   p-value: ", 
        if (out3$pvalue < 1e-04) {
          "<0.0001"
        } else {
          format(round(out3$pvalue, 4), nsmall = 4)}, "\n\n", sep = "")
    cat("Type C test is based on a t-distribution (one-sided),", 
        "\nall other tests are based on mixture of F distributions.\n\n")
  }
  else {
    cat("All tests are based on mixture of F distributions", 
        "\n(Type C test is not applicable because of equality restriktions)\n\n")
  }
  
}