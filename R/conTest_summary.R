# this function is largely based on the summary.orlm() function 
# from the ic.infer package.

summary.conTest.conLM <- function(object, test = "F", ...) {
  
  x <- object
  if (!("conLM" %in% class(x))) {
    stop("x must be of class \"conLM\" or \"conRLM\"")
  }
  
  vnames <- names(x$b.unrestr)
  Amat <- x$constraints
  meq <- x$neq
  bvec <- x$rhs
  
  colnames(Amat) <- vnames
  out.rest <- cbind(Amat, c(rep("   ==", meq), rep("   >=", nrow(Amat) -
                                                       meq)), bvec)
  
  rownames(out.rest) <- paste(1:nrow(out.rest), ":", sep = "")
  
  colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs")
  out.rest <- cbind(rep(" ", nrow(out.rest)), out.rest)
  out.rest[x$iact, 1] <- "A"
  if (nrow(Amat) == meq) {
    out.rest[1:nrow(Amat), 1] <- "A"
  }  
  out.rest <- as.data.frame(out.rest)
    names(out.rest)[1] <- ""
  
  ldots <- list(...)
  ldots$type <- NULL
  CALL <- c(list(object = x, test = test), ldots)
  
  # fit hypothesis tests
  CALL$type <- "global"
  out0 <- do.call("conTest", CALL)
  CALL$type <- "A"
  out1 <- do.call("conTest", CALL)
  CALL$type <- "B"
  out2 <- do.call("conTest", CALL)
  
  OUT <- list()
  
  cat("\nRestriktor: hypothesis tests (", x$df.residual, "error degrees of freedom ):\n")
  if (x$R2.org == x$R2.reduced) {
    cat("\nMultiple R-squared remains", round(x$R2.org, 4),"\n")
  } else {
    cat("\nMultiple R-squared reduced from", round(x$R2.org, 4), "to", round(x$R2.reduced, 4),"\n")  
  }
  cat("\n(rows indicated with an \"A\" are active (=) restriktions)\n")
  print(out.rest, quote = FALSE, scientific = FALSE)
  cat("\nOverview of all available hypothesis tests:\n")
  cat("\nGlobal test: H0: all parameters are restrikted to be equal", "\n", 
      "        vs. HA: at least one restriktion strictly true", "\n")
  cat("       Test statistic: ", out0$Ts, ",   p-value: ", 
      if (out0$pvalue < 1e-04) {
        "<0.0001"
      } else {
        format(out0$pvalue, digits = 4)
      }, "\n\n", sep = "")
  
  OUT$Global$Ts <- out0$Ts 
  OUT$Global$pvalue <- out0$pvalue[1]
  
  ###
  cat("Type A test: H0: all restriktions active (=)", "\n", 
      "        vs. HA: at least one inequality restriktion strictly true", "\n")
  cat("       Test statistic: ", out1$Ts, ",   p-value: ", 
      if (out1$pvalue < 1e-04) {
        "<0.0001"
      } else {
        format(out1$pvalue, digits = 4)
      }, "\n\n", sep = "")
  
  OUT$TPA$Ts <- out1$Ts 
  OUT$TPA$pvalue <- out1$pvalue[1]
  ###
  cat("Type B test: H0: all restriktions true", "\n", 
      "        vs. HA: at least one restriktion violated ", "\n")
  cat("       Test statistic: ", out2$Ts, ",   p-value: ", 
      if (out2$pvalue < 1e-04) {
        "<0.0001"
      } else {
        format(out2$pvalue, digits = 4)
      }, "\n\n", sep = "")
  
  OUT$TPB$Ts <- out2$Ts 
  OUT$TPB$pvalue <- out2$pvalue[1]
  
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
          format(out3$pvalue, digits = 4)
        }, "\n\n", sep = "")
    cat("Note: Type C test is based on a t-distribution (one-sided),", 
        "\n      all other tests are based on mixture of F-distributions.\n\n")
  
    OUT$TPC$Ts <- out3$Ts 
    OUT$TPC$pvalue <- out3$pvalue[1]
  }
  else {
    cat("Note: All tests are based on mixture of F-distributions", 
        "\n      (Type C test is not applicable because of equality restriktions)\n\n")
  }
  
  
  invisible(OUT)
}