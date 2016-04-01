summary.conTest <- function(x, test = "F", digits = max(3, getOption("digits") - 4), ...) {
  
  if (!("conTest" %in% class(x))) {
    stop("x must be of class \"conTest\"")
  }
  
  if (nrow(x$Amat) != x$meq) {
    cat("\nRestriktor: inequality constrained hypothesis tests\n\n")
    
    if (length(x$b.constr) > 0L) {
      cat("Constrained coefficients:\n")
      print(x$b.constr, digits = digits, scientific = FALSE, print.gap = 2L,
            quote = FALSE)
    } else {
      cat("No coefficients\n")
    }
    
    vnames <- names(x$b.constr)
    if (is.null(vnames)) { 
      vnames <- paste("m", 1:length(x$b.constr), sep = "") 
    }
    Amat <- x$Amat
    colnames(Amat) <- vnames
    out.rest <- cbind(Amat, c(rep("   ==", x$meq), rep("   >=", nrow(x$Amat) -
                                                         x$meq)), x$bvec)
    
    rownames(out.rest) <- paste(1:nrow(out.rest), ":", sep = "")
    
    colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs")
    out.rest <- cbind(rep(" ", nrow(out.rest)), out.rest)
    out.rest[x$iact, 1] <- "A"
    out.rest <- as.data.frame(out.rest)
    names(out.rest)[1] <- ""
    cat("\n(rows indicated with an \"A\" are active constraints)\n")
    print(out.rest, quote = FALSE, scientific = FALSE)
    
    cat("\n\nGlobal model test under the inequality constraints:\n")
    ### !!! Exceptions for models without intercept to be implemented
    cat("       Test statistic: ", round(out0$Ts, 9), ",   p-value: ", 
      if (out0$pvalue < 1e-04) 
            "<0.0001"
          else format(round(out0$pvalue, 4), nsmall = 4), 
          "\n\n", sep = "")
    
      cat("Type A test: H0: all constraints active(=)", "\n", 
          "            HA: at least one restriction strictly true (>)", 
          "\n")
      ## weights are used from previous test
      cat("       Test statistic: ", round(outA$Ts, 9), ",   p-value: ", 
          if (outA$pvalue < 1e-04) 
            "<0.0001"
          else format(round(outA$pvalue, 4), nsmall = 4), 
          "\n", sep = "")
      
      cat("Type B test: H0: all constraints true", "\n", 
          "        vs. HA: at least one restriction false", 
          "\n")
      ## weights are used from previous test
      cat("       Test statistic: ", round(outB$T, 9), ",   p-value: ", 
          if (outB$pvalue < 1e-04) 
            "<0.0001"
          else format(round(outB$pvalue, 4), nsmall = 4), 
          "\n", sep = "")
      
      
      if (!x$meq > 0) {
        cat("Type C test: H0: at least one restriction false or active (=)", "\n", 
            "        vs. HA: all constraints strictly true (>)", 
            "\n")
        
        cat("       Test statistic: ", round(outC$Ts, 9), 
            ",   p-value: ", if (outC$pvalue < 1e-04) 
              "<0.0001"
            else format(round(outC$pvalue, 4), nsmall = 4), 
            "\n", sep = "")
        cat("\nType C test based on a one-sided t-distribution \n\n")
      } else {
        cat("\n(Type C test not applicable because of equality constraints)\n\n")
      }

  }
}  
  
  
  
  
  
  
  
#  if (nrow(x$Amat) != x$meq) {
#    cat("\nConstrained hypothesis test type", x$type, "\n\n")
#    out.test <- c(sprintf("%.4f", x$Ts), 
#                  if (x$pvalue < 1e-04) { "<0.0001" } 
#                  else { sprintf("%.4f", x$pvalue) })
#      names(out.test) <- c(" Test statistic", "p-value")
#  } else { #equality constraints only
#    if (x$test == "F") {
#      cat("\nConstrained hypothesis test\n\n")
  #    out.test <- c(x$Ts, x$df, x$df.residual, if (x$pvalue < 1e-04) { "<0.0001" } else { x$pvalue})
#      out.test <- c(sprintf("%.4f", x$Ts), sprintf("%.4f", x$df), sprintf("%.4f", x$df.residual), 
#                               if (x$pvalue < 1e-04) { "<0.0001" } 
#                               else { sprintf("%.4f", x$pvalue) })
#      names(out.test) <- c(" Test statistic", "df", "df.residual", "p-value")
#    } else if (x$test == "Wald") {
#      cat("\nConstrained hypothesis test\n\n")
  #    out.test <- c(x$Ts, x$df, if (x$pvalue < 1e-04) { "<0.0001" } else { x$pvalue})
#      out.test <- c(sprintf("%.4f", x$Ts), sprintf("%.4f", x$df), 
#                               if (x$pvalue < 1e-04) { "<0.0001" } 
#                               else { sprintf("%.4f", x$pvalue) })
#      names(out.test) <- c(" Test statistic", "df", "p-value")
#    }  
#  } 
#  print(out.test, quote = FALSE, digits = digits)
#}

