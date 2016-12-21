print.conTest <- function(x, digits = max(3, getOption("digits") - 2), ...) {

  if (!(inherits(x, "conTest"))) {
    stop("x must be of class \"conTest\"")
  }
  
  Amat <- x[[1]]$Amat
  meq  <- x[[1]]$meq
  bvec <- x[[1]]$bvec
  
  if (!is.null(attr(x[[1]]$pvalue, "boot_type")) || length(x) > 1 || "ceq" %in% names(x)) {
  cat("\nRestriktor: restricted hypothesis tests (", x[[1]]$df.residual, "residual degrees of freedom ):\n")
  } else {
    cat("\nRestriktor: restricted hypothesis tests:\n")
  }
  if (length(x) == 1 && !(names(x) %in% c("C"))) {
    if (x[[1]]$boot %in% c("parametric", "model.based")) {
      cat("( Number of successful bootstrap draws:", attr(x[[1]]$pvalue, "R"),")\n")
    }
  }
  if ((x[[1]]$R2.org - x[[1]]$R2.reduced) < 1e-08) {
    cat("\nMultiple R-squared remains", sprintf("%5.3f", x[[1]]$R2.org),"\n\n")
  } else {
    cat("\nMultiple R-squared reduced from", sprintf("%5.3f", x[[1]]$R2.org), "to", 
        sprintf("%5.3f", x[[1]]$R2.reduced),"\n\n")  
  }
  
  vnames <- names(x[[1]]$b.unrestr)
  colnames(Amat) <- vnames
  
  out.rest <- cbind(round(Amat, 4), c(rep("   ==", meq), rep("   >=", nrow(Amat) - 
                                                     meq)), bvec, " ")
  
  rownames(out.rest) <- paste(1:nrow(out.rest), ":", sep = "")
  
  colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs", "active")
  idx <- ncol(out.rest)
  out.rest[x[[1]]$iact, idx] <- TRUE
  # in case of equality constraints only all constraints are active (==)
  if (nrow(Amat) == meq) {
    out.rest[1:nrow(Amat), idx] <- TRUE
  }  
  out.rest <- as.data.frame(out.rest)
  
  
  if (length(x) > 1L) {
    cat("\nConstraint matrix and information about which constraint is active:\n\n")
    print(out.rest, quote = FALSE, scientific = FALSE)
    
    cat("\n\nOverview of all available hypothesis tests:\n")
    
    cat("\nGlobal test: H0: all parameters are restricted to be equal", "\n", 
        "        vs. HA: at least one restriktion strictly is true", "\n")
    cat("       Test statistic: ", x$global$Ts, ",   p-value: ", 
        if (!is.na(x$global$pvalue) && x$global$pvalue < 1e-04) { 
          "<0.0001"
        } else if (!is.na(x$global$pvalue)) { 
          format(x$global$pvalue, digits = 4)
        } else {
          as.numeric(NA)
        }, "\n\n", sep = "")
    ###
    cat("Type A test: H0: all restriktions are active (==)", "\n", 
        "        vs. HA: at least one inequality restriktion is strictly true", "\n")
    cat("       Test statistic: ", x$A$Ts, ",   p-value: ", 
        if (!is.na(x$A$pvalue) && x$A$pvalue < 1e-04) { 
          "<0.0001"
        } else if (!is.na(x$A$pvalue)) { 
          format(x$A$pvalue, digits = 4)
        } else {
          as.numeric(NA)
        }, "\n\n", sep = "")
    ###
    cat("Type B test: H0: all restriktions are true", "\n", 
        "        vs. HA: at least one restriktion is violated ", "\n")
    cat("       Test statistic: ", x$B$Ts, ",   p-value: ", 
        if (!is.na(x$B$pvalue) && x$B$pvalue < 1e-04) { 
          "<0.0001"
        } else if (!is.na(x$B$pvalue)) { 
          format(x$B$pvalue, digits = 4)
        } else {
          as.numeric(NA)
        }, "\n\n", sep = "")
    ###
    if (length(x) == 4 && meq == 0) {
      cat("Type C test: H0: at least one restriktion is false or active (==)", 
          "\n", "        vs. HA: all restriktions are strictly true (>)", "\n")
      cat("       Test statistic: ", x$C$Ts, ",   p-value: ", 
          if (!is.na(x$C$pvalue) && x$C$pvalue < 1e-04) { 
            "<0.0001"
          } else if (!is.na(x$C$pvalue)) { 
            format(x$C$pvalue, digits = 4)
          } else {
            as.numeric(NA)
          }, "\n\n", sep = "")
      cat("Note: Type C test is based on a t-distribution (one-sided),", 
          "\n      all other tests are based on a mixture of F-distributions.\n\n")
    }
    else {
      cat("Note: All tests are based on a mixture of F-distributions", 
          "\n      (Type C test is not applicable because of equality restriktions)\n\n")
    }
  } else {
    x <- x[[1]]
    df.bar <- attr(x$pvalue, "df.bar")
    
    if (!is.na(x$pvalue)) {
      out.test <- c(sprintf("%.4f", x$Ts), 
                    if (x$pvalue < 1e-04) { "<0.0001" } 
                    else { sprintf("%.4f", x$pvalue) })  
    } else {
      out.test <- c(sprintf("%.4f", x$Ts), as.numeric(NA)) 
    }
    names(out.test)[1] <- sprintf(" %s%s", x$test,"-test statistic")
    names(out.test)[2] <- sprintf("%s", "p-value")
    
    if (nrow(x$Amat) > x$meq) {
      if (x$type == "global") {
        cat("\nGlobal test: H0: all parameters are restricted to be equal", "\n", 
            "        vs. HA: at least one restriktion is strictly true", "\n\n")
        print(out.test, quote = FALSE, scientific = FALSE)
        if (!is.null(df.bar)) {
          cat("\nThis test is based on a mixture of F-distributions on", df.bar, 
              "\ndegrees of freedom and", x$df.residual, "residual degrees of freedom.\n\n")
        }
        cat("\nConstraint matrix and information about which constraint is active:\n\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestricted estimate under H0:\n")
        print.default(format(x$b.eqrestr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nrestricted estimate under union of H0 and HA:\n")
        print.default(format(x$b.restr, digits = digits),
                      print.gap = 2, quote = FALSE)
      } else if (x$type == "A") {
        cat("\nType A test: H0: all restriktions are active (==)", "\n", 
            "        vs. HA: at least one inequality restriktion is strictly true", "\n\n")
        print(out.test, quote = FALSE, scientific = FALSE)        
        if (!is.null(df.bar)) {
        cat("\nThis test is based on a mixture of F-distributions on", df.bar, 
            "\ndegrees of freedom and", x$df.residual, "residual degrees of freedom.\n\n")
        } 
        cat("\nConstraint matrix and information about which constraint is active:\n\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestricted estimate under H0:\n")
        print.default(format(x$b.eqrestr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nrestricted estimate under union of H0 and HA:\n")
        print.default(format(x$b.restr, digits = digits),
                      print.gap = 2, quote = FALSE)
      } else if (x$type == "B" && x$meq.alt == 0L) {
        cat("\nType B test: H0: all restriktions are true", "\n", 
            "        vs. HA: at least one restriktion is violated", "\n\n")
        print(out.test, quote = FALSE)
        if (!is.null(df.bar)) {
          cat("\nThis test is based on a mixture of F-distributions on", df.bar, 
              "\ndegrees of freedom and", x$df.residual, "residual degrees of freedom.\n\n")
        }
        cat("\nConstraint matrix and information about which constraint is active:\n\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestricted estimate under H0:\n")
        print.default(format(x$b.restr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nUnrestricted estimate:\n")
        print.default(format(x$b.unrestr, digits = digits),
                      print.gap = 2, quote = FALSE)
      } else if (x$type == "B" && x$meq.alt > 0L) {
        cat("\nType B test: H0: all restriktions are true (>= or =)", "\n",
            "        vs. HA: at least one restriktion is violated (<), some =-restriktions maintained",
              "\n\n")
        print(out.test, quote = FALSE)
        if (!is.null(df.bar)) {
          cat("\nThis test is based on a mixture of F-distributions on", df.bar, 
              "\ndegrees of freedom and", x$df.residual, "residual degrees of freedom.\n\n")
        }
        cat("\nConstraint matrix and information about which constraint is active:\n\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestricted estimate under H0:\n")
        print.default(format(x$b.restr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nrestricted estimate under HA:\n")
        print.default(format(x$b.restr.alt, digits = digits),
                      print.gap = 2, quote = FALSE)
        } else if (x$type == "C") {
          cat("\nType C test: H0: at least one restriktion is false or active (==)", 
              "\n", "        vs. HA: all restriktions are strictly true (>)", "\n\n")
          print(out.test, quote = FALSE)
          cat("\nThis test is based on a one-sided t-distributions on", x$df.residual, 
              "residual \ndegrees of freedom.\n\n")
          cat("\nConstraint matrix and information about which constraint is active:\n\n")
          print(out.rest, quote = FALSE, scientific = FALSE)
          cat("\nunrestricted estimate:\n")
          print.default(format(x$b.unrestr, digits = digits),
                        print.gap = 2, quote = FALSE)
        }
    } else { #equality constraints only
      cat("\n","classical test: H0: all restriktions are active (==)", 
          "\n","            vs. HA: at least one equality restriktion is violated", "\n\n")
      print(out.test, quote = FALSE)
      cat("\n\n(all rows are active restriktions under H0, H1 is unrestricted!)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
      cat("\nrestricted estimate under H0:\n")
      print.default(format(x$b.restr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\nunrestricted estimate:\n")
      print.default(format(x$b.unrestr, digits = digits),
                    print.gap = 2, quote = FALSE)
    }
  }
}

