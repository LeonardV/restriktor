print.conTest <- function(x, digits = max(3, getOption("digits") - 2), ...) {

  if (!(inherits(x, "conTest"))) {
    stop("x must be of class \"conTest\"")
  }
  
  if (length(x) > 5) {
    type <- x$type  
    x <- list(x)
    names(x) <- type
  }
  
  Amat <- x[[1]]$Amat
  meq  <- x[[1]]$meq
  bvec <- x[[1]]$bvec
  #rdf  <- x[[1]]$df.residual
  boot <- x[[1]]$boot
  model.org  <- x[[1]]$model.org
  b.unrestr  <- x[[1]]$b.unrestr
  iact <- x[[1]]$iact
  R    <- attr(x[[1]]$pvalue, "R")
  
  #cat("\nRestriktor: restricted hypothesis tests (", rdf, "residual degrees of freedom ):\n")
  cat("\nRestriktor: restricted hypothesis tests: \n")
  
  if (!("C" %in% names(x))) {
    if (boot %in% c("parametric", "model.based")) {
      cat("( Number of successful bootstrap draws:", R,")\n")
    }
  } else {
    cat("\n")
  }
  
  if (!inherits(model.org, "glm")) {
    R2.reduced <- x[[1]]$R2.reduced
    R2.org     <- x[[1]]$R2.org
    if (all((R2.org - R2.reduced) < 1e-08)) {
      cat("\nMultiple R-squared remains", sprintf("%5.3f", R2.org),"\n")
    } else {
      cat("\nMultiple R-squared reduced from", sprintf("%5.3f", R2.org), "to", 
          sprintf("%5.3f", R2.reduced),"\n")  
    }
  }
  
  colnames(Amat) <- names(b.unrestr)
  out.rest <- cbind(round(Amat, 4), c(rep("   ==", meq), rep("   >=", nrow(Amat) - 
                                                     meq)), bvec, " ")
  
  rownames(out.rest) <- paste(seq_len(nrow(out.rest)), ":", sep = "")
  
  colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs", "active")
  idx <- ncol(out.rest)
  out.rest[, idx] <- "no"
  out.rest[iact, idx] <- "yes"
  # in case of equality constraints only all constraints are active (==)
  if (nrow(Amat) == meq) {
    out.rest[seq_len(nrow(Amat)), idx] <- "yes"
  }  
  out.rest <- as.data.frame(out.rest)
  
  
  if (length(x) > 1L) {
    cat("\nConstraint matrix:\n")
    print(out.rest, quote = FALSE, scientific = FALSE)
    
    cat("\n\nOverview of all available hypothesis tests:\n")
    
    if (!is.null(x$global)) {
      cat("\nGlobal test: H0: all parameters are restricted to be equal (==)\n", 
          "        vs. HA: at least one inequality restriction is strictly true (>)\n")
      cat("       Test statistic: ", sprintf("%.4f", x$global$Ts), ",   p-value: ", 
          if (!is.na(x$global$pvalue) && x$global$pvalue < 1e-04) { 
            "<0.0001"
          } else if (!is.na(x$global$pvalue)) { 
            format(x$global$pvalue, digits = 4)
          } else {
            as.numeric(NA)
          }, "\n\n", sep = "")
    }
    ###
    if (!is.null(x$A)) {
      cat("Type A test: H0: all restrictions are equalities (==)", "\n", 
          "        vs. HA: at least one inequality restriction is strictly true (>)\n")
      cat("       Test statistic: ", sprintf("%.4f", x$A$Ts), ",   p-value: ", 
          if (!is.na(x$A$pvalue) && x$A$pvalue < 1e-04) { 
            "<0.0001"
          } else if (!is.na(x$A$pvalue)) { 
            format(x$A$pvalue, digits = 4)
          } else {
            as.numeric(NA)
          }, "\n\n", sep = "")
    }
    ###
    if (!is.null(x$B)) {
      if (x$B$meq.alt == 0L) {
        cat("Type B test: H0: all restrictions hold in the population\n", 
            "        vs. HA: at least one restriction is violated\n")
      } else if (x$B$meq.alt > 0L) {
        cat("Type B test: H0: all restrictions hold in the population\n", 
            "        vs. HA: at least one restriction is violated (<),", 
            "\n                  some equality restrictions are maintained\n")
      }
      cat("       Test statistic: ", sprintf("%.4f", x$B$Ts), ",   p-value: ", 
          if (!is.na(x$B$pvalue) && x$B$pvalue < 1e-04) { 
            "<0.0001"
          } else if (!is.na(x$B$pvalue)) { 
            format(x$B$pvalue, digits = 4)
          } else {
            as.numeric(NA)
          }, "\n\n", sep = "")
    }
    ###
    if (!is.null(x$C)) {
      cat("Type C test: H0: at least one restriction is false or active (==)", 
          "\n", "        vs. HA: all restrictions are strictly true (>)\n")
      cat("       Test statistic: ", sprintf("%.4f", x$C$Ts), ",   p-value: ", 
          if (!is.na(x$C$pvalue) && x$C$pvalue < 1e-04) { 
            "<0.0001"
          } else if (!is.na(x$C$pvalue)) { 
            format(x$C$pvalue, digits = 4)
          } else {
            as.numeric(NA)
          }, "\n\n", sep = "")
      cat("Note: Type C test is based on a t-distribution (one-sided),", 
          "\n      all other tests are based on a mixture of F-distributions.\n\n")
    } else {
      if (inherits(x[[1]]$model.org, "mlm")) {
        cat("Note: All tests are based on a mixture of F-distributions", 
            "\n      (Type C test is not (yet) available for object of class mlm.)\n\n")  
      } else {
        cat("Note: All tests are based on a mixture of F-distributions", 
            "\n      (Type C test is not applicable because of equality restrictions.)\n\n")
      }
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
        cat("\n\nGlobal test: H0: all parameters are restricted to be equal (==)", "\n", 
            "        vs. HA: at least one inequality restriction is strictly true (>)\n\n")
        print(out.test, quote = FALSE, scientific = FALSE)
        if (!is.null(df.bar)) {
          cat("\nThis test is based on a mixture of F-distributions on", df.bar, 
              "\ndegrees of freedom and", x$df.residual, "residual degrees of freedom.\n\n")
        }
        cat("\nConstraint matrix:\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestricted estimate under H0:\n")
        print.default(format(x$b.eqrestr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nrestricted estimate under HA:\n")
        print.default(format(x$b.restr, digits = digits),
                      print.gap = 2, quote = FALSE)
      } else if (x$type == "A") {
        cat("\nType A test: H0: all restrictions are equalities (==)", "\n", 
            "        vs. HA: at least one inequality restriction is strictly true (>)\n\n")
        print(out.test, quote = FALSE, scientific = FALSE)        
        if (!is.null(df.bar)) {
        cat("\nThis test is based on a mixture of F-distributions on", df.bar, 
            "\ndegrees of freedom and", x$df.residual, "residual degrees of freedom.\n\n")
        } 
        cat("\nConstraint matrix:\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestricted estimate under H0:\n")
        print.default(format(x$b.eqrestr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nrestricted estimate under HA:\n")
        print.default(format(x$b.restr, digits = digits),
                      print.gap = 2, quote = FALSE)
      } else if (x$type == "B" && x$meq.alt == 0L) {
        cat("\nType B test: H0: all restrictions hold in the population", "\n", 
            "        vs. HA: at least one restriction is violated\n\n")
        print(out.test, quote = FALSE)
        if (!is.null(df.bar)) {
          cat("\nThis test is based on a mixture of F-distributions on", df.bar, 
              "\ndegrees of freedom and", x$df.residual, "residual degrees of freedom.\n\n")
        }
        cat("\nConstraint matrix:\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestricted estimate under H0:\n")
        print.default(format(x$b.restr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nUnrestricted estimate:\n")
        print.default(format(x$b.unrestr, digits = digits),
                      print.gap = 2, quote = FALSE)
      } else if (x$type == "B" && x$meq.alt > 0L) {
        cat("\nType B test: H0: all restrictions hold in the population", "\n", 
            "        vs. HA: at least one restriction is violated (<),", 
            "\n                  some equality restrictions are maintained\n\n")
        print(out.test, quote = FALSE)
        if (!is.null(df.bar)) {
          cat("\nThis test is based on a mixture of F-distributions on", df.bar, 
              "\ndegrees of freedom and", x$df.residual, "residual degrees of freedom.\n\n")
        }
        cat("\nConstraint matrix:\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestricted estimate under H0:\n")
        print.default(format(x$b.restr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nrestricted estimate under HA:\n")
        print.default(format(x$b.restr.alt, digits = digits),
                      print.gap = 2, quote = FALSE)
        } else if (x$type == "C") {
          cat("\nType C test: H0: at least one restriction is false or active (==)", 
              "\n", "        vs. HA: all restrictions are strictly true (>)\n\n")
          print(out.test, quote = FALSE)
          cat("\nThis test is based on a one-sided t-distributions on", x$df.residual, 
              "residual \ndegrees of freedom.\n\n")
          cat("\nConstraint matrix:\n")
          print(out.rest, quote = FALSE, scientific = FALSE)
          cat("\nunrestricted estimate:\n")
          print.default(format(x$b.unrestr, digits = digits),
                        print.gap = 2, quote = FALSE)
        }
    } else { #equality constraints only
      cat("\n","classical test: H0: all restrictions are active (==)", 
          "\n","            vs. HA: at least one equality restriction is violated\n\n")
      print(out.test, quote = FALSE)
      cat("\n\n(all rows are active restrictions under H0, H1 is unrestricted!)\n")
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






print.conTestLavaan <- function(x, digits = max(3, getOption("digits") - 2), ...) {
  object <- x
  cat("\nRestriktor: restricted hypothesis tests:\n\n")
  cat("  Variable names in model         :", unlist(object$fit.h2@Data@ov.names[1]), "\n")
  cat("  Number of variables             :", object$fit.h2@Model@nvar[1], "\n")
  cat("  Number of groups                :", object$fit.h2@Data@ngroups, "\n")
  cat("  Used sample size per group      :", unlist(object$fit.h2@Data@nobs), "\n")
  cat("  Used sample size                :", sum(unlist(object$fit.h2@Data@nobs)), "\n")
  cat("  Total sample size               :", sum(unlist(object$fit.h2@Data@norig)), "\n\n")
  cat("  Estimator                       :", object$fit.h2@Options$estimator, "\n")
  cat("  Missing data                    :", object$fit.h2@Options$missing, "\n")
  cat("  Bootstrap method                :", object$type, "\n")
  cat("  Double bootstrap method         :", object$double.bootstrap, "\n")
  
  dbtype <- object$double.bootstrap
  # original LRT for hypothesis test Type A
  TsA <- attr(object$bootA, "D.original")
  # original LRT for hypothesis test Type B
  TsB <- attr(object$bootB, "D.original")
  # unadjusted pvalues for Ts
  pvalueA <- object$bootA[1]
  pvalueB <- object$bootB[1]
  alpha <- object$double.bootstrap.alpha
  
  ###
  if (dbtype == "no") {
    if (!is.null(TsA)) {
    cat("\n\n  Type A test: H0: all restriktions active (=)", "\n",
        "          vs. H1: at least one restriktion strictly true (>)", "\n")
    cat("         Test statistic: ", format(round(TsA, digits), nsmall = digits), ", unadjusted p-value: ",
        if (pvalueA < 1e-04) {
          "<0.0001"
        } else {
          format(round(pvalueA, digits), nsmall = digits)}, " (alpha = ", alpha, ") ", "\n\n", sep = "")
    }
    if (!is.null(TsB)) {
    cat("\n  Type B test: H0: all restriktions true", "\n",
        "          vs. H1: at least one restriktion false", "\n")
    cat("         Test statistic: ", format(round(TsB, digits), nsmall = digits), ", unadjusted p-value: ",
        if (pvalueB < 1e-04) {
          "<0.0001"
        } else {
          format(round(pvalueB, digits), nsmall = digits)}, " (alpha = ", alpha, ") ", "\n", sep = "")
    }
  } else if (dbtype == "standard") {
    # adjusted nominal levels
    adj.alphaA  <- attr(object$bootA, "adj.alpha")
    adj.alphaB  <- attr(object$bootB, "adj.alpha")
    # adjusted pvalues for Ts
    adj.pvalueA <- attr(object$bootA, "adj.pvalue")
    adj.pvalueB <- attr(object$bootB, "adj.pvalue")
    if (!is.null(TsA)) {
      cat("\n\n  Type A test: H0: all restriktions active (=)", "\n",
          "          vs. H1: at least one restriktion strictly true (>)", "\n")
      cat("         Test statistic: ", format(round(TsA, digits),
                                              nsmall = digits), ", adjusted p-value: ",
          if (adj.pvalueA < 1e-04) {
            "<0.0001"
          } else {
            format(round(adj.pvalueA, digits), nsmall = digits)}, " (alpha = ", alpha, ") ", "\n", sep = "")
      cat("                                ", "  unadjusted p-value: ",
          if (pvalueA < 1e-04) {
            "<0.0001"
          } else {
            format(round(pvalueA, digits), nsmall = digits)}, " (alpha = ",
          format(round(adj.alphaA, digits), nsmall = digits), ") ", "\n\n", sep = "")
    }
    if (!is.null(TsB)) {    
      cat("\n  Type B test: H0: all restriktions true", "\n",
          "          vs. H1: at least one restriktion false", "\n")
      cat("         Test statistic: ", format(round(TsB, digits), nsmall = digits), ", adjusted p-value: ",
          if (adj.pvalueB < 1e-04) {
            "<0.0001"
          } else {
            format(round(adj.pvalueB, digits), nsmall = digits)}, " (alpha = ", alpha, ") ", "\n", sep = "")
      cat("                               ", "   unadjusted p-value: ",
          if (pvalueB < 1e-04) {
            "<0.0001"
          } else {
            format(round(pvalueB, digits), nsmall = digits)}, " (alpha = ",
          format(round(adj.alphaB, digits), nsmall = digits), ") ", "\n\n", sep = "")
    }
  } else if (dbtype == "FDB") {
    # adjusted pvalues for Ts
    adj.pvalueA <- attr(object$bootA, "adj.pvalue")
    adj.pvalueB <- attr(object$bootB, "adj.pvalue")
    if (!is.null(TsA)) {
      cat("\n\n  Type A test: H0: all restriktions active (=)", "\n",
          "          vs. H1: at least one restriktion strictly true (>)", "\n")
      cat("         Test statistic: ", format(round(TsA, digits),
                                              nsmall = digits), ", adjusted p-value: ",
          if (adj.pvalueA < 1e-04) {
            "<0.0001"
          } else {
            format(round(adj.pvalueA, digits), nsmall = digits)}, " (alpha = ", alpha, ") ", "\n", sep = "")
    }
    if (!is.null(TsB)) {    
      cat("\n  Type B test: H0: all restriktions true", "\n",
          "          vs. H1: at least one restriktion false", "\n")
      cat("         Test statistic: ", format(round(TsB, digits), nsmall = digits), ", adjusted p-value: ",
          if (adj.pvalueB < 1e-04) {
            "<0.0001"
          } else {
            format(round(adj.pvalueB, digits), nsmall = digits)}, " (alpha = ", alpha, ") ", "\n", sep = "")
    }
  }
  
  if (dbtype == "no") {
    cat("\n  No double bootstrap method is set. The results may be spurious.\n\n")
  }
  
}
