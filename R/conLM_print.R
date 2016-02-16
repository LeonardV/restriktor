print.conLM <- function(x, digits = max(3, getOption("digits") - 3),
                        bootCIs = TRUE, bty = "basic", level = 0.95, ...) {

  # bty = "stud" needs bootstrap variances
  if (bootCIs & !is.null(x$bootout) & !bty %in% c("norm", "basic", "perc", "bca")) {
    stop("bty is invalid.")
  }
  if (bootCIs & !is.null(x$bootout) & (level < 0.5 | level > 1)) {
    stop("invalid confidence level")
  }

  cat("\nRestriktor: constrained linear model:\n\n")

  cat("Residuals:\n")
  sr <- summary(c(x$residuals))
  print(sr[c(1,2,3,5,6)], digits = digits, scientific = FALSE, print.gap = 2L,
                quote = FALSE)
  cat("\n")

  if (length(x$b.constr) && is.null(x$bootout)) {
    cat("Coefficients:\n")
    tval <- ifelse(x$se != 0, x$b.constr/x$se, 0L)
    coefficients <- cbind(x$b.constr, x$se, tval, 2 * pt(abs(tval),
                                                    x$df.residual, lower.tail = FALSE))
    dimnames(coefficients) <- list(names(x$model.org$coefficients),
                                   c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    print(coefficients, digits = digits, scientific = FALSE, print.gap = 2L, quote = FALSE)
    cat("\nConstrained model: R2 reduced from", x$R2.org, "to", x$R2.reduced,"\n")
  } else if (length(x$b.constr) && !is.null(x$bootout)) {
    if (bootCIs) {
      cis <- matrix(0, length(x$b.constr), 2)
      colnames(cis) <- c("lower", "upper")
      for (i in 1:length(x$b.constr)) {
        if (!bty %in% c("norm", "perc"))
          cis[i, ] <- boot.ci(x$bootout, conf = level,
                              type = bty, index = i)[[bty]][4:5]
        if (bty == "perc")
          cis[i, ] <- boot.ci(x$bootout, conf = level,
                              type = bty, index = i)[["percent"]][4:5]
        if (bty == "norm")
          cis[i, ] <- boot.ci(x$bootout, conf = level,
                              type = bty, index = i)[["normal"]][2:3]
      }
      cat("\nCoefficients from constrained model\nwith",
          100 * level, "pct bootstrap confidence intervals (",bty,"):", "\n")
    }

    se <- apply(x$bootout$t, 2, sd)
    est <- round(coef(x), 9)
    icc <- cbind(est, se, round(cis, 9))
    colnames(icc) <- c("Estimate", "Std. Error", "Lower", "Upper")
    print(icc, quote = FALSE, digits = digits)

#    print(x$bootout, digits = digits, scientific = FALSE)
    cat("\nConstrained model: R2 reduced from", x$R2.org, "to", x$R2.reduced,"\n")
  } else {
    cat("No coefficients\n")
  }  
  cat("\n")
  invisible(x)
}
