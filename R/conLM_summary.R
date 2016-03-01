summary.conLM <- function(x, digits = max(3, getOption("digits") - 3),
                          bootCIs = TRUE, bty = "basic", level = 0.95, 
                          signif.stars = getOption("show.signif.stars"), ...) {

  # bty = "stud" needs bootstrap variances
  if (bootCIs & !is.null(x$bootout) & !bty %in% c("norm", "basic", "perc", "bca")) {
    stop("bty is invalid.")
  }
  if (bootCIs & !is.null(x$bootout) & (level < 0.5 | level > 1)) {
    stop("invalid confidence level")
  }

  cat("\nRestriktor: constrained linear model:\n\n")

  cat("Residuals:\n")
  resid <- x$residuals
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- quantile(resid)
   names(rq) <- nam
  print(rq, digits, ...)
  cat("\n")

  if (length(x$b.constr) && is.null(x$bootout)) {
    cat("Coefficients:\n")
    se <- attr(x$information, "se")
    information <- MASS:::ginv(x$information)
    std.error <- sqrt(diag(information))
    tval <- ifelse(std.error != 0, x$b.constr/std.error, 0L)
    coefficients <- cbind(x$b.constr, std.error, tval, 2 * pt(abs(tval),
                                                    x$df.residual, lower.tail = FALSE))
    dimnames(coefficients) <- list(names(x$model.org$coefficients),
                                   c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    coefficients[,4][coefficients[,4] < 2e-16] <- 2e-16
#    print(coefficients, digits = digits, scientific = FALSE, print.gap = 2L, quote = FALSE)
    printCoefmat(coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    cat("\n")
    if (se == "defaults") {
      cat("Homoskedastic standard errors\n")
    } else {
      cat("Heteroskedastic standard errors:", se ,"\n")
    }
    cat("Constrained model: R2 reduced from", round(x$R2.org,3), "to", round(x$R2.reduced, 3),"\n")
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
    cat("\nConstrained model: R2 reduced from", round(x$R2.org,3), "to", round(x$R2.reduced, 3),"\n")
  } else {
    cat("No coefficients\n")
  }  
  cat("\n")
  invisible(x)
}
