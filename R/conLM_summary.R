summary.conLM <- function(x, digits = max(3, getOption("digits") - 2),
                          bootCIs = TRUE, bty = "basic", level = 0.95, 
                          signif.stars = getOption("show.signif.stars")) {

  # bty = "stud" needs bootstrap variances
  if (bootCIs & !is.null(x$bootout) & !bty %in% c("norm", "basic", "perc", "bca")) {
    stop("bty is invalid.")
  }
  if (bootCIs & !is.null(x$bootout) & (level < 0.5 | level > 1)) {
    stop("invalid confidence level")
  }

  if (class(x)[1] == "conLM") {
    cat("\nRestriktor: constrained linear model:\n\n")
  } else if (class(x)[1] == "conRLM") {
    cat("\nRestriktor: constrained robust linear model:\n\n")
  }
  
  cat("Residuals:\n")
  resid <- x$residuals
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- quantile(resid)
   names(rq) <- nam
  print(rq, digits)
  cat("\n")

  se <- x$se
  if (length(x$b.constr) && is.null(x$bootout) && !(se == "no")) {
    cat("Coefficients:\n")
    vcovHC <- sandwich(x, bread.=bread(x), meat.=meatHC(x, type = se))
    SE <-
#    if (se == "const" | se == "default") {
#      sqrt(diag(x$information.inverted))
#    } else {
      sqrt(diag(vcovHC))  
#    }  
    
    ##########
    tval <- ifelse(SE != 0, x$b.constr/SE, 0L)
    coefficients <- cbind(x$b.constr, SE, tval, 2 * pt(abs(tval),
                                                    x$df.residual, lower.tail = FALSE))
    dimnames(coefficients) <- list(names(x$model.org$coefficients),
                                   c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    
    ######################### new defined parameters ############################
    if (!(is.null(x$partable)) && any(x$partable$op == ":=")) {
      b.def <- x$CON$def.function(coef(x))
      JAC <- lav_func_jacobian_simple(func = x$CON$def.function, x = coef(x))
      def.cov <- JAC %*% vcovHC %*% t(JAC)
      diag.def.cov <- diag(def.cov)
      diag.def.cov[ diag.def.cov < 0 ] <- as.numeric(NA)
      SE.def <- sqrt(diag.def.cov)
      tval.def <- ifelse(SE.def != 0, b.def/SE.def, 0L)
      coefficients <- rbind(coefficients, cbind(b.def, SE.def, tval.def, 2 * pt(abs(tval.def),
                            x$df.residual, lower.tail = FALSE)))
    }  
    ############################################################################
      coefficients[,4][coefficients[,4] < 2e-16] <- 2e-16
    printCoefmat(coefficients, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA")
     
    cat("\n")
    if (se == "const") {
      cat("Homoskedastic standard errors\n")
    } else {
      cat("Heteroskedastic robust standard errors:", se ,"\n")
    }
    
    if (round(x$R2.org,3) != round(x$R2.reduced, 3)) {
      cat("Constrained model: R2 reduced from", round(x$R2.org,3), "to", round(x$R2.reduced, 3),"\n")  
    } else {
      cat("Constrained model: R2 remains", round(x$R2.org,3),"\n")
    }
    
#    if (class(x)[1] == "conRLM") {
#      cat("Convergence in ")
#    }

  } else if (length(x$b.constr) && !is.null(x$bootout)) {
    if (bootCIs) {
      cis <- matrix(0, length(x$b.constr), 2)
      colnames(cis) <- c("lower", "upper")
      for (i in 1:length(x$b.constr)) {
        if (!bty %in% c("norm", "perc")) {
          cis[i, ] <- boot.ci(x$bootout, conf = level,
                              type = bty, index = i)[[bty]][4:5]
        } else if (bty == "perc") {
          cis[i, ] <- boot.ci(x$bootout, conf = level,
                              type = bty, index = i)[["percent"]][4:5]
        } else if (bty == "norm") {
          cis[i, ] <- boot.ci(x$bootout, conf = level,
                              type = bty, index = i)[["normal"]][2:3]
        }  
      }
      cat("\nCoefficients from constrained model\nwith",
          100 * level, "pct bootstrap confidence intervals (",bty,"):", "\n")
    }

    se <- apply(x$bootout$t, 2, sd)
    est <- round(coef(x), 9)
    icc <- cbind(est, se, round(cis, 9))
    colnames(icc) <- c("Estimate", "Std. Error", "Lower", "Upper")
    print(icc, quote = FALSE, digits = digits)

    cat("\nConstrained model: R2 reduced from", round(x$R2.org,3), "to", round(x$R2.reduced, 3),"\n")
  } else if (se == "no") {
    cat("Coefficients:\n")
    print(coef(x), digits = digits, scientific = FALSE, print.gap = 2L,
          quote = FALSE)
    cat("\nConstrained model: R2 reduced from", round(x$R2.org,3), "to", round(x$R2.reduced, 3),"\n")
  } else {
    cat("No coefficients\n")
  }  
  cat("\n")

#  cat("\n")
  # acknowledgment: code taken from the goric package.
  # The goric_penalty() function uses a simulation approach for calculating the
  # level probabilities, while these weights can be calculated using the 
  # multivariate normal distribution function.
  s2unc <- x$s2unc.ml
  X <- model.matrix(x$model.org)[,,drop=FALSE]
  invW <- kronecker(solve(s2unc), t(X) %*% X)
  W <- solve(invW)
  Amat <- fit.con$Amat

  if (!(nrow(x$Amat) == x$meq)) {
    LP <- rev(con_wt(fit.con$Amat%*%W%*%t(fit.con$Amat), meq = x$meq))
    penalty <- 1 + sum( (1:ncol(W)) * c(LP, rep(0, ncol(W)-length(LP))) )
    goric <- -2*(fit.con$loglik - penalty)
    delta <- goric - min(goric)
    goric_weights <- exp(-delta/2) / sum(exp(-delta/2))
    result_goric <- c(fit.con$loglik, penalty, goric = goric, 
                          goric_weights = round(goric_weights,3))
      names(result_goric) <- c("Loglik", "Penalty", "Goric", "Weights")
    cat("Generalized Order-Restricted Information Criterion:\n")
    print(result_goric, digits = digits)
  }
  invisible(x)
}
