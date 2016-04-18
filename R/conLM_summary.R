#<FIXME> rewrite in order to extract coefficients from summary <\FIXME>

summary.conLM <- function(x, digits = max(3, getOption("digits") - 2),
                          bootCIs = TRUE, bty = "basic", level = 0.95, 
                          signif.stars = getOption("show.signif.stars"), 
                          ICtype = "GORIC", ...) {

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
  if (is.null(weights(x))) {
    resid <- x$residuals  
  } else {
    resid <- sqrt(x$weights) * x$residuals
  }
  
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- quantile(resid)
   names(rq) <- nam
  print(rq, digits)
  cat("\n")

  se <- x$se
  if (length(x$b.constr) && is.null(x$bootout) && !(se == "no")) {
    cat("Coefficients:\n")
    vcovHC <- sandwich(x, bread.=bread(x), meat.=meatHC(x, type = se))
    SE <- sqrt(diag(vcovHC))  
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
      def.cov <- JAC %*% vcovHC %*% t(JAC) #JAC, Amat[meq]
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
      cat("Homoskedastic standard errors.\n")
    } else {
      cat("Heteroskedastic robust standard errors:", se ,"\n")
    }
    
    if (round(x$R2.org,3) != round(x$R2.reduced, 3)) {
      cat("Constrained model: R2 reduced from", round(x$R2.org,3), "to", round(x$R2.reduced, 3),"\n")  
    } else {
      cat("Constrained model: R2 remains", round(x$R2.org,3),"\n")
    }

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
  
  if (class(x)[1] == "conRLM") {
    s2ml.unc <- summary_rlm(x$model.org, ml = TRUE)$stddev^2
  } else if (class(x)[1] == "conLM") {
    s2ml.unc <- x$s2.unc.ml     
  }
  
  # compute goric
  # REF: Kuiper, R.M.; Hoijtink, H.J.A.; Silvapulle, M. J. (2012) 
  # Journal of statistical planning and inference, volume 142, pp. 2454 - 2463
  ## TO DO: add small samples correction 
  Amat <- x$Amat
  meq <- x$meq
  X <- model.matrix(x$model.org)[,,drop=FALSE]
  Y <- as.matrix(x$model.org$model[, attr(x$model.org$terms, "response")])
  invW <- kronecker(solve(s2ml.unc), t(X) %*% X)
  W <- solve(invW)
  
  if (meq < nrow(Amat)) {
    wt <- con_wt(Amat %*% W %*% t(Amat), meq = meq)
  } else {
    wt <- 1
  }
  # construct weight vector
  wtExt <- rep(0L, ncol(W))
  if (meq > 0L) {
    wtExt[(meq+1):(meq + length(wt))] <- wt
  } else {
    wtExt[1:length(wt)] <- wt
  }
  # penalty term
  PT <- 1 + sum( (1:ncol(W)) * rev(wtExt))
  goric <- -2*(x$loglik - PT)
  result_goric <- c(x$loglik, PT, goric)
      names(result_goric) <- c("Loglik", "Penalty", "Goric")
  cat("\nGeneralized Order-Restricted Information Criterion:\n")
  print(result_goric, digits = digits)
  ####
  
  invisible(x)
}
