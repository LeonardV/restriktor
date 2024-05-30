print.summary.restriktor <- function(x, digits = max(3, getOption("digits") - 2),
                                     signif.stars = getOption("show.signif.stars"), ...) {
  
  resid <- x$residuals
  rdf <- x$rdf
  se.type <- x$se.type[1]
  bootCIs <- attr(x$se.type, "bootCIs")
  bty <- attr(x$se.type, "bty")
  level <- attr(x$se.type, "level")
  
  cat("\nCall:\n", gsub("\\n", "\n", paste(deparse(x$call), collapse="\n"), 
                        fixed = TRUE), "\n\n", sep = "")
  
  if (inherits(x, "summary.conRLM")) {
    cat("Restriktor: restricted robust linear model:\n\n")
  } else if (inherits(x, "summary.conGLM")) {
    cat("Restriktor: restricted generalized linear model:\n\n")
  } else if (inherits(x, "summary.conLM")) {
    cat("Restriktor: restricted linear model:\n\n")
  } else if (inherits(x, "summary.conMLM")) {
    cat("Restriktor: restricted multivariate linear model:\n\n")
  }  
  
  cat(if (!is.null(x$weights) && diff(range(x$weights))) 
    "Weighted ", "Residuals:\n", sep = "")
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(c(resid))) == 2L) {
      structure(apply(t(resid), 1L, quantile), dimnames = list(nam, 
                                                               dimnames(resid)[[2L]]))
    } else {
      zz <- zapsmall(quantile(resid), digits + 1L)
      structure(zz, names = nam)
    }
    print(rq, digits = digits, ...)
  } else if (rdf > 0L) {
    print(resid, digits = digits, ...)
  }
  
  coefs <- x$coefficients
  if (se.type %in% c("boot.model.based", "boot.standard") && bootCIs) {
    cat("\nCoefficients from restricted model\nwith",
      100 * level, "pct bootstrap confidence intervals (",bty,"):\n ")  
  } else {
    cat("\nCoefficients:\n")
  } 
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
               na.print = "NA")
  
  if (inherits(x, c("summary.conLM", "summary.conRLM"))) {
    if (inherits(x, "summary.conLM")) {
      stddev <- sqrt(x$s2)
    } else {
      stddev <- x$stddev
    }
    cat("\nResidual standard error:", format(signif(stddev, 
                                                    digits)), "on", rdf, "degrees of freedom")
  }
  
  if (se.type == "standard") {
    cat("\nStandard errors:", se.type ,"\n")
  } else if (se.type == "const") {
    cat("\nHomoskedastic standard errors.\n")
  } else if (se.type %in% c("boot.model.based", "boot.standard")) {
    if (se.type == "boot.model.based") {
      se.type <- "model-based"
    } else if (se.type == "boot.standard") {
      se.type <- "standard"
    }
    cat("\nBootstrapped standard errors:", se.type ,"\n")
  } else {
    cat("\nHeteroskedastic robust standard errors:", se.type ,"\n")
  }
  
  if (!(inherits(x, c("summary.conGLM", "summary.conMLM"))) && !is.na(x$R2.reduced)) {
    if (all((x$R2.org - x$R2.reduced) < 1e-08)) {
     cat("Multiple R-squared remains", sprintf("%5.3f", x$R2.org),"\n")
    } else {
     cat("Multiple R-squared reduced from", sprintf("%5.3f", x$R2.org), "to", 
         sprintf("%5.3f", x$R2.reduced),"\n")  
    }
  }
  
  if (inherits(x, "summary.conGLM")) {
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
      format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null", 
                                                          "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("deviance.null", 
                                                                                                                           "deviance")]), digits = max(5L, digits + 1L)), " on", 
                                           format(unlist(x[c("df.residual.null", "rdf")])), " degrees of freedom\n"), 
                                     1L, paste, collapse = " "), sep = "")
  }
  
  
  if (inherits(x, "summary.conRLM")) {
    cat("\n")
    robWeights(x$wgt)
    cat("\n")
  }
  
  goric <- x$goric
  if (!is.null(goric)) {
    ll <- attr(goric, "loglik") 
    PT <- attr(goric, "penalty")
    result.goric <- c(ll, PT, goric)
      names(result.goric) <- c("Loglik", "Penalty", paste0(attr(goric, "type"))) 
    if (attr(goric, "type") == "goric") {
      cat("\nGeneralized order-restricted information criterion: \n")
    } else if (attr(goric, "type") == "goricc") {
      cat("\nSmall sample generalized order-restricted information criterion: \n")
    } else if (attr(goric, "type") == "gorica") {
      cat("\nGeneralized order-restricted information criteron approximation: \n")
    } else if (attr(goric, "type") == "goricac") {
      cat("\nSmall sample generalized order-restricted information criteron approximation: \n")
    } 
    print(result.goric, digits = digits)
  }  
  cat("\n")
  message(x$messages$mix_weights)
  cat("\n")
  
  invisible(x)
}
