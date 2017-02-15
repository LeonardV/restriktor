print.summary.restriktor <- function(x, digits = max(3, getOption("digits") - 2),
                                     signif.stars = getOption("show.signif.stars"), ...) {
  
  resid <- x$residuals
  rdf <- x$rdf
  se.type <- x$se.type[1]
  bootCIs <- attr(x$se.type, "bootCIs")
  bty <- attr(x$se.type, "bty")
  level <- attr(x$se.type, "level")
  
  #cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("\nCall:\n", gsub("\\n", "\n", paste(deparse(x$call), collapse="\n"), 
                        fixed = TRUE), "\n\n", sep = "")
  
  if (inherits(x, "summary.conRLM")) {
    cat("Restriktor: restricted robust linear model:\n\n")
  } else if (inherits(x, "summary.conGLM")) {
    cat("Restriktor: restricted generalized linear model:\n\n")
  } else if (inherits(x, "summary.conLM")) {
    cat("Restriktor: restricted linear model:\n\n")
  } 
  
  cat(if (!is.null(x$weights) && diff(range(x$weights))) 
    "Weighted ", "Residuals:\n", sep = "")
  if (rdf > 5L) {
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- if (length(dim(resid)) == 2L) {
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
  if (se.type %in% c("boot.model.based", "boot.standard") & bootCIs) {
    cat("\nCoefficients from restricted model\nwith",
      100 * level, "pct bootstrap confidence intervals (",bty,"):\n ")  
  } else {
    cat("\nCoefficients:\n")
  } 
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
               na.print = "NA")
  
  if (inherits(x, c("summary.conLM", "summary.conRLM"))) {
    cat("\nResidual standard error:", format(signif(sqrt(x$s2.restr), 
                                                    digits)), "on", rdf, "degrees of freedom")
  } 
  
  #cat("\n")
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
  
  if (!(inherits(x, "summary.conGLM"))) { 
    if ((x$R2.org - x$R2.reduced) < 1e-08) {
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
  
  goric <- x$goric
  if (!is.null(goric)) {
    ll <- attr(goric, "loglik") 
    PT <- attr(goric, "penalty")
    result.goric <- c(ll, PT, goric)
    names(result.goric) <- c("Loglik", "Penalty", "Goric")
    cat("\nGeneralized Order-Restricted Information Criterion:\n")
    print(result.goric, digits = digits)
  }  
  cat("\n")
  invisible(x)
}
