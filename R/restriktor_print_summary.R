print.summary.restriktor <- function(x, digits = max(3, getOption("digits") - 2),
                                     signif.stars = getOption("show.signif.stars"), ...) {
  
  resid <- x$residuals
  rdf <- x$rdf
  se_type <- x$se_type[1]
  bootCIs <- attr(x$se_type, "bootCIs")
  bty <- attr(x$se_type, "bty")
  level <- attr(x$se_type, "level")
  
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  if (inherits(x, "summary.conRLM")) {
    cat("Restriktor: restricted robust linear model:\n\n")
  } else if (inherits(x, "summary.conLM")) {
    cat("Restriktor: restricted linear model:\n\n")
  } else if (inherits(x, "summary.conGLM")) {
    cat("Restriktor: restricted generalized linear model:\n\n")
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
  if (se_type %in% c("boot.model.based", "boot.standard") & bootCIs) {
    cat("\nCoefficients from restricted model\nwith",
      100 * level, "pct bootstrap confidence intervals (",bty,"):\n ")  
  } else {
    cat("\nCoefficients:\n")
  } 
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
               na.print = "NA")
  
  if (inherits(x, c("summary.conLM", "summary.conRLM"))) {
    cat("\nResidual standard error:", format(signif(sqrt(x$s2_restr), 
                                                    digits)), "on", rdf, "degrees of freedom")
  } 
  
  #cat("\n")
  if (se_type == "standard") {
    cat("\nStandard errors:", se_type ,"\n")
  } else if (se_type == "const") {
    cat("\nHomoskedastic standard errors.\n")
  } else if (se_type %in% c("boot.model.based", "boot.standard")) {
    if (se_type == "boot.model.based") {
      se_type <- "model-based"
    } else if (se_type == "boot.standard") {
      se_type <- "standard"
    }
    cat("\nBootstrapped standard errors:", se_type ,"\n")
  } else {
    cat("\nHeteroskedastic robust standard errors:", se_type ,"\n")
  }
  
  if (!(inherits(x, "summary.conGLM"))) { 
    if ((x$R2_org - x$R2_reduced) < 1e-08) {
     cat("Multiple R-squared remains", sprintf("%5.3f", x$R2_org),"\n")
    } else {
     cat("Multiple R-squared reduced from", sprintf("%5.3f", x$R2_org), "to", 
         sprintf("%5.3f", x$R2_reduced),"\n")  
    }
  }
  
  if (inherits(x, "summary.conGLM")) {
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
      format(x$dispersion_restr), ")\n\n", apply(cbind(paste(format(c("Null", 
                                                                "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("deviance_null", 
                                                                                                                                 "deviance")]), digits = max(5L, digits + 1L)), " on", 
                                                 format(unlist(x[c("df.residual_null", "rdf")])), " degrees of freedom\n"), 
                                           1L, paste, collapse = " "), sep = "")
  }
  
  goric <- x$goric
  if (!is.null(goric)) {
    ll <- attr(goric, "loglik") 
    PT <- attr(goric, "penalty")
    result_goric <- c(ll, PT, goric)
    names(result_goric) <- c("Loglik", "Penalty", "Goric")
    cat("\nGeneralized Order-Restricted Information Criterion:\n")
    print(result_goric, digits = digits)
  }  
  cat("\n")
  invisible(x)
}
