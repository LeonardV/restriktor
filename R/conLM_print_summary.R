print.summary.conLM <- function(object, digits = max(3, getOption("digits") - 2),
                          signif.stars = getOption("show.signif.stars"), ...) {
  
  x <- object
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  resid <- x$residuals
  rdf <- x$rdf
  se.type <- x$se.type[1]
  bootCIs <- attr(x$se.type, "bootCIs")
  bty <- attr(x$se.type, "bty")
  level <- attr(x$se.type, "level")
  
  if (class(x)[1] == "summary.conLM") {
    cat("Restriktor: constrained linear model:\n\n")
  } else if (class(x)[1] == "summary.conRLM") {
    cat("Restriktor: constrained robust linear model:\n\n")
  }
  
  cat(if (!is.null(x$weights) && diff(range(x$weights))) 
    "Weighted ", "Residuals:\n\n", sep = "")
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
    cat("\nCoefficients from constrained model\nwith",
      100 * level, "pct bootstrap confidence intervals (",bty,"):\n ")  
  } else {
    cat("\nCoefficients:\n")
  } 
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
               na.print = "NA")
  
  cat("\n")  
  if (se.type == "const") {
   cat("Homoskedastic standard errors.\n")
  } else if (se.type %in% c("boot.model.based", "boot.standard")) {
    cat("Bootstrapped standard errors:", se.type ,"\n")
  } else {
    cat("Heteroskedastic robust standard errors:", se.type ,"\n")
  }
   
  if (x$R2.org == x$R2.reduced) {
   cat("Constrained model: R2 remains", round(x$R2.org, 4),"\n")
  } else {
   cat("Constrained model: R2 reduced from", round(x$R2.org, 4), "to", round(x$R2.reduced, 4),"\n")  
  }

  goric <- x$goric
  ll <- attr(goric, "loglik") 
  PT <- attr(goric, "penalty")
  result_goric <- c(ll, PT, goric)
  names(result_goric) <- c("Loglik", "Penalty", "Goric")
  cat("\nGeneralized Order-Restricted Information Criterion:\n")
  print(result_goric, digits = digits)
  
  cat("\n")
  invisible(x)
}
