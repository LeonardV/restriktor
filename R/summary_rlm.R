# acknowledgement: code taken from MASS package
# adapted to correct rdf for equality constraints by LV.
summary_rlm <- function (object, method = c("XtX", "XtWX"), 
                         correlation = FALSE, 
                         Amat = NULL, meq = NULL, ...) {
  method <- match.arg(method)
  s <- object$s
  coef <- object$coefficients
  ptotal <- length(coef)
  wresid <- object$wresid
  res <- object$residuals
  n <- length(wresid)
  if (any(na <- is.na(coef))) 
    coef <- coef[!na]
  cnames <- names(coef)
  p <- length(coef)
  rinv <- diag(p)
  dimnames(rinv) <- list(cnames, cnames)
  wts <- if (length(object$weights)) {
    object$weights
  } else {
    rep(1, n)
  }
  if (length(object$call$wt.method) && object$call$wt.method == "case") {
    rdf <- sum(wts) - p
    if (!is.null(Amat)) {
      rdf <- sum(wts) - (p - qr(Amat[0:meq,])$rank)
    } 
    w <- object$psi(wresid/s)
    S <- sum(wts * (wresid * w)^2)/rdf
    psiprime <- object$psi(wresid/s, deriv = 1)
    m1 <- sum(wts * psiprime)
    m2 <- sum(wts * psiprime^2)
    nn <- sum(wts)
    mn <- m1/nn
    kappa <- 1 + p * (m2 - m1^2/nn)/(nn - 1)/(nn * mn^2)
    stddev <- sqrt(S) * (kappa/mn)
  } else {
    res <- res * sqrt(wts)
    rdf <- n - p
    if (!is.null(Amat)) {
      rdf <- n - (p - qr(Amat[0:meq,])$rank)
    }
    w <- object$psi(wresid/s)
    S <- sum((wresid * w)^2)/rdf
    psiprime <- object$psi(wresid/s, deriv = 1)
    mn <- mean(psiprime)
    kappa <- 1 + p * var(psiprime)/(n * mn^2)
    stddev <- sqrt(S) * (kappa/mn)
  }
  X <- if (length(object$weights)) 
    object$x * sqrt(object$weights)
  else object$x
  if (method == "XtWX") {
    mn <- sum(wts * w)/sum(wts)
    X <- X * sqrt(w/mn)
  }
  R <- qr(X)$qr
  R <- R[1L:p, 1L:p, drop = FALSE]
  R[lower.tri(R)] <- 0
  rinv <- solve(R, rinv)
  dimnames(rinv) <- list(cnames, cnames)
  rowlen <- (rinv^2 %*% rep(1, p))^0.5
  names(rowlen) <- cnames
  if (correlation) {
    correl <- rinv * array(1/rowlen, c(p, p))
    correl <- correl %*% t(correl)
  }
  else correl <- NULL
  coef <- array(coef, c(p, 3L))
  dimnames(coef) <- list(cnames, c("Value", "Std. Error", "t value"))
  coef[, 2] <- rowlen %o% stddev
  coef[, 3] <- coef[, 1]/coef[, 2]
  object <- object[c("call", "na.action")]
  object$residuals <- res
  object$coefficients <- coef
  object$sigma <- s
  object$stddev <- stddev
  object$df <- c(p, rdf, ptotal)
  object$r.squared <- NA
  object$cov.unscaled <- rinv %*% t(rinv)
  object$correlation <- correl
  object$terms <- NA
  
  object
}