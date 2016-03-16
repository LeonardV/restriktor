#adjusted functions from the sandwich package.
sandwich <- function(x, bread. = bread, meat. = meat, ...) {
  if(is.list(x$model.org) && !is.null(x$model.org$na.action)) class(x$model.org$na.action) <- "omit"
  if(is.function(bread.)) bread. <- bread.(x)
  if(is.function(meat.)) meat. <- meat.(x, ...)
  n <- NROW(estfun(x))
  return(1/n * (bread. %*% meat. %*% bread.))
}


bread <- function(x, ...) {
  UseMethod("bread")
}

estfun <- function(x, ...) {
  UseMethod("estfun")
}

bread.lm <- function(x, ...)
{
  if(!is.null(x$na.action)) class(x$na.action) <- "omit"
#  X <- model.matrix(x$model.org)[,,drop=FALSE]
#  s2 <- x$s2
  cov <- 1/x$s2 * x$information #solve(t(X) %*% X)
  return(cov * as.vector(sum(summary(x$model.org)$df[1:2])))
}


bread.rlm <- function(x, ...) {
  if(!is.null(x$model.org$na.action)) class(x$na.action) <- "omit"
  xmat <- model.matrix(x$model.org)
  xmat <- naresid(x$model.org$na.action, xmat)
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  res <- residuals(x)
  psi_deriv <- function(z) tukeyChi(z, deriv = 1)
  rval <- sqrt(abs(as.vector(psi_deriv(res/x$model.org$s)/x$model.org$s))) * wts * xmat    
  rval <- chol2inv(qr.R(qr(rval))) * nrow(xmat)
  return(rval)
}


estfun.lm <- function(x, ...) {
  xmat <- model.matrix(x$model.org)
  xmat <- naresid(x$model.org$na.action, xmat)
  if(any(alias <- is.na(coef(x)))) xmat <- xmat[, !alias, drop = FALSE]
  wts <- weights(x$model.org)
  if(is.null(wts)) wts <- 1
  res <- residuals(x)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(zoo:::is.zoo(res)) rval <- zoo:::zoo(rval, index(res), attr(res, "frequency"))
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  return(rval)
}


estfun.rlm <- function(x, ...) {
  xmat <- model.matrix(x$model.org)
  xmat <- naresid(x$model.org$na.action, xmat)
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  res <- residuals(x)
  psi <- function(z) tukeyChi(z) * z
  rval <- as.vector(psi(res/x$model.org$s)) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  if(zoo:::is.zoo(res)) rval <- zoo:::zoo(rval, index(res), attr(res, "frequency"))
  return(rval)
}


meatHC <- function(x, 
                   type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
                   omega = NULL) {
  ## ensure that NAs are omitted
  if(is.list(x$model.org) && !is.null(x$model.org$na.action)) class(x$model.org$na.action) <- "omit"
  
  ## extract X
  X <- model.matrix(x$model.org)
  if(any(alias <- is.na(coef(x)))) X <- X[, !alias, drop = FALSE]
  attr(X, "assign") <- NULL
  n <- NROW(X)
  
  ## get hat values and residual degrees of freedom
  #diaghat <- try(hatvalues(x$model.org), silent = TRUE)
  diaghat <- diag(X%*%(1/x$s2 * x$information)%*%t(X))
  df <- n - NCOL(X)
  
  ## the following might work, but "intercept" is also claimed for "coxph"
  ## res <- if(attr(terms(x), "intercept") > 0) estfun(x)[,1] else rowMeans(estfun(x)/X, na.rm = TRUE)
  ## hence better rely on
  ef <- estfun(x)
  res <- rowMeans(ef/X, na.rm = TRUE)
  ## handle rows with just zeros
  res[apply(abs(ef) < .Machine$double.eps, 1L, all)] <- 0
  
  ## if omega not specified, set up using type
  if(is.null(omega)) {
    type <- match.arg(type)
    if(type == "HC") type <- "HC0"
    switch(type,
           "const" = { omega <- function(residuals, diaghat, df) rep(1, length(residuals)) * sum(residuals^2)/df },
           "HC0"   = { omega <- function(residuals, diaghat, df) residuals^2 },
           "HC1"   = { omega <- function(residuals, diaghat, df) residuals^2 * length(residuals)/df },
           "HC2"   = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat) },
           "HC3"   = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat)^2 },
           "HC4"   = { omega <- function(residuals, diaghat, df) {
             n <- length(residuals)
             p <- as.integer(round(sum(diaghat),  digits = 0))
             delta <- pmin(4, n * diaghat/p)
             residuals^2 / (1 - diaghat)^delta
           }},
           "HC4m"  = { omega <- function(residuals, diaghat, df) {
             gamma <- c(1.0, 1.5) ## as recommended by Cribari-Neto & Da Silva
             n <- length(residuals)
             p <- as.integer(round(sum(diaghat),  digits = 0))
             delta <- pmin(gamma[1], n * diaghat/p) + pmin(gamma[2], n * diaghat/p)
             residuals^2 / (1 - diaghat)^delta
           }},
           "HC5"   = { omega <- function(residuals, diaghat, df) {
             k <- 0.7 ## as recommended by Cribari-Neto et al.
             n <- length(residuals)
             p <- as.integer(round(sum(diaghat),  digits = 0))
             delta <- pmin(n * diaghat/p, pmax(4, n * k * max(diaghat)/p))
             residuals^2 / sqrt((1 - diaghat)^delta)
           }}
    )
  }
  
  ## process omega
  if(is.function(omega)) omega <- omega(res, diaghat, df)
  rval <- sqrt(omega) * X
  rval <- crossprod(rval)/n
  
  return(rval)
}
