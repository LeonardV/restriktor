# adjusted functions from the sandwich package.
# adapted by LV

sandwich <- function(x, bread. = bread, meat. = meatHC, ...) {
  if (is.function(bread.)) { bread. <- bread.(x) }
  if (is.function(meat.)) { meat. <- meat.(x, ...) }
  n <- NROW(estfun(x))
  return(1/n * (bread. %*% meat. %*% bread.))
}


# meat <- function(x, adjust = FALSE, ...)
# {
#   if (is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
#   psi <- estfun(x, ...)
#   k <- NCOL(psi)
#   n <- NROW(psi)
#   rval <- crossprod(as.matrix(psi))/n
#   if (adjust) rval <- n/(n-k) * rval
#   rownames(rval) <- colnames(rval) <- colnames(psi)
#   return(rval)
# }


bread <- function(x, ...) {
  UseMethod("bread")
}

bread.conLM <- function(x, ...) {
  I.inv <- attr(x$information, "inverted")
  cov.unscaled <- 1/x$s2.restr * I.inv
  return(cov.unscaled * as.vector(sum(summary(x$model.org)$df[1:2])))
}


bread.conRLM <- function(x, ...) {
    xmat <- model.matrix(x)
    wts <- weights(x)
    if (is.null(wts)) { wts <- 1 }
    res <- residuals(x)
    psi_deriv <- function(z) x$model.org$psi(z, deriv = 1)
    rval <- sqrt(abs(as.vector(psi_deriv(res / x$model.org$s) / x$model.org$s))) * wts * xmat    
    rval <- chol2inv(qr.R(qr(rval))) * nrow(xmat)
    rval <- solve(rval) #/ x$model.org$s * nrow(xmat)                           #solve(rval)
    
    is.augmented <- TRUE
    if (all(c(x$constraints) == 0)) { is.augmented <- FALSE }
    
    rval <- con_augmented_information(information  = rval, 
                                      is.augmented = is.augmented,
                                      X            = xmat, 
                                      b.unrestr    = x$b.unrestr, 
                                      b.restr      = x$b.restr, 
                                      Amat         = x$constraints, 
                                      bvec         = x$rhs, 
                                      meq          = x$neq)
    return(rval$information)
}



estfun <- function(x, ...) {
  UseMethod("estfun")
}

estfun.conLM <- function(x, ...) {
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if (any(alias <- is.na(coef(x)))) { 
    xmat <- xmat[, !alias, drop = FALSE] 
  }
  wts <- weights(x)
  if (is.null(wts)) { wts <- 1 }
  res <- residuals(x)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
#  if (is.zoo(res)) {
#    rval <- zoo(rval, index(res), attr(res, "frequency"))
#  }
#  if (is.ts(res)) {
#    rval <- ts(rval, start = start(res), frequency = frequency(res))
#  }
  return(rval)
}

estfun.conRLM <- function(x, ...) {
  xmat <- model.matrix(x)
  wts <- weights(x)
  if (is.null(wts)) { wts <- 1 }
  res <- residuals(x)
  psi <- function(z) x$model.org$psi(z) * z
  rval <- as.vector(psi(res/x$model.org$s)) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
#  if (is.ts(res)) {
#    rval <- ts(rval, start = start(res), frequency = frequency(res))
#  }
#  if (is.zoo(res)) {
#    rval <- zoo(rval, index(res), attr(res, "frequency"))
#  }
  
  return(rval)
}


meatHC <- function(x, 
                   type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
                   omega = NULL) {
  
  X <- model.matrix(x$model.org)
  if (any(alias <- is.na(coef(x)))) X <- X[, !alias, drop = FALSE]
  attr(X, "assign") <- NULL
  n <- NROW(X)
  
  ### get hat values and residual degrees of freedom ###
  if (inherits(x, "conRLM")) {
    #W <- sqrt(x$model.org$w)                                     #constrained or unconstrained?
    # sqrt of the weights
    rW <- sqrt(x$w)
  } else if (class(x)[1] == "conLM") {
      if (!is.null(x$weights)) {
        rW <- sqrt(x$weights)
      } else {
        rW <- rep(1, n)
      }
  } 
  
  ### compute hat matrix ###
  ### code added by LV ###
  
  ## matrix form ##
  # diaghat <- diag(X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W)
  
  # QR factorization #
  # rescaled X
  rWX <- rW * X
  # QR factorization
  QR <- qr.default(rWX)
  Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))
  # for weighted least squares
  Q1 <- (1 / rW) * Q
  Q2 <- rW * Q
  ## diagonal
  diaghat <- rowSums(Q1 * Q2)  
  
  
  p <- NCOL(X)
  # correct df for equality constraints
  df <- n - (p - qr(x$constraints[0:x$neq,])$rank)
  ## the following might work, but "intercept" is also claimed for "coxph"
  ## res <- if(attr(terms(x), "intercept") > 0) estfun(x)[,1] else rowMeans(estfun(x)/X, na.rm = TRUE)
  ## hence better rely on
  ef <- estfun(x)
  res <- rowMeans(ef/X, na.rm = TRUE)
  ## handle rows with just zeros
  res[apply(abs(ef) < .Machine$double.eps, 1L, all)] <- 0
  
  ## if omega not specified, set up using type
  if (is.null(omega)) {
    type <- match.arg(type)
    if (type == "HC") type <- "HC0"
    switch (type,
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
             p <- as.integer(round(sum(diaghat), digits = 0))
             delta <- pmin(n * diaghat/p, pmax(4, n * k * max(diaghat)/p))
             residuals^2 / sqrt((1 - diaghat)^delta)
           }}
    )
  }
  
  ## process omega
  if (is.function(omega)) { omega <- omega(res, diaghat, df) }
  rval <- sqrt(omega) * X
  rval <- crossprod(rval)/n
  
  return(rval)
}
