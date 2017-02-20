# adjusted functions from the sandwich package.
# adapted by LV.

sandwich <- function(x, bread. = bread, meat. = meatHC, ...) {
  if (is.function(bread.)) { bread. <- bread.(x) }
  if (is.function(meat.)) { meat. <- meat.(x, ...) }
  n <- NROW(estfun(x))
  return(1/n * (bread. %*% meat. %*% bread.))
}


bread <- function(x, ...) {
  UseMethod("bread")
}

bread.conLM <- function(x, ...) {
  cov <- attr(x$information, "inverted")
  cov.unscaled <- 1/x$s2.restr * cov
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
    rval <- solve(rval) 
    
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


bread.conGLM <- function(x, ...) {
  sx <- summary(x$model.org)
  wres <- as.vector(residuals(x, "working")) * weights(x, "working")
  dispersion.restr <- if (substr(x$model.org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) {
    1
  } else {
    sum(wres^2) / sum(weights(x, "working"))
  }
  
  cov <- attr(x$information, "inverted") #* x$df.residual / as.vector(sum(sx$df[1:2]))
  cov.unscaled <- cov / x$dispersion
  
  #cov.unscaled <- 1 / dispersion.restr * I.inv
  
  return(cov.unscaled * as.vector(sum(sx$df[1:2])) * dispersion.restr)
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

  return(rval)
}


estfun.conGLM <- function(x, ...) {
  xmat <- model.matrix(x)
  if (any(alias <- is.na(coef(x)))) {
    xmat <- xmat[, !alias, drop = FALSE]
  }
  wres <- as.vector(residuals(x, "working")) * weights(x, "working")
  dispersion <- if (substr(x$model.org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) {
    1
  } else {
    sum(wres^2, na.rm = TRUE) / sum(weights(x, "working"), na.rm = TRUE)
  }
  
  rval <- wres * xmat / dispersion
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
#  res <- residuals(x, type = "pearson")
#  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
#  if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
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
  if (type %in% c("HC2", "HC3", "HC4", "HC4m", "HC5")) {
    if (inherits(x, "conRLM")) {
      rW <- sqrt(x$w)
    } else if (inherits(x, "conLM")) {
        if (!is.null(x$weights)) {
          rW <- sqrt(x$weights)
        } else {
          rW <- rep(1, n)
        }
    } else if (inherits(x, "conGLM")) {
      rW <- sqrt(x$weights)
      }
    
    # it may happen that a weight equals 0.
    # to avoid computational issues, we replace the zero with 1e-08
    idx.rW <- which(rW == 0)
    rW[idx.rW] <- 1e-08
    
    if (any(rW == 0)) {
      warning("Restriktor WARNING: weights ", idx.rW, " are exactly zero. This causes problems for", 
              "\ncomputing the hat-matrix. To compute the ", sQuote(type), " standard errors, we added a small",
              "\nnumber 1e-08. Note that the results may not be trustworthy.")
    }
    
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
  }
  
  
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
  rval <- crossprod(rval) / n
  
  return(rval)
}
