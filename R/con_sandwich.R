# adjusted functions from the sandwich package (based on sandwich version 3.1-0).
# adapted for restricted estimation by LV.
#
# the bread functions replace the classical unscaled covariance matrix by the
# inverted augmented information matrix (Schoenberg, 1997), so that the
# (in)equality restrictions are taken into account. without restrictions the
# results are identical to sandwich::vcovHC().

sandwich <- function(x, bread. = bread, meat. = meatHC, ...) {
  if (is.list(x) && !is.null(x$na.action)) { class(x$na.action) <- "omit" }
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
  cov.unscaled <- 1/x$s2 * cov
  return(cov.unscaled * as.vector(sum(summary(x$model.org)$df[1:2])))
}


bread.conRLM <- function(x, ...) {
  xmat <- model.matrix(x)
  wts <- weights(x)
  if (is.null(wts)) { wts <- 1 }
  res <- residuals(x)
  # x$scale is re-estimated under the restrictions; it equals x$model.org$s
  # when no restrictions are violated.
  psi_deriv <- function(z) x$model.org$psi(z, deriv = 1)
  rval <- sqrt(abs(as.vector(psi_deriv(res / x$scale) / x$scale))) * wts * xmat
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

  # x$information is scaled by 1/x$dispersion, hence the inverted information
  # is unscaled here and rescaled with the ML-type dispersion, exactly as in
  # sandwich::bread.glm().
  cov <- attr(x$information, "inverted")
  cov.unscaled <- cov / x$dispersion

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
  xmat <- naresid(x$na.action, xmat)
  wts <- weights(x)
  if (is.null(wts)) { wts <- 1 }
  res <- residuals(x)
  psi <- function(z) x$model.org$psi(z) * z
  rval <- as.vector(psi(res / x$scale)) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL

  return(rval)
}


estfun.conGLM <- function(x, ...) {
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
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

  return(rval)
}


meatHC <- function(x,
                   type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5"),
                   omega = NULL) {

  if (is.list(x) && !is.null(x$na.action)) { class(x$na.action) <- "omit" }

  type <- match.arg(type)
  if (type == "HC") { type <- "HC0" }

  X <- model.matrix(x$model.org)
  if (any(alias <- is.na(coef(x)))) X <- X[, !alias, drop = FALSE]
  attr(X, "assign") <- NULL
  n <- NROW(X)

  # correct the residual degrees of freedom for equality constraints
  df <- n - (NCOL(X) - qr(x$constraints[seq_len(x$neq), , drop = FALSE])$rank)

  ### get hat values ###
  # hatvalues() cannot be applied to restricted model objects, hence the
  # diagonal of the hat matrix is computed directly via a QR factorization
  # of the (weighted) model matrix.
  if (type %in% c("HC2", "HC3", "HC4", "HC4m", "HC5")) {
    if (inherits(x, "conRLM")) {
      # weights from the IWLS process (psi(resid/scale))
      rW <- sqrt(x$wgt)
    } else if (inherits(x, "conGLM")) {
      # working weights
      rW <- sqrt(x$weights)
    } else if (!is.null(x$weights)) {
      # prior weights (weighted least squares)
      rW <- sqrt(x$weights)
    } else {
      rW <- rep(1, n)
    }

    Xh <- rW * X
    # under equality constraints the fitted values live in the column space
    # of X %*% N, with N a basis for the null space of the equality rows;
    # the hat values are therefore based on that reduced space. inequality
    # constraints are not counted, consistent with the df correction above.
    if (x$neq > 0L) {
      Aeq <- x$constraints[seq_len(x$neq), , drop = FALSE]
      QRt <- qr(t(Aeq))
      N <- qr.Q(QRt, complete = TRUE)[, -seq_len(QRt$rank), drop = FALSE]
      Xh <- Xh %*% N
    }

    QR <- qr.default(Xh)
    Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))
    # diagonal of the hat matrix; observations with a zero weight simply get
    # a zero hat value.
    diaghat <- rowSums(Q^2)
  }

  ## extract the empirical estimating functions and residuals
  ef <- estfun(x)
  res <- rowMeans(ef / X, na.rm = TRUE)
  ## handle rows with just zeros
  all0 <- apply(abs(ef) < .Machine$double.eps, 1L, all)
  res[all0] <- 0
  ## in case of just zeros the residuals cannot be recovered from the
  ## estimating functions; for the "const" estimator fall back to the model
  ## residuals (e.g., observations fully downweighted by conRLM)
  if (any(all0) && type == "const") {
    if (inherits(x, "conGLM")) {
      res <- as.vector(residuals(x, "working")) * weights(x, "working")
      if (!(substr(x$model.org$family$family, 1L, 17L) %in% c("poisson", "binomial", "Negative Binomial"))) {
        res <- res * sum(weights(x, "working"), na.rm = TRUE) / sum(res^2, na.rm = TRUE)
      }
    } else {
      res <- as.vector(residuals(x))
      if (!is.null(weights(x))) { res <- res * weights(x) }
    }
  }

  ## if omega not specified, set up using type
  if (is.null(omega)) {
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
