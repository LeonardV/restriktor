conRLM_fit <- function (x, y, weights, ..., w = rep(1, nrow(x)), init = "ls", 
                        psi = psi.bisquare, scale.est = c("MAD", "Huber", "proposal 2"), 
                        k2 = 1.345, method = c("M", "MM"), 
                        wt.method = c("inv.var", "case"), 
                        maxit = 20, acc = 1e-04, test.vec = "resid", 
                        lqs.control = NULL, Amat, bvec, meq,
                        tol = sqrt(.Machine$double.eps)) {

  irls.delta <- function(old, new) sqrt(sum((old - new)^2)/max(1e-20, 
                                                               sum(old^2)))
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(r * w, 1L, length(r)) %*% x)/sqrt(matrix(w, 
                                                             1L, length(r)) %*% (x^2))))/sqrt(sum(w * r^2))
  }
  wmad <- function(x, w) {
    o <- sort.list(abs(x))
    x <- abs(x)[o]
    w <- w[o]
    p <- cumsum(w)/sum(w)
    n <- sum(p < 0.5)
    if (p[n + 1L] > 0.5) 
      x[n + 1L]/0.6745
    else (x[n + 1L] + x[n + 2L])/(2 * 0.6745)
  }
  method <- match.arg(method)
  wt.method <- match.arg(wt.method)
  nmx <- deparse(substitute(x))
  if (is.null(dim(x))) {
    x <- as.matrix(x)
    colnames(x) <- nmx
  }
  else x <- as.matrix(x)
  if (is.null(colnames(x))) 
    colnames(x) <- paste("X", seq(ncol(x)), sep = "")
  if (qr(x)$rank < ncol(x)) 
    stop("'x' is singular: singular fits are not implemented in 'rlm'")
  if (!(any(test.vec == c("resid", "coef", "w", "NULL")) || 
        is.null(test.vec))) 
    stop("invalid 'test.vec'")
  xx <- x
  yy <- y
  if (!missing(weights)) {
    if (length(weights) != nrow(x)) 
      stop("length of 'weights' must equal number of observations")
    if (any(weights < 0)) 
      stop("negative 'weights' value")
    if (wt.method == "inv.var") {
      fac <- sqrt(weights)
      y <- y * fac
      x <- x * fac
      wt <- NULL
    }
    else {
      w <- w * weights
      wt <- weights
    }
  }
  else wt <- NULL
  if (method == "M") {
    scale.est <- match.arg(scale.est)
    if (!is.function(psi)) 
      psi <- get(psi, mode = "function")
    arguments <- list(...)
    if (length(arguments)) {
      pm <- pmatch(names(arguments), names(formals(psi)), 
                   nomatch = 0L)
      if (any(pm == 0L)) 
        warning("some of ... do not match")
      pm <- names(arguments)[pm > 0L]
      formals(psi)[pm] <- unlist(arguments[pm])
    }
    if (is.character(init)) {
      temp <- if (init == "ls") 
        lm.wfit(x, y, w, method = "qr")
      else if (init == "lts") {
        if (is.null(lqs.control)) 
          lqs.control <- list(nsamp = 200L)
        do.call("lqs", c(list(x, y, intercept = FALSE), 
                         lqs.control))
      }
      else stop("'init' method is unknown")
      coef <- temp$coefficient
      resid <- temp$residuals
    }
    else {
      if (is.list(init)) 
        coef <- init$coef
      else coef <- init
      resid <- drop(y - x %*% coef)
    }
  }
  else if (method == "MM") {
    scale.est <- "MM"
    psi <- psi.bisquare
    temp <- do.call("lqs", c(list(x, y, intercept = FALSE, 
                                    method = "S", k0 = 1.54764), lqs.control))
      coef <- temp$coefficients
      resid <- temp$residuals
      if (length(arguments <- list(...))) 
        if (match("c", names(arguments), nomatch = 0L)) {
          c0 <- arguments$c
          if (c0 > 1.54764) 
            formals(psi)$c <- c0
          else warning("'c' must be at least 1.548 and has been ignored")
        }
      scale <- temp$scale
  } else {
    stop("'method' is unknown")
  }
    
  done <- FALSE
  conv <- NULL
  n1 <- (if (is.null(wt)) 
    nrow(x)
    else sum(wt)) - ncol(x)
  theta <- 2 * pnorm(k2) - 1
  gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
  if (scale.est != "MM") 
    scale <- if (is.null(wt)) 
      mad(resid, 0)
  else wmad(resid, wt)
  for (iiter in 1L:maxit) {
    if (!is.null(test.vec)) 
      testpv <- get(test.vec)
    if (scale.est != "MM") {
      scale <- if (scale.est == "MAD") 
        if (is.null(wt)) 
          median(abs(resid))/0.6745
      else wmad(resid, wt)
      else if (is.null(wt)) 
        sqrt(sum(pmin(resid^2, (k2 * scale)^2))/(n1 * 
                                                   gamma))
      else sqrt(sum(wt * pmin(resid^2, (k2 * scale)^2))/(n1 * 
                                                           gamma))
      if (scale == 0) {
        done <- TRUE
        break
      }
    }
    w <- psi(resid/scale)
    if (!is.null(wt)) 
      w <- w * weights
    
    #############################################################
    W <- sqrt(w)
    fit <- con_solver_rlm(X = x * W, y = y * W,
                          Amat = Amat, bvec = bvec,
                          meq = meq, maxit = maxit, 
                          absval = tol)
    
    coef <- fit$solution
    resid <- drop(y - x %*% coef)
    iact <- fit$iact
    #############################################################
    
    #temp <- lm.wfit(x, y, w, method = "qr")
    #coef <- temp$coefficients
    #resid <- temp$residuals
    if (!is.null(test.vec)) 
      convi <- irls.delta(testpv, get(test.vec))
    else convi <- irls.rrxwr(x, w, resid)
    conv <- c(conv, convi)
    done <- (convi <= acc)
    if (done) 
      break
  }
  if (!done) 
    warning(gettextf("'rlm' failed to converge in %d steps", 
                     maxit), domain = NA)
  fitted <- drop(xx %*% coef)
  cl <- match.call()
  cl[[1L]] <- as.name("rlm")
  fit <- list(coefficients = coef, residuals = yy - fitted, 
              wresid = resid, effects = temp$effects, rank = temp$rank, 
              fitted.values = fitted, assign = temp$assign, qr = temp$qr, 
              df.residual = NA, w = w, s = scale, psi = psi, k2 = k2, 
              weights = if (!missing(weights)) weights, conv = conv, 
              converged = done, x = xx, call = cl, 
              Amat = Amat, bvec = bvec, meq = meq, iact = iact)
  class(fit) <- c("conRLM", "rlm", "lm")
  fit
}