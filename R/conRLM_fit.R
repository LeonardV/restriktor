# Fit an inequality contrained robust linear model.
# If method = "M", the original coef and scale are used as starting value.
# For method = "MM", this can only be when there are no equality constraints.
conRLM_fit <- function(model, maxit = 5000,
                       acc = 1e-14, lqs.control= NULL, 
                       Amat = NULL, bvec = NULL, meq = 0L, ...) {
    
 # acknowledgement: the irls.delta(); irls.rrxwr(); wmad() functions are taken 
 # from the rlm.default function from the MASS package.

 mf <- match.call(expand.dots = FALSE)
 mf$method <- mf$wt.method <- mf$model <- mf$x.ret <- mf$y.ret <- mf$contrasts <- mf$... <- NULL
  
# if (method == "model.frame") {
#   return(mf)
# }
  
  irls.delta <- function(old, new) {
    sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  }

  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(r * w, 1L, length(r)) %*% x)/
              sqrt(matrix(w, 1L, length(r)) %*% (x^2))))/sqrt(sum(w * r^2))
  }
  wmad <- function(x, w) {
    o <- sort.list(abs(x)); x <- abs(x)[o]; w <- w[o]
    p <- cumsum(w)/sum(w)
    n <- sum(p < 0.5)
    if (p[n + 1L] > 0.5) x[n + 1L]/0.6745 else (x[n + 1L] + x[n + 2L])/(2*0.6745)
  }
  
  # response variable
  y <- as.matrix(model$model[, attr(model$terms, "response")])
  # model matrix
  x <- model.matrix(model)[,,drop = FALSE]
  # scale
  scale <- model$s
  # weights
  weights <- model$weights
  
  # original model function call
  call <- as.list(model$call)
  # which method M or MM estimation
  method <- call[["method"]]
  if (is.null(method)) { method <- "M" } 
  # psi function (only bisquare)
  psi <- call[["psi"]]
  if (is.null(psi)) { psi <- "psi.huber" }
  psi <- get(as.character(psi))
  # weights
  w <- model$weights
  if (any(w != 1L)) { stop("weights are not implemented (yet).") }
  # weights method, inverse of the variance or case 
  wt.method <- call[["wt.method"]]
  if (is.null(wt.method)) { wt.method <- "inv.var" } 
  # down weights 
  w <- call[["w"]]
  if (is.null(w)) { w <- rep(1, nrow(x)) }
  # scale estimator - depends on method
  scale.est <- call[["scale.est"]]
  if (is.null(scale.est)) { scale.est <- "MAD" }
  # tuning constant used for Huber proposal 2 scale estimation.
  k2 <- call[["k2"]]
  if (is.null(k2)) { k2 <- 1.345 }
  # the stopping criterion is based on changes in this vector
  test.vec <- call[["test.vec"]]
  if (is.null(test.vec)) { test.vec <- "resid" }
  
  # handling weights
  if (!missing(weights)) {
    if (wt.method == "inv.var") {
      fac <- sqrt(weights)
      y <- y * fac
      x <- x * fac
      wt <- NULL
    } else {
      w <- w * weights
      wt <- weights
    }
  } else { 
    wt <- NULL
  }

  # M-estimation
  if (method == "M") {
    coef <- model$coefficient
    resid <- model$residuals
  # MM-estimation  
  } else if (method == "MM") {
    scale.est <- "MM"
    # Which columns of X are constraint to zero by equality constraints on the paramters.
    # These columns must be removed before computing the robust scale.
    # Can we check this without QP?
    if (meq > 0L) {
      Dmat <- diag(ncol(Amat))
      dvec <- cbind(rep(1, ncol(Amat)))
      QP <- solve.QP(Dmat, dvec, t(Amat[1:meq,,drop = FALSE]), bvec[1:meq], 
                     meq = nrow(Amat[1:meq,,drop = FALSE]))$solution
      x.idx <- QP %in% 0
      temp <- do.call("lqs",
                      c(list(x = x[,!x.idx, drop = FALSE], y, intercept = FALSE, 
                             method = "S", k0 = 1.54764), lqs.control)) 
      coef <- temp$coefficients
      resid <- temp$residuals
      #resid0 <- resid
      scale <- temp$scale
      psi <- psi.bisquare
    } else {
      coef <- model$coefficients
      resid <- model$residuals
      #resid0 <- resid
      scale <- model$s
      psi <- psi.bisquare
    }
  } else {
    stop("'method' is unknown")
  }  
  
  done <- FALSE
  conv <- NULL
  # acknowledgement: code taken from the rlm.default() function from the 
  # MASS package.
  n1 <- ( if (is.null(wt)) { nrow(x) } else { sum(wt) } ) - ncol(x) 
  theta <- 2 * pnorm(k2) - 1
  gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
  
#  if (scale.est != "MM") {
#    scale <- if (is.null(wt)) mad(resid, 0) else wmad(resid, wt)
#  }
  
  for(iiter in 1L:maxit) {
    if (!is.null(test.vec)) { 
      testpv <- get(test.vec)
    }
    if (scale.est != "MM") {
      # acknowledgement: code taken from the rlm.default() function from the 
      # MASS package.
      scale <- 
        if (scale.est == "MAD") {
          if (is.null(wt)) { median(abs(resid)) / 0.6745 } else { wmad(resid, wt) }
        } else if (is.null(wt)) {
          sqrt(sum(pmin(resid^2, (k2 * scale)^2)) / (n1*gamma))
        } else { 
          sqrt(sum(wt*pmin(resid^2, (k2 * scale)^2)) / (n1*gamma))
        }
      
      if (scale == 0) {
        done <- TRUE
        break
      }
    }
    w <- psi(resid / scale)
    if (!is.null(wt)) {
      w <- w * weights
    }
    W <- diag(c(w))
    XX <- t(x) %*% W %*% x
    Xy <- t(x) %*% W %*% y
    QP <- solve.QP(Dmat = XX, dvec
                   = Xy, Amat = t(Amat), 
                   bvec = bvec, meq = meq)
    
    coef <- QP$solution
    resid <- drop(y - x %*% coef)
    iact <- QP$iact

    if (!is.null(test.vec)) {
      convi <- irls.delta(testpv, get(test.vec))
    } else {
      convi <- irls.rrxwr(x, w, resid)
    }
    conv <- c(conv, convi)
    done <- (convi <= acc)
    if (done) {
      break
    }
  }
  
  if (!done) {
    warning(gettextf("'conRLM' failed to converge in %d steps", maxit),
            domain = NA)
  }
  fitted <- drop(x %*% coef)
  names(coef) <- colnames(x)

  cl <- match.call()
  cl[[1L]] <- as.name("conRLM, rlm")
  
  fit <- list(coefficients = coef, 
              residuals = c(y - fitted), 
              wresid = resid,
              #resid0 = if (method == "MM") { resid0 },
              fitted.values = fitted,
              df.residual = NA, w = w,
              scale = scale, psi = psi, k2 = k2,
              weights = if (!missing(weights)) weights,
              conv = conv, converged = done, iter = iiter, x = x, call = cl, 
              Amat = Amat, bvec = bvec, meq = meq, iact = iact)
  class(fit) <- c("conRLM", "rlm", "conLM", "lm")
  
  OUT <- fit
  
  OUT
}
