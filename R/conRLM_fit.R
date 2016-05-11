# Fit an inequality contrained robust linear model.
# If method = "M", the original coef and scale are used as starting value.
# For method = "MM", the original coef and scale can only be used in case of 
# no equality constraints.
conRLM_fit <- function(model,  
                       Amat = NULL, bvec = NULL, meq = 0L, 
                       tol = sqrt(.Machine$double.eps), ...) {
    
 # acknowledgement: the irls.delta(); irls.rrxwr(); wmad() functions are taken 
 # from the rlm.default function from the MASS package.

 mf <- match.call(expand.dots = FALSE)
 mf$method <- mf$wt.method <- mf$model <- mf$x.ret <- mf$y.ret <- mf$contrasts <- mf$... <- NULL
  
 irls.delta <- function(old, new) {
  sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
 }

 irls.rrxwr <- function(x, w, r) {
  w <- sqrt(w)
  max(abs((matrix(r * w, 1L, length(r)) %*% x) / 
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
 if (any(weights != 1L)) { stop("Restriktor ERROR: weights are not implemented (yet).") }
 # original model function call
 call <- as.list(model$call)
 # which method M or MM estimation
 method <- call[["method"]]
 # tukey's tuning constant
 cc <- call[["c"]]
 # weights method, inverse of the variance or case 
 wt.method <- call[["wt.method"]]
 # down weights 
 w <- call[["w"]]
 # scale estimator - depends on method
 scale.est <- call[["scale.est"]]
 # tuning constant used for Huber proposal 2 scale estimation.
 k2 <- call[["k2"]]
 # the stopping criterion is based on changes in this vector
 test.vec <- call[["test.vec"]]
 # an optional list of control values for lqs.
 lqs.control <- call[["lqs.control"]]
 # the limit on the number of IWLS iterations.
 maxit <- call[["maxit"]]
 # the accuracy for the stopping criterion.
 acc <- call[["acc"]]
 
 xx <- x
 yy <- y
 
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

  psi <- psi.bisquare
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
    # Is this fool proof?
    if (meq > 0L) {
      Dmat <- diag(ncol(Amat))
      dvec <- cbind(rep(1, ncol(Amat)))
      QP <- solve.QP(Dmat, dvec, t(Amat[1:meq,,drop = FALSE]), bvec[1:meq], 
                     meq = nrow(Amat[1:meq,,drop = FALSE]))$solution
      QP[abs(QP) < tol] <- 0L
      x.idx <- QP %in% 0
      temp <- do.call("lqs",
                      c(list(x = x[,!x.idx, drop = FALSE], y, intercept = FALSE, 
                             method = "S", k0 = 1.54764), lqs.control)) 
      coef  <- temp$coefficients
      resid <- temp$residuals
      scale <- temp$scale  
    } else {
      coef  <- model$coefficients
      resid <- model$residuals
      scale <- model$s
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
    w <- psi(resid / scale, c = cc)
    if (!is.null(wt)) {
      w <- w * weights
    }
    # constained estimation part.
    W <- diag(c(w))
    XX <- t(x) %*% W %*% x
    Xy <- t(x) %*% W %*% y
    QP <- solve.QP(Dmat = XX, dvec = Xy, Amat = t(Amat), 
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
  fitted <- drop(xx %*% coef)
  names(coef) <- colnames(x)

  cl <- match.call()
  cl[[1L]] <- as.name("conRLM, rlm")
  
  fit <- list(coefficients = coef, 
              residuals = c(yy - fitted), 
              wresid = resid,
              #resid0 = if (method == "MM") { resid0 },
              fitted.values = fitted,
              df.residual = NA, w = w,
              scale = scale, psi = psi, k2 = k2,
              weights = if (!missing(weights)) weights,
              conv = conv, converged = done, iter = iiter, x = x, call = cl, 
              Amat = Amat, bvec = bvec, meq = meq, iact = iact)
  class(fit) <- c("conRLM","rlm")
  
  OUT <- fit
  
  OUT
}
