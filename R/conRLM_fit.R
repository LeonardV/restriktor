#adjusted code from the rlm() function (MASS package).
conRLM_fit <- function(x, y, weights, w = rep(1, nrow(x)),
                      init = "ls", 
                      scale.est = c("MAD", "Huber", "proposal 2"), k2 = 1.345,
                      method = c("M", "MM"), wt.method = c("inv.var", "case"),
                      maxit = 500, acc = 1e-09, test.vec = "resid", lqs.control= NULL,
                      Amat = NULL, bvec = NULL, meq = 0L, ...)
  {
    irls.delta <- function(old, new)
      sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
    irls.rrxwr <- function(x, w, r)
    {
      w <- sqrt(w)
      max(abs((matrix(r * w, 1L, length(r)) %*% x)/
                sqrt(matrix(w, 1L, length(r)) %*% (x^2))))/sqrt(sum(w * r^2))
    }
    wmad <- function(x, w)
    {
      o <- sort.list(abs(x)); x <- abs(x)[o]; w <- w[o]
      p <- cumsum(w)/sum(w)
      n <- sum(p < 0.5)
      if (p[n + 1L] > 0.5) x[n + 1L]/0.6745 else (x[n + 1L] + x[n + 2L])/(2*0.6745)
    }
    
    method <- match.arg(method)
    wt.method <- match.arg(wt.method)
    nmx <- deparse(substitute(x))
    if(is.null(dim(x))) {
      x <- as.matrix(x)
      colnames(x) <- nmx
    } else x <- as.matrix(x)
    if(is.null(colnames(x)))
      colnames(x) <- paste("X", seq(ncol(x)), sep="")
    if(qr(x)$rank < ncol(x))
      stop("'x' is singular: singular fits are not implemented in 'rlm'")
    
    if(!(any(test.vec == c("resid", "coef", "w", "NULL"))
         || is.null(test.vec))) stop("invalid 'test.vec'")
    ## deal with weights
    xx <- x
    yy <- y
    if(!missing(weights)) {
      if(length(weights) != nrow(x))
        stop("length of 'weights' must equal number of observations")
      if(any(weights < 0)) stop("negative 'weights' value")
      if(wt.method == "inv.var") {
        fac <- sqrt(weights)
        y <- y*fac; x <- x* fac
        wt <- NULL
      } else {
        w <- w * weights
        wt <- weights
      }
    } else wt <- NULL
    
    if(method == "M") {
      scale.est <- match.arg(scale.est)
      psi <- psi.bisquare
      if(!is.function(psi)) psi <- get(psi, mode="function")
      ## match any ... args to those of psi.
#      arguments <- list(...)
#      if(length(arguments)) {
#        pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0L)
#        if(any(pm == 0L)) warning("some of ... do not match")
#        pm <- names(arguments)[pm> 0L]
#        formals(psi)[pm] <- unlist(arguments[pm])
#      }
      if(is.character(init)) {
        if(init == "ls") {
            temp <- lm.wfit(x, y, w, method="qr") 
        }
        else if(init == "lts") {                                                
          if(is.null(lqs.control)) lqs.control <- list(nsamp=200L)
          temp <- do.call("lqs", c(list(x, y, intercept = FALSE), lqs.control))
        } else stop("'init' method is unknown")
          coef <- temp$coefficient
          resid <- temp$residuals
      } else {
        if(is.list(init)) coef <- init$coef
        else coef <- init
        resid <- drop(y - x %*% coef)
      }
    } else if(method == "MM") {
      scale.est <- "MM"
        temp <- do.call("lqs",
                        c(list(x, y, intercept = FALSE, method = "S",
                               k0 = 1.54764), lqs.control))
        coef <- temp$coefficients
        resid <- temp$residuals
        resid0 <- resid
        scale <- temp$scale
      psi <- psi.bisquare
  #    if(length(arguments <- list(...)))
  #      if(match("c", names(arguments), nomatch = 0L)) {
  #        c0 <- arguments$c
  #        if (c0 > 1.54764) formals(psi)$c <- c0
  #        else
  #          warning("'c' must be at least 1.548 and has been ignored")
  #      }
    } else {
      stop("'method' is unknown")
    }
    done <- FALSE
    conv <- NULL
    n1 <- (if(is.null(wt)) nrow(x) else sum(wt)) - ncol(x) 
    theta <- 2*pnorm(k2) - 1
    gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
    ## At this point the residuals are weighted for inv.var and
    ## unweighted for case weights.  Only Huber handles case weights
    ## correctly.
    if(scale.est != "MM")
      scale <- if(is.null(wt)) mad(resid, 0) else wmad(resid, wt)

    for(iiter in 1L:maxit) {
      if(!is.null(test.vec)) testpv <- get(test.vec)
      if(scale.est != "MM") {
        scale <- if(scale.est == "MAD")
          if(is.null(wt)) median(abs(resid))/0.6745 else wmad(resid, wt)
        else if(is.null(wt))
          sqrt(sum(pmin(resid^2, (k2 * scale)^2))/(n1*gamma))
        else sqrt(sum(wt*pmin(resid^2, (k2 * scale)^2))/(n1*gamma))
        if(scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid/scale)
      if(!is.null(wt)) w <- w * weights
      if(is.null(Amat)) {
        temp <- lm.wfit(x, y, w, method="qr")
        coef <- temp$coefficients
        resid <- temp$residuals
        iact <- 0L
      }
      else {
        w <- diag(c(w))
        XX <- t(x) %*% w %*% x
        Xy <- t(x) %*% w %*% y
        out <- con_my_solve_QP(Dmat = XX, dvec = Xy, Amat = t(Amat), bvec = bvec, 
                               meq = meq)
        
        if(out$status == -1 || out$status == -2) {
          done <- FALSE
          break
        }
        coef <- out$solution 
        resid <- drop(y - x %*% coef)
        iact <- out$iact
      }
      
      
      if(!is.null(test.vec)) convi <- irls.delta(testpv, get(test.vec))
      else convi <- irls.rrxwr(x, w, resid)                                      
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if(done) break
    }

    if(!done)
      warning(gettextf("'rlm' failed to converge in %d steps", maxit),
              domain = NA)
    fitted <- drop(xx %*% coef)
      names(coef) <- colnames(x)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1L]] <- as.name("rlm")
    fit <- list(coefficients = coef, residuals = yy - fitted, wresid = resid,
                resid0 = if(method=="MM") { resid0 },
                effects = temp$effects,
                rank = temp$rank, fitted.values = fitted,
                assign = temp$assign, qr = temp$qr, 
                df.residual = NA, w = w,
                s = scale, psi = psi, k2 = k2,
                weights = if(!missing(weights)) weights,
                conv = conv, converged = done, x = xx, call = cl, 
                iact = iact)
    class(fit) <- c("conRLM", "rlm", "lm")
    fit
  }
