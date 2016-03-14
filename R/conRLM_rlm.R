#compute constraint estimates
conRLM.rlm <- function(model, constraints, debug = FALSE,
                       se = "default", control = NULL,
                       bvec = NULL, meq = 0L,
                       tol = sqrt(.Machine$double.eps), ...) { 
  
  if (!("rlm" %in% class(model))) {
    stop("ERROR: model must be of class rlm.")
  }
  
  # If "default", conventional standard errors are computed based on inverting
  # the expected augmented information matrix.
  if(!(se %in% c("default","standard","boot.residual","boot.model.based","boot.standard")))
    stop("ERROR: se must be \"standard\", \"boot.model.based\" or \"boot.standard\"")
  if (se == "default") {
    se <- "standard"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  
  Amat <- Amatw
  bvec <- bvecw
  meq <- meqw
  
  if(is.null(bvec)) { bvec <- rep(0L, nrow(Amat)) }
  
  # acknowledgement: code taken from ic.infer package (Ulrike Groemping)
  if (nrow(Amat) - meq - 2 > 2) {
    if (!is.numeric(try(matrix(0, floor((nrow(Amat) - meq -
                                         2)/2), choose(nrow(Amat) - meq, floor((nrow(Amat) - meq -
                                                                                2)/2))), silent = TRUE)))
      stop(paste("test does not work, too many inequality constraints, \n",
                 "interim matrix with ", floor((nrow(Amat) - meq)/2) *
                   choose(nrow(Amat) - meq, floor((nrow(Amat) - meq)/2)),
                 " elements cannot be created", sep = ""))
  }
  
  if(is.null(object$model)) {
    Y <- object$residuals + object$fitted
  } else {
    Data <- object$model
    mfit <- object$model
    Y <- model.response(mfit)  
  }  
  X  <- model.matrix(object)[,,drop = FALSE]
  n <- dim(X)[1]
  beta.unconstr <- coef(object)
  # number of parameters
  p <- length(beta.unconstr)
  df.residual <- nrow(X)-ncol(X)
  
  if(ncol(Amat) != length(beta.unconstr)) {
    stop("length coefficients and ncol(Amat) must be identical")
  }
    
  if (all(Amat %*% c(beta.unconstr) - bvec >= 0 * bvec) & meq == 0) {
    #R2 and adjusted R2, code taken from lmrob() function.
    resid <- object$residuals
    pred <- object$fitted.values
    resp <- pred + resid 
    wgt <- object$w
    if(is.null(object$model)) {
      df.int <- if (any(colnames(X) == "(Intercept)")) 1L else 0L
    } else {
        df.int <- if (attr(object$terms, "intercept")) 1L else 0L
    } 
    resp.mean <- if (df.int == 1L) sum(wgt * resp)/sum(wgt) else 0
    yMy <- sum(wgt * (resp - resp.mean)^2)
    rMr <- sum(wgt * resid^2)
    # bi-square correction
    correc <- 1.207617 
    R2 <- r2correc <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    R2.adjusted <- 1 - (1 - r2correc) * ((n - df.int) / df.residual)
    
    OUT <- list(object.rlm = object,  
                beta.constr = beta.unconstr, 
                residuals = resid(object), 
                fitted.values = fitted(object), 
                weights = object$weights, R2 = R2, R2.adjusted = R2.adjusted,
                df.residual = df.residual, scale = object$s, 
                VCOV = NULL, 
                Amat = Amat, bvec = bvec, meq = meq, iact = 0L, 
                converged = object$converged)
  }
  else {
    call.rlm <- as.list(object$call)
    call.rlm <- call.rlm[-1]
    
    if (call.rlm[["method"]] == "M") {
      if (is.null(call.rlm[["psi"]])) {
        stop("test only applicable with \"psi=psi.bisquare\".")
      } 
      else if (!(call.rlm[["psi"]] == "psi.bisquare")) {
        stop("test only applicable with \"psi=psi.bisquare\".")  
      }
    }
    #fit inequality constrained model
    if (is.null(call.rlm[["formula"]])) {
      call.rlm[["data"]] <- NULL
      call.rlm[["x"]] <- NULL
      call.rlm[["y"]] <- NULL
      call.my <- list(x = X, y = Y, Amat = Amat, meq = meq, bvec = bvec)      
      CALL <- c(call.rlm, call.my)
      rfit <- do.call("icrlm.fit", CALL)
    }
    else {
      call.my <- list(Amat = Amat, meq = meq, bvec = bvec)  
      call.rlm[["data"]] <- as.name("Data")
      CALL <- c(call.rlm, call.my)
      rfit <- do.call("icrlm.formula", CALL)
    }
    
    beta.constr <- coef(rfit)
      beta.constr[abs(beta.constr) < sqrt(.Machine$double.eps)] <- 0L
    
    # vector with the indices of the active constraints
    iact <- rfit$iact
    
    #scale parameter
    sigma <- object$s
    #standard deviation
    tau.hat <- MASS:::summary.rlm(object)$stddev  
    
    #R2 and adjusted R2, code taken from lmrob() function.
    resid <- rfit$residuals
    pred <- rfit$fitted.values
    resp <- pred + resid 
    wgt <- rfit$w
    if (is.null(object$model)) {
      df.int <- if (any(colnames(X) == "(Intercept)")) 1L else 0L
    } else {
      df.int <- if (attr(object$terms, "intercept")) 1L else 0L
    }
    resp.mean <- if (df.int == 1L) sum(wgt * resp)/sum(wgt) else 0
    yMy <- sum(wgt * (resp - resp.mean)^2)
    rMr <- sum(wgt * resid^2)
    # bi-square correction
    correc <- 1.207617 
    R2 <- r2correc <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    R2.adjusted <- 1 - (1 - r2correc) * ((n - df.int) / df.residual)
  
    OUT <- list(object.rlm = object,  
                beta.constr = beta.constr, 
                residuals = resid, 
                fitted.values = pred, 
                weights = rfit$weights, R2 = R2, R2.adjusted = R2.adjusted,
                df.residual = df.residual, scale = rfit$s, 
                VCOV = NULL, 
                Amat = Amat, bvec = bvec, meq = meq, iact = iact, 
                converged = object$converged)
  }
  class(OUT) <- c("conRLM")
    
  OUT
}  
  
