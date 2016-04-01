#compute constraint estimates
conRLM.rlm <- function(model, constraints, debug = FALSE,
                       se = "default", control = NULL,
                       bvec = NULL, meq = 0L,
                       tol = sqrt(.Machine$double.eps), ...) { 
  
  Amat <- Amatw; bvec <- bvecw; meq <- meqw
  
  # check for psi.bisquare loss function
  call.model <- as.list(model$call)
  if (!is.null(call.model$method)) {
    if (call.model$method == "MM") {
      if (is.null(call.model$psi)) {
        call.model$psi <- "psi.bisquare"  
      }
    } else if (call.model$method == "M") {
      if (is.null(call.model$psi)) {
        call.model$psi <- "psi.huber"  
      }
    }
  }
  if (is.null(call.model$psi)) {
    stop("Restriktor ERROR: psi must be a bisquare function.")
  } else if (!is.null(call.model$psi)) {
      if (!(call.model$psi == "psi.bisquare")) {
        stop("Restriktor ERROR: psi must be a bisquare function.")
      }
  }
  
  # check class
  if (!("rlm" %in% class(model))) {
    stop("Restriktor ERROR: model must be of class rlm.")
  }
  if (se == "default" | se == "standard") {
    se <- "const"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (missing(constraints) && is.null(bvec)) { 
    bvec <- rep(0L, nrow(Amat)) 
  }
  
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
  
  
  if (is.null(model$model)) {
    Y <- model$residuals + model$fitted
  } else {
    Y <- as.matrix(model$model[, attr(model$terms, "response")])
  }  
  X  <- model.matrix(model)[,,drop = FALSE]
  b.unconstr <- coef(model)  
  so <- MASS:::summary.rlm(model)
  df.residual <- so$df[2]

  if (ncol(Amat) != length(b.unconstr)) {
    stop("length coefficients and ncol(Amat) must be identical")
  }
  
  #R2 and adjusted R2, code taken from lmrob() function.
  resid <- model$residuals
  pred <- model$fitted.values
  resp <- pred + resid 
  wgt <- model$w
  if (is.null(model$model)) {
    df.int <- if (any(colnames(X) == "(Intercept)")) 1L else 0L
  } else {
    df.int <- if (attr(model$terms, "intercept")) 1L else 0L
  }
  resp.mean <- if (df.int == 1L) sum(wgt * resp)/sum(wgt) else 0
  yMy <- sum(wgt * (resp - resp.mean)^2)
  rMr <- sum(wgt * resid^2)
  # bi-square correction
  correc <- 1.207617 
  R2.org <- r2correc <- (yMy - rMr) / (yMy + rMr * (correc - 1))
  #R2.adjusted <- 1 - (1 - r2correc) * ((n - df.int) / df.residual)

  temp <- Sestimator(x = X, y = Y, lqs.control= NULL)
  resid0 <- temp$resid  
  scale <- temp$scale
  
  ll.out <- con_loglik_lm(X = X, y = Y, b = b.unconstr, detU = 1)
  ll <- ll.out$loglik
  s2.ml <- ll.out$Sigma
  tau.hat <- so$stddev  
  s2 <- tau.hat^2
  
  if (all(Amat %*% c(b.unconstr) - bvec >= 0 * bvec) & meq == 0) {
    
    b.constr <- b.unconstr
        
    OUT <- list(CON = NULL,
                partable = NULL,
                constraints = constraints,
                b.unconstr = b.unconstr,
                b.constr = b.unconstr,
                residuals = resid(model),
                init.resid = resid0,
                fitted.values = fitted(model),
                weights = weights(model),
                R2.org = R2.org,
                R2.reduced = R2.org,
                df.residual = df.residual,
                scale = scale, 
                s2 = s2, s2.ml = s2.ml, 
                s2.unc = s2, s2.unc.ml = s2.ml,
                loglik = ll, 
                Sigma = vcov(model),                                            #probably not so robust!
                Amat = Amat, bvec = bvec, meq = meq, iact = 0L,
                converged = model$converged, iter = NULL,
                bootout = NULL)
  } else {
    call.rlm <- as.list(model$call)
    call.rlm <- call.rlm[-1]
    call.rlm[["weights"]] <- weights(model) 
    
    #fit inequality constrained robust model
#    if (is.null(call.rlm[["formula"]])) {
      call.rlm[["data"]] <- NULL
      call.rlm[["x"]] <- NULL
      call.rlm[["y"]] <- NULL
      call.my <- list(x = X, y = Y, Amat = Amat, meq = meq, bvec = bvec)      
      CALL <- c(call.rlm, call.my)
      rfit <- do.call("conRLM_fit", CALL)
#    }
#    else {
#      call.my <- list(Amat = Amat, meq = meq, bvec = bvec)  
      #call.rlm[["data"]] <- as.name("DATA")
#      CALL <- c(call.rlm, call.my)
#      rfit <- do.call("conRLM.formula", CALL)
#    }
    
    b.constr <- coef(rfit)
    b.constr[abs(b.constr) < sqrt(.Machine$double.eps)] <- 0L
    resid0 <- rfit$resid0
    resid  <- rfit$residuals
    scale  <- rfit$s
          
    ll.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
    ll <- ll.out$loglik
    cons2.ml <- ll.out$Sigma    
    
    tau.hat <- MASS:::summary.rlm(rfit)$stddev                                  # <FIXME> in case of "x == 0 or 0 == x" p should be adjusted. </FIXME>
    
    cons2 <- tau.hat^2
    iact <- rfit$iact
    
    #R2 and adjusted R2, code taken from lmrob() function.
    pred <- rfit$fitted.values
    resp <- pred + resid 
    wgt <- diag(rfit$w)
    if (is.null(model$model)) {
      df.int <- if (any(colnames(X) == "(Intercept)")) 1L else 0L
    } else {
      df.int <- if (attr(model$terms, "intercept")) 1L else 0L
    }
    resp.mean <- if (df.int == 1L) sum(wgt * resp)/sum(wgt) else 0
    yMy <- sum(wgt * (resp - resp.mean)^2)
    rMr <- sum(wgt * resid^2)
    # bi-square correction
    correc <- 1.207617 
    R2.reduced <- r2correc <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    #R2.adjusted <- 1 - (1 - r2correc) * ((n - df.int) / df.residual)
  
    OUT <- list(CON = NULL,
                partable = NULL,
                constraints = constraints,
                b.unconstr = b.unconstr,
                b.constr = b.constr,
                residuals = resid,
                init.resid = resid0,
                wresid = rfit$wresid,
                fitted.values = pred,
                weights = rfit$weights,
                R2.org = R2.org,
                R2.reduced = R2.reduced,
                df.residual = df.residual,
                scale = rfit$s,                                                               
                s2 = cons2, s2.ml = cons2.ml,
                s2.unc = s2, s2.unc.ml = s2.ml,
                loglik = ll, 
                Sigma = vcov(model),                                             #probably not robust???
                Amat = Amat, bvec = bvec, meq = meq, iact = iact,
                converged = rfit$converged, iter = rfit$iter,
                bootout = NULL)
  }
  
  OUT$model.org <- model
  OUT$CON <- if (is.character(constraints)) { CON }
  OUT$partable <- if (is.character(constraints)) { partable }
  OUT$se <- se 
  
  if (se != "no") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      information <- 1/s2 * crossprod(X)  
      
      #only for MM-estimation? 
      #information <- solve(vcovMM(X, resid0, resid, scale))                     
      information.inv <- con_augmented_information(information = information,
                                                   X = X, b.unconstr = b.unconstr, 
                                                   b.constr = b.constr,
                                                   Amat = Amat, 
                                                   bvec = bvec, meq = meq)
      OUT$information.inverted <- information.inv
    } else if (se == "boot.model.based") { 
      OUT$bootout <- con_boot_lm(model, B = ifelse(is.null(control$B),
                                                   999, control$B), 
                                 fixed = TRUE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "no")
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(model, B = ifelse(is.null(control$B),
                                                   999, control$B),
                                 fixed = FALSE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "no")
    }
  }
  
  class(OUT) <- c("conRLM","conLM","rlm")
  
  return(OUT)
}
