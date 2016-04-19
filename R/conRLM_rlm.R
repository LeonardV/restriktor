#compute constraint robust estimates
conRLM.rlm <- function(model, constraints, debug = FALSE,
                       se = "default", control = NULL,
                       bvec = NULL, meq = 0L,
                       tol = sqrt(.Machine$double.eps), ...) { 
  
  cl <- match.call()
  Amat <- Amatw; bvec <- bvecw; meq <- meqw
  
  # check class
  if (!("rlm" %in% class(model))) {
    stop("Restriktor ERROR: model must be of class rlm.")
  }
  if (se == "default" | se == "standard") {
    se <- "const"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (!(se %in% c("none","const","boot.model.based","boot.standard","HC","HC0",
                  "HC1","HC2","HC3","HC4","HC4m","HC5"))) {
    stop("standard error method ", sQuote(se), " unknown.")
  }
  if (se == "boot.model.based" & any(Amat[,1] == 1)) { 
    stop("no constraints on intercept possible for model based bootstrap.")
  }
  if (missing(constraints) && is.null(bvec)) { 
    bvec <- rep(0L, nrow(Amat)) 
  }
  
  Y <- as.matrix(model$model[, attr(model$terms, "response")])
  X  <- model.matrix(model)[,,drop = FALSE]
  # unconstrained coefficients
  b.unconstr <- coef(model)  
  
  # original model function call
  call <- as.list(model$call)
  # which method M or MM estimation
  method <- call[["method"]]
  if (is.null(method)) {
    model$call[["method"]] <- "M"
  } else if (method == "MM") {
    call[["psi"]] <- "psi.bisquare"
  } else if (is.name(method)) {
    stop("define ", sQuote(method), " inside function.")
  }
  # psi function (only bisquare, will be checked later!)
  psi <- as.character(call[["psi"]])
  if (is.null(psi)) {
    psi <- model$call[["psi"]] <- "psi.huber"
  } 
  # weights
  weights <- model$weights
  if (any(weights != 1L)) {
    stop("weights are not implemented (yet).")
  }
  # weights method, inverse of the variance or case 
  wt.method <- call[["wt.method"]]
  if (is.null(wt.method)) {
    model$call[["wt.method"]] <- "inv.var"
  } else if (is.name(wt.method)) {
    stop("define ", sQuote(wt.method), " inside function.")
  }
  # down weights 
  w <- call[["w"]]
  if (is.null(w)) {
    model$call[["w"]] <- rep(1, nrow(X))
  } else if (is.name(w)) {
    stop("define ", sQuote(w), " inside function.")
  }
  # scale estimator - depends on method
  scale.est <- call[["scale.est"]]
  if (is.null(scale.est)) {
    model$call[["scale.est"]] <- "MAD"
  } else if (is.name(scale.est)) {
    stop("define ", sQuote(scale.est), " inside function.")
  }
  # tuning constant used for Huber proposal 2 scale estimation.
  k2 <- call[["k2"]]
  if (is.null(k2)) {
    model$call[["k2"]] <- 1.345
  } else if (is.name(k2)) {
    stop("define ", sQuote(k2), " inside function.")
  }
  # the stopping criterion is based on changes in this vector
  test.vec <- call[["test.vec"]]
  if (is.null(test.vec)) {
    model$call[["test.vec"]] <- "resid"
  } else if (is.name(resid)) {
    stop("define ", sQuote(resid), " inside function.")
  }
  
  # Restriktor only supports bisquare psi loss function (for now).
  if (method == "M") {
    if (psi != "psi.bisquare") {
      stop("restriktor only supports the bisquare loss function (for now).")
    }
  }
  
  # acknowledgement: check for too many inequality constraints is taken from 
  # the ic.infer package. 
  if (nrow(Amat) - meq - 2 > 2) {
    if (!is.numeric(try(matrix(0, floor((nrow(Amat) - meq -
                                         2)/2), choose(nrow(Amat) - meq, floor((nrow(Amat) - meq -
                                                                                2)/2))), silent = TRUE)))
      stop(paste("test does not work, too many inequality constraints, \n",
                 "interim matrix with ", floor((nrow(Amat) - meq)/2) *
                   choose(nrow(Amat) - meq, floor((nrow(Amat) - meq)/2)),
                 " elements cannot be created", sep = ""))
  }
  
  # Adjusted summary function from MASS:::summary.rlm().
  # The df needs to be corrected for equality constraints and
  # we need the df for ML. 
  so <- summary_rlm(model, Amat = Amat, meq = meq)
  # unconstrained scale estimate for the standard errors
  tau <- so$stddev
  # # unconstrained scale estimate for the variance
  s2.con <- s2.unc <- tau^2
  
  #R^2 
  # acknowledment: code taken from the lmrob() function from the robustbase 
  # package
  resid <- model$resid
  pred <- model$fitted.values
  resp <- pred + resid 
  wgt <- weights 
  if (is.null(model$model)) {
    df.int <- if (any(colnames(X) == "(Intercept)")) { 1L } else { 0L }
  } else {
    df.int <- if (attr(model$terms, "intercept")) { 1L } else { 0L }
  }
  resp.mean <- if (df.int == 1L) { sum(wgt * resp) / sum(wgt) } else { 0L }
  yMy <- sum(wgt * (resp - resp.mean)^2)
  rMr <- sum(wgt * resid^2)
  # bi-square correction
  correc <- 1.207617 
  R2.org <- (yMy - rMr) / (yMy + rMr * (correc - 1))
  
  # compute the loglikelihood
  ll.unc <- con_loglik_lm(X = X, y = Y, b = b.unconstr, w = w)
  LL.unc <- ll.unc$loglik
  
  # # check if the constraints are in line with the data    
  if (all(Amat %*% c(b.unconstr) - bvec >= 0 * bvec) & meq == 0) {  
    b.constr <- b.unconstr
        
    OUT <- list(CON = NULL,
                parTable = parTable,
                constraints = constraints,
                b.unconstr = b.unconstr,
                b.constr = b.unconstr,
                residuals = model$resid,
                wresid = model$wresid,
                #init.resid = resid0,
                fitted.values = fitted(model),
                weights = w,
                R2.org = R2.org,
                R2.reduced = R2.org,
                df.residual = so$df[2],
                scale = model$s, 
                s2.unc = s2.unc, 
                s2.con = s2.unc, 
                loglik = LL.unc, 
                Sigma = vcov(model),                                            #probably not so robust!
                Amat = Amat, bvec = bvec, meq = meq, iact = 0L,
                converged = model$converged, iter = model$iter,
                bootout = NULL, call = cl)
  } else {
    # update model
    call.my <- list(Amat = Amat, meq = meq, bvec = bvec)
    CALL <- c(list(model), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    b.constr <- rfit$coefficients
      b.constr[abs(b.constr) < sqrt(.Machine$double.eps)] <- 0L
    iact <- rfit$iact
    
    ll.con <- con_loglik_lm(X = X, y = Y, b = b.unconstr, w = w)
    LL.con <- ll.con$loglik
        
    so.rlm <- summary_rlm(rfit, Amat = Amat, meq = meq)
    tau <- so.rlm$stddev
    s2.con <- tau^2
    
    #R^2
    pred <- rfit$fitted.values
    resid <- rfit$residuals
    resp <- pred + resid 
    wgt <- rfit$weights
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
    R2.reduced <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    #R2.adjusted <- 1 - (1 - r2correc) * ((n - df.int) / df.residual)
  
    OUT <- list(CON = NULL,
                parTable = parTable,
                constraints = constraints,
                b.unconstr = b.unconstr,
                b.constr = b.constr,
                residuals = resid,
                #init.resid = resid0,
                wresid = rfit$wresid,
                fitted.values = pred,
                weights = w,
                R2.org = R2.org,
                R2.reduced = R2.reduced,
                df.residual = so$df[2], # correction for equalities constraints is included.
                scale = rfit$s,                                                               
                s2.unc = s2.unc, 
                s2.con = s2.con, 
                loglik = LL.con, 
                Sigma = vcov(model),                                             #probably not robust???
                Amat = Amat, bvec = bvec, meq = meq, iact = iact,
                converged = rfit$converged, iter = rfit$iter,
                bootout = NULL, call = cl)
  }
  
  OUT$model.org <- model
  OUT$CON <- if (is.character(constraints)) { CON }
  OUT$se <- se 
  
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      information <- 1/s2.con * crossprod(X) 
      OUT$information <- information      
      #only for MM-estimation? 
      #information <- solve(vcovMM(X, resid0, resid, scale))                     
      inverted.information <- con_augmented_information(information = information,
                                                        X = X, 
                                                        b.unconstr = b.unconstr, 
                                                        b.constr = b.constr,
                                                        Amat = Amat, 
                                                        bvec = bvec, meq = meq)
      attr(OUT$information, "inverted.information") <- inverted.information        
    } else if (se == "boot.model.based") { 
      OUT$bootout <- con_boot_lm(model, B = ifelse(is.null(control$B),
                                                   999, control$B), 
                                 fixed = TRUE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(model, B = ifelse(is.null(control$B),
                                                   999, control$B),
                                 fixed = FALSE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    }
  }
  
  class(OUT) <- c("conRLM","rlm")
  
  return(OUT)
}
