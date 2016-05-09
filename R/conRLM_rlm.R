#compute constrained robust estimates
conRLM.rlm <- function(model, constraints, debug = FALSE,
                       se = "default", B = 999, 
                       bvec = NULL, meq = 0L, 
                       tol = sqrt(.Machine$double.eps), ...) { 
  
  cl <- match.call()
  Amat <- Amatw; bvec <- bvecw; meq <- meqw
  
  # check class
  if (!("rlm" %in% class(model))) {
    stop("Restriktor ERROR: model must be of class rlm.")
  }
  if (se == "default") {
    se <- "standard"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (!(se %in% c("none","standard","const","boot.model.based","boot.standard","HC","HC0",
                  "HC1","HC2","HC3","HC4","HC4m","HC5"))) {
    stop("standard error method ", sQuote(se), " unknown.")
  }
  if (se == "boot.model.based" & any(Amat[,1] == 1)) { 
    stop("no restriktions on intercept possible for 'se = boot.model.based' bootstrap method.")
  }
  if (missing(constraints) && is.null(bvec)) { 
    bvec <- rep(0L, nrow(Amat)) 
  }
  
  y <- as.matrix(model$model[, attr(model$terms, "response")])
  X  <- model.matrix(model)[,,drop = FALSE]
  # unconstrained coefficients
  b.unrestr <- coef(model)  
  
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
  # prior weights for each case
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
  # initial down-weighting for each case
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
  
  # check
  if (method == "M") {
    if (psi != "psi.bisquare") {
      stop("restriktor only supports the bisquare loss function (for now).")
    }
  }
  if(ncol(Amat) != length(b.unrestr)) {
    stop("length coefficients and ncol(constraints) must be identical")
  }
  
  # Adjusted summary function from MASS:::summary.rlm().
  # The df needs to be corrected for equality constraints and
  # we need the df for ML. 
  so.org <- summary_rlm(model)
  # unconstrained scale estimate for the standard errors
  tau.unc <- tau <- so.org$stddev
  
  #R^2 
  # acknowledment: code taken from the lmrob() function from the robustbase 
  # package
  residuals <- model$residuals
  pred <- model$fitted.values
  resp <- pred + residuals 
  wgt <- weights 
  if (is.null(model$model)) {
    df.int <- if (any(colnames(X) == "(Intercept)")) { 1L } else { 0L }
  } else {
    df.int <- if (attr(model$terms, "intercept")) { 1L } else { 0L }
  }
  resp.mean <- if (df.int == 1L) { sum(wgt * resp) / sum(wgt) } else { 0L }
  yMy <- sum(wgt * (resp - resp.mean)^2)
  rMr <- sum(wgt * residuals^2)
  # bi-square correction
  correc <- 1.207617 
  R2.org <- (yMy - rMr) / (yMy + rMr * (correc - 1))
  
  # compute the loglikelihood
  ll.unc <- con_loglik_lm(X = X, y = y, b = b.unrestr, w = weights)
  LL.unc <- ll.unc$loglik
  # # check if the constraints are in line with the data    
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) & meq == 0) {  
    b.restr <- b.unrestr
    # To compute the vcovMM, we need the initial residuals from the S-estimator.
    # These are not available in the original rlm object.
    #S <- Sestimator(x = X, y = y)
    #resid0 <- residuals(S)
    OUT <- list(CON = NULL,
                parTable = parTable,
                #constraints = constraints,
                b.unrestr = b.unrestr,
                b.restr = b.unrestr,
                residuals = model$residuals,
                wresid = model$wresid,
                #init.residuals = resid0,
                fitted.values = model$fitted,
                weights = model$weights,
                w = model$w, 
                scale = model$s, 
                psi = model$psi,
                R2.org = R2.org,
                R2.reduced = R2.org,
                df.residual = so.org$df[2],
                s2.unc = tau^2, 
                s2.restr = tau^2, 
                loglik = LL.unc, 
                Sigma = vcov(model),                                            #probably not so robust!
                Amat = Amat, bvec = bvec, meq = meq, iact = 0L,
                converged = model$converged, iter = model$iter,
                bootout = NULL, call = cl)
  } else {
    # add constraints to model
    call.my <- list(Amat = Amat, meq = meq, bvec = bvec)
    # collect all original model arguments and add constraints
    CALL <- c(list(model), call.my)
    # fit constraint robust liner model
    rfit <- do.call("conRLM_fit", CALL)
    # constrained coefs
    b.restr <- rfit$coefficients
    b.restr[abs(b.restr) < sqrt(.Machine$double.eps)] <- 0L
    # number of active inequality constraints
    iact <- rfit$iact
    # weights
    weights <- rfit$weights
    # psi(resid/scale) these are the weights used for downweighting the cases.
    w <- rfit$w
    #resid0 <- rfit$resid0
    # compute loglik
    ll.restr <- con_loglik_lm(X = X, y = y, b = b.unrestr, w = weights)
    LL.restr <- ll.restr$loglik
    # we need the std. deviation adjusted for equality constraints (meq > 0)
    so.restr <- summary_rlm(rfit, Amat = Amat, meq = meq)
    tau <- so.restr$stddev
    
    #R^2 under the constrained model
    # predicted values under constraints
    pred <- rfit$fitted.values
    # residuals under constraints
    residuals <- rfit$residuals
    # response under constraints
    resp <- pred + residuals 
    # weights do not change under constraints
    wgt <- rfit$weights
    df.int <- if (attr(model$terms, "intercept")) { 1L } else { 0L }
    resp.mean <- if (df.int == 1L) { sum(wgt * resp)/sum(wgt) } else { 0 }
    yMy <- sum(wgt * (resp - resp.mean)^2)
    rMr <- sum(wgt * residuals^2)
    # bi-square correction
    correc <- 1.207617 
    R2.reduced <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    
    OUT <- list(CON = NULL,
                parTable = parTable,
                #constraints = constraints,
                b.unrestr = b.unrestr,
                b.restr = b.restr,
                residuals = residuals,
                #init.residuals = resid0,
                wresid = rfit$wresid,
                fitted.values = pred,
                weights = weights,
                w = w, 
                scale = model$s,
                R2.org = R2.org,
                R2.reduced = R2.reduced,
                df.residual = so.restr$df[2], # correction for equalities constraints is included.
                #df.residual = so.org$df[2], 
                s2.unc = tau.unc^2, 
                s2.restr = tau^2, 
                loglik = LL.restr, 
                Sigma = vcov(model),                                             #probably not so robust???
                Amat = Amat, bvec = bvec, meq = meq, iact = iact,
                converged = rfit$converged, iter = rfit$iter,
                bootout = NULL, call = cl)
  }
  
  OUT$model.org <- model
  OUT$CON <- if (is.character(constraints)) { CON }
  OUT$se <- se 
  
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      #if (model$call[["method"]] == "MM") {
      #  V <- vcovMM(X = X, resid0 = resid0, residuals = residuals, scale = model$s)  
      #  OUT$information <- solve(V)
      #} else {
      #OUT$information <- tau^2 * cov.unscaled 
      OUT$information <- 1/tau^2 * crossprod(X)
      #}
      inverted.information <- con_augmented_information(information = OUT$information,
                                                        X = X, 
                                                        b.unrestr = b.unrestr, 
                                                        b.restr = b.restr,
                                                        Amat = Amat, 
                                                        bvec = bvec, meq = meq)
      
      attr(OUT$information, "inverted.information") <- inverted.information        
    } else if (se == "boot.model.based") { 
      OUT$bootout <- con_boot_lm(model, B = B, 
                                 fixed = TRUE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(model, B = B,
                                 fixed = FALSE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    }
  }
  
  class(OUT) <- c("conRLM","conLM","rlm","lm")
  
  return(OUT)
}
