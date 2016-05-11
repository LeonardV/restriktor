#compute constrained robust estimates
conRLM.rlm <- function(model, constraints, se = "default", B = 999, 
                       bvec = NULL, meq = 0L, control = NULL, ...) { 
  
  cl <- match.call()
  Amat <- Amatw; bvec <- bvecw; meq <- meqw
  
  # check class
  if (!("rlm" %in% class(model))) {
    stop("Restriktor ERROR: model must be of class rlm.")
  }
  # standard error stuff
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
  b.unrestr <- coef(model)  
  
  ## extract all available settings from the original model and else
  ## asign default settings from rlm. Except for the psi function, because
  ## only the tukey's bisquare loss function is supported by restriktor.
  # original model function call
  call <- as.list(model$call)
  # which method M or MM estimation
  method <- call[["method"]]
  if (is.null(method)) model$call[["method"]] <- "M"
  
  # check (only tukey's bisquare)
  # psi function 
  psi <- call[["psi"]]
  if (is.null(psi)) {
    if (method == "M") {
      stop("only tukey's bisquare loss function is supported.")
    } 
  } else {
    if (psi != "psi.bisquare") {
      stop("only tukey's bisquare loss function is supported.")
    }
  }
  # tukeys bisquare tuning constant
  cc <- call[["c"]]
  if (is.null(cc)) model$call[["c"]] <- 4.685061
  # prior weights for each case
  weights <- model$weights
  if (any(weights != 1L)) stop("weights are not implemented (yet).")
  # weights method, inverse of the variance or case 
  wt.method <- call[["wt.method"]]
  if (is.null(wt.method)) model$call[["wt.method"]] <- "inv.var"
  # initial down-weighting for each case
  w <- call[["w"]]
  if (is.null(w)) model$call[["w"]] <- rep(1, nrow(X))
  # scale estimator - depends on method
  scale.est <- call[["scale.est"]]
  if (is.null(scale.est)) model$call[["scale.est"]] <- "MAD"
  # tuning constant used for Huber proposal 2 scale estimation.
  k2 <- call[["k2"]]
  if (is.null(k2)) model$call[["k2"]] <- 1.345
  # the stopping criterion is based on changes in this vector
  test.vec <- call[["test.vec"]]
  if (is.null(test.vec)) model$call[["test.vec"]] <- "resid"
  # the limit on the number of IWLS iterations.
  maxit <- call[["maxit"]]
  if (is.null(maxit)) model$call[["maxit"]] <- 20
  # the accuracy for the stopping criterion.
  acc <- call[["acc"]]
  if (is.null(acc)) model$call[["acc"]] <- 1e-4
  
  # check
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
                constraints = Amat, rhs = bvec, neq = meq, iact = 0L,
                converged = model$converged, iter = model$iter,
                bootout = NULL, call = cl)
  } else {
    # add constraints to model 
    call.my <- list(Amat = Amat, meq = meq, bvec = bvec, 
                    tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol))
    # collect all original model arguments and add constraints
    #CALL <- c(list(model), call.my)
    CALL <- c(list(model), call.my)
    CALL <- CALL[!duplicated(CALL)]
    # fit constraint robust liner model
    rfit <- do.call("conRLM_fit", CALL)
    # constrained coefs
    b.restr <- rfit$coefficients
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
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
                constraints = Amat, rhs = bvec, neq = meq, iact = iact,
                converged = rfit$converged, iter = rfit$iter,
                bootout = NULL, call = cl)
  }
  
  OUT$model.org <- model
  OUT$CON <- if (is.character(constraints)) { CON }
  OUT$se <- se 
  
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      #  V <- vcovMM(X = X, resid0 = resid0, residuals = residuals, scale = model$s)  
      OUT$information <- 1/tau^2 * crossprod(X)
      information <- con_augmented_information(information = OUT$information,
                                               X = X, 
                                               b.unrestr = b.unrestr, 
                                               b.restr = b.restr,
                                               Amat = Amat, 
                                               bvec = bvec, meq = meq)
      attr(OUT$information, "inverted.information")  <- information$inverted.information
      attr(OUT$information, "augmented.information") <- information$augmented.information        
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
