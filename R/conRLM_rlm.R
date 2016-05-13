#compute restrikted robust estimates
conRLM.rlm <- function(model, constraints, se = "default", B = 999, 
                       rhs = NULL, neq = 0L, bootWt = FALSE, R = 99999,
                       parallel = "no", ncpus = 1L, cl = NULL, seed = NULL, 
                       control = NULL, verbose = FALSE, debug = FALSE, ...) { 
  
  # check class
  if (!(class(model)[1] == "rlm")) {
    stop("Restriktor ERROR: model must be of class rlm.")
  }
  
  cl <- match.call()
  # rename for internal use
  bvec <- rhs; meq <- neq
  # construct constraint matrix/vector.
  restr_OUT <- con_constraints(model, constraints = constraints, 
                               bvec = bvec, meq = meq, debug = debug)  
  # a list with useful information about the restriktions.}
  CON <- restr_OUT$CON
  # a parameter table with information about the observed variables in the model 
  # and the imposed restriktions.}
  parTable <- restr_OUT$parTable
  # constraints matrix
  Amat <- restr_OUT$Amat
  # rhs 
  bvec <- restr_OUT$bvec
  # neq
  meq <- restr_OUT$meq
  
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
  
  ## extract all available options from the original model, else
  ## asign default settings from rlm. Except for the psi function, because
  ## only the tukey's bisquare loss function is supported by restriktor.
  
  # original model function call
  call <- as.list(model$call)
  # which method M or MM estimation
  method <- call[["method"]]
  if (is.null(method)) model$call[["method"]] <- "M"
  # check (only tukey's bisquare supported)
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
  if (any(weights != 1L)) stop("weights are not implemented (yet).")            # FIXME
  # weights method, inverse of the variance or case 
  wt.method <- call[["wt.method"]]
  if (is.null(wt.method)) model$call[["wt.method"]] <- "inv.var"
  # initial down-weighting for each case
  w <- call[["w"]]
  if (is.null(w)) model$call[["w"]] <- rep(1, length(weights))
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
  
  # response varialbe
  y <- as.matrix(model$model[, attr(model$terms, "response")])
  # model matrix
  X  <- model.matrix(model)[ , ,drop = FALSE]
  # model summary
  so.org <- summary(model)
  # unrestrikted coefficients
  b.unrestr <- coef(model)  
  # unrestrikted scale estimate for the standard deviation
  tau <- so.org$stddev
  # residual degrees of freedom
  #rdf <- so.org$df[2]
  # residuals
  residuals <- model$residuals
  # sampel size
  n <- length(residuals)
  
  # check
  if(ncol(Amat) != length(b.unrestr)) {
    stop("restriktor ERROR: length coefficients and ncol(Amat) must be identical")
  }
  
  # compute R-squared 
  # acknowledment: code taken from the lmrob() function from the robustbase package
  df.int <- ifelse(attr(model$terms, "intercept"), 1L, 0L)
  y.mean <- if (df.int == 1L) { 
    sum(weights * y) / sum(weights) 
    } else { 0L }
  yMy <- sum(weights * (y - y.mean)^2)
  rMr <- sum(weights * residuals^2)
  # tukey's bi-square correction
  correc <- 1.207617 
  R2.org <- (yMy - rMr) / (yMy + rMr * (correc - 1))
  
  # compute the loglikelihood
  ll.unc <- con_loglik_lm(X = X, y = y, b = b.unrestr, w = weights)
  LL.unc <- ll.unc$loglik
  
  # ML unconstrained MSE
  tau2ml <- (tau^2 * so.org$df[2]) / n
  invW <- kronecker(solve(tau2ml), t(X) %*% X)
  W <- solve(invW)
  # compute mixing weights
  if (bootWt) { # compute mixing weights based on simulation
    wt <- mix.boot(VCOV = W,
                   Amat = Amat, 
                   meq = meq, 
                   R = R,
                   parallel = parallel,
                   ncpus = ncpus,
                   cl = cl,
                   seed = seed,
                   verbose = verbose)
  } else if (!bootWt & (meq < nrow(Amat))) { # compute mixing weights based on mvnorm
    wt <- rev(con_wt(Amat %*% W %*% t(Amat), meq = meq))
  } else if (!bootWt & (meq == nrow(Amat))) { # only equality constraints
    wt <- rep(0L, ncol(W) + 1)
    wt.idx <- ncol(W) - meq + 1
    wt[wt.idx] <- 1
  }
  attr(wt, "bootWt") <- bootWt
  
  
  # # check if the constraints are in line with the data    
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) & meq == 0) {  
    b.restr <- b.unrestr
    OUT <- list(CON = NULL,
                parTable = parTable,
                wt = wt,
                b.unrestr = b.unrestr,
                b.restr = b.unrestr,
                residuals = model$residuals,
                wresid = model$wresid,
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
    CALL <- c(list(model), call.my)
    CALL <- CALL[!duplicated(CALL)]
    # fit constraint robust liner model
    rfit <- do.call("conRLM_fit", CALL)
    # constrained coefs
    b.restr <- rfit$coefficients
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    residuals <- rfit$residuals
    # weights
    weights <- rfit$weights
    # psi(resid/scale) these are the weights used for downweighting the cases.
    w <- rfit$w
    # compute loglik
    ll.restr <- con_loglik_lm(X = X, y = y, b = b.unrestr, w = weights)
    LL.restr <- ll.restr$loglik
    
    # adjusted summary.rlm function to get the right residual degrees of freedom
    # and standard deviation in case of equality restriktions.
    so.restr <- summary_rlm(rfit, Amat = Amat, meq = meq)
    #rdf <- so.restr$df[2]
    tau.restr <- so.restr$stddev
    
    #R^2 under the restrikted model
    df.int <- if (attr(model$terms, "intercept")) { 1L } else { 0L }
    resp.mean <- if (df.int == 1L) { sum(weights * y) / sum(weights) } else { 0 }
    yMy <- sum(weights * (y - resp.mean)^2)
    rMr <- sum(weights * residuals^2)
    # tukey's bi-square correction
    correc <- 1.207617 
    R2.reduced <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    
    OUT <- list(CON = NULL,
                parTable = parTable,
                wt = wt,
                b.unrestr = b.unrestr,
                b.restr = b.restr,
                residuals = residuals,
                wresid = rfit$wresid,
                fitted.values = fitted,
                weights = weights,
                w = w, 
                scale = model$s,
                R2.org = R2.org,
                R2.reduced = R2.reduced,
                df.residual = so.org$df[2], 
                s2.unc = tau^2, 
                s2.restr = tau.restr^2, 
                loglik = LL.restr, 
                Sigma = vcov(model),                                             #probably not so robust???
                constraints = Amat, rhs = bvec, neq = meq, iact = rfit$iact,
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
      OUT$bootout <- con_boot_lm(model, B = B, fixed = TRUE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(model, B = B, fixed = FALSE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    }
  }
  
  class(OUT) <- c("conRLM","conLM","rlm","lm")
  
  return(OUT)
}
