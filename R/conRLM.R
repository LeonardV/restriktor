#compute restrikted robust estimates
conRLM.rlm <- function(object, constraints = NULL, se = "standard", 
                       B = 999, rhs = NULL, neq = 0L, Wt = TRUE,
                       bootWt = FALSE, bootWt.R = 99999,
                       parallel = "no", ncpus = 1L, cl = NULL, 
                       seed = NULL, control = NULL, verbose = FALSE, 
                       debug = FALSE, ...) { 
  
  # check class
  if (!(class(object)[1] == "rlm")) {
    stop("Restriktor ERROR: object must be of class rlm.")
  }
  
  
  # timing
  start.time0 <- start.time <- proc.time()[3]; timing <- list()
  
  # store call
  mc <- match.call()
  
  # rename for internal use
  Amat <- constraints
  bvec <- rhs
  meq <- neq
  
  # response varialbe
  y <- as.matrix(object$model[, attr(object$terms, "response")])
  # model matrix
  X  <- model.matrix(object)[ , ,drop = FALSE]
  # model summary
  so.org <- summary(object)
  # unrestrikted coefficients
  b.unrestr <- coef(object)
  # vcov
  Sigma <- vcov(object) #solve(invW)
  # unrestrikted scale estimate for the standard deviation
  tau <- so.org$stddev
  # residual degrees of freedom
  #rdf <- so.org$df[2]
  # residuals
  residuals <- object$residuals
  # sampel size
  n <- dim(X)[1]
  # number of parameters
  p <- dim(X)[2]
  
  timing$preparation <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # deal with constraints
  if (!is.null(constraints)) {
    restr_OUT <- con_constraints(object, 
                                 constraints = Amat, 
                                 bvec        = bvec, 
                                 meq         = meq, 
                                 debug       = debug)  
    # a list with useful information about the restriktions.}
    CON <- restr_OUT$CON
    # a parameter table with information about the observed variables in the object 
    # and the imposed restriktions.}
    parTable <- restr_OUT$parTable
    # constraints matrix
    Amat <- restr_OUT$Amat
    # rhs 
    bvec <- restr_OUT$bvec
    # neq
    meq <- restr_OUT$meq
  } else if (is.null(constraints)) { # no constraints specified - needed for GORIC to include unconstrained model
    CON <- NULL
    parTable <- NULL
    Amat <- rbind(rep(0L, p))
    bvec <- rep(0L, nrow(Amat))
    meq  <- 0L
  } 
  
  # compute the reduced row-echelon form of the constraints matrix
  rAmat <- GaussianElimination(t(Amat))
  if (Wt && !bootWt) {
    if (rAmat$rank < nrow(Amat) && rAmat$rank != 0L) {
      stop(paste("Restriktor ERROR: The constraint matrix must have full row-rank ( choose e.g. rows", 
                 paste(rAmat$pivot, collapse = " "), ", or try to set bootWt = TRUE)"))
    }
  } else {
    if (rAmat$rank < nrow(Amat) && 
        !(se %in% c("none", "boot.model.based", "boot.standard")) && 
        rAmat$rank != 0L) {
      se <- "none"
      warning(paste("Restriktor Warning: No standard errors could be computed. 
                      The constraint matrix must have full row-rank ( choose e.g. rows", 
                    paste(rAmat$pivot, collapse = " "), "). 
                      Try to set se = \"boot.model.based\" or \"boot.standard\"."))  
    }
  }
  
  timing$constraints <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # standard error stuff
  if (se == "default") {
    se <- "standard"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (!(se %in% c("none","standard","const","boot.model.based","boot.standard","HC","HC0",
                  "HC1","HC2","HC3","HC4","HC4m","HC5"))) {
    stop("Restriktor ERROR: standard error method ", sQuote(se), " unknown.")
  }
  if (se == "boot.model.based" & any(Amat[,1] == 1)) { 
    stop("Restriktor ERROR: no restriktions on intercept possible for 'se = boot.model.based' bootstrap method.")
  }
  if(ncol(Amat) != length(b.unrestr)) {
    stop("Restriktor ERROR: the columns of \"constraints\" does not match with the number of parameters.")
  }
  
  
  ## extract all available options from the original model, else
  ## asign default settings from rlm. Except for the psi function, because
  ## only the tukey's bisquare loss function is supported by restriktor.
  
  # original model function call
  call <- as.list(object$call)
  # which method M or MM estimation
  method <- call[["method"]]
  if (is.null(method)) object$call[["method"]] <- "M"
  # check (only tukey's bisquare supported)
  psi <- call[["psi"]]
  if (is.null(psi)) {
    if (object$call[["method"]] == "M") {
      stop("Restriktor ERROR: only tukey's bisquare loss function is supported.")
    }
  } else {
    if (psi != "psi.bisquare") {
      stop("Restriktor ERROR: only tukey's bisquare loss function is supported.")
    }
  }
  # tukeys bisquare tuning constant
  cc <- call[["c"]]
  if (is.null(cc)) { object$call[["c"]] <- 4.685061 }
  # prior weights for each case
  weights <- object$weights
  if (any(weights != 1L)) stop("Restriktor ERROR: prior weights are not implemented (yet).")            # FIXME
  # weights method, inverse of the variance or case 
  wt.method <- call[["wt.method"]]
  if (is.null(wt.method)) { object$call[["wt.method"]] <- "inv.var" }
  # initial down-weighting for each case
  w <- call[["w"]]
  if (is.null(w)) { object$call[["w"]] <- rep(1, length(weights)) }
  # scale estimator - depends on method
  scale.est <- call[["scale.est"]]
  if (is.null(scale.est)) { object$call[["scale.est"]] <- "MAD" }
  # tuning constant used for Huber proposal 2 scale estimation.
  k2 <- call[["k2"]]
  if (is.null(k2)) { object$call[["k2"]] <- 1.345 }
  # the stopping criterion is based on changes in this vector
  test.vec <- call[["test.vec"]]
  if (is.null(test.vec)) { object$call[["test.vec"]] <- "resid" }
  # the limit on the number of IWLS iterations.
  maxit <- call[["maxit"]]
  if (is.null(maxit)) { object$call[["maxit"]] <- 50 }
  # the accuracy for the stopping criterion.
  acc <- call[["acc"]]
  if (is.null(acc)) { object$call[["acc"]] <- 1e-4 }
  
  timing$mc_org <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
 
  
  # compute R-squared 
  # acknowledment: code taken from the lmrob() function from the robustbase package
  df.int <- ifelse(attr(object$terms, "intercept"), 1L, 0L)
  y.mean <- if (df.int == 1L) { 
    sum(weights * y) / sum(weights) 
    } else { 0L }
  yMy <- sum(weights * (y - y.mean)^2)
  rMr <- sum(weights * residuals^2)
  # tukey's bi-square correction
  correc <- 1.207617 
  R2.org <- (yMy - rMr) / (yMy + rMr * (correc - 1))
  
  timing$R2 <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # compute the loglikelihood
  ll.unc <- con_loglik_lm(X = X, 
                          y = y, 
                          b = b.unrestr, 
                          w = weights)
  LL.unc <- ll.unc$loglik
  
  timing$LLik <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  # ML unconstrained MSE
  #tau2ml <- (tau^2 * so.org$df[2]) / n
  #invW <- kronecker(solve(tau2ml), t(X) %*% X)
  
  ## compute mixing weights
  is.augmented <- TRUE
  if (Wt) {
    ## compute mixing weights
    if (all(c(Amat) == 0)) { # unrestrikted case
      wt <- c(rep(0L, p), 1)
      is.augmented <- FALSE
    } else if (bootWt) { # compute mixing weights based on simulation
      wt <- con_weightsBoot(VCOV     = Sigma,
                            Amat     = Amat, 
                            meq      = meq, 
                            R        = bootWt.R,
                            parallel = parallel,
                            ncpus    = ncpus,
                            cl       = cl,
                            seed     = seed,
                            verbose  = verbose)
    } else if (!bootWt & (meq < nrow(Amat))) { # compute mixing weights based on mvnorm
      #check
    #  if ((qr(Amat)$rank < nrow(Amat))) {
    #    stop("Restriktor ERROR: constraint matrix must have full row-rank. try set bootWt = TRUE.")
    #  }
      wt <- rev(con_weights(Amat %*% Sigma %*% t(Amat), meq = meq))
    } else if (!bootWt & (meq == nrow(Amat))) { # only equality constraints
      wt <- rep(0L, ncol(Sigma) + 1)
      wt.idx <- ncol(Sigma) - meq + 1
      wt[wt.idx] <- 1
    }
    attr(wt, "bootWt") <- bootWt
    if (bootWt) { attr(wt, "bootWt.R") <- bootWt.R }
  } else {
    wt <- NULL
  }
  
  timing$wt <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # # check if the constraints are in line with the data    
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) & meq == 0) {  
    b.restr <- b.unrestr
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                wt          = wt,
                b.unrestr   = b.unrestr,
                b.restr     = b.unrestr,
                residuals   = object$residuals,
                wresid      = object$wresid,
                fitted      = object$fitted,
                weights     = object$weights,
                w           = object$w, 
                scale       = object$s, 
                scale.restr = object$s,
                psi         = object$psi,
                R2.org      = R2.org,
                R2.reduced  = R2.org,
                df.residual = so.org$df[2],
                s2.unc      = tau^2, 
                s2.restr    = tau^2, 
                loglik      = LL.unc, 
                Sigma       = Sigma,                                            #probably not so robust!
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                iact        = 0L,
                converged   = object$converged, 
                iter        = object$iter,
                bootout     = NULL, 
                control     = control)
  } else {
    # add constraints to model 
    call.my <- list(Amat = Amat, 
                    meq  = meq, 
                    bvec = bvec, 
                    tol  = ifelse (is.null(control$tol), 
                                   sqrt(.Machine$double.eps), control$tol))
    # collect all original model arguments and add constraints
    CALL <- c(list(object), call.my)
    CALL <- CALL[!duplicated(CALL)]
    # fit constraint robust liner model
    rfit <- do.call("conRLM_fit", CALL)
    
    timing$conRLM_fit <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    # constrained coefs
    b.restr <- rfit$coefficients
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    fitted <- X %*% b.restr
    residuals <- rfit$residuals
    # weights
    weights <- rfit$weights
    # psi(resid/scale) these are the weights used for downweighting the cases.
    w <- rfit$w
    
    # compute loglik
    ll.restr <- con_loglik_lm(X = X, 
                              y = y, 
                              b = b.unrestr, 
                              w = weights)
    LL.restr <- ll.restr$loglik
    
    timing$LLik <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    # in case of equality constraints we need to correct the residual df 
    if (length(object$call[["wt.method"]]) && object$call[["wt.method"]] == "case") {
      rdf <- sum(weights) - p
    #  w   <- rfit$psi(rfit$wresid/rfit$scale)
      S   <- sum(weights * (rfit$wresid * w)^2) / rdf
      std <- summary(rfit)$stddev / sqrt(S)
      
      if (!is.null(Amat)) {
        rdf <- sum(weights) - (p - qr(Amat[0:meq,])$rank)
        S.new <- sum(weights * (rfit$wresid * w)^2) / rdf
        tau.restr <- (summary(rfit)$stddev / sqrt(S)) * sqrt(S.new)
      } else {
        tau.restr <- std * sqrt(S)
      }
    } else {
      rdf <- n - p
     # w   <- rfit$psi(rfit$wresid/rfit$scale)
      S   <- sum((rfit$wresid * w)^2) / rdf
      std <- summary(rfit)$stddev / sqrt(S)
      
      if (!is.null(Amat)) {
        rdf <- n - (p - qr(Amat[0:meq,])$rank)
        S.new <- sum((rfit$wresid * w)^2) / rdf
        tau.restr <- (summary(rfit)$stddev / sqrt(S)) * sqrt(S.new)
      } else {
        tau.restr <- std * sqrt(S)
      }
    }
    
    timing$df_correction <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
      
    #R^2 under the restrikted object
    df.int <- if (attr(object$terms, "intercept")) { 1L } else { 0L }
    resp.mean <- if (df.int == 1L) { sum(weights * y) / sum(weights) } else { 0 }
    yMy <- sum(weights * (y - resp.mean)^2)
    rMr <- sum(weights * residuals^2)
    # tukey's bi-square correction
    correc <- 1.207617 
    R2.reduced <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    
    timing$R2 <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                wt          = wt,
                b.unrestr   = b.unrestr,
                b.restr     = b.restr,
                residuals   = residuals,
                wresid      = rfit$wresid,
                fitted      = fitted,
                weights     = weights,
                w           = w, 
                scale       = object$s,
                scale.restr = rfit$scale,
                R2.org      = R2.org,
                R2.reduced  = R2.reduced,
                df.residual = so.org$df[2], 
                s2.unc      = tau^2, 
                s2.restr    = tau.restr^2, 
                loglik      = LL.restr, 
                Sigma       = Sigma,                                             #probably not so robust???
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                iact        = rfit$iact,
                converged   = rfit$converged, 
                iter        = rfit$iter,
                bootout     = NULL, 
                control     = control)
  }
  
  OUT$model.org <- object
  OUT$se <- se 
  
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      #  V <- vcovMM(X = X, resid0 = resid0, residuals = residuals, scale = model$s)  
      OUT$information <- 1/tau^2 * crossprod(X)
      information.inv <- con_augmented_information(information  = OUT$information,
                                                   is.augmented = is.augmented,
                                                   X            = X, 
                                                   b.unrestr    = b.unrestr, 
                                                   b.restr      = b.restr,
                                                   Amat         = Amat, 
                                                   bvec         = bvec, 
                                                   meq          = meq)
          
      attr(OUT$information, "inverted")  <- information.inv$information
      attr(OUT$information, "inverted.augmented") <- information.inv$information.augmented
      
      timing$inv_aug_information <- (proc.time()[3] - start.time)
      start.time <- proc.time()[3]
    } else if (se == "boot.model.based") { 
      OUT$bootout <- con_boot_lm(object      = object, 
                                 B           = B, 
                                 fixed       = TRUE,
                                 constraints = Amat,
                                 rhs         = bvec, 
                                 neq         = meq, 
                                 se          = "none",
                                 bootWt      = bootWt,
                                 bootWt.R    = bootWt.R,
                                 parallel    = parallel, 
                                 ncpus       = ncpus, 
                                 cl          = cl)
      timing$boot_model_based <- (proc.time()[3] - start.time)
      start.time <- proc.time()[3]
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(object      = object, 
                                 B           = B, 
                                 fixed       = FALSE,
                                 constraints = Amat,
                                 rhs         = bvec, 
                                 neq         = meq, 
                                 se          = "none",
                                 bootWt      = bootWt,
                                 bootWt.R    = bootWt.R,
                                 parallel    = parallel, 
                                 ncpus       = ncpus, 
                                 cl          = cl)
      timing$boot_standard <- (proc.time()[3] - start.time)
      start.time <- proc.time()[3]
    }
  }
  
  OUT$timing$total <- (proc.time()[3] - start.time0)
  
  class(OUT) <- c("conRLM","conLM","rlm","lm")
  
  OUT

}
