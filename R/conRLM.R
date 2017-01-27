#compute restrikted robust estimates
conRLM.rlm <- function(object, constraints = NULL, se = "standard", 
                       B = 999, rhs = NULL, neq = 0L, 
                       Wt = "mvnorm", bootWt.R = 99999,
                       parallel = "no", ncpus = 1L, cl = NULL, 
                       seed = NULL, control = list(), verbose = FALSE, 
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
  meq  <- neq
  
  # response varialbe
  y <- as.matrix(object$model[, attr(object$terms, "response")])
  # model matrix
  X  <- model.matrix(object)[ , ,drop = FALSE]
  # model summary
  so <- summary(object)
  # weights
  weights <- weights(object)
  # unrestrikted coefficients
  b_unrestr <- coef(object)
  # vcov
  Sigma <- vcov(object) 
  # unrestrikted scale estimate for the standard deviation: 
  # tau^2 * solve(t(X)%*%X) equals vcov(object)
  tau <- so$stddev
  # residuals
  residuals <- residuals(object) # NOT working residual
  
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
    # a list with useful information about the restrictions.}
    CON <- restr_OUT$CON
    # a parameter table with information about the observed variables in the object 
    # and the imposed restrictions.}
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
  
  
  if (!(Wt %in% c("mvnorm", "boot", "none"))) {
    stop("Restriktor ERROR: ", sQuote(Wt), " method unknown. Choose from \"mvnorm\", \"boot\", or \"none\"")
  }
  
  # compute the reduced row-echelon form of the constraints matrix
  rAmat <- GaussianElimination(t(Amat))
  if (Wt == "mvnorm") {
    if (rAmat$rank < nrow(Amat) && rAmat$rank != 0L) {
      stop(paste("Restriktor ERROR: The constraint matrix must have full row-rank ( choose e.g. rows", 
                 paste(rAmat$pivot, collapse = " "), ", or try to set Wt = \"boot\""))
    }
  } else {
    if (rAmat$rank < nrow(Amat) && 
        !(se %in% c("none", "boot.model.based", "boot.standard")) && 
        rAmat$rank != 0L) {
      se <- "none"
      warning(paste("Restriktor Warning: No standard errors could be computed.", 
                    "\nThe constraint matrix must have full row-rank ( choose e.g. rows", 
                    paste(rAmat$pivot, collapse = " "), ").",
                    "\nTry to set se = \"boot.model.based\" or \"boot.standard\"."))  
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
  
  if(ncol(Amat) != length(b_unrestr)) {
    stop("Restriktor ERROR: the columns of constraints does not", 
         "\nmatch with the number of parameters.")
  }
  
  
  # original model function call
  call_org <- as.list(object$call)
  # which method M or MM estimation
  method <- call_org[["method"]]
  if (is.null(method)) object$call[["method"]] <- "M"
  # check (only tukey's bisquare supported)
  psi <- call_org[["psi"]]
  if (is.null(psi)) {
    if (object$call[["method"]] == "M") {
      stop("Restriktor ERROR: only tukey's bisquare loss function is supported.")
    }
  } else {
    if (psi != "psi.bisquare") {
      stop("Restriktor ERROR: only tukey's bisquare loss function is supported.")
    }
  }
  
  # compute R-squared 
  # acknowledment: code taken from the lmrob() function from the robustbase package
  df_int <- ifelse(attr(object$terms, "intercept"), 1L, 0L)
  y_mean <- if (df_int == 1L) { 
    sum(weights * y) / sum(weights) 
    } else { 0L }
  yMy <- sum(weights * (y - y_mean)^2)
  rMr <- sum(weights * residuals^2)
  # tukey's bi-square correction
  correc <- 1.207617 
  R2_org <- (yMy - rMr) / (yMy + rMr * (correc - 1))
  
  # compute the loglikelihood
  ll_unc <- con_loglik_lm(object)
  
  # ML unconstrained MSE
  #tau2ml <- (tau^2 * so$df[2]) / n
  #invW <- kronecker(solve(tau2ml), t(X) %*% X)
  
  ## compute mixing weights
  is.augmented <- TRUE
  if (Wt != "none") {
    ## compute mixing weights
    if (all(c(Amat) == 0)) { # unrestrikted case
      wt <- c(rep(0L, p), 1)
      is.augmented <- FALSE
    } else if (Wt == "boot") { # compute mixing weights based on simulation
      wt <- con_weights_boot(VCOV     = Sigma,
                             Amat     = Amat, 
                             meq      = meq, 
                             R        = bootWt.R,
                             parallel = parallel,
                             ncpus    = ncpus,
                             cl       = cl,
                             seed     = seed,
                             verbose  = verbose)
      attr(wt, "bootWt.R") <- bootWt.R 
    } else if (Wt == "mvnorm" && (meq < nrow(Amat))) { # compute mixing weights based on mvnorm
      wt <- rev(con_weights(Amat %*% Sigma %*% t(Amat), meq = meq))
    } else if (Wt == "mvnorm" && (meq == nrow(Amat))) { # only equality constraints
      wt <- rep(0L, ncol(Sigma) + 1)
      wt.idx <- ncol(Sigma) - meq + 1
      wt[wt.idx] <- 1
    }
  } else {
    wt <- NA
  }
  attr(wt, "method") <- Wt
  
  timing$wt <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # # check if the constraints are in line with the data    
  if (all(Amat %*% c(b_unrestr) - bvec >= 0 * bvec) & meq == 0) {  
    b_restr <- b_unrestr
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                wt          = wt,
                b_unrestr   = b_unrestr,
                b_restr     = b_unrestr,
                residuals   = residuals,
                wresid      = object$wresid,
                fitted      = fitted(object),
                weights     = weights,  # prior weights
                w           = object$w, # weights used in the IWLS process
                scale       = object$s, 
                psi         = object$psi,
                R2_org      = R2_org,
                R2_reduced  = R2_org,
                df.residual = so$df[2],
                s2_unrestr  = tau, 
                s2_restr    = tau, 
                loglik      = ll_unc, 
                Sigma       = Sigma, # probably not so robust!
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                iact        = 0L,
                converged   = object$converged, 
                #iter        = object$iter,
                bootout     = NULL, 
                control     = control)
  } else {
    # collect all original model arguments and add constraints
    call_org$formula <- call_org$subset <- call_org$na.action <- 
      call_org$model <- call_org$x.ret <- call_org$y.ret <- 
      call_org$contrasts <- call_org$data <- call_org$weights <- 
      call_org$psi <- call_org$X <- call_org$y <- NULL
    
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = Amat, bvec = bvec, meq = meq))
    
    # fit constraint robust liner model
    rfit <- do.call("conRLM_fit", CALL)
    
    timing$conRLM_fit <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    # constrained coefs
    b_restr <- coefficients(rfit)
      names(b_restr) <- names(b_unrestr)
    b_restr[abs(b_restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    fitted <- fitted(rfit)
    residuals <- residuals(rfit)
    
    # psi(resid/scale) these are the weights used for downweighting the cases.
    w <- rfit$w
    
    # compute loglik
    ll_restr <- con_loglik_lm(rfit)
    
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
        tau_restr <- (summary(rfit)$stddev / sqrt(S)) * sqrt(S.new)
      } else {
        tau_restr <- std * sqrt(S)
      }
    } else {
      rdf <- n - p
     # w   <- rfit$psi(rfit$wresid/rfit$scale)
      S   <- sum((rfit$wresid * w)^2) / rdf
      std <- summary(rfit)$stddev / sqrt(S)
      
      if (!is.null(Amat)) {
        rdf <- n - (p - qr(Amat[0:meq,])$rank)
        S.new <- sum((rfit$wresid * w)^2) / rdf
        tau_restr <- (summary(rfit)$stddev / sqrt(S)) * sqrt(S.new)
      } else {
        tau_restr <- std * sqrt(S)
      }
    }
    
      
    #R^2 under the restrikted object
    df_int <- if (attr(object$terms, "intercept")) { 1L } else { 0L }
    resp_mean <- if (df_int == 1L) { sum(weights * y) / sum(weights) } else { 0 }
    yMy <- sum(weights * (y - resp_mean)^2)
    rMr <- sum(weights * residuals^2)
    # tukey's bi-square correction
    correc <- 1.207617 
    R2_reduced <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                wt          = wt,
                b_unrestr   = b_unrestr,
                b_restr     = b_restr,
                residuals   = residuals,
                wresid      = rfit$wresid,
                fitted      = fitted,
                weights     = weights,
                w           = w, 
                scale       = rfit$s,
                R2_org      = R2_org,
                R2_reduced  = R2_reduced,
                df.residual = so$df[2], 
                s2_unrestr  = tau, 
                s2_restr    = tau_restr, 
                loglik      = ll_restr, 
                Sigma       = Sigma,                                             #probably not so robust???
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                iact        = rfit$iact,
                converged   = rfit$converged, 
                bootout     = NULL, 
                control     = control)
  }
  
  OUT$model_org <- object
  OUT$se <- se 
  OUT$information <- 1 / tau^2 * crossprod(X)
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      #  V <- vcovMM(X = X, resid0 = resid0, residuals = residuals, scale = model$s)  
      information.inv <- con_augmented_information(information  = OUT$information,
                                                   is.augmented = is.augmented,
                                                   X            = X, 
                                                   b_unrestr    = b_unrestr, 
                                                   b_restr      = b_restr,
                                                   Amat         = Amat, 
                                                   bvec         = bvec, 
                                                   meq          = meq)
          
      attr(OUT$information, "inverted")  <- information.inv$information
      attr(OUT$information, "augmented") <- information.inv$information.augmented
      
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
                                 Wt          = "none",
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
                                 Wt          = "none",
                                 parallel    = parallel, 
                                 ncpus       = ncpus, 
                                 cl          = cl)
      
      timing$boot_standard <- (proc.time()[3] - start.time)
      start.time <- proc.time()[3]
    }
  }
  
  OUT$timing$total <- (proc.time()[3] - start.time0)
  
  class(OUT) <- c("restriktor", "conRLM", "conLM")
  
  OUT

}
