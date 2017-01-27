conGLM.glm <- function(object, constraints = NULL, se = "standard", 
                       B = 999, rhs = NULL, neq = 0L, 
                       Wt = "mvnorm", bootWt.R = 99999,
                       parallel = "no", ncpus = 1L, cl = NULL, 
                       seed = NULL, control = list(), verbose = FALSE, 
                       debug = FALSE, ...) {
    
  # check class
  if (!(class(object)[1] %in% c("glm"))) {
    stop("Restriktor ERROR: object must be of class glm.")
  }
  
  # timing
  start.time0 <- start.time <- proc.time()[3]; timing <- list()
  
  # store call
  mc <- match.call()
  
  # rename for internal use
  Amat <- constraints
  bvec <- rhs 
  meq  <- neq
  
  # response variable
  y <- as.matrix(object$model[, attr(object$terms, "response")])
  # model matrix
  X <- model.matrix(object)[,,drop = FALSE]
  # variance-covariance matrix
  Sigma <- vcov(object)
  # familiy and link function
  fam <- object$family
  # model summary
  so <- summary(object)
  # dispersion
  dispersion <- so$dispersion
  # prior weigths
  prior.weights <- weights(object, "prior")
  # working weights
  weights <- weights(object, "working")
  # unconstrained estimates
  b_unrestr <- coef(object)
  # number of parameters
  p <- length(coef(object))
  # sample size
  n <- dim(X)[1]
  
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
  } else if (is.null(constraints)) { 
    # no constraints specified - needed for GORIC to include unconstrained model
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
                 paste(rAmat$pivot, collapse = " "), ", or try to set Wt = \"boot\")"))
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
  # parallel housekeeping
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") {
      have_snow <- TRUE
    }  
    if (!have_mc && !have_snow) {
      ncpus <- 1L
    }  
  }
  
  # standard error methods
  if (se == "default") {
    se <- "standard"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (!(se %in% c("none","standard","const","boot.model.based","boot.standard",
                  "HC","HC0","HC1","HC2","HC3","HC4","HC4m","HC5"))) {
    stop("Restriktor ERROR: standard error method ", sQuote(se), " unknown.")
  }
  
  if(ncol(Amat) != length(b_unrestr)) {
    stop("Restriktor ERROR: length coefficients and the number of",
         "\ncolumns constraints-matrix must be identical")
  }
  
  # unrestricted (weighted) log-likelihood
  LL_unc <- con_loglik_glm(object)
  
  is.augmented <- TRUE
  # compute mixing weights
  if (Wt != "none") {
    # unrestrikted case
    if (all(c(Amat) == 0)) { 
      wt <- c(rep(0L, p), 1)
      is.augmented <- FALSE
    } else if (Wt == "boot") { 
      # compute mixing weights based on simulation
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
    } else if (Wt == "mvnorm" && (meq < nrow(Amat))) {
      # compute mixing weights based on mvtnorm
      wt <- rev(con_weights(Amat %*% Sigma %*% t(Amat), meq = meq))
    } else if (Wt == "mvnorm" && (meq == nrow(Amat))) {
      # only equality constraints
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
  
  # check if the constraints are not in line with the data, else skip optimization
  if (all(Amat %*% c(b_unrestr) - bvec >= 0 * bvec) & meq == 0) {
    b_restr <- b_unrestr
    
    OUT <- list(CON               = CON,
                call              = mc,
                timing            = timing,
                parTable          = parTable,
                family            = object$family,
                wt                = wt,
                b_unrestr         = b_unrestr,
                b_restr           = b_unrestr,
                residuals         = residuals(object, "working"), # unweighted residuals
                fitted            = fitted(object),
                prior.weights     = prior.weights, #prior weights
                weights           = weights, #working weights, weights final iteration
                df.residual       = object$df.residual,
                df.residual_null  = object$df.null,
                dispersion        = so$dispersion, 
                loglik            = LL_unc,
                aic               = object$aic,
                deviance_null     = object$null.deviance,
                deviance          = object$deviance,
                Sigma             = Sigma,
                constraints       = Amat, 
                rhs               = bvec, 
                neq               = meq, 
                iact              = 0L, 
                converged         = object$converged,
                bootout           = NULL, 
                control           = control)  
  } else {
    ## compute constrained estimates for glm()
    
    # collect all original model arguments and add constraints
    call_org <- as.list(object$call)
    
    CALL <- c(list(x = X, y = y, weights = prior.weights, 
                   start = coef(object), etastart = call_org[["etastart"]],
                   mustart = call_org[["mustart"]], offset = object$offset,
                   family = fam, control = call_org[["control"]],
                   intercept = attr(object$terms, "intercept"),
                   Amat = Amat, bvec = bvec, meq = meq))
    
    if (is.null(CALL$control)) { CALL$control <- list() }
    
    # fit restricted generalized liner model
    fit_glmc <- do.call("conGLM_fit", CALL)
    
    # restricted estimates
    b_restr <- fit_glmc$coefficients
    b_restr[abs(b_restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    
    b_restr <- as.vector(b_restr)
      names(b_restr) <- names(b_unrestr)
    fitted <- fitted(fit_glmc)
    residuals <- residuals(fit_glmc, "working")
    
    # weights
    weights <- weights(fit_glmc, "working")
    
    timing$optim <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    # compute restricted dispersion
    df.residual <- n - (p - qr(Amat[0:meq,])$rank)
    
    if (fit_glmc$family$family %in% c("poisson", "binomial")) {
      dispersion <- 1
    } else {
      if (any(weights(fit_glmc) == 0)) {
        warning("observations with zero weight not used for calculating dispersion")
      }
      dispersion <- sum((weights * residuals^2)[weights > 0]) / df.residual
    }
    
    # restricted log-likelihood
    LL_restr <- con_loglik_glm(fit_glmc)

    OUT <- list(CON               = CON,
                call              = mc,
                timing            = timing,
                parTable          = parTable,
                family            = fit_glmc$family,
                wt                = wt,
                b_unrestr         = b_unrestr,
                b_restr           = b_restr,
                residuals         = residuals, # unweighted residuals
                fitted            = fitted,
                prior.weights     = prior.weights, #prior weights
                weights           = weights, #working weights, weights final iteration
                df.residual       = df.residual(object),
                df.residual_null  = fit_glmc$df.null,
                dispersion        = dispersion,
                loglik            = LL_restr,
                aic               = fit_glmc$aic,
                deviance_null     = fit_glmc$null.deviance,
                deviance          = fit_glmc$deviance,
                Sigma             = Sigma,
                constraints       = Amat, 
                rhs               = bvec, 
                neq               = meq, 
                iact              = fit_glmc$iact,
                converged         = fit_glmc$converged,
                bootout           = NULL, 
                control           = control)
  }
  
  # original object
  OUT$model_org <- object
  # type standard error
  OUT$se <- se
  # compute standard errors based on the augmented inverted information matrix or
  # based on the standard bootstrap or model.based bootstrap
  W <- diag(weights) # working
  OUT$information <- t(X) %*% W %*% X / dispersion
  
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
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
      
      if (attr(object$terms, "intercept") && any(Amat[,1] == 1)) {
          stop("Restriktor ERROR: no restriktions on intercept possible",
               "\n       for 'se = boot.model.based' bootstrap method.")
      }
      
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
  
  class(OUT) <- c("restriktor", "conGLM")
  
  OUT
}