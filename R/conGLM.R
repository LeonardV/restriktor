conGLM.glm <- function(object, constraints = NULL, se = "standard", 
                       B = 999L, rhs = NULL, neq = 0L, mix_weights = "pmvnorm", 
                       parallel = "no", ncpus = 1L, cl = NULL, seed = NULL, 
                       control = list(), verbose = FALSE, debug = FALSE, ...) {
    
  # check class
  if (!(class(object)[1] == "glm")) {
    stop("Restriktor ERROR: object must be of class glm.")
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
  # check method to compute chi-square-bar weights
  if (!(mix_weights %in% c("pmvnorm", "boot", "none"))) {
    stop("Restriktor ERROR: ", sQuote(mix_weights), " method unknown. Choose from \"pmvnorm\", \"boot\", or \"none\"")
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
  b.unrestr <- coef(object)
  b.unrestr[abs(b.unrestr) < ifelse(is.null(control$tol), 
                                    sqrt(.Machine$double.eps), 
                                    control$tol)] <- 0L
  # number of parameters
  p <- length(coef(object))
  # sample size
  n <- dim(X)[1]
  # unrestricted (weighted) log-likelihood
  ll.unrestr <- con_loglik_glm(object)
  
  if (debug) {
    print(list(loglik.unrestr = ll.unrestr))
  }  
  
  timing$preparation <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # deal with constraints
  if (!is.null(constraints)) {
    restr.OUT <- con_constraints(object, 
                                 VCOV        = Sigma,
                                 est         = b.unrestr,
                                 constraints = Amat, 
                                 bvec        = bvec, 
                                 meq         = meq, 
                                 debug       = debug)  
    # a list with useful information about the restriktions.}
    CON <- restr.OUT$CON
    # a parameter table with information about the observed variables in the object 
    # and the imposed restriktions.}
    parTable <- restr.OUT$parTable
    # constraints matrix
    Amat <- restr.OUT$Amat
    # rhs 
    bvec <- restr.OUT$bvec
    # neq
    meq <- restr.OUT$meq
  } else if (is.null(constraints)) { 
    # no constraints specified - needed for GORIC to include unconstrained model
    CON <- NULL
    parTable <- NULL
    Amat <- rbind(rep(0L, p))
    bvec <- rep(0L, nrow(Amat))
    meq  <- 0L
  } 
  
  if (length(Amat) == 0L) {
    Amat <- rbind(rep(0L, p))
    bvec <- rep(0L, nrow(Amat))
    meq  <- 0L
  }
  
  ## create list for warning messages
  messages <- list()
  
  ## check if constraint matrix is of full-row rank. 
  # rAmat <- GaussianElimination(t(Amat))
  Amat_meq_PT <- PT_Amat_meq(Amat, meq)
  rAmat <- Amat_meq_PT$RREF

  if (rAmat$rank < nrow(Amat) &&
             !(se %in% c("none", "boot.model.based", "boot.standard")) &&
             rAmat$rank != 0L) {
    se <- "none"
    warning(paste("\nRestriktor Warning: No standard errors could be computed, because",
                  "the constraint matrix must be full row-rank.",
                  "Try se = \"boot.model.based\" or \"boot.standard\"."))
  }
  
  
  ## some checks
  if(ncol(Amat) != length(b.unrestr)) {
    stop(paste("Restriktor ERROR: length coefficients and the number of",
         "columns constraints-matrix must be identical"))
  }
  
  if (!(nrow(Amat) == length(bvec))) {
    stop("nrow(Amat) != length(bvec)")
  }
  
  
  timing$constraints <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # check if the constraints are not in line with the data, else skip optimization
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) && meq == 0) {
    b.restr <- b.unrestr
    
    OUT <- list(CON               = CON,
                call              = mc,
                timing            = timing,
                parTable          = parTable,
                family            = object$family,
                b.unrestr         = b.unrestr,
                b.restr           = b.unrestr,
                residuals         = residuals(object, "working"), # unweighted residuals
                fitted            = object$fitted.values,
                prior.weights     = prior.weights, #prior weights
                weights           = weights, #working weights, weights final iteration
                df.residual       = object$df.residual,
                df.residual.null  = object$df.null,
                dispersion        = so$dispersion, 
                loglik            = ll.unrestr,
                aic               = object$aic,
                deviance.null     = object$null.deviance,
                deviance          = object$deviance,
                Sigma             = Sigma,
                constraints       = Amat, 
                rhs               = bvec, 
                neq               = meq, 
                wt.bar            = NULL,
                iact              = 0L, 
                converged         = object$converged,
                iter              = object$iter,
                bootout           = NULL, 
                control           = control)  
  } else {
    ## compute constrained estimates for glm()
    # collect all original model arguments and add constraints
    call.org <- as.list(object$call)
    
    CALL <- c(list(x = X, y = y, weights = prior.weights, 
                   start = coef(object), etastart = call.org[["etastart"]],
                   mustart = call.org[["mustart"]], offset = object$offset,
                   family = fam, control = call.org[["control"]],
                   intercept = attr(object$terms, "intercept"),
                   Amat = Amat, bvec = bvec, meq = meq))
    
    if (is.null(CALL$control)) { CALL$control <- list() }
    
    # fit restricted generalized liner model
    fit.glmc <- do.call("conGLM_fit", CALL)
    
    timing$optim <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    # restricted estimates
    b.restr <- fit.glmc$coefficients
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    b.restr <- as.vector(b.restr)
      names(b.restr) <- names(b.unrestr)
    fitted <- fit.glmc$fitted
    residuals <- residuals(fit.glmc, "working")
    # weights
    weights <- weights(fit.glmc, "working")
    # compute restricted dispersion
    df.residual <- n - (p - qr(Amat[0:meq,])$rank)
    
    if (fit.glmc$family$family %in% c("poisson", "binomial")) {
      dispersion <- 1
    } else {
      if (any(weights(fit.glmc) == 0)) {
        warning("observations with zero weight not used for calculating dispersion")
      }
      dispersion <- sum((weights * residuals^2)[weights > 0]) / df.residual
    }
    
    # restricted log-likelihood
    ll.restr <- con_loglik_glm(fit.glmc)

    if (debug) {
      print(list(loglik.restr = ll.restr))
    }
    
    OUT <- list(CON               = CON,
                call              = mc,
                timing            = timing,
                parTable          = parTable,
                family            = fit.glmc$family,
                b.unrestr         = b.unrestr,
                b.restr           = b.restr,
                residuals         = residuals, # unweighted residuals
                fitted            = fitted,
                prior.weights     = prior.weights, #prior weights
                weights           = weights, #working weights, weights final iteration
                df.residual       = df.residual(object),
                df.residual.null  = fit.glmc$df.null,
                dispersion        = dispersion,
                loglik            = ll.restr,
                aic               = fit.glmc$aic,
                deviance.null     = fit.glmc$null.deviance,
                deviance          = fit.glmc$deviance,
                Sigma             = Sigma,
                constraints       = Amat, 
                rhs               = bvec, 
                neq               = meq, 
                wt.bar            = NULL,
                iact              = fit.glmc$iact,
                converged         = fit.glmc$converged,
                iter              = fit.glmc$iter,
                bootout           = NULL, 
                control           = control)
  }
  
  # original object
  OUT$model.org <- object
  # type standard error
  OUT$se <- se
  # compute standard errors based on the augmented inverted information matrix or
  # based on the standard bootstrap or model.based bootstrap
  W <- diag(weights) # working
  OUT$information <- t(X) %*% W %*% X / dispersion
  
  if (se != "none") {
    is.augmented <- TRUE
    if (all(c(Amat) == 0)) { 
      # unrestricted case
      is.augmented <- FALSE
    } 
    
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      information.inv <- try(con_augmented_information(information  = OUT$information,
                                                       is.augmented = is.augmented,
                                                       X            = X, 
                                                       b.unrestr    = b.unrestr, 
                                                       b.restr      = b.restr,
                                                       Amat         = Amat, 
                                                       bvec         = bvec, 
                                                       meq          = meq), silent = TRUE)
      
      if (inherits(information.inv, "try-error")) {
        stop(paste("Restriktor Warning: No standard errors could be computed.",
                   "Try to set se = \"none\", \"boot.model.based\" or \"boot.standard\"."))
      }
      
      
      attr(OUT$information, "inverted")  <- information.inv$information
      attr(OUT$information, "augmented") <- information.inv$information.augmented
      
      if (debug) {
        print(list(information = OUT$information))
      }
    } else if (se == "boot.model.based") {
      if (attr(object$terms, "intercept") && any(Amat[,1] == 1)) {
          stop(paste("Restriktor ERROR: no restrictions on intercept possible",
               "for 'se = boot.model.based' bootstrap method."))
      }
      OUT$bootout <- con_boot_lm(object      = object, 
                                 B           = B, 
                                 fixed       = TRUE, 
                                 Amat        = Amat,
                                 bvec        = bvec, 
                                 meq         = meq, 
                                 se          = "none",
                                 mix_weights = "none",
                                 parallel    = parallel, 
                                 ncpus       = ncpus, 
                                 cl          = cl)
      if (debug) {
        print(list(bootout = OUT$bootout))
      }
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(object      = object, 
                                 B           = B, 
                                 fixed       = FALSE, 
                                 Amat        = Amat,
                                 bvec        = bvec, 
                                 meq         = meq, 
                                 se          = "none",
                                 mix_weights = "none",
                                 parallel    = parallel, 
                                 ncpus       = ncpus, 
                                 cl          = cl)
      if (debug) {
        print(list(bootout = OUT$bootout))
      }
    }
    timing$standard.error <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
  }
  
  start.time <- proc.time()[3]
  
  ## determine level probabilies
  if (mix_weights != "none" && inherits(object, "goric")) {
    RREF <- Amat_meq_PT$RREF
    OUT$PT_Amat <- PT_Amat <- Amat_meq_PT$PT_Amat
    OUT$PT_meq <- PT_meq <- Amat_meq_PT$PT_meq
    
    if (mix_weights == "pmvnorm") {
      if (RREF$rank < nrow(PT_Amat) && RREF$rank != 0L) {
        messages$mix_weights_rank <- paste(
          "Restriktor message: Since the constraint matrix is not full row-rank, the level probabilities", 
          "are calculated using mix_weights = \"boot\" (the default is mix_weights = \"pmvnorm\").",
          "For more information see ?restriktor.\n"
        )
        mix_weights <- "boot"
      }
    } 
    
    wt.bar <- calculate_weight_bar(Amat = PT_Amat, meq = PT_meq, VCOV = Sigma, 
                                   mix_weights = mix_weights, seed = seed, 
                                   control = control, verbose = verbose, ...) 
  } else {
    if (mix_weights == "pmvnorm") {
      if (rAmat$rank < nrow(Amat) && rAmat$rank != 0L) {
        messages$mix_weights_rank <- paste(
          "Restriktor message: Since the constraint matrix is not full row-rank, the level probabilities", 
          "are calculated using mix_weights = \"boot\" (the default is mix_weights = \"pmvnorm\").",
          "For more information see ?restriktor.\n"
        )
        mix_weights <- "boot"
      }
    }
    wt.bar <- calculate_weight_bar(Amat = Amat, meq = meq, VCOV = Sigma, 
                                   mix_weights = mix_weights, seed = seed, 
                                   control = control, verbose = verbose, ...) 
  }
  attr(wt.bar, "method") <- mix_weights
  OUT$wt.bar <- wt.bar
  
  if (debug) {
    print(list(mix.weigths = wt.bar))
  }
  
  timing$mix_weights <- (proc.time()[3] - start.time)
  OUT$messages <- messages
  OUT$timing$total <- (proc.time()[3] - start.time0)
  
  class(OUT) <- c("restriktor", "conGLM")
  
  OUT
}