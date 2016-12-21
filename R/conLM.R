conLM.lm <- function(object, constraints = NULL, se = "standard", 
                     B = 999, rhs = NULL, neq = 0L, Wt = TRUE, 
                     bootWt = FALSE, bootWt.R = 99999,
                     parallel = "no", ncpus = 1L, cl = NULL, 
                     seed = NULL, control = NULL, verbose = FALSE, 
                     debug = FALSE, ...) {
  
  # check class
  if (!(class(object)[1] %in% c("lm"))) {
    stop("Restriktor ERROR: object must be of class lm.")
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
  # model summary
  so <- summary(object)
  # unconstrained residual variance (weighted)
  s2.unc <- so$sigma^2
  # weigths
  weights <- weights(object)
  # unconstrained estimates
  b.unrestr <- coef(object)
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
  if (attr(object$terms, "intercept")) {
    if (se == "boot.model.based" & any(Amat[,1] == 1)) { 
      stop("Restriktor ERROR: no restriktions on intercept possible for 'se = boot.model.based' bootstrap method.")
    }
  }
  if(ncol(Amat) != length(b.unrestr)) {
    stop("Restriktor ERROR: length coefficients and the number of columns constraints-matrix must be identical")
  }

  # compute (weighted) log-likelihood
  ll.unc <- con_loglik_lm(X = X, 
                          y = y, 
                          b = b.unrestr, 
                          w = weights)
  LL.unc <- ll.unc$loglik
  
  timing$LLik <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  # ML unconstrained MSE
  #s2ml.unc <- s2.unc * object$df.residual / n
  #invW <- kronecker(solve(s2.unc), t(X) %*% X)
  W <- vcov(object)#solve(invW)
  
  is.augmented <- TRUE
  # compute mixing weights
  if (Wt) {
    ## compute mixing weights
    if (all(c(Amat) == 0)) { # unrestrikted case
      wt <- c(rep(0L, p), 1)
      is.augmented <- FALSE
    } else if (bootWt) { # compute mixing weights based on simulation
      wt <- con_weights_boot(VCOV     = W,
                             Amat     = Amat, 
                             meq      = meq, 
                             R        = bootWt.R,
                             parallel = parallel,
                             ncpus    = ncpus,
                             cl       = cl,
                             seed     = seed,
                             verbose  = verbose)
      # compute mixing weights based on mvtnorm
    } else if (!bootWt & (meq < nrow(Amat))) {
      wt <- rev(con_weights(Amat %*% W %*% t(Amat), meq = meq))
    # only equality constraints
    } else if (!bootWt & (meq == nrow(Amat))) { 
      wt <- rep(0L, ncol(W) + 1)
      wt.idx <- ncol(W) - meq + 1
      wt[wt.idx] <- 1
    }
    attr(wt, "bootWt") <- bootWt
    if (bootWt) { attr(wt, "bootWt.R") <- bootWt.R }
  } else {
    wt <- NULL
  }
  
  timing$wt <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # check if the constraints are not in line with the data, else skip optimization
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) & meq == 0) {
    b.restr <- b.unrestr
    s2.restr <- s2.unc
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                wt          = wt,
                b.unrestr   = b.unrestr,
                b.restr     = b.unrestr,
                residuals   = object$residuals, # unweighted residuals
                fitted      = object$fitted,
                weights     = object$weights,
                df.residual = object$df.residual,
                R2.org      = so$r.squared, 
                R2.reduced  = so$r.squared,
                s2.unc      = s2.unc, 
                s2.restr    = s2.unc, 
                loglik      = LL.unc, 
                Sigma       = W,
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                iact        = 0L, 
                bootout     = NULL, 
                control     = control)  
  } else {
    # compute constrained estimates for lm()
    out.QP <- con_solver(X         = X, 
                         y         = y, 
                         b.unrestr = b.unrestr, 
                         w         = weights, 
                         Amat      = Amat,
                         bvec      = bvec, 
                         meq       = meq, 
                         absval    = ifelse(is.null(control$absval), 
                                         sqrt(.Machine$double.eps), 
                                         control$absval),
                         maxit     = ifelse(is.null(control$maxit), 1e04, 
                                         control$maxit))
    b.restr <- matrix(out.QP$solution, ncol = ncol(y))
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    
    timing$optim <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    ll.restr <- con_loglik_lm(X = X, 
                              y = y, 
                              b = b.restr, 
                              w = weights)
    LL.restr <- ll.restr$loglik
    
    timing$conLLik <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    # lm
    if (ncol(y) == 1L) {
      b.restr <- as.vector(b.restr)
        names(b.restr) <- names(b.unrestr)
      fitted <- X %*% b.restr
      residuals <- y - fitted
      
      # compute R^2
      if (is.null(weights)) {
        mss <- if (attr(object$terms, "intercept")) {
          sum((fitted - mean(fitted))^2)
        } else {
          sum(fitted^2)
        }
        rss <- sum(residuals^2)
      } else {
        mss <- if (attr(object$terms, "intercept")) {
          m <- sum(weights * fitted / sum(weights))
          sum(weights * (fitted - m)^2)
        } else {
          sum(weights * fitted^2)
        }
        rss <- sum(weights * residuals^2)
      }
      R2.reduced <- mss / (mss + rss)
      
      timing$R2 <- (proc.time()[3] - start.time)
      start.time <- proc.time()[3]
      
      # compute residual degreees of freedom, corrected for equality constraints.
      # only for the estimation, not for testing. 
      df.residual <- n - (p - qr(Amat[0:meq,])$rank)
      # compute weighted residuals
      if (is.null(weights)) {
        s2.restr <- sum(residuals^2) / df.residual  
      } else {
        s2.restr <- sum(weights * residuals^2) / df.residual
      }
    } else { #mlm <FIXME>
      stop("mlm not implemented (yet)")
      # residuals <- y - fitted
      #   rownames(b.restr) <- rownames(b.unrestr)
      #   colnames(b.restr) <- colnames(b.unrestr)
      # fitted <- X %*% b.restr
      # df <- object$df.residual
      # df.residual <- df + qr(Amat[0:meq,])$rank
      # se <- "none"
    }

    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                wt          = wt,
                b.unrestr   = b.unrestr,
                b.restr     = b.restr,
                residuals   = residuals, # unweighted residuals
                fitted      = fitted,
                weights     = weights,
                df.residual = object$df.residual,
                R2.org      = so$r.squared, 
                R2.reduced  = R2.reduced,
                s2.unc      = s2.unc, 
                s2.restr    = s2.restr, 
                loglik      = LL.restr, 
                Sigma       = W,
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                iact        = out.QP$iact, 
                bootout     = NULL, 
                control     = control)
  }
  
  # original object
  OUT$model.org <- object
  # type standard error
  OUT$se <- se
  # compute standard errors based on the augmented inverted information matrix or
  # based on the standard bootstrap or model.based bootstrap
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      OUT$information <- 1/s2.restr * crossprod(X)
      
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
  
  if (ncol(y) == 1L) {
    class(OUT) <- c("restriktor", "conLM", "lm")
  } else if (ncol(y) > 1L) {
    class(OUT) <- c("restriktor", "conMLM", "mlm")
  }
    
    OUT
}