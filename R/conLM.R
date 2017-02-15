conLM.lm <- function(object, constraints = NULL, se = "standard", 
                     B = 999, rhs = NULL, neq = 0L, mix.weights = "pmvnorm", 
                     mix.bootstrap = 99999L, parallel = "no", ncpus = 1L, cl = NULL, 
                     seed = NULL, control = list(), verbose = FALSE, 
                     debug = FALSE, ...) {
  
  # check class
  if (!(class(object)[1] %in% c("lm"))) {
    stop("Restriktor ERROR: object must be of class lm.")
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
  if (!(mix.weights %in% c("pmvnorm", "boot", "none"))) {
    stop("Restriktor ERROR: ", sQuote(mix.weights), " method unknown. Choose from \"pmvnorm\", \"boot\", or \"none\"")
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
  ## model matrix
  X <- model.matrix(object)[,,drop = FALSE]
  # model summary
  so <- summary(object)
  # unconstrained residual variance (weighted)
  s2.unrestr <- so$sigma^2
  # weigths
  weights <- weights(object)
  # unconstrained estimates
  b.unrestr <- coef(object)
  # ML unconstrained MSE
  #invW <- kronecker(solve(s2.unrestr), t(X) %*% X)
  Sigma <- vcov(object)#solve(invW)
  # number of parameters
  p <- length(coef(object))
  # sample size
  n <- dim(X)[1]
  # compute (weighted) log-likelihood
  ll.unrestr <- con_loglik_lm(object)
  
  if (debug) {
    print(list(loglik.unc = ll.unrestr))
  }
  
  timing$preparation <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # deal with constraints
  if (!is.null(constraints)) {
    restr.OUT <- con_constraints(object, 
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
  
  # compute the reduced row-echelon form of the constraints matrix
  rAmat <- GaussianElimination(t(Amat))
  if (mix.weights == "pmvnorm") {
    if (rAmat$rank < nrow(Amat) && rAmat$rank != 0L) {
      stop(paste("Restriktor ERROR: The constraint matrix must have full row-rank", 
                 "\n  (choose e.g. rows", paste(rAmat$pivot, collapse = " "), ",or try to set mix.weights = \"boot\")"))
    }
  } 
  else if (rAmat$rank < nrow(Amat) && 
             !(se %in% c("none", "boot.model.based", "boot.standard")) && 
             rAmat$rank != 0L) {
      se <- "none"
      warning(paste("Restriktor Warning: No standard errors could be computed. 
                      The constraint matrix must have full row-rank ( choose e.g. rows", 
                    paste(rAmat$pivot, collapse = " "), "), 
                      Try se = \"boot.model.based\" or \"boot.standard\"."))  
  }
  
  timing$constraints <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  if (ncol(Amat) != length(b.unrestr)) {
    stop("Restriktor ERROR: length coefficients and the number of",
         "\n       columns constraints-matrix must be identical")
  }

  is.augmented <- TRUE
  # compute chi-square-bar weights
  if (mix.weights != "none") {
    if (nrow(Amat) == meq) {
      # equality constraints only
      wt <- rep(0L, ncol(Sigma) + 1)
      wt.idx <- ncol(Sigma) - meq + 1
      wt[wt.idx] <- 1
    } else if (all(c(Amat) == 0)) { 
      # unrestricted case
      wt <- c(rep(0L, p), 1)
      is.augmented <- FALSE
    } else if (mix.weights == "boot") { 
      # compute chi-square-bar weights based on Monte Carlo simulation
      wt <- con_weights_boot(VCOV     = Sigma,
                             Amat     = Amat, 
                             meq      = meq, 
                             R        = mix.bootstrap,
                             parallel = parallel,
                             ncpus    = ncpus,
                             cl       = cl,
                             seed     = seed,
                             verbose  = verbose)
      attr(wt, "mix.bootstrap") <- mix.bootstrap 
    } else if (mix.weights == "pmvnorm" && (meq < nrow(Amat))) {
      # compute chi-square-bar weights based on pmvnorm
      wt <- rev(con_weights(Amat %*% Sigma %*% t(Amat), meq = meq))
    } 
  } else {
    wt <- NA
  }
  attr(wt, "method") <- mix.weights
  
  if (debug) {
    print(list(wt = wt))
  }
  
  timing$wt <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # check if the constraints are not in line with the data, else skip optimization
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) & meq == 0) {
    b.restr  <- b.unrestr
    s2.restr <- s2.unrestr
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.unrestr,
                residuals   = object$residuals, # unweighted residuals
                fitted      = object$fitted,
                weights     = weights,
                df.residual = object$df.residual,
                R2.org      = so$r.squared, 
                R2.reduced  = so$r.squared,
                s2.unrestr  = s2.unrestr, 
                s2.restr    = s2.unrestr, 
                loglik      = ll.unrestr, 
                Sigma       = Sigma,
                constraints = Amat, 
                wt          = wt,
                rhs         = bvec, 
                neq         = meq, 
                iact        = 0L, 
                bootout     = NULL, 
                control     = control)  
  } else {
    # compute constrained estimates using quadprog
    out.QP <- con_solver_lm(X         = X, 
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
    
    # lm
    if (ncol(y) == 1L) {
      b.restr <- as.vector(b.restr)
        names(b.restr) <- names(b.unrestr)
      fitted <- X %*% b.restr
      residuals <- y - fitted
      
      # compute log-likelihood
      object.restr <- list(residuals = residuals, weights = weights)
      ll.restr <- con_loglik_lm(object.restr)
      
      if (debug) {
        print(list(loglik.restr = ll.restr))
      }  
      
      # compute R^2
      if (is.null(weights)) {
        mss <- if (attr(object$terms, "intercept")) {
          sum((fitted - mean(fitted))^2)
        } else { sum(fitted^2) }
        rss <- sum(residuals^2)
      } else {
        mss <- if (attr(object$terms, "intercept")) {
          m <- sum(weights * fitted / sum(weights))
          sum(weights * (fitted - m)^2)
        } else { sum(weights * fitted^2) }
        rss <- sum(weights * residuals^2)
      }
      R2.reduced <- mss / (mss + rss)
      
      # compute residual degreees of freedom, corrected for equality constraints.
      df.residual <- n - (p - qr(Amat[0:meq,])$rank)
      # compute weighted residuals
      if (is.null(weights)) {
        s2.restr <- sum(residuals^2) / df.residual  
      } else {
        s2.restr <- sum(weights * residuals^2) / df.residual
      }
    } else { #mlm <FIXME>
      stop("mlm not implemented (yet)")
    }

    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.restr,
                residuals   = residuals, # unweighted residuals
                fitted      = fitted,
                weights     = weights,
                df.residual = object$df.residual,
                R2.org      = so$r.squared, 
                R2.reduced  = R2.reduced,
                s2.unrestr  = s2.unrestr, 
                s2.restr    = s2.restr, 
                loglik      = ll.restr, 
                Sigma       = Sigma,
                constraints = Amat, 
                wt          = wt,
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
  OUT$information <- 1/s2.restr * crossprod(X)
  # compute standard errors based on the augmented inverted information matrix or
  # based on the standard bootstrap or model.based bootstrap
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      information.inv <- con_augmented_information(information  = OUT$information,
                                                   is.augmented = is.augmented,
                                                   X            = X, 
                                                   b.unrestr    = b.unrestr, 
                                                   b.restr      = b.restr,
                                                   Amat         = Amat, 
                                                   bvec         = bvec, 
                                                   meq          = meq) 
          
      attr(OUT$information, "inverted")  <- information.inv$information
      attr(OUT$information, "augmented") <- information.inv$information.augmented
      if (debug) {
        print(list(information = OUT$information))
      }
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
                                 mix.weights          = "none",
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
                                 constraints = Amat,
                                 rhs         = bvec, 
                                 neq         = meq, 
                                 se          = "none",
                                 mix.weights          = "none",
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
  
  OUT$timing$total <- (proc.time()[3] - start.time0)
  
  if (ncol(y) == 1L) {
    class(OUT) <- c("restriktor", "conLM")
  } else if (ncol(y) > 1L) {
    class(OUT) <- c("restriktor", "conMLM")
  }
    
  OUT
}