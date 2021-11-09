conMLM.mlm <- function(object, constraints = NULL, se = "none", 
                       B = 999, rhs = NULL, neq = 0L, mix.weights = "pmvnorm", 
                       mix.bootstrap = 99999L, parallel = "no", ncpus = 1L, cl = NULL, 
                       seed = NULL, control = list(), verbose = FALSE, 
                       debug = FALSE, ...) {
  
  # check class
  if (!(class(object)[1] == "mlm")) {
    stop("Restriktor ERROR: object must be of class mlm.")
  }
  
  # not available yet.
  se <- "none"
  
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
  # sample size
  n <- dim(X)[1]
  # number of parameters
  p <- length(coef(object))
  # unconstrained estimates
  b.unrestr <- coef(object)
  b.unrestr[abs(b.unrestr) < ifelse(is.null(control$tol), 
                                    sqrt(.Machine$double.eps), 
                                    control$tol)] <- 0L
  # unconstrained residual variance (weighted)
  residuals <- object$residuals
  # residual degree of freedom
  df.residual <- object$df.residual
  s2 <- (t(residuals) %*% residuals) / df.residual #sigma(object)
  # vcov
  Sigma <- vcov(object)
  # model summary
  so <- summary(object)
  # original R2
  R2.org <- unlist(lapply(so, function(x) x$r.squared))
  # compute log-likelihood
  object.restr <- list(residuals = residuals)
  ll.unrestr <- con_loglik_lm(object.restr)
  
  if (debug) {
    print(list(loglik.unc = ll.unrestr))
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
  
  # if only new parameters are defined and no constraints
  if (length(Amat) == 0L) {
    Amat <- rbind(rep(0L, p))
    bvec <- rep(0L, nrow(Amat))
    meq  <- 0L
  }
  
  ## create list for warning messages
  messages <- list()
  
  ## check if constraint matrix is of full-row rank. 
  rAmat <- GaussianElimination(t(Amat))
  if (mix.weights == "pmvnorm") {
    if (rAmat$rank < nrow(Amat) && rAmat$rank != 0L) {
      messages$mix_weights <- paste(
        "Restriktor message: Since the constraint matrix is not full row-rank, the level probabilities 
 are calculated using mix.weights = \"boot\" (the default is mix.weights = \"pmvnorm\").
 For more information see ?restriktor.\n"
      )
      mix.weights <- "boot"
    }
  } else if (rAmat$rank < nrow(Amat) &&
             !(se %in% c("none", "boot.model.based", "boot.standard")) &&
             rAmat$rank != 0L) {
    se <- "none"
    warning(paste("\nRestriktor Warning: No standard errors could be computed.
                    The constraint matrix must be full row-rank.
                    Try se = \"boot.model.based\" or \"boot.standard\"."))
  }
  

  ## some checks
  if (ncol(Amat) != length(b.unrestr)) {
    stop("Restriktor ERROR: length coefficients and the number of",
         "\n       columns constraints-matrix must be identical")
  }
  
  if (!(nrow(Amat) == length(bvec))) {
    stop("nrow(Amat) != length(bvec)")
  }
  
  timing$constraints <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # check if the constraints are not in line with the data, else skip optimization
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) && meq == 0) {
    b.restr  <- b.unrestr
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.unrestr,
                residuals   = residuals, 
                fitted      = object$fitted.values,
                df.residual = df.residual,
                R2.org      = R2.org, 
                R2.reduced  = R2.org,
                s2          = s2, 
                loglik      = ll.unrestr, 
                Sigma       = Sigma,
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                wt.bar      = NULL,
                iact        = 0L,
                Niter       = 2L,
                bootout     = NULL, 
                control     = control)  
  } else {
    # compute constrained estimates using quadprog
    out.solver <- con_solver_lm(X         = X, 
                                y         = y, 
                                w         = NULL, 
                                Amat      = Amat,
                                bvec      = bvec, 
                                meq       = meq, 
                                absval    = ifelse(is.null(control$absval), 
                                                   sqrt(.Machine$double.eps), 
                                                   control$absval),
                                maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                   control$maxit))
    out.QP <- out.solver$qp
    b.restr <- matrix(out.QP$solution, ncol = ncol(y))
      colnames(b.restr) <- colnames(b.unrestr)
      rownames(b.restr) <- rownames(b.unrestr)
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    # number of iterations
    Niter <- out.solver$Niter
    
    timing$optim <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
     
    # compute log-likelihood
    fitted <- X %*% b.restr
    residuals <- y - fitted
    s2 <- out.solver$s2
    object.restr <- list(residuals = residuals)
    ll.restr <- con_loglik_lm(object.restr)
    
    if (debug) {
      print(list(loglik.restr = ll.restr))
    }  
    
    # compute R^2
    mss <- if (attr(object$terms, "intercept")) {
      colSums((fitted - colMeans(fitted))^2)
    } else { colSums(fitted^2) }
    rss <- colSums(residuals^2)
    R2.reduced <- mss / (mss + rss)
    
    # compute residual degreees of freedom, corrected for equality constraints.
    df.residual <- df.residual + qr(Amat[0:meq,])$rank
   
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.restr,
                residuals   = residuals, # unweighted residuals
                fitted      = fitted,
                df.residual = object$df.residual,
                R2.org      = R2.org, 
                R2.reduced  = R2.reduced,
                s2          = s2,
                loglik      = ll.restr, 
                Sigma       = Sigma,
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                wt.bar      = NULL,
                iact        = out.QP$iact, 
                Niter       = Niter,
                bootout     = NULL, 
                control     = control)
  }
  
  # original object
  OUT$model.org <- object
  # # type standard error
  OUT$se <- se
  # OUT$information <- 1/s2 * crossprod(X)
  OUT$information <- NULL
  # # compute standard errors based on the augmented inverted information matrix or
  # # based on the standard bootstrap or model.based bootstrap
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
        stop(paste("Restriktor Warning: No standard errors could be computed.
                      Try to set se = \"none\", \"boot.model.based\" or \"boot.standard\"."))
      }

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
                                 Amat        = Amat,
                                 bvec        = bvec,
                                 meq         = meq,
                                 se          = "none",
                                 mix.weights = "none",
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
                                 mix.weights = "none",
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
  
  ## determine level probabilities
  if (mix.weights != "none") {
    if (nrow(Amat) == meq) {
      # equality constraints only
      wt.bar <- rep(0L, ncol(Sigma) + 1)
      wt.bar.idx <- ncol(Sigma) - qr(Amat)$rank + 1
      wt.bar[wt.bar.idx] <- 1
    } else if (all(c(Amat) == 0)) { 
      # unrestricted case
      wt.bar <- c(rep(0L, p), 1)
    } else if (mix.weights == "boot") { 
      # compute chi-square-bar weights based on Monte Carlo simulation
      wt.bar <- con_weights_boot(VCOV     = Sigma,
                                 Amat     = Amat, 
                                 meq      = meq, 
                                 R        = mix.bootstrap,
                                 parallel = parallel, 
                                 ncpus    = ncpus, 
                                 cl       = cl,
                                 seed     = seed,
                                 verbose  = verbose)
      attr(wt.bar, "mix.bootstrap") <- mix.bootstrap 
    } else if (mix.weights == "pmvnorm" && (meq < nrow(Amat))) {
      # compute chi-square-bar weights based on pmvnorm
      wt.bar <- rev(con_weights(Amat %*% Sigma %*% t(Amat), meq = meq))
    } 
  } else {
    wt.bar <- NA
  }
  attr(wt.bar, "method") <- mix.weights
  OUT$wt.bar <- wt.bar
  
  if (debug) {
    print(list(mix.weights = wt.bar))
  }
  
  timing$mix.weights <- (proc.time()[3] - start.time)
  OUT$messages <- messages
  OUT$timing$total <- (proc.time()[3] - start.time0)
  
  class(OUT) <- c("restriktor", "conMLM")
  
  OUT
}