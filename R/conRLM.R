#compute restrikted robust estimates
conRLM.rlm <- function(object, constraints = NULL, se = "standard", 
                       B = 999, rhs = NULL, neq = 0L, mix.weights = "pmvnorm", 
                       mix.bootstrap = 99999L, parallel = "no", ncpus = 1L, cl = NULL, 
                       seed = NULL, control = list(), verbose = FALSE, 
                       debug = FALSE, ...) { 
  
  # check class
  if (!(class(object)[1] == "rlm")) {
    stop("Restriktor ERROR: object must be of class rlm.")
  }
  # standard error methods
  if (se == "default") {
    se <- "standard"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (!(se %in% c("none","standard","const","boot.model.based","boot.standard","HC","HC0",
                  "HC1","HC2","HC3","HC4","HC4m","HC5"))) {
    stop("Restriktor ERROR: standard error method ", sQuote(se), " unknown.")
  }
  # check method to compute chi-square-bar weights
  if (!(mix.weights %in% c("pmvnorm", "boot", "none"))) {
    stop("Restriktor ERROR: ", sQuote(mix.weights), " method unknown. Choose from \"pmvnorm\", \"boot\", or \"none\"")
  }
  
  # original model function call
  call.org <- as.list(object$call)
  # M or MM estimation?
  method <- call.org[["method"]]
  # loss function
  psi <- call.org[["psi"]]
  if (is.null(method)) { 
    object$call[["method"]] <- "M"
  }
  # check (only tukey's bisquare supported)
  if (is.null(psi)) {
    if (object$call[["method"]] == "M") {
      stop("Restriktor ERROR: only tukey's bisquare loss function is supported.")
    }
  } else {
    if (psi != "psi.bisquare") {
      stop("Restriktor ERROR: only tukey's bisquare loss function is supported.")
    }
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
  b.unrestr <- coef(object)
  b.unrestr[abs(b.unrestr) < ifelse(is.null(control$tol), 
                                    sqrt(.Machine$double.eps), 
                                    control$tol)] <- 0L
  # unrestrikted estimate of scale
  scale <- object$s
  # vcov(object) incorrect! The rlm() function uses a robust sum of squares. 
  ## Yohai (1987, p. 648). Eq. 4.2, 4.3 and 4.4.
  # (for homoscedastic regression)
  cc <- ifelse(is.null(call.org$c), 4.685061, call.org$c)
  res <- y - X %*% b.unrestr
  rstar <- res / scale
  a <- mean(tukeyChi(rstar, cc, deriv = 1)^2)
  b <- mean(tukeyChi(rstar, cc, deriv = 2))
  tau2 <- scale^2 * a/b^2
  Sigma <- tau2 * solve(crossprod(X))
  # a scale estimate used for the standard errors
  stddev <- sqrt(tau2) #so$stddev
  # residuals
  residuals <- residuals(object) # NOT working residual
  # sample size
  n <- dim(X)[1]
  # number of parameters
  p <- dim(X)[2]
  # compute the loglikelihood
  ll.unrestr <- con_loglik_lm(object)
  
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
    # a list with useful information about the restrictions.}
    CON <- restr.OUT$CON
    # a parameter table with information about the observed variables in the object 
    # and the imposed restrictions.}
    parTable <- restr.OUT$parTable
    # constraints matrix
    Amat <- restr.OUT$Amat
    # rhs 
    bvec <- restr.OUT$bvec
    # neq
    meq <- restr.OUT$meq
  } else if (is.null(constraints)) { # no constraints specified - needed for GORIC to include unconstrained model
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
  if(ncol(Amat) != length(b.unrestr)) {
    stop("Restriktor ERROR: the columns of constraints does not", 
         "\n       match with the number of parameters.")
  }
  
  if (!(nrow(Amat) == length(bvec))) {
    stop("nrow(Amat) != length(bvec)")
  }
  
  timing$constraints <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # compute R-squared 
  # acknowledgement: code taken from the lmrob() function from the robustbase package
  wgt <- object$w
  df.int <- ifelse(attr(object$terms, "intercept"), 1L, 0L)
  y.mean <- if (df.int == 1L) { 
    sum(wgt * y) / sum(wgt) 
    } else { 0L }
  yMy <- sum(wgt * (y - y.mean)^2)
  rMr <- sum(wgt * residuals^2)
  # tukey's bi-square correction
  correc <- 1.207617 
  R2.org <- (yMy - rMr) / (yMy + rMr * (correc - 1))
  
  start.time <- proc.time()[3]
  
  # # check if the constraints are in line with the data    
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) & meq == 0) {  
    b.restr   <- b.unrestr
    #scale.restr <- scale
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.unrestr,
                residuals   = residuals,
                wresid      = object$wresid,
                fitted      = object$fitted.values,
                weights     = weights,  # prior weights
                wgt         = object$w, # weights used in the IWLS process
                scale       = object$s, 
                stddev      = stddev,
                psi         = object$psi,
                R2.org      = R2.org,
                R2.reduced  = R2.org,
                df.residual = so$df[2],
                loglik      = ll.unrestr, 
                Sigma       = Sigma, 
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                wt.bar      = NULL,
                iact        = 0L,
                converged   = object$converged, 
                iter        = object$conv,
                bootout     = NULL, 
                control     = control)
  } else {
    # collect all original model arguments and add constraints
    call.org$formula <- call.org$subset <- call.org$na.action <- 
      call.org$model <- call.org$x.ret <- call.org$y.ret <- 
      call.org$contrasts <- call.org$data <- call.org$weights <- 
      call.org$psi <- call.org$X <- call.org$y <- NULL
    
    CALL <- c(call.org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = Amat, bvec = bvec, meq = meq))
    
    # fit constrained robust liner model
    rfit <- do.call("conRLM_fit", CALL)
    
    timing$conRLM.fit <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    # constrained coefs
    b.restr <- coefficients(rfit)
      names(b.restr) <- names(b.unrestr)
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    fitted <- rfit$fitted
    residuals <- rfit$residuals
    # psi(resid/scale) these are the weights used for downweighting the cases.
    wgt <- rfit$w
    # compute loglik
    ll.restr <- con_loglik_lm(rfit)
    
    if (debug) {
      print(list(loglik.restr = ll.restr))
    }
    
    # compute correct stddev
    scale <- rfit$s
    res <- y - X %*% b.restr
    rstar <- res / scale
    a <- mean(tukeyChi(rstar, cc, deriv = 1)^2)
    b <- mean(tukeyChi(rstar, cc, deriv = 2))
    tau2 <- scale^2 * a/b^2
    stddev <- sqrt(tau2)
    
    
    #R^2 under the restricted object
    df.int <- if (attr(object$terms, "intercept")) { 1L } else { 0L }
    resp.mean <- if (df.int == 1L) { sum(wgt * y) / sum(wgt) } else { 0 }
    yMy <- sum(wgt * (y - resp.mean)^2)
    rMr <- sum(wgt * residuals^2)
    # tukey's bi-square correction
    correc <- 1.207617 
    R2.reduced <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.restr,
                residuals   = residuals,
                wresid      = rfit$wresid,
                fitted      = fitted,
                weights     = weights,
                wgt         = wgt, 
                scale       = scale,
                stddev      = stddev,
                R2.org      = R2.org,
                R2.reduced  = R2.reduced,
                df.residual = so$df[2], 
                loglik      = ll.restr, 
                Sigma       = Sigma,                             
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                wt.bar      = NULL,
                iact        = rfit$iact,
                converged   = rfit$converged, 
                iter        = rfit$conv,
                bootout     = NULL, 
                control     = control)
  }
  
  OUT$model.org <- object
  OUT$se <- se 
  OUT$information <- 1 / stddev^2 * crossprod(X)
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
  
  class(OUT) <- c("restriktor", "conRLM", "conLM")
  
  OUT

}
