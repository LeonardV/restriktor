### computes the F-bar, LRT-bar and score-bar test statistic ####
## REF: Silvapulle and Sen (2005). Constrained statistical inference. Chapter 4.
conTestF.conGLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                            p.distr = rnorm, parallel = "no", ncpus = 1L,
                            cl = NULL, seed = 1234, verbose = FALSE,
                            control = NULL, ...) {

  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!inherits(object, "conGLM")) {
    stop("Restriktor ERROR: object must be of class conGLM.")
  }
  if (type != "global") {
    type <- toupper(type)
  }    
  if(!(type %in% c("A","B","global"))) {
    stop("Restriktor ERROR: type must be \"A\", \"B\", or \"global\"")
  }
  if(!(boot %in% c("no","residual","model.based","parametric","mix.weights"))) {
    stop("Restriktor ERROR: boot method unknown.")
  }
  if (boot == "residual") {
    boot <- "model.based"
  }
  
  # original model
  model.org <- object$model.org
  # model matrix
  #X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  #y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # sample size
#  n <- dim(X)[1]
  # weights
  #w <- weights(model.org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- vcov(model.org) 
  # parameter estimates
  b.unrestr <- object$b.unrestr
  b.restr <- object$b.restr
  b.eqrestr <- NULL
  b.restr.alt <- NULL
  # length parameter vector
  p <- length(b.unrestr)
  # variable names
  #vnames <- names(b.unrestr)
  # constraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  control <- c(object$control, control)
  # remove duplicated elements from control list
  control <- control[!duplicated(control)]
  # get tolerance for control if exists
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  }
  
  # check for intercept                                          
  intercept <- any(attr(terms(model.org), "intercept"))
  if (type == "global") {
    if (intercept) { 
      AmatG <- cbind(rep(0, (p - 1)), diag(rep(1, p - 1))) 
    } else {
      AmatG <- diag(1, p)
      for (i in 1:p) {
        AmatG[i,i-1] <- -1
      }
      AmatG <- AmatG[-1,]
    }
    AmatX <- AmatG %*% (diag(rep(1, p)) - t(Amat) %*%            
                          MASS::ginv(Amat %*% t(Amat)) %*% Amat)
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "type.org") <- "global"
    } else {
      # remove all rows with only zeros
      AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
      rAmatX <- GaussianElimination(t(AmatX), tol = tol)
      AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
    }
    AmatG <- rbind(AmatX, Amat)
    bvecG <- c(rep(0, nrow(AmatX)), bvec)
    attr(Amat, "Amat.global") <- AmatG
    attr(bvec, "bvec.global") <- bvecG
  }
  
  if (type == "global") {
    CALL <- list(object = model.org, constraints = AmatG, rhs = bvecG, 
                 neq = nrow(AmatG), se = "none", mix.weights = "none")
    glm.fit <- do.call("restriktor", CALL)
    
    b.eqrestr <- coef(glm.fit)
    # compute global test statistic
    Ts <- c( t(b.restr - b.eqrestr) %*% solve(Sigma, b.restr - b.eqrestr) ) 
  } else if (type == "A") {
    CALL <- list(object = model.org, constraints = Amat, rhs = bvec, 
                 neq = nrow(Amat), se = "none", mix.weights = "none")
    glm.fit <- do.call("restriktor", CALL)
    
    b.eqrestr <- coef(glm.fit)
    # compute test statistic for hypothesis test type A
    Ts <- c( t(b.restr - b.eqrestr) %*% solve(Sigma, b.restr - b.eqrestr) ) 
  } else if (type == "B") {
    if (meq.alt == 0L) {
      # compute test statistic for hypothesis test type B when no equalities are
      # preserved in the alternative hypothesis.
      Ts <- c( t(b.unrestr - b.restr) %*% solve(Sigma, b.unrestr - b.restr) ) 
    } else {
      if (meq.alt > 0L && meq.alt <= meq) {
        # compute test statistic for hypothesis test type B when some equalities may 
        # be preserved in the alternative hypothesis.
        CALL <- list(object = model.org, constraints = Amat[1:meq.alt,,drop = FALSE], 
                     rhs = bvec[1:meq.alt], neq = meq.alt, 
                     se = "none", mix.weights = "none")
        glm.fit <- do.call("restriktor", CALL)
        
        b.restr.alt <- coef(glm.fit)
        
        Ts <- c( t(b.restr - b.restr.alt) %*% solve(Sigma, b.restr - b.restr.alt) ) 
      } else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  # The test statistics based on inequality constraints are often
  # distributed as mixtures of chi-squares. These mixing weights can be computed
  # using the multivariate normal distribution with additional Monte Carlo steps
  # or via bootstrapping. The pvalue can also be computed directly via 
  # the parametric bootstrap or model based bootstrap, without fist computing 
  # the mixing weights.
  if (!(attr(object$wt.bar, "method") == "none") && boot == "no") {
    wt.bar <- object$wt.bar
    pvalue <- con_pvalue_Fbar(wt.bar      = wt.bar, 
                              Ts.org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq.alt     = meq.alt)
    #wt.mix <- wt.bar
    #attributes(wt.mix) <- NULL
    attr(pvalue, "wt.bar") <- wt.bar
    attr(pvalue, "wt.bar.method") <- attr(wt.bar, "method")
  } else if (boot == "parametric") {
    if (!is.function(p.distr)) {
      p.distr <- get(p.distr, mode = "function")
    }
    arguments <- list(...)
    pnames <- names(formals(p.distr))
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    pm <- names(arguments)[pm > 0L]
    formals(p.distr)[pm] <- unlist(arguments[pm])
    
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts.org   = Ts, 
                                         type     = type, 
                                         test     = "F", 
                                         meq.alt  = meq.alt, 
                                         R        = R, 
                                         p.distr  = p.distr,
                                         parallel = parallel,
                                         ncpus    = ncpus, 
                                         cl       = cl,
                                         seed     = seed, 
                                         verbose  = verbose)
  } else if (boot == "model.based") {
    
    if (!(family(model.org)$family %in% c("gaussian","gamma"))) {
      stop("Restriktor ERROR: the model.based bootstrap is not available for the ", family(model.org)$family, " family.")
    }
    
    pvalue <- con_pvalue_boot_model_based(object, 
                                          Ts.org   = Ts, 
                                          type     = type, 
                                          test     = "F", 
                                          meq.alt  = meq.alt,
                                          R        = R,
                                          parallel = parallel, 
                                          ncpus    = ncpus,
                                          cl       = cl, 
                                          seed     = seed, 
                                          verbose  = verbose)
  } else {
    pvalue <- as.numeric(NA)
  } 
  
  # necessary for the print function
  if (!is.null(attr(type, "type.org"))) {
    type <- "global"
  }
  
  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq.alt     = meq.alt,
              iact        = object$iact,
              type        = type,
              test        = "F",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b.eqrestr   = b.eqrestr,
              b.unrestr   = b.unrestr,
              b.restr     = b.restr,
              b.restr.alt = b.restr.alt,
              boot        = boot,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}


# REF: Silvapulle and Sen (2005). Constrained statistical inference. Chapter 4.
conTestLRT.conGLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                              p.distr = rnorm, parallel = "no", ncpus = 1L,
                              cl = NULL, seed = 1234, verbose = FALSE,
                              control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!inherits(object, "conGLM")) {
    stop("Restriktor ERROR: object must be of class conGLM.")
  }
  if (type != "global") {
    type <- toupper(type)
  }  
  if(!(type %in% c("A","B","global"))) {
    stop("Restriktor ERROR: type must be \"A\", \"B\", or \"global\"")
  }
  if(!(boot %in% c("no", "residual", "model.based", "parametric", "mix.weights"))) {
    stop("Restriktor ERROR: boot method unknown.")
  }
  if (boot == "residual") {
    boot <- "model.based"
  }
  
  # original model
  model.org <- object$model.org
  # model matrix
  X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # weights
  #w <- weights(model.org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- vcov(model.org) 
  # parameter estimates
  b.unrestr <- object$b.unrestr
  b.restr <- object$b.restr
  b.eqrestr <- NULL
  b.restr.alt <- NULL
  # length parameter vector
  p <- length(b.unrestr)
  # variable names
  #vnames <- names(b.unrestr)
  # constraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  control <- c(object$control, control)
  # remove duplicated elements from control list
  control <- control[!duplicated(control)]
  # get tolerance for control if exists
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  }
  
  # check for intercept                                          
  intercept <- any(attr(terms(model.org), "intercept"))
  if (type == "global") {
    if (intercept) { 
      AmatG <- cbind(rep(0, (p - 1)), diag(rep(1, p - 1))) 
    } else {
      AmatG <- diag(1, p)
      for (i in 1:p) {
        AmatG[i,i-1] <- -1
      }
      AmatG <- AmatG[-1,]
    }
    AmatX <- AmatG %*% (diag(rep(1, p)) - t(Amat) %*%            
                          MASS::ginv(Amat %*% t(Amat)) %*% Amat)
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "type.org") <- "global"
    } else {
      # remove all rows with only zeros
      AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
      rAmatX <- GaussianElimination(t(AmatX), tol = tol)
      AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
    }
    AmatG <- rbind(AmatX, Amat)
    bvecG <- c(rep(0, nrow(AmatX)), bvec)
    attr(Amat, "Amat.global") <- AmatG
    attr(bvec, "bvec.global") <- bvecG
  }
  
  if (type == "global") {  
    CALL <- list(object = model.org, constraints = AmatG, 
                 rhs = bvecG, neq = nrow(AmatG), se = "none", mix.weights = "none")
    glm.fit <- do.call("restriktor", CALL)
    
    ll0 <- glm.fit$loglik
    ll1 <- object$loglik
    
    Ts <- -2*(ll0 - ll1)
  } else if (type == "A") {
    CALL <- list(object = model.org, constraints = Amat, 
                 rhs = bvec, neq = nrow(Amat), se = "none", mix.weights = "none")
    glm.fit <- do.call("restriktor", CALL)
    
    ll0 <- glm.fit$loglik
    ll1 <- object$loglik
    
    Ts <- -2*(ll0 - ll1)
  } else if (type == "B") {
    if (meq.alt == 0L) {
      ll0 <- object$loglik 
      ll1 <- logLik(model.org)
      
      Ts <- -2*(ll0 - ll1)
    }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        CALL <- list(object = model.org, constraints = Amat[1:meq.alt,,drop=FALSE], 
                     rhs = bvec[1:meq.alt], neq = meq.alt, se = "none", mix.weights = "none")
        glm.fit <- do.call("restriktor", CALL)
        
        ll0 <- glm.fit$loglik
        ll1 <- object$loglik
        
        Ts <- -2*(ll0 - ll1)
      }
      else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  if (!(attr(object$wt.bar, "method") == "none") && boot == "no") { 
    wt.bar <- object$wt.bar
    pvalue <- con_pvalue_Fbar(wt.bar      = wt.bar, 
                              Ts.org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq.alt     = meq.alt)
    attr(pvalue, "wt.bar") <- wt.bar
  } else if (boot == "parametric") {
    if (!is.function(p.distr)) {
      p.distr <- get(p.distr, mode = "function")
    }
    arguments <- list(...)
    pnames <- names(formals(p.distr))
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    pm <- names(arguments)[pm > 0L]
    formals(p.distr)[pm] <- unlist(arguments[pm])
    
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts.org   = Ts, 
                                         type     = type, 
                                         test     = "LRT", 
                                         meq.alt  = meq.alt,
                                         R        = R, 
                                         p.distr  = p.distr,
                                         parallel = parallel,
                                         ncpus    = ncpus, 
                                         cl       = cl,
                                         seed     = seed, 
                                         verbose  = verbose)
  } else if (boot == "model.based") {
    
    if (!(family(model.org)$family %in% c("gaussian","gamma"))) {
      stop("Restriktor ERROR: the model.based bootstrap is not available for the ", family(model.org)$family, " family.")
    }
    
    pvalue <- con_pvalue_boot_model_based(object, 
                                          Ts.org   = Ts, 
                                          type     = type, 
                                          test     = "LRT",
                                          meq.alt  = meq.alt,
                                          R        = R, 
                                          parallel = parallel,
                                          ncpus    = ncpus, 
                                          cl       = cl,
                                          seed     = seed, 
                                          verbose  = verbose)
  } else {
    pvalue <- as.numeric(NA)
  }  
  
  # necessary for the print function
  if (!is.null(attr(type, "type.org"))) {
    type <- "global"
  }
  
  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq.alt     = meq.alt,
              iact        = object$iact,
              type        = type,
              test        = "LRT",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b.eqrestr   = b.eqrestr,
              b.unrestr   = b.unrestr,
              b.restr     = b.restr,
              b.restr.alt = b.restr.alt,
              boot        = boot,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}



# REF: Robertson, Silvapulle and Sen (2005, chapter 4)
conTestScore.conGLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                               p.distr = rnorm, parallel = "no", ncpus = 1L,
                               cl = NULL, seed = 1234, verbose = FALSE,
                               control = NULL, ...) {
  
  #rename for internal use
  meq.alt <- neq.alt

  # checks
  if (!("conGLM" %in% class(object))) {
    stop("Restriktor ERROR: object must be of class conGLM.")
  }
  if (type != "global") {
    type <- toupper(type)
  }
  if(!(type %in% c("A","B","global"))) {
    stop("Restriktor ERROR: type must be \"A\", \"B\", or \"global\"")
  }
  if(!(boot %in% c("no", "residual", "model.based", "parametric", "mix.weights"))) {
    stop("Restriktor ERROR: boot method unknown.")
  }
  if (boot == "residual") {
    boot <- "model.based"
  }

  # original model
  model.org <- object$model.org
  # family
  fam <- family(model.org)
  # model matrix
  X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  #y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- vcov(model.org)
  # sample size
  n <- dim(X)[1]
  # number of parameters
  p <- dim(X)[2]
  # weights
  w <- weights(model.org)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  #W <- diag(w)
  # parameter estimates
  b.unrestr <- object$b.unrestr
  b.restr <- object$b.restr
  b.eqrestr <- NULL
  b.restr.alt <- NULL
  # length parameter vector
  p <- length(b.unrestr)
  # variable names
  #vnames <- names(b.unrestr)
  # restraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  control <- c(object$control, control)
  # remove duplicated elements from control list
  control <- control[!duplicated(control)]
  # get tolerance for control if exists
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  }

  # check for intercept
  intercept <- any(attr(terms(model.org), "intercept"))
  if (type == "global") {
    if (intercept) {
      AmatG <- cbind(rep(0, (p - 1)), diag(rep(1, p - 1)))
    } else {
      AmatG <- diag(1, p)
      for (i in 1:p) {
        AmatG[i,i-1] <- -1
      }
      AmatG <- AmatG[-1,]
    }
    AmatX <- AmatG %*% (diag(rep(1, p)) - t(Amat) %*%
                          MASS::ginv(Amat %*% t(Amat)) %*% Amat)

    if (all(abs(AmatX) < tol)) {
      type <- "A"
      attr(type, "type.org") <- "global"
    } else {
      # remove all rows with only zeros
      AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
      rAmatX <- GaussianElimination(t(AmatX), tol = tol)
      AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
    }
    AmatG <- rbind(AmatX, Amat)
    bvecG <- c(rep(0, nrow(AmatX)), bvec)
    attr(Amat, "Amat.global") <- AmatG
    attr(bvec, "bvec.global") <- bvecG
  }

  if (type == "global") {
    CALL <- list(object = model.org, constraints = AmatG, 
                 rhs = bvecG, neq = nrow(AmatG), se = "none", mix.weights = "none")
    glm.fit <- do.call("restriktor", CALL)
    
    b.eqrestr <- coef(glm.fit)
    
    res0 <- residuals(glm.fit, "working")
    res1 <- residuals(object, "working")
    wres0 <- as.vector(res0) * weights(glm.fit, "working")
    wres1 <- as.vector(res1) * weights(object, "working")
    
    df0 <- n - (p - nrow(AmatG)) 
    df1 <- n - (p - qr(Amat[0:meq,])$rank)
    
    dispersion0 <- if (fam$family %in% c("poisson", "binomial")) 1
      else sum(wres0^2) / df0
    dispersion1 <- if (fam$family %in% c("poisson", "binomial")) 1
      else sum(wres1^2) / df1
    
    # score vector
    G0 <- colSums(wres0 * X / dispersion0)
    G1 <- colSums(wres1 * X / dispersion1)
    
    # information matrix under the null-hypothesis
    W <- diag(weights(glm.fit))
    I0 <- t(X) %*% W %*% X / dispersion0
    
    # score test-statistic
    Ts <- c((G0 - G1) %*% solve(I0 , (G0 - G1)))
    ###############################################
  } else if (type == "A") {
    CALL <- list(object = model.org, constraints = Amat, 
                 rhs = bvec, neq = nrow(Amat), se = "none", mix.weights = "none")
    glm.fit <- do.call("restriktor", CALL)
    
    b.eqrestr <- coef(glm.fit)
    
    res0 <- residuals(glm.fit, "working")
    res1 <- residuals(object, "working")
    wres0 <- as.vector(res0) * weights(glm.fit, "working")
    wres1 <- as.vector(res1) * weights(object, "working")

    df0 <- n - (p - nrow(Amat)) 
    df1 <- n - (p - qr(Amat[0:meq,])$rank)
    
    dispersion0 <- if (fam$family %in% c("poisson", "binomial")) 1
      else sum(wres0^2) / df0
    dispersion1 <- if (fam$family %in% c("poisson", "binomial")) 1
      else sum(wres1^2) / df1
    
    # score vector
    G0 <- colSums(wres0 * X / dispersion0)
    G1 <- colSums(wres1 * X / dispersion1)
    
    # information matrix under the null-hypothesis
    W <- diag(weights(glm.fit))
    I0 <- t(X) %*% W %*% X / dispersion0
    
    # score test-statistic
    Ts <- c((G0 - G1) %*% solve(I0 , (G0 - G1)))
    ###############################################
  } else if (type == "B") {
    if (meq.alt == 0L) {
      res0 <- residuals(object, "working")
      res1 <- residuals(model.org, "working")
      wres0 <- as.vector(res0) * weights(object, "working")
      wres1 <- as.vector(res1) * weights(model.org, "working")
      
      df0 <- n - (p - qr(Amat[0:meq,])$rank)
      df1 <- n - p
      
      dispersion0 <- if (fam$family %in% c("poisson", "binomial")) 1
        else sum(wres0^2) / df0
      dispersion1 <- if (fam$family %in% c("poisson", "binomial")) 1
        else sum(wres1^2) / df1
      
      # score vector
      G0 <- colSums(wres0 * X / dispersion0)
      G1 <- colSums(wres1 * X / dispersion1)
      
      # information matrix under the null-hypothesis
      W <- diag(weights(object))
      I0 <- t(X) %*% W %*% X / dispersion0
      
      # score test-statistic
      Ts <- c((G0 - G1) %*% solve(I0 , (G0 - G1)))
      ###############################################
    }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        CALL <- list(object = model.org, constraints = Amat[1:meq.alt,,drop=FALSE], 
                     rhs = bvec[1:meq.alt], neq = meq.alt, se = "none", mix.weights = "none")
        glm.fit <- do.call("restriktor", CALL)
        
        b.restr.alt <- coef(glm.fit)
        
        res0 <- residuals(object, "working")
        res1 <- residuals(glm.fit, "working")
        wres0 <- as.vector(res0) * weights(glm.fit, "working")
        wres1 <- as.vector(res1) * weights(object, "working")
        
        df0 <- n - (p - qr(Amat[0:meq,])$rank)
        df1 <- n - (p - qr(Amat[0:meq,])$rank)
        
        dispersion0 <- if (fam$family %in% c("poisson", "binomial")) 1
          else sum(wres0^2) / df0
        dispersion1 <- if (fam$family %in% c("poisson", "binomial")) 1
          else sum(wres1^2) / df1
        
        # score vector
        G0 <- colSums(wres0 * X / dispersion0)
        G1 <- colSums(wres1 * X / dispersion1)
        
        # information matrix under the null-hypothesis
        W <- diag(weights(object))
        I0 <- t(X) %*% W %*% X / dispersion0
        
        # score test-statistic
        Ts <- c((G0 - G1) %*% solve(I0 , (G0 - G1)))
      }
      else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  }

  if (!(attr(object$wt.bar, "method") == "none") && boot == "no") {
    wt.bar <- object$wt.bar
    pvalue <- con_pvalue_Fbar(wt.bar      = wt.bar,
                              Ts.org      = Ts,
                              df.residual = df.residual,
                              type        = type,
                              Amat        = Amat,
                              bvec        = bvec,
                              meq         = meq,
                              meq.alt     = meq.alt)
    attr(pvalue, "wt.bar") <- wt.bar
  } else if (boot == "parametric") {
    if (!is.function(p.distr)) {
      p.distr <- get(p.distr, mode = "function")
    }
    arguments <- list(...)
    pnames <- names(formals(p.distr))
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    pm <- names(arguments)[pm > 0L]
    formals(p.distr)[pm] <- unlist(arguments[pm])
    
    pvalue <- con_pvalue_boot_parametric(object,
                                         Ts.org   = Ts,
                                         type     = type,
                                         test     = "score",
                                         meq.alt  = meq.alt,
                                         R        = R,
                                         p.distr  = p.distr,
                                         parallel = parallel,
                                         ncpus    = ncpus,
                                         cl       = cl,
                                         seed     = seed,
                                         verbose  = verbose)
  } else if (boot == "model.based") {
    stop("Restriktor ERROR: model.based bootstrap is not implemented yet.")
    pvalue <- con_pvalue_boot_model_based(object,
                                          Ts.org   = Ts,
                                          type     = type,
                                          test     = "score",
                                          meq.alt  = meq.alt,
                                          R        = R,
                                          parallel = parallel,
                                          ncpus    = ncpus,
                                          cl       = cl,
                                          seed     = seed,
                                          verbose  = verbose)
  } else {
    pvalue <- as.numeric(NA)
  }

  # necessary for the print function
  if (!is.null(attr(type, "type.org"))) {
    type <- "global"
  }

  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq.alt     = meq.alt,
              iact        = object$iact,
              type        = type,
              test        = "Score",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b.eqrestr   = b.eqrestr,
              b.unrestr   = b.unrestr,
              b.restr     = b.restr,
              b.restr.alt = b.restr.alt,
              boot        = boot,
              model.org   = model.org)

  class(OUT) <- "conTest"

  OUT
}
