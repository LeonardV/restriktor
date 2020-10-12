### computes the F-bar, LRT-bar and score-bar test statistic ####
# REFs: 
# Silvapulle and Sen (2005). Constrained statistical inference. Chapter 2.
# Wolak, F. An exact test for multiple inequality and equality 
# constraints in the linear regression model Journal of the American 
# statistical association, 1987, 82, 782-793
conTestF.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                           p.distr = rnorm, parallel = "no", ncpus = 1L,
                           cl = NULL, seed = 1234, verbose = FALSE,
                           control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!inherits(object, "conLM")) {
    stop("Restriktor ERROR: object must be of class conLM.")
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
  X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # sample size
  n <- dim(X)[1]
  # weights
  w <- weights(model.org)
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
  vnames <- names(b.unrestr)
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
    # call quadprog
    b.eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               #b.unrestr = b.unrestr,
                               w         = w, 
                               Amat      = AmatG,
                               bvec      = bvecG, 
                               meq       = nrow(AmatG),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$qp$solution
    # fix estimates < tol to zero 
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    # compute global test statistic
    Ts <- c( t(b.restr - b.eqrestr) %*% solve(Sigma, b.restr - b.eqrestr) ) 
  } else if (type == "A") {
    b.eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               #b.unrestr = b.unrestr,
                               w         = w, 
                               Amat      = Amat,
                               bvec      = bvec, 
                               meq       = nrow(Amat),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$qp$solution
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    # compute test statistic for hypothesis test type A
    Ts <- c( t(b.restr - b.eqrestr) %*% solve(Sigma, b.restr - b.eqrestr) ) 
  } else if (type == "B") {
    if (meq.alt == 0L) {
      # compute test statistic for hypothesis test type B when no equalities are
      # preserved in the alternative hypothesis.
      Ts <- c( t(b.unrestr - b.restr) %*% solve(Sigma, b.unrestr - b.restr) ) 
    } else {
      if (meq.alt > 0L && meq.alt <= meq) {
        b.restr.alt <- con_solver_lm(X         = X, 
                                     y         = y, 
                                     #b.unrestr = b.unrestr,
                                     w         = w, 
                                     Amat      = Amat[1:meq.alt,,drop = FALSE],
                                     bvec      = bvec[1:meq.alt], 
                                     meq       = meq.alt,
                                     absval    = ifelse(is.null(control$absval), 1e-09, 
                                                        control$absval),
                                     maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                        control$maxit))$qp$solution
        b.restr.alt[abs(b.restr.alt) < tol] <- 0L
        names(b.restr.alt) <- vnames
        # compute test statistic for hypothesis test type B when some equalities may 
        # be preserved in the alternative hypothesis.
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
                                          control  = control,
                                          verbose  = verbose)
   } else if (boot == "model.based") {
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
                                           control  = control,
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
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              boot        = boot,
              model.org   = model.org)

  class(OUT) <- "conTest"

  OUT
}


# REF: Silvapulle and Sen (2005). Constrained statistical inference. Chapter 3.
# Wolak, F. An exact test for multiple inequality and equality constraints in the 
# linear regression model Journal of the American statistical association, 1987, 82, 782-793
conTestLRT.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                             p.distr = rnorm, parallel = "no", ncpus = 1L,
                             cl = NULL, seed = 1234, verbose = FALSE,
                             control = NULL, ...) {

  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!inherits(object, "conLM")) {
    stop("Restriktor ERROR: object must be of class conLM.")
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
  w <- weights(model.org)
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
  vnames <- names(b.unrestr)
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
    b.eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               #b.unrestr = b.unrestr, 
                               w         = w, 
                               Amat      = AmatG,
                               bvec      = bvecG, 
                               meq       = nrow(AmatG),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                                 control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                 control$maxit))$qp$solution
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    b.eqrestr <- as.vector(b.eqrestr)
    
    fitted0 <- X %*% b.eqrestr
    residuals0 <- y - fitted0
    object.eqrestr <- list(residuals = residuals0)
    object.eqrestr$weights   <- object$weights
    
    ll0 <- con_loglik_lm(object.eqrestr)
    ll1 <- object$loglik
    
    Ts <- -2*(ll0 - ll1)
  } else if (type == "A") {
    b.eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               #b.unrestr = b.unrestr,
                               w         = w, 
                               Amat      = Amat,
                               bvec      = bvec, 
                               meq       = nrow(Amat),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$qp$solution
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    b.eqrestr <- as.vector(b.eqrestr)
    
    fitted0 <- X %*% b.eqrestr
    residuals0 <- y - fitted0
    object.eqrestr <- list(residuals = residuals0)
    object.eqrestr$weights <- object$weights
    
    ll0 <- con_loglik_lm(object.eqrestr)
    ll1 <- object$loglik
    
    Ts <- -2*(ll0 - ll1)
  } else if (type == "B") {
      if (meq.alt == 0L) {
        
        ll0 <- object$loglik
        
        fitted1 <- X %*% object$b.unrestr
        residuals1 <- y - fitted1
        object.unrestr <- list(residuals = residuals1)
        object.unrestr$weights   <- object$weights
        ll1 <- con_loglik_lm(object.unrestr)
        
        Ts <- -2*(ll0 - ll1)
      }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        b.restr.alt <- con_solver_lm(X         = X, 
                                     y         = y, 
                                     #b.unrestr = b.unrestr,
                                     w         = w,
                                     Amat      = Amat[1:meq.alt,,drop=FALSE],
                                     bvec      = bvec[1:meq.alt], 
                                     meq       = meq.alt,
                                     absval    = ifelse(is.null(control$absval), 1e-09, 
                                                        control$absval),
                                     maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                        control$maxit))$qp$solution
        b.restr.alt[abs(b.restr.alt) < tol] <- 0L
        names(b.restr.alt) <- vnames

        ll0 <- object$loglik
        
        fitted1 <- X %*% b.restr.alt
        residuals1 <- y - fitted1
        object.restr.alt <- list(residuals = residuals1)
        object.restr.alt$weights <- object$weights
        ll1 <- con_loglik_lm(object.restr.alt)
        
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
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              boot        = boot,
              model.org   = model.org)

  class(OUT) <- "conTest"

  OUT
}


# REF: Robertson, Wright and Dykstra (1988, p. 321). Order constrained statistical inference.
# Global score 
conTestScore.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                               p.distr = rnorm, parallel = "no", ncpus = 1L,
                               cl = NULL, seed = 1234, verbose = FALSE,
                               control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!inherits(object, "conLM")) {
    stop("Restriktor ERROR: object must be of class conLM.")
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
  vnames <- names(b.unrestr)
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
    b.eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               #b.unrestr = b.unrestr,
                               w         = w, 
                               Amat      = AmatG, 
                               bvec      = bvecG, 
                               meq       = nrow(AmatG),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$qp$solution
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    
    res0 <- y - X %*% b.eqrestr
    res1 <- y - X %*% object$b.restr
    # score vector
    df0 <- n - (p - nrow(AmatG))  
    s20 <- sum(res0^2) / df0
    df1 <- n - (p - qr(Amat[0:meq,])$rank)
    s21 <- sum(res1^2) / df1
    G0 <- colSums(as.vector(res0) * w * X) / s20
    G1 <- colSums(as.vector(res1) * w * X) / s21
    # information matrix under the null-hypothesis
    I0 <- 1 / s20 * crossprod(X)#(t(X) %*% X)
    # score test-statistic
    Ts <- c((G0 - G1) %*% solve(I0, (G0 - G1)))
    ###############################################
    # df0 <- n - (p - nrow(AmatG))
    # df1 <- n - (p - qr(Amat[0:meq,])$rank)
    # s20 <- sum((y - X %*% b.eqrestr)^2) / df0
    # s21 <- sum((y - X %*% b.restr)^2) / df1
    # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.eqrestr)))
    # d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b.restr)))
    # I0 <- 1/s20 * (t(X) %*% X)
    # Ts2 <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
    ###############################################
  } else if (type == "A") {
    b.eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               #b.unrestr = b.unrestr,
                               w         = w, 
                               Amat      = Amat,
                               bvec      = bvec, 
                               meq       = nrow(Amat),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$qp$solution
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    
    res0 <- y - X %*% b.eqrestr
    res1 <- y - X %*% object$b.restr
    # score vector
    df0 <- n - (p - nrow(Amat))
    s20 <- sum(res0^2) / df0
    df1 <- n - (p - qr(Amat[0:meq,])$rank)
    s21 <- sum(res1^2) / df1
    G0 <- colSums(as.vector(res0) * w * X) / s20
    G1 <- colSums(as.vector(res1) * w * X) / s21
    # information matrix under the null-hypothesis
    I0 <- 1 / s20 * crossprod(X)#(t(X) %*% X)
    # score test-statistic
    Ts <- c((G0 - G1) %*% solve(I0, (G0 - G1))) 
    #############################################
    # df0 <- n - (p - nrow(Amat))
    # df1 <- n - (p - qr(Amat[0:meq,])$rank)
    # s20 <- sum((y - X %*% b.eqrestr)^2) / df0
    # s21 <- sum((y - X %*% b.restr)^2) / df1
    # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.eqrestr)))
    # d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b.restr)))
    # I0 <- 1/s20 * (t(X) %*% X)
    # Ts2 <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
    ###############################################
  } else if (type == "B") {
    if (meq.alt == 0L) {
      res0 <- y - X %*% object$b.restr
      res1 <- y - X %*% object$b.unrestr
      
      #res0 <- residuals(object)
      #res1 <- residuals(model.org)
      # score vector
      df0 <- n - (p - qr(Amat[0:meq,])$rank)
      s20 <- sum(res0^2) / df0
      df1 <- n - p
      s21 <- sum(res1^2) / df1
      G0 <- colSums(as.vector(res0) * w * X) / s20
      G1 <- colSums(as.vector(res1) * w * X) / s21
      # information matrix under the null-hypothesis
      I0 <- 1 / s20 * crossprod(X)#(t(X) %*% X)
      # score test-statistic
      Ts <- c((G0 - G1) %*% solve(I0, (G0 - G1)))
      #########
      # df0 <- n - (p - qr(Amat[0:meq,])$rank)
      # df1 <- n - p
      # s20 <- sum((y - X %*% b.restr)^2) / df0
      # s21 <- sum((y - X %*% b.unrestr)^2) / df1
      # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.restr)))
      # d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b.unrestr)))
      # I0 <- 1/s20 * (t(X) %*% X)
      # Ts2 <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
      ###############################################
      } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        b.restr.alt <- con_solver_lm(X         = X, 
                                     y         = y, 
                                     #b.unrestr = b.unrestr,
                                     w         = w,
                                     Amat      = Amat[1:meq.alt,,drop=FALSE],
                                     bvec      = bvec[1:meq.alt], 
                                     meq       = meq.alt,
                                     absval    = ifelse(is.null(control$absval), 1e-09,
                                                        control$absval),
                                     maxit     = ifelse(is.null(control$maxit), 1e04,
                                                        control$maxit))$qp$solution
        b.restr.alt[abs(b.restr.alt) < tol] <- 0L
        names(b.restr.alt) <- vnames
        
        res0 <- y - X %*% object$b.restr
        res1 <- y - X %*% b.restr.alt
        
        # res0 <- residuals(object)
        # res1 <- y - X %*% b.restr.alt
        # score vector
        df0 <- n - (p - qr(Amat[0:meq,])$rank)
        s20 <- sum(res0^2) / df0
        df1 <- n - (p - qr(Amat[0:meq,])$rank)
        s21 <- sum(res1^2) / df1
        G0 <- colSums(as.vector(res0) * w * X) / s20
        G1 <- colSums(as.vector(res1) * w * X) / s21
        # information matrix
        I0 <- 1 / s20 * crossprod(X)#(t(X) %*% X)
        # score test-statistic
        Ts <- c((G0 - G1) %*% solve(I0, (G0 - G1)))
        ########################
        # df1 <- df0 <- n - (p - qr(Amat[0:meq,])$rank)
        # s20 <- sum((y - X %*% b.restr)^2) / df0
        # s21 <- sum((y - X %*% b.restr.alt)^2) / df1
        # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.restr)))
        # d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b.restr.alt)))
        # I0 <- 1/s20 * (t(X) %*% X)
        # Ts2 <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
        ###############################################
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
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              boot        = boot,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}



## hypothesis test type C is based on a t-distribution.
## intersection-union test 
## REF: S. Sasabuchi (1980). A Test of a Multivariate Normal Mean with 
## Composite Hypotheses Determined by Linear Inequalities. Biometrika Trust, 67 (2), 429-439.
conTestC.restriktor <- function(object, ...) {
  
  if (!(inherits(object, "restriktor"))) {
    stop("Restriktor ERROR: object must be of class restriktor.")
  }
  
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  
  if (meq > 0L) {
    stop("Restriktor ERROR: test not applicable with equality restrictions.")
  }
  
  b.unrestr   <- object$b.unrestr
  Sigma       <- object$Sigma
  df.residual <- object$df.residual
    
  Ts <- as.vector(min((Amat %*% b.unrestr - bvec) / 
                        sqrt(diag(Amat %*% Sigma %*% t(Amat))))) 
  
  pvalue <- 1 - pt(Ts, df.residual)

  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq.alt     = 0L,
              iact        = object$iact,
              type        = "C",
              test        = "t",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b.unrestr   = b.unrestr,
              b.restr     = object$b.restr,
              Sigma       = Sigma,
              R2.org      = ifelse(!is.null(object$R2.org), object$R2.org, as.numeric(NA)),
              R2.reduced  = ifelse(!is.null(object$R2.reduced), object$R2.reduced, as.numeric(NA)),
              boot        = "no",
              model.org   = object$model.org)
  
#  OUT <- list(OUT)
#    names(OUT) <- "C"
  
  class(OUT) <- "conTest"
  
  OUT
}
