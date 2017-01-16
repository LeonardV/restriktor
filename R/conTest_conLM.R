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
  meq_alt <- neq.alt
  
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
  model_org <- object$model_org
  # model matrix
  X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # sample size
  n <- dim(X)[1]
  # weights
  w <- weights(model_org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- vcov(model_org) 
  # parameter estimates
  b_unrestr <- object$b_unrestr
  b_restr <- object$b_restr
  b_eqrestr <- NULL
  b_restr_alt <- NULL
  # length parameter vector
  p <- length(b_unrestr)
  # variable names
  vnames <- names(b_unrestr)
  # constraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  #control
  control <- object$control
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  }
  
  if (is.null(control$tol)) {
    tol <- sqrt(.Machine$double.eps)
  } else {
    tol <- control$tol
  }
  
  # check for intercept                                          
  intercept <- any(attr(terms(model_org), "intercept"))
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
        attr(type, "type_org") <- "global"
    } else {
      # remove all rows with only zeros
      AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
      rAmatX <- GaussianElimination(t(AmatX), tol = tol)
      AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
    }
    AmatG <- rbind(AmatX, Amat)
    bvecG <- c(rep(0, nrow(AmatX)), bvec)
      attr(Amat, "Amat_global") <- AmatG
      attr(bvec, "bvec_global") <- bvecG
  }
  
  if (type == "global") {
    # call quadprog
    b_eqrestr <- con_solver_lm(X         = X, 
                            y         = y, 
                            b_unrestr = b_unrestr,
                            w         = w, 
                            Amat      = AmatG,
                            bvec      = bvecG, 
                            meq       = nrow(AmatG),
                            absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                            maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
  
    # fix estimates < tol to zero 
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    # compute global test statistic
    Ts <- c(t(b_restr - b_eqrestr) %*% solve(Sigma, b_restr - b_eqrestr))
  } else if (type == "A") {
    b_eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               b_unrestr = b_unrestr,
                               w         = w, 
                               Amat      = Amat,
                               bvec      = bvec, 
                               meq       = nrow(Amat),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    # compute test statistic for hypothesis test type A
    Ts <- c(t(b_restr - b_eqrestr) %*% solve(Sigma, b_restr - b_eqrestr))
  } else if (type == "B") {
    if (meq_alt == 0L) {
      # compute test statistic for hypothesis test type B when no equalities are
      # preserved in the alternative hypothesis.
      Ts <- c(t(b_unrestr - b_restr) %*% solve(Sigma, b_unrestr - b_restr))
    } else {
      if (meq_alt > 0L && meq_alt <= meq) {
        b_restr_alt <- con_solver_lm(X         = X, 
                                  y         = y, 
                                  b_unrestr = b_unrestr,
                                  w         = w, 
                                  Amat      = Amat[1:meq_alt,,drop = FALSE],
                                  bvec      = bvec[1:meq_alt], 
                                  meq       = meq_alt,
                                  absval    = ifelse(is.null(control$absval), 1e-09, 
                                                     control$absval),
                                  maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                     control$maxit))$solution
        b_restr_alt[abs(b_restr_alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps),                                        
                                              control$tol)] <- 0L
        names(b_restr_alt) <- vnames
        # compute test statistic for hypothesis test type B when some equalities may 
        # be preserved in the alternative hypothesis.
        Ts <- c(t(b_restr - b_restr_alt) %*% solve(Sigma, b_restr - b_restr_alt))
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
  
  if (!is.null(object$wt) && boot == "no") {
    wt <- object$wt
    # is this fool proof? 
    # The number of bootstrap samples must be large enough to avoid spurious results.
    wt <- rev(wt)
    if (attr(object$wt, "bootWt")) {
      if (attr(object$wt, "bootWt.R") < 999) {
        stop("Restriktor ERROR: increase the number of bootstrap draws. Preferably to a large number e.g., bootWt.R = 99999")
      }
      wt.idx <- which(wt == 0)
      wt <- wt[-wt.idx]
    }
    
    pvalue <- con_pvalue_Fbar(wt          = wt, 
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
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
                                          Ts_org   = Ts, 
                                          type     = type, 
                                          test     = "F", 
                                          meq_alt  = meq_alt, 
                                          R        = R, 
                                          p.distr  = p.distr,
                                          parallel = parallel,
                                          ncpus    = ncpus, 
                                          cl       = cl,
                                          seed     = seed, 
                                          verbose  = verbose)
   } else if (boot == "model.based") {
     pvalue <- con_pvalue_boot_model_based(object, 
                                           Ts_org   = Ts, 
                                           type     = type, 
                                           test     = "F", 
                                           meq_alt  = meq_alt,
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
  if (!is.null(attr(type, "type_org"))) {
    type <- "global"
  }
  
  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq_alt     = meq_alt,
              iact        = object$iact,
              type        = type,
              test        = "F",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b_eqrestr   = b_eqrestr,
              b_unrestr   = b_unrestr,
              b_restr     = b_restr,
              b_restr_alt = b_restr_alt,
              Sigma       = Sigma,
              R2_org      = object$R2_org,
              R2_reduced  = object$R2_reduced,
              boot        = boot,
              model_org   = model_org)

  OUT <- list(OUT)
    names(OUT) <- type
  
  class(OUT) <- "conTest"

  OUT
}


# REF: Silvapulle and Sen (2005). Constrained statistical inference. Chapter 3.
# Wolak, F. An exact test for multiple inequality and equality constraints in the linear regression model Journal of the American statistical association, 1987, 82, 782-793
conTestLRT.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                             p.distr = rnorm, parallel = "no", ncpus = 1L,
                             cl = NULL, seed = 1234, verbose = FALSE,
                             control = NULL, ...) {

  # rename for internal use
  meq_alt <- neq.alt
  
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
  model_org <- object$model_org
  # model matrix
  X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # weights
  w <- weights(model_org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- vcov(model_org) 
  # parameter estimates
  b_unrestr <- object$b_unrestr
  b_restr <- object$b_restr
  b_eqrestr <- NULL
  b_restr_alt <- NULL
  # length parameter vector
  p <- length(b_unrestr)
  # variable names
  vnames <- names(b_unrestr)
  # constraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  #control
  control <- object$control
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  }
  
  if (is.null(control$tol)) {
    tol <- sqrt(.Machine$double.eps)
  } else {
    tol <- control$tol
  }
  
  # check for intercept                                          
  intercept <- any(attr(terms(model_org), "intercept"))
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
        attr(type, "type_org") <- "global"
    } else {
      # remove all rows with only zeros
      AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
      rAmatX <- GaussianElimination(t(AmatX), tol = tol)
      AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
    }
    AmatG <- rbind(AmatX, Amat)
    bvecG <- c(rep(0, nrow(AmatX)), bvec)
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
  }
  
  if (type == "global") {  
    b_eqrestr <- con_solver_lm(X         = X, 
                            y         = y, 
                            b_unrestr = b_unrestr, 
                            w         = w, 
                            Amat      = AmatG,
                            bvec      = bvecG, 
                            meq       = nrow(AmatG),
                            absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                            maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    ll0.out <- con_loglik_lm(X = X, 
                             y = y, 
                             b = b_eqrestr, 
                             w = w)
    ll0 <- ll0.out$loglik
    
    ll1 <- object$loglik
    Ts <- -2*(ll0 - ll1)
  } else if (type == "A") {
    b_eqrestr <- con_solver_lm(X         = X, 
                            y         = y, 
                            b_unrestr = b_unrestr,
                            w         = w, 
                            Amat      = Amat,
                            bvec      = bvec, 
                            meq       = nrow(Amat),
                            absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                            maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    ll0.out <- con_loglik_lm(X = X, 
                             y = y, 
                             b = b_eqrestr, 
                             w = w)
    ll0 <- ll0.out$loglik
    ll1 <- object$loglik
    
    Ts <- -2*(ll0 - ll1)
  } else if (type == "B") {
      if (meq_alt == 0L) {
        ll0 <- object$loglik
        ll1.out <- con_loglik_lm(X = X, 
                                 y = y, 
                                 b = b_unrestr, 
                                 w = w)
        ll1 <- ll1.out$loglik
        Ts <- -2*(ll0 - ll1)
      }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq_alt > 0L && meq_alt <= meq) {
        b_restr_alt <- con_solver_lm(X         = X, 
                                  y         = y, 
                                  b_unrestr = b_unrestr,
                                  w         = w,
                                  Amat      = Amat[1:meq_alt,,drop=FALSE],
                                  bvec      = bvec[1:meq_alt], 
                                  meq       = meq_alt,
                                  absval    = ifelse(is.null(control$absval), 1e-09, 
                                                     control$absval),
                                  maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                     control$maxit))$solution
        b_restr_alt[abs(b_restr_alt) < ifelse(is.null(control$tol),                                        
                                              sqrt(.Machine$double.eps),                                        
                                              control$tol)] <- 0L
        names(b_restr_alt) <- vnames

        ll0 <- object$loglik
        ll1.out <- con_loglik_lm(X = X, 
                                 y = y, 
                                 b = b_restr_alt, 
                                 w = w)
        ll1 <- ll1.out$loglik
        Ts <- -2*(ll0 - ll1)
      }
      else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  if (!is.null(object$wt) && boot == "no") { 
    wt <- object$wt
    # is this fool proof? 
    # The number of bootstrap samples must be large enough to avoid spurious results.
    wt <- rev(wt)
    if (attr(object$wt, "bootWt")) {
      if (attr(object$wt, "bootWt.R") < 999) {
        stop("Restriktor ERROR: increase the number of bootstrap draws. Preferably to a large number e.g., bootWt.R = 99999")
      }
      wt.idx <- which(wt == 0)
      wt <- wt[-wt.idx]
    }
  
    pvalue <- con_pvalue_Fbar(wt          = wt, 
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
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
                                         Ts_org   = Ts, 
                                         type     = type, 
                                         test     = "LRT", 
                                         meq_alt  = meq_alt,
                                         R        = R, 
                                         p.distr  = p.distr,
                                         parallel = parallel,
                                         ncpus    = ncpus, 
                                         cl       = cl,
                                         seed     = seed, 
                                         verbose  = verbose)
  } else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based(object, 
                                          Ts_org   = Ts, 
                                          type     = type, 
                                          test     = "LRT",
                                          meq_alt  = meq_alt,
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
  if (!is.null(attr(type, "type_org"))) {
    type <- "global"
  }
  
  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq_alt     = meq_alt,
              iact        = object$iact,
              type        = type,
              test        = "LRT",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b_eqrestr   = b_eqrestr,
              b_unrestr   = b_unrestr,
              b_restr     = b_restr,
              b_restr_alt = b_restr_alt,
              Sigma       = Sigma,
              R2_org      = object$R2_org,
              R2_reduced  = object$R2_reduced,
              boot        = boot,
              model_org   = model_org)

  OUT <- list(OUT)
  names(OUT) <- type
  
  class(OUT) <- "conTest"

  OUT
}


# REF: Robertson, Wright and Dykstra (1988, p. 321). Order constrained statistical inference.
conTestScore.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                               p.distr = rnorm, parallel = "no", ncpus = 1L,
                               cl = NULL, seed = 1234, verbose = FALSE,
                               control = NULL, ...) {
  
  # rename for internal use
  meq_alt <- neq.alt
  
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
  model_org <- object$model_org
  # model matrix
  X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- vcov(model_org) 
  # sample size
  n <- dim(X)[1]
  # number of parameters
  p <- dim(X)[2]
  # weights
  w <- weights(model_org)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  #W <- diag(w)
  # parameter estimates
  b_unrestr <- object$b_unrestr
  b_restr <- object$b_restr
  b_eqrestr <- NULL
  b_restr_alt <- NULL
  # length parameter vector
  p <- length(b_unrestr)
  # variable names
  vnames <- names(b_unrestr)
  # restraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  #control
  control <- object$control
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  }
  
  if (is.null(control$tol)) {
    tol <- sqrt(.Machine$double.eps)
  } else {
    tol <- control$tol
  }
  
  # check for intercept                                          
  intercept <- any(attr(terms(model_org), "intercept"))
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
        attr(type, "type_org") <- "global"
    } else {
      # remove all rows with only zeros
      AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
      rAmatX <- GaussianElimination(t(AmatX), tol = tol)
      AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
    }
    AmatG <- rbind(AmatX, Amat)
    bvecG <- c(rep(0, nrow(AmatX)), bvec)
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
  }

  if (type == "global") {
    b_eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               b_unrestr = b_unrestr,
                               w         = w, 
                               Amat      = AmatG, 
                               bvec      = bvecG, 
                               meq       = nrow(AmatG),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$solution
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    
    res0 <- y - X %*% b_eqrestr
    res1 <- residuals(object)
    
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
    Ts <- (G0 - G1) %*% solve(I0, (G0 - G1))
    ###############################################
    # df0 <- n - (p - nrow(AmatG))
    # df1 <- n - (p - qr(Amat[0:meq,])$rank)
    # s20 <- sum((y - X %*% b_eqrestr)^2) / df0
    # s21 <- sum((y - X %*% b_restr)^2) / df1
    # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b_eqrestr)))
    # d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b_restr)))
    # I0 <- 1/s20 * (t(X) %*% X)
    # Ts2 <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
    ###############################################
  } else if (type == "A") {
    b_eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               b_unrestr = b_unrestr,
                               w         = w, 
                               Amat      = Amat,
                               bvec      = bvec, 
                               meq       = nrow(Amat),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$solution
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    
    res0 <- y - X %*% b_eqrestr
    res1 <- residuals(object)
    
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
    Ts <- (G0 - G1) %*% solve(I0, (G0 - G1))
    #############################################
    # df0 <- n - (p - nrow(Amat))
    # df1 <- n - (p - qr(Amat[0:meq,])$rank)
    # s20 <- sum((y - X %*% b_eqrestr)^2) / df0
    # s21 <- sum((y - X %*% b_restr)^2) / df1
    # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b_eqrestr)))
    # d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b_restr)))
    # I0 <- 1/s20 * (t(X) %*% X)
    # Ts2 <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
    ###############################################
  } else if (type == "B") {
    if (meq_alt == 0L) {
      res0 <- residuals(object)
      res1 <- residuals(model_org)
      
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
      Ts <- (G0 - G1) %*% solve(I0, (G0 - G1))
      
      #########
      # df0 <- n - (p - qr(Amat[0:meq,])$rank)
      # df1 <- n - p
      # s20 <- sum((y - X %*% b_restr)^2) / df0
      # s21 <- sum((y - X %*% b_unrestr)^2) / df1
      # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b_restr)))
      # d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b_unrestr)))
      # I0 <- 1/s20 * (t(X) %*% X)
      # Ts2 <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
      ###############################################
      } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq_alt > 0L && meq_alt <= meq) {
        b_restr_alt <- con_solver_lm(X         = X, 
                                  y         = y, 
                                  b_unrestr = b_unrestr,
                                  w         = w,
                                  Amat      = Amat[1:meq_alt,,drop=FALSE],
                                  bvec      = bvec[1:meq_alt], 
                                  meq       = meq_alt,
                                  absval    = ifelse(is.null(control$absval), 1e-09,
                                                     control$absval),
                                  maxit     = ifelse(is.null(control$maxit), 1e04,
                                                     control$maxit))$solution
        b_restr_alt[abs(b_restr_alt) < ifelse(is.null(control$tol),                                        
                                              sqrt(.Machine$double.eps),                                        
                                              control$tol)] <- 0L
        names(b_restr_alt) <- vnames
        
        res0 <- residuals(object)
        res1 <- y - X %*% b_restr_alt
        
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
        Ts <- (G0 - G1) %*% solve(I0, (G0 - G1))
        ########################
        # df1 <- df0 <- n - (p - qr(Amat[0:meq,])$rank)
        # s20 <- sum((y - X %*% b_restr)^2) / df0
        # s21 <- sum((y - X %*% b_restr_alt)^2) / df1
        # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b_restr)))
        # d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b_restr_alt)))
        # I0 <- 1/s20 * (t(X) %*% X)
        # Ts2 <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
        ###############################################
      }
      else {
      stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  if (!is.null(object$wt) && boot == "no") {
    wt <- object$wt
    wt <- rev(wt)
    # is this fool proof? 
    # The number of bootstrap samples must be large enough to avoid spurious results.
    if (attr(object$wt, "bootWt")) {
      if (attr(object$wt, "bootWt.R") < 999) {
        stop("Restriktor ERROR: increase the number of bootstrap draws. Preferably to a large number e.g., bootWt.R = 99999")
      }
      wt.idx <- which(wt == 0)
      wt <- wt[-wt.idx]
    }
  
    pvalue <- con_pvalue_Fbar(wt          = wt, 
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
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
                                         Ts_org   = Ts, 
                                         type     = type, 
                                         test     = "score",
                                         meq_alt  = meq_alt, 
                                         R        = R, 
                                         p.distr  = p.distr, 
                                         parallel = parallel,
                                         ncpus    = ncpus, 
                                         cl       = cl,
                                         seed     = seed, 
                                         verbose  = verbose)
  } else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based(object, 
                                          Ts_org   = Ts, 
                                          type     = type, 
                                          test     = "score",
                                          meq_alt  = meq_alt,
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
  if (!is.null(attr(type, "type_org"))) {
    type <- "global"
  }
  
  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq_alt     = meq_alt,
              iact        = object$iact,
              type        = type,
              test        = "Score",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b_eqrestr   = b_eqrestr,
              b_unrestr   = b_unrestr,
              b_restr     = b_restr,
              b_restr_alt = b_restr_alt,
              Sigma       = Sigma,
              R2_org      = object$R2_org,
              R2_reduced  = object$R2_reduced,
              boot        = boot,
              model_org   = model_org)
  
  OUT <- list(OUT)
  names(OUT) <- type
  
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
  
  b_unrestr   <- object$b_unrestr
  Sigma       <- object$Sigma
  df.residual <- object$df.residual
    
  Ts <- as.vector(min((Amat %*% b_unrestr - bvec) / 
                        sqrt(diag(Amat %*% Sigma %*% t(Amat))))) 
  
  pvalue <- 1 - pt(Ts, df.residual)

  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq_alt     = 0L,
              iact        = object$iact,
              type        = "C",
              test        = "t",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b_unrestr   = b_unrestr,
              b_restr     = object$b_restr,
              Sigma       = Sigma,
              R2_org      = ifelse(!is.null(object$R2_org), object$R2_org, as.numeric(NA)),
              R2_reduced  = ifelse(!is.null(object$R2_reduced), object$R2_reduced, as.numeric(NA)),
              boot        = "no",
              model_org   = object$model_org)
  
  OUT <- list(OUT)
    names(OUT) <- "C"
  
  class(OUT) <- "conTest"
  
  OUT
}
