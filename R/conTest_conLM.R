## computes the F, LRT and score test statistic ##
# REFs: 
# Silvapulle and Sen (2005). Constrained statistical inference. Chapter 2.
# Wolak, F. An exact test for multiple inequality and equality constraints in the linear regression model Journal of the American statistical association, 1987, 82, 782-793
conTestF.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                           p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                           cl = NULL, seed = 1234, verbose = FALSE,
                           control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conLM" %in% class(object))) {
    stop("object must be of class conLM.")
  }
  if (type != "global") {
    type <- toupper(type)
  }    
  if(!(type %in% c("A","B","global"))) {
    stop("type must be \"A\", \"B\", or \"global\"")
  }
  if(!(boot %in% c("no","residual","model.based","parametric","mix.weights"))) {
    stop("ERROR: boot method unknown.")
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
  #control
  control <- object$control
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality restriktions only.")
  }
  
  if (is.null(control$tol)) {
    tol <- sqrt(.Machine$double.eps)
  } else {
    tol <- control$tol
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
                          solve(Amat %*% t(Amat), Amat))
    
    bvecG <- rep(0L, nrow(AmatG))
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "org_global") <- "org_global"
    }
    
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
    
    # call quadprog
    b.eqrestr <- con_solver(X         = X, 
                            y         = y, 
                            b.unrestr = b.unrestr,
                            w         = w, 
                            Amat      = AmatG,
                            bvec      = bvecG, 
                            meq       = nrow(AmatG),
                            absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                            maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
  
    # fix estimates < tol to zero 
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    # compute global test statistic
    Ts <- c(t(b.restr - b.eqrestr) %*% solve(Sigma, b.restr - b.eqrestr))
  } else if (type == "A") {
    b.eqrestr <- con_solver(X         = X, 
                            y         = y, 
                            b.unrestr = b.unrestr,
                            w         = w, 
                            Amat      = Amat,
                            bvec      = bvec, 
                            meq       = nrow(Amat),
                            absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                            maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    # compute test statistic for hypothesis test type A
    Ts <- c(t(b.restr - b.eqrestr) %*% solve(Sigma, b.restr - b.eqrestr))
  } else if (type == "B") {
    if (meq.alt == 0L) {
      # compute test statistic for hypothesis test type B when no equalities are
      # preserved in the alternative hypothesis.
      Ts <- as.vector(t(b.unrestr - b.restr) %*% solve(Sigma, b.unrestr - b.restr))
    } else {
      if (meq.alt != 0L && meq.alt <= meq) {
        b.restr.alt <- con_solver(X         = X, 
                                  y         = y, 
                                  b.unrestr = b.unrestr,
                                  w         = w, 
                                  Amat      = Amat[1:meq.alt,,drop = FALSE],
                                  bvec      = bvec[1:meq.alt], 
                                  meq       = meq.alt,
                                  absval    = ifelse(is.null(control$absval), 1e-09, 
                                                     control$absval),
                                  maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                     control$maxit))$solution
        b.restr.alt[abs(b.restr.alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps),                                        
                                              control$tol)] <- 0L
        names(b.restr.alt) <- vnames
        # compute test statistic for hypothesis test type B when some equalities may 
        # be preserved in the alternative hypothesis.
        Ts <- as.vector(t(b.restr - b.restr.alt) %*% solve(Sigma, b.restr - b.restr.alt))
      } else {
        stop("neq.alt must not be larger than neq.")
      }
    }
  } 

  # The test statistics based on inequality constraints are often
  # distributed as mixtures of chi-squares. These mixing weights can be computed
  # using the multivariate normal distribution or via bootstrapping. The
  # pvalue can also be computed directly via the parametric bootstrap or model
  # based bootstrap, without fist computing the mixing weights.
  
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
  
  if (boot == "no") {
    pvalue <- con_pvalue_Fbar(wt          = wt, 
                              Ts.org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq.alt     = meq.alt)
   } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts.org   = Ts, 
                                         type     = type, 
                                         test     = "F", 
                                         meq.alt  = meq.alt, 
                                         R        = R, 
                                         p.distr  = p.distr,
                                         df       = df, 
                                         parallel = parallel,
                                         ncpus    = ncpus, 
                                         cl       = cl,
                                         seed     = seed, 
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
                                          verbose  = verbose)
  } 
  
  OUT <- list(CON = object$CON,
              type = type,
              boot = boot,
              b.eqrestr = b.eqrestr,
              b.unrestr = b.unrestr,
              b.restr = b.restr,
              b.restr.alt = b.restr.alt,
              Amat = Amat,
              bvec = bvec,
              meq = meq,
              meq.alt = meq.alt,
              iact = object$iact,
              df.residual = df.residual,
              Sigma = Sigma,
              Ts = Ts,
              pvalue = pvalue,
              model.org = model.org)

  class(OUT) <- "conTest"

  OUT
}


# REF: Silvapulle and Sen (2005). Constrained statistical inference. Chapter 3.
# Wolak, F. An exact test for multiple inequality and equality constraints in the linear regression model Journal of the American statistical association, 1987, 82, 782-793
conTestLRT.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                             p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                             cl = NULL, seed = 1234, verbose = FALSE,
                             control = NULL, ...) {

  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conLM" %in% class(object))) {
    stop("object must be of class conLM.")
  }
  if (type != "global") {
    type <- toupper(type)
  }  
  if(!(type %in% c("A","B","global"))) {
    stop("type must be \"A\", \"B\", or \"global\"")
  }
  if(!(boot %in% c("no", "residual", "model.based", "parametric", "mix.weights"))) {
    stop("ERROR: boot method unknown.")
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
  #control
  control <- object$control
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality restriktions only.")
  }
  
  if (is.null(control$tol)) {
    tol <- sqrt(.Machine$double.eps)
  } else {
    tol <- control$tol
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
                          solve(Amat %*% t(Amat), Amat))
    
    bvecG <- rep(0L, nrow(AmatG))
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "org_global") <- "org_global"
    }
    
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
    
    b.eqrestr <- con_solver(X         = X, 
                            y         = y, 
                            b.unrestr = b.unrestr, 
                            w         = w, 
                            Amat      = AmatG,
                            bvec      = bvecG, 
                            meq       = nrow(AmatG),
                            absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                            maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    ll0.out <- con_loglik_lm(X = X, 
                             y = y, 
                             b = b.eqrestr, 
                             w = w)
    ll0 <- ll0.out$loglik
    
    ll1 <- object$loglik
    Ts <- -2*(ll0 - ll1)
  } else if (type == "A" | type == "Ax") {
    b.eqrestr <- con_solver(X         = X, 
                            y         = y, 
                            b.unrestr = b.unrestr,
                            w         = w, 
                            Amat      = Amat,
                            bvec      = bvec, 
                            meq       = nrow(Amat),
                            absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                            maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    ll0.out <- con_loglik_lm(X = X, 
                             y = y, 
                             b = b.eqrestr, 
                             w = w)
    ll0 <- ll0.out$loglik
    
    ll1 <- object$loglik
    Ts <- -2*(ll0 - ll1)
  } else if (type == "B") {
      if (meq.alt == 0L) {
        ll0 <- object$loglik
        ll1.out <- con_loglik_lm(X = X, 
                                 y = y, 
                                 b = b.unrestr, 
                                 w = w)
        ll1 <- ll1.out$loglik
        Ts <- -2*(ll0 - ll1)
      }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        b.restr.alt <- con_solver(X         = X, 
                                  y         = y, 
                                  b.unrestr = b.unrestr,
                                  w         = w,
                                  Amat      = Amat[1:meq.alt,,drop=FALSE],
                                  bvec      = bvec[1:meq.alt], 
                                  meq       = meq.alt,
                                  absval    = ifelse(is.null(control$absval), 1e-09, 
                                                     control$absval),
                                  maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                     control$maxit))$solution
        b.restr.alt[abs(b.restr.alt) < ifelse(is.null(control$tol),                                        
                                              sqrt(.Machine$double.eps),                                        
                                              control$tol)] <- 0L
        names(b.restr.alt) <- vnames

        ll0 <- object$loglik
        ll1.out <- con_loglik_lm(X = X, 
                                 y = y, 
                                 b = b.restr.alt, 
                                 w = w)
        ll1 <- ll1.out$loglik
        Ts <- -2*(ll0 - ll1)
      }
      else {
        stop("neq.alt must not be larger than neq.")
      }
    }
  } 
  
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
  
  if (boot == "no") {
    pvalue <- con_pvalue_Fbar(wt          = wt, 
                              Ts.org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq.alt     = meq.alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts.org   = Ts, 
                                         type     = type, 
                                         test     = "LRT", 
                                         meq.alt  = meq.alt,
                                         R        = R, 
                                         p.distr  = p.distr,
                                         df       = df, 
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
  } 
  
  OUT <- list(CON = object$CON,
              type = type,
              boot = boot,
              b.eqrestr = b.eqrestr,
              b.unrestr = b.unrestr,
              b.restr = b.restr,
              b.restr.alt = b.restr.alt,
              Amat = Amat,
              bvec = bvec,
              meq = meq,
              meq.alt = meq.alt,
              iact = object$iact,
              df.residual = df.residual,
              Sigma = Sigma,
              Ts = Ts,
              pvalue = pvalue,
              model.org = object$model.org)

  class(OUT) <- "conTest"

  OUT
}


# REF: Robertson, Wright and Dykstra (1988, p. 321). Order constrained statistical inference.
conTestScore.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                               p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                               cl = NULL, seed = 1234, verbose = FALSE,
                               control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conLM" %in% class(object))) {
    stop("object must be of class conLM.")
  }
  if (type != "global") {
    type <- toupper(type)
  }  
  if(!(type %in% c("A","B","global"))) {
    stop("type must be \"A\", \"B\", or \"global\"")
  }
  if(!(boot %in% c("no", "residual", "model.based", "parametric", "mix.weights"))) {
    stop("ERROR: boot method unknown.")
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
  #control
  control <- object$control
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality restriktions only.")
  }
  
  if (is.null(control$tol)) {
    tol <- sqrt(.Machine$double.eps)
  } else {
    tol <- control$tol
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
                          solve(Amat %*% t(Amat), Amat))
    
    bvecG <- rep(0L, nrow(AmatG))
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "org_global") <- "org_global"
    }
    
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG

    b.eqrestr <- con_solver(X         = X, 
                            y         = y, 
                            b.unrestr = b.unrestr,
                            w         = w, 
                            Amat      = AmatG, 
                            bvec      = bvecG, 
                            meq       = nrow(AmatG),
                            absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                            maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
      names(b.eqrestr) <- vnames
    
    # df0 <- n-(p-nrow(AmatG)) 
    # s20 <- sum(w*(y - X %*% b.eqrestr)^2) / df0
    # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.eqrestr)))
    # i <- 1/s20 * (t(X) %*% W %*% X) / n
    # U <- 1/sqrt(n) * solve(i) %*% d0
    # UI <- t(U) %*% i
    # D <- i
    # b <- solve.QP(Dmat = D, 
    #               dvec = UI, 
    #               Amat = t(Amat), 
    #               bvec = bvec, 
    #               meq  = meq)$solution
    # Ts <- t(U) %*% i %*% U - ( t(U-b) %*% i %*% (U-b) ) 
    ###############################################
    df0 <- n - (p - nrow(AmatG))   
    df1 <- n - (p - qr(Amat[0:meq,])$rank)
    s20 <- sum((y - X %*% b.eqrestr)^2) / df0
    s21 <- sum((y - X %*% b.restr)^2) / df1
    d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.eqrestr)))
    d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b.restr)))
    I0 <- 1/s20 * (t(X) %*% X)
    Ts <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
    ###############################################
  } else if (type == "A" | type == "Ax") {
    b.eqrestr <- con_solver(X         = X, 
                            y         = y, 
                            b.unrestr = b.unrestr,
                            w         = w, 
                            Amat      = Amat,
                            bvec      = bvec, 
                            meq       = nrow(Amat),
                            absval    = ifelse(is.null(control$absval), 1e-09, 
                                               control$absval),
                            maxit     = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    
    ######################################################
    # df0 <- n - (p - nrow(Amat))                                                  
    # s20 <- sum(w*(y - X %*% b.eqrestr)^2) / df0
    #d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.eqrestr)))
    # i <- 1/s20 * (t(X) %*% W %*% X) / n
    # U <- 1/sqrt(n) * solve(i) %*% d0
    # UI <- t(U) %*% i
    # D <- i
    # b <- solve.QP(Dmat = D,
    #               dvec = UI,
    #               Amat = t(Amat),
    #               bvec = bvec,
    #               meq  = meq)$solution
    # Ts <- t(U) %*% i %*% U - ( t(U-b) %*% i %*% (U-b) )
    ###############################################
    df0 <- n - (p - nrow(Amat))   
    df1 <- n - (p - qr(Amat[0:meq,])$rank)
    s20 <- sum((y - X %*% b.eqrestr)^2) / df0
    s21 <- sum((y - X %*% b.restr)^2) / df1
    d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.eqrestr)))
    d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b.restr)))
    I0 <- 1/s20 * (t(X) %*% X)
    Ts <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
    ###############################################
  } else if (type == "B") {
      if (meq.alt == 0L) {
        # df0 <- n - (p - qr(Amat[0:meq,])$rank)
        # s20 <- sum(w*(y - X %*% b.restr)^2) / df0
        # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.restr)))
        # i <- 1/s20 * (t(X) %*% W %*% X) / n
        # U <- 1/sqrt(n) * solve(i) %*% d0
        # UI <- t(U) %*% i
        # Ts <- t(U) %*% i %*% U 
        ###############################################
        df0 <- n - (p - qr(Amat[0:meq,])$rank)
        df1 <- n - p   
        s20 <- sum((y - X %*% b.restr)^2) / df0
        s21 <- sum((y - X %*% b.unrestr)^2) / df1
        d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.restr)))
        d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b.unrestr)))
        I0 <- 1/s20 * (t(X) %*% X)
        Ts <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
        ###############################################
      }
      else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        b.restr.alt <- con_solver(X         = X, 
                                  y         = y, 
                                  b.unrestr = b.unrestr,
                                  w         = w,
                                  Amat      = Amat[1:meq.alt,,drop=FALSE],
                                  bvec      = bvec[1:meq.alt], 
                                  meq       = meq.alt,
                                  absval    = ifelse(is.null(control$absval), 1e-09,
                                                     control$absval),
                                  maxit     = ifelse(is.null(control$maxit), 1e04,
                                                     control$maxit))$solution
        b.restr.alt[abs(b.restr.alt) < ifelse(is.null(control$tol),                                        
                                              sqrt(.Machine$double.eps),                                        
                                              control$tol)] <- 0L
        names(b.restr.alt) <- vnames
        
        # df0 <- n - (p - qr(Amat[0:meq,])$rank)
        # s20 <- sum(w*(y - X %*% b.restr)^2) / df0
        # d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.restr)))
        # i <- 1/s20 * (t(X) %*% W %*% X) / n
        # U <- 1/sqrt(n) * solve(i) %*% d0
        # UI <- t(U) %*% i
        # D <- i
        # b <- solve.QP(Dmat = D, 
        #               dvec = UI, 
        #               Amat = t(Amat[1:meq.alt, ,drop = FALSE]), 
        #               bvec = bvec[1:meq.alt], 
        #               meq  = meq.alt)$solution
        # Ts <- t(U) %*% i %*% U - ( t(U-b) %*% i %*% (U-b) ) 
        ###############################################
        df1 <- df0 <- n - (p - qr(Amat[0:meq,])$rank)
        s20 <- sum((y - X %*% b.restr)^2) / df0
        s21 <- sum((y - X %*% b.restr.alt)^2) / df1
        d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.restr)))
        d1 <- 1/s21 * (t(X) %*% (w*(y - X %*% b.restr.alt)))
        I0 <- 1/s20 * (t(X) %*% X)
        Ts <- t(c(d0 - d1)) %*% solve(I0) %*% c(d0 - d1)
        ###############################################
      }
      else {
      stop("neq.alt must not be larger than neq.")
      }
    }
  } 
  
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
  
  if (boot == "no") {
    pvalue <- con_pvalue_Fbar(wt          = wt, 
                              Ts.org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq.alt     = meq.alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts.org   = Ts, 
                                         type     = type, 
                                         test     = "score",
                                         meq.alt  = meq.alt, 
                                         R        = R, 
                                         p.distr  = p.distr, 
                                         df       = df,
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
  } 
  
  OUT <- list(CON = object$CON,
              type = type,
              boot = boot,
              b.eqrestr = b.eqrestr,
              b.unrestr = b.unrestr,
              b.restr = b.restr,
              b.restr.alt = b.restr.alt,
              Amat = Amat,
              bvec = bvec,
              meq = meq,
              meq.alt = meq.alt,
              iact = object$iact,
              df.residual = df.residual,
              Sigma = Sigma,
              Ts = Ts,
              pvalue = pvalue,
              model.org = object$model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}



# hypothesis test type C is based on a t-distribution.
# intersection-union test 
# REF: S. Sasabuchi (1980). A Test of a Multivariate Normal Mean with Composite Hypotheses Determined by Linear Inequalities. Biometrika Trust, 67 (2), 429-439.
conTestC.conLM <- function(object, type = "C", ...) {
  
  if (!("conLM" %in% class(object))) {
    stop("object must be of class conLM.")
  }
  
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  model.org <- object$model.org
  Sigma <- vcov(model.org) 
  df.residual <- object$df.residual
  b.unrestr <- object$b.unrestr
    
  if (meq == 0L) {
    Ts <- as.vector(min((Amat %*% b.unrestr - bvec) / 
                          sqrt(diag(Amat %*% Sigma %*% t(Amat)))))
    pvalue <- 1 - pt(Ts, df.residual)
  } else {
    stop("test not applicable with equality restriktions.")
  }
  
  OUT <- list(type = "C",
              Ts = Ts, 
              pvalue = pvalue,
              b.unrestr = b.unrestr,
              Amat = Amat,
              bvec = bvec,
              meq = meq)
  
  class(OUT) <- "conTest"
  
  OUT
}
