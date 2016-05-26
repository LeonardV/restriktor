## computes the F, LRT and score test statistic ##
##
# REF: Silvapulle and Sen (2005). Constrained statistical inference. Chapter 2.
conTestF.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", B = 9999, 
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
  # variable names
  vnames <- names(b.unrestr)
  # constraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq

  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality restriktions only.")
  }
  
  if (type == "global") {
    # check for intercept
    intercept <- any(attr(terms(model.org), "intercept"))
    l <- length(b.restr)
    if (intercept) {
      bvecg <- bvec
      Amatg <- cbind(rep(0, (l - 1)), diag(rep(1, l - 1))) 
      Amatx <- Amatg %*% (diag(rep(1, l)) - t(Amat) %*% 
                                solve(Amat %*% t(Amat), Amat))
      if (!all(abs(Amatx) == 0)) {
        Amatx <- Amatx[!rowSums(abs(Amatx) == 0) == l, , drop = FALSE]
        if (nrow(Amatx) > 1) {
          Amat.rref <- GaussianElimination(Amatx)
          if (Amat.rref$rank == 1) {
            Amatx <- matrix(Amatx[1, ], 1, ncol(Amatx))
          } else {
            if (Amat.rref$rank < nrow(Amatx)) {
              Amatx <- Amatx[Amat.rref$pivot, , drop = FALSE]
            }
          }
        }
        Amatg <- rbind(Amatx, Amat)
        bvecg <- c(rep(0, nrow(Amatx)), bvec)
      } else {
        Amatg <- Amat
        bvecg <- bvec
      }
    } else {
        stop("Restriktor ERROR: test not ready yet for models without intercept.")      
    } 
    # call quadprog
    b.eqrestr <- con_solver(b.unrestr, 
                            X      = X, 
                            y      = y, 
                            w      = w, 
                            Amat   = Amatg,
                            bvec   = bvecg, 
                            meq    = nrow(Amatg),
                            absval = ifelse(is.null(control$absval), 1e-09, 
                                            control$absval),
                            maxit  = ifelse(is.null(control$maxit), 1e04, 
                                            control$maxit))$solution
    # fix estimates to zero 
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
      names(b.eqrestr) <- vnames
    # compute global test statistic
    Ts <- c(t(b.restr - b.eqrestr) %*% solve(Sigma, b.restr - b.eqrestr))
  } else if (type == "A") {
    b.eqrestr <- con_solver(b.unrestr, 
                            X      = X, 
                            y      = y, 
                            w      = w, 
                            Amat   = Amat,
                            bvec   = bvec, 
                            meq    = nrow(Amat),
                            absval = ifelse(is.null(control$absval), 1e-09, 
                                            control$absval),
                            maxit  = ifelse(is.null(control$maxit), 1e04, 
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
        b.restr.alt <- con_solver(b.unrestr, 
                                  X      = X, 
                                  y      = y, 
                                  w      = w, 
                                  Amat   = Amat[1:meq.alt,,drop = FALSE],
                                  bvec   = bvec[1:meq.alt], 
                                  meq    = meq.alt,
                                  absval = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                                  maxit  = ifelse(is.null(control$maxit), 1e04, 
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
                                         R        = B, 
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
conTestLRT.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", B = 9999, 
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
  # variable names
  vnames <- names(b.unrestr)
  # constraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality restriktions only.")
  }
  
  if (type == "global") {
    # check for intercept
    intercept <- any(attr(terms(model.org), "intercept"))
    l <- length(b.restr)
    if (intercept) {
      bvecg <- bvec
      Amatg <- cbind(rep(0, (l - 1)), diag(rep(1, l - 1))) 
      Amatx <- Amatg %*% (diag(rep(1, l)) - t(Amat) %*% 
                            solve(Amat %*% t(Amat), Amat))
      if (!all(abs(Amatx) == 0)) {
        Amatx <- Amatx[!rowSums(abs(Amatx) == 0) == l, , drop = FALSE]
        if (nrow(Amatx) > 1) {
          Amat.rref <- GaussianElimination(Amatx)
          if (Amat.rref$rank == 1) {
            Amatx <- matrix(Amatx[1, ], 1, ncol(Amatx))
          } else {
            if (Amat.rref$rank < nrow(Amatx)) {
              Amatx <- Amatx[Amat.rref$pivot, , drop = FALSE]
            }
          }
        }
        Amatg <- rbind(Amatx, Amat)
        bvecg <- c(rep(0, nrow(Amatx)), bvec)
      } else {
        Amatg <- Amat
        bvecg <- bvec
      }
    } else {
      stop("Restriktor ERROR: test not ready yet for models without intercept.")      
    } 
    
    b.eqrestr <- con_solver(b.unrestr, 
                            X      = X, 
                            y      = y, 
                            w      = w, 
                            Amat   = Amatg,
                            bvec   = bvecg, 
                            meq    = nrow(Amatg),
                            absval = ifelse(is.null(control$absval), 1e-09, 
                                            control$absval),
                            maxit  = ifelse(is.null(control$maxit), 1e04, 
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
  } else if (type == "A") {
    b.eqrestr <- con_solver(b.unrestr, 
                            X      = X, 
                            y      = y, 
                            w      = w, 
                            Amat   = Amat,
                            bvec   = bvec, 
                            meq    = nrow(Amat),
                            absval = ifelse(is.null(control$absval), 1e-09, 
                                            control$absval),
                            maxit  = ifelse(is.null(control$maxit), 1e04, 
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
        b.restr.alt <- con_solver(b.unrestr, 
                                  X      = X, 
                                  y      = y, 
                                  w      = w,
                                  Amat   = Amat[1:meq.alt,,drop=FALSE],
                                  bvec   = bvec[1:meq.alt], 
                                  meq    = meq.alt,
                                  absval = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                                  maxit  = ifelse(is.null(control$maxit), 1e04, 
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
                                         R        = B, 
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
                                          R        = B, 
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


# REF: Silvapulle, M.J. and Silvapulle, P. (1995). A score Test Against One-Sided Alternatives
# Journal of the American Statistical Association, Vol. 90, No. 429 (Mar., 1995), pp. 342-349
conTestScore.conLM <- function(object, type = "A", neq.alt = 0, boot = "no", B = 9999, 
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
  W <- diag(w)
  # parameter estimates
  b.unrestr <- object$b.unrestr
  b.restr <- object$b.restr
  b.eqrestr <- NULL
  b.restr.alt <- NULL
  # variable names
  vnames <- names(b.unrestr)
  # restraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality restriktions only.")
  }

  if (type == "global") {
    # check for intercept
    intercept <- any(attr(terms(model.org), "intercept"))
    l <- length(b.restr)
    if (intercept) {
      bvecg <- bvec
      Amatg <- cbind(rep(0, (l - 1)), diag(rep(1, l - 1))) 
      Amatx <- Amatg %*% (diag(rep(1, l)) - t(Amat) %*% 
                            solve(Amat %*% t(Amat), Amat))
      if (!all(abs(Amatx) == 0)) {
        Amatx <- Amatx[!rowSums(abs(Amatx) == 0) == l, , drop = FALSE]
        if (nrow(Amatx) > 1) {
          Amat.rref <- GaussianElimination(Amatx)
          if (Amat.rref$rank == 1) {
            Amatx <- matrix(Amatx[1, ], 1, ncol(Amatx))
          } else {
            if (Amat.rref$rank < nrow(Amatx)) {
              Amatx <- Amatx[Amat.rref$pivot, , drop = FALSE]
            }
          }
        }
        Amatg <- rbind(Amatx, Amat)
        bvecg <- c(rep(0, nrow(Amatx)), bvec)
      } else {
        Amatg <- Amat
        bvecg <- bvec
      }
    } else {
      stop("Restriktor ERROR: test not ready yet for models without intercept.")      
    } 
    
    b.eqrestr <- con_solver(b.unrestr, 
                            X      = X, 
                            y      = y, 
                            w      = w, 
                            Amat   = Amatg, 
                            bvec   = bvecg, 
                            meq    = nrow(Amatg),
                            absval = ifelse(is.null(control$absval), 1e-09, 
                                            control$absval),
                            maxit  = ifelse(is.null(control$maxit), 1e04, 
                                            control$maxit))$solution
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
      names(b.eqrestr) <- vnames
    
    df <- n-(p-nrow(Amatg)) 
    s20 <- sum(w*(y - X %*% b.eqrestr)^2) / df
    d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.eqrestr)))
    i <- 1/s20 * (t(X) %*% W %*% X)
    U <- 1/sqrt(n) * solve(i) %*% d0
    UI <- t(U) %*% i
    D <- i
    b <- solve.QP(Dmat = D, 
                  dvec = UI, 
                  Amat = t(Amat), 
                  bvec = bvec, 
                  meq  = meq)$solution
    Ts <- t(U) %*% i %*% U - ( t(U-b) %*% i %*% (U-b) ) 
    Ts <- as.numeric(n * Ts)
  } else if (type == "A") {
    b.eqrestr <- con_solver(b.unrestr, 
                            X      = X, 
                            y      = y, 
                            w      = w, 
                            Amat   = Amat,
                            bvec   = bvec, 
                            meq    = nrow(Amat),
                            absval = ifelse(is.null(control$absval), 1e-09, 
                                            control$absval),
                            maxit  = ifelse(is.null(control$maxit), 1e04, 
                                            control$maxit))$solution
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol),                                        
                                      sqrt(.Machine$double.eps),                                        
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    
    df <- n - (p - nrow(Amat))
    s20 <- sum(w*(y - X %*% b.eqrestr)^2) / df
    d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.eqrestr)))
    i <- 1/s20 * (t(X) %*% W %*% X)
    U <- 1/sqrt(n) * solve(i) %*% d0
    UI <- t(U) %*% i
    D <- i
    b <- solve.QP(Dmat = D, 
                  dvec = UI, 
                  Amat = t(Amat), 
                  bvec = bvec, 
                  meq  = meq)$solution
    Ts <- t(U) %*% i %*% U - ( t(U-b) %*% i %*% (U-b) ) 
    Ts <- as.numeric(n*Ts)
  } else if (type == "B") {
      if (meq.alt == 0L) {
        df <- n - (p - qr(Amat[0:meq,])$rank)
        s20 <- sum(w*(y - X %*% b.restr)^2) / df
        d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.restr)))
        i <- 1/s20 * (t(X) %*% W %*% X)
        U <- 1/sqrt(n) * solve(i) %*% d0
        UI <- t(U) %*% i
        Ts <- t(U) %*% i %*% U 
        Ts <- as.numeric(n*Ts)
      }
      else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        b.restr.alt <- con_solver(b.unrestr, 
                                  X      = X, 
                                  y      = y, 
                                  w      = w,
                                  Amat   = Amat[1:meq.alt,,drop=FALSE],
                                  bvec   = bvec[1:meq.alt], 
                                  meq    = meq.alt,
                                  absval = ifelse(is.null(control$absval), 1e-09,
                                                  control$absval),
                                  maxit  = ifelse(is.null(control$maxit), 1e04,
                                                  control$maxit))$solution
        b.restr.alt[abs(b.restr.alt) < ifelse(is.null(control$tol),                                        
                                              sqrt(.Machine$double.eps),                                        
                                              control$tol)] <- 0L
        names(b.restr.alt) <- vnames
        
        df <- n - (p - qr(Amat[0:meq,])$rank)
        s20 <- sum(w*(y - X %*% b.restr)^2) / df
        d0 <- 1/s20 * (t(X) %*% (w*(y - X %*% b.restr)))
        i <- 1/s20 * (t(X) %*% W %*% X)
        U <- 1/sqrt(n) * solve(i) %*% d0
        UI <- t(U) %*% i
        D <- i
        b <- solve.QP(Dmat = D, 
                      dvec = UI, 
                      Amat = t(Amat[1:meq.alt, ,drop = FALSE]), 
                      bvec = bvec[1:meq.alt], 
                      meq  = meq.alt)$solution
        Ts <- t(U) %*% i %*% U - ( t(U-b) %*% i %*% (U-b) ) 
        Ts <- as.numeric(n*Ts)
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
                                         R        = B, 
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
                                          R        = B, 
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
# REF: S. Sasabuchi (1980). A Test of a Multivariate Normal Mean with Composite 
# Hypotheses Determined by Linear Inequalities. Biometrika Trust, 67 (2), 429-439.
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
