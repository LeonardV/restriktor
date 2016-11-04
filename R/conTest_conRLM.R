conTestF.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                            p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                            cl = NULL, seed = 1234, verbose = FALSE,
                            control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conRLM" %in% class(object))) {
    stop("Restriktor ERROR: object must be of class conRLM.")
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
  # tukey's bisquare tuning constant
  cc <- model.org$call[["c"]]
  # model matrix
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- object$Sigma
  # unconstrained scale
  scale <- model.org$s
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
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
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
                          MASS::ginv(Amat %*% t(Amat)) %*% Amat)
    
    bvecG <- rep(0L, nrow(AmatG))
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "org_global") <- "org_global"
    }
    
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
    
    #fit inequality constrained robust model
    call.my <- list(Amat = AmatG, meq = nrow(AmatG), bvec = bvecG,
                    tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol))
    # collect all original model arguments and add constraints
    CALL <- c(list(model.org), call.my)
    CALL <- CALL[!duplicated(CALL)]
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    Ts <- robustFm(x         = X, 
                   y         = y, 
                   b.unrestr = b.unrestr,
                   b.eqrestr = b.eqrestr, 
                   b.restr   = b.restr, 
                   scale     = scale, 
                   cc        = ifelse(is.null(cc), 4.685061, cc))
  } else if (type == "A" | type == "Ax") {
    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec,
                    tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol))
    # collect all original model arguments and add constraints
    CALL <- c(list(model.org), call.my)
    CALL <- CALL[!duplicated(CALL)]
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    Ts <- robustFm(x         = X, 
                   y         = y, 
                   b.unrestr = b.unrestr,
                   b.eqrestr = b.eqrestr, 
                   b.restr   = b.restr, 
                   scale     = scale, 
                   cc        = ifelse(is.null(cc), 4.685061, cc))
  } else if (type == "B") {
    if (meq.alt == 0L) {
      Ts <- robustFm(x         = X, 
                     y         = y, 
                     b.unrestr = b.unrestr,
                     b.eqrestr = b.restr, 
                     b.restr   = b.unrestr, 
                     scale     = scale, 
                     cc        = ifelse(is.null(cc), 4.685061, cc))
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        call.my <- list(Amat = Amat[1:meq.alt, , drop = FALSE], 
                        meq = meq.alt, bvec = bvec[1:meq.alt],
                        tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                      control$tol))
        # collect all original model arguments and add constraints
        CALL <- c(list(model.org), call.my)
        CALL <- CALL[!duplicated(CALL)]
        rfit <- do.call("conRLM_fit", CALL)
        b.restr.alt <- rfit$coefficients
        b.restr.alt[abs(b.restr.alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps), 
                                              control$tol)] <- 0L
        names(b.restr.alt) <- vnames
        Ts <- robustFm(x         = X, 
                       y         = y, 
                       b.unrestr = b.unrestr,
                       b.eqrestr = b.restr, 
                       b.restr   = b.restr.alt, 
                       scale     = scale, 
                       cc        = ifelse(is.null(cc), 4.685061, cc))
      } else {
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
  } else {
    pvalue <- as.numeric(NA)
  }  
  
  OUT <- list(CON         = object$CON,
              type        = type,
              boot        = boot,
              b.eqrestr   = b.eqrestr,
              b.unrestr   = b.unrestr,
              b.restr     = b.restr,
              b.restr.alt = b.restr.alt,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq.alt     = meq.alt,
              iact        = object$iact,
              df.residual = df.residual,
              Sigma       = Sigma,
              Ts          = Ts,
              pvalue      = pvalue,
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}



conTestWald.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                               p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                               cl = NULL, seed = 1234, verbose = FALSE,
                               control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conRLM" %in% class(object))) {
    stop("Restriktor ERROR: object must be of class conRLM.")
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
  # tukey's bisquare tuning constant
  cc <- model.org$call[["c"]]
  # model matrix
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # weights
  #w <- weights(model.org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- object$Sigma
  # unconstrained scale
  scale <- model.org$s
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
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
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
                          MASS::ginv(Amat %*% t(Amat)) %*% Amat)
    
    bvecG <- rep(0L, nrow(AmatG))
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "org_global") <- "org_global"
    }
    
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
    
    #fit inequality constrained robust model
    call.my <- list(Amat = AmatG, meq = nrow(AmatG), bvec = bvecG,
                    tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol))
    # collect all original model arguments and add constraints
    CALL <- c(list(model.org), call.my)
    CALL <- CALL[!duplicated(CALL)]
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol), 
                                       sqrt(.Machine$double.eps), 
                                       control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    out0 <- robustWaldScores(x         = X, 
                             y         = y,  
                             b.eqrestr = b.eqrestr, 
                             b.restr   = b.restr, 
                             b.unrestr = b.unrestr,
                             scale     = scale, 
                             test      = "Wald", 
                             cc        = ifelse(is.null(cc), 4.685061, cc))
    Ts <- out0$Ts
    Sigma <- out0$V
  } else if (type == "A" | type == "Ax") {
    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec,
                    tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol))
    # collect all original model arguments and add constraints
    CALL <- c(list(model.org), call.my)
    CALL <- CALL[!duplicated(CALL)]
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    
    out1 <- robustWaldScores(x         = X, 
                             y         = y,  
                             b.eqrestr = b.eqrestr, 
                             b.restr   = b.restr, 
                             b.unrestr = b.unrestr,
                             scale     = scale, 
                             test      = "Wald", 
                             cc        = ifelse(is.null(cc), 4.685061, cc))
    Ts <- out1$Ts
    Sigma <- out1$V
  }
  else if (type == "B") {
    if (meq.alt == 0L) {
      out2 <- robustWaldScores(x         = X, 
                               y         = y, 
                               b.eqrestr = b.restr, 
                               b.restr   = b.unrestr,
                               b.unrestr = b.unrestr,
                               scale     = scale, 
                               test      = "Wald", 
                               cc        = ifelse(is.null(cc), 4.685061, cc))
      Ts <- out2$Ts
      Sigma <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        call.my <- list(Amat = Amat[1:meq.alt,,drop=FALSE], meq = meq.alt, 
                        bvec = bvec[1:meq.alt],
                        tol  = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                      control$tol))
        # collect all original model arguments and add constraints
        CALL <- c(list(model.org), call.my)
        CALL <- CALL[!duplicated(CALL)]
        rfit <- do.call("conRLM_fit", CALL)
        b.restr.alt <- rfit$coefficients
        b.restr.alt[abs(b.restr.alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps), 
                                              control$tol)] <- 0L
        names(b.restr.alt) <- vnames
        out3 <- robustWaldScores(x         = X, 
                                 y         = y,  
                                 b.eqrestr = b.restr, 
                                 b.restr   = b.restr.alt,
                                 b.unrestr = b.unrestr,
                                 scale     = scale, 
                                 test      = "Wald", 
                                 cc        = ifelse(is.null(cc), 4.685061, cc))
        Ts <- out3$Ts
        Sigma <- out3$V
      } else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  # We need to recalculate the weights based on V_hat = Sigma instead on solve(t(X)%*%X)
  # Do we have to? The differences look very small.
  ## compute mixing weights
  if (!is.null(object$wt)) {
    if (all(c(Amat) == 0)) { # unrestrikted case
      wt <- c(rep(0L, p), 1)
    } else if (attr(object$wt, "bootWt")) { # compute mixing weights based on simulation
      wt <- con_weightsBoot(VCOV     = Sigma,
                            Amat     = Amat,
                            meq      = meq,
                            R        = attr(object$wt, "bootWt.R"),
                            parallel = parallel,
                            ncpus    = ncpus,
                            cl       = cl,
                            seed     = seed,
                            verbose  = verbose)
    } else if (!attr(object$wt, "bootWt") && (meq < nrow(Amat))) { # compute mixing weights based on mvnorm
      wt <- rev(con_weights(Amat %*% Sigma %*% t(Amat), meq = meq))
    } else if (!attr(object$wt, "bootWt") && (meq == nrow(Amat))) { # only equality constraints
      wt <- rep(0L, ncol(Sigma) + 1)
      wt.idx <- ncol(Sigma) - meq + 1
      wt[wt.idx] <- 1
    }
  }
  
  
  #wt <- object$wt
  # is this fool proof? 
  # The number of bootstrap samples must be large enough to avoid spurious results.
  if (!is.null(object$wt) && boot == "no") {
    wt <- rev(wt)
    if (attr(object$wt, "bootWt")) {
      if (attr(object$wt, "bootWt.R") < 999) {
        stop("Restriktor ERROR: increase the number of bootstrap draws. Preferably to a large number e.g., bootWt.R = 99999")
      }
      wt.idx <- which(wt == 0)
      wt <- wt[-wt.idx]
    }
  
    # compute pvalue based on F-distribution
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
                                         test     = "Wald",
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
                                          test     = "Wald",
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
  
  OUT <- list(CON         = object$CON,
              type        = type,
              boot        = boot,
              b.eqrestr   = b.eqrestr,
              b.unrestr   = b.unrestr,
              b.restr     = b.restr,
              b.restr.alt = b.restr.alt,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq.alt     = meq.alt,
              iact        = object$iact,
              df.residual = df.residual,
              Sigma       = Sigma,
              Ts          = Ts,
              pvalue      = pvalue,
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}


#REF: Silvapulle, M.J. (1992) Robust Wald-Type Test of One-Sided Hypotheses in the linear model.
#American Statistical Association, 87, 417.
conTestWald2.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                                p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                                cl = NULL, seed = 1234, verbose = FALSE,
                                control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conRLM" %in% class(object))) {
    stop("Restriktor ERROR: object must be of class conRLM.")
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
  # tukey's bisquare tuning constant
#  cc <- model.org$call[["c"]]
  # model matrix
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
#  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- object$Sigma
  # unconstrained scale
#  scale <- model.org$s
  tau <- object$s2.unc
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
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
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
                          MASS::ginv(Amat %*% t(Amat)) %*% Amat)
    
    bvecG <- rep(0L, nrow(AmatG))
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "org_global") <- "org_global"
    }
    
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
    
    #fit inequality constrained robust model
    call.my <- list(Amat = AmatG, meq = nrow(AmatG), bvec = bvecG,
                    tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol))
    # collect all original model arguments and add constraints
    CALL <- c(list(model.org), call.my)
    CALL <- CALL[!duplicated(CALL)]
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    Ts <- robustWaldXX(x         = X, 
                       b.eqrestr = b.eqrestr, 
                       b.restr   = b.restr, 
                       b.unrestr = b.unrestr, 
                       tau       = tau)
  } else if (type == "A" | type == "Ax") {
    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec,
                    tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol))
    # collect all original model arguments and add constraints
    CALL <- c(list(model.org), call.my)
    CALL <- CALL[!duplicated(CALL)]
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    Ts <- robustWaldXX(x         = X, 
                       b.eqrestr = b.eqrestr, 
                       b.restr   = b.restr, 
                       b.unrestr = b.unrestr, 
                       tau       = tau)
  } else if (type == "B") {
    if (meq.alt == 0L) {
      Ts <- robustWaldXX(x         = X, 
                         b.eqrestr = b.restr, 
                         b.restr   = b.unrestr, 
                         b.unrestr = b.unrestr, 
                         tau       = tau)
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        call.my <- list(Amat = Amat[1:meq.alt, , drop = FALSE], 
                        meq = meq.alt, bvec = bvec[1:meq.alt],
                        tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                      control$tol))
        # collect all original model arguments and add constraints
        CALL <- c(list(model.org), call.my)
        CALL <- CALL[!duplicated(CALL)]
        rfit <- do.call("conRLM_fit", CALL)
        b.restr.alt <- rfit$coefficients
        b.restr.alt[abs(b.restr.alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps), 
                                              control$tol)] <- 0L
        names(b.restr.alt) <- vnames
        Ts <- robustWaldXX(x         = X, 
                           b.eqrestr = b.restr, 
                           b.restr   = b.restr.alt, 
                           b.unrestr = b.restr.alt, 
                           tau       = tau)
      } else {
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
  } else {
    pvalue <- as.numeric(NA)
  } 
  
  OUT <- list(CON         = object$CON,
              type        = type,
              boot        = boot,
              b.eqrestr   = b.eqrestr,
              b.unrestr   = b.unrestr,
              b.restr     = b.restr,
              b.restr.alt = b.restr.alt,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq.alt     = meq.alt,
              iact        = object$iact,
              df.residual = df.residual,
              Sigma       = Sigma,
              Ts          = Ts,
              pvalue      = pvalue,
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}



conTestScore.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                                p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                                cl = NULL, seed = 1234, verbose = FALSE,
                                control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conRLM" %in% class(object))) {
    stop("Restriktor ERROR: object must be of class conRLM.")
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
  # tukey's bisquare tuning constant
  cc <- model.org$call[["c"]]
  # model matrix
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # weights
  #w <- weights(model.org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- object$Sigma
  # unconstrained scale
  scale <- model.org$s
  # parameter estimates
  b.unrestr <- object$b.unrestr
  b.restr <- object$b.restr
  b.eqrestr <- NULL
  b.restr.alt <- NULL
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
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
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
                          MASS::ginv(Amat %*% t(Amat)) %*% Amat)
    
    bvecG <- rep(0L, nrow(AmatG))
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "org_global") <- "org_global"
    }
    
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
    
    #fit inequality constrained robust model
    call.my <- list(Amat = AmatG, meq = nrow(AmatG), bvec = bvecG,
                    tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol))
    # collect all original model arguments and add constraints
    CALL <- c(list(model.org), call.my)
    CALL <- CALL[!duplicated(CALL)]
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    out0 <- robustWaldScores(x         = X, 
                             y         = y,  
                             b.eqrestr = b.eqrestr, 
                             b.restr   = b.restr, 
                             b.unrestr = b.unrestr,
                             scale     = scale, 
                             test      = "score", 
                             cc        = ifelse(is.null(cc), 4.685061, cc))    
    Ts <- out0$Ts
    Sigma <- out0$V
  } else if (type == "A" | type == "Ax") {
    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec,
                    tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol))
    # collect all original model arguments and add constraints
    CALL <- c(list(model.org), call.my)
    CALL <- CALL[!duplicated(CALL)]
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b.eqrestr) <- vnames
    
    out1 <- robustWaldScores(x         = X, 
                             y         = y,  
                             b.eqrestr = b.eqrestr, 
                             b.restr   = b.restr, 
                             b.unrestr = b.unrestr,
                             scale     = scale, 
                             test      = "score", 
                             cc        = ifelse(is.null(cc), 4.685061, cc))
    Ts <- out1$Ts
    Sigma <- out1$V
  } else if (type == "B") {
    if (meq.alt == 0L) {
      out2 <- robustWaldScores(x         = X, 
                               y         = y,  
                               b.eqrestr = b.restr, 
                               b.restr   = b.unrestr,
                               b.unrestr = b.unrestr,
                               scale     = scale, 
                               test      = "score", 
                               cc        = ifelse(is.null(cc), 4.685061, cc))
      Ts <- out2$Ts
      Sigma <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        call.my <- list(Amat = Amat[1:meq.alt,,drop=FALSE], meq = meq.alt, 
                        bvec = bvec[1:meq.alt],
                        tol = ifelse (is.null(control$tol), sqrt(.Machine$double.eps), 
                                      control$tol))
        # collect all original model arguments and add constraints
        CALL <- c(list(model.org), call.my)
        CALL <- CALL[!duplicated(CALL)]
        rfit <- do.call("conRLM_fit", CALL)
        b.restr.alt <- rfit$coefficients
        b.restr.alt[abs(b.restr.alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps), 
                                              control$tol)] <- 0L
        names(b.restr.alt) <- vnames
        out3 <- robustWaldScores(x         = X, 
                                 y         = y,  
                                 b.eqrestr = b.restr, 
                                 b.restr   = b.restr.alt,
                                 b.unrestr = b.unrestr,
                                 scale     = scale, 
                                 test      = "score", 
                                 cc        = ifelse(is.null(cc), 4.685061, cc))
        Ts <- out3$Ts
        Sigma <- out3$V
      } else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  # We need to recalculate the weights based on V_hat = Sigam instead on solve(t(X)%*%X)
  # Do we have to? The differences are very small.
  ## compute mixing weights
  if (all(c(Amat) == 0)) { # unrestrikted case
    wt <- c(rep(0L, p), 1)
  } else if (attr(object$wt, "bootWt")) { # compute mixing weights based on simulation
    wt <- con_weightsBoot(VCOV     = Sigma,
                          Amat     = Amat,
                          meq      = meq,
                          R        = attr(object$wt, "bootWt.R"),
                          parallel = parallel,
                          ncpus    = ncpus,
                          cl       = cl,
                          seed     = seed,
                          verbose  = verbose)
  } else if (!attr(object$wt, "bootWt") & (meq < nrow(Amat))) { # compute mixing weights based on mvnorm
    wt <- rev(con_weights(Amat %*% Sigma %*% t(Amat), meq = meq))
  } else if (!attr(object$wt, "bootWt") & (meq == nrow(Amat))) { # only equality constraints
    wt <- rep(0L, ncol(Sigma) + 1)
    wt.idx <- ncol(Sigma) - meq + 1
    wt[wt.idx] <- 1
  }
  
  #wt <- object$wt
  # is this fool proof? 
  # The number of bootstrap samples must be large enough to avoid spurious results.
  if (!is.null(object$wt) && boot == "no") {
    wt <- rev(wt)
    if (attr(object$wt, "bootWt")) {
      if (attr(object$wt, "bootWt.R") < 999) {
        stop("Restriktor ERROR: increase the number of bootstrap draws. Preferably to a large number e.g., bootWt.R = 99999")
      }
      wt.idx <- which(wt == 0)
      wt <- wt[-wt.idx]
    }
    
    # compute pvalue based on F-distribution
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
  } else {
    pvalue <- as.numeric(NA)
  }  
  
  OUT <- list(CON         = object$CON,
              type        = type,
              boot        = boot,
              b.eqrestr   = b.eqrestr,
              b.unrestr   = b.unrestr,
              b.restr     = b.restr,
              b.restr.alt = b.restr.alt,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq.alt     = meq.alt,
              iact        = object$iact,
              df.residual = df.residual,
              Sigma       = Sigma,
              Ts          = Ts,
              pvalue      = pvalue,
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}
