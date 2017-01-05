conTestF.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                            p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                            cl = NULL, seed = 1234, verbose = FALSE,
                            control = NULL, ...) {
  
  # rename for internal use
  meq_alt <- neq.alt
  
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
  model_org <- object$model_org
  # orignal model call
  call_org <- as.list(model_org$call)
  # model matrix
  X <- model.matrix(model_org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # weights
  weights <- model_org$weights
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- object$Sigma
  # unconstrained scale
  scale <- model_org$s
  # parameter estimates
  b_unrestr <- object$b_unrestr
  b_restr <- object$b_restr
  b_eqrestr <- NULL
  b_restr.alt <- NULL
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
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
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
    call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = AmatG, bvec = bvecG, 
                             meq = nrow(AmatG), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b_eqrestr <- rfit$coefficients
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    Ts <- robustFm(x         = X, 
                   y         = y, 
                   b_unrestr = b_unrestr,
                   b_eqrestr = b_eqrestr, 
                   b_restr   = b_restr, 
                   scale     = scale, 
                   cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))$Ts
  } else if (type == "A") {
    call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = Amat, bvec = bvec, 
                             meq = nrow(Amat), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b_eqrestr <- rfit$coefficients
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    Ts <- robustFm(x         = X, 
                   y         = y, 
                   b_unrestr = b_unrestr,
                   b_eqrestr = b_eqrestr, 
                   b_restr   = b_restr, 
                   scale     = scale, 
                   cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))$Ts
  } else if (type == "B") {
    if (meq_alt == 0L) {
      Ts <- robustFm(x         = X, 
                     y         = y, 
                     b_unrestr = b_unrestr,
                     b_eqrestr = b_restr, 
                     b_restr   = b_unrestr, 
                     scale     = scale, 
                     cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))$Ts
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq_alt > 0L && meq_alt <= meq) {
        call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
        CALL <- c(call_org, list(x = X, y = y, weights = weights,
                                 psi = psi.bisquare, 
                                 Amat = Amat[1:meq_alt, , drop = FALSE], 
                                 meq = meq_alt, bvec = bvec[1:meq_alt],
                                 tol = tol))
        
        rfit <- do.call("conRLM_fit", CALL)
        b_restr.alt <- rfit$coefficients
        b_restr.alt[abs(b_restr.alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps), 
                                              control$tol)] <- 0L
        names(b_restr.alt) <- vnames
        Ts <- robustFm(x         = X, 
                       y         = y, 
                       b_unrestr = b_unrestr,
                       b_eqrestr = b_restr, 
                       b_restr   = b_restr.alt, 
                       scale     = scale, 
                       cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))$Ts
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
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts_org   = Ts, 
                                         type     = type, 
                                         test     = "F",
                                         meq_alt  = meq_alt,
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
              b_restr.alt = b_restr.alt,
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



conTestWald.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                               p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                               cl = NULL, seed = 1234, verbose = FALSE,
                               control = NULL, ...) {
  
  # rename for internal use
  meq_alt <- neq.alt
  
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
  model_org <- object$model_org
  # original model call
  call_org <- as.list(model_org$call)
  # model matrix
  X <- model.matrix(model_org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # weights
  weights <- model_org$weights
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- object$Sigma
  # unconstrained scale
  scale <- model_org$s
  # parameter estimates
  b_unrestr <- object$b_unrestr
  b_restr <- object$b_restr
  b_eqrestr <- NULL
  b_restr.alt <- NULL
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
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
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
    call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = AmatG, bvec = bvecG, 
                             meq = nrow(AmatG), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b_eqrestr <- rfit$coefficients
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol), 
                                       sqrt(.Machine$double.eps), 
                                       control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    out0 <- robustWaldScores(x         = X, 
                             y         = y,  
                             b_eqrestr = b_eqrestr, 
                             b_restr   = b_restr, 
                             b_unrestr = b_unrestr,
                             scale     = scale, 
                             test      = "Wald", 
                             cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))
    Ts <- out0$Ts
    Sigma <- out0$V
  } else if (type == "A") {
    call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = Amat, bvec = bvec, 
                             meq = nrow(Amat), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b_eqrestr <- rfit$coefficients
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    
    out1 <- robustWaldScores(x         = X, 
                             y         = y,  
                             b_eqrestr = b_eqrestr, 
                             b_restr   = b_restr, 
                             b_unrestr = b_unrestr,
                             scale     = scale, 
                             test      = "Wald", 
                             cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))
    Ts <- out1$Ts
    Sigma <- out1$V
  }
  else if (type == "B") {
    if (meq_alt == 0L) {
      out2 <- robustWaldScores(x         = X, 
                               y         = y, 
                               b_eqrestr = b_restr, 
                               b_restr   = b_unrestr,
                               b_unrestr = b_unrestr,
                               scale     = scale, 
                               test      = "Wald", 
                               cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))
      Ts <- out2$Ts
      Sigma <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq_alt > 0L && meq_alt <= meq) {
        call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
        CALL <- c(call_org, list(x = X, y = y, weights = weights,
                                 psi = psi.bisquare, 
                                 Amat = Amat[1:meq_alt,,drop=FALSE], 
                                 meq = meq_alt, 
                                 bvec = bvec[1:meq_alt], tol = tol))
        
        rfit <- do.call("conRLM_fit", CALL)
        b_restr.alt <- rfit$coefficients
        b_restr.alt[abs(b_restr.alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps), 
                                              control$tol)] <- 0L
        names(b_restr.alt) <- vnames
        out3 <- robustWaldScores(x         = X, 
                                 y         = y,  
                                 b_eqrestr = b_restr, 
                                 b_restr   = b_restr.alt,
                                 b_unrestr = b_unrestr,
                                 scale     = scale, 
                                 test      = "Wald", 
                                 cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))
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
      wt <- con_weights_boot(VCOV     = Sigma,
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
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts_org   = Ts, 
                                         type     = type, 
                                         test     = "Wald",
                                         meq_alt  = meq_alt,
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
                                          Ts_org   = Ts, 
                                          type     = type, 
                                          test     = "Wald",
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
              test        = "Wald",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b_eqrestr   = b_eqrestr,
              b_unrestr   = b_unrestr,
              b_restr     = b_restr,
              b_restr.alt = b_restr.alt,
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


#REF: Silvapulle, M.J. (1992) Robust Wald-Type Test of One-Sided Hypotheses in the linear model.
#American Statistical Association, 87, 417.
conTestWald2.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                                p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                                cl = NULL, seed = 1234, verbose = FALSE,
                                control = NULL, ...) {
  
  # rename for internal use
  meq_alt <- neq.alt
  
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
  model_org <- object$model_org
  # origianl model call
  call_org <- as.list(model_org$call)
  # model matrix
  X <- model.matrix(model_org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # weights
  weights <- model_org$weights
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- object$Sigma
  # unconstrained scale
  tau <- object$s2_unrestr
  # parameter estimates
  b_unrestr <- object$b_unrestr
  b_restr <- object$b_restr
  b_eqrestr <- NULL
  b_restr.alt <- NULL
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
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
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
    call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = AmatG, bvec = bvecG, 
                             meq = nrow(AmatG), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b_eqrestr <- rfit$coefficients
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    Ts <- robustWaldXX(x         = X, 
                       b_eqrestr = b_eqrestr, 
                       b_restr   = b_restr, 
                       b_unrestr = b_unrestr, 
                       tau       = tau)$Ts
  } else if (type == "A") {
    call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = Amat, bvec = bvec, 
                             meq = nrow(Amat), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b_eqrestr <- rfit$coefficients
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    Ts <- robustWaldXX(x         = X, 
                       b_eqrestr = b_eqrestr, 
                       b_restr   = b_restr, 
                       b_unrestr = b_unrestr, 
                       tau       = tau)$Ts
  } else if (type == "B") {
    if (meq_alt == 0L) {
      Ts <- robustWaldXX(x         = X, 
                         b_eqrestr = b_restr, 
                         b_restr   = b_unrestr, 
                         b_unrestr = b_unrestr, 
                         tau       = tau)$Ts
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq_alt > 0L && meq_alt <= meq) {
        call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
        CALL <- c(call_org, list(x = X, y = y, weights = weights,
                                 psi = psi.bisquare, 
                                 Amat = Amat[1:meq_alt, , drop = FALSE], 
                                 meq = meq_alt, bvec = bvec[1:meq_alt], 
                                 tol = tol))
        
        rfit <- do.call("conRLM_fit", CALL)
        b_restr.alt <- rfit$coefficients
        b_restr.alt[abs(b_restr.alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps), 
                                              control$tol)] <- 0L
        names(b_restr.alt) <- vnames
        Ts <- robustWaldXX(x         = X, 
                           b_eqrestr = b_restr, 
                           b_restr   = b_restr.alt, 
                           b_unrestr = b_restr.alt, 
                           tau       = tau)$Ts
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
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts_org   = Ts, 
                                         type     = type, 
                                         test     = "Wald2",
                                         meq_alt  = meq_alt,
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
                                          Ts_org   = Ts, 
                                          type     = type, 
                                          test     = "Wald2",
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
              test        = "Wald2",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b_eqrestr   = b_eqrestr,
              b_unrestr   = b_unrestr,
              b_restr     = b_restr,
              b_restr.alt = b_restr.alt,
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



conTestScore.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                                p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                                cl = NULL, seed = 1234, verbose = FALSE,
                                control = NULL, ...) {
  
  # rename for internal use
  meq_alt <- neq.alt
  
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
  model_org <- object$model_org
  # original model call
  call_org <- as.list(model_org$call)
  # model matrix
  X <- model.matrix(model_org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # weights
  weights <- model_org$weights
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- object$Sigma
  # unconstrained scale
  scale <- model_org$s
  # parameter estimates
  b_unrestr <- object$b_unrestr
  b_restr <- object$b_restr
  b_eqrestr <- NULL
  b_restr.alt <- NULL
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
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
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
    call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = AmatG, bvec = bvecG, 
                             meq = nrow(AmatG), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b_eqrestr <- rfit$coefficients
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    out0 <- robustWaldScores(x         = X, 
                             y         = y,  
                             b_eqrestr = b_eqrestr, 
                             b_restr   = b_restr, 
                             b_unrestr = b_unrestr,
                             scale     = scale, 
                             test      = "score", 
                             cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))
    Ts <- out0$Ts
    Sigma <- out0$V
  } else if (type == "A") {
    call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
                             psi = psi.bisquare, 
                             Amat = Amat, bvec = bvec, 
                             meq = nrow(Amat), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b_eqrestr <- rfit$coefficients
    b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol)] <- 0L
    names(b_eqrestr) <- vnames
    
    out1 <- robustWaldScores(x         = X, 
                             y         = y,  
                             b_eqrestr = b_eqrestr, 
                             b_restr   = b_restr, 
                             b_unrestr = b_unrestr,
                             scale     = scale, 
                             test      = "score", 
                             cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))
    Ts <- out1$Ts
    Sigma <- out1$V
  } else if (type == "B") {
    if (meq_alt == 0L) {
      out2 <- robustWaldScores(x         = X, 
                               y         = y,  
                               b_eqrestr = b_restr, 
                               b_restr   = b_unrestr,
                               b_unrestr = b_unrestr,
                               scale     = scale, 
                               test      = "score", 
                               cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))
      Ts <- out2$Ts
      Sigma <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq_alt > 0L && meq_alt <= meq) {
        call_org$weights <- call_org$psi <- call_org$X <- call_org$y <- NULL
        CALL <- c(call_org, list(x = X, y = y, weights = weights,
                                 psi = psi.bisquare, 
                                 Amat = Amat[1:meq_alt,,drop=FALSE], 
                                 meq = meq_alt, 
                                 bvec = bvec[1:meq_alt], tol = tol))
        
        rfit <- do.call("conRLM_fit", CALL)
        b_restr.alt <- rfit$coefficients
        b_restr.alt[abs(b_restr.alt) < ifelse(is.null(control$tol), 
                                              sqrt(.Machine$double.eps), 
                                              control$tol)] <- 0L
        names(b_restr.alt) <- vnames
        out3 <- robustWaldScores(x         = X, 
                                 y         = y,  
                                 b_eqrestr = b_restr, 
                                 b_restr   = b_restr.alt,
                                 b_unrestr = b_unrestr,
                                 scale     = scale, 
                                 test      = "score", 
                                 cc        = ifelse(is.null(call_org$c), 4.685061, call_org$c))
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
    wt <- con_weights_boot(VCOV     = Sigma,
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
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts_org   = Ts, 
                                         type     = type, 
                                         test     = "score",
                                         meq_alt  = meq_alt,
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
              b_restr.alt = b_restr.alt,
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
