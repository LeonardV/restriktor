### computes the F-bar, Wald(2)-bar and score-bar test statistic ####
conTestF.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                            p.distr = rnorm, parallel = "no", ncpus = 1L,
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
  # tolerance
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
  }
  
  call_org$formula <- call_org$data <- call_org$weights <- 
    call_org$subset <- call_org$na.action <- call_org$model <- 
    call_org$x.ret <- call_org$y.ret <- call_org$contrasts <- 
    call_org$X <- call_org$y <- NULL
    
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
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
        CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
  
  if (!(attr(object$wt, "method") == "none") && boot == "no") {
    wt <- object$wt
    pvalue <- con_pvalue_Fbar(wt          = rev(wt), 
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
    attr(pvalue, "wt") <- wt
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
                               p.distr = rnorm, parallel = "no", ncpus = 1L,
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
  # tolerance
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
  }
  
  call_org$formula <- call_org$data <- call_org$weights <- 
    call_org$subset <- call_org$na.action <- call_org$model <- 
    call_org$x.ret <- call_org$y.ret <- call_org$contrasts <- 
    call_org$X <- call_org$y <- NULL
  
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
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
    V <- out0$V
  } else if (type == "A") {
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
    V <- out1$V
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
      V <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq_alt > 0L && meq_alt <= meq) {
        CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
        V <- out3$V
      } else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  ##############################################################################
  # We need to recalculate the weights based on V_hat = V instead on solve(t(X)%*%X)
  # Do we have to? The differences look very small.
  ## compute mixing weights
  if (!(attr(object$wt, "method") == "none")) {
    if (nrow(Amat) == meq) {
    # equality constraints only
      wt <- rep(0L, ncol(V) + 1)
      wt.idx <- ncol(V) - meq + 1
      wt[wt.idx] <- 1
    } else if (attr(object$wt, "method") == "boot") { 
      # compute chi-square-bar weights based on Monte Carlo simulation
      wt <- con_weights_boot(VCOV     = V,
                             Amat     = Amat, 
                             meq      = meq, 
                             R        = attr(object$wt, "bootWt.R"),
                             parallel = parallel,
                             ncpus    = ncpus,
                             cl       = cl,
                             seed     = seed,
                             verbose  = verbose)
      attr(wt, "bootWt.R") <- attr(object$wt, "bootWt.R") 
    } else if (attr(object$wt, "method") == "pmvnorm" && (meq < nrow(Amat))) {
      # compute chi-square-bar weights based on pmvnorm
      wt <- rev(con_weights(Amat %*% V %*% t(Amat), meq = meq))
    } 
  } 
  attr(wt, "method") <- attr(object$wt, "method")
  ##############################################################################
  
  if (!(attr(object$wt, "method") == "none") && boot == "no") {
    # compute pvalue based on F-distribution
    pvalue <- con_pvalue_Fbar(wt          = rev(wt), 
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
    attr(pvalue, "wt") <- wt
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
                                         test     = "Wald",
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
              V           = V,
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
                                p.distr = rnorm, parallel = "no", ncpus = 1L,
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
  # unrestrikted scale estimate for the standard deviation: 
  tau <- sqrt(object$s2_unrestr)
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
  # tolerance
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
  }
  
  call_org$formula <- call_org$data <- call_org$weights <- 
    call_org$subset <- call_org$na.action <- call_org$model <- 
    call_org$x.ret <- call_org$y.ret <- call_org$contrasts <- 
    call_org$X <- call_org$y <- NULL
  
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
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
        CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
  
  if (!(attr(object$wt, "method") == "none") && boot == "no") {
    wt <- object$wt
    pvalue <- con_pvalue_Fbar(wt          = rev(wt), 
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
    attr(pvalue, "wt") <- wt
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
                                         test     = "Wald2",
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
                                p.distr = rnorm, parallel = "no", ncpus = 1L,
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
  # tolerance
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
  }
  
  call_org$formula <- call_org$data <- call_org$weights <- 
    call_org$subset <- call_org$na.action <- call_org$model <- 
    call_org$x.ret <- call_org$y.ret <- call_org$contrasts <- 
    call_org$X <- call_org$y <- NULL
  
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
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
    V <- out0$V
  } else if (type == "A") {
    CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
    V <- out1$V
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
      V <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq_alt > 0L && meq_alt <= meq) {
        CALL <- c(call_org, list(x = X, y = y, weights = weights,
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
        V <- out3$V
      } else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  # We need to recalculate the weights based on V_hat = Sigam instead on solve(t(X)%*%X)
  # Do we have to? The differences are very small.
  ## compute mixing weights
  if (!(attr(object$wt, "method") == "none")) {
    if (nrow(Amat) == meq) {
      # equality constraints only
      wt <- rep(0L, ncol(V) + 1)
      wt.idx <- ncol(V) - meq + 1
      wt[wt.idx] <- 1
    } else if (attr(object$wt, "method") == "boot") { 
      # compute chi-square-bar weights based on Monte Carlo simulation
      wt <- con_weights_boot(VCOV     = V,
                             Amat     = Amat, 
                             meq      = meq, 
                             R        = attr(object$wt, "bootWt.R"),
                             parallel = parallel,
                             ncpus    = ncpus,
                             cl       = cl,
                             seed     = seed,
                             verbose  = verbose)
      attr(wt, "bootWt.R") <- attr(object$wt, "bootWt.R")
    } else if (attr(object$wt, "method") == "pmvnorm" && (meq < nrow(Amat))) {
      # compute chi-square-bar weights based on pmvnorm
      wt <- rev(con_weights(Amat %*% V %*% t(Amat), meq = meq))
    } 
  } 
  attr(wt, "method") <- attr(object$wt, "method")
  
  if (!attr(object$wt, "method") == "none" && boot == "no") {
    # compute pvalue based on F-distribution
    pvalue <- con_pvalue_Fbar(wt          = rev(wt), 
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
    attr(pvalue, "wt") <- wt
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
              b_restr.alt = b_restr.alt,
              Sigma       = Sigma,
              V           = V,
              R2_org      = object$R2_org,
              R2_reduced  = object$R2_reduced,
              boot        = boot,
              model_org   = model_org)
  
  OUT <- list(OUT)
  names(OUT) <- type
  
  class(OUT) <- "conTest"
  
  OUT
}
