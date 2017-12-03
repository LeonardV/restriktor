### computes the F-bar, Wald(2)-bar and score-bar test statistic ####
conTestF.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                            p.distr = rnorm, parallel = "no", ncpus = 1L,
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
  # orignal model call
  call.org <- as.list(model.org$call)
  # model matrix
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # weights
  weights <- model.org$weights
  # unconstrained df
  df.residual <- object$df.residual
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
  control <- c(object$control, control)
  # remove duplicated elements from control list
  control <- control[!duplicated(control)]
  # get tolerance for control if exists
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
  }
  
  call.org$formula <- call.org$data <- call.org$weights <- 
    call.org$subset <- call.org$na.action <- call.org$model <- 
    call.org$x.ret <- call.org$y.ret <- call.org$contrasts <- 
    call.org$X <- call.org$y <- NULL
    
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
    CALL <- c(call.org, list(x = X, y = y, weights = weights,
                             Amat = AmatG, bvec = bvecG, 
                             meq = nrow(AmatG), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    Ts <- robustFm(x         = X, 
                   y         = y, 
                   b.unrestr = b.unrestr,
                   b.eqrestr = b.eqrestr, 
                   b.restr   = b.restr, 
                   scale     = scale, 
                   cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))$Ts
  } else if (type == "A") {
    CALL <- c(call.org, list(x = X, y = y, weights = weights,
                             Amat = Amat, bvec = bvec, 
                             meq = nrow(Amat), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    Ts <- robustFm(x         = X, 
                   y         = y, 
                   b.unrestr = b.unrestr,
                   b.eqrestr = b.eqrestr, 
                   b.restr   = b.restr, 
                   scale     = scale, 
                   cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))$Ts
  } else if (type == "B") {
    if (meq.alt == 0L) {
      Ts <- robustFm(x         = X, 
                     y         = y, 
                     b.unrestr = b.unrestr,
                     b.eqrestr = b.restr, 
                     b.restr   = b.unrestr, 
                     scale     = scale, 
                     cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))$Ts
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        CALL <- c(call.org, list(x = X, y = y, weights = weights,
                                 Amat = Amat[1:meq.alt, , drop = FALSE], 
                                 meq = meq.alt, bvec = bvec[1:meq.alt],
                                 tol = tol))
        
        rfit <- do.call("conRLM_fit", CALL)
        b.restr.alt <- rfit$coefficients
        b.restr.alt[abs(b.restr.alt) < tol] <- 0L
        names(b.restr.alt) <- vnames
        Ts <- robustFm(x         = X, 
                       y         = y, 
                       b.unrestr = b.unrestr,
                       b.eqrestr = b.restr, 
                       b.restr   = b.restr.alt, 
                       scale     = scale, 
                       cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))$Ts
      } else {
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
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              boot        = boot,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}



#REF: Silvapulle, M.J. (1992) Robust Wald-Type Test of One-Sided Hypotheses in the linear model.
#American Statistical Association, 87, 417.
conTestWald.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                               p.distr = rnorm, parallel = "no", ncpus = 1L,
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
  # origianl model call
  call.org <- as.list(model.org$call)
  # model matrix
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # weights
  weights <- model.org$weights
  # unconstrained df
  df.residual <- object$df.residual
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
  control <- c(object$control, control)
  # remove duplicated elements from control list
  control <- control[!duplicated(control)]
  # get tolerance for control if exists
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
  }
  
  call.org$formula <- call.org$data <- call.org$weights <- 
    call.org$subset <- call.org$na.action <- call.org$model <- 
    call.org$x.ret <- call.org$y.ret <- call.org$contrasts <- 
    call.org$X <- call.org$y <- NULL
  
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
    CALL <- c(call.org, list(x = X, y = y, weights = weights,
                             Amat = AmatG, bvec = bvecG, 
                             meq = nrow(AmatG), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    Ts <- robustWald(x         = X,
                     y         = y,
                     b.eqrestr = b.eqrestr, 
                     b.restr   = b.restr, 
                     b.unrestr = b.unrestr,
                     Amat      = AmatG,
                     scale     = scale,
                     cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))$Ts
    
  } else if (type == "A") {
    CALL <- c(call.org, list(x = X, y = y, weights = weights,
                             Amat = Amat, bvec = bvec, 
                             meq = nrow(Amat), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    Ts <- robustWald(x         = X,
                     y         = y,
                     b.eqrestr = b.eqrestr, 
                     b.restr   = b.restr, 
                     b.unrestr = b.unrestr, 
                     Amat      = Amat,
                     scale     = scale,
                     cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))$Ts
  } else if (type == "B") {
    if (meq.alt == 0L) {
      Ts <- robustWald(x         = X,
                       y         = y,
                       b.eqrestr = b.restr, 
                       b.restr   = b.unrestr, 
                       b.unrestr = b.unrestr, 
                       Amat      = Amat,
                       scale     = scale,
                       cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))$Ts
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        CALL <- c(call.org, list(x = X, y = y, weights = weights,
                                 Amat = Amat[1:meq.alt, , drop = FALSE], 
                                 meq = meq.alt, bvec = bvec[1:meq.alt], 
                                 tol = tol))
        
        rfit <- do.call("conRLM_fit", CALL)
        b.restr.alt <- rfit$coefficients
        b.restr.alt[abs(b.restr.alt) < tol] <- 0L
        names(b.restr.alt) <- vnames
        Ts <- robustWald(x         = X,
                         y         = y,
                         b.eqrestr = b.restr, 
                         b.restr   = b.restr.alt, 
                         b.unrestr = b.restr.alt,
                         Amat      = Amat,
                         scale     = scale,
                         cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))$Ts
      } else {
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
                                         test     = "Wald",
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
              test        = "Wald",
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



conTestScore.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                                p.distr = rnorm, parallel = "no", ncpus = 1L,
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
  # original model call
  call.org <- as.list(model.org$call)
  # model matrix
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # weights
  weights <- model.org$weights
  # unconstrained df
  df.residual <- object$df.residual
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
  control <- c(object$control, control)
  # remove duplicated elements from control list
  control <- control[!duplicated(control)]
  # get tolerance for control if exists
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restriktions only.")
  }
  
  call.org$formula <- call.org$data <- call.org$weights <- 
    call.org$subset <- call.org$na.action <- call.org$model <- 
    call.org$x.ret <- call.org$y.ret <- call.org$contrasts <- 
    call.org$X <- call.org$y <- NULL
  
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
    CALL <- c(call.org, list(x = X, y = y, weights = weights,
                             Amat = AmatG, bvec = bvecG, 
                             meq = nrow(AmatG), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    out0 <- robustScores(x         = X, 
                         y         = y,  
                         b.eqrestr = b.eqrestr, 
                         b.restr   = b.restr, 
                         b.unrestr = b.unrestr,
                         Amat      = AmatG,
                         scale     = scale, 
                         test      = "score", 
                         cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))
    Ts <- out0$Ts
    V <- out0$V
  } else if (type == "A") {
    CALL <- c(call.org, list(x = X, y = y, weights = weights,
                             Amat = Amat, bvec = bvec, 
                             meq = nrow(Amat), tol = tol))
    
    rfit <- do.call("conRLM_fit", CALL)
    b.eqrestr <- rfit$coefficients
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    names(b.eqrestr) <- vnames
    
    out1 <- robustScores(x         = X, 
                         y         = y,
                         b.eqrestr = b.eqrestr, 
                         b.restr   = b.restr, 
                         b.unrestr = b.unrestr,
                         Amat      = Amat,
                         scale     = scale, 
                         test      = "score", 
                         cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))
    Ts <- out1$Ts
    V <- out1$V
  } else if (type == "B") {
    if (meq.alt == 0L) {
      out2 <- robustScores(x         = X, 
                           y         = y,  
                           b.eqrestr = b.restr, 
                           b.restr   = b.unrestr,
                           b.unrestr = b.unrestr,
                           scale     = scale,
                           Amat      = Amat,
                           test      = "score", 
                           cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))
      Ts <- out2$Ts
      V <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        CALL <- c(call.org, list(x = X, y = y, weights = weights,
                                 Amat = Amat[1:meq.alt,,drop=FALSE], 
                                 meq = meq.alt, 
                                 bvec = bvec[1:meq.alt], tol = tol))
        
        rfit <- do.call("conRLM_fit", CALL)
        b.restr.alt <- rfit$coefficients
        b.restr.alt[abs(b.restr.alt) < tol] <- 0L
        names(b.restr.alt) <- vnames
        out3 <- robustScores(x         = X, 
                             y         = y,  
                             b.eqrestr = b.restr, 
                             b.restr   = b.restr.alt,
                             b.unrestr = b.restr.alt,
                             scale     = scale, 
                             Amat      = Amat,
                             test      = "score", 
                             cc        = ifelse(is.null(call.org$c), 4.685061, call.org$c))
        Ts <- out3$Ts
        V <- out3$V
      } else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  # We need to recalculate the weights based on V.hat = Sigma instead on solve(t(X)%*%X)
  # Is this really needed? The differences look very small.
  ## compute mixing weights
  if (!(attr(object$wt.bar, "method") == "none")) {
    if (nrow(Amat) == meq) {
      # equality constraints only
      wt.bar <- rep(0L, ncol(V) + 1)
      wt.bar.idx <- ncol(V) - meq + 1
      wt.bar[wt.bar.idx] <- 1
    } else if (attr(object$wt.bar, "method") == "boot") { 
      # compute chi-square-bar weights based on Monte Carlo simulation
      wt.bar <- con_weights_boot(VCOV     = V,
                                 Amat     = Amat, 
                                 meq      = meq, 
                                 R        = attr(object$wt.bar, "mix.bootstrap"),
                                 parallel = parallel,
                                 ncpus    = ncpus,
                                 cl       = cl,
                                 seed     = seed,
                                 verbose  = verbose)
      attr(wt.bar, "mix.bootstrap") <- attr(object$wt.bar, "mix.bootstrap")
    } else if (attr(object$wt.bar, "method") == "pmvnorm" && (meq < nrow(Amat))) {
      # compute chi-square-bar weights based on pmvnorm
      wt.bar <- rev(con_weights(Amat %*% V %*% t(Amat), meq = meq))
    } 
  } else {
    wt.bar <- as.numeric(NA)
  } 
  attr(wt.bar, "method") <- attr(object$wt.bar, "method")
  
  if (!attr(object$wt.bar, "method") == "none" && boot == "no") {
    # compute pvalue based on F-distribution
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
              V           = V,
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              boot        = boot,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}
