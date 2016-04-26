conTestF.rlm <- function(object, type = "A", boot = "no", meq.alt = 0,
                         control = NULL, tol = sqrt(.Machine$double.eps), ...) {
  
  if (type != "global") {
    type <- toupper(type)
  }  
  if (!("conRLM" %in% class(object))) {
    stop("object must be of class conRLM.")
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
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # weights
  #w <- weights(model.org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  COV <- object$Sigma
  # unconstrained scale
  scale <- model.org$s
  # parameter estimates
  b.unconstr <- object$b.unconstr
  b.constr <- object$b.constr
  b.eqconstr <- NULL
  b.constr.alt <- NULL
  # variable names
  vnames <- names(b.unconstr)
  # constraints stuff
  Amat <- object$Amat
  bvec <- object$bvec
  meq  <- object$meq
  
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality restriktions only.")
  }
  
  if (type == "global") {
    # check for intercept
    intercept <- any(attr(terms(model.org), "intercept"))
    g <- length(object$b.constr)
    if (intercept) {
      Amatg <- cbind(rep(0, (g - 1)), diag(rep(1, g - 1))) 
      bvecg <- rep(0, g - 1) 
    } else {
      stop("Restriktor ERROR: test not applicable for models without intercept.")      
    } 
    #fit inequality constrained robust model
    call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)
        CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    Ts <- robustFm(x = X, y = y,  beta0 = b.eqconstr, betaA = b.constr, 
                   scale = scale, cc = 4.685061)
  } else if (type == "A") {
    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec)
        CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    Ts <- robustFm(x = X, y = y,  beta0 = b.eqconstr, betaA = b.constr, 
                   scale = scale, cc = 4.685061)
  } else if (type == "B") {
    if (meq.alt == 0L) {
      Ts <- robustFm(x = X, y = y,  beta0 = b.constr, betaA = b.unconstr, 
                     scale = scale, cc = 4.685061)
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        call.my <- list(Amat = Amat[1:meq.alt, , drop = FALSE], 
                        meq = meq.alt, bvec = bvec[1:meq.alt])
        CALL <- c(list(model.org), call.my)
        rfit <- do.call("conRLM_fit", CALL)
        
        b.constr.alt <- rfit$coefficients
        b.constr.alt[abs(b.constr.alt) < tol] <- 0L
        names(b.constr.alt) <- vnames
        Ts <- robustFm(x = X, y = y,  beta0 = b.constr, betaA = b.constr.alt, 
                       scale = scale, cc = 4.685061)
      } else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  } 
  
  if (boot == "no") {
    # compute mixing weights
    wt <- con_wt(Amat %*% COV %*% t(Amat), meq = meq)
    pvalue <- con_pvalue_Fbar(wt = wt, Ts.org = Ts, 
                              df.residual = df.residual, type = type,
                              Amat = Amat, bvec = bvec, meq = meq, 
                              meq.alt = meq.alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, Ts.org = Ts, type = type, test = "F",
                                         meq.alt = meq.alt,
                                         R = ifelse(is.null(control$B), 9999, control$B),
                                         p.distr = ifelse(is.null(control$p.distr), "N", control$p.distr),
                                         df = ifelse(is.null(control$df), 7, control$df),
                                         parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                         ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                         cl = ifelse(is.null(control$cl), NULL, control$cl),
                                         seed = ifelse(is.null(control$seed), 1234, control$seed),
                                         verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  } else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based(object, Ts.org = Ts, type = type, test = "F",
                                          meq.alt = meq.alt,
                                          R = ifelse(is.null(control$B), 9999, control$B),
                                          parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                          ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                          cl = ifelse(is.null(control$cl), NULL, control$cl),
                                          seed = ifelse(is.null(control$seed), 1234, control$seed),
                                          verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  } else if (boot == "mix.weights") {
    # compute weights based on simulation
    wt <- mix.boot(object, type = type, meq.alt = meq.alt, 
                       R = ifelse(is.null(control$B), 9999, control$B),
                       parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                       ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                       cl = ifelse(is.null(control$cl), NULL, control$cl),
                       seed = ifelse(is.null(control$seed), 1234, control$seed),
                       verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
    
    wt <- rev(wt)
    wt.idx <- which(wt == 0)
    wt <- wt[-wt.idx]
    
    pvalue <- con_pvalue_Fbar(wt = wt, Ts.org = Ts, 
                              df.residual = df.residual, type = type,
                              Amat = Amat, bvec = bvec, meq = meq, 
                              meq.alt = meq.alt)
  }
  
  OUT <- list(CON = object$CON,
              type = type,
              boot = boot,
              b.eqconstr = NULL,
              b.unconstr = b.unconstr,
              b.constr = b.constr,
              b.constr.alt = b.constr.alt,
              Amat = Amat,
              bvec = bvec,
              meq = meq,
              meq.alt = meq.alt,
              iact = object$iact,
              df.residual = df.residual,
              COV = COV,
              Ts = Ts,
              pvalue = pvalue,
              model.org = model.org)
  
  if (type == "A" | type == "global") { 
    OUT$b.eqconstr <- b.eqconstr 
  }
  
  class(OUT) <- "conTest"
  
  OUT
  
}





conTestWald.rlm <- function(object, type = "A", boot = "no", meq.alt = 0,
                             control = NULL, tol = sqrt(.Machine$double.eps), ...) {
  
  if (type != "global") {
    type <- toupper(type)
  }  
  
  if (!("conRLM" %in% class(object))) {
    stop("object must be of class conRLM.")
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
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # weights
  #w <- weights(model.org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  COV <- object$Sigma
  # unconstrained scale
  scale <- model.org$s
  # parameter estimates
  b.unconstr <- object$b.unconstr
  b.constr <- object$b.constr
  b.eqconstr <- NULL
  b.constr.alt <- NULL
  # variable names
  vnames <- names(b.unconstr)
  # constraints stuff
  Amat <- object$Amat
  bvec <- object$bvec
  meq  <- object$meq
  
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality restriktions only.")
  }
  
  if (type == "global") {
    # check for intercept
    intercept <- any(attr(terms(model.org), "intercept"))
    g <- length(object$b.constr)
    if (intercept) {
      Amatg <- cbind(rep(0, (g - 1)), diag(rep(1, g - 1))) 
      bvecg <- rep(0, g - 1) 
    } else {
      stop("Restriktor ERROR: test not applicable for models without intercept.")      
    } 
    #fit inequality constrained robust model
    call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)
    CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    out0 <- robustWaldScores(x = X, y = y,  beta0 = b.eqconstr, 
                             betaA = b.constr, scale = scale, test = "Wald")
    Ts <- out0$Ts
    COV <- out0$V
  } else if (type == "A") {
    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec)
    CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    
    out1 <- robustWaldScores(x = X, y = y,  beta0 = b.eqconstr, 
                             betaA = b.constr, scale = scale, test = "Wald")
    Ts <- out1$Ts
    COV <- out1$V
  }
  else if (type == "B") {
    if (meq.alt == 0L) {
      out2 <- robustWaldScores(x = X, y = y,  beta0 = b.constr, 
                               betaA = b.unconstr, scale = scale, test = "Wald")
      Ts <- out2$Ts
      COV <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        call.my <- list(Amat = Amat[1:meq.alt,,drop=FALSE], meq = meq.alt, 
                        bvec = bvec[1:meq.alt])
        CALL <- c(list(model.org), call.my)
        rfit <- do.call("conRLM_fit", CALL)
    
        b.constr.alt <- rfit$coefficients
        b.constr.alt[abs(b.constr.alt) < tol] <- 0L
        names(b.constr.alt) <- vnames
        out3 <- robustWaldScores(x = X, y = y,  beta0 = b.constr, 
                                 betaA = b.constr.alt, scale = scale, 
                                 test = "Wald")
        Ts <- out3$Ts
        COV <- out3$V
      } else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  } 
  
  if (boot == "no") {
    # compute weights
    wt <- con_wt(Amat %*% COV %*% t(Amat), meq = meq)
    # compute pvalue based on F-distribution
    pvalue <- con_pvalue_Fbar(wt = wt, Ts.org = Ts, 
                              df.residual = df.residual, type = type,
                              Amat = Amat, bvec = bvec, meq = meq, 
                              meq.alt = meq.alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, Ts.org = Ts, type = type, test = "Wald",
                                         meq.alt = meq.alt,
                                         R = ifelse(is.null(control$B), 9999, control$B),
                                         p.distr = ifelse(is.null(control$p.distr), "N", control$p.distr),
                                         df = ifelse(is.null(control$df), 7, control$df),
                                         parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                         ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                         cl = ifelse(is.null(control$cl), NULL, control$cl),
                                         seed = ifelse(is.null(control$seed), 1234, control$seed),
                                         verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  } else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based(object, Ts.org = Ts, type = type, test = "Wald",
                                          meq.alt = meq.alt,
                                          R = ifelse(is.null(control$B), 9999, control$B),
                                          parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                          ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                          cl = ifelse(is.null(control$cl), NULL, control$cl),
                                          seed = ifelse(is.null(control$seed), 1234, control$seed),
                                          verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  } else if (boot == "mix.weights") {
    # compute weights based on simulation
    wt <- mix.boot(object, type = type, meq.alt = meq.alt, 
                       R = ifelse(is.null(control$B), 9999, control$B),
                       parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                       ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                       cl = ifelse(is.null(control$cl), NULL, control$cl),
                       seed = ifelse(is.null(control$seed), 1234, control$seed),
                       verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
    
    wt <- rev(wt)
    wt.idx <- which(wt == 0)
    wt <- wt[-wt.idx]
    
    pvalue <- con_pvalue_Fbar(wt = wt, Ts.org = Ts, 
                              df.residual = df.residual, type = type,
                              Amat = Amat, bvec = bvec, meq = meq, 
                              meq.alt = meq.alt)
  }
  
  OUT <- list(CON = object$CON,
              type = type,
              boot = boot,
              b.eqconstr = NULL,
              b.unconstr = b.unconstr,
              b.constr = b.constr,
              b.constr.alt = b.constr.alt,
              Amat = Amat,
              bvec = bvec,
              meq = meq,
              meq.alt = meq.alt,
              iact = object$iact,
              df.residual = df.residual,
              COV = COV,
              Ts = Ts,
              pvalue = pvalue,
              model.org = model.org)
  
  
  if (type == "A" | type == "global") { 
    OUT$b.eqconstr <- b.eqconstr 
  }
  
  class(OUT) <- "conTest"
  
  OUT
  
}




conTestScore.rlm <- function(object, type = "A", boot = "no", meq.alt = 0,
                             control = NULL, tol = sqrt(.Machine$double.eps), ...) {
  
  if (type != "global") {
    type <- toupper(type)
  }  
  if (!("conRLM" %in% class(object))) {
    stop("object must be of class conRLM.")
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
  X <- model.matrix(model.org)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # weights
  #w <- weights(model.org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  COV <- object$Sigma
  # unconstrained scale
  scale <- model.org$s
  # parameter estimates
  b.unconstr <- object$b.unconstr
  b.constr <- object$b.constr
  b.eqconstr <- NULL
  b.constr.alt <- NULL
  # variable names
  vnames <- names(b.unconstr)
  # constraints stuff
  Amat <- object$Amat
  bvec <- object$bvec
  meq  <- object$meq
  
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality restriktions only.")
  }
  
  if (type == "global") {
    # check for intercept
    intercept <- any(attr(terms(model.org), "intercept"))
    g <- length(object$b.constr)
    if (intercept) {
      Amatg <- cbind(rep(0, (g - 1)), diag(rep(1, g - 1))) 
      bvecg <- rep(0, g - 1) 
    } else {
      stop("Restriktor ERROR: test not applicable for models without intercept.")      
    } 
    #fit inequality constrained robust model
    call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)
        CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    out0 <- robustWaldScores(x = X, y = y,  beta0 = b.eqconstr, 
                             betaA = b.constr, scale = scale, test = "score")    
    Ts <- out0$Ts
    COV <- out0$V
  } else if (type == "A") {
    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec)
            CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    
    out1 <- robustWaldScores(x = X, y = y,  beta0 = b.eqconstr, 
                             betaA = b.constr, scale = scale, test = "score")
    Ts <- out1$Ts
    COV <- out1$V
  } else if (type == "B") {
    if (meq.alt == 0L) {
      out2 <- robustWaldScores(x = X, y = y,  beta0 = b.constr, 
                               betaA = b.unconstr, scale = scale, test = "score")
      Ts <- out2$Ts
      COV <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        call.my <- list(Amat = Amat[1:meq.alt,,drop=FALSE], meq = meq.alt, 
                        bvec = bvec[1:meq.alt])
                CALL <- c(list(model.org), call.my)
        rfit <- do.call("conRLM_fit", CALL)
        
        b.constr.alt <- rfit$coefficients
        b.constr.alt[abs(b.constr.alt) < tol] <- 0L
        names(b.constr.alt) <- vnames
        out3 <- robustWaldScores(x = X, y = y,  beta0 = b.constr, 
                                 betaA = b.constr.alt, scale = scale, 
                                 test = "score")
        Ts <- out3$Ts
        COV <- out3$V
      } else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  } 
  
  if (boot == "no") {
    # compute weights
    wt <- con_wt(Amat %*% COV %*% t(Amat), meq = meq)
    # compute pvalue based on F-distribution
    pvalue <- con_pvalue_Fbar(wt = wt, Ts.org = Ts, 
                              df.residual = df.residual, type = type,
                              Amat = Amat, bvec = bvec, meq = meq, 
                              meq.alt = meq.alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, Ts.org = Ts, type = type, test = "score",
                                         meq.alt = meq.alt,
                                         R = ifelse(is.null(control$B), 9999, control$B),
                                         p.distr = ifelse(is.null(control$p.distr), "N", control$p.distr),
                                         df = ifelse(is.null(control$df), 7, control$df),
                                         parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                         ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                         cl = ifelse(is.null(control$cl), NULL, control$cl),
                                         seed = ifelse(is.null(control$seed), 1234, control$seed),
                                         verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  } else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based(object, Ts.org = Ts, type = type, test = "score",
                                          meq.alt = meq.alt,
                                          R = ifelse(is.null(control$B), 9999, control$B),
                                          parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                          ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                          cl = ifelse(is.null(control$cl), NULL, control$cl),
                                          seed = ifelse(is.null(control$seed), 1234, control$seed),
                                          verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  } else if (boot == "mix.weights") {
    # compute weights based on simulation
    wt <- mix.boot(object, type = type, meq.alt = meq.alt, 
                       R = ifelse(is.null(control$B), 9999, control$B),
                       parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                       ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                       cl = ifelse(is.null(control$cl), NULL, control$cl),
                       seed = ifelse(is.null(control$seed), 1234, control$seed),
                       verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
    
    wt <- rev(wt)
    wt.idx <- which(wt == 0)
    wt <- wt[-wt.idx]
    
    pvalue <- con_pvalue_Fbar(wt = wt, Ts.org = Ts, 
                              df.residual = df.residual, type = type,
                              Amat = Amat, bvec = bvec, meq = meq, 
                              meq.alt = meq.alt)
  }
  
  OUT <- list(CON = object$CON,
              type = type,
              boot = boot,
              b.eqconstr = NULL,
              b.unconstr = b.unconstr,
              b.constr = b.constr,
              b.constr.alt = b.constr.alt,
              Amat = Amat,
              bvec = bvec,
              meq = meq,
              meq.alt = meq.alt,
              iact = object$iact,
              df.residual = df.residual,
              COV = COV,
              Ts = Ts,
              pvalue = pvalue,
              model.org = model.org)
  
  
  if (type == "A" | type == "global") { 
    OUT$b.eqconstr <- b.eqconstr 
  }
  
  class(OUT) <- "conTest"
  
  OUT
  
}
