conTestF.rlm <- function(object, type = "A", boot = "no", meq.alt = 0,
                         control = NULL, tol = sqrt(.Machine$double.eps), ...) {
  
  if (type != "global") {
    type <- toupper(type)
  }  
  
  if (!("conRLM" %in% class(object))) {
    stop("object must be of class conRLM().")
  }
  if(!(type %in% c("A","B","C","global"))) {
    stop("type must be \"A\", \"B\", \"C\" or \"global\"")
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
  Y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
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
    stop("test not applicable for object with equality constraints only.")
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
    # call.rlm <- as.list(model.org$call)
    # call.rlm <- call.rlm[-1]
    # call.rlm[["weights"]] <- w
    # call.rlm[["data"]] <- NULL
    # call.rlm[["x"]] <- NULL
    # call.rlm[["y"]] <- NULL
    # call.my <- list(x = X, y = Y, Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)      
    # CALL <- c(call.rlm, call.my)
    # rfit <- do.call("conRLM_fit", CALL)

    call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)
        CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    Ts <- robustFm(x = X, y = Y,  beta0 = b.eqconstr, betaA = b.constr, 
                   scale = scale, cc = 4.685061)
  } else if (type == "A") {
    #call.rlm <- as.list(model.org$call)
    #call.rlm <- call.rlm[-1]
    #call.rlm[["weights"]] <- w
    #call.rlm[["data"]] <- NULL
    #call.rlm[["x"]] <- NULL
    #call.rlm[["y"]] <- NULL
    #call.my <- list(x = X, y = Y, Amat = Amat, meq = nrow(Amat), bvec = bvec)      
    #CALL <- c(call.rlm, call.my)
    #rfit <- do.call("conRLM_fit", CALL)

    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec)
        CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    Ts <- robustFm(x = X, y = Y,  beta0 = b.eqconstr, betaA = b.constr, 
                   scale = scale, cc = 4.685061)
  } else if (type == "B") {
    if (meq.alt == 0L) {
      Ts <- robustFm(x = X, y = Y,  beta0 = b.constr, betaA = b.unconstr, 
                     scale = scale, cc = 4.685061)
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        # call.rlm <- as.list(model.org$call)
        # call.rlm <- call.rlm[-1]
        # call.rlm[["weights"]] <- w
        # call.rlm[["data"]] <- NULL
        # call.rlm[["x"]] <- NULL
        # call.rlm[["y"]] <- NULL
        # call.my <- list(x = X, y = Y, Amat = Amat[1:meq.alt, , drop = FALSE], 
        #                 meq = meq.alt, bvec = bvec[1:meq.alt])      
        # CALL <- c(call.rlm, call.my)
        # rfit <- do.call("conRLM_fit", CALL)

        call.my <- list(Amat = Amat[1:meq.alt, , drop = FALSE], 
                        meq = meq.alt, bvec = bvec[1:meq.alt])
        CALL <- c(list(model.org), call.my)
        rfit <- do.call("conRLM_fit", CALL)
        
        b.constr.alt <- rfit$coefficients
        b.constr.alt[abs(b.constr.alt) < tol] <- 0L
        names(b.constr.alt) <- vnames
        Ts <- robustFm(x = X, y = Y,  beta0 = b.constr, betaA = b.constr.alt, 
                       scale = scale, cc = 4.685061)
      } else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  } else if (type == "C") { # intersection-union test (Sasabuchi, 1980)
    if (meq == 0L) {
      Ts <- as.vector(min((Amat %*% b.unconstr - bvec) /
                            sqrt(diag(Amat %*% COV %*% t(Amat)))))
      pvalue <- 1-pt(Ts, df.residual)
      names(pvalue) <- "pt.value"
    } else {
      stop("test not applicable with equality constraints.")
    }
  }
  
  if (!(type == "C")) {
    if (boot == "no") {
      # compute mixing weights
      wt.bar <- con_wt(Amat %*% COV %*% t(Amat), meq = meq)
      pvalue <- con_pvalue_Fbar(wt = wt.bar, Ts.org = Ts, 
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
      wt.bar <- mix.boot(object, type = type, meq.alt = meq.alt, 
                         R = ifelse(is.null(control$B), 9999, control$B),
                         parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                         ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                         cl = ifelse(is.null(control$cl), NULL, control$cl),
                         seed = ifelse(is.null(control$seed), 1234, control$seed),
                         verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
      
      pvalue <- con_pvalue_Fbar(wt = wt.bar, Ts.org = Ts, 
                                df.residual = df.residual, type = type,
                                Amat = Amat, bvec = bvec, meq = meq, 
                                meq.alt = meq.alt)
    }
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
    stop("object must be of class conRLM().")
  }
  if(!(type %in% c("A","B","C","global"))) {
    stop("type must be \"A\", \"B\", \"C\" or \"global\"")
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
  Y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
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
    stop("test not applicable for object with equality constraints only.")
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
    # call.rlm <- as.list(model.org$call)
    # call.rlm <- call.rlm[-1]
    # call.rlm[["weights"]] <- weights(model.org)
    # call.rlm[["data"]] <- NULL
    # call.rlm[["x"]] <- NULL
    # call.rlm[["y"]] <- NULL
    # call.my <- list(x = X, y = Y, Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)      
    # CALL <- c(call.rlm, call.my)
    # rfit <- do.call("conRLM_fit", CALL)
    
    call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)
    CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    out0 <- robustWaldScores(x = X, y = Y,  beta0 = b.eqconstr, 
                             betaA = b.constr, scale = scale)#, Amat = Amatg, 
                             #bvec = bvecg, meq = meq)
    Ts <- out0$RWald
    COV <- out0$V
  } else if (type == "A") {
    # call.rlm <- as.list(model.org$call)
    # call.rlm <- call.rlm[-1]
    # call.rlm[["weights"]] <- weights(model.org)
    # call.rlm[["data"]] <- NULL
    # call.rlm[["x"]] <- NULL
    # call.rlm[["y"]] <- NULL
    # call.my <- list(x = X, y = Y, Amat = Amat, meq = nrow(Amat), bvec = bvec)      
    # CALL <- c(call.rlm, call.my)
    # rfit <- do.call("conRLM_fit", CALL)
    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec)
    CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    
    out1 <- robustWaldScores(x = X, y = Y,  beta0 = b.eqconstr, 
                             betaA = b.constr, scale = scale)#, Amat = Amat, 
                             #bvec = bvec, meq = meq)
    Ts <- out1$RWald
    COV <- out1$V
  }
  else if (type == "B") {
    if (meq.alt == 0L) {
      out2 <- robustWaldScores(x = X, y = Y,  beta0 = b.constr, 
                               betaA = b.unconstr, scale = scale)#, Amat = Amat, 
                               #bvec = bvec, meq = Amat)
      Ts <- out2$RWald
      COV <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        # call.rlm <- as.list(model.org$call)
        # call.rlm <- call.rlm[-1]
        # call.rlm[["weights"]] <- weights(model.org)
        # call.rlm[["data"]] <- NULL
        # call.rlm[["x"]] <- NULL
        # call.rlm[["y"]] <- NULL
        # call.my <- list(x = X, y = Y, Amat = Amat[1:meq.alt,,drop=FALSE], 
        #                 meq = meq.alt, bvec = bvec[1:meq.alt])      
        # CALL <- c(call.rlm, call.my)
        # rfit <- do.call("conRLM_fit", CALL)
        
        call.my <- list(Amat = Amat[1:meq.alt,,drop=FALSE], meq = meq.alt, 
                        bvec = bvec[1:meq.alt])
        CALL <- c(list(model.org), call.my)
        rfit <- do.call("conRLM_fit", CALL)
    
        b.constr.alt <- rfit$coefficients
        b.constr.alt[abs(b.constr.alt) < tol] <- 0L
        names(b.constr.alt) <- vnames
        out3 <- robustWaldScores(x = X, y = Y,  beta0 = b.constr, 
                                 betaA = b.constr.alt, scale = scale)#, 
                                 #Amat = Amat, meq = meq.alt,                    
                                 #bvec = bvec)
        Ts <- out3$RWald
        COV <- out3$V
      } else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  } else if (type == "C") { # intersection-union test (Sasabuchi, 1980)
    if (meq == 0L) {
      Ts <- as.vector(min((Amat %*% b.unconstr - bvec) /
                            sqrt(diag(Amat %*% COV %*% t(Amat)))))
      pvalue <- 1-pt(Ts, df.residual)
      names(pvalue) <- "pt.value"
    } else {
      stop("test not applicable with equality constraints.")
    }
  }
  
  if (!(type == "C")) {
    if (boot == "no") {
      # compute weights
      wt.bar <- con_wt(Amat %*% COV %*% t(Amat), meq = meq)
      # compute pvalue based on F-distribution
      pvalue <- con_pvalue_Fbar(wt = wt.bar, Ts.org = Ts, 
                                df.residual = df.residual, type = type,
                                Amat = Amat, bvec = bvec, meq = meq, 
                                meq.alt = meq.alt)
    } else if (boot == "parametric") {
      pvalue <- con_pvalue_boot_parametric(object, Ts.org = Ts, type = type, test = "wald",
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
      pvalue <- con_pvalue_boot_model_based(object, Ts.org = Ts, type = type, test = "wald",
                                            meq.alt = meq.alt,
                                            R = ifelse(is.null(control$B), 9999, control$B),
                                            parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                            ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                            cl = ifelse(is.null(control$cl), NULL, control$cl),
                                            seed = ifelse(is.null(control$seed), 1234, control$seed),
                                            verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
    } else if (boot == "mix.weights") {
      # compute weights based on simulation
      wt.bar <- mix.boot(object, type = type, meq.alt = meq.alt, 
                         R = ifelse(is.null(control$B), 9999, control$B),
                         parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                         ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                         cl = ifelse(is.null(control$cl), NULL, control$cl),
                         seed = ifelse(is.null(control$seed), 1234, control$seed),
                         verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
      
      pvalue <- con_pvalue_Fbar(wt = wt.bar, Ts.org = Ts, 
                                df.residual = df.residual, type = type,
                                Amat = Amat, bvec = bvec, meq = meq, 
                                meq.alt = meq.alt)
    }
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
    stop("object must be of class conRLM().")
  }
  if(!(type %in% c("A","B","C","global"))) {
    stop("type must be \"A\", \"B\", \"C\" or \"global\"")
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
  Y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
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
    stop("test not applicable for object with equality constraints only.")
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
    # call.rlm <- as.list(model.org$call)
    # call.rlm <- call.rlm[-1]
    # call.rlm[["weights"]] <- w
    # #    if (is.null(call.rlm[["formula"]])) {
    # call.rlm[["data"]] <- NULL
    # call.rlm[["x"]] <- NULL
    # call.rlm[["y"]] <- NULL
    # call.my <- list(x = X, y = Y, Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)      
    # CALL <- c(call.rlm, call.my)
    # rfit <- do.call("conRLM_fit", CALL)
    call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)
        CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    out0 <- robustWaldScores(x = X, y = Y,  beta0 = b.eqconstr, 
                             betaA = b.constr, scale = scale)#, Amat = Amatg, 
    #bvec = bvecg, meq = meq)
    Ts <- out0$Rscore
    COV <- out0$V
  } else if (type == "A") {
    # call.rlm <- as.list(model.org$call)
    # call.rlm <- call.rlm[-1]
    # call.rlm[["weights"]] <- w
    # call.rlm[["data"]] <- NULL
    # call.rlm[["x"]] <- NULL
    # call.rlm[["y"]] <- NULL
    # call.my <- list(x = X, y = Y, Amat = Amat, meq = nrow(Amat), bvec = bvec)      
    # CALL <- c(call.rlm, call.my)
    # rfit <- do.call("conRLM_fit", CALL)
    
    call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec)
            CALL <- c(list(model.org), call.my)
    rfit <- do.call("conRLM_fit", CALL)
    
    b.eqconstr <- rfit$coefficients
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    
    out1 <- robustWaldScores(x = X, y = Y,  beta0 = b.eqconstr, 
                             betaA = b.constr, scale = scale)#, Amat = Amat, 
    #bvec = bvec, meq = meq)
    Ts <- out1$Rscore
    COV <- out1$V
  } else if (type == "B") {
    if (meq.alt == 0L) {
      out2 <- robustWaldScores(x = X, y = Y,  beta0 = b.constr, 
                               betaA = b.unconstr, scale = scale)#, Amat = Amat, 
      #bvec = bvec, meq = meq)
      Ts <- out2$Rscore
      COV <- out2$V
    } else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt != 0L && meq.alt <= meq) {
        # call.rlm <- as.list(model.org$call)
        # call.rlm <- call.rlm[-1]
        # call.rlm[["weights"]] <- w
        # call.rlm[["data"]] <- NULL
        # call.rlm[["x"]] <- NULL
        # call.rlm[["y"]] <- NULL
        # call.my <- list(x = X, y = Y, Amat = Amat[1:meq.alt,,drop=FALSE], 
        #                 meq = meq.alt, bvec = bvec[1:meq.alt])      
        # CALL <- c(call.rlm, call.my)
        # rfit <- do.call("conRLM_fit", CALL)
        
        call.my <- list(Amat = Amat[1:meq.alt,,drop=FALSE], meq = meq.alt, 
                        bvec = bvec[1:meq.alt])
                CALL <- c(list(model.org), call.my)
        rfit <- do.call("conRLM_fit", CALL)
        
        b.constr.alt <- rfit$coefficients
        b.constr.alt[abs(b.constr.alt) < tol] <- 0L
        names(b.constr.alt) <- vnames
        out3 <- robustWaldScores(x = X, y = Y,  beta0 = b.constr, 
                                 betaA = b.constr.alt, scale = scale)#,
        #Amat = Amat, meq = meq.alt,                    #[1:meq.alt,,drop=FALSE]
        #bvec = bvec)
        Ts <- out3$Rscore
        COV <- out3$V
      } else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  } else if (type == "C") { # intersection-union test (Sasabuchi, 1980)
    if (meq == 0L) {
      Ts <- as.vector(min((Amat %*% b.unconstr - bvec) /
                            sqrt(diag(Amat %*% COV %*% t(Amat)))))
      pvalue <- 1-pt(Ts, df.residual)
      names(pvalue) <- "pt.value"
    } else {
      stop("test not applicable with equality constraints.")
    }
  }
  
  if (!(type == "C")) {
    if (boot == "no") {
      # compute weights
      wt.bar <- con_wt(Amat %*% COV %*% t(Amat), meq = meq)
      # compute pvalue based on F-distribution
      pvalue <- con_pvalue_Fbar(wt = wt.bar, Ts.org = Ts, 
                                df.residual = df.residual, type = type,
                                Amat = Amat, bvec = bvec, meq = meq, 
                                meq.alt = meq.alt)
    } else if (boot == "parametric") {
      pvalue <- con_pvalue_boot_parametric(object, Ts.org = Ts, type = type, test = "Score",
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
      pvalue <- con_pvalue_boot_model_based(object, Ts.org = Ts, type = type, test = "Score",
                                            meq.alt = meq.alt,
                                            R = ifelse(is.null(control$B), 9999, control$B),
                                            parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                            ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                            cl = ifelse(is.null(control$cl), NULL, control$cl),
                                            seed = ifelse(is.null(control$seed), 1234, control$seed),
                                            verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
    } else if (boot == "mix.weights") {
      # compute weights based on simulation
      wt.bar <- mix.boot(object, type = type, meq.alt = meq.alt, 
                         R = ifelse(is.null(control$B), 9999, control$B),
                         parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                         ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                         cl = ifelse(is.null(control$cl), NULL, control$cl),
                         seed = ifelse(is.null(control$seed), 1234, control$seed),
                         verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
      
      pvalue <- con_pvalue_Fbar(wt = wt.bar, Ts.org = Ts, 
                                df.residual = df.residual, type = type,
                                Amat = Amat, bvec = bvec, meq = meq, 
                                meq.alt = meq.alt)
    }
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



# #<FIXME> remove intercept X[,-1] for assymetric error distributions.</FIXME>
# conTestWald2.rlm <- function(object, type = "A", boot = "no", meq.alt = 0,
#                              control = NULL, tol = sqrt(.Machine$double.eps), ...) {
#   
#   if (type != "global") {
#     type <- toupper(type)
#   }  
#   
#   if (!("conRLM" %in% class(object))) {
#     stop("object must be of class conRLM().")
#   }
#   if(!(type %in% c("A","B","C","global"))) {
#     stop("type must be \"A\", \"B\", \"C\" or \"global\"")
#   }
#   
#   if(!(boot %in% c("no", "residual", "model.based", "parametric", "mix.weights"))) {
#     stop("ERROR: boot method unknown.")
#   }
#   
#   if (boot == "residual") {
#     boot <- "model.based"
#   }
#   
#   
#   # original model
#   model.org <- object$model.org
#   X <- model.matrix(model.org)[,,drop=FALSE]
#   Y <- model.org$model[, attr(model.org$terms, "response")]
#   # unconstrained vcov
#   COV <- object$Sigma                                                           #use the robust vcovMM?
#   # unconstrained scale estimate for the standard errors
#   tau <- MASS:::summary.rlm(model.org)$stddev  
#   # unconstrained parameters
#   b.unconstr <- object$b.unconstr
#   # constrained parameters
#   b.constr <- object$b.constr
#   vnames <- names(b.unconstr)
#   # constraint syntax
#   Amat <- object$Amat
#   bvec <- object$bvec
#   meq  <- object$meq
#   # depends on type
#   b.eqconstr <- NULL
#   b.constr.alt <- NULL
#   
#   
#   if (meq == nrow(Amat)) {
#     stop("test not applicable for object with equality constraints only.")
#   }
#   if (!(length(b.unconstr) == nrow(COV))) {
#     stop("length b.unconstr and nrow(COV) must be identical.")
#   }
#   
#   if (type == "global") {
#     # check for intercept
#     intercept <- any(attr(terms(model.org), "intercept"))
#     g <- length(object$b.constr)
#     if (intercept) {
#       Amatg <- cbind(rep(0, (g - 1)), diag(rep(1, g - 1))) 
#       bvecg <- rep(0, g - 1) 
#     } else {
#       stop("Restriktor ERROR: test not applicable for models without intercept.")      
#     } 
#     #fit inequality constrained robust model
#     call.rlm <- as.list(model.org$call)
#     call.rlm <- call.rlm[-1]
#     call.rlm[["weights"]] <- weights(model.org)
#     #    if (is.null(call.rlm[["formula"]])) {
#     call.rlm[["data"]] <- NULL
#     call.rlm[["x"]] <- NULL
#     call.rlm[["y"]] <- NULL
#     call.my <- list(x = X, y = Y, Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)      
#     CALL <- c(call.rlm, call.my)
#     rfit <- do.call("conRLM_fit", CALL)
#     #    }
#     #    else {
#     #      call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg)  
#     #      call.rlm[["data"]] <- as.name("Data")
#     #      CALL <- c(call.rlm, call.my)
#     #      rfit <- do.call("conRLM.formula", CALL)
#     #    }
#     b.eqconstr <- rfit$coefficients
#     b.eqconstr[abs(b.eqconstr) < tol] <- 0L
#     names(b.eqconstr) <- vnames
#     
#     Ts <- robustWaldXX(x = X, beta0 = b.eqconstr, beta1 = b.constr, 
#                        beta2 = b.unconstr, tau = tau) 
#   } else if (type == "A") {
#     call.rlm <- as.list(model.org$call)
#     call.rlm <- call.rlm[-1]
#     call.rlm[["weights"]] <- weights(model.org)
#     #   if (is.null(call.rlm[["formula"]])) {
#     call.rlm[["data"]] <- NULL
#     call.rlm[["x"]] <- NULL
#     call.rlm[["y"]] <- NULL
#     call.my <- list(x = X, y = Y, Amat = Amat, meq = nrow(Amat), bvec = bvec)      
#     CALL <- c(call.rlm, call.my)
#     rfit <- do.call("conRLM_fit", CALL)
#     #    }
#     #    else {
#     #      call.my <- list(Amat = Amat, meq = nrow(Amat), bvec = bvec)  
#     #      call.rlm[["data"]] <- as.name("Data")
#     #      CALL <- c(call.rlm, call.my)
#     #      rfit <- do.call("conRLM.formula", CALL)
#     #    }
#     b.eqconstr <- rfit$coefficients
#     b.eqconstr[abs(b.eqconstr) < tol] <- 0L
#     names(b.eqconstr) <- vnames
#     
#     Ts <- robustWaldXX(x = X, beta0 = b.eqconstr, beta1 = b.constr, 
#                        beta2 = b.unconstr, tau = tau)   
#   }
#   else if (type == "B") {
#     if (meq.alt == 0L) {
#       Ts <- robustWaldXX(x = X, beta0 = b.constr, beta1 = b.unconstr, 
#                          beta2 = b.unconstr, tau = tau)   
#     }
#     else {
#       # some equality may be preserved in the alternative hypothesis.
#       if(meq.alt != 0L && meq.alt <= meq) {
#         call.rlm <- as.list(model.org$call)
#         call.rlm <- call.rlm[-1]
#         call.rlm[["weights"]] <- weights(model.org)
#         #        if (is.null(call.rlm[["formula"]])) {
#         call.rlm[["data"]] <- NULL
#         call.rlm[["x"]] <- NULL
#         call.rlm[["y"]] <- NULL
#         call.my <- list(x = X, y = Y, Amat = Amat[1:meq.alt,,drop=FALSE], 
#                         meq = meq.alt, bvec = bvec[1:meq.alt])      
#         CALL <- c(call.rlm, call.my)
#         rfit <- do.call("conRLM_fit", CALL)
#         #       }
#         #        else {
#         #          call.my <- list(Amat = Amat[1:meq.alt,,drop=FALSE], meq = meq.alt, 
#         #                          bvec = bvec[1:meq.alt])  
#         #          call.rlm[["data"]] <- as.name("Data")
#         #          CALL <- c(call.rlm, call.my)
#         #          rfit <- do.call("conRLM.formula", CALL)
#         #        }
#         b.constr.alt <- rfit$coefficients
#         b.constr.alt[abs(b.constr.alt) < tol] <- 0L
#         names(b.constr.alt) <- vnames
#         Ts <- robustWaldXX(x = X, beta0 = b.constr, beta1 = b.constr.alt, 
#                            beta2 = b.unconstr, tau = tau)   
#       }
#       else {
#         stop("meq.alt must not be larger than meq.")
#       }
#     }
#   } else if (type == "C") { # intersection-union test (Sasabuchi, 1980)
#     if (meq == 0L) {
#       Ts <- as.vector(min((Amat %*% b.unconstr - bvec) /
#                             sqrt(diag(Amat %*%COV%*% t(Amat)))))
#       #      names(Ts) <- "Tbar"
#     } else {
#       stop("test not applicable with equality constraints.")
#     }
#   }
#   
#   if (boot == "no") {
#     pvalue <- con_pvalue_Fbar(COV = crossprod(X), Ts.org = Ts, df.residual, type = type,
#                               Amat, bvec, meq, meq.alt)
#   } else if (boot == "parametric") {
#     pvalue <- con_pvalue_boot_parametric(object, Ts.org = Ts, type = type, test = "wald2",
#                                          meq.alt = meq.alt,
#                                          R = ifelse(is.null(control$B), 9999, control$B),
#                                          p.distr = ifelse(is.null(control$p.distr), "N", control$p.distr),
#                                          df = ifelse(is.null(control$df), 7, control$df),
#                                          parallel = ifelse(is.null(control$parallel), "no", control$parallel),
#                                          ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
#                                          cl = ifelse(is.null(control$cl), NULL, control$cl),
#                                          seed = ifelse(is.null(control$seed), 1234, control$seed),
#                                          verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
#   } else if (boot == "model.based") {
#     pvalue <- con_pvalue_boot_model_based(object, Ts.org = Ts, type = type, test = "wald2",
#                                           meq.alt = meq.alt,
#                                           R = ifelse(is.null(control$B), 9999, control$B),
#                                           parallel = ifelse(is.null(control$parallel), "no", control$parallel),
#                                           ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
#                                           cl = ifelse(is.null(control$cl), NULL, control$cl),
#                                           seed = ifelse(is.null(control$seed), 1234, control$seed),
#                                           verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
#   } else if (boot == "mix.weights") {
#     pvalue <- con_pvalue_boot_weights(object, pbar = "pfbar", Ts.org = Ts, df.residual = df.residual, 
#                                       type = type, meq.alt = meq.alt, 
#                                       R = ifelse(is.null(control$B), 9999, control$B),
#                                       parallel = ifelse(is.null(control$parallel), "no", control$parallel),
#                                       ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
#                                       cl = ifelse(is.null(control$cl), NULL, control$cl),
#                                       seed = ifelse(is.null(control$seed), 1234, control$seed),
#                                       verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
#   }
#   
#   OUT <- list(CON = object$CON,
#               type = type,
#               boot = boot,
#               b.eqconstr = NULL,
#               b.unconstr = b.unconstr,
#               b.constr = b.constr,
#               b.constr.alt = b.constr.alt,
#               Amat = Amat,
#               bvec = bvec,
#               meq = meq,
#               meq.alt = meq.alt,
#               iact = object$iact,
#               df.residual = df.residual,
#               COV = COV,
#               Ts = Ts,
#               pvalue = pvalue,
#               model.org = model.org)
#   
#   
#   if (type == "A" | type == "global") { 
#     OUT$b.eqconstr <- b.eqconstr 
#   }
#   
#   class(OUT) <- "conTest"
#   
#   OUT
#   
# }