conTestF.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", B = 9999, 
                            p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                            cl = NULL, seed = 1234, verbose = FALSE,
                            control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conRLM" %in% class(object))) {
    stop("object must be of class conRLM.")
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
    #fit inequality constrained robust model
    call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg,
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
                   b.eqrestr = b.eqrestr, 
                   b.restr   = b.restr, 
                   scale     = scale, 
                   cc        = ifelse(is.null(cc), 4.685061, cc))
  } else if (type == "A") {
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
                   b.eqrestr = b.eqrestr, 
                   b.restr   = b.restr, 
                   scale     = scale, 
                   cc        = ifelse(is.null(cc), 4.685061, cc))
  } else if (type == "B") {
    if (meq.alt == 0L) {
      Ts <- robustFm(x         = X, 
                     y         = y, 
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
                       b.eqrestr = b.restr, 
                       b.restr   = b.restr.alt, 
                       scale     = scale, 
                       cc        = ifelse(is.null(cc), 4.685061, cc))
      } else {
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
              model.org = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}



conTestWald.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", B = 9999, 
                               p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                               cl = NULL, seed = 1234, verbose = FALSE,
                               control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conRLM" %in% class(object))) {
    stop("object must be of class conRLM.")
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
    #fit inequality constrained robust model
    call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg,
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
  } else if (type == "A") {
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
                                          test     = "Wald",
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
              model.org = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}




conTestScore.conRLM <- function(object, type = "A", neq.alt = 0, boot = "no", B = 9999, 
                                p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                                cl = NULL, seed = 1234, verbose = FALSE,
                                control = NULL, ...) {
  
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!("conRLM" %in% class(object))) {
    stop("object must be of class conRLM.")
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
    #fit inequality constrained robust model
    call.my <- list(Amat = Amatg, meq = nrow(Amatg), bvec = bvecg,
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
  } else if (type == "A") {
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
              model.org = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}