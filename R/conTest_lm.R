# computes the F, LRT and Score test statistic
# To do: implement E-bar test statistic.

# REF: Silvapulle and Sen (2005). Constrained statistical inference. Chapter 2.
conTestF.lm <- function(object, type = "A", boot = "no", meq.alt = 0,
                        control = NULL, tol = sqrt(.Machine$double.eps), ...) {

  if (type != "global") {
    type <- toupper(type)
  }  
  
  if (!("conLM" %in% class(object))) {
    stop("object must be of class conLM().")
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

  model.org <- object$model.org
  X <- model.matrix(object)[,,drop=FALSE]
  Y <- model.org$model[, attr(model.org$terms, "response")]
  cov <- vcov(model.org) 
  b.unconstr <- object$b.unconstr
  vnames <- names(b.unconstr)
  b.constr <- object$b.constr
  b.eqconstr <- NULL
  b.constr.alt <- NULL
  Ts <- as.numeric(NA)
  Amat <- object$Amat
  bvec <- object$bvec
  meq <- object$meq

  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality constraints only.")
  }
  if (!(length(b.unconstr) == nrow(cov))) {
    stop("length b.unconstr and nrow(cov) must be identical.")
  }

  if (type == "global") {
    # check for intercept
    intercept <- any(attr(terms(model.org), "intercept"))
    g <- length(object$b.constr)
    if (intercept) {
      Amatg <- cbind(rep(0, (g - 1)), diag(rep(1, g - 1))) 
      bvecg <- rep(0, g - 1) 
    } else {
        stop("Restriktor ERROR: test not ready yet for models without intercept.")      
    } 
    b.eqconstr <- con_solver(b.unconstr, X = X, y = Y, Amat = Amatg,
                                bvec = bvecg, meq = nrow(Amatg),
                                tol = ifelse(is.null(control$tol), 1e-09, 
                                             control$tol),
                                maxit = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    Ts <- c(t(b.constr - b.eqconstr) %*% solve(cov, b.constr - b.eqconstr))
  } else if (type == "A") {
    # optimizer
    b.eqconstr <- con_solver(b.unconstr, X = X, y = Y, Amat = Amat,
                                bvec = bvec, meq = nrow(Amat),
                                tol = ifelse(is.null(control$tol), 1e-09, 
                                             control$tol),
                                maxit = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
      names(b.eqconstr) <- vnames
    Ts <- c(t(b.constr - b.eqconstr) %*% solve(cov, b.constr - b.eqconstr))
  }
  else if (type == "B") {
    if (meq.alt == 0L) {
      # Fbar test statistic
      Ts <- as.vector(t(b.unconstr - b.constr) %*% solve(cov, b.unconstr - b.constr))
    }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if(meq.alt != 0L && meq.alt <= meq) {
        b.constr.alt <- con_solver(b.unconstr, X = X, y = Y,
                                   Amat = Amat[1:meq.alt,,drop=FALSE],
                                   bvec = bvec[1:meq.alt], meq = meq.alt,
                                   tol = ifelse(is.null(control$tol), 1e-09, 
                                                control$tol),
                                   maxit = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$solution
        names(b.constr.alt) <- vnames
        Ts <- as.vector(t(b.constr - b.constr.alt) %*% solve(cov, b.constr - b.constr.alt))
      }
      else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  } else if (type == "C") { # intersection-union test (Sasabuchi, 1980)
    if (meq == 0L) {
      Ts <- as.vector(min((Amat %*% b.unconstr - bvec) /
                            sqrt(diag(Amat %*%cov%*% t(Amat)))))
    } else {
      stop("test not applicable with equality constraints.")
    }
  }

  if (boot == "no") {
    pvalue <- con_pvalue_Fbar(cov, Ts.org = Ts, object$df.residual, type = type,
                              Amat = Amat, bvec = bvec, meq = meq, meq.alt = meq.alt)
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
    pvalue <- con_pvalue_boot_weights(object, type = type, pbar = "pfbar", Ts.org = Ts, 
                                      df.residual = object$df.residual, 
                                      meq.alt = meq.alt,
                                      R = ifelse(is.null(control$B), 9999, control$B),
                                      parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                      ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                      cl = ifelse(is.null(control$cl), NULL, control$cl),
                                      seed = ifelse(is.null(control$seed), 1234, control$seed),
                                      verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
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
              df.residual = object$df.residual,
              cov = cov,
              Ts = Ts,
              pvalue = pvalue,
              model.org = model.org)


  if (type == "A" | type == "global") { 
    OUT$b.eqconstr <- b.eqconstr 
  }

  class(OUT) <- "conTest"

  OUT

}


# REF: Silvapulle and Sen (2005). Constrained statistical inference. Chapter 3.
conTestLRT.lm <- function(object, type = "A", boot = "no", meq.alt = 0,
                           control = NULL, tol = sqrt(.Machine$double.eps), ...) {

  if (type != "global") {
    type <- toupper(type)
  }  
  
  if (!("conLM" %in% class(object))) {
    stop("object must be of class conLM().")
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

  model.org <- object$model.org
  X <- model.matrix(model.org)[,,drop=FALSE]
  Y <- cbind(model.org$model[, attr(model.org$terms, "response")])
  cov <- vcov(model.org)   
  b.unconstr <- object$b.unconstr
    vnames <- names(b.unconstr)
  b.constr <- object$b.constr
  b.eqconstr <- NULL
  b.constr.alt <- NULL
  Ts <- as.numeric(NA)
    names(Ts) <- "LRT"
  Amat <- object$Amat
  bvec <- object$bvec
  meq <- object$meq

  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality constraints only.")
  }
  if (!(length(b.unconstr) == nrow(cov))) {
    stop("length b.unconstr and nrow(cov) must be identical.")
  }

  if (type == "global") {
    # check for intercept
    intercept <- object$model.org$assign[1] == 0L
    g <- length(object$b.constr)
    if (intercept) {
      Amatg <- cbind(rep(0, (g - 1)), diag(rep(1, g - 1))) 
      bvecg <- rep(0, g - 1) 
    } else {
      Amatg <- diag(rep(1, g))
      bvecg <- rep(0, g) 
    } 
    b.eqconstr <- con_solver(b.unconstr, X = X, y = Y, Amat = Amatg,
                             bvec = bvecg, meq = nrow(Amatg),
                             tol = ifelse(is.null(control$tol), 1e-09, 
                                          control$tol),
                             maxit = ifelse(is.null(control$maxit), 1e04, 
                                            control$maxit))$solution
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    # test statistics
    ll0.out <- con_loglik_lm(X = X, y = Y, b = b.eqconstr, detU = 1)
    ll0 <- ll0.out$loglik
    
    ll1.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
    ll1 <- ll1.out$loglik
    Ts <- -2*(ll0 - ll1)
  } else if (type == "A") {
    # optimizer
    b.eqconstr <- con_solver(b.unconstr, X = X, y = Y, Amat = Amat,
                                bvec = bvec, meq = nrow(Amat),
                                tol = ifelse(is.null(control$tol), 1e-09, 
                                             control$tol),
                                maxit = ifelse(is.null(control$maxit), 1e04, 
                                               control$maxit))$solution
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
      names(b.eqconstr) <- vnames
    # test statistics
    ll0.out <- con_loglik_lm(X = X, y = Y, b = b.eqconstr, detU = 1)
    ll0 <- ll0.out$loglik
    
    ll1.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
    ll1 <- ll1.out$loglik
    Ts <- -2*(ll0 - ll1)
  }
  else if (type == "B") {
    if (meq.alt == 0L) {
      # test statistic
      ll0.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
      ll0 <- ll0.out$loglik
    
      ll1.out <- con_loglik_lm(X = X, y = Y, b = b.unconstr, detU = 1)
      ll1 <- ll1.out$loglik
      Ts <- -2*(ll0 - ll1)
    }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        b.constr.alt <- con_solver(b.unconstr, X = X, y = Y,
                                   Amat = Amat[1:meq.alt,,drop=FALSE],
                                   bvec=bvec[1:meq.alt], meq = meq.alt,
                                   tol = ifelse(is.null(control$tol), 1e-09, 
                                                control$tol),
                                   maxit = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$solution
        names(b.constr.alt) <- vnames

        ll0.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
        ll0 <- ll0.out$loglik
        
        ll1.out <- con_loglik_lm(X = X, y = Y, b = b.constr.alt, detU = 1)
        ll1 <- ll1.out$loglik
        Ts <- -2*(ll0 - ll1)
      }
      else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  }
  # intersection-union test (Sasabuchi, 1980)
  else if (type == "C") {
    if (meq == 0L) {
      Ts <- as.vector(min((Amat %*% b.unconstr - bvec) /
                            sqrt(diag(Amat %*%cov%*% t(Amat)))))
      names(Ts) <- "Tbar"
    }
    else {
      stop("test not applicable with equality constraints.")
    }
  }

  if (boot == "no") {
    pvalue <- con_pvalue_Fbar(cov, Ts.org = Ts, object$df.residual, type = type,
                              Amat, bvec, meq, meq.alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, Ts.org = Ts, type = type, test = "LRT",
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
    pvalue <- con_pvalue_boot_model_based(object, Ts.org = Ts, type = type, test = "LRT",
                                          meq.alt = meq.alt,
                                          R = ifelse(is.null(control$B), 9999, control$B),
                                          parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                          ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                          cl = ifelse(is.null(control$cl), NULL, control$cl),
                                          seed = ifelse(is.null(control$seed), 1234, control$seed),
                                          verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  } else if (boot == "mix.weights") {
    pvalue <- con_pvalue_boot_weights(object, type = type, pbar = "pfbar", Ts.org = Ts, df.residual = object$df.residual, 
                                      meq.alt = meq.alt,
                                      R = ifelse(is.null(control$B), 9999, control$B),
                                      parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                      ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                      cl = ifelse(is.null(control$cl), NULL, control$cl),
                                      seed = ifelse(is.null(control$seed), 1234, control$seed),
                                      verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
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
              df.residual = object$df.residual,
              cov = cov,
              Ts = Ts,
              pvalue = pvalue,
              model.org = object$model.org)

  if(type == "A" | type == "global") { 
    OUT$b.eqconstr <- b.eqconstr 
  }

  class(OUT) <- "conTest"

  OUT

}


# REF: Silvapulle, M.J. and Silvapulle, P. (1995). A Score Test Against One-Sided Alternatives
# Journal of the American Statistical Association, Vol. 90, No. 429 (Mar., 1995), pp. 342-349
conTestScore.lm <- function(object, type = "A", boot = "no", meq.alt = 0,
                            control = NULL, tol = sqrt(.Machine$double.eps), ...) {
  
  if (type != "global") {
    type <- toupper(type)
  }  
  
  if (!("conLM" %in% class(object))) {
    stop("object must be of class conLM().")
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
  
  model.org <- object$model.org
  X <- model.matrix(object)[,,drop=FALSE]
  Y <- cbind(model.org$model[, attr(model.org$terms, "response")])
  n <- dim(X)[1]
  p <- dim(X)[2]
  cov <- object$Sigma
  b.unconstr <- object$b.unconstr
  vnames <- names(b.unconstr)
  b.constr <- object$b.constr
  b.eqconstr <- NULL
  b.constr.alt <- NULL
  Ts <- as.numeric(NA)
  names(Ts) <- "Score"
  
  Amat <- object$Amat
  bvec <- object$bvec
  meq <- object$meq
  
  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality constraints only.")
  }

  if (type == "global") {
    # check for intercept
    intercept <- object$model.org$assign[1] == 0L
    g <- length(object$b.constr)
    if (intercept) {
      Amatg <- cbind(rep(0, (g - 1)), diag(rep(1, g - 1))) 
      bvecg <- rep(0, g - 1) 
    } else {
      Amatg <- diag(rep(1, g))
      bvecg <- rep(0, g) 
    }
    
    b.eqconstr <- con_solver(b.unconstr, X = X, y = Y, Amat = Amatg,
                             bvec = bvecg, meq = nrow(Amatg),
                             tol = ifelse(is.null(control$tol), 1e-09, 
                                          control$tol),
                             maxit = ifelse(is.null(control$maxit), 1e04, 
                                            control$maxit))$solution
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    
    s20  <- sum((Y - X%*%b.eqconstr)^2) / (n-p) #no df adjustment necessary for global test.
    d0   <- 1/s20 * t(X)%*%(Y - X%*%b.eqconstr)
    i <- 1/s20 * (t(X)%*%X)
    U <- 1/sqrt(n) * solve(i) %*% d0
    UI <- t(U) %*% i
    D <- i
    b <- quadprog:::solve.QP(Dmat=D , dvec=UI, Amat=t(Amat), bvec=bvec, meq=meq)$solution
    Ts <- t(U) %*% i %*% U - ( t(U-b) %*% i %*% (U-b) ) 
    Ts <- as.numeric(n*Ts)
  } else if (type == "A") {
    
    b.eqconstr <- con_solver(b.unconstr, X = X, y = Y, Amat = Amat,
                             bvec = bvec, meq = nrow(Amat),
                             tol = ifelse(is.null(control$tol), 1e-09, 
                                          control$tol),
                             maxit = ifelse(is.null(control$maxit), 1e04, 
                                            control$maxit))$solution
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
    names(b.eqconstr) <- vnames
    
    # if parameters are constraint to a number, df shoulde be adjusted.
    con.idx <- which(object$partable$op == "<" | object$partable$op == ">")
    user.equal  <- object$partable
      user.equal$op[con.idx] <- "=="
    pEq.corr <- dfEq_correction(user.equal)
    p <- NCOL(X)
    df <- (n - p + pEq.corr)
    s20 <- sum((Y - X%*%b.eqconstr)^2) / df
    
    cat("CHECK DF S^2! =", s20, "...df = ", df, "\n")
    
    d0 <- 1/s20 * t(X)%*%(Y - X%*%b.eqconstr)
    i <- 1/s20 * (t(X)%*%X)
    U <- 1/sqrt(n) * solve(i) %*% d0
    UI <- t(U) %*% i
    D <- i
    b <- quadprog:::solve.QP(Dmat = D , dvec = UI, Amat = t(Amat), bvec = bvec, 
                             meq = meq)$solution
    Ts <- t(U) %*% i %*% U - ( t(U-b) %*% i %*% (U-b) ) 
    Ts <- as.numeric(n*Ts)
  }
  else if (type == "B") {
    if (meq.alt == 0L) {
      pEq.corr <- dfEq_correction(object$partable)
      p <- NCOL(X)
      df <- (n - p + pEq.corr)
      s20 <- sum((Y - X%*%b.constr)^2) / df
      
      cat("CHECK DF S^2! =", s20, "...df = ", df, "\n")
      
      d0   <- 1/s20 * t(X)%*%(Y - X%*%b.constr)
      i <- 1/s20 * (t(X)%*%X)
      U <- 1/sqrt(n) * solve(i) %*% d0
      UI <- t(U) %*% i
      Ts <- t(U) %*% i %*% U 
      Ts <- as.numeric(n*Ts)
    }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if(meq.alt != 0L && meq.alt <= meq) {
        b.constr.alt <- con_solver(b.unconstr, X = X, y = Y,
                                   Amat = Amat[1:meq.alt,,drop=FALSE],
                                   bvec = bvec[1:meq.alt], meq = meq.alt,
                                   tol = ifelse(is.null(control$tol), 1e-09, 
                                                control$tol),
                                   maxit = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$solution
        names(b.constr.alt) <- vnames
        
        pEq.corr <- dfEq_correction(object$partable)
        p <- NCOL(X)
        df <- (n - p + pEq.corr)
        s20 <- sum((Y - X%*%b.constr)^2) / df
        
        cat("CHECK DF S^2! =", s20, "...df = ", df, "\n")
        
        d0   <- 1/s20 * t(X)%*%(Y - X%*%b.constr)
        i <- 1/s20 * (t(X)%*%X)
        U <- 1/sqrt(n) * solve(i) %*% d0
        UI <- t(U) %*% i
        D <- i
        b <- quadprog:::solve.QP(Dmat = D, dvec = UI, 
                                 Amat = t(Amat[1:meq.alt, , drop = FALSE]), 
                                 bvec = bvec[1:meq.alt], meq = meq.alt)$solution
        Ts <- t(U) %*% i %*% U - ( t(U-b) %*% i %*% (U-b) ) 
        Ts <- as.numeric(n*Ts)
      }
      else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  }
  # intersection-union test (Sasabuchi, 1980)
  else if (type == "C") {
    if (meq == 0L) {
      Ts <- as.vector(min((Amat %*% b.unconstr - bvec) /
                            sqrt(diag(Amat %*%cov%*% t(Amat)))))
      names(Ts) <- "Tbar"
    }
    else {
      stop("test not applicable with equality constraints.")
    }
  }
  
  if (boot == "no") {
    pvalue <- con_pvalue_Chibar(cov, Ts.org = Ts, object$df.residual, type = type,
                                Amat, bvec, meq, meq.alt)
  }
  else if (boot == "parametric") {
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
  }
  else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based(object, Ts.org = Ts, type = type, test = "Score",
                                          meq.alt = meq.alt,
                                          R = ifelse(is.null(control$B), 9999, control$B),
                                          parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                          ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                          cl = ifelse(is.null(control$cl), NULL, control$cl),
                                          seed = ifelse(is.null(control$seed), 1234, control$seed),
                                          verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  } else if (boot == "mix.weights") {
    pvalue <- con_pvalue_boot_weights(object, type = type, pbar = "pfbar", Ts.org = Ts, df.residual = object$df.residual,
                                      meq.alt = meq.alt,
                                      R = ifelse(is.null(control$B), 9999, control$B),
                                      parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                      ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                      cl = ifelse(is.null(control$cl), NULL, control$cl),
                                      seed = ifelse(is.null(control$seed), 1234, control$seed),
                                      verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
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
              df.residual = object$df.residual,
              cov = cov,
              Ts = Ts,
              pvalue = pvalue,
              model.org = object$model.org)
  
  if(type == "A" | type == "global") { 
    OUT$b.eqconstr <- b.eqconstr 
  }
  
  class(OUT) <- "conTest"
  
  OUT
  
}
