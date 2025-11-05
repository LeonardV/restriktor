conTest_ceq.conLM <- function(object, test = "F", boot = "no", 
                              R = 9999, p.distr = rnorm, 
                              parallel = "no", ncpus = 1L, cl = NULL, 
                              seed = 1234, verbose = FALSE, ...) {
  
  if (!inherits(object, "conLM")) {
    stop("object must be of class conLM.")
  }
  if(!is.null(weights(object))) {
    stop("weights not supported (yet).")
  }
  
  test <- tolower(test)
  stopifnot(test %in% c("f","wald","lrt","score"))
  
  CON  <- object$CON
  CON$Amat <- Amat <- object$constraints
  CON$bvec <- object$rhs
  CON$meq  <- meq <- object$neq
  
  if (#length(CON$ceq.linear.idx)     > 0  && # some linear eq. constraints
    #length(CON$ceq.nonlinear.idx) == 0L && # no nonlinear eq. constraints
    #length(CON$cin.linear.idx)    == 0L && # no inequality constraints
    #length(CON$cin.nonlinear.idx) == 0L
    nrow(Amat) == meq) {
    
    if (test == "default") {
      test <- "f"
    }
    
    # here we perform the usual Wald/F test...
    if (test == "wald") {
      #theta.r <- object$b.unrestr
      Wald.out <- con_test_Wald(Sigma   = object$Sigma,
                                JAC     = Amat,         
                                theta.r = Amat %*% object$b.unrestr - object$rhs) 
      
      OUT <- append(CON, Wald.out)
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "f") {
      Wald.out <- con_test_Wald(Sigma   = object$Sigma,
                                JAC     = Amat,
                                theta.r = Amat%*%object$b.unrestr - object$rhs) 
      # convert Wald to F
      OUT <- append(CON, Wald.out)
      OUT$test <- "F"
      OUT$Ts <- Wald.out$Ts / Wald.out$df
      OUT$df.residual <- object$df.residual
      OUT$pvalue <- 1 - pf(OUT$Ts, OUT$df, OUT$df.residual)
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "score") {
      OUT <- CON
      OUT$test <- "Score"
      # model matrix
      X <- model.matrix(object)[,,drop = FALSE]
      n <- dim(X)[1]
      #p <- dim(X)[2]
      
      # weights
      if (is.null(weights(object))) {
        w <- rep(1, n)
      } else {
        w <- weights(object)  
      }
      
      # residuals under the null-hypothesis
      res0 <- residuals(object, "working")
      # degrees-of-freedom under the null-hypothesis
      #df0 <- (n - (p - qr(Amat[0:meq,,drop = FALSE])$rank))
      # sigma^2
      s20 <- object$s2#sum(res0^2) / df0
      
      # information matrix under the null-hypothesis
      I0 <- object$information
      # score vector
      G0 <- colSums(as.vector(res0) * w * X) / s20
      # score test-statistic
      OUT$Ts <- G0 %*% solve(I0, G0)
      # df
      OUT$df <- nrow(Amat)
      # p-value based on chisq
      OUT$pvalue    <- 1 - pchisq(OUT$Ts, df = OUT$df)
      OUT$b.restr   <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "lrt") {
      OUT <- CON
      OUT$test <- "LRT"
      ll0 <- object$loglik
      ll1 <- logLik(object$model.org)
      OUT$Ts <- -2*(ll0 - ll1)
      OUT$df <- nrow(Amat)
      OUT$pvalue    <- 1 - pchisq(OUT$Ts, df = OUT$df)
      OUT$b.restr   <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
      } else {
      stop("restriktor ERROR: test ", sQuote(test), " not (yet) implemented.")
    }
  } else if (#length(CON$ceq.nonlinear.idx) == 0L &&
    #length(CON$cin.linear.idx)     > 0L && # some inequalities restr.
    #length(CON$cin.nonlinear.idx) == 0L
    nrow(Amat != meq)) {
    stop("test not applicable with inequality constraints.")
  } else if (length(CON$ceq.nonlinear.idx) > 0L || length(CON$cin.nonlinear.idx) > 0L) {
    stop("ERROR: can not handle (yet) nonlinear (in)equality constraints")
  }

  OUT$boot <- boot
  
  if (boot == "parametric") {
    if (!is.function(p.distr)) {
      p.distr <- get(p.distr, mode = "function")
    }
    arguments <- list(...)
    pnames <- names(formals(p.distr))
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    pm <- names(arguments)[pm > 0L]
    formals(p.distr)[pm] <- unlist(arguments[pm])
    
    OUT$pvalue <- con_pvalue_boot_parametric(object, 
                                             Ts.org   = OUT$Ts, 
                                             type     = "A",
                                             test     = test, 
                                             R        = R, 
                                             p.distr  = p.distr,
                                             parallel = parallel,
                                             ncpus    = ncpus, 
                                             cl       = cl,
                                             seed     = seed, 
                                             verbose  = verbose)
  } else if (boot == "model.based") {
    OUT$pvalue <- con_pvalue_boot_model_based(object, 
                                              Ts.org   = OUT$Ts, 
                                              type     = "A", 
                                              test     = test,
                                              R        = R,
                                              parallel = parallel, 
                                              ncpus    = ncpus, 
                                              cl       = cl, 
                                              seed     = seed, 
                                              verbose  = verbose)
  } 
  
  OUT$R2.org <- object$R2.org
  OUT$R2.reduced  <- object$R2.reduced
  OUT$model.org <- object$model.org
  
  class(OUT) <- "conTest"

  OUT
}



con_test_Wald <- function(Sigma, JAC, theta.r) {
  
  # remove redundant rows from JAC *and* theta_r
  npar <- ncol(JAC)
  
  JAC.aug <- cbind(JAC, theta.r)
  Q <- qr(t(JAC.aug))
  JAC.full <- t(qr.X(Q, ncol = Q$rank))
  JAC <- JAC.full[,seq_len(npar),drop = FALSE]
  theta.r <- as.numeric(JAC.full[,(npar + 1L)])
  
  # restricted vcov
  info.r  <- JAC %*% Sigma %*% t(JAC)
  
  # Wald test statistic
  Wald <- as.numeric(t(theta.r) %*% solve(info.r) %*% theta.r)
  
  # df
  Wald.df <- nrow(JAC)
  
  # p-value based on chisq
  Wald.pvalue <- 1 - pchisq(Wald, df = Wald.df)
  
  OUT <- list(test   = "Wald",
              Ts     = Wald,
              df     = Wald.df,
              pvalue = Wald.pvalue)
  
  OUT
}


