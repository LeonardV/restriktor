conTest_ceq.conRLM <- function(object, test = "F", boot = "no", 
                               R = 9999, p.distr = rnorm,  
                               parallel = "no", ncpus = 1L, cl = NULL, 
                               seed = 1234, verbose = FALSE, ...) {
  
  if (!inherits(object, "conRLM")) {
    stop("object must be of class conRLM.")
  }
  
  test <- tolower(test)
  stopifnot(test %in% c("f","wald","score"))
  
  model.org <- object$model.org
  y <- as.matrix(object$model.org$model[, attr(object$model.org$terms, "response")])
  # model matrix
  X <- model.matrix(object)[,,drop = FALSE]
  # tukey's bisquare tuning constant
  cc <- model.org$call[["c"]]
  # coefficients under equality restriktions
  b.eqrestr <- object$b.restr
  # unrestrikted coefficients
  b.unrestr <- coef(model.org)
  # unrestrikted scale estimate
  scale <- model.org$s
  # scale estimate used for the standard errors
  #scale <- summary(model.org)$stddev
  # restriktion stuff
  CON  <- object$CON
  
  CON$Amat <- Amat <- object$constraints
  CON$bvec <- object$rhs
  CON$meq  <- meq <- object$neq
  
  if (nrow(Amat) == meq) {
    
    if (test == "default") {
      test <- "f"
    }
    
    # here we perform the usual Wald/F test...
    if (test == "f") {
      F.out <- robustFm(x         = X, 
                        y         = y, 
                        b.unrestr = b.unrestr, 
                        b.eqrestr = b.eqrestr, 
                        b.restr   = b.unrestr, 
                        scale     = scale, 
                        cc        = ifelse(is.null(cc), 4.685061, cc))
          
      OUT <- append(CON, F.out)
      # df
      OUT$df <- nrow(Amat)
      OUT$Ts <- OUT$Ts / OUT$df
      # rdf
      OUT$df.residual <- df.residual(object) 
      # p-value based on chisq
      OUT$pvalue <- 1 - pf(OUT$Ts, OUT$df, OUT$df.residual)
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "score") {
      Score.out <- robustScores(x         = X, 
                                y         = y, 
                                b.eqrestr = b.eqrestr, 
                                b.restr   = b.unrestr, 
                                b.unrestr = b.unrestr,
                                Amat      = Amat,
                                scale     = scale, 
                                test      = test, 
                                cc        = ifelse(is.null(cc), 4.685061, cc))
      OUT <- append(CON, Score.out)
      OUT$df <- nrow(Amat)
      OUT$df.residual <- df.residual(object) 
      # p-value based on chisq
      OUT$pvalue <- 1 - pchisq(OUT$Ts, df = OUT$df) 
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "wald") {  
      Wald.out <- robustWald(x         = X, 
                             y         = y,
                             b.eqrestr = b.eqrestr, 
                             b.restr   = b.unrestr, 
                             b.unrestr = b.unrestr,
                             scale     = scale,
                             Amat      = Amat,
                             cc        = ifelse(is.null(cc), 4.685061, cc))
      OUT <- append(CON, Wald.out)
      
      OUT$df <- nrow(Amat)
      OUT$pvalue <- 1 - pchisq(OUT$Ts, df = OUT$df) 
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr 
    } else {
      stop("restriktor ERROR: test ", sQuote(test), " not (yet) implemented.")
    }
  } else if (nrow(Amat != meq)) {
    stop("test not applicable with inequality constraints.")
  } else if (length(CON$ceq.nonlinear.idx) > 0L ||
             length(CON$cin.nonlinear.idx) > 0L) {
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
