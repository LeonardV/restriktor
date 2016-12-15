conTestEq.conRLM <- function(object, test = "F", boot = "no", 
                             R = 9999, p.distr = "N", df = 7, 
                             parallel = "no", ncpus = 1L, cl = NULL, 
                             seed = 1234, verbose = FALSE, ...) {
  
  if (!("conRLM" %in% class(object))) {
    stop("object must be of class conRLM.")
  }
  if(!any(weights(object) == 1)) {
    stop("weights not supported (yet).")
  }
  
  test <- tolower(test)
  stopifnot(test %in% c("f","wald","wald2","score"))
  
  model.org <- object$model.org
  y <- as.matrix(object$model.org$model[, attr(object$model.org$terms, "response")])
  # model matrix
  X <- model.matrix(object)[,,drop = FALSE]
  # tukey's bisquare tuning constant
  cc <- model.org$call[["c"]]
  # coefficients under equality restriktions
  b.eqrestr <- coef(object)
  # unrestrikted coefficients
  b.unrestr <- coef(model.org)
  # unrestrikted scale estimate
  scale <- model.org$s
  # scale estimate used for the standard errors
  tau <- object$s2.unc
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
      F_out <- robustFm(x         = X, 
                        y         = y, 
                        b.unrestr = b.unrestr, 
                        b.eqrestr = b.eqrestr, 
                        b.restr   = b.unrestr, 
                        scale     = scale, 
                        cc        = ifelse(is.null(cc), 4.685061, cc))
          
      OUT <- append(CON, F_out)
      # df
      OUT$df <- nrow(Amat)
      OUT$Ts <- OUT$Ts / OUT$df
      # rdf
      OUT$df.residual <- df.residual(object) 
      # p-value based on chisq
      OUT$pvalue <- 1 - pf(OUT$Ts, OUT$df, OUT$df.residual)
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "wald" || test == "score") {
      WaldScore_out <- robustWaldScores(x         = X, 
                                        y         = y, 
                                        b.eqrestr = b.eqrestr, 
                                        b.restr   = b.unrestr, 
                                        b.unrestr = b.unrestr,
                                        scale     = scale, 
                                        test      = test, 
                                        cc        = ifelse(is.null(cc), 4.685061, cc))
      OUT <- append(CON, WaldScore_out)
      OUT$df <- nrow(Amat)
      OUT$df.residual <- df.residual(object) 
      # p-value based on chisq
      OUT$pvalue <- 1 - pchisq(OUT$Ts, df = OUT$df) 
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "wald2") {  
      Wald2_out <- robustWaldXX(x         = X, 
                                b.eqrestr = b.eqrestr, 
                                b.restr   = b.unrestr, 
                                b.unrestr = b.unrestr, 
                                tau       = tau)
      OUT <- append(CON, Wald2_out)
      
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
    OUT$pvalue <- con_pvalue_boot_parametric(object, 
                                             Ts.org   = OUT$Ts, 
                                             type     = "A",
                                             test     = test, 
                                             R        = R, 
                                             p.distr  = p.distr,
                                             df       = df, 
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
  
  OUT$R2.org      <- object$R2.org
  OUT$R2.reduced  <- object$R2.reduced
  
  OUT <- list(OUT)
  names(OUT) <- "ceq"
  
  class(OUT) <- "conTest"
  
  OUT
}
