conTest_ceq.conRLM <- function(object, test = "F", boot = "no", 
                             R = 9999, p.distr = rnorm,  
                             parallel = "no", ncpus = 1L, cl = NULL, 
                             seed = 1234, verbose = FALSE, ...) {
  
  if (!inherits(object, "conRLM")) {
    stop("object must be of class conRLM.")
  }
  
  test <- tolower(test)
  stopifnot(test %in% c("f","wald","wald2","score"))
  
  model_org <- object$model_org
  y <- as.matrix(object$model_org$model[, attr(object$model_org$terms, "response")])
  # model matrix
  X <- model.matrix(object)[,,drop = FALSE]
  # tukey's bisquare tuning constant
  cc <- model_org$call[["c"]]
  # coefficients under equality restriktions
  b_eqrestr <- coef(object)
  # unrestrikted coefficients
  b_unrestr <- coef(model_org)
  # unrestrikted scale estimate
  scale <- model_org$s
  # scale estimate used for the standard errors
  tau <- sqrt(object$s2_unrestr)
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
                        b_unrestr = b_unrestr, 
                        b_eqrestr = b_eqrestr, 
                        b_restr   = b_unrestr, 
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
      OUT$b_restr <- object$b_restr
      OUT$b_unrestr <- object$b_unrestr
    } else if (test == "wald" || test == "score") {
      WaldScore_out <- robustWaldScores(x         = X, 
                                        y         = y, 
                                        b_eqrestr = b_eqrestr, 
                                        b_restr   = b_unrestr, 
                                        b_unrestr = b_unrestr,
                                        scale     = scale, 
                                        test      = test, 
                                        cc        = ifelse(is.null(cc), 4.685061, cc))
      OUT <- append(CON, WaldScore_out)
      OUT$df <- nrow(Amat)
      OUT$df.residual <- df.residual(object) 
      # p-value based on chisq
      OUT$pvalue <- 1 - pchisq(OUT$Ts, df = OUT$df) 
      OUT$b_restr <- object$b_restr
      OUT$b_unrestr <- object$b_unrestr
    } else if (test == "wald2") {  
      Wald2_out <- robustWaldXX(x         = X, 
                                b_eqrestr = b_eqrestr, 
                                b_restr   = b_unrestr, 
                                b_unrestr = b_unrestr, 
                                tau       = tau)
      OUT <- append(CON, Wald2_out)
      
      OUT$df <- nrow(Amat)
      OUT$pvalue <- 1 - pchisq(OUT$Ts, df = OUT$df) 
      OUT$b_restr <- object$b_restr
      OUT$b_unrestr <- object$b_unrestr 
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
                                             Ts_org   = OUT$Ts, 
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
                                              Ts_org   = OUT$Ts, 
                                              type     = "A", 
                                              test     = test, 
                                              R        = R, 
                                              parallel = parallel, 
                                              ncpus    = ncpus, 
                                              cl       = cl, 
                                              seed     = seed, 
                                              verbose  = verbose)
  } 
  
  OUT$R2_org      <- object$R2_org
  OUT$R2_reduced  <- object$R2_reduced
  OUT$model_org <- object
  
  OUT <- list(OUT)
  names(OUT) <- "ceq"
  
  class(OUT) <- "conTest"
  
  OUT
}
