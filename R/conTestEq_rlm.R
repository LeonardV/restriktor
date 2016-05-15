conTestEq.rlm <- function(object, test = "default", ...) {
  
  if(!any(weights(object) == 1)) {
    stop("weights not supported (yet).")
  }
  
  model.org <- object$model.org
  y <- as.matrix(object$model.org$model[, attr(object$model.org$terms, "response")])
  # model matrix
  X <- model.matrix(object)[,,drop = FALSE]
  # tukey's bisquare tuning constant
  cc <- model.org$call[["c"]]
  # coefficients under equality restriktions
  beta0 <- coef(object)
  # unrestrikted coefficients
  betaA <- coef(model.org)
  # unrestrikted scale estimate
  scale <- model.org$s
  # restriktion stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  
  test <- tolower(test)
  if (nrow(Amat) == meq) {
    
    if (test == "default") {
      test <- "f"
    }
    
    # here we perform the usual Wald/F test...
    if (test == "wald" || test == "x2" || test == "chisq" || test == "score") {
      OUT <- robustWaldScores(x = X, y = y, beta0 = beta0, betaA = betaA, 
                              scale = scale, test = test, 
                              cc = ifelse(is.null(cc), 4.685061, cc))
      # df
      OUT$df <- nrow(Amat)
      OUT$Ts <- OUT$Ts / OUT$df
      # rdf
      OUT$df.residual <- df.residual(object) 
      # p-value based on chisq
      OUT$pvalue <- 1 - pf(OUT$Ts, OUT$df, OUT$df.residual)
      OUT$Amat <- Amat
      OUT$bvec <- bvec
      OUT$meq  <- meq
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "f" || test == "ftest") {
      Fmm <- robustFm(x = X, y = y, beta0 = beta0, betaA = betaA, 
                      scale = scale, cc = ifelse(is.null(cc), 4.685061, cc))
      
      OUT <- list()
      OUT$test <- "F"
      # df
      OUT$df <- nrow(Amat)
      OUT$Ts <- Fmm / OUT$df
      # rdf
      OUT$df.residual <- df.residual(object) 
      # p-value based on chisq
      OUT$pvalue <- 1 - pf(OUT$Ts, OUT$df, OUT$df.residual)
      OUT$Amat <- Amat
      OUT$bvec <- bvec
      OUT$meq  <- meq
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
  
  class(OUT) <- "conTest"
  
  OUT
}
