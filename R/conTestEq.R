conTestEq.lm <- function(object, test = "default", ...) {
  
  #  CON <- object$CON
  #  if (is.null(CON)) { stop("not (yet) possible when constraints input is numerical.") }
  
  test <- tolower(test)
  if (#length(CON$ceq.linear.idx)     > 0  && # some linear eq. constraints
    #length(CON$ceq.nonlinear.idx) == 0L && # no nonlinear eq. constraints
    #length(CON$cin.linear.idx)    == 0L && # no inequality constraints
    #length(CON$cin.nonlinear.idx) == 0L
    nrow(object$Amat) == object$meq) {
    
    if (test == "default") {
      test <- "wald"
    }
    
    # here we perform the usual Wald/F test...
    if (test == "wald" || test == "x2" || test == "chisq") {
      OUT <- con_test_Wald(Sigma    = object$Sigma,
                           JAC     = object$Amat, #CON$ceq.JAC,
                           theta.r = c(object$Amat %*% object$b.unrestr)) #CON$ceq.theta)
    OUT$Amat <- object$Amat
    OUT$bvec <- object$bvec
    OUT$meq <- object$meq
    OUT$b.restr <- object$b.restr
    OUT$b.unrestr <- object$b.unrestr
    } else if (test == "f" || test == "ftest") {
      Wald <- con_test_Wald(Sigma    = object$Sigma,
                            JAC     = object$Amat, #CON$ceq.JAC,
                            theta.r = c(object$Amat %*% object$b.unrestr)) #CON$ceq.theta))
      # convert Wald to F
      OUT <- list()
      OUT$test <- "F"
      OUT$Ts <- Wald$Ts / Wald$df
      OUT$df   <- Wald$df
      #X <- model.matrix(object$model.org)[,,drop=FALSE]
      OUT$df.residual <- df.residual(object) #nrow(X)-ncol(X)
      OUT$pvalue <- 1 - pf(OUT$Ts, OUT$df, OUT$df.residual)
      OUT$Amat <- object$Amat
      OUT$bvec <- object$bvec
      OUT$meq <- object$meq
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else {
      stop("restriktor ERROR: test ", sQuote(test), " not (yet) implemented. Choose \"F\" or \"Wald\"")
    }
  } else if (#length(CON$ceq.nonlinear.idx) == 0L &&
    #length(CON$cin.linear.idx)     > 0L && # some inequalities restr.
    #length(CON$cin.nonlinear.idx) == 0L
    nrow(object$Amat != object$meq)) {
    
    stop("test not applicable with inequality constraints.")
  } else if (length(CON$ceq.nonlinear.idx) > 0L ||
               length(CON$cin.nonlinear.idx) > 0L) {
    stop("ERROR: can not handle (yet) nonlinear (in)equality constraints")
  }

  class(OUT) <- "conTest"

  OUT
}



