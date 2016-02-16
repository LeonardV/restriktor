conTestEq <- function(object, test = "default", ...) {

#  CON <- object$CON
#  if (is.null(CON)) { stop("not (yet) possible when constraints input is numerical.") }

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
      OUT <- con_test_Wald(VCOV    = vcov(object$model.org, labels = FALSE),
                           JAC     = object$Amat, #CON$ceq.JAC,
                           theta.r = c(object$Amat %*% coef(object$model.org)) )#CON$ceq.theta)
    } else if (test == "f" || test == "ftest") {
      Wald <- con_test_Wald(VCOV    = vcov(object$model.org, labels = FALSE),
                            JAC     = object$Amat, #CON$ceq.JAC,
                            theta.r = c(object$Amat %*% coef(object$model.org)) )#CON$ceq.theta))
      # convert Wald to F
      OUT <- list()
      OUT$test <- "F"
      OUT$stat <- Wald$stat / Wald$df
      OUT$df   <- Wald$df
      X <- model.matrix(object$model.org)[,,drop=FALSE]
      OUT$df.residual <- nrow(X)-ncol(X)#df.residual(object)
      OUT$p.value <- 1 - pf(OUT$stat, OUT$df, OUT$df.residual)
    } else {
      stop("restriktor ERROR: test ", sQuote(test), " not implemented.")
    }
  } else if (#length(CON$ceq.nonlinear.idx) == 0L &&
            #length(CON$cin.linear.idx)     > 0L && # some inequalities constr.
            #length(CON$cin.nonlinear.idx) == 0L
            nrow(object$Amat != object$meq)) {

    stop("test not applicable with inequality constraints.")
  } else if (length(CON$ceq.nonlinear.idx) > 0L ||
            length(CON$cin.nonlinear.idx) > 0L) {
    stop("ERROR: can not handle (yet) nonlinear (in)equality constraints")
  }

  OUT
}
