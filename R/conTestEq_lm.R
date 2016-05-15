conTestEq.lm <- function(object, test = "default", ...) {
  
  if(!is.null(weights(object))) {
    stop("weights not supported (yet).")
  }
  
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  
  test <- tolower(test)
  if (#length(CON$ceq.linear.idx)     > 0  && # some linear eq. constraints
    #length(CON$ceq.nonlinear.idx) == 0L && # no nonlinear eq. constraints
    #length(CON$cin.linear.idx)    == 0L && # no inequality constraints
    #length(CON$cin.nonlinear.idx) == 0L
    nrow(Amat) == meq) {
    
    if (test == "default") {
      test <- "f"
    }
    
    # here we perform the usual Wald/F test...
    if (test == "wald" || test == "x2" || test == "chisq") {
      OUT <- con_test_Wald(Sigma   = object$Sigma,
                           JAC     = Amat, #CON$ceq.JAC,
                           theta.r = c(Amat %*% object$b.unrestr)) #CON$ceq.theta)
      OUT$Amat <- Amat
      OUT$bvec <- bvec
      OUT$meq  <- meq
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "f" || test == "ftest") {
      Wald <- con_test_Wald(Sigma   = object$Sigma,
                            JAC     = Amat, #CON$ceq.JAC,
                            theta.r = c(Amat %*% object$b.unrestr)) #CON$ceq.theta))
      # convert Wald to F
      OUT <- list()
      OUT$test <- "F"
      OUT$Ts <- Wald$Ts / Wald$df
      OUT$df <- Wald$df
      OUT$df.residual <- df.residual(object) 
      OUT$pvalue <- 1 - pf(OUT$Ts, OUT$df, OUT$df.residual)
      OUT$Amat <- Amat
      OUT$bvec <- bvec
      OUT$meq  <- meq
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else if (test == "score") {
      # response variable
      y <- as.matrix(object$model.org$model[, attr(object$model.org$terms, "response")])
      # model matrix
      X <- model.matrix(object)[,,drop = FALSE]
      n <- dim(X)[1]
      p <- dim(X)[2]
      # MSE 
      s20 <- sum((y - X%*%object$b.restr)^2) / (n - (p - qr(Amat[0:meq,,drop = FALSE])$rank))
      # information matrix
      I  <- 1/s20 * t(X) %*% X
      # score vector
      d0 <- as.numeric(1/s20 * t(X) %*% (y - X %*% object$b.restr))
      OUT <- list()
      # score test statistic
      OUT$Ts <- as.numeric(d0 %*% solve(I) %*% d0)
      # df
      OUT$df <- nrow(Amat)
      # p-value based on chisq
      OUT$pvalue <- 1 - pchisq(OUT$Ts, df = OUT$df)
      OUT$Amat <- Amat
      OUT$bvec <- bvec
      OUT$meq  <- meq
      OUT$b.restr <- object$b.restr
      OUT$b.unrestr <- object$b.unrestr
    } else {
      stop("restriktor ERROR: test ", sQuote(test), " not (yet) implemented.")
    }
  } else if (#length(CON$ceq.nonlinear.idx) == 0L &&
    #length(CON$cin.linear.idx)     > 0L && # some inequalities restr.
    #length(CON$cin.nonlinear.idx) == 0L
    nrow(Amat != meq)) {
    stop("test not applicable with inequality constraints.")
  } else if (length(CON$ceq.nonlinear.idx) > 0L ||
               length(CON$cin.nonlinear.idx) > 0L) {
    stop("ERROR: can not handle (yet) nonlinear (in)equality constraints")
  }

  class(OUT) <- "conTest"

  OUT
}



