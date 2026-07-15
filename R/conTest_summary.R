conTest_summary.restriktor <- function(object, test = "F", ...) {
  
  if (!(inherits(object, "restriktor"))) {
    stop("restriktor ERROR: object must be of class restriktor")
  }
  
  Amat <- object$constraints
  meq  <- object$neq
  if (meq == nrow(Amat)) {
    stop("restriktor ERROR: test not applicable for object with equality restrictions only.")
  }
  
  ldots <- list(...)
  CALL <- c(list(object = object, test = test), ldots)
  
  OUT <- list()
  # fit all available hypothesis tests
  CALL$type <- "global"
  OUT$global <- do.call("conTest", CALL)
  CALL$type <- "A"
  OUT$A <- do.call("conTest", CALL)
  CALL$type <- "B"
  OUT$B <- do.call("conTest", CALL)
  
  #OUT <- list(global = OUT1, A = OUT2, B = OUT3)
  
  if (meq == 0 && !inherits(object, "conMLM")) {
    CALL$type <- "C"
    OUT$C <- do.call("conTest", CALL)
    #OUT <- c(OUT, OUT4)
  }

  # joint tail probability of the type B and type A statistics under the
  # shared (least favorable) null Rb = r (Wolak, 1989, Theorem 1). Only
  # available when the marginal p-values are based on the chi-bar-square
  # weights (no bootstrap) and under the full null (neq.alt = 0).
  wt.bar <- object$wt.bar
  if ((is.null(ldots$boot) || ldots$boot == "no") &&
      (is.null(ldots$neq.alt) || ldots$neq.alt == 0L) &&
      inherits(object, c("conLM", "conGLM", "conRLM")) &&
      !is.null(wt.bar) && !is.null(attr(wt.bar, "method")) &&
      attr(wt.bar, "method") != "none") {
    OUT$joint <- con_pvalue_joint(object, Ts.A = OUT$A$Ts, Ts.B = OUT$B$Ts)
  }

  class(OUT) <- c("conTest")

  OUT
}
