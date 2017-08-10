conTest_summary.restriktor <- function(object, test = "F", ...) {
  
  if (!(inherits(object, "restriktor"))) {
    stop("Restriktor ERROR: object must be of class restriktor")
  }
  
  Amat <- object$constraints
  meq  <- object$neq
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
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
  
  class(OUT) <- c("conTest")
  
  OUT
}