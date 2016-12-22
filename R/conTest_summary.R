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
  OUT1 <- do.call("conTest", CALL)
  CALL$type <- "A"
  OUT2 <- do.call("conTest", CALL)
  CALL$type <- "B"
  OUT3 <- do.call("conTest", CALL)
  
  OUT <- c(OUT1, OUT2, OUT3)
  
  if (meq == 0) {
    CALL$type <- "C"
    OUT4 <- do.call("conTest", CALL)
    OUT <- c(OUT, OUT4)
  }
  
  class(OUT) <- c("conTest")
  
  OUT
}