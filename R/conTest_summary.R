conTest_summary.conLM <- function(object, test = "F", ...) {
  
  x <- object
  if (!(inherits(object, c("conLM","conRLM")))) {
    stop("Restriktor ERROR: object must be of class \"conLM\" or \"conRLM\"")
  }
  
  ldots <- list(...)
  CALL <- c(list(object = x, test = test), ldots)
  
  OUT <- list()
  # fit all available hypothesis tests
  CALL$type <- "global"
  OUT1 <- do.call("conTest", CALL)
  CALL$type <- "A"
  OUT2 <- do.call("conTest", CALL)
  CALL$type <- "B"
  OUT3 <- do.call("conTest", CALL)
  
  OUT <- c(OUT1, OUT2, OUT3)
  
  if (x$neq == 0) {
    CALL$type <- "C"
    OUT4 <- do.call("conTest", CALL)
    OUT <- c(OUT, OUT4)
  }
  
  class(OUT) <- c("conTest")
  
  OUT
}