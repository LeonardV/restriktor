# standard bootstrap
con_bootdata_lm <- function(data, wt, indices, ...) {
  l <- list(...)
  CALL <- as.list(l$CALL)
  CALL[[1]] <- NULL
  # the data obtained via object$model may contains an extra column
  # with the weights and must be removed first, before we can continue.
  if (!is.null(CALL$weights)) {
    data <- data[-ncol(data)]
  }
  CALL$data <- data[indices,]
  CALL$weights <- wt 
  l$object <- do.call("lm", CALL)
  
  do.call("restriktor", l)$b.restr
}


con_bootdata_rlm <- function(data, wt, indices, ...) {
  l <- list(...)
  CALL <- as.list(l$CALL)
  CALL[[1]] <- NULL
  if (!is.null(CALL$weights)) {
    data <- data[-ncol(data)]
  }
  CALL$data <- data[indices,]
  l$object <- do.call("rlm", CALL)
  
  do.call("restriktor", l)$b.restr
}
