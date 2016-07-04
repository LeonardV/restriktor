con_bootdata_lm <- function(data, indices, ...) {
  l <- list(...)
  CALL <- as.list(l$CALL)
  CALL[[1]] <- NULL
  CALL$data <- data[indices,]
  l$object <- do.call("lm", CALL)
  OUT <- do.call("restriktor", l)$b.restr

  OUT
}


con_bootdata_rlm <- function(data, indices, ...) {
  l <- list(...)
  CALL <- as.list(l$CALL)
  CALL[[1]] <- NULL
  CALL$data <- data[indices,]
  l$object <- do.call("rlm", CALL)
  OUT <- do.call("restriktor", l)$b.restr
  
  OUT
}
