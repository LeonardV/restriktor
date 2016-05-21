con_bootdata_lm <- function(data, indices, ...) {
  l <- list(...)
  CALL <- as.list(l$CALL)
  CALL[[1]] <- NULL
  CALL$data <- data[indices,]
  l$model <- do.call("lm", CALL)
  #l$model <- lm(form, data = dat)
  OUT <- do.call("restriktor", l)$b.restr

  OUT
}


con_bootdata_rlm <- function(data, indices, ...) {
  l <- list(...)
  CALL <- as.list(l$CALL)
  CALL[[1]] <- NULL
  CALL$data <- data[indices,]
  l$model <- do.call("rlm", CALL)
  #l$model <- rlm(form, data = dat)
  OUT <- do.call("restriktor", l)$b.restr
  
  OUT
}
