# standard bootstrap
con_bootdata_lm <- function(data, indices, ...) {
  CALL <- list(...)
  object <- CALL$object.org
  # original model formula
  form <- formula(object)
  p <- length(attr(object$terms, "term.labels")) + 1
  DATA <- data[indices,1:p]
  CALL$object <- update(object, formula = form, data = DATA)
  #CALL$object <- do.call("lm", CALL.lm)
  
  # call restriktor
  do.call("restriktor", CALL)$b.restr
}

