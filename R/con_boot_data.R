# standard bootstrap
con_bootdata_lm <- function(data, indices, ...) {
  l <- list(...)
  object <- l$object
  # original model formula
  form <- formula(object)
  # in case of weights, the data.frame includes the weights
  p <- length(attr(object$terms, "term.labels")) + 1
  # resample data
  DATA <- data[indices,1:p]
  # update lm object with boot data
  l$object <- update(object, formula = form, data = DATA)
  # calll restriktor
  do.call("restriktor", l)$b.restr
}

