# standard bootstrap
con_boot_data <- function(data, indices, ...) {
  ldots <- list(...)
  z <- ldots$object
  # original model formula
  form <- formula(z)
  # in case of weights, the data.frame includes the weights
  # resample data
  DATA <- data[indices,1:ncol(data)]
  # update lm object with boot data
  ldots$object <- update(z, formula = form, data = DATA)
  # calll restriktor
  OUT <- do.call("restriktor", ldots)$b_restr
  
  OUT
}

