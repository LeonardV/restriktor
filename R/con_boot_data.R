# standard bootstrap
con_boot_data <- function(data, indices, ...) {
  ldots <- list(...)
  z <- ldots$object
  # in case of weights, the data.frame includes the weights
  # resample data
  DATA <- data[indices, ]
  # update lm object with boot data
  ldots$object <- update(z, formula = formula(z), data = DATA)
  # calll restriktor
  OUT <- do.call("restriktor", ldots)$b.restr
  
  OUT
}

