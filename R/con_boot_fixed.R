# model based bootstrap
con_boot_fixed <- function(data, indices, ...) {
  ldots <- list(...)
  z <- ldots$object
  res <- data$res[indices]
  DATA <- as.data.frame(cbind(data$fit + res, data[,-1]))
    colnames(DATA) <- colnames(data)
  DATA <- DATA[ ,1:(ncol(data)-2), drop = FALSE]  
  ldots$object <- update(z, formula = formula(z), data = DATA)
  
  OUT <- do.call("restriktor", ldots)$b.restr  
  
  OUT
}
