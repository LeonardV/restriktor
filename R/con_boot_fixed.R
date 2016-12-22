# model based bootstrap
con_boot_fixed <- function(data, indices, ...) {
    ldots <- list(...)
    z <- ldots$object
    res <- data$res[indices]
    DATA <- cbind(data$fit + res, data[ ,2:(ncol(data)-2), drop = FALSE])
    DATA <- as.data.frame(DATA)
      colnames(DATA) <- colnames(data[ ,1:(ncol(data)-2)])
    form <- formula(z)
    p <- length(attr(z$terms, "term.labels")) + 1
    DATA <- data[indices,1:p]
    ldots$object <- update(z, formula = form, data = DATA)
    
    OUT <- do.call("restriktor", ldots)$b.restr  
    
    OUT
}
