# model based bootstrap
con_boot_fixed_lm <- function(data, indices, ...) {
    CALL <- list(...)
    object <- CALL$object.org
    res <- data$res[indices]
    DATA <- cbind(data$fit + res, data[ ,2:(ncol(data)-2), drop = FALSE])
    DATA <- as.data.frame(DATA)
      colnames(DATA) <- colnames(data[ ,1:(ncol(data)-2)])
    form <- formula(object)
    p <- length(attr(object$terms, "term.labels")) + 1
    DATA <- data[indices,1:p]
    CALL$object <- update(object, formula = form, data = DATA)
    
    OUT <- do.call("restriktor", CALL)$b.restr  
    
    OUT
}

