con_boot_fixed_lm <- function(data, indices, ...) {
    # because we need the original model formula, the original model is parsed via the ...
    l <- list(...)
    # model formula
    form <- l$form
    res <- data$res[indices]
    boot_data <- cbind(data$fit + res, data[ ,2:(ncol(data)-2), drop = FALSE])
    boot_data <- as.data.frame(boot_data)
    colnames(boot_data) <- colnames(data[ ,1:(ncol(data)-2)])
    l$model <- lm(form, data = boot_data)
    out <- do.call("restriktor", l)$b.restr  
    
    OUT <- out
    
      OUT
}



