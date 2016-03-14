con_boot_fixed_lm <- function(data, indices, ...) {
    l <- list(...)
    form <- l$form
    e <- data$e[indices]
    dat <- cbind(data$fit+e, data[,2:(ncol(data)-2), drop=FALSE])
    dat <- as.data.frame(dat)
    colnames(dat) <- colnames(data[,1:(ncol(data)-2)])
    l$model <- lm(form, data = dat)
    out <- do.call("restriktor", l)$b.constr  
    out
}



