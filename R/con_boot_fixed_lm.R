# acknowledgement: code taken from ic.infer package
con_boot_fixed_lm <- function(data, indices, ...) {

    e <- data$e[indices]
    dat <- cbind(data$fit+e, data[,2:(ncol(data)-2), drop=FALSE])
    dat <- as.data.frame(dat)
      colnames(dat) <- colnames(data[,1:(ncol(data)-2)])
    out <- restriktor(lm(dat), ...)$b.constr
    return(out)
}



