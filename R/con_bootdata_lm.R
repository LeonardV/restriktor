con_bootdata_lm <- function(data, indices, ...){
    ## df.error chosen arbitrary, since it is irrelevant for b.constr
#    data <- as.data.frame(sapply(data, as.numeric))
    dat <- data[indices,]
    out <- restriktor(lm(dat), ...)$b.constr
    return(out)
#    bx <- restriktor(cov.wt(dat[,-1], wt = dat$wt)$cov, df.error = 10, ...)$b.constr
    ## intercept
#    bi <- mean(dat[,2]) - sum(colMeans(dat[,3:ncol(dat)])*bx)

#    return(c(bi,bx))
}
