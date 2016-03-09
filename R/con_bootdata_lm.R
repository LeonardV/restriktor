con_bootdata_lm <- function(data, indices, ...){
    l <- list(...)
    model <- l$model  
    dat <- data[indices,]
    l$model <- lm(formula(model), data = dat)
    out <- do.call("restriktor", l)  
    out <- out$b.constr
  
    return(out)
}
