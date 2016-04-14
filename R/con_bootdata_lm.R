con_bootdata_lm <- function(data, indices, ...){
    l <- list(...)
    form <- l$form
    dat <- data[indices,]
    l$model <- lm(form, data = dat)
    out <- do.call("restriktor", l)$b.constr
    
    OUT <- out
    
      OUT
}
