#acknowledgement: code taken from ic.infer package
#slightly adapted by LV
con_boot_lm <- function(model, B = 1000, fixed = FALSE, ...) { #constraints, bvec, meq
    ## check for admissible model
    if (!("lm" %in% class(model)))
        stop("ERROR: model must be of class lm.")

    ## prepare data for bootstrap sampling
    resp <- attr(model$terms, "response")
#    resp.names <- as.character(formula(model))[[2]]
#    y <- as.matrix(model$model[,resp])
#      colnames(y) <- resp.names
    xcol <- which(rowSums(attr(model$terms, "factors")) > 0)
    DATA <- as.data.frame(model$model[, c(resp, xcol)])
    wt <- weights(model)
    if (is.null(wt))
        wt <- rep(1/nrow(DATA), nrow(DATA))
    if (!fixed)
      bootout <- boot(cbind(DATA), con_bootdata_lm, B,
                      model = model, ...)
    else {
        e <- model$residuals
        fit <- model$fitted.values
        bootout <- boot(data.frame(DATA, fit = fit, e = e), con_boot_fixed_lm,
                        B, model = model, ...)
    }
    
    bootout
}


