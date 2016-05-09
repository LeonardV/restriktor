con_boot_lm <- function(model, B = 1000, fixed = FALSE, ...) { 
    
  if (!(c("lm","rlm") %in% class(model))) {
    stop("ERROR: model must be of class lm or rlm.")
  }
  # original model formula
  form <- formula(model)
  DATA <- as.data.frame(model$model)
  # weights
  wt <- weights(model)
  if (is.null(wt)) {
    wt <- rep(1/nrow(DATA), nrow(DATA))
  } # standard bootstrap   
  if (!fixed) {
    bootout <- boot(cbind(wt = wt, DATA), con_bootdata_lm, B,
                    form = form, ...)
  } else { # model based bootstrap
    res <- model$residuals
    fit <- model$fitted.values
    bootout <- boot(data.frame(DATA, fit = fit, res = res), con_boot_fixed_lm,
                    B, form = form, ...)
  }
  
  OUT <- bootout
  
    OUT
}


