con_boot_lm <- function(model, B = 999, fixed = FALSE, ...) { 
    
  if (class(model)[1] != "lm") {
    stop("ERROR: model must be of class lm")
  }
  # original model formula
  #form <- formula(model)
  DATA <- as.data.frame(model$model)
  # weights
  wt <- weights(model)
  if (is.null(wt)) {
    wt <- rep(1, nrow(DATA))
  } # standard bootstrap   
  if (!fixed) {
    bootout <- boot(cbind(wt = wt, DATA), 
                    con_bootdata_lm, 
                    B = B, ...)
                    #form = form, ...)
  } else { # model based bootstrap
    res <- model$residuals
    fit <- model$fitted.values
    bootout <- boot(data.frame(DATA, fit = fit, res = res), 
                    con_boot_fixed_lm, B, ...)
  }
  
  OUT <- bootout
  
    OUT
}


con_boot_rlm <- function(model, B = 999, fixed = FALSE, ...) { 
  
  if (class(model)[1] != "rlm") {
    stop("ERROR: model must be of class rlm")
  }
  # original model formula
  #form <- formula(model)
  DATA <- as.data.frame(model$model)
  # weights
  wt <- weights(model)
  if (is.null(wt)) {
    wt <- rep(1/nrow(DATA), nrow(DATA))
  } # standard bootstrap   
  if (!fixed) {
    bootout <- boot(cbind(wt = wt, DATA), con_bootdata_rlm, 
                    B, ...)
  } else { # model based bootstrap
    res <- model$residuals
    fit <- model$fitted.values
    bootout <- boot(data.frame(DATA, fit = fit, res = res), 
                    con_boot_fixed_rlm, B, ...)
  }
  
  OUT <- bootout
  
  OUT
}

