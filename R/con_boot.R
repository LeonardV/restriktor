con_boot_lm <- function(object, B = 999, fixed = FALSE, ...) { 
    
  if (class(object)[1] != "lm") {
    stop("ERROR: object must be of class lm")
  }
  # original object formula
  #form <- formula(object)
  DATA <- as.data.frame(object$model)
  # weights
  wt <- weights(object)
  if (is.null(wt)) {
    wt <- rep(1, nrow(DATA))
  } # standard bootstrap   
  if (!fixed) {
    bootout <- boot(DATA,  
                    con_bootdata_lm, 
                    R = B, wt = wt, ...)
                    #form = form, ...)
  } else { # model based bootstrap
    res <- object$residuals
    fit <- object$fitted.values
    bootout <- boot(data.frame(DATA, fit = fit, res = res), 
                    con_boot_fixed_lm, R = B, ...)
  }
  
  OUT <- bootout
  
    OUT
}


con_boot_rlm <- function(object, B = 999, fixed = FALSE, ...) { 
  
  if (class(object)[1] != "rlm") {
    stop("ERROR: object must be of class rlm")
  }
  # original model formula
  #form <- formula(object)
  DATA <- as.data.frame(object$model)
  # weights
  wt <- weights(object)
  if (is.null(wt)) {
    wt <- rep(1/nrow(DATA), nrow(DATA))
  } # standard bootstrap   
  if (!fixed) {
    bootout <- boot(cbind(wt = wt, DATA), con_bootdata_rlm, 
                    B, ...)
  } else { # model based bootstrap
    res <- object$residuals
    fit <- object$fitted.values
    bootout <- boot(data.frame(DATA, fit = fit, res = res), 
                    con_boot_fixed_rlm, B, ...)
  }
  
  OUT <- bootout
  
  OUT
}

