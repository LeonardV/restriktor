con_boot_lm <- function(object, B = 999, fixed = FALSE, ...) { 
    
  if (!any(class(object) %in% c("lm", "rlm"))) {
    stop("ERROR: object must be of class lm or rlm")
  }
  DATA <- as.data.frame(object$model)
  if (!fixed) { # standard bootstrap
    bootout <- boot(DATA,  
                    con_bootdata_lm, 
                    R = B, ...)
  } else { # model based bootstrap
    res <- object$residuals
    fit <- object$fitted.values
    bootout <- boot(data.frame(DATA, fit = fit, res = res), 
                    con_boot_fixed_lm, R = B, ...)
  }
  
  OUT <- bootout
  
    OUT
}
