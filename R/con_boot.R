con_boot_lm <- function(B = 999, fixed = FALSE, ...) { 
    
  ldots <- list(...)
  z <- ldots$object
  DATA <- as.data.frame(z$model)
  if (!fixed) { 
    # standard bootstrap
    bootout <- boot(DATA,  
                    con_boot_data, 
                    R = B, ...)
  } else { 
    # model based bootstrap
    res <- z$residuals
    fit <- z$fitted.values
    bootout <- boot(data.frame(DATA, fit = fit, res = res), 
                    con_boot_fixed, R = B, ...)
  }
  
  bootout
  
}
