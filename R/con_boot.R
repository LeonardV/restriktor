con_boot_lm <- function(object, B = 999, fixed = FALSE, Amat, 
                        bvec, meq, se = "none", mix_weights = "none", 
                        parallel = parallel, ncpus = ncpus, cl = cl) { 
    
  DATA <- as.data.frame(object$model)
  if (!fixed) { 
    # standard bootstrap
    bootout <- boot(DATA, con_boot_data, R = B, object = object, 
                    constraints = Amat, rhs = bvec, neq = meq, 
                    se = "none", mix_weights = "none", 
                    parallel = parallel, ncpus = ncpus, cl = cl)
  } else { 
    # model based bootstrap
    res <- object$residuals
    fit <- object$fitted.values
    bootout <- boot(data.frame(DATA, fit = fit, res = res), 
                    con_boot_fixed, R = B, object = object, 
                    constraints = Amat, rhs = bvec, neq = meq, 
                    se = "none", mix_weights = "none", 
                    parallel = parallel, ncpus = ncpus, cl = cl)
  }
  
  bootout
}
