# model based bootstrap
con_boot_fixed_lm <- function(data, indices, ...) {
    l <- list(...)
    CALL <- as.list(l$CALL)
    CALL[[1]] <- NULL  
    res <- data$res[indices]
    boot_data <- cbind(data$fit + res, data[ ,2:(ncol(data)-2), drop = FALSE])
    boot_data <- as.data.frame(boot_data)
      colnames(boot_data) <- colnames(data[ ,1:(ncol(data)-2)])
    CALL$data <- boot_data
    l$object <- do.call("lm", CALL)
    OUT <- do.call("restriktor", l)$b.restr  
    
    OUT
}


con_boot_fixed_rlm <- function(data, indices, ...) {
  l <- list(...)
  CALL <- as.list(l$CALL)
  CALL[[1]] <- NULL
  res <- data$res[indices]
  boot_data <- cbind(data$fit + res, data[ ,2:(ncol(data)-2), drop = FALSE])
  boot_data <- as.data.frame(boot_data)
    colnames(boot_data) <- colnames(data[ ,1:(ncol(data)-2)])
  CALL$data <- boot_data
  l$object <- do.call("rlm", CALL)
  OUT <- do.call("restriktor", l)$b.restr  
  
  OUT
}

