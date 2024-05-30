conTestD <- function(model = NULL, data = NULL, constraints = NULL, 
                     type = c("A","B"), R = 1000L, bootstrap.type = "bollen.stine", 
                     return.test = TRUE, neq.alt = 0, 
                     double.bootstrap = "standard", double.bootstrap.R = 249, 
                     double.bootstrap.alpha = 0.05, 
                     parallel = c("no", "multicore", "snow"), 
                     ncpus = 1L, cl = NULL, verbose = FALSE, ...) {

  # fit unrestricted model
  fit.h2 <- sem(model, ..., data = data, test = "standard") #se = "none"
  
  stopifnot(!is.null(constraints))
  
  ## fit null-model
  # add constraints to parameter table
  CON <- attr(lavParseModelString(constraints), "constraints")
  parTable <- fit.h2@ParTable
  for (con in seq_len(length(CON))) {
    parTable <- lav_partable_add(parTable, CON[[con]])
  }
  
  # replace <, > with ==
  user.equal <- parTable
  for (con in seq_len(length(CON))) {
    if (CON[[con]]$op %in% c("<", ">")) {
      this.lhs <- CON[[con]]$lhs
      this.op  <- CON[[con]]$op
      this.rhs <- CON[[con]]$rhs
      
      # find this line in user.equal@ParTable
      idx <- which(user.equal$lhs == this.lhs,
                   user.equal$op  == this.op,
                   user.equal$rhs == this.rhs)
      if (length(idx) == 0L) { # not found, give warning?
        stop("lavaan ERROR: no inequality constraints (<, >) found.")
      }
      
      # change operator to ==
      user.equal$op[idx] <- "=="
    }
  }
    
  fit.h0 <- sem(user.equal, ..., data = data, test = "standard") # se = "none"
  
  if ("A" %in% type) {
    bootA <- bootstrapD(h0 = fit.h0, h1 = fit.h2, 
                        constraints = constraints, 
                        R = R, bootstrap.type = bootstrap.type, type = "A", 
                        verbose = verbose, return.D = return.test, 
                        double.bootstrap = double.bootstrap, 
                        double.bootstrap.R = double.bootstrap.R, 
                        double.bootstrap.alpha = double.bootstrap.alpha, 
                        parallel = parallel, ncpus = ncpus, cl = cl)
  }
  
  if ("B" %in% type) {
    bootB <- bootstrapD(h0 = fit.h0, h1 = fit.h2, 
                        constraints = constraints, 
                        R = R, bootstrap.type = bootstrap.type, type = "B", 
                        verbose = verbose, return.D = return.test, 
                        double.bootstrap = double.bootstrap, 
                        double.bootstrap.R = double.bootstrap.R, 
                        double.bootstrap.alpha = double.bootstrap.alpha, 
                        parallel = parallel, ncpus = ncpus, cl = cl)
  }
  
  
  output <- list(fit.h0 = fit.h0, fit.h2 = fit.h2,
                 double.bootstrap = double.bootstrap, 
                 double.bootstrap.alpha = double.bootstrap.alpha, 
                 return.test = return.test, 
                 type = bootstrap.type)
  
  if ("A" %in% type) {
    output$bootA <- bootA
  }
  if ("B" %in% type) {
    output$bootB <- bootB
  }
  
  class(output) <- "conTestLavaan"
  return(output)
}