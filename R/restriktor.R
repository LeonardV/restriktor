restriktor <- function(object, constraints = NULL, ...) {
  
  # check the class of object
  if (!inherits(object, c("lm","rlm","glm","mlm"))) {
    stop("restriktor only works for lm(), mlm(), rlm() and glm().")
  }
  
  arguments <- list(...)
  if (length(arguments)) {
    pnames <- c("se", "B", "rhs", "neq", "mix.weights", "mix.bootstrap", "parallel", 
                "ncpus", "cl", "seed", "control", "verbose", "debug")
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    if (any(pm == 0L)) { 
      pm.idx <- which(pm == 0L)
      stop("Restriktor ERROR: ", names(arguments[pm.idx]), " invalid argument(s).")
    }
  }
  
  if (class(object)[1] == "lm") {
    UseMethod("conLM")
  } 
  else if (class(object)[1] == "rlm") {
    UseMethod("conRLM")
  }
  else if (class(object)[1] == "glm") {
    UseMethod("conGLM")
  }
  else if (class(object)[1] == "mlm") {
    UseMethod("conMLM")
  }
  
}
