conLM  <- function(object, constraints, ...) UseMethod("conLM")
conRLM <- function(object, constraints, ...) UseMethod("conRLM")
conGLM <- function(object, constraints, ...) UseMethod("conGLM")
conMLM <- function(object, constraints, ...) UseMethod("conMLM")


restriktor <- function(object, constraints = NULL, ...) {
  
  # check the class of object
  if (!inherits(object, c("lm","rlm","glm","mlm"))) {
    stop("Restriktor only works for lm(), mlm(), rlm() and glm().")
  }
  
  arguments <- list(...)
  if (length(arguments)) {
    pnames <- c("se", "B", "rhs", "neq", "mix_weights",  
                "auxilliary", "emControl", "parallel", "ncpus", "cl", "seed", 
                "control", "verbose", "debug", "auto_bound",
                # for rtmvnorm() function
                "lower", "upper", "algorithm",
                "burn.in.samples", "start.values", "thinning")
    
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    if (any(pm == 0L)) { 
      pm.idx <- which(pm == 0L)
      stop("Restriktor Error: ", names(arguments[pm.idx]), " invalid argument(s).")
    }
  }

  if (class(object)[1] %in% c("lm", "aov")) {
    conLM(object, constraints, ...)
  } else if (class(object)[1] == "rlm") {
    conRLM(object, constraints, ...)
  } else if (class(object)[1] == "glm") {
    conGLM(object, constraints, ...)
  } else if (class(object)[1] == "mlm") {
    conMLM(object, constraints, ...)
  } else {
    stop("Restriktor Error: I don't know how to handles objects of class", 
         paste(class(object), sep = ","), call. = FALSE)
  } 
}
