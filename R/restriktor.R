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
    pnames <- c("se", "B", "rhs", "neq", "mix.weights", "mix.bootstrap", "missing",
                "auxilliary", "emControl", "parallel", "ncpus", "cl", "seed", "control",  
                "verbose", "debug", "auto_bound")
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    if (any(pm == 0L)) { 
      pm.idx <- which(pm == 0L)
      stop("Restriktor ERROR: ", names(arguments[pm.idx]), " invalid argument(s).")
    }
  }
  
  missing <- arguments$missing
  
  if (is.null(missing)) { 
    missing <- "none"
  } else if (missing %in% c("em", "EM", "two.stage", "twostage")) {
    missing <- "two.stage" 
  } else if (!(missing %in% c("none", "two.stage"))) {
      stop("Restriktor ERROR: missing method ", sQuote(missing), " unknown.", call. = FALSE)
  }  
  
  if (missing == "two.stage") {
    # EM algorithm for incomplete multivariate normal data 
    data_imp <- two_stage(object     = object, 
                          emControl  = arguments$emControl, 
                          auxilliary = arguments$auxilliary)
    
    object <- update(object, data = data_imp)
  }
  
  arguments$missing <- NULL
  
  if (class(object)[1] == "lm") {
    conLM(object, constraints, ...)
  } 
  else if (class(object)[1] == "rlm") {
    conRLM(object, constraints, ...)
  }
  else if (class(object)[1] == "glm") {
    conGLM(object, constraints, ...)
  }
  else if (class(object)[1] == "mlm") {
    conMLM(object, constraints, ...)
  }  

}
