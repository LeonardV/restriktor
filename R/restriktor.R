restriktor <- function(object, constraints = NULL, ...) {
  
  # check the class of object
  if (!inherits(object, c("lm","rlm","glm"))) {
    stop("restriktor only works for lm(), rlm() and glm().")
  }
  
  if (class(object)[1] %in% c("lm","mlm")) {
    UseMethod("conLM")
  } 
  else if (class(object)[1] %in% "rlm") {
    UseMethod("conRLM")
  }
  else if (class(object)[1] %in% "glm") {
    UseMethod("conGLM")
  }
  
}
