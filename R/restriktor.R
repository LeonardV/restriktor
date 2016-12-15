restriktor <- function(object, constraints = NULL, ...) {
  
  # check the class of object
  if (!any(class(object) %in% c("lm", "rlm"))) {
    stop("restriktor only works for lm(), rlm()")
  }
  
  if (class(object)[1] %in% c("lm","mlm")) {
    UseMethod("conLM")
  } 
  else if (class(object)[1] %in% "rlm") {
    UseMethod("conRLM")
  }
  
}
