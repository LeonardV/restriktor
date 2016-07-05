################################### TO DO ######################################
# - add weights + down weights to RLM.
# - implement E-bar test statistic for lm.
# - add support for no intercept objects.
# - add small sample correction for GORIC
# - add support for mlm and glm.
# - check for right function arguments
###############################################################################
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
