################################### TO DO ######################################
# - add weights + down weights to RLM.
# - implement E-bar test statistic for lm.
# - add support for no intercept models.
# - add small sample correction for GORIC, see REF. above.
# - add support for mlm and glm.
# - add test-statistic for class of rlm with equality restriktions only.
# - add weights to conTestEq
# - check for right function arguments
###############################################################################
restriktor <- function(model, constraints = NULL, ...) {
  
  # check the class of object
  if (!any(class(model) %in% c("lm", "rlm"))) {
    stop("restriktor only works for lm(), rlm()")
  }
  
  if (class(model)[1] %in% c("lm","mlm")) {
    UseMethod("conLM")
  } 
  else if (class(model)[1] %in% "rlm") {
    UseMethod("conRLM")
  }
  
}
