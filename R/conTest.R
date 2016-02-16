conTest <- function(object, test = "F", type = "A", boot = "none",
                    meq.alt = 0L, control = NULL, 
                    tol = sqrt(.Machine$double.eps), ...) {
  
  out <- 
  if ("lm" %in% class(object)) {
    if (test == "F") {
      conTestF_lm(object, type = type, boot = boot, meq.alt = meq.alt,
                  control = control, tol = tol)
    } else if (test == "LRT") {
      conTestLRT_lm(object, type = type, boot = boot, meq.alt = meq.alt,
                    control = control, tol = tol) 
    }
    else if ("lm" %in% class(object) && !(test %in% c("F","LRT"))) {
      stop("test ", sQuote(test), " not (yet) implemented.")      
    }
  } else if (!("lm" %in% class(object))) {
      stop("object class ", sQuote(class(object))," not yet implemented.")
  }
  
  out
}