conTest <- function(object, type = "A", ...) {
  
  l <- list(...)
  if (nrow(object$Amat) > object$meq) {
    if (!("test" %in% names(l))) {
      test <- "F"
    } else {
      test <- l$test
    }
    
    if (test == "F" | test == "f") {
      UseMethod("conTestF")
    } else if (test == "LRT" | test == "lrt") {
      UseMethod("conTestLRT")
    } else if (test == "score" | test == "Score") {
      UseMethod("conTestScore")
    } else if (test == "Wald" | test == "wald") {
       if ("lm" %in% class(object)) {
         UseMethod("conTestF")
       } else {
         UseMethod("conTestWald")
       }
#    } else if (test == "Wald2" | test == "wald2") {
#      UseMethod("conTestWald2")  
    } else {
      stop("restriktor ERROR: test ", sQuote(test), " not implemented.")
    }
  } else if (nrow(object$Amat) == object$meq) {
    UseMethod("conTestEq")
  } else {
    stop("Restriktor ERROR: Amat and meq do not match.")
  }
}  
