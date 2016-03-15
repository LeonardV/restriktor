conTest <- function(object, type = "A", ...) {
  
  l <- list(...)
  if (nrow(object$Amat) > object$meq) {
    if (!("test" %in% names(l))) {
      test <- "F"
    } else {
      test <- l$test
    } 
    if (test == "F") {
      UseMethod("conTestF")
    } else if (test == "LRT") {
      UseMethod("conTestLRT")  
    } else if (test == "score") {
      UseMethod("conTestScore")  
    } else if (test == "wald") {
      UseMethod("conTestWald")  
    } else if (test == "wald2") {
      UseMethod("conTestWald2")  
    } else {
      stop("restriktor ERROR: test ", sQuote(test), " not (yet) implemented.")
    }
  } else if (nrow(object$Amat) == object$meq) {
    UseMethod("conTestEq")    
  } else {
    stop("Restriktor ERROR: Amat and meq do not match.")
  }
}  
