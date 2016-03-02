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
      UseMethod("conTestscore")  
    } else {
      stop("restriktor ERROR: test ", sQuote(test), " not (yet) implemented. Choose \"F\", \"score\", or \"LRT\"")
    }
  } else if (nrow(object$Amat) == object$meq) {
    UseMethod("conTestEq")    
  }
}  
