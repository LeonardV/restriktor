conTest <- function(object, type = "A", ...) {
  
  if (!inherits(object, "conLM")) {
    stop("object must be of class \"conLM\" or \"conRLM\"")
  }
  
  l <- list(...)
  if (nrow(object$Amat) > object$meq) {
    if (!("test" %in% names(l))) {
      test <- "F"
    } else {
      test <- l$test
    }
    if (class(object)[1] == "conLM" && !(test %in% c("F","Wald","LRT","score"))) {
      stop("restriktor ERROR: test ", sQuote(test), " unknown. Choose F, LRT or score.")  
    } else if (class(object)[1] == "conRLM" && !(test %in% c("F","Wald","score"))) {
      stop("restriktor ERROR: test ", sQuote(test), " unknown. Choose F, Wald or score.")  
    } 
    if (type %in% c("A","B","global")) {
      if (test == "F") {
        UseMethod("conTestF")
      } else if (test == "LRT") {
        UseMethod("conTestLRT")
      } else if (test == "score") {
        UseMethod("conTestScore")
      } else if (test == "Wald") {
         if (class(object)[1] == "conLM") {
           UseMethod("conTestF")
         } else if (class(object)[1] == "conRLM") {
           UseMethod("conTestWald")
         }
      } 
    } else if (type == "C") {
      UseMethod("conTestC")
    } else if (type == "summary") {
      UseMethod("summary.conTest")     
    } else {
      stop("type ", sQuote, " unknown.")
    }
  } else if (nrow(object$Amat) == object$meq) {
    UseMethod("conTestEq")
  } else {
    stop("Restriktor ERROR: Amat and meq do not match.")
  }
}  
