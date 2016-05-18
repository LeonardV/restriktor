conTest <- function(object, type = "summary", ...) {
  
  if (!inherits(object, "conLM")) {
    stop("object must be of class \"conLM\" or \"conRLM\"")
  }
  
  Amat <- object$constraints
  meq <- object$neq
  
  l <- list(...)
  if (nrow(Amat) > meq) {
    if (!("test" %in% names(l))) {
      test <- "F"
    } else {
      test <- l$test
    }
    if (!("boot" %in% names(l))) {
      boot <- "no"
    } else {
      boot <- l$boot
    }
    # check
    if (boot == "no") {
      # acknowledgement: check for too many inequality constraints is taken from 
      # the ic.infer package. Prevent long runs with late aborts because of too 
      # low memory.
      if (nrow(Amat) - meq - 2 > 2) {
        if (!is.numeric(try(matrix(0, floor((nrow(Amat) - meq -
                                               2)/2), choose(nrow(Amat) - meq, floor((nrow(Amat) - meq -
                                                                                        2)/2))), silent = TRUE)))
          stop(paste("test does not work, too many inequality restriktions, \n",
                     "interim matrix with ", floor((nrow(Amat) - meq)/2) *
                       choose(nrow(Amat) - meq, floor((nrow(Amat) - meq)/2)),
                     " elements cannot be created", sep = ""))
      }
    } 
    # check 
    if (!is.null(l$control)) {
      if ("B" %in% names(l$control)) {
        if (l$control$B < 9999) {
          warning("The number of bootstrap samples for computing the mixing weights is too low and may cause spurious results (by default B = 99999).")
        }
      }
    }
    
    # check
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
  } else if (nrow(Amat) == meq) {
      UseMethod("conTestEq") # Wald, F and score test
  } else {
    stop("Restriktor ERROR: constraints and neq do not match.")
  }
}  
