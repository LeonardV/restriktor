conTest <- function(object, type = "summary", ...) {
  
  if (!inherits(object, "conLM")) {
    stop("object must be of class \"conLM\" or \"conRLM\"")
  }
  
  stopifnot(type %in% c("A","B","C","global","summary"))
  
  l <- list(...)
  Amat <- object$constraints
  meq  <- object$neq
  
  if (nrow(Amat) == 0) {
    stop("no constraints were specified.")
  } else if (nrow(Amat) > meq) {
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
    if (type %in% c("A","B","global")) {
      if (class(object)[1] == "conLM") {
        if (!(test %in% c("F","LRT","score"))) {
          stop("restriktor ERROR: test ", sQuote(test), " unknown. Choose F, LRT or score.")  
        } 
        if (test == "F") {
          UseMethod("conTestF")
        } else if (test == "LRT") {
          UseMethod("conTestLRT")
        } else if (test == "score") {
          UseMethod("conTestScore")
        }
      } else if (class(object)[1] == "conRLM") { 
        if (!(test %in% c("F","Wald","Wald2","score"))) {
          stop("restriktor ERROR: test ", sQuote(test), " unknown. Choose F, Wald, Wald2 or score.")  
        } 
        if (test == "F") {
          UseMethod("conTestF")
        } else if (test == "Wald") {
          UseMethod("conTestWald")
        } else if (test == "Wald2") {
          UseMethod("conTestWald2")
        } else if (test == "score") {
          UseMethod("conTestScore")
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
