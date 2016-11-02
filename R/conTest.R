conTest <- function(object, type = "summary", test = "F", ...) {
  
  # some checks
  if (!inherits(object, "conLM")) {
    stop("object must be of class \"conLM\" or \"conRLM\"")
  }
  if (type != "global" && type != "summary") {
    type <- toupper(type)
  }
  if ( !(type %in% c("A","B","C","global","summary")) ) {
    stop("restriktor ERROR: hypothesis type ", sQuote(type), " unknown. \nPlease, choose from \"A\", \"B\", \"C\", \"global\", or \"summary\".")
  }
  if (is.null(object$wt)) {
    stop("restriktor ERROR: no chi-square-bar weights computed. Set Wt = TRUE in the restriktor() function.")
  } 
  
  l <- list(...)
  Amat <- object$constraints
  meq  <- object$neq
  
  if (all(Amat == 0)) {
    stop("Restriktor ERROR: no constraints specified!")
  } else if (nrow(Amat) > meq) {
    boot <- l$boot
    if (!is.null(boot) && boot == "no") {
      # Prevents long runs with late aborts because of too low memory.
      # acknowledgement: taken from the ic.infer package (Gromping, 2010)
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
      UseMethod("conTestEq") # classical Wald, F and score test
  } else {
    stop("Restriktor ERROR: constraints and neq do not match.")
  }
}  
