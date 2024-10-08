conTestF        <- function(object, type = type, ...) UseMethod("conTestF")
conTestLRT      <- function(object, type = type, ...) UseMethod("conTestLRT")
conTestScore    <- function(object, type = type, ...) UseMethod("conTestScore")
conTestWald     <- function(object, type = type, ...) UseMethod("conTestWald")
conTest_summary <- function(object, ...) UseMethod("conTest_summary")
conTest_ceq     <- function(object, ...) UseMethod("conTest_ceq")
conTestC        <- function(object, ...) UseMethod("conTestC")
#conTestD        <- function(object, type = type, ...) UseMethod("conTestD")


conTest <- function(object, constraints = NULL, type = "summary", test = "F", 
                    rhs = NULL, neq = 0L, ...) {
  
  
  if (inherits(object, "restriktor")) {
    test <- tolower(test)
    type <- tolower(type)
    stopifnot(type %in% c("global","a","b","c","summary"))
    
    ldots <- list(...)
    Amat <- object$constraints
    meq  <- object$neq
    
    if (all(Amat == 0)) {
      stop("Restriktor ERROR: no constraints specified!")
    } else if (nrow(Amat) > meq) {
      boot <- ldots$boot
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
      if (type %in% c("a","b","global")) {
        if (class(object)[2] %in% "conLM") {
          if (!(test %in% c("f","lrt","score"))) {
            stop("Restriktor ERROR: test ", sQuote(test), " unknown. Choose F, LRT or score.")  
          } 
          if (test == "f") {
            conTestF(object, type = type, ...)
          } else if (test == "lrt") {
            conTestLRT(object, type = type, ...)
          } else if (test == "score") {
            conTestScore(object, type = type, ...)
          } 
        } else if (class(object)[2] %in% "conRLM") { 
          if (!(test %in% c("f","wald","score"))) {
            stop("Restriktor ERROR: test ", sQuote(test), " unknown. Choose F, Wald or score.")  
          } 
          if (test == "f") {
            conTestF(object, type = type, ...)
          } else if (test == "wald") {
            conTestWald(object, type = type, ...)
          } else if (test == "score") {
            conTestScore(object, type = type, ...)
          } 
        } else if (class(object)[2] %in% "conGLM") {
          if (!(test %in% c("f","lrt","score"))) {
            stop("Restriktor ERROR: test ", sQuote(test), " unknown. Choose Wald, LRT or score.")  
          } 
          if (test == "f") {
            conTestF(object, type = type, ...)
          } else if (test == "lrt") {
            conTestLRT(object, type = type, ...)
          } else if (test == "score") {
            conTestScore(object, type = type, ...)
          } 
        } else if (class(object)[2] %in% "conMLM") {
          if (!(test %in% c("lrt"))) {
            stop("Restriktor ERROR: test ", sQuote(test), " unknown. Only LRT available for now.")  
          }
          if (test == "lrt") {
            conTestLRT(object, type = type, ...)
          } 
        }
      } else if (type == "c" && (inherits(object, c("conLM","conRLM","conGLM","conMLM")))) {
        if (inherits(object, "conMLM")) {
          stop("Restriktor ERROR: hypothesis test type c not available (yet).")
        }
        conTestC(object, ...)
      } else if (type == "summary" && (inherits(object, c("conLM","conRLM","conGLM","conMLM")))) {
        conTest_summary(object, test = test, ...)
      } else {
        stop("type ", sQuote(type), " unknown.") 
      }
    } else if (nrow(Amat) == meq) {
      conTest_ceq(object, ...) 
    } else {
      stop("Restriktor ERROR: constraints and neq do not match.")
    }
  } else if (inherits(object, c("lm","mlm","rlm","glm")) && !is.null(constraints)) {
    ldots <- list(...)
    if (is.null(ldots$se)) { 
      ldots$se <- "none"
    }
    
    m.restr <- match(names(ldots), c("se", "B", "mix_weights", 
                                     "parallel", "ncpus", "cl", "seed", "control", 
                                     "verbose", "debug"), 0L)
    
    CALL.restr <- c(list(object = object, constraints = constraints, rhs = rhs, 
                         neq = neq), ldots[m.restr > 0L])
    fit.restr <- do.call("restriktor", CALL.restr) 
    
    m.test <- match(names(ldots), c("neq.alt", "boot", "R", "p.distr", "df", 
                                    "parallel", "ncpus", "cl", "seed", "control", 
                                    "verbose"), 0L)
    CALL.test <- c(list(object = fit.restr, test = test, type = type), 
                   ldots[m.test > 0L])
    do.call("conTest", CALL.test)
  } 
  # else if (inherits(object, "character") && !is.null(constraints)) {
  #   
  #   class(object) <- "conLavaan"
  #   conTestD(object, type = type, ...) 
  #   
  # } 
  else if (is.null(constraints)) {
    if (class(object)[1] == "conTest") { 
      stop("Restriktor ERROR: object is already of class conTest.")
    } 
    else {
      stop("Restriktor ERROR: no constraints found.")
    }
  }
}  
