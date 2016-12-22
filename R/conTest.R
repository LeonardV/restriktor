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
        if (!inherits(object, "conRLM")) {
          if (!(test %in% c("f","lrt","score"))) {
            stop("restriktor ERROR: test ", sQuote(test), " unknown. Choose F, LRT or score.")  
          } 
          if (test == "f") {
            UseMethod("conTestF")
          } else if (test == "lrt") {
            UseMethod("conTestLRT")
          } else if (test == "score") {
            UseMethod("conTestScore")
          } 
        } else if (inherits(object, "conRLM")) { 
          if (!(test %in% c("f","wald","wald2","score"))) {
            stop("restriktor ERROR: test ", sQuote(test), " unknown. Choose F, Wald, Wald2 or score.")  
          } 
          if (test == "f") {
            UseMethod("conTestF")
          } else if (test == "wald") {
            UseMethod("conTestWald")
          } else if (test == "wald2") {
            UseMethod("conTestWald2")
          } else if (test == "score") {
            UseMethod("conTestScore")
          } 
        }
      } else if (type == "c") {
        UseMethod("conTestC")
      } else if (type == "summary") {
        UseMethod("conTest_summary")     
      } else {
        stop("type ", sQuote, " unknown.")
      }
    } else if (nrow(Amat) == meq) {
        UseMethod("conTest_ceq") # classical Wald, F and score test with equality constraints only
    } else {
      stop("Restriktor ERROR: constraints and neq do not match.")
    }
  } else if (inherits(object, c("lm","rlm")) && !is.null(constraints)) {
    ldots <- list(...)
    if (is.null(ldots$se)) { 
      ldots$se <- "none"
    }
    
    m_restr <- match(names(ldots), c("se", "B", "Wt", "bootWt", "bootWt.R", 
                                     "parallel", "ncpus", "cl", "seed", "control", 
                                     "verbose", "debug"), 0L)
    
    CALL_restr <- c(list(object = object, constraints = constraints, rhs = rhs, 
                         neq = neq), ldots[m_restr > 0L])
    fit_restr <- do.call("restriktor", CALL_restr) 
    
    m_test <- match(names(ldots), c("neq.alt", "boot", "R", "p.distr", "df", 
                                    "parallel", "ncpus", "cl", "seed", "control", 
                                    "verbose"), 0L)
    CALL_test <- c(list(object = fit_restr, test = test, type = type), 
                   ldots[m_test > 0L])
    do.call("conTest", CALL_test)
  } else if (is.null(constraints)) {
    stop("Restriktor ERROR: no constraints specified!")
  }
}  
