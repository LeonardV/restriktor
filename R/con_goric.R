goric <- function(object, ...) {
  UseMethod("goric")
}


goric.default <- function(object, ..., 
                          comparison = c("unconstrained", "complement", "none"), 
                          VCOV = NULL, sample.nobs = NULL,
                          type = "goric", bound = NULL, debug = FALSE) {

  mc <- match.call()
  CALL <- as.list(mc[-1])
  
  # some checks
  comparison <- tolower(comparison)
  comparison <- match.arg(comparison)
  stopifnot(comparison %in% c("unconstrained", "complement", "none"))
  
  type <- tolower(type)
  stopifnot(type %in% c("goric", "goricc", "gorica", "goricac"))
  
  ldots <- list(...)
  arguments <- c("B", "mix.weights", "mix.bootstrap", "parallel", "ncpus", 
                 "cl", "seed", "control", "verbose", "debug")
  m.restr <- pmatch(names(ldots), arguments, 0L)
  if (length(m.restr) == 0L) {
    ldots2 <- ldots
  } else {
    ldots2 <- ldots[m.restr == 0L]
  }  

  if (inherits(object, "matrix")) {
    if (dim(object)[2] > 1L) {
      stop("restriktor ERROR: object must be a vector.")
    }
    object <- as.vector(object)
  }
  
  # what is the constraint input type
  objectList  <- c(list(object), ldots2)  
  if (is.name(CALL$object)) {
    names(objectList)[1] <- as.character(CALL$object)
  } else {
    names(objectList)[1] <- as.character("object")
  }
    
  isRestr     <- sapply(objectList, function(x) inherits(x, "restriktor"))
  isCharacter <- sapply(ldots2,     function(x) inherits(x, "character"))
  isList      <- sapply(ldots2,     function(x) inherits(x, "list"))

  ## catching user errors
  # only ldots2
  isRestr2 <- sapply(ldots2, function(x) inherits(x, "restriktor"))
  col.check <- rbind(isRestr2, isCharacter, isList)
  if (length(col.check) > 0) {
    unknown.idx <- colSums(col.check) == 0
  } 
  
  if (length(ldots2) > 0 && any(unknown.idx)) {  
    stop("restriktor ERROR: argument ", sQuote(names(ldots2)[unknown.idx]), " unknown.", call. = FALSE)
  }
    
  # create output list
  ans <- list()
  
  ## how to deal with constraints
  # if all objects are of class restriktor
  if (all(isRestr)) {   
    conList   <- objectList
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     goric       = type,
                                                     sample.nobs = sample.nobs))
    ans$model.org <- object$model.org
    sample.nobs   <-  nrow(model.frame(object$model.org))
    # unrestricted VCOV
    VCOV <- vcov(ans$model.org)
    constraints <- conList
  # if the constraints syntax is of class character, e.g., x1 < x2; x2 < x3
  } else if (inherits(object, "lm") && all(isCharacter)) {
    # check if constraints exists
    if (length(ldots2) == 0) {
      stop("restriktor ERROR: no constraints found!", call. = FALSE)
    }
    constraints <- ldots2[isCharacter]
    # standard errors are not needed
    ldots3 <- ldots[m.restr > 0L]
    if (is.null(ldots$se)) { 
      ldots3$se <- "none"
    }
    # fit restriktor object for each hypothesis
    conList <- list()
    for (c in 1:length(constraints)) {
      CALL.restr <- c(list(object = object,
                           constraints = constraints[[c]]), ldots3)
      conList[[c]] <- do.call("restriktor", CALL.restr) 
    }
    # compute summary for each restriktor object. 
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     goric       = type,
                                                     sample.nobs = sample.nobs))
    # add unrestricted object to output
    ans$model.org <- object
    # unrestricted VCOV
    VCOV <- vcov(ans$model.org)
    sample.nobs <- nrow(model.frame(object))
    idx <- length(conList) 
    objectnames <- vector("character", idx)
    CALL$object <- NULL
  # if the constraints are a list with constraints, rhs and neq for each hypothesis   
  } else if (inherits(object, "lm") && all(isList)) {
      # create lists
      # constraints <- list(); rhs <- list(); neq <- list()
      # extract constraints, rhs and neq
      #constraints <- list()
      constraints <- lapply(ldots2, FUN = function(x) { x$constraints } )
      constraints.check <- sapply(constraints, FUN = function(x) { is.null(x) } )
      if (any(constraints.check)) {
        stop("restriktor ERROR: the constraints must be specified as a list. \nFor example, h1 <- list(constraints = 'x1 > 0')", call. = FALSE)
      }
      rhs <- lapply(ldots2, FUN = function(x) {x$rhs} )
      neq <- lapply(ldots2, FUN = function(x) {x$neq} )
      # standard errors are not needed
      ldots3 <- ldots[m.restr > 0L]
      if (is.null(ldots$se)) { 
        ldots3$se <- "none"
      }
      # fit restriktor object for each hypothesis
      conList <- list()
      for (c in 1:length(constraints)) {
        CALL.restr <- c(list(object = object,
                             constraints = constraints[[c]],
                             rhs = rhs[[c]],
                             neq = neq[[c]]), ldots3)
        
        conList[[c]] <- do.call("restriktor", CALL.restr) 
      }
      # compute symmary for each restriktor object. Here is the goric value 
      # computed. Note, not the gorica value
      isSummary <- lapply(conList, function(x) summary(x, 
                                                       goric       = type,
                                                       sample.nobs = sample.nobs))
      # add unrestricted object to output
      ans$model.org <- object
      # unrestricted VCOV
      VCOV <- vcov(ans$model.org) 
      sample.nobs <- nrow(model.frame(object))
      idx <- length(conList) 
      objectnames <- vector("character", idx)
      CALL$object <- NULL
    } else if (inherits(object, "numeric")) {
      if (type %in% c("goric", "goricc")) {
        stop("restriktor ERROR: object of class numeric is only supported for type = 'gorica(c)'.")
      }
      if (is.null(VCOV)) {
        stop("restriktor ERROR: the argument VCOV is not found.")
      }
      if (is.null(sample.nobs) && type %in% c("goricc", "goricac")) {
        stop("restriktor ERROR: the argument sample.nobs is not found.")
      }
      
      isCharacter <- sapply(ldots2, function(x) inherits(x, "character"))    
      isList      <- sapply(ldots2, function(x) inherits(x, "list"))
      
      ldots3 <- ldots[m.restr > 0L]
      if (all(isList)) { 
        # create lists
        # constraints <- list(); rhs <- list(); neq <- list()
        # extract constraints, rhs and neq
        #constraints <- list()
        constraints <- lapply(ldots2, FUN = function(x) { x$constraints } )
        constraints.check <- sapply(constraints, FUN = function(x) { is.null(x)} )
        if (any(constraints.check) || length(constraints) == 0) {
          stop("restriktor ERROR: no constraints found! The constraints must be \nspecified as a list. For example, h1 <- list(constraints = 'x1 > 0')")
        }
        rhs <- lapply(ldots2, FUN = function(x) {x$rhs} )
        neq <- lapply(ldots2, FUN = function(x) {x$neq} )
        
        # fit restriktor object for each hypothesis
        conList <- list()
        for (c in 1:length(constraints)) {
          CALL.restr <- c(list(object      = c(object),
                               constraints = constraints[[c]],
                               rhs         = rhs[[c]],
                               neq         = neq[[c]],
                               VCOV        = VCOV),
                          ldots3)
          conList[[c]] <- do.call("con_gorica_est", CALL.restr) 
        }
        isSummary <- lapply(conList, function(x) summary(x, 
                                                         type        = type,
                                                         sample.nobs = sample.nobs)) 
      } else if (all(isCharacter)) {
        # if constraints is character 
        constraints <- ldots2[isCharacter]
        # fit restriktor object for each hypothesis
        conList <- list()
        for (c in 1:length(constraints)) {
          CALL.restr <- c(list(object      = c(object),
                               constraints = constraints[[c]],
                               VCOV        = VCOV),
                          ldots3)
          conList[[c]] <- do.call("con_gorica_est", CALL.restr)  
        }
        isSummary <- lapply(conList, function(x) summary(x, 
                                                         type        = type,
                                                         sample.nobs = sample.nobs)) 
      }
      CALL$object <- NULL; CALL$comparison <- NULL; CALL$type <- NULL
      CALL$VCOV <- NULL; CALL$bound <- NULL; 
      CALL$debug <- NULL; CALL$B <- NULL; CALL$mig.weights <- NULL; 
      CALL$mix.bootstrap <- NULL; CALL$parallel <- NULL; CALL$ncpus <- NULL; 
      CALL$cl <- NULL; CALL$seed <- NULL; CALL$control <- NULL; 
      CALL$verbose <- NULL; 
      
      if (!exists("conList")) {
        stop("restriktor ERROR: no constraints found!", call. = FALSE)
      }
    } else {
      stop("restriktor ERROR: I don't know how to handle an object of class ", paste0(class(object)[1]))
    }

  
  ## add objectnames if not available
  # constraints must be a list
  if (!is.list(constraints)) {
    constraints <- list(constraints)
  }
  
  if (any(is.null(names(constraints))) || all(names(constraints) == "")) {  
    objectnames <- paste0("H", 1:length(constraints))
  } else {
    objectnames <- names(constraints) 
  }

  
  if (comparison == "complement" && length(conList) > 1L) {
    comparison <- "unconstrained"
    warning("restriktor WARNING: if comparison = 'complement', only one order-restricted hypothesis\n",
            " is allowed (for now). Therefore, comparison is set to 'unconstrained'.",
            call. = FALSE)
  } 

  

# compute complement ------------------------------------------------------
  df.c <- NULL
  if (comparison == "complement") {
      # unrestricted estimates
      if (inherits(object, "numeric")) {
        b.unrestr <- object
      } else {
        b.unrestr <- coef(ans$model.org)
      }
      # restricted estimates
      b.restr <- conList[[1]]$b.restr
      # number of parameters
      p <- length(b.unrestr)
      # level probabilities
      wt.bar <- conList[[1]]$wt.bar
      # constraints matrix
      Amat <- conList[[1]]$constraints
      # number of equalities
      meq <- conList[[1]]$neq
      # rhs
      bvec <- conList[[1]]$rhs
    # extract equalities and inequalities
    if (meq > 0) {
      Amat.ceq <- Amat[1:meq, , drop = FALSE]
      bvec.ceq <- bvec[1:meq]
      Amat.ciq <- Amat[-c(1:meq), , drop = FALSE]
      bvec.ciq <- bvec[-c(1:meq)]
    } else {
      Amat.ceq <- matrix( , nrow = 0, ncol = ncol(Amat))
      bvec.ceq <- rep(0, 0)
      Amat.ciq <- Amat[ , , drop = FALSE]
      bvec.ciq <- bvec
    }
    
    if (!is.null(bound) && meq == 0L) {
      warning("restriktor WARNING: bounds are only available for equality restrictions \n",
              " and are therefore ignored.", call. = FALSE)
      bound <- NULL
    } 
    
    # for now
    #bound <- NULL
    

# about equality bounds ---------------------------------------------------
    if (!is.null(bound)) {

      betasc <- NA
      
      # check if any equality constraints
      if (nrow(Amat.ceq) != length(bound)) {
        stop("Restriktor ERROR: the number of bounds must be equal to the number of equality restrictions.")
      }
      
      # upper-bound = rhs + bound
      ub <- bvec.ceq + bound
      # lower-bound = rhs - bound
      lb <- bvec.ceq - bound
      
      # correct user error
      if ( length(ub) == 1L || length(lb) == 1L && meq >= 1L ) {
        ub <- rep(ub, meq)
        lb <- rep(lb, meq)
      }
      
      nr.meq <- 1:meq
      
      ## If no constraints are violated, then the complement is on the boundary
      ## of the bounds.
      
      # check if unconstrained mle are violated
      if (nrow(Amat.ciq) > 0) {
        check.ciq <- all(Amat.ciq %*% b.unrestr - bvec.ciq >= 0)
      } else {
        check.ciq <- TRUE
      }
      # check if unconstrained mle lay between the bounds
      check.ub <- all(Amat.ceq %*% b.unrestr - bvec[1:meq] - ub <= 0)
      check.lb <- all(Amat.ceq %*% b.unrestr - bvec[1:meq] - lb >= 0)
     
      # check if unrestricted mle lay in boundary area
      if (check.ciq & check.ub & check.lb) {
      # log-likelihood model
        if (type == "goric") {
          llm <- logLik(ans$model.org)
        } else if (type == "gorica") {
          llm <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
        }

        ## determine log-likelihood complement
        # upper-bound                                                           # enkel voor unconstrained mle die niet tussen de bounds liggen?
        llc.ub <- list()
        for (l in nr.meq) {
          Amat.ub <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
          bvec.ub <- c(ub[l], bvec.ciq)
          Hc.ub <- restriktor(ans$model.org, constraints = Amat.ub,
                              neq = 1, rhs = bvec.ub,
                              mix.weights = "none", se = "none")
          if (type == "goric") {  
           llc.ub[[l]] <- logLik(Hc.ub)
          } else if (type == "gorica") {
           llc.ub[[l]] <- dmvnorm(c(b.unrestr - Hc.ub$b.restr), 
                                  sigma = VCOV, log = TRUE)
          }
        }
       # lower-bound                                                            # enkel voor unconstrained mle die niet tussen de bounds liggen?
        llc.lb <- list()
        for (l in nr.meq) {
          Amat.lb <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
          bvec.lb <- c(lb[l], bvec.ciq)
          Hc.lb <- restriktor(ans$model.org, constraints = Amat.lb,
                              neq = 1, rhs = bvec.lb,
                              mix.weights = "none", se = "none")
          if (type == "goric") {
            llc.lb[[l]] <- logLik(Hc.lb)
          } else if (type == "gorica") {
            llc.lb[[l]] <- dmvnorm(c(b.unrestr - Hc.lb$b.restr),
                                   sigma = VCOV, log = TRUE)
          }
        }
       
       # if any inequalty constraints, then all equalities are released and the
       # inequality constraints are treated as equalities, separately.
        if (nrow(Amat.ciq) > 0L) {
          llc.s <- list()
          nr.ciq <- 1:nrow(Amat.ciq)
          for (l in nr.ciq) {
            Hc.s <- restriktor(ans$model.org,
                               constraints = Amat.ciq[l, , drop = FALSE],
                               neq = 1, rhs = bvec.ciq[l],
                               mix.weights = "none", se = "none")
            b.s <- coef(Hc.s)
            b.s[abs(b.s) < sqrt(.Machine$double.eps)] <- 0L
            
            Amatx <- rbind(Amat.ciq, Amat.ceq, -Amat.ceq)                       # correct?
            bvecx <- c(bvec.ciq, lb, lb)                                        # correct?
            if (all(Amatx %*% b.s - bvecx >= 0)) {                              # correct check?
              # log-likelihood model
              if (type == "goric") {
                llc.s[[l]] <- logLik(Hc.s)
              } else if (type == "gorica") {
                llc.s[[l]] <- dmvnorm(c(b.unrestr - Hc.s$b.restr),
                                      sigma = VCOV, log = TRUE)
              }
            }
          }
        } else {
          llc.s <- NULL
        }
        llc <- unlist(c(llc.ub, llc.lb, llc.s))
        llc <- max(llc)

        if (debug) {
          cat(llc)
        }
      } else {
        # determine the log-likelihood_m in case the unconstrained mle
        # lay outside the range restrictions.
        # 2^q2 combinations.
        if (type == "goric") {
          llc <- logLik(ans$model.org)
        } else if (type == "gorica") {
          llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
        }
       
        # if all bounds are zero, then llm = logLik(object)
        if (!all(bound == 0L)) {
        
        
        # upper-bound                                                           # enkel voor unconstrained mle die niet tussen de bounds liggen?
        llm.ub <- list()
        for (l in nr.meq) {
          Amat.ub <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
          bvec.ub <- c(ub[l], bvec.ciq)
          Hm.ub <- restriktor(ans$model.org, constraints = Amat.ub,
                              neq = 1, rhs = bvec.ub,
                              mix.weights = "none", se = "none")
          if (type == "goric") {  
            llm.ub[[l]] <- logLik(Hm.ub)
          } else if (type == "gorica") {
            llm.ub[[l]] <- dmvnorm(c(b.unrestr - Hm.ub$b.restr), 
                                   sigma = VCOV, log = TRUE)
          }
        }
        # lower-bound                                                            # enkel voor unconstrained mle die niet tussen de bounds liggen?
        llm.lb <- list()
        for (l in nr.meq) {
          Amat.lb <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
          bvec.lb <- c(lb[l], bvec.ciq)
          Hm.lb <- restriktor(ans$model.org, constraints = Amat.lb,
                              neq = 1, rhs = bvec.lb,
                              mix.weights = "none", se = "none")
          if (type == "goric") {
            llm.lb[[l]] <- logLik(Hm.lb)
          } else if (type == "gorica") {
            llm.lb[[l]] <- dmvnorm(c(b.unrestr - Hm.lb$b.restr),
                                   sigma = VCOV, log = TRUE)
          }
        }
        
        #### 
        # wat doen we met de inqualities ??
        ####
        
        llm.s <- NULL
        llm <- unlist(c(llm.ub, llm.lb, llm.s))
        if (debug) {
          cat("log-likelihood_m =", llm, "\n")
        }
 
        llm <- max(llm)
        } else {
         if (type == "goric") {
           llm <- logLik(conList[[1]])
         } else if (type == "gorica") {
           llm <- dmvnorm(c(b.unrestr - b.restr), sigma = VCOV, log = TRUE)
         }
        }
        llm
      
      }
#################################### no bounds ################################    
    } else if (is.null(bound)) {
      # check if any equality constraint is violated
      check.ceq <- !(all(Amat.ceq %*% c(b.unrestr) - bvec.ceq == 0))
      if (nrow(Amat) > meq) {
        # check if any inequality constraint is violated
        check.ciq <- !(all(Amat.ciq %*% c(b.unrestr) - bvec.ciq >= 0))
      } else {
        check.ciq <- FALSE
      }
      ## compute log-likelihood for complement
      # check if any constraint is violated
      if (check.ciq || check.ceq) {    
        if (type %in% c("goric", "goricc")) { 
          llc <- logLik(ans$model.org)
          betasc <- b.unrestr
        } else if (type %in% c("gorica", "goricac")) {
          llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
          betasc <- b.unrestr
        }
        if (debug) {
          cat("log-likelihood_c value =", llc, "\n")
        }
        # if any constraints is violated LL_c = LL_u
      } else if (nrow(Amat) > meq && !(all(c(Amat) == 0L))) {
        ll <- list()
        betas <- list()
        # number of rows
        nr <- 1:nrow(Amat)
        # remove rows corresponding to equality constraints
        if (meq > 0L) {
          nr <- nr[-c(0:meq)]
        }
        # treat each row of Amat as an equality constraint
        for (l in 1:length(nr)) {
          idx <- c(nr[l], nr[-l])
          Amatx <- Amat[idx, , drop = FALSE]
          if (type %in% c("goric", "goricc")) {          
            Hc.restr <- restriktor(ans$model.org, constraints = Amatx, 
                                   neq = 1, rhs = bvec[idx], 
                                   mix.weights = "none", se = "none")
            betas[[l]] <- coef(Hc.restr)
            ll[[l]]    <- logLik(Hc.restr)
          } else if (type %in% c("gorica", "goricac")) {
            ldots4 <- ldots[m.restr > 0L]
            ldots4$mix.weights <- "none"
            CALL.restr <- c(list(object      = b.unrestr,
                                 constraints = Amatx,
                                 rhs         = bvec[idx],
                                 neq         = 1,
                                 VCOV        = VCOV),
                            ldots4)
            Hc.restr   <- do.call("con_gorica_est", CALL.restr) 
            betas[[l]] <- Hc.restr$b.restr
            ll[[l]]    <- dmvnorm(c(b.unrestr - Hc.restr$b.restr), 
                                  sigma = VCOV, log = TRUE)            
          }
        }
        if (debug) {
          cat("log-likelihood value =", ll[[l]], "\n")
        }
        # take the highest log-likelihood value as a substitute for the complement
        ll.unlist <- unlist(ll)
        ll.idx <- which(ll.unlist == max(ll.unlist))
        llc <- max(ll.unlist)
        betasc <- betas[[ll.idx]]
      } else if (nrow(Amat) == meq) { 
        # redundant, this will be catched by the first statement. In case of equality 
        # constraints only, the complement is equal to the unconstrained log-likelihood
        if (type %in% c("goric", "goricc")) {
          llc <- logLik(ans$model.org)
          betasc <- b.unrestr
        } else if (type %in% c("gorica", "goricac")) {
          llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
          betasc <- b.unrestr
        }
        if (debug) {
          cat("log-likelihood_c value =", llc, "\n")
        }
      } else if (all(c(Amat) == 0L)) {
        # unconstrained setting
        stop("restriktor ERROR: For an unconstrained hypothesis no complement exists.")
      } else {
        stop("restriktor ERROR: you might have found a bug, please contact me at: lgf.vanbrabant@gmail.com!")
      }
      if (type %in% c("goric", "goricc")) {
        llm <- logLik(conList[[1]])
      } else if (type %in% c("gorica", "goricac")) {
        llm <- dmvnorm(c(b.unrestr - b.restr), sigma = VCOV, log = TRUE)
      }
    } 
    
    # compute the number of free parameters f in the complement
    p   <- ncol(VCOV)
    # rank q1
    lq1 <- qr(Amat.ciq)$rank
    # rank q2
    lq2 <- qr(Amat.ceq)$rank
    # free parameters. Note that Amat includes q1 and q2 constraints
    f   <- p - qr(Amat)$rank # p - q1 - q2
    if (debug) { cat("number of free parameters =", (f + lq2), "\n") }
    
    # compute penalty term value PTc
    if (type %in% c("goric", "gorica")) {
      idx <- length(wt.bar)
      if (attr(wt.bar, "method") == "boot") {
        PTc <- as.numeric(1 + p - wt.bar[idx-meq] * lq1)
      } else if (attr(wt.bar, "method") == "pmvnorm") {
        # here, the q2 equalities are not included in wt.bar. Hence, they do not 
        # have to be subtracted.
        PTc <- as.numeric(1 + p - wt.bar[idx] * lq1) - length(bound)            # Dit moet nog voor de overige PTc gefixed worden!                             
      } else {
        stop("restriktor ERROR: no level probabilities (chi-bar-square weights) found.")
      }
    } else if (type %in% c("goricc", "goricac")) {
      idx <- length(wt.bar) 
      if (is.null(sample.nobs)) {
        stop("restriktor ERROR: the argument sample.nobs is not found.")
      }
      N <- sample.nobs
      # small sample correction
     if (attr(wt.bar, "method") == "boot") {
       PTc <- 1 + wt.bar[idx-meq] * (N * (p - lq1) + (p - lq1) + 2) / (N - (p - lq1) - 2) + 
         (1 - wt.bar[idx-meq]) * (N * p + p + 2) / (N - p - 2)
     } else if (attr(wt.bar, "method") == "pmvnorm") {
        PTc <- 1 + wt.bar[idx] * (N * (p - lq1) + (p - lq1) + 2) / (N - (p - lq1) - 2) + 
          (1 - wt.bar[idx]) * (N * p + p + 2) / (N - p - 2) 
     }
    }
    if (debug) {
      cat("penalty term value =", PTc, "\n")
    }
    
    # correct PTc for gorica(c)
    if (type %in% c("gorica", "goricac")) {
     PTc <- PTc - 1
    }
  } 
    
  ## compute loglik-value, goric(a)-values, and PT-values if comparison = unconstrained
  if (comparison == "unconstrained") { 
    PTm <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    if (type %in% c("goric", "goricc")) {
      # model
      llm <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
      # unrestricted
      llu <- logLik(ans$model.org)
    } else if (type %in% c("gorica", "goricac")) {
      # model 
      if (inherits(object, "numeric")) {
        llm <- unlist(lapply(conList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
                                                            sigma = VCOV, log = TRUE) )) 
        # unrestricted
        llu <- dmvnorm(rep(0, ncol(VCOV)), sigma = VCOV, log = TRUE) 
      } else { 
        llm <- unlist(lapply(conList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
                                                      sigma = VCOV, log = TRUE) )) 
        # unrestricted
        llu <- dmvnorm(rep(0, ncol(VCOV)), sigma = VCOV, log = TRUE) 
      }
    }
    
    if (type %in% c("goric", "gorica")) {
      PTu <- 1 + ncol(VCOV)
    } else if (type %in% c("goricc", "goricac")) {
      if (is.null(sample.nobs)) {
        stop("restriktor ERROR: the argument sample.nobs is not found.")
      }
      N   <- sample.nobs
      PTu <- ( (N * (ncol(VCOV) + 1) / (N - ncol(VCOV) - 2) ) ) 
    }
    
    ## correct PT for gorica(c)
    if (type %in% c("gorica", "goricac")) {
      PTu <- PTu - 1
    }
    
    goric.Hm <- -2*(llm - PTm)
    goric.Hu <- -2*(llu - PTu)
    df.Hm <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                        goric = goric.Hm)
    df.Hm$model <- as.character(df.Hm$model)
    df.u <- data.frame(model = "unconstrained", loglik = llu, penalty = PTu, 
                       goric = goric.Hu)
    df <- rbind(df.Hm, df.u)
    names(df)[4] <- type
  } else if (comparison == "none") {
    ## compute loglik-value, goric(a)-values, and PT-values if comparison = none
    PT <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    if (type %in% c("goric", "goricc")) {
      ll <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
      goric.Hm <- unlist(lapply(isSummary, function(x) x$goric[1]))
      df <- data.frame(model = objectnames, loglik = ll, penalty = PT, 
                       goric = goric.Hm)
      df$model <- as.character(df$model)
    } else if (type %in% c("gorica", "goricac")) {
      ll <- unlist(lapply(conList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
                                                         sigma = VCOV, log = TRUE) )) 
      goric.Hm <- -2*(ll - PT)
      df <- data.frame(model = objectnames, loglik = ll, penalty = PT, 
                       gorica = goric.Hm)
      df$model <- as.character(df$model)
    }
  } else if (comparison == "complement") {
    ## compute loglik-value, goric(a)-values, and PT-values if comparison = complement
    PTm <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    
    if (debug) {
     print(PTm)
    }
    
    # The PTm should not be corrected here! The PTm is already corrected in the 
    # restriktor.summary() function.
    if (type %in% c("goric", "goricc")) {
      # model
      goric.Hm <- -2*(llm - PTm)
      df.Hm  <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                           goric = goric.Hm)
      df.Hm$model <- as.character(df.Hm$model)
      # complement
      goric.Hc <- -2*(llc - PTc)
      df.c  <- data.frame(model = "complement", loglik = llc, penalty = PTc, 
                          goric = goric.Hc)
      df <- rbind(df.Hm, df.c)
    } else if (type %in% c("gorica", "goricac")) {
      # model
      gorica.Hm <- -2*(llm - PTm)
      df.Hm  <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                           gorica = gorica.Hm)
      df.Hm$model <- as.character(df.Hm$model)
      # complement
      gorica.Hc <- -2*(llc - PTc)
      df.c  <- data.frame(model = "complement", loglik = llc, penalty = PTc, 
                          gorica = gorica.Hc)
      df <- rbind(df.Hm, df.c)
    }
  } else {
    stop("restriktor ERROR: I cannot compute goric-values.")
  }

  ans$objectList  <- conList
  ans$objectNames <- objectnames
  # compute goric weights and relative weights
  delta <- df$goric - min(df$goric)
  goric.weights <- exp(-delta / 2) / sum(exp(-delta / 2))
  df$goric.weights <- goric.weights
    names(df)[5] <- paste0(type, ".weights")
  ans$result <- df

  # compute relative weights
  modelnames <- as.character(df$model)
  if (length(modelnames) > 1) {
    rw <- goric.weights %*% t(1/goric.weights)
    # it might happen that a diagonal value results in NaN.
    diag(rw) <- 1
    rownames(rw) <- modelnames
    colnames(rw) <- paste0("vs. ", modelnames)
    ans$relative.gw <- rw
  }
  
  # if (comparison == "complement" && is.null(bound)) {
  #   ans$ormle$b.restr.complement <- betasc
  # }
  
  # list all object estimates
  coefs <- lapply(conList, FUN = function(x) { coef(x) } )
  max.length <- max(sapply(coefs, length))
  coefs <- lapply(coefs, function(v) { c(v, rep(NA, max.length-length(v)))})
  coefs <- as.data.frame(do.call("rbind", coefs))
  if (comparison == "complement") {
    coefs <- rbind(coefs, betasc)
    rownames(coefs) <- c(objectnames, "complement")
  } else if (comparison == "unconstrained") {
    coefs <- rbind(coefs, conList[[1]]$b.unrestr)
    rownames(coefs) <- c(objectnames, "unconstrained")
  } else {
    rownames(coefs) <- objectnames    
  }
  
  # list of constraints
  ans$constraints <- lapply(conList, FUN = function(x) { x$constraints } )
    names(ans$constraints) <- ans$objectNames
  ans$rhs <- lapply(conList, FUN = function(x) { x$rhs } )
    names(ans$rhs) <- ans$objectNames
  ans$neq  <- lapply(conList, FUN = function(x) { x$neq } )
    names(ans$neq) <- ans$objectNames
  
  ans$ormle$b.restr <- coefs  
  ans$comparison <- comparison
  ans$type <- type
  ans$messages$mix_weights <- do.call("rbind", lapply(isSummary, FUN = function(x) { x$messages$mix_weights }))
  ans$messages$mix_weights <- ans$messages$mix_weights [!duplicated(ans$messages$mix_weights )]
  
  if (type %in% c("goric", "goricc")) {
    class(ans) <- "con_goric"
  } else if (type %in% c("gorica", "goricac")) {
    class(ans) <- c("con_gorica", "con_goric")
  }
  
  ans
}



goric.lavaan <- function(object, ...,
                         comparison = "unconstrained",
                         type = "gorica",
                         standardized = FALSE,
                         bound = NULL, debug = FALSE) {
  
  if (!inherits(object, "lavaan")) {
    stop("restriktor ERROR: the object must be of class lavaan.")
  }
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop("restriktor ERROR: object of class lavaan is only supported for type = 'gorica(c)'.")
  }
  
  objectList <- list(...)
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  mcList$object       <- NULL
  mcList$comparison   <- NULL
  mcList$type         <- NULL
  mcList$standardized <- NULL
  
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  
  est <- con_gorica_est_lav(object, standardized)
  objectList$object       <- est$estimate
  objectList$VCOV         <- est$VCOV
  objectList$comparison   <- comparison
  objectList$type         <- type
  objectList$bound        <- bound
  objectList$debug <- debug
  if (type == "goricac") {
    objectList$sample.nobs <- lavInspect(object, what = "ntotal")
  }
  res <- do.call(goric.default, objectList)
  
  res
}


goric.lm <- function(object, ...,
                     comparison = "unconstrained",
                     type = "goric",
                     bound = NULL, debug = FALSE) {
  
  if (!inherits(object, "lm")) {
    stop("restriktor ERROR: the object must be of class lm, glm, mlm, rlm.")
  }
  
  
  objectList <- list(...)
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  mcList$object       <- NULL
  mcList$comparison   <- NULL
  mcList$type         <- NULL
  
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  
  objectList$object       <- object
  objectList$comparison   <- comparison
  objectList$type         <- type
  objectList$bound        <- bound
  objectList$debug        <- debug
  if (type == "goricac") {
    objectList$sample.nobs <- length(residuals(object))
  }
  res <- do.call(goric.default, objectList)
  
  res
}


goric.restriktor <- function(object, ...,
                             comparison = "unconstrained",
                             type = "goric",
                             bound = NULL, debug = FALSE) {
  
  if (!inherits(object, "restriktor")) {
    stop("restriktor ERROR: the object must be of class restriktor.")
  }
  
  objectList <- list(...)
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  mcList$comparison   <- NULL
  mcList$type         <- NULL
  
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  
  objectList$comparison   <- comparison
  objectList$type         <- type
  objectList$bound        <- bound
  objectList$debug        <- debug
  if (type == "goricac") {
    objectList$sample.nobs <- length(residuals(object))
  }
  res <- do.call(goric.default, objectList)
  
  res
}


goric.numeric <- function(object, ...,
                          VCOV = NULL,
                          comparison = "unconstrained",
                          type = "gorica", sample.nobs = NULL,
                          bound = NULL, debug = FALSE) {
  
  if (!inherits(object, "numeric")) {
    stop("restriktor ERROR: the object must be of class numeric.")
  }
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop("restriktor ERROR: object of class numeric is only supported for type = 'gorica(c)'.")
  }
  
  objectList <- list(...)
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  mcList$object       <- NULL
  mcList$comparison   <- NULL
  mcList$type         <- NULL
  mcList$VCOV         <- NULL
  
  #<FIXME>
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  #</FIXME>
  
  objectList$object       <- object
  objectList$VCOV         <- VCOV
  objectList$comparison   <- comparison
  objectList$type         <- type
  objectList$bound        <- bound
  objectList$debug        <- debug
  if (type == "goricac") {
    objectList$sample.nobs <- sample.nobs
  }
  res <- do.call(goric.default, objectList)
  
  res
}


goric.CTmeta <- function(object, ...,
                         VCOV = NULL,
                         comparison = "unconstrained",
                         type = "gorica", sample.nobs = NULL,
                         bound = NULL, debug = FALSE) {
  
  if (!inherits(object, "CTmeta")) {
    stop("restriktor ERROR: the object must be of class CTmeta.")
  }
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop("restriktor ERROR: object of class CTmeta is only supported for type = 'gorica(c)'.")
  }
  
  objectList <- list(...)
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  mcList$object       <- NULL
  mcList$comparison   <- NULL
  mcList$type         <- NULL
  mcList$VCOV         <- NULL
  
  #<FIXME>
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  #</FIXME>
  
  objectList$object       <- coef(object)
  objectList$VCOV         <- vcov(object)
  objectList$comparison   <- comparison
  objectList$type         <- type
  objectList$bound        <- bound
  objectList$debug        <- debug
  if (type == "goricac") {
    objectList$sample.nobs <- sample.nobs
  }
  res <- do.call(goric.default, objectList)
  
  res
}



print.con_goric <- function(x, digits = max(3, getOption("digits") - 4), ...) {

  type <- x$type
  comparison <- x$comparison
  dig <- paste0("%6.", digits, "f")
  x2 <- lapply(x$result[-1], sprintf, fmt = dig)
  x2$model <- x$result$model
  df <- as.data.frame(x2)
  df <- df[, c(5,1,2,3,4)]
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))

  if (type == "goric") {
    cat("generalized order-restriced information criterion: \n")
  } else if (type == "gorica") {
    cat("generalized order-restriced information criterion approximation:\n")
  } else if (type == "goricc") {
    cat("small sample generalized order-restriced information criterion:\n")
  } else if (type == "goricac") {
    cat("small sample generalized order-restriced information criterion approximation:\n")
  }
  
  cat("\nResults:\n")
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  cat("---")
  if (comparison == "complement") {
    objectnames <- as.character(df$model)
    relative.gw <- apply(x$relative.gw, 2, sprintf, fmt = dig)
    cat("\nThe order-restricted hypothesis", sQuote(objectnames[1]), "has", 
        sprintf("%s", relative.gw[1,2]), "times more support than its complement.\n")
  } 
  
  cat("\n")
  cat(x$messages$mix_weights)
  
  invisible(x)
}




summary.con_goric <- function(object, brief = TRUE, 
                          digits = max(3, getOption("digits") - 4), ...) {
  
  x <- object
  type <- x$type
  comparison <- x$comparison
  dig <- paste0("%6.", digits, "f")
  x2 <- lapply(x$result[-1], sprintf, fmt = dig)
  x2$model <- x$result$model
  df <- as.data.frame(x2)
  df <- df[, c(5,1,2,3,4)]
  
  objectnames   <- as.character(df$model)
  #goric.weights <- as.numeric(as.character(df$goric.weights))

  Amat <- x$constraints
  meq  <- x$neq
  bvec <- x$rhs
  iact <- lapply(x$objectList, FUN = function(x) { x$iact } )

  if (type == "goric") {
    cat("\nrestriktor: generalized order-restriced information criterion: \n")
  } else if (type == "gorica") {
    cat("\nrestriktor: generalized order-restriced information criterion approximation:\n")
  } else if (type == "goricc") {
    cat("\nrestriktor: small sample generalized order-restriced information criterion:\n")
  } else if (type == "goricac") {
    cat("\nrestriktor: small sample generalized order-restriced information criterion approximation:\n")
  }
  
  cat("\nResults:\n")  
  
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  cat("---\n")
  
  if (!is.null(x$relative.gw)) {
    if (type == "goric") {
      cat("\nRelative GORIC-weights:\n")
    } else if (type == "gorica") {
      cat("\nRelative GORICA-weights:\n")
    } else if (type == "goricc") {
      cat("\nRelative GORICC-weights:\n")
    } else if (type == "goricac") {
      cat("\nRelative GORICAC-weights:\n")
    }
    
    relative.gw <- apply(x$relative.gw, 2, sprintf, fmt = dig)
    
    rownames(relative.gw) <- rownames(x$relative.gw)
    if (max(as.numeric(relative.gw)) >= 1e4) {
      print(format(x$relative.gw, digits = digits, scientific = TRUE), 
            print.gap = 2, quote = FALSE)
    } else {
      print(format(trimws(relative.gw), digits = digits, scientific = FALSE), 
            print.gap = 2, quote = FALSE)
    }
    cat("---\n")
    if (length(unique(df$loglik)) != length(df$loglik)) {
      cat("Note: In case of equal log-likelihood (loglik) values, the 
      relative weights are solely based on the difference in penalty values.\n")
    }
  }
  
  if (comparison == "complement") {
    cat("The order-restricted hypothesis", sQuote(objectnames[1]), "has", 
        sprintf("%s", relative.gw[1,2]), "times more support than its complement.\n")
  } 
  
  if (!brief) {
    cat("\n\nOrder/Inequality restricted coefficients:\n")
    coefs <- trimws(apply(x$ormle$b.restr, 2, sprintf, fmt = dig))
    rownames(coefs) <- rownames(x$ormle$b.restr)
    print(format(coefs, digits = digits, scientific = FALSE), print.gap = 2L,
          quote = FALSE)
    cat("---\n")
    
    vnames <- names(x$ormle$b.restr)
    fn <- function(Amat, bvec, meq, iact, vnames) {
      colnames(Amat) <- vnames
      out.rest <- cbind(round(Amat, 4), c(rep("   ==", meq), rep("   >=", nrow(Amat) - 
                                                                   meq)), bvec, " ")
      rownames(out.rest) <- paste(1:nrow(out.rest), ":", sep = "")
      colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs", "active")
      idx <- ncol(out.rest)
      out.rest[, idx] <- "no"
      out.rest[iact, idx] <- "yes"
      if (nrow(Amat) == meq) {
        out.rest[1:nrow(Amat), idx] <- "yes"
      }  
      out.rest <- as.data.frame(out.rest)
      
      out.rest
    }
    
    conMat <- list()
    for (i in 1:length(x$objectList)) {
      conMat[[i]] <- fn(Amat = Amat[[i]], bvec = bvec[[i]], meq = meq[[i]], 
                        iact = iact[[i]], vnames = vnames)  
    }
    names(conMat) <- x$objectNames
    
    if (comparison == "complement") {
     conMat$complement <- paste("not", x$objectNames) 
    }
    
    cat("\nConstraint matrices:\n")
    print(conMat, quote = FALSE, scientific = FALSE)
   
    invisible(x)
  }
  cat("\n")
  cat(x$messages$mix_weights)
}


