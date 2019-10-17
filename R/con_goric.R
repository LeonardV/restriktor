goric <- function(object, ..., 
                  comparison = c("unconstrained", "complement", "none"), 
                  VCOV = NULL, sample.nobs = NULL,
                  type = "goric", bound = NULL, debug = FALSE) {
  
  mc <- match.call()
  CALL <- as.list(mc)
  CALL[[1]] <- NULL
  
  # some checks
  comparison <- tolower(comparison)
  comparison <- match.arg(comparison)
  stopifnot(comparison %in% c("unconstrained", "complement", "none"))
  
  type <- tolower(type)
  stopifnot(type %in% c("goric", "goricc", "gorica", "goricca"))
  
  ldots <- list(...)
  m.restr <- match(names(ldots), c("B", "mix.weights", "mix.bootstrap", 
                                   "parallel", "ncpus", "cl", "seed", "control", 
                                   "verbose", "debug"), 0L)
  if (length(m.restr) == 0L) {
    ldots2 <- ldots
  } else {
    ldots2 <- ldots[m.restr == 0L]
  }  

  if (inherits(object, "matrix")) {
    if (dim(object)[2] > 1L) {
      stop("Restriktor ERROR: object must be a vector.")
    }
    object <- as.vector(object)
  }
  
  # what is the constraint input type
  objectList  <- c(list(object), ldots2)
  isRestr     <- sapply(objectList, function(x) inherits(x, "restriktor"))
  isCharacter <- sapply(ldots2,     function(x) inherits(x, "character"))
  isList      <- sapply(ldots2,     function(x) inherits(x, "list"))
  
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
    idx  <- length(conList) 
    objectnames <- vector("character", idx)
    for (i in 1:idx) {
      if (length(as.character(CALL[[i]])) > 1) {
        CALL[[i]] <- paste0("H", i)
      }
      objectnames[i] <- as.character(CALL[[i]])
    }
  # if the constraints syntax is of class character, e.g., x1 < x2; x2 < x3
  } else if ( inherits(object, "lm") && all(isCharacter) ) {
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
    for (i in 1:idx) {
      if (length(as.character(CALL[[i]])) > 1) {
        CALL[[i]] <- paste0("H", i)
      }
      objectnames[i] <- as.character(CALL[[i]])
    }
  # if the constraints are a list with constraints, rhs and neq for each hypothesis   
  } else if (inherits(object, "lm") && (all(isList))) {
      # create lists
      #constraints <- list(); rhs <- list(); neq <- list()
      # extract constraints, rhs and neq
      constraints <- lapply(ldots2, FUN = function(x) {x$constraints} )
      constraints.check <- sapply(constraints, FUN = function(x) { is.null(x) } )
      if (any(constraints.check)) {
        stop("Restriktor ERROR: the constraints must be specified as a list. E.g., h1 <- list(constraints = 'x1 > 0')")
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
      for (i in 1:idx) {
        if (length(as.character(CALL[[i]])) > 1) {
          CALL[[i]] <- paste0("H", i)
        }
        objectnames[i] <- as.character(CALL[[i]])
      }
    } else if (inherits(object, "numeric")) {
      if (type %in% c("goric", "goricc")) {
        stop("Restriktor ERROR: object of class numeric is only supported for type = 'goric(c)a'.")
      }
      if (is.null(VCOV)) {
        stop("Restriktor ERROR: the argument VCOV is not found.")
      }
      if (is.null(sample.nobs) && type %in% c("goricc", "goricca")) {
        stop("Restriktor ERROR: the argument sample.nobs is not found.")
      }
      
      isCharacter <- sapply(ldots2, function(x) inherits(x, "character"))    
      isList      <- sapply(ldots2, function(x) inherits(x, "list"))
      
      ldots3 <- ldots[m.restr > 0L]
      if (all(isList)) { 
        # create lists
        #constraints <- list(); rhs <- list(); neq <- list()
        # extract constraints, rhs and neq
        constraints       <- lapply(ldots2,      FUN = function(x) { x$constraints } )
        constraints.check <- sapply(constraints, FUN = function(x) { is.null(x)} )
        if (any(constraints.check) || length(constraints) == 0) {
          stop("Restriktor ERROR: no constraints found! The constraints must be specified as a list. E.g., h1 <- list(constraints = 'x1 > 0')")
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
      
      idx <- length(conList) 
      objectnames <- vector("character", idx)
      for (i in 1:idx) { 
        if (length(as.character(CALL[[i]])) > 1) {
          CALL[[i]] <- paste0("H", i)
        }
        objectnames[i] <- as.character(CALL[[i]])
      }
    } else {
      stop("Restriktor ERROR: I don't know how to handle an object of class ", paste0(class(object)[1]))
    }


  if (comparison == "complement" && length(conList) > 1L) {
    comparison <- "unconstrained"
    warning("Restriktor WARNING: if comparison = 'complement', only one order-restricted hypothesis\n",
            "                      is allowed (for now). Therefore, comparison is set to 'unconstrained'.")
  } 


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
              "                      and are therefore ignored.")
      bound <- NULL
    } 
    
    # for now
    bound <- NULL
    if (!is.null(bound)) {
      # if (meq > 0) {
      #   Amat.ceq <- Amat[1:meq, ,drop = FALSE]
      #   bvec.ceq <- bvec[1:meq]
      #   Amat.ciq <- Amat[-c(1:meq), , drop = FALSE]
      #   bvec.ciq <- bvec[-c(1:meq)]
      # } else {
      #   Amat.ceq <- matrix( , nrow = 0, ncol = ncol(Amat))
      #   bvec.ceq <- rep(0, 0)
      #   Amat.ciq <- Amat[ , , drop = FALSE]
      #   bvec.ciq <- bvec
      # }
      # 
      # # upper-bound = rhs + bound
      # ub <- bvec.ceq + bound
      # # lower-bound = rhs - bound
      # lb <- bvec.ceq - bound
      # 
      # ## check if any equality constraints
      # if (meq > 0L) {
      #   # correct user error
      #   if ( ((length(ub) == 1L) || (length(lb) == 1L)) && meq >= 1L ) {
      #     ub <- rep(ub, meq)
      #     lb <- rep(lb, meq)
      #   }
      #   # check
      #   if ( (length(ub) != meq) || (length(lb) != meq) ) {
      #     stop("restriktor ERROR: the number of bounds is not equal to number of equality constraints (neq).")
      #   }
      # 
      #   # check if unconstrained mle are violated
      #   check.ciq <- all(Amat.ciq %*% c(b.unrestr) - bvec.ciq >= 0)
      #   # check if unconstrained mle lay between the bounds
      #   check.ub <- all(Amat.ceq %*% c(b.unrestr) <= (ub + .Machine$double.eps))
      #   check.lb <- all(Amat.ceq %*% c(b.unrestr) >= (lb - .Machine$double.eps))
      # 
      #   # check if unrestricted mle lay in boundary area
      #   if (check.ciq && check.ub && check.lb) {
      #     # log-likelihood model
      #     if (type == "goric") {
      #       llm <- logLik(ans$model.org)
      #     } else if (type == "gorica") {
      #       llm <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
      #     }
      #     # determine log-likelihood_c
      #     nr.meq <- 1:meq
      #     # upper-bound
      #     llc.ub <- list()
      #     for (l in nr.meq) {
      #       Amat.ub <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
      #       bvec.ub <- c(ub[l], bvec.ciq)
      #       Hc.ub <- restriktor(ans$model.org, constraints = Amat.ub,
      #                           neq = 1, rhs = bvec.ub,
      #                           mix.weights = "none", se = "none")
      #       if (type == "goric") {  
      #         llc.ub[[l]] <- logLik(Hc.ub)
      #       } else if (type == "gorica") {
      #         llc.ub[[l]] <- dmvnorm(c(b.unrestr - Hc.ub$b.restr), 
      #                                sigma = VCOV, log = TRUE)
      #       }
      #     }
      #     # lower-bound
      #     llc.lb <- list()
      #     for (l in nr.meq) {
      #       Amat.lb <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
      #       bvec.lb <- c(lb[l], bvec.ciq)
      #         Hc.lb <- restriktor(ans$model.org, constraints = Amat.lb,
      #                             neq = 1, rhs = bvec.lb,
      #                             mix.weights = "none", se = "none")
      #       if (type == "goric") {
      #         llc.lb[[l]] <- logLik(Hc.lb)
      #       } else if (type == "gorica") {
      #         llc.lb[[l]] <- dmvnorm(c(b.unrestr - Hc.lb$b.restr), 
      #                                sigma = VCOV, log = TRUE)
      #       }
      #     }
      # 
      #     ## only if any ciq
      #     if (length(Amat.ciq) > 0L) {
      #       Amatx <- rbind(Amat.ciq, Amat.ceq, -Amat.ceq)        # correct?
      #       bvecx <- c(bvec.ciq, lb, lb)                         # correct?
      # 
      #       llc.s <- list()
      #       nr.ciq <- 1:nrow(Amat.ciq)
      #       for (l in nr.ciq) {
      #         Hc.s <- restriktor(ans$model.org, 
      #                            constraints = Amat.ciq[l, , drop = FALSE],
      #                            neq = 1, rhs = bvec.ciq[l],
      #                            mix.weights = "none", se = "none")
      #         b.s <- coef(Hc.s)
      #         b.s[abs(b.s) < sqrt(.Machine$double.eps)] <- 0L
      #         if (all(Amatx %*% c(b.s) >= bvecx)) {               # correct check?
      #           # log-likelihood model
      #           if (type == "goric") {
      #             llc.s[[l]] <- logLik(Hc.s)
      #           } else if (type == "gorica") {
      #             llc.s[[l]] <- dmvnorm(c(b.unrestr - Hc.s$b.restr), 
      #                                   sigma = VCOV, log = TRUE)
      #           }
      #         }
      #       }
      #     } else {
      #       llc.s <- NULL
      #     }
      # 
      #     llc <- unlist(c(llc.ub, llc.lb, llc.s))
      #     llc <- max(llc)
      # 
      #     if (debug) {
      #       print(llc)
      #     }
      #   } else {
      #     # determine the log-likelihood_m in case the unconstrained mle
      #     # lay outside the range restrictions.
      #     # 2^q2 combinations.
      #     if (type == "goric") {
      #       llc <- logLik(ans$model.org)
      #     } else if (type == "gorica") {
      #       llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
      #     }
      # 
      #     # if all bounds are zero, then llm = logLik(object)
      #     if (!all(bound == 0L)) {
      #       # get all combinations for ub and lb
      # 
      #       len.bvec.ceq <- length(ub)
      #       perm <- rep(list(rep(c(-1,1), (2^len.bvec.ceq)/2)), len.bvec.ceq)
      #       perm.grid <- unique(as.matrix(expand.grid(perm)))
      #       nr.perm <- 1:(2^len.bvec.ceq)
      #       # if ub = lb, then the bounds must be zero.
      #       # Otherwise, +ub, -lb
      #       bound.zero.idx <- which(ub == lb)
      # 
      #       llm <- list()
      #       for (m in nr.perm) {
      #         #perm.vec <- perm.grid[m, ]
      #         perm.vec <- apply(Amat.ceq, 2, function(x) perm.grid[m, ] * x)
      #         ub.idx <- which(perm.grid[m, ] == -1)
      #         lb.idx <- which(perm.grid[m, ] ==  1)
      #         order.idx <- c(ub.idx, lb.idx)
      #         bounds.new <- c(-ub[ub.idx], lb[lb.idx])
      #         bounds.new <- bounds.new[order(order.idx, decreasing = FALSE)]
      # 
      #         Amatx <- rbind(Amat.ciq, Amat.ceq, -Amat.ceq)
      #         bvecx <- c(bvec.ciq, bounds.new, bounds.new)
      # 
      #         #Amat.ceq.perm <- t( t(Amat.ceq) * perm.vec )
      #         Amat.ceq.perm <- perm.vec
      #         Amat.new <- rbind(Amat.ceq.perm, Amat.ciq)
      #         bvec.new <- c(bounds.new, bvec.ciq)
      # 
      #         if (length(bound.zero.idx) == 0L) {
      #           Amat.new.sort <- Amat.new
      #           bvec.new.sort <- bvec.new
      #         } else {
      #           # constraint rows for which bound == 0 are placed on top
      #           Amat.new.sort <- rbind(Amat.new[bound.zero.idx, , drop = FALSE],
      #                                  Amat.new[-bound.zero.idx, , drop = FALSE])
      #           bvec.new.sort <- c(bvec.new[bound.zero.idx],
      #                              bvec.new[-bound.zero.idx])
      #         }
      # 
      #         Hm <- restriktor(ans$model.org, constraints = Amat.new.sort,
      #                          neq = length(bound.zero.idx), rhs = bvec.new.sort,
      #                          mix.weights = "none", se = "none")
      #         beta.Hm <- coef(Hm)
      #         beta.Hm[abs(beta.Hm) < sqrt(.Machine$double.eps)] <- 0L
      #         if (all( Amatx %*% c(beta.Hm) - bvecx + sqrt(.Machine$double.eps) >= 0 )) {
      #           if (type == "goric") {
      #             llm[[m]] <- logLik(Hm)
      #           } else if (type == "gorica") {
      #             llm[[m]] <- dmvnorm(c(Hm$b.unrestr - Hm$b.restr),
      #                                 sigma = VCOV, log = TRUE)
      #           }
      #         }
      #       }
      #       llm <- unlist(llm)
      #       if (debug) {
      #         cat("log-likelihood_m =", llm, "\n")
      #       }
      # 
      #       llm <- max(llm)
      #     } else {
      #       if (type == "goric") {
      #         llm <- logLik(conList[[1]])
      #       } else if (type == "gorica") {
      #         llm <- dmvnorm(c(b.unrestr - b.restr),
      #                        sigma = VCOV, log = TRUE)
      #       }
      #     }
      #     llm
      #   }
      # }
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
        } else if (type %in% c("gorica", "goricca")) {
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
          } else if (type %in% c("gorica", "goricca")) {
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
        } else if (type %in% c("gorica", "goricca")) {
          llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
          betasc <- b.unrestr
        }
        if (debug) {
          cat("log-likelihood_c value =", llc, "\n")
        }
      } else if (all(c(Amat) == 0L)) {
        # unconstrained setting
        stop("Restriktor ERROR: no complement exists for the unconstrained hypothesis.")
      } else {
        stop("Restriktor ERROR: you might have found a bug, please contact me!")
      }
      if (type %in% c("goric", "goricc")) {
        llm <- logLik(conList[[1]])
      } else if (type %in% c("gorica", "goricca")) {
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
        PTc <- as.numeric(1 + p - wt.bar[idx] * lq1)
      } else {
        stop("Restriktor ERROR: no level probabilities (chi-bar-square weights) found.")
      }
    } else if (type %in% c("goricc", "goricca")) {
      idx <- length(wt.bar) 
      if (is.null(sample.nobs)) {
        stop("Restriktor ERROR: the argument sample.nobs is not found.")
      }
      N <- sample.nobs
      # small sample correction
      if (attr(wt.bar, "method") == "boot") {
        lPT <- 0 : ncol(Amat)
        PTc  <- sum( ( (N * (lPT + 1) / (N - lPT - 2) ) ) * wt.bar[idx-meq])
      } else if (attr(wt.bar, "method") == "pmvnorm") {
        # min.col <- ncol(Amat) - nrow(Amat) # p - q1 - q2
        # max.col <- ncol(Amat) - meq        # p - q2
        # lPT <- min.col : max.col
        # PTc <- 1 + p - sum( ((N * (lPT + 1) / (N - lPT - 2))) * wt.bar[idx])   
        
        PTu <- ( (N * (ncol(VCOV) + 1) / (N - ncol(VCOV) - 2) ) ) 
        PTm <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
        PTc <- 1 + p + log(1 - exp(PTm - PTu))
      }
    }
    if (debug) {
      cat("penalty term value =", PTc, "\n")
    }
  } 
    
  # compute loglik-value, goric(a)-values, and PT-values
  if (comparison == "unconstrained") { 
    PTm <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    if (type %in% c("goric", "goricc")) {
      # model
      llm <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
      # unrestricted
      llu <- logLik(ans$model.org)
    } else if (type %in% c("gorica", "goricca")) {
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
    } else if (type %in% c("goricc", "goricca")) {
      lPT <- ncol(VCOV)
      if (is.null(sample.nobs)) {
        stop("Restriktor ERROR: the argument sample.nobs is not found.")
      }
      N   <- sample.nobs
      PTu <- sum( ( (N * (lPT + 1) / (N - lPT - 2) ) ) * 1L)
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
    PT <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    if (type %in% c("goric", "goricc")) {
      ll <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
      goric.Hm <- unlist(lapply(isSummary, function(x) x$goric[1]))
      df <- data.frame(model = objectnames, loglik = ll, penalty = PT, 
                       goric = goric.Hm)
      df$model <- as.character(df$model)
    } else if (type %in% c("gorica", "goricca")) {
#      if (inherits(object, "numeric")) {
        ll <- unlist(lapply(conList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
                                                           sigma = VCOV, log = TRUE) )) 
      # } else { 
      #   ll <- unlist(lapply(conList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
      #                                                      sigma = VCOV, log = TRUE) )) 
      # }
      goric.Hm <- -2*(ll - PT)
      df <- data.frame(model = objectnames, loglik = ll, penalty = PT, 
                       gorica = goric.Hm)
      df$model <- as.character(df$model)
    }
  } else if (comparison == "complement") {
    PTm <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
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
    } else if (type %in% c("gorica", "goricca")) {
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
    stop("Restriktor ERROR: I cannot compute goric-values.")
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
    colnames(rw) <- paste0("vs ", modelnames)
    ans$relative.gw <- rw
  }
  
  if (comparison == "complement" && is.null(bound)) {
    ans$ormle$b.restr.complement <- betasc
  }
  
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
  
  if (type %in% c("goric", "goricc")) {
    class(ans) <- "goric"
  } else if (type %in% c("gorica", "goricca")) {
    class(ans) <- c("gorica", "goric")
  }
  
  ans
}



print.goric <- function(x, digits = max(3, getOption("digits") - 4), ...) {

  type <- x$type
  comparison <- x$comparison
  dig <- paste0("%6.", digits, "f")
  x2 <- lapply(x$result[-1], sprintf, fmt = dig)
  x2$model <- x$result$model
  df <- as.data.frame(x2)
  df <- df[, c(5,1,2,3,4)]
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))

  if (type == "goric") {
    cat("\nRestriktor: generalized order-restriced information criterion: \n")
  } else if (type == "gorica") {
    cat("\nRestriktor: generalized order-restriced information criterion approximation:\n")
  } else if (type == "goricc") {
    cat("\nRestriktor: small sample generalized order-restriced information criterion:\n")
  } else if (type == "goricca") {
    cat("\nRestriktor: small sample generalized order-restriced information criterion approximation:\n")
  }
  
  cat("\nResults:\n")
  
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  cat("---\n")
  if (comparison == "complement") {
    objectnames <- as.character(df$model)
    relative.gw <- apply(x$relative.gw, 2, sprintf, fmt = dig)
    
    cat("The order-restricted hypothesis", sQuote(objectnames[1]), "has", 
        sprintf("%s", relative.gw[1,2]), "times more support than its complement.\n")
  } 
  invisible(x)
}



summary.goric <- function(object, brief = TRUE, 
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
    cat("\nRestriktor: generalized order-restriced information criterion: \n")
  } else if (type == "gorica") {
    cat("\nRestriktor: generalized order-restriced information criterion approximation:\n")
  } else if (type == "goricc") {
    cat("\nRestriktor: small sample generalized order-restriced information criterion:\n")
  } else if (type == "goricca") {
    cat("\nRestriktor: small sample generalized order-restriced information criterion approximation:\n")
  }
  
  cat("\nResults:\n")  
  
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  cat("---\n")
  
  if (!is.null(x$relative.gw)) {
    if (type == "goric") {
      cat("\n\nRelative GORIC-weights:\n")
    } else if (type == "gorica") {
      cat("\n\nRelative GORICA-weights:\n")
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
}



coef.goric <- function(object, ...)  {
  return(object$ormle$b.restr)
}

coef.gorica_est <- function(object, ...)  {
  return(object$b.restr)
}
