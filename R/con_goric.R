goric <- function(object, ...) { UseMethod("goric") }


goric.default <- function(object, ..., hypotheses = NULL,
                          comparison = c("unconstrained", "complement", "none"), 
                          VCOV = NULL, sample.nobs = NULL,
                          type = "goric", 
                          auto.bound = FALSE, debug = FALSE) {

  constraints <- hypotheses
  # class objects
  object_class <- unlist(lapply(object, class))
  
  # is object class restriktor
  if ("restriktor" %in% object_class && !is.null(constraints)) {
    warning("restriktor Warning: hypotheses are inherited from the restriktor object and are therefore ignored.", 
            call. = FALSE)
    constraints <- NULL
  }
  
  # auto.bounds are ignored (for now)
  auto.bound <- NULL
  
  mc <- match.call()
  CALL <- as.list(mc[-1])
  
  # some checks
  comparison <- tolower(comparison)
  comparison <- match.arg(comparison)
  stopifnot(comparison %in% c("unconstrained", "complement", "none"))
  
  type <- tolower(type)
  stopifnot(type %in% c("goric", "goricc", "gorica", "goricac"))
  
  ldots <- list(...)
  ldots$missing <- NULL

  conChar <- sapply(constraints, function(x) inherits(x, "character"))
  isConChar <- all(conChar)

  if (!"restriktor" %in% object_class) { 
    if (!is.list(constraints)) { 
      stop("Restriktor ERROR: The 'hypotheses' argument must be a named list. 
           Please provide hypotheses in the following format: 'list(H1 = H1)' or 
           'list(S1 = list(H11, H12), S2 = list(H21, H22))'", call. = FALSE)
    }
    
    # give constraints list a name if null
    if (is.null(names(constraints))) {
      names(constraints) <- paste0("H", 1:length(constraints)) 
    }

  }
# -------------------------------------------------------------------------
  # create output list
  ans <- list()
  
  ## deal with objects of different classes
  if ("restriktor" %in% object_class) {   
    # if all objects are of class restriktor
    conList   <- object
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     goric       = type,
                                                     sample.nobs = sample.nobs))
    ans$hypotheses_usr <- lapply(conList, function(x) x$CON$constraints)
    ans$model.org <- object[[1]]$model.org
    sample.nobs   <- nrow(model.frame(object[[1]]$model.org))
    # unrestricted VCOV
    VCOV <- vcov(ans$model.org)
    
    #constraints <- conList
  } else if (any(object_class %in% c("lm","rlm","glm","mlm")) && isConChar) { 
    # standard errors are not needed
    ldots$se <- "none"
    # fit restriktor object for each hypothesis
    conList <- list()
    for (con in 1:length(constraints)) {
      CALL.restr <- c(list(object      = object$object,
                           constraints = constraints[[con]]), ldots)
      conList[[con]] <- do.call("restriktor", CALL.restr)
    }
    names(conList) <- names(constraints)
    # compute summary for each restriktor object. 
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     goric       = type,
                                                     sample.nobs = sample.nobs))
    ans$hypotheses_usr <- lapply(conList, function(x) x$CON$constraints)
    # add unrestricted object to output
    ans$model.org <- object[[1]]
    # unrestricted VCOV
    VCOV <- vcov(ans$model.org)
    sample.nobs <- nrow(model.frame(object[[1]]))
    idx <- length(conList) 
    objectnames <- vector("character", idx)
    CALL$object <- NULL
  } else if (any(object_class %in% c("lm","rlm","glm","mlm")) && !isConChar) {
    # tolower names Amat and rhs
    for (i in seq_along(constraints)) { 
      names(constraints[[i]]) <- tolower(names(constraints[[i]]))
      if (any(!names(constraints[[i]]) %in% c("constraints", "rhs", "neq"))) { 
        stop("Restriktor ERROR: The list objects must be named 'constraints', 'rhs' and 'neq', e.g.:
              h1 <- list(constraints = c(0,1,0))
              h2 <- list(constraints = rbind(c(0,1,0), c(0,0,1)), rhs = c(0.5, 1), neq = 0)
              hypotheses = list(H1 = h1, H2 = h2).", 
             call. = FALSE)
      }
    }

    # standard errors are not needed
    ldots$se <- "none"
    # fit restriktor object for each hypothesis
    conList <- list()
    for (con in 1:length(constraints)) {
      CALL.restr <- append(list(object      = object$object,
                                constraints = constraints[[con]]$constraints,
                                rhs         = constraints[[con]]$rhs,
                                neq         = constraints[[con]]$neq), ldots)
      
      conList[[con]] <- do.call("restriktor", CALL.restr) 
    }
    names(conList) <- names(constraints)
    # compute symmary for each restriktor object. Here is the goric value 
    # computed. Note: not the gorica value
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     goric       = type,
                                                     sample.nobs = sample.nobs))
    # add unrestricted object to output
    ans$model.org <- object[[1]]
    # unrestricted VCOV
    VCOV <- vcov(ans$model.org) 
    sample.nobs <- nrow(model.frame(object[[1]]))
    idx <- length(conList) 
    objectnames <- vector("character", idx)
    CALL$object <- NULL
  } else if ("numeric" %in% object_class & isConChar) {
    # fit restriktor object for each hypothesis
    conList <- list()
    for (con in 1:length(constraints)) {
      CALL.restr <- append(list(object      = object$object, 
                                constraints = constraints[[con]],
                                VCOV        = as.matrix(VCOV)), ldots)
      conList[[con]] <- do.call("con_gorica_est", CALL.restr)  
    }
    names(conList) <- names(constraints)
    
    ans$hypotheses_usr <- lapply(conList, function(x) x$CON$constraints)
    
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     type        = type,
                                                     sample.nobs = sample.nobs)) 
    } else if ("numeric" %in% object_class & !isConChar) {
    # tolower names Amat and rhs
    for (i in seq_along(constraints)) { 
      names(constraints[[i]]) <- tolower(names(constraints[[i]]))
      if (any(!names(constraints[[i]]) %in% c("constraints", "rhs", "neq"))) {
        stop("Restriktor ERROR: The list objects must be named 'constraints', 'rhs' and 'neq', e.g.:
              h1 <- list(constraints = c(0,1,0))
              h2 <- list(constraints = rbind(c(0,1,0), c(0,0,1)), rhs = c(0.5, 1), neq = 0)
              hypotheses = list(H1 = h1, H2 = h2).", 
             call. = FALSE)
      }
    }
    
      conList <- list()
      for (con in 1:length(constraints)) {
        CALL.restr <- append(list(object      = object$object,
                                  VCOV        = as.matrix(VCOV),
                                  constraints = constraints[[con]]$constraints,
                                  rhs         = constraints[[con]]$rhs,
                                  neq         = constraints[[con]]$neq), ldots)
        conList[[con]] <- do.call("con_gorica_est", CALL.restr) 
      }
      names(conList) <- names(constraints)
      
      isSummary <- lapply(conList, function(x) summary(x, 
                                                       type        = type,
                                                       sample.nobs = sample.nobs)) 
    } else {
      stop("restriktor ERROR: I don't know how to handle an object of class ", paste0(class(object)[1]))
    }

  
  ## add objectnames if not available
  # constraints must be a list
  if (!is.list(constraints)) {
    constraints <- list(constraints)
  }
  
  
  objectnames <- names(constraints) 
  # constraints are inherited
  if ("restriktor" %in% object_class) {
    objectnames <- paste0("H", 1:length(object))
  } else if (any(is.null(names(constraints))) || all(names(constraints) == "")) {  
    objectnames <- paste0("H", 1:length(constraints))
  } 

  
  if (comparison == "complement" && length(conList) > 1L) {
    warning("Restriktor Warning: Only one hypothesis is allowed (for now) when comparison = 'complement'.",
            " Setting comparison to 'unconstrained' instead.", call. = FALSE)
    comparison <- "unconstrained"
  } 

  

# compute complement ------------------------------------------------------
  df.c <- NULL
  if (comparison == "complement") {
      # unrestricted estimates
      if (inherits(object$object, "numeric")) {
        b.unrestr <- object$object
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
      # remove all zero rows
      Amat <- Amat[apply(Amat, 1, function(x) !all(x == 0)), , drop = FALSE]
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
    
    if (!is.null(auto.bound) && meq == 0L) {
      warning("restriktor Warning: auto.bounds are only available for equality restrictions \n",
              " and are therefore ignored.", call. = FALSE)
      auto.bound <- NULL
    } 
    
# about equality auto.bounds ---------------------------------------------------
#     if (!is.null(auto.bound)) {
# 
#       betasc <- NA
#       
#       if (auto.bound <= 0L) {
#         stop("Restriktor ERROR: the auto.bounds must be positive!")
#       }
#       
#       # upper-bound = rhs + bound
#       ub <- bvec.ceq + bound
#       # lower-bound = rhs - bound
#       lb <- bvec.ceq - bound
#       
#       nr.meq <- 1:meq
#       
#       ## If no constraints are violated, then the complement is on the boundary
#       ## of the bounds.
#       
#       Amatx <- rbind(Amat.ceq, -Amat.ceq, Amat.ciq)                  
#       bvecx <- c(lb, -ub, bvec.ciq)                                  
#       
#       # constraints are in agreement with the unconstrained estimates
#       if (all(Amatx %*% b.unrestr - bvecx >= 0)) {
#         if (type == "goric") { 
#           llm <- logLik(ans$model.org)
#         } else if (type == "gorica") {
#           llm <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE) 
#         }
#         # check all bounds
#         
#         llc <- list()
#         for (i in 1:meq) {
#          bvecx2 <- bvecx
#          
#          ## this is the same as specifying a range restriction, e.g., x1 < 0.05; x1 > -0.05
#          bvecx2[i] <- ub[i]
#          Hc.ub <- restriktor(ans$model.org, constraints = Amatx,
#                              neq = i, rhs = bvecx2,
#                              mix.weights = "none", se = "none")
#          
#          bvecx2[i] <- lb[i]
#          Hc.lb <- restriktor(ans$model.org, constraints = Amatx,
#                              neq = i, rhs = bvecx2,
#                              mix.weights = "none", se = "none")
#          
#          bvecx2[i + meq] <- ub[i]
#          Hc.ub2 <- restriktor(ans$model.org, constraints = Amatx,
#                               neq = i, rhs = bvecx2,
#                               mix.weights = "none", se = "none")
#          
#          bvecx2[i + meq] <- lb[i]
#          Hc.lb2 <- restriktor(ans$model.org, constraints = Amatx,
#                               neq = i, rhs = bvecx2,
#                               mix.weights = "none", se = "none")
#          
#          
#          llc[[i]] <- c(logLik(Hc.ub), logLik(Hc.lb), logLik(Hc.ub2), logLik(Hc.lb2))
#          
#          # hoe doet Muthen het met about equality
#          # hoe zit het nu met de true hypothesis rate
#         }
#        
#         llc <- max(unlist(llc))
#       } else {
#         ## If one constraint is violated, then the complement is equal to the 
#         ## unconstrained model
#         if (type == "goric") { 
#           llc <- logLik(ans$model.org)
#         } else if (type == "gorica") {
#           llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE) 
#         } 
#       }
# 
# #################################### no bounds ################################    
#     } else 
#       
    if (is.null(auto.bound)) {
      # check if any equality constraint is violated
      check.ceq <- !(all(c(Amat.ceq %*% c(b.unrestr)) - bvec.ceq == 0))
      if (nrow(Amat) > meq) {
        # check if any inequality constraint is violated
        check.ciq <- !(all(c(Amat.ciq %*% c(b.unrestr)) - bvec.ciq >= 0))
      } else {
        check.ciq <- FALSE
      }
      # compute log-likelihood for complement
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
            #ldots4 <- ldots[m.restr > 0L]
            ldots$mix.weights <- "none"
            CALL.restr <- append(list(object      = b.unrestr,
                                      constraints = Amatx,
                                      rhs         = bvec[idx],
                                      neq         = 1,
                                      VCOV        = VCOV),
                                 ldots)
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
        stop("restriktor ERROR: no complement exists for an unconstrained hypothesis.")
      } else {
        stop("restriktor ERROR: you might have found a bug, please contact me at: info@restriktor.org!")
      }
      if (type %in% c("goric", "goricc")) {
        llm <- logLik(conList[[1]])
      } else if (type %in% c("gorica", "goricac")) {
        llm <- dmvnorm(c(b.unrestr - b.restr), sigma = VCOV, log = TRUE)
      }
    } 
    
    # compute the number of free parameters f in the complement
    p <- ncol(VCOV)
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
      # if (inherits(object$object, "numeric")) {
      #   llm <- unlist(lapply(conList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
      #                                                       sigma = VCOV, log = TRUE) ))  
      #   # unrestricted
      #   llu <- dmvnorm(rep(0, ncol(VCOV)), sigma = VCOV, log = TRUE) 
      # } else { 
        llm <- unlist(lapply(conList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
                                                            sigma = VCOV, log = TRUE) )) 
        # unrestricted
        llu <- dmvnorm(rep(0, ncol(VCOV)), sigma = VCOV, log = TRUE) 
      # }
    }
    
    if (type %in% c("goric", "gorica")) {
      PTu <- 1 + ncol(VCOV)
    } else if (type %in% c("goricc", "goricac")) {
      if (is.null(sample.nobs)) {
        stop("restriktor ERROR: if type = \'goric(a)c\' the argument \'sample.nobs\' needs to be provided.",
             call. = FALSE)
      }
      N <- sample.nobs
      # unconstrained penalty
      PTu <- ( (N * (ncol(VCOV) + 1) / (N - ncol(VCOV) - 2) ) ) 
    }
    
    ## correct PT for gorica(c)
    if (type %in% c("gorica", "goricac")) {
      PTu <- PTu - 1
    }
    
    if (is.null(PTm)) {
      stop("restriktor ERROR: no chi-bar-square weights are found. Use mix.weights = 'pmvnorm' (default) or 'boot'.", call. = FALSE)
    }
    
    goric.Hm <- -2*(llm - PTm)
    goric.Hu <- -2*(llu - PTu)
    df.Hm <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                        goric = goric.Hm)
    df.Hm$model <- as.character(df.Hm$model)
    df.u <- data.frame(model = "unconstrained", loglik = llu, penalty = PTu, 
                       goric = goric.Hu)
    df <- rbind(df.Hm, df.u)
    rownames(df) <- NULL
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
    }
    df <- rbind(df.Hm, df.c)
    names(df)[4] <- type
  } else {
    stop("restriktor ERROR: I don't know how to compute the goric-values.")
  }

  LL <- -2*df$loglik
  delta_LL <- LL - min(LL)
  loglik.weights <- exp(0.5 * -delta_LL) / sum(exp(0.5 * -delta_LL))
  penalty.weights <- exp(-df$penalty) / sum(exp(-df$penalty))
  df$loglik.weights <- loglik.weights
  df$penalty.weights <- penalty.weights
  
  ans$objectList  <- conList
  ans$objectNames <- objectnames
  # compute goric weights and ratio weights
  delta <- df$goric - min(df$goric)
  goric.weights <- exp(0.5 * -delta) / sum(exp(0.5 * -delta))
  df$goric.weights <- goric.weights
  names(df)[7] <- paste0(type, ".weights")
  rownames(df) <- NULL
  ans$result <- df

  # compute ratio weights
  modelnames <- as.character(df$model)
  if (length(modelnames) > 1) {
    goric_rw   <- goric.weights %*% t(1/goric.weights)
    penalty_rw <- penalty.weights %*% t(1/penalty.weights)
    loglik_rw  <- loglik.weights %*% t(1/loglik.weights)
    # it might happen that a diagonal value results in NaN.
    diag(goric_rw) <- 1
    diag(penalty_rw) <- 1
    diag(loglik_rw) <- 1
    rownames(goric_rw) <- modelnames
    rownames(penalty_rw) <- modelnames
    rownames(loglik_rw) <- modelnames
    colnames(goric_rw) <- paste0("vs. ", modelnames)
    colnames(penalty_rw) <- paste0("vs. ", modelnames)
    colnames(loglik_rw) <- paste0("vs. ", modelnames)
    ans$ratio.gw <- goric_rw
    ans$ratio.pw <- penalty_rw
    ans$ratio.lw <- loglik_rw
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
  ans$messages$mix_weights <- do.call("rbind", lapply(isSummary, FUN = function(x) { x$messages$mix_weights }))
  ans$messages$mix_weights <- ans$messages$mix_weights [!duplicated(ans$messages$mix_weights )]
  
  if (type %in% c("goric", "goricc")) {
    class(ans) <- "con_goric"
  } else if (type %in% c("gorica", "goricac")) {
    class(ans) <- c("con_gorica", "con_goric")
  }
  
  ans
}

# -------------------------------------------------------------------------



# object of class lm ------------------------------------------------------
goric.lm <- function(object, ..., hypotheses = NULL,
                     comparison = "unconstrained",
                     type = "goric",
                     missing = "none", auxiliary = c(), emControl = list(),
                     auto.bound = NULL, debug = FALSE) {
  
  if (!inherits(object, "lm")) {
    stop("restriktor ERROR: the object must be of class lm, glm, mlm, rlm.")
  }
  
  if (is.null(hypotheses)) {
   stop("restriktor ERROR: The 'hypotheses' argument is missing. Please make sure 
        to provide a valid set of hypotheses, for example, hypotheses = list(h1 = 'x1 > x2 > x3').", call. = FALSE) 
  }
  
  if (!is.list(hypotheses)) {
    stop("restriktor ERROR: the hypotheses must be specified as a list. 
         For example, hypotheses = list(h1 = 'x1 > x2 > x3')", call. = FALSE)
  }
  
  if (missing %in% c("em", "EM", "two.stage", "twostage")) {
    missing <- "two.stage" 
  } 
  
  objectList <- list(...)
  
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  
  # only one object of class lm is allowed
  isLm <- unlist(lapply(objectList, function(x) class(x)[1] %in% c("lm", "glm", "mlm", "rlm")))
  if (sum(isLm) > 1L) {
    stop("restriktor ERROR: multiple objects of class lm found, only 1 is allowed.") 
  }
  
  if (missing == "two.stage") {
    if (family(object)$family != "gaussian") {
      stop("Restriktor ERROR: \"two.stage\" is not available in the categorical setting")
    }
    tsm  <- two_stage_matrices(object, auxiliary = auxiliary, emControl = emControl)
    vcov <- two_stage_sandwich(tsm)
    est <- coef(tsm$fitTarget)
    N <- tsm$N
    
    if (!type %in% c("gorica", "goricac")) {
      stop("restriktor EROR: missing = \"two.stage\" is only (for now) available for type = 'gorica(c)'\n", call. = FALSE)
    }
    objectList$object <- est
    objectList$VCOV   <- vcov
  } else {
    objectList$object <- object
    objectList$VCOV        <- NULL
  } 
    
    objectList$hypotheses  <- hypotheses
    objectList$comparison  <- comparison
    objectList$type        <- type
    objectList$auto.bound  <- auto.bound
    objectList$debug       <- debug
    objectList$missing     <- NULL
    objectList$auxiliary   <- NULL
    objectList$emControl   <- NULL
    
    if (type == "goricac") {
      objectList$sample.nobs <- length(residuals(object))
    }
    
    objectList <- sapply(objectList, function(x) eval(x))
    
    if (missing == "none") {
      object_idx <- grepl("object", names(objectList))
      objectList <- append(list(object = objectList[object_idx]), objectList[!object_idx])
    }
    
    # which arguments are allowed
    arguments <- c("B", "mix.weights", "mix.bootstrap", "parallel", "ncpus", 
                   "cl", "seed", "control", "verbose", "debug", "comparison",
                   "type", "auto.bound", "hypotheses", "missing", "auxiliary",
                   "VCOV", "sample.nobs", "object")
    
    # check for unkown arguments
    pm <- pmatch(names(objectList), arguments, 0L)
    
    if (any(pm == 0)) {
      stop("restriktor ERROR: argument ", sQuote(names(objectList)[pm == 0]), " unknown.", call. = FALSE)
    }
    
    if (missing == "two.stage") {
      res <- do.call(goric.numeric, objectList)
    } else {    
      #res <- do.call(goric.default, c(objectList[object_idx], objectList[!object_idx]))
      res <- do.call(goric.default, objectList)
    }
    # res <- do.call(goric.numeric, args = list(object = est, VCOV = vcov, sample.nobs = N,
    #                                           constraints = constraints,
    #                                           comparison = comparison, type = type))    
  
  res
}



# object of class restriktor ----------------------------------------------
goric.restriktor <- function(object, ..., hypotheses = NULL,
                             comparison = "unconstrained",
                             type = "goric",
                             auto.bound = NULL, debug = FALSE) {
  
  # hypotheses are inherited from the restriktor object
  
  if (!inherits(object, "restriktor")) {
    stop("restriktor ERROR: the object must be of class restriktor.")
  }
  
  objectList <- list(...)
  
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  
  objectList$hypotheses  <- hypotheses
  objectList$comparison  <- comparison
  objectList$type        <- type
  objectList$auto.bound  <- auto.bound
  objectList$debug       <- debug
  objectList$VCOV        <- NULL
  
  if (type == "goricac") {
    objectList$sample.nobs <- length(residuals(object))
  }
  
  objectList <- sapply(objectList, function(x) eval(x))
  
  # multiple objects of class restriktor are allowed
  isRestr <- unlist(lapply(objectList, function(x) class(x)[1] == "restriktor"))

  # put all objects of class restriktor in one list
  objectList <- append(list(object = objectList[isRestr]), objectList[!isRestr])
  
  # which arguments are allowed
  arguments <- c("B", "mix.weights", "mix.bootstrap", "parallel", "ncpus", 
                 "cl", "seed", "control", "verbose", "debug", "comparison",
                 "type", "auto.bound", "hypotheses", "missing", "auxiliary",
                 "VCOV", "sample.nobs", "object")
  
  # check for unkown arguments
  pm <- pmatch(names(objectList), arguments, 0L)
  if (any(pm == 0)) {
    stop("restriktor ERROR: argument ", sQuote(names(objectList)[pm == 0]), " unknown.", call. = FALSE)
  }
  
  res <- do.call(goric.default, c(objectList[isRestr], objectList[!isRestr]))
  
  res
}



# object of class numeric -------------------------------------------------
goric.numeric <- function(object, ..., hypotheses = NULL,
                          VCOV = NULL,
                          comparison = "unconstrained",
                          type = "gorica", sample.nobs = NULL,
                          auto.bound = NULL, debug = FALSE) {
  
  if (!inherits(object, "numeric")) {
    stop("restriktor ERROR: the object must be of class numeric.")
  }
  
  if (is.null(hypotheses)) {
    stop("restriktor ERROR: The 'hypotheses' argument is missing. Please make sure 
         to provide a valid set of hypotheses, for example, hypotheses = list(h1 = 'x1 > x2 > x3').", call. = FALSE)   }
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop("restriktor ERROR: object of class numeric is only supported for type = 'gorica(c)'.")
  }
  
  if (is.null(VCOV)) {
    stop("restriktor ERROR: the argument VCOV is not found.")
  } else {
    # check if it is a matrix
    if (!is.matrix(VCOV)) {
      # used in lme4
      if (inherits(VCOV, "dpoMatrix")) {
        VCOV <- as.matrix(VCOV)
      }
    }
  }
  
  if (is.null(sample.nobs) && type %in% c("goricac")) {
    stop("restriktor ERROR: the argument sample.nobs is not found.")
  }
  
  if (!is.list(hypotheses)) {
    stop("restriktor ERROR: the hypotheses must be specified as a list. 
         For example, hypotheses = list(h1 = 'x1 > x2 > x3')", call. = FALSE)
  }
  
  
  objectList <- list(...)
  
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  
  objectList$object      <- object
  objectList$hypotheses  <- hypotheses
  objectList$comparison  <- comparison
  objectList$type        <- type
  objectList$auto.bound  <- auto.bound
  objectList$debug       <- debug
  objectList$VCOV        <- VCOV
  objectList$sample.nobs <- sample.nobs
  
  object_idx <- grepl("object", names(objectList))
  objectList <- append(list(object = objectList[object_idx]), objectList[!object_idx])
  
  # which arguments are allowed
  arguments <- c("B", "mix.weights", "mix.bootstrap", "parallel", "ncpus", 
                 "cl", "seed", "control", "verbose", "debug", "comparison",
                 "type", "auto.bound", "hypotheses", "missing", "auxiliary",
                 "VCOV", "sample.nobs", "object")
  
  # check for unkown arguments
  pm <- pmatch(names(objectList), arguments, 0L)
  if (any(pm == 0)) {
    stop("restriktor ERROR: argument ", sQuote(names(objectList)[pm == 0]), " unknown.", call. = FALSE)
  }
  
  res <- do.call(goric.default, c(objectList[object_idx], objectList[!object_idx]))
  
  
  res
}


# object of class lavaan --------------------------------------------------
goric.lavaan <- function(object, ..., hypotheses = NULL,
                         comparison = "unconstrained",
                         type = "gorica",
                         standardized = FALSE,
                         auto.bound = NULL, debug = FALSE) {
  
  if (!inherits(object, "lavaan")) {
    stop("restriktor ERROR: the object must be of class lavaan.")
  }
  
  if (is.null(hypotheses)) {
    stop("restriktor ERROR: The 'hypotheses' argument is missing. Please make sure 
         to provide a valid set of hypotheses, for example, hypotheses = list(h1 = 'x1 > x2 > x3').", call. = FALSE) 
  }
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop("restriktor ERROR: object of class lavaan is only supported for type = 'gorica(c)'.")
  }
  
  
  if (!is.list(hypotheses)) {
    stop("restriktor ERROR: the hypotheses must be specified as a list. 
         For example, hypotheses = list(h1 = 'x1 > x2 > x3')", call. = FALSE)
  }
  
  objectList <- list(...)
  
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  mcList$standardized <- NULL
  
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  
  est <- con_gorica_est_lav(object, standardized)
  objectList$object     <- est$estimate
  objectList$VCOV       <- est$VCOV
  objectList$comparison <- comparison
  objectList$type       <- type
  objectList$auto.bound <- auto.bound
  objectList$debug      <- debug
  
  
  if (type == "goricac") {
    objectList$sample.nobs <- lavInspect(object, what = "ntotal")
  }
  
  object_idx <- grepl("object", names(objectList))
  objectList <- append(list(object = objectList[object_idx]), objectList[!object_idx])
  
  
  # which arguments are allowed
  arguments <- c("B", "mix.weights", "mix.bootstrap", "parallel", "ncpus", 
                 "cl", "seed", "control", "verbose", "debug", "comparison",
                 "type", "auto.bound", "hypotheses", "missing", "auxiliary",
                 "VCOV", "sample.nobs", "object")
  
  # check for unkown arguments
  pm <- pmatch(names(objectList), arguments, 0L)
  if (any(pm == 0)) {
    stop("restriktor ERROR: argument ", sQuote(names(objectList)[pm == 0]), " unknown.", call. = FALSE)
  }
  
  res <- do.call(goric.default, c(objectList[object_idx], objectList[!object_idx]))
  
  res
}


# object of class CTmeta --------------------------------------------------
goric.CTmeta <- function(object, ..., hypotheses = NULL,
                         comparison = "unconstrained",
                         type = "gorica", sample.nobs = NULL,
                         auto.bound = NULL, debug = FALSE) {
  
  if (!inherits(object, "CTmeta")) {
    stop("restriktor ERROR: the object must be of class CTmeta.")
  }
  
  if (is.null(hypotheses)) {
    stop("restriktor ERROR: The 'hypotheses' argument is missing. Please make sure to 
         provide a valid set of hypotheses, for example, hypotheses = list(h1 = 'x1 > x2 > x3').", call. = FALSE) 
  }
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop("restriktor ERROR: object of class CTmeta is only supported for type = 'gorica(c)'.")
  }
  
  if (!is.list(hypotheses)) {
    stop("restriktor ERROR: the hypotheses must be specified as a list. 
         For example, hypotheses = list(h1 = 'x1 > x2 > x3')", call. = FALSE)
  }
  
  objectList <- list(...)
  
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList
  objectList$object <- NULL
  
  objectList$object <- coef(object) 
  objectList$VCOV   <- vcov(object)
  objectList$hypotheses  <- hypotheses
  objectList$comparison  <- comparison
  objectList$type        <- type
  objectList$auto.bound  <- auto.bound
  objectList$debug       <- debug
  
  
  if (type == "goricac") {
    objectList$sample.nobs <- sample.nobs
  }
  
  object_idx <- grepl("object", names(objectList))
  objectList <- append(list(object = objectList[object_idx]), objectList[!object_idx])  
  
  # which arguments are allowed
  arguments <- c("B", "mix.weights", "mix.bootstrap", "parallel", "ncpus", 
                 "cl", "seed", "control", "verbose", "debug", "comparison",
                 "type", "auto.bound", "hypotheses", "missing", "auxiliary",
                 "VCOV", "sample.nobs", "object")
  
  # check for unkown arguments
  pm <- pmatch(names(objectList), arguments, 0L)
  if (any(pm == 0)) {
    stop("restriktor ERROR: argument ", sQuote(names(objectList)[pm == 0]), " unknown.", call. = FALSE)
  }
  
  res <- do.call(goric.default, c(objectList[object_idx], objectList[!object_idx]))
  
  res
}




# object of class rma -----------------------------------------------------
goric.rma <- function(object, ..., hypotheses = NULL,
                      VCOV = NULL,
                      comparison = "unconstrained",
                      type = "gorica", sample.nobs = NULL,
                      auto.bound = NULL, debug = FALSE) {
  
  if (!inherits(object, c("rma"))) {
    stop("restriktor ERROR: the object must be of class lm, glm, mlm, rlm, rma (only 'rma.uni').")
  }
  
  if (is.null(hypotheses)) {
    stop("restriktor ERROR: The 'hypotheses' argument is missing. Please make sure to provide a valid set of hypotheses, for example, hypotheses = list(h1 = 'x1 > x2 > x3').", call. = FALSE) 
  }
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop("restriktor ERROR: object of class rma is only supported for type = 'gorica(c)'.")
  }
  
  if (!is.list(hypotheses)) {
    stop("restriktor ERROR: the hypotheses must be specified as a list. \nFor example, hypotheses = list(h1 = 'x1 > x2 > x3')", call. = FALSE)
  }
  
  objectList <- list(...)
  
  mcList <- as.list(match.call())
  mcList <- mcList[-c(1)]
  
  mcnames <- names(mcList) == ""
  lnames <- as.character(mcList[mcnames])
  names(mcList)[mcnames] <- lnames
  objectList <- mcList  
  objectList$object <- NULL
  
  objectList$object <- coef(object) 
  objectList$VCOV   <- vcov(object)
  objectList$comparison   <- comparison
  objectList$type         <- type
  objectList$auto.bound   <- auto.bound
  objectList$debug        <- debug
  
  if (type == "goricac") {
    objectList$sample.nobs <- sample.nobs
  }
  
  object_idx <- grepl("object", names(objectList))
  objectList <- append(list(object = objectList[object_idx]), objectList[!object_idx])  
  
  # which arguments are allowed
  arguments <- c("B", "mix.weights", "mix.bootstrap", "parallel", "ncpus", 
                 "cl", "seed", "control", "verbose", "debug", "comparison",
                 "type", "auto.auto.bound", "hypotheses", "missing", "auxiliary",
                 "VCOV", "sample.nobs", "object")
  
  # check for unkown arguments
  pm <- pmatch(names(objectList), arguments, 0L)
  if (any(pm == 0)) {
    stop("restriktor ERROR: argument ", sQuote(names(objectList)[pm == 0]), " unknown.", call. = FALSE)
  }
  
  res <- do.call(goric.default, c(objectList[object_idx], objectList[!object_idx]))
  
  res
}





print.con_goric <- function(x, digits = max(3, getOption("digits") - 4), ...) {

  type <- x$type
  comparison <- x$comparison
  dig <- paste0("%6.", digits, "f")
  x2 <- lapply(x$result[-1], sprintf, fmt = dig)
  df <- data.frame(model = x$result$model, x2)

  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))

  if (type == "goric") {
    cat("generalized order-restricted information criterion: \n")
  } else if (type == "gorica") {
    cat("generalized order-restricted information criterion approximation:\n")
  } else if (type == "goricc") {
    cat("small sample generalized order-restricted information criterion:\n")
  } else if (type == "goricac") {
    cat("small sample generalized order-restricted information criterion approximation:\n")
  }

  wt_bar <- sapply(x$objectList, function(x) attr(x$wt.bar, "method") == "boot")
  # we need a check if the hypothesis is equalities only
  ceq_only <- sapply(x$objectList, function(x) nrow(x$constraints) == x$neq)
  wt_bar <- as.logical(wt_bar * !ceq_only)
  
  if (sum(wt_bar) > 0) {
    wt_method_boot <- x$objectList[wt_bar]
    wt_bootstrap_draws  <- sapply(wt_method_boot, function(x) attr(x$wt.bar, "mix.bootstrap"))
    wt_bootstrap_errors <- lapply(wt_method_boot, function(x) attr(x$wt.bar, "error.idx"))
    max_nchar <- max(nchar(names(wt_method_boot)))
    
    len <- length(wt_method_boot)
    if (len > 0) { 
      cat("\n")
      cat("Level probabilities:\n")
      cat("  Number of requested bootstrap draws", wt_bootstrap_draws[1], "\n")
      for (i in 1:len) {
        #cat("Number of successful bootstrap draws for", names(wt_method_boot)[i], ":", (wt_bootstrap_draws[1] - length(wt_bootstrap_errors[[i]])), "\n")
        cat(paste0("  Number of successful bootstrap draws for ", sprintf(paste0("%", max_nchar, "s"), names(wt_method_boot)[i]), ": ", (wt_bootstrap_draws[1] - length(wt_bootstrap_errors[[i]])), "\n"))
      }
    }
  }
    
  cat("\nResults:\n")
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  
  if (comparison == "complement") {
    objectnames <- as.character(df$model)
    class(x$ratio.gw) <- "numeric"
    cat("---", "\nThe order-restricted hypothesis", sQuote(objectnames[1]), "has", sprintf("%.3f", x$ratio.gw[1,2]), "times more support than its complement.\n\n")
  } else {
    cat("---\n")  
  }
  
  message(x$messages$mix_weights)
  
  invisible(x)
}




summary.con_goric <- function(object, brief = TRUE, 
                          digits = max(3, getOption("digits") - 4), ...) {
  
  x <- object
  type <- x$type
  comparison <- x$comparison
  dig <- paste0("%6.", digits, "f")
  
  ratio.gw <- x$ratio.gw
  rownames(ratio.gw) <- rownames(x$ratio.gw)
  
  x2 <- lapply(x$result[-1], sprintf, fmt = dig)
  df <- data.frame(model = x$result$model, x2)
  
  objectnames <- as.character(df$model)
  
  Amat <- x$constraints
  meq  <- x$neq
  bvec <- x$rhs
  iact <- lapply(x$objectList, FUN = function(x) { x$iact } )

  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  
  if (type == "goric") {
    cat("generalized order-restricted information criterion: \n")
  } else if (type == "gorica") {
    cat("generalized order-restricted information criterion approximation:\n")
  } else if (type == "goricc") {
    cat("small sample generalized order-restricted information criterion:\n")
  } else if (type == "goricac") {
    cat("small sample generalized order-restricted information criterion approximation:\n")
  }
  
  wt_bar <- sapply(x$objectList, function(x) attr(x$wt.bar, "method") == "boot")
  # we need a check if the hypothesis is equalities only
  ceq_only <- sapply(x$objectList, function(x) nrow(x$constraints) == x$neq)
  wt_bar <- as.logical(wt_bar * !ceq_only)
  
  if (sum(wt_bar) > 0) {
    wt_method_boot <- x$objectList[wt_bar]
    wt_bootstrap_draws  <- sapply(wt_method_boot, function(x) attr(x$wt.bar, "mix.bootstrap"))
    wt_bootstrap_errors <- lapply(wt_method_boot, function(x) attr(x$wt.bar, "error.idx"))
    max_nchar <- max(nchar(names(wt_method_boot)))
    
    len <- length(wt_method_boot)
    if (len > 0) { 
      cat("\n")
      cat("Level probabilities:\n")
      cat("  Number of requested bootstrap draws", wt_bootstrap_draws[1], "\n")
      for (i in 1:len) {
        #cat("Number of successful bootstrap draws for", names(wt_method_boot)[i], ":", (wt_bootstrap_draws[1] - length(wt_bootstrap_errors[[i]])), "\n")
        cat(paste0("  Number of successful bootstrap draws for ", sprintf(paste0("%", max_nchar, "s"), names(wt_method_boot)[i]), ": ", (wt_bootstrap_draws[1] - length(wt_bootstrap_errors[[i]])), "\n"))
      }
    }
  }
  
  cat("\nResults:\n")  
  
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE, right = TRUE)
  cat("---\n")
  
  
  if (any(combn(as.numeric(df$loglik), 2, FUN = function(x) abs(diff(x))) <= 1e-5)) { 
    cat("Note: If log-likelihood values are equal or close to each other, the goric ratio weights are determined", 
        "only by the difference in penalty values. Please check the ratio penalty-weights.\n\n")
  }
  if (comparison == "complement") {
    cat("The order-restricted hypothesis", sQuote(objectnames[1]), "has", 
        sprintf("%.3f", as.numeric(ratio.gw[1,2])), "times more support than its complement.\n\n")
  } 
  
  if (!is.null(x$ratio.gw)) {
    if (type == "goric") {
      cat("\nRatio GORIC-weights:\n")
    } else if (type == "gorica") {
      cat("\nRatio GORICA-weights:\n")
    } else if (type == "goricc") {
      cat("\nRatio GORICC-weights:\n")
    } else if (type == "goricac") {
      cat("\nRatio GORICAC-weights:\n") 
    }

    ratio.gw <- apply(x$ratio.gw, 2, sprintf, fmt = dig)
    rownames(ratio.gw) <- rownames(x$ratio.gw)
    class(ratio.gw) <- "numeric"
    
    if (max(ratio.gw, na.rm = TRUE) >= 1e4) {
      print(format(ratio.gw, digits = digits, scientific = TRUE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = FALSE) 
    } else {
      print(format(ratio.gw, digits = digits, scientific = FALSE, trim = TRUE), 
            print.gap = 5, quote = FALSE, right = FALSE)
    }
    cat("---\n")
  }
  
  if (!is.null(x$ratio.pw)) {
    cat("\nRatio penalty-weights:\n")
    ratio.pw <- apply(x$ratio.pw, 2, sprintf, fmt = dig)
    rownames(ratio.pw) <- rownames(x$ratio.pw)
    class(ratio.pw) <- "numeric"
    
    if (max(ratio.pw, na.rm = TRUE) >= 1e4) {
      print(format(x$ratio.pw, digits = digits, scientific = TRUE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = FALSE)
    } else {
      print(format(ratio.pw, digits = digits, scientific = FALSE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = FALSE)
    }
    cat("---\n")
  }
  
  if (!is.null(x$ratio.lw)) {
    cat("\nRatio loglik-weights:\n")
    ratio.lw <- apply(x$ratio.lw, 2, sprintf, fmt = dig)
    rownames(ratio.lw) <- rownames(x$ratio.lw) 
    class(ratio.lw) <- "numeric"
    
    if (max(ratio.lw, na.rm = TRUE) >= 1e4) {
      print(format(ratio.lw, digits = digits, scientific = TRUE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = FALSE)
    } else {
      print(format(ratio.lw, digits = digits, scientific = FALSE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = FALSE)
    }
    cat("---\n")
  }
  

  if (!brief) {
    cat("\n\nOrder-restricted coefficients:\n")
    coefs <- trimws(apply(x$ormle$b.restr, 2, sprintf, fmt = dig))
    coefs[coefs == "NA"] <- ""
    rownames(coefs) <- rownames(x$ormle$b.restr)
    print(format(coefs, digits = digits, scientific = FALSE), print.gap = 2,
          quote = FALSE)
    cat("---\n")
    
    vnames <- names(x$ormle$b.restr)
    vnames_len <- length(x$objectList)
    first_na <- apply(x$ormle$b.restr[1:vnames_len, ], 1, function(x) { which(is.na(x))[1] }) -1
    first_na[is.na(first_na)] <- 0
    
    selected_names <- list()
    for (i in 1:length(first_na)) {
      if (first_na[i] == 0) {
        selected_names[[i]] <- vnames
      } else {
        selected_names[[i]] <- vnames[1:first_na[i]]
      }
    }
    
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
    for (i in 1:vnames_len) {
      conMat[[i]] <- fn(Amat = Amat[[i]], bvec = bvec[[i]], meq = meq[[i]], 
                        iact = iact[[i]], vnames = selected_names[[i]])  
    }
    names(conMat) <- x$objectNames
    
    if (comparison == "complement") {
     conMat$complement <- paste("not", x$objectNames) 
    }
    
    cat("\nRestriction matrices:\n")
    print(conMat, quote = FALSE, scientific = FALSE)
   
    invisible(x)
  } else {
    if (!is.null(object$hypotheses_usr)) {
      cat("\norder-restricted hypotheses:\n\n")
      hypotheses_usr <- object$hypotheses_usr
      #hypotheses_usr <- lapply(hypotheses_usr, function(x) unlist(strsplit(x, "\n")))
      # remove user specified paramters
      #hypotheses_usr <- lapply(hypotheses_usr, function(x) x[!grepl(":=", x)])
      
      for (i in 1:length(x$objectList)) {
        text <- gsub("(\\n\\s+)+", "\n", hypotheses_usr[[i]])
        cat(paste0(objectnames[i],":\n", trimws(gsub("\\h+", " ", text, perl = TRUE))), "\n\n")
      }
      
      #calculate max length of vectors
      #max_length <- max(sapply(hypotheses_usr, function(x) length(x)))
      # for (j in 1:length(hypotheses_usr)) {
      #   length(hypotheses_usr[[j]]) <- max_length
      # }
      # hypotheses_usr <- do.call(rbind, hypotheses_usr)
      # row.names(hypotheses_usr) <- paste0(objectnames[1:length(x$objectList)], ":")
      # hypotheses_usr[is.na(hypotheses_usr)] <- ""
      # name.width <- max(sapply(hypotheses_usr, nchar))
      # hypotheses_usr <- format(hypotheses_usr, width = name.width, justify = "left")
      # hypotheses_usr <- as.data.frame(hypotheses_usr)
      # names(hypotheses_usr) <- NULL
      # print(hypotheses_usr)
    }
  }
  cat("\n")
  message(x$messages$mix_weights)
}


