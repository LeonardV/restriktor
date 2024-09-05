goric <- function(object, ...) { UseMethod("goric") }


goric.default <- function(object, ..., hypotheses = NULL,
                          comparison = NULL, 
                          VCOV = NULL, sample_nobs = NULL,
                          type = "goric", control = list(),
                          debug = FALSE) {
  
  # the following classes are allowed (for now)
  obj_class <- class(object)
  classes <- c("aov", "lm", "glm", "mlm", "rlm", "numeric", "lavaan", "CTmeta", 
               "rma.uni", "nlmerMod", "glmerMod", "merMod")
  check_class <- obj_class %in% classes
  if (!any(check_class)) {
    stop(paste("Objects of class", paste(obj_class, collapse = ", "), "are not supported. Supported classes are:", paste(classes, collapse = ", "), "."))
  }
  
  if (is.null(hypotheses)) {
    stop(paste("Restriktor ERROR: The 'hypotheses' argument is missing. Please make sure",
         "to provide a valid set of hypotheses, for example, hypotheses =",
         "list(h1 = 'x1 > x2 > x3')."), call. = FALSE)
  } else {
    if (!is.list(hypotheses)) {
      stop(paste("Restriktor ERROR: the hypotheses must be specified as a list.",
      "For example, hypotheses = list(h1 = 'x1 > x2 > x3')", call. = FALSE))
    }
    
    if (length(hypotheses) == 1 && is.null(comparison)) {
      comparison <- "complement"
    }
    
  }
  
  if (is.null(sample_nobs) && type %in% c("goricac")) {
    stop(paste("Restriktor ERROR: the argument sample_nobs is not found."))
  }
  
  if (!is.null(VCOV)) {
    # check if it is of class matrix
    if (inherits(VCOV, "dpoMatrix")) {
      VCOV <- as.matrix(VCOV)
    }
  }
  
  ldots <- list(...)
  ldots$missing <- NULL
  ldots$control <- control
  
  # which arguments are allowed
  goric_arguments <- c("B", "mix_weights", "parallel", 
                       "ncpus", "cl", "seed", "control", "verbose", "debug", 
                       "comparison", "type", "hypotheses", "auxiliary",
                       "VCOV", "sample_nobs", "object",
                       # for rtmvnorm() function
                       "lower", "upper", "algorithm",
                       "burn.in.samples", "start.values", "thinning")
  
  # check for unkown arguments
  pm <- pmatch(names(ldots), goric_arguments, 0L)
  if (any(pm == 0)) {
    stop("Restriktor ERROR: argument ", 
         sQuote(paste(names(ldots)[pm == 0], collapse = ", ")), 
         " unknown.", call. = FALSE)
  }
  
  if (any(c("lower", "upper", "algorithm", "burn.in.samples", "start.values", 
            "thinning") %in% names(ldots))) {
    ldots$mix_weights <- "boot"
  }
  
  constraints <- hypotheses
  # class objects
  object_class <- obj_class

  # some checks
  if (!is.null(comparison)) {
    comparison <- tolower(comparison)
  }
  comparison <- match.arg(comparison, c("unconstrained", "complement", "none"))
  type <- tolower(type)
  type <- match.arg(type, c("goric", "goricc", "gorica", "goricac"))
  
  conChar <- sapply(constraints, function(x) inherits(x, "character"))
  isConChar <- all(conChar)

# -------------------------------------------------------------------------
  # create output list
  ans <- list()
  
  ## deal with objects of different classes
  # if ("restriktor" %in% object_class) {   
  #   # if all objects are of class restriktor
  #   conList   <- object
  #   isSummary <- lapply(conList, function(x) summary(x, 
  #                                                    goric       = type,
  #                                                    sample.nobs = sample_nobs))
  #   
  #   PT_Amat <- lapply(isSummary, function(x) x$PT_Amat)
  #   PT_meq  <- lapply(isSummary, function(x) x$PT_meq)
  #   
  #   for (lnames in names(conList)) {
  #     conList[[lnames]]$PT_Amat <- PT_Amat[[lnames]]
  #     conList[[lnames]]$PT_meq <- PT_meq[[lnames]]
  #   }
  #   
  #   ans$hypotheses_usr <- lapply(conList, function(x) x$CON$constraints)
  #   ans$model.org <- object[[1]]$model.org
  #   sample_nobs   <- nrow(model.frame(object[[1]]$model.org))
  #   # unrestricted VCOV
  #   VCOV <- vcov(ans$model.org)
  # } else 
  if (any(object_class %in% c("lm","rlm","glm","mlm")) && isConChar) { 
    # standard errors are not needed
    ldots$se <- "none"
    
    # the mixing-weights (LPs) are computed differently for object with class goric,
    # then for objects of class restriktor. In case of range restrictions 
    # (e.g., 0 < x < 1), the restrictions are treated as equality constraints for
    # computing PT. 
    # 
    class(object) <- append(class(object), "goric")
    
    # fit restriktor object for each hypothesis
    conList <- lapply(constraints, function(constraint) {
      CALL.restr <- append(list(object      = object, 
                                constraints = constraint), ldots)
      do.call("restriktor", CALL.restr)
    })
    names(conList) <- names(constraints)
    # compute summary for each restriktor object. 
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     goric       = type,
                                                     sample.nobs = sample_nobs))
    
    PT_Amat <- lapply(isSummary, function(x) x$PT_Amat)
    PT_meq  <- lapply(isSummary, function(x) x$PT_meq)
    
    for (lnames in names(conList)) {
      conList[[lnames]]$PT_Amat <- PT_Amat[[lnames]]
      conList[[lnames]]$PT_meq <- PT_meq[[lnames]]
    }
    
    ans$hypotheses_usr <- lapply(conList, function(x) x$CON$constraints)
    # add unrestricted object to output
    ans$model.org <- object
    # unrestricted VCOV
    VCOV <- vcov(ans$model.org)
    sample_nobs <- nrow(model.frame(object))
    idx <- length(conList) 
    objectnames <- vector("character", idx)
  } else if (any(object_class %in% c("aov", "lm","rlm","glm","mlm")) && !isConChar) {
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
    class(object) <- append(class(object), "goric")
    # fit restriktor object for each hypothesis
    conList <- lapply(constraints, function(constraint) {
      CALL.restr <- append(list(object      = object,
                                constraints = constraint$constraints,
                                rhs         = constraint$rhs,
                                neq         = constraint$neq), ldots)
      do.call("restriktor", CALL.restr)
    })
    
    names(conList) <- names(constraints)
    # compute symmary for each restriktor object. Here is the goric value 
    # computed. Note: not the gorica value
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     goric       = type,
                                                     sample.nobs = sample_nobs))
    
    PT_Amat <- lapply(isSummary, function(x) x$PT_Amat)
    PT_meq  <- lapply(isSummary, function(x) x$PT_meq)
    
    for (lnames in names(conList)) {
      conList[[lnames]]$PT_Amat <- PT_Amat[[lnames]]
      conList[[lnames]]$PT_meq <- PT_meq[[lnames]]
    }
    
    # add unrestricted object to output
    ans$model.org <- object
    # unrestricted VCOV
    VCOV <- vcov(ans$model.org) 
    sample_nobs <- nrow(model.frame(object))
    idx <- length(conList) 
    objectnames <- vector("character", idx)
    #CALL$object <- NULL
  } else if ("numeric" %in% object_class && isConChar) {
    # fit restriktor object for each hypothesis
    conList <- lapply(constraints, function(constraint) {
      CALL.restr <- append(list(object      = object, 
                                constraints = constraint,
                                VCOV        = as.matrix(VCOV)), ldots)
      do.call("con_gorica_est", CALL.restr)
    })
    
    names(conList) <- names(constraints)
    
    ans$hypotheses_usr <- lapply(conList, function(x) x$CON$constraints)
    
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     type        = type,
                                                     sample.nobs = sample_nobs)) 
    } else if ("numeric" %in% object_class && !isConChar) {
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
      conList <- lapply(constraints, function(constraint) {
        CALL.restr <- append(list(object      = object,
                                  VCOV        = as.matrix(VCOV),
                                  constraints = constraint$constraints,
                                  rhs         = constraint$rhs,
                                  neq         = constraint$neq), ldots)
        do.call("con_gorica_est", CALL.restr)
      })
      names(conList) <- names(constraints)
      
      isSummary <- lapply(conList, function(x) summary(x, 
                                                       type        = type,
                                                       sample.nobs = sample_nobs)) 
    } else {
      stop("Restriktor ERROR: I don't know how to handle an object of class ", paste0(class(object)[1]))
    }

  
  ## add objectnames if not available
  # constraints must be a list
  if (!is.list(constraints)) {
    constraints <- list(constraints)
  }
  
  objectnames <- names(constraints) 
  # constraints are inherited
  if (any(is.null(names(constraints)))) {  
    objectnames <- paste0("H", seq_len(length(constraints)))
  } else if (any(names(constraints) == "")) {
    objectnames <- ifelse(names(constraints) == "", paste0("H", seq_along(names(constraints))), 
                          names(constraints))
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
      if (inherits(object, "numeric")) {
        b.unrestr <- object
      } else {
        b.unrestr <- coef(ans$model.org)
      }
      # restricted estimates
      b.restr <- conList[[1]]$b.restr
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
      Amat.ceq <- matrix(numeric(0), nrow = 0, ncol = ncol(Amat))
      bvec.ceq <- rep(0, 0)
      Amat.ciq <- Amat[, , drop = FALSE]
      bvec.ciq <- bvec
    }
    
    # compute log-likelihood for complement
    # moet dit obv PT_Amat en PT_meq?
    LL_c <- compute_complement_likelihood(ans$model.org, VCOV, 
                                          Amat, Amat.ciq, Amat.ceq, 
                                          bvec, bvec.ciq, bvec.ceq, 
                                          meq, b.unrestr, type, ldots,
                                          debug = debug)
    llc <- LL_c$llc
    betasc <- LL_c$betasc
    
    # compute log-likelihood model
    if (type %in% c("goric", "goricc")) {
      llm <- logLik(conList[[1]])
    } else if (type %in% c("gorica", "goricac")) {
      llm <- dmvnorm(c(b.unrestr - b.restr), sigma = VCOV, log = TRUE)
    }

    # compute complement penalty term value 
    PTc <- penalty_complement_goric(Amat = conList[[1]]$PT_Amat, 
                                    meq  = conList[[1]]$PT_meq, 
                                    type, wt.bar, 
                                    debug = debug, 
                                    sample.nobs = sample_nobs)   
  } 
   
  ## for complement compute loglik-value, goric(a)-values, and PT-values if comparison = unconstrained
  switch(comparison,
          "unconstrained" = { 
            PTm <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
            if (type %in% c("goric", "goricc")) {
              # model
              llm <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
              # unrestricted
              llu <- logLik(ans$model.org)
            } else if (type %in% c("gorica", "goricac")) {
                llm <- unlist(lapply(conList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
                                                                    sigma = VCOV, log = TRUE) )) 
                # unrestricted
                llu <- dmvnorm(rep(0, ncol(VCOV)), sigma = VCOV, log = TRUE) 
            }
            
            if (type %in% c("goric", "gorica")) {
              PTu <- 1 + ncol(VCOV)
            } else if (type %in% c("goricc", "goricac")) {
              if (is.null(sample_nobs)) {
                stop("Restriktor ERROR: if type = \'goric(a)c\' the argument \'sample_nobs\' needs to be provided.",
                     call. = FALSE)
              }
              N <- sample_nobs
              # unconstrained penalty
              PTu <- ( (N * (ncol(VCOV) + 1) / (N - ncol(VCOV) - 2) ) ) 
            }
            
            ## correct PT for gorica(c)
            if (type %in% c("gorica", "goricac")) {
              PTu <- PTu - 1
            }
            
            if (is.null(PTm)) {
              stop("Restriktor ERROR: no chi-bar-square weights (a.k.a. level probabilities) are found. Use mix_weights = 'pmvnorm' (default) or 'boot'.", call. = FALSE)
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
           },
          "complement" = {
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
            },
         "none" = {
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
        },
          stop("Restriktor ERROR: I don't know how to compute the goric-values.")
  )

  ans$objectList  <- conList
  ans$objectNames <- objectnames
  
  # calculate LL, PT, goric weights and ratios
  model_comparison_metrics <- calculate_model_comparison_metrics(df)
  
  df$loglik.weights  <- model_comparison_metrics$loglik_weights
  df$penalty.weights <- model_comparison_metrics$penalty_weights
  df$goric.weights   <- model_comparison_metrics$goric_weights
  df$goric.weights_without_unc <- model_comparison_metrics$goric_weights_without_unc
  names(df)[7] <- paste0(type, ".weights")
  if (!is.null(df$goric.weights_without_unc)) {
    names(df)[8] <- paste0(type, ".weights_without_unc")
  }
  rownames(df) <- NULL

  ans$result <- df
  ans$ratio.gw <- model_comparison_metrics$goric_rw
  ans$ratio.pw <- model_comparison_metrics$penalty_rw
  ans$ratio.lw <- model_comparison_metrics$loglik_rw

  
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
  
  # Extracting and naming in one step for each attribute
  attributes <- c("constraints", "rhs", "neq")
  for (attr in attributes) {
    extracted <- lapply(conList, function(x) x[[attr]])
    names(extracted) <- ans$objectNames
    ans[[attr]] <- extracted
  }
  
  ans$Sigma <- VCOV
  ans$b.unrestr <- conList[[1]]$b.unrestr
  ans$ormle$b.restr <- coefs  
  ans$comparison <- comparison
  ans$type <- type

  # Assign class based on type\
  classMappings <- list(
    "goric"   = "con_goric",
    "goricc"  = "con_goric",
    "gorica"  = c("con_gorica", "con_goric"),
    "goricac" = c("con_gorica", "con_goric")
  )
  
  class(ans) <- classMappings[[type]]
  
  return(ans)
}



# object of class restriktor ----------------------------------------------
# goric.restriktor <- function(object, ..., hypotheses = NULL,
#                              comparison = "unconstrained",
#                              type = "goric",
#                              debug = FALSE) {
#   
#   # hypotheses are inherited from the restriktor object
#   
#   if (!inherits(object, "restriktor")) {
#     stop("Restriktor ERROR: the object must be of class restriktor.")
#   }
#   
#   objectList <- list(...)
#   
#   mcList <- as.list(match.call())
#   mcList <- mcList[-c(1)]
#   
#   mcnames <- names(mcList) == ""
#   lnames <- as.character(mcList[mcnames])
#   names(mcList)[mcnames] <- lnames
#   objectList <- mcList  
#   
#   objectList$hypotheses  <- hypotheses
#   objectList$comparison  <- comparison
#   objectList$type        <- type
#   objectList$debug       <- debug
#   objectList$VCOV        <- NULL
#   
#   if (type == "goricac") {
#     objectList$sample.nobs <- length(residuals(object))
#   }
#   
#   objectList <- sapply(objectList, function(x) eval(x))
#   
#   # multiple objects of class restriktor are allowed
#   isRestr <- unlist(lapply(objectList, function(x) class(x)[1] == "restriktor"))
# 
#   restr_objectList <- list(object = objectList[isRestr])
#   arguments_objectList <- objectList[!isRestr]
#   
#   # put all objects of class restriktor in one list
#   #objectList <- append(list(object = objectList[isRestr]), objectList[!isRestr])
#   
#   res <- do.call(goric.default, c(restr_objectList, arguments_objectList)) 
#   
#   res
# }



# object of class lm ------------------------------------------------------
goric.lm <- function(object, ..., hypotheses = NULL,
                     comparison = NULL,
                     type = "goric",
                     missing = "none", auxiliary = c(), emControl = list(),
                     debug = FALSE) {
  
  if (missing %in% c("em", "EM", "two.stage", "twostage")) {
    missing <- "two.stage" 
  } 
  
  objectList <- list(...)
  
  #mcList <- as.list(match.call())
  #mcList <- mcList[-c(1)]
  
  #mcnames <- names(mcList) == ""
  #lnames <- as.character(mcList[mcnames])
  #names(mcList)[mcnames] <- lnames
  #objectList <- mcList  
  
  # only one object of class lm is allowed
  isLm <- unlist(lapply(objectList, function(x) class(x)[1] %in% c("lm", "glm", "mlm", "rlm")))
  if (sum(isLm) > 1L) {
    stop(paste("Restriktor ERROR: multiple objects of class lm found, only 1 is allowed."), call. = FALSE)
  }
  
  if (missing == "two.stage") {
    if (family(object)$family != "gaussian") {
      stop(paste("Restriktor ERROR: \"two.stage\" is not available in the categorical setting"), call. = FALSE)
    }
    tsm  <- two_stage_matrices(object, auxiliary = auxiliary, emControl = emControl)
    vcov <- two_stage_sandwich(tsm)
    est  <- coef(tsm$fitTarget)
    #N <- tsm$N
    
    if (!type %in% c("gorica", "goricac")) {
      stop(paste("restriktor EROR: missing = \"two.stage\" is only (for now)", 
           "available for type = 'gorica(c)'"), call. = FALSE)
    }
    objectList$object <- est
    objectList$VCOV   <- vcov
  } else {
    objectList$object <- object
    objectList$VCOV   <- NULL
  } 
    
    objectList$hypotheses  <- hypotheses
    objectList$comparison  <- comparison
    objectList$type        <- type
    objectList$debug       <- debug
    objectList$missing     <- NULL
    objectList$auxiliary   <- NULL
    objectList$emControl   <- NULL
    
    if (type == "goricac") {
      objectList$sample_nobs <- length(residuals(object))
    }
    
    # if (missing == "none") {
    #   object_idx <- grepl("object", names(objectList))
    #   objectList <- append(list(object = objectList[object_idx]), objectList[!object_idx])
    # }
    
    if (missing == "two.stage") {
      res <- do.call(goric.numeric, objectList)
    } else {    
      res <- do.call(goric.default, objectList)
    }
  res
}


# object of class numeric -------------------------------------------------
goric.numeric <- function(object, ..., hypotheses = NULL,
                          VCOV = NULL,
                          comparison = NULL,
                          type = "gorica", sample_nobs = NULL,
                          debug = FALSE) {
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop(paste("Restriktor ERROR: object of class numeric is only supported for", 
         "type = 'gorica(c)'."), call. = FALSE)
  }
  
  if (is.null(VCOV)) {
    stop(paste("Restriktor ERROR: the argument VCOV is not found."), call. = FALSE)
  } 
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = object,
    VCOV = VCOV,
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    sample_nobs = sample_nobs,
    debug = debug
  )
  
  # Voeg extra argumenten toe aan de objectList
  extraArgs <- list(...)
  objectList <- c(objectList, extraArgs)
  
  # Roep de goric.default functie aan met de samengestelde lijst
  res <- do.call(goric.default, objectList)
  
  res
}


# object of class lavaan --------------------------------------------------
goric.lavaan <- function(object, ..., hypotheses = NULL,
                         comparison = NULL,
                         type = "gorica",
                         standardized = FALSE,
                         debug = FALSE) {
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop(paste("Restriktor ERROR: object of class lavaan is only supported for", 
         "type = 'gorica(c)'."), call. = FALSE)
  }
  
  est <- con_gorica_est_lav(object, standardized)
  objectList <- list(
    object = est$estimate,
    VCOV = est$VCOV,
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    debug = debug
  )
  
  if (type == "goricac") {
    objectList$sample_nobs <- lavInspect(object, what = "ntotal")
  }
  
  # Voeg extra argumenten toe aan de objectList
  extraArgs <- list(...)
  objectList <- c(objectList, extraArgs)
  
  # Roep de goric.default functie aan met de samengestelde lijst
  res <- do.call(goric.default, objectList)
  
  res
}


# object of class CTmeta --------------------------------------------------
goric.CTmeta <- function(object, ..., hypotheses = NULL,
                         comparison = NULL,
                         type = "gorica", sample_nobs = NULL,
                         debug = FALSE) {
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop(paste("Restriktor ERROR: object of class CTmeta is only supported for",
               "type = 'gorica(c)'."), call. = FALSE)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = coef(object),
    VCOV = vcov(object),
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    debug = debug
  )
  
  # Voeg extra argumenten toe aan de objectList
  extraArgs <- list(...)
  objectList <- c(objectList, extraArgs)
  
  # Roep de goric.default functie aan met de samengestelde lijst
  res <- do.call(goric.default, objectList)
  
  res
}


# object of class rma -----------------------------------------------------
goric.rma <- function(object, ..., hypotheses = NULL,
                      VCOV = NULL,
                      comparison = NULL,
                      type = "gorica", sample_nobs = NULL,
                      debug = FALSE) {
  
  if (!inherits(object, c("rma.uni"))) {
    stop(paste("Restriktor ERROR: the object must be of class 'rma.uni'."), call. = FALSE)
  }
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop(paste("Restriktor ERROR: object of class rma is only supported for type = 'gorica(c)'."), call. = FALSE)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = coef(object),
    VCOV = vcov(object),
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    debug = debug
  )
  
  # Voeg extra argumenten toe aan de objectList
  extraArgs <- list(...)
  objectList <- c(objectList, extraArgs)
  
  # Roep de goric.default functie aan met de samengestelde lijst
  res <- do.call(goric.default, objectList)

  res
}



## lme4
# glmerMod, lmerMod, nlmerMod


# object of class nlmerMod -----------------------------------------------------
goric.nlmerMod <- function(object, ..., hypotheses = NULL,
                           VCOV = NULL, comparison = NULL,
                           type = "gorica", sample_nobs = NULL,
                           debug = FALSE) {
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop(paste("Restriktor ERROR: object of class nlmerMod is only supported for", 
               "type = 'gorica(c)'."), call. = FALSE)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = object@beta,
    VCOV = suppressWarnings(vcov(object)),
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    sample_nobs = sample_nobs,
    debug = debug
  )
  
  names(objectList$object) <- colnames(vcov(object))
  
  # Voeg extra argumenten toe aan de objectList
  extraArgs <- list(...)
  objectList <- c(objectList, extraArgs)
  
  # Roep de goric.default functie aan met de samengestelde lijst
  res <- do.call(goric.default, objectList)
  
  res
}


# object of class glmerMod -----------------------------------------------------
goric.glmerMod <- function(object, ..., hypotheses = NULL,
                           VCOV = NULL, comparison = NULL,
                           type = "gorica", sample_nobs = NULL,
                           debug = FALSE) {
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop(paste("Restriktor ERROR: object of class glmerMod is only supported for", 
               "type = 'gorica(c)'."), call. = FALSE)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = object@beta,
    VCOV = suppressWarnings(vcov(object)),
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    sample_nobs = sample_nobs,
    debug = debug
  )
  
  names(objectList$object) <- colnames(vcov(object))
  
  # Voeg extra argumenten toe aan de objectList
  extraArgs <- list(...)
  objectList <- c(objectList, extraArgs)
  
  # Roep de goric.default functie aan met de samengestelde lijst
  res <- do.call(goric.default, objectList)
  
  res
}


# object of class glmerMod -----------------------------------------------------
goric.lmerMod <- function(object, ..., hypotheses = NULL,
                          VCOV = NULL, comparison = NULL,
                          type = "gorica", sample_nobs = NULL,
                          debug = FALSE) {
  
  if (!c(type %in% c("gorica", "goricac"))) {
    stop(paste("Restriktor ERROR: object of class lmerMod is only supported for", 
               "type = 'gorica(c)'."), call. = FALSE)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = object@beta,
    VCOV = suppressWarnings(vcov(object)),
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    sample_nobs = sample_nobs,
    debug = debug
  )
  
  names(objectList$object) <- colnames(vcov(object))
  
  # Voeg extra argumenten toe aan de objectList
  extraArgs <- list(...)
  objectList <- c(objectList, extraArgs)
  
  # Roep de goric.default functie aan met de samengestelde lijst
  res <- do.call(goric.default, objectList)
  
  res
}


# goric.nlmerMod <- function(object, ..., hypotheses = NULL,
#                       VCOV = NULL,
#                       comparison = NULL,
#                       type = "gorica", sample_nobs = NULL,
#                       debug = FALSE) {
#   
#   if (!c(type %in% c("gorica", "goricac"))) {
#     stop(paste("Restriktor ERROR: object of class nlmerMod is only supported for", 
#          "type = 'gorica(c)'."), call. = FALSE)
#   }
#   
#   objectList <- list(...)
#   objectList$VCOV <- suppressWarnings(vcov(object))
#   cnames <- colnames(objectList$VCOV)
#   objectList$object <- object@beta
#   names(objectList$object) <- cnames
#   objectList$hypotheses <- hypotheses
#   objectList$comparison <- comparison
#   objectList$type <- type
#   objectList$debug <- debug
#   
#   if (type == "goricac") {
#     objectList$sample_nobs <- sample_nobs
#   }
#   
#   object_idx <- grepl("object", names(objectList))
#   objectList <- append(list(object = objectList[object_idx]), objectList[!object_idx])  
#   
#   res <- do.call(goric.default, c(objectList[object_idx], objectList[!object_idx]))
#   
#   res
# }