goric <- function(object, ...) { UseMethod("goric") }


goric.default <- function(object, ..., hypotheses = NULL,
                          comparison = NULL, type = "goric", 
                          VCOV = NULL, sample_nobs = NULL,
                          penalty_factor = 2,
                          Heq = FALSE, control = list(), debug = FALSE) {
  
  # the following classes are allowed (for now)
  obj_class <- class(object)
  classes <- c("aov", "lm", "glm", "mlm", "rlm", "numeric", "lavaan", "CTmeta", 
               "rma.uni", "rma.mv", "lmerMod", "nlmerMod", "glmerMod", "merMod")
  check_class <- obj_class %in% classes
  if (!any(check_class)) {
    stop(paste0("\nrestriktor ERROR: Objects of class", paste(obj_class, collapse = ", "), 
               "are not supported. Supported classes are:", paste(classes, collapse = ", "),"."), call. = FALSE)
  }
  
  if (!is.list(hypotheses) || is.null(hypotheses)) {
    stop(paste("\nrestriktor ERROR: The 'hypotheses' argument is missing or not a list.",
         "Please make sure to provide a valid set of hypotheses, for example, hypotheses =",
         "list(h1 = 'x1 > x2 > x3')."), call. = FALSE)
  } 
  
  if (Heq == TRUE && !is.null(hypotheses$Heq)) { # can happen if it comes from benchmarks function
    hypotheses <- hypotheses[-1]
  }
  num_hypotheses <- length(hypotheses)
  
  # Set default comparison if needed
  if (is.null(comparison)) {
    if (num_hypotheses == 1) {
      comparison <- "complement"
    } else {
      comparison <- "unconstrained"
    }
  }
 
  comparison <- match.arg(comparison, c("unconstrained", "complement", "none"))
  
  # Validate Heq
  if (Heq && num_hypotheses > 1) {
    stop(paste0("\nrestriktor ERROR: Heq = TRUE is only allowed when there is one order-restricted hypothesis. ",
                "Now, there are ", num_hypotheses, " order-restricted hypotheses." ), call. = FALSE)
  }

  # Adjust comparison if necessary based on hypotheses
  if (comparison == "complement" && num_hypotheses > 1) {
    warning("\nrestriktor WARNING: More than one hypothesis provided. Therefore, 'comparison' set to 'unconstrained'.", call. = FALSE)
    comparison <- "unconstrained"
  }
  
  # Ignore Heq for other comparisons
  if (comparison %in% c("unconstrained", "none") && Heq) {
      warning("\nrestriktor WARNING: The 'Heq' argument is ignored. The specified", 
              " hypothesis is only valid when the order-restricted hypothesis is compared",
              " to its complement.", call. = FALSE)
      Heq <- FALSE
  }
  # TO DO kan deze weg of gaat het dan fout in de code?

  
  if (!is.null(VCOV)) {
    # check if scalar
    if (length(VCOV) == 1) {
      VCOV <- as.matrix(VCOV)
    }
    # check if matrix and whether p times p matrix
    if (!is.matrix(VCOV)) {
      VCOV <- as.matrix(VCOV)
    }
    if (dim(VCOV)[1] != dim(VCOV)[2]) {
      stop(paste("\nrestriktor ERROR: The covariance matrix (VCOV) should be a square matrix."), call. = FALSE)
    }
    # check if it is of class matrix
    if (inherits(VCOV, "dpoMatrix")) {
      VCOV <- as.matrix(VCOV)
      # as.matrix(suppressWarnings(vcov(object)))
    }
    if (any(is.na(VCOV))) {
      stop(paste("\nrestriktor ERROR: The covariance matrix (VCOV) contains NA or NaN values.", 
           "Please check your data or model specification."), call. = FALSE)
    }
  }
  
  if (penalty_factor < 0) {
    stop(paste("\nrestriktor ERROR: the penalty factor must be >= 0."), call. = FALSE)
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
    stop("\nrestriktor ERROR: argument ", 
         sQuote(paste(names(ldots)[pm == 0], collapse = ", ")), 
         " unknown.", call. = FALSE)
  }
  
  if (any(c("lower", "upper", "algorithm", "burn.in.samples", "start.values", 
            "thinning") %in% names(ldots))) {
    ldots$mix_weights <- "boot"
  }

  # some checks
  if (!is.null(comparison)) {
    comparison <- tolower(comparison)
  }
  
  #if (length(hypotheses) == 1 & Heq & comparison == "complement") {
  if (length(hypotheses) == 1 & Heq & comparison != "none") {
    Hceq <- gsub("<|>", "=", hypotheses[[1]])
    hypotheses <- append(list(Heq = Hceq), hypotheses)
  } 
    
  constraints <- hypotheses
  # class objects
  object_class <- obj_class

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
  #   #VCOV <- vcov(ans$model.org) 
  #   VCOV <- VCOV.unbiased(ans$model.org, sample_nobs)
  # } else 
  if (any(object_class %in% c("aov","lm","rlm","glm","mlm")) && isConChar) { 
    # TO DO in mlm geeft coef() een matrix, je moet dan as.vector doen (coef.named.vector(...)).
    #       Maakt dat in dit deel nog uit? Ik zie het niet, graag ff checken of dat klopt.
    #       Dan ook in volgende 'if' evt.
    # TO DO voor mlm krijg je  namen die niet werken in goric... Hoe dat aanpassen?
    
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
                                                     goric          = type,
                                                     sample.nobs    = sample_nobs,
                                                     penalty_factor = penalty_factor))
    
    PT_Amat <- lapply(isSummary, function(x) x$PT_Amat)
    PT_meq  <- lapply(isSummary, function(x) x$PT_meq)
    
    for (lnames in names(conList)) {
      conList[[lnames]]$PT_Amat <- PT_Amat[[lnames]]
      conList[[lnames]]$PT_meq <- PT_meq[[lnames]]
    }
    
    ans$hypotheses_usr <- lapply(conList, function(x) x$CON$constraints)
    ans$model.org <- object # add unrestricted object to output
    if (is.null(sample_nobs)) {
      sample_nobs <- nrow(model.frame(object)) # This is sample size N, not group sizes
    } 
    VCOV <- VCOV.unbiased(ans$model.org, sample_nobs) # unrestricted VCOV + check sample_nobs
    idx <- length(conList) 
    objectnames <- vector("character", idx)
  } else if (any(object_class %in% c("aov","lm","rlm","glm","mlm")) && !isConChar) {
    # tolower names Amat and rhs
    for (i in seq_along(constraints)) { 
      names(constraints[[i]]) <- tolower(names(constraints[[i]])) 
      if (any(!names(constraints[[i]]) %in% c("constraints", "rhs", "neq"))) { 
        stop("\nrestriktor ERROR: The list objects must be named 'constraints', 'rhs' and 'neq', e.g.:
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
                                                     goric          = type,
                                                     sample.nobs    = sample_nobs,
                                                     penalty_factor = penalty_factor))
    
    PT_Amat <- lapply(isSummary, function(x) x$PT_Amat)
    PT_meq  <- lapply(isSummary, function(x) x$PT_meq)
    
    for (lnames in names(conList)) {
      conList[[lnames]]$PT_Amat <- PT_Amat[[lnames]]
      conList[[lnames]]$PT_meq <- PT_meq[[lnames]]
    }
    
    ans$model.org <- object # add unrestricted object to output
    if (is.null(sample_nobs)) {
      sample_nobs <- nrow(model.frame(object)) # This is sample size N, not group sizes
    } 
    VCOV <- VCOV.unbiased(ans$model.org, sample_nobs) # unrestricted VCOV + check sample_nobs
    idx <- length(conList) 
    objectnames <- vector("character", idx)
    #CALL$object <- NULL
  } else if ("numeric" %in% object_class && isConChar) {
    
    if (is.null(sample_nobs) && type %in% c("goricac")) {
      stop(paste("\nrestriktor ERROR: the argument sample_nobs is not found,",
                 "which is needed when type = 'goricac'."), call. = FALSE)
    } else if (type %in% c("goricac")) {
      sample_nobs <- check_sample_nobs(sample_nobs) # take sum in case of group sizes
    }
    
    # fit restriktor object for each hypothesis
    conList <- lapply(constraints, function(constraint) {
      CALL.restr <- append(list(object      = object, 
                                constraints = constraint,
                                VCOV        = as.matrix(VCOV)), 
                                ldots)
      do.call("con_gorica_est", CALL.restr)
    })
    
    names(conList) <- names(constraints)
    
    ans$hypotheses_usr <- lapply(conList, function(x) x$CON$constraints)
    
    isSummary <- lapply(conList, function(x) summary(x, 
                                                     type           = type,
                                                     sample.nobs    = sample_nobs,
                                                     penalty_factor = penalty_factor)) 
  } else if ("numeric" %in% object_class && !isConChar) {
    
    if (is.null(sample_nobs) && type %in% c("goricac")) {
      stop(paste("\nrestriktor ERROR: the argument sample_nobs is not found,",
                 "which is needed when type = 'goricac'."), call. = FALSE)
    } else if (type %in% c("goricac")) {
      sample_nobs <- check_sample_nobs(sample_nobs) # take sum in case of group sizes
    }
    
    # tolower names Amat and rhs
    for (i in seq_along(constraints)) { 
      names(constraints[[i]]) <- tolower(names(constraints[[i]]))
      if (any(!names(constraints[[i]]) %in% c("constraints", "rhs", "neq"))) {
        stop("\nrestriktor ERROR: The list objects must be named 'constraints', 'rhs' and 'neq', e.g.:
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
                                                     type           = type,
                                                     sample.nobs    = sample_nobs,
                                                     penalty_factor = penalty_factor)) 
  } else {
    stop("\nrestriktor ERROR: I don't know how to handle an object of class ", paste0(class(object)[1]))
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
    if (Heq) {
      objectnames <- ifelse(names(constraints) == "", paste0("H", seq_along(names(constraints)) - 1), 
                            names(constraints))
    } else {
      objectnames <- ifelse(names(constraints) == "", paste0("H", seq_along(names(constraints))), 
                            names(constraints))
    }
  }

 
  if (comparison == "complement" && length(conList) == 1L && 
      nrow(conList[[1]]$constraints) == conList[[1]]$neq) {
    comparison  <- "unconstrained"
    message("\nrestriktor Message: The complement of a hypothesis with only equality", 
    " constraints is the unconstrained model. Comparison set to 'unconstrained' instead.")
  }
  
  names(conList) <- objectnames
# compute complement ------------------------------------------------------
  df.c <- NULL
   if (comparison == "complement") { 
    Hm <- setdiff(names(conList), "Heq")  
     
    # unrestricted estimates
    #if (inherits(object, "numeric")) {
    if (inherits(conList[[Hm]]$b.unrestr, "numeric")) {
      #b.unrestr <- object
      # TO DO hier gaat het fout als selectie param in hypo
      b.unrestr <- conList[[Hm]]$b.unrestr
      #
      # TO DO Wat als object niet een vector is?
      #       NB Als: b.unrestr <- as.vector(object)
      #          dan heb je geen namen meer (die juist wel in hypo terug moeten komen)
      #          Bij een matrix dan die col of row names uitlezen en gebruiken, of obv vcov!
    } else {
      b.unrestr <- coef.named.vector(ans$model.org)
    }
    # VCOV
    VCOV <- conList[[Hm]]$VCOV
    # restricted estimates
    b.restr <- conList[[Hm]]$b.restr
    # level probabilities
    wt.bar <- conList[[Hm]]$wt.bar
    # constraints matrix
    Amat <- conList[[Hm]]$constraints
    # remove all zero rows
    Amat <- Amat[apply(Amat, 1, function(x) !all(x == 0)), , drop = FALSE]
    # number of equalities
    meq <- conList[[Hm]]$neq
    # rhs
    bvec <- conList[[Hm]]$rhs
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
  # moet dit obv PT_Amat en PT_meq? TO DO
  LL_c <- compute_complement_likelihood(ans$model.org, VCOV,
                                        Amat, Amat.ciq, Amat.ceq, 
                                        bvec, bvec.ciq, bvec.ceq, 
                                        meq, b.unrestr, type, ldots,
                                        debug = debug)
  llc <- LL_c$llc
  betasc <- LL_c$betasc
  
  # compute log-likelihood model
  if (type %in% c("goric", "goricc")) {
    llm <- logLik(conList[[Hm]])
  } else if (type %in% c("gorica", "goricac")) {
    llm <- dmvnorm(c(b.unrestr - b.restr), sigma = VCOV, log = TRUE)
    #llm <- mnormt::dmnorm(c(b.unrestr - b.restr), varcov = VCOV, log = FALSE) 
  }

  # compute complement penalty term value 
  PTc <- penalty_complement_goric(Amat = conList[[Hm]]$PT_Amat, 
                                  meq  = conList[[Hm]]$PT_meq, 
                                  type, wt.bar, 
                                  debug = debug, 
                                  sample.nobs = sample_nobs)   
  } 
   
  
  
  if (comparison == "complement" && Heq) {
    # restricted estimates
    b.restr <- conList$Heq$b.restr
    # compute log-likelihood model
    if (type %in% c("goric", "goricc")) {
      llceq <- logLik(conList$Heq)                                             
    } else if (type %in% c("gorica", "goricac")) {
      llceq <- dmvnorm(c(b.unrestr - b.restr), sigma = VCOV, log = TRUE)
      #llceq <- mnormt::dmnorm(c(b.unrestr - b.restr), varcov = VCOV, log = FALSE) 
    }
    llm <- c(llceq, llm)
    names(llm) <- c("Heq", Hm)
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
                llm <- unlist(lapply(conList, function(x) { 
                  dmvnorm(c(x$b.unrestr - x$b.restr), sigma = VCOV, log = TRUE) 
                  # mnormt::dmnorm(c(x$b.unrestr - x$b.restr), varcov = VCOV, log = FALSE) 
                  } 
                )) 
                # unrestricted
                llu <- dmvnorm(rep(0, ncol(VCOV)), sigma = VCOV, log = TRUE) 
                #llu <- mnormt::dmnorm(rep(0, ncol(VCOV)), varcov = VCOV, log = FALSE) 
            }
            
            if (type %in% c("goric", "gorica")) {
              PTu <- 1 + ncol(VCOV)
            } else if (type %in% c("goricc", "goricac")) {
              
              if (is.null(sample_nobs)) {
                stop(paste("\nrestriktor ERROR: the argument sample_nobs is not found,",
                           "which is needed when type = type = \'goric(a)c\'."), call. = FALSE)
              } else {
                sample_nobs <- check_sample_nobs(sample_nobs) # take sum in case of group sizes
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
              stop("\nrestriktor ERROR: no chi-bar-square weights (a.k.a. level probabilities) are found. Use mix_weights = 'pmvnorm' (default) or 'boot'.", call. = FALSE)
            }
            
            goric.Hm <- -2*llm + penalty_factor*PTm #-2*(llm - PTm)
            goric.Hu <- -2*llu + penalty_factor*PTu #-2*(llu - PTu)
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
              goric.Hm <- -2*llm + penalty_factor*PTm #-2*(llm - PTm)
              df.Hm  <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                                   goric = goric.Hm)
              df.Hm$model <- as.character(df.Hm$model)
              # complement
              goric.Hc <- -2*llc + penalty_factor*PTc #-2*(llc - PTc)
              df.c  <- data.frame(model = "complement", loglik = llc, penalty = PTc, 
                                  goric = goric.Hc)
            } else if (type %in% c("gorica", "goricac")) {
              # model
              gorica.Hm <- -2*llm + penalty_factor*PTm #-2*(llm - PTm)
              df.Hm  <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                                   gorica = gorica.Hm)
              df.Hm$model <- as.character(df.Hm$model)
              # complement
              gorica.Hc <- -2*llc + penalty_factor*PTc #-2*(llc - PTc)
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
            ll <- unlist(lapply(conList, function(x) { 
              dmvnorm(c(x$b.unrestr - x$b.restr), sigma = VCOV, log = TRUE)
              #mnormt::dmnorm(c(x$b.unrestr - x$b.restr), varcov = VCOV, log = FALSE) 
              }
            )) 
            goric.Hm <- -2*ll + penalty_factor*PT #-2*(ll - PT)
            df <- data.frame(model = objectnames, loglik = ll, penalty = PT, 
                             gorica = goric.Hm)
            df$model <- as.character(df$model)
          }
        },
          stop("\nrestriktor ERROR: I don't know how to compute the GORIC values.")
  )

  ans$objectList  <- conList
  ans$objectNames <- objectnames
  
  # calculate LL, PT, goric weights and ratios
  model_comparison_metrics <- calculate_model_comparison_metrics(df)
  
  df$loglik.weights  <- model_comparison_metrics$loglik_weights
  df$penalty.weights <- model_comparison_metrics$penalty_weights
  df$goric.weights   <- model_comparison_metrics$goric_weights
  df$goric.weights_without_unc <- model_comparison_metrics$goric_weights_without_unc
  df$goric.weights_without_heq <- model_comparison_metrics$goric_weights_without_heq
  
  names(df)[4] <- type
  names(df)[7] <- paste0(type, ".weights")
  if (!is.null(df$goric.weights_without_unc)) {
    names(df)[8] <- paste0(type, ".weights_without_unc")
  }
  if (!is.null(df$goric.weights_without_heq)) {
    names(df)[8] <- paste0(type, ".weights_without_heq")
  }
  
  rownames(df) <- NULL

  ans$result <- df
  ans$ratio.gw <- model_comparison_metrics$goric_rw
  ans$ratio.pw <- model_comparison_metrics$penalty_rw
  ans$ratio.lw <- model_comparison_metrics$loglik_rw

  
  # list all object estimates
  coefs <- lapply(conList, FUN = function(x) { coef(x) } )
  # TO DO hier ook aanpassing nodig (nb ws wel namen nodig dan ook)?
  max.length <- max(sapply(coefs, length))
  coefs <- lapply(coefs, function(v) { c(v, rep(NA, max.length - length(v)))})
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
  
  ans$sample_nobs <- sample_nobs
  ans$VCOV <- VCOV 
  ans$b.unrestr <- conList[[1]]$b.unrestr
  ans$ormle$b.restr <- coefs  
  ans$comparison <- comparison
  ans$type <- type
  ans$penalty_factor <- penalty_factor
  ans$Heq <- Heq

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
#     stop("\nrestriktor ERROR: the object must be of class restriktor.")
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
  isLm <- unlist(lapply(objectList, function(x) class(x)[1] %in% c("aov", "lm", "glm", "mlm", "rlm")))
  # TO DO werkte voor rlm niet, zie evt TO DOs in conRLM (wordt via restriktor functie in LL bepaling aangeroepen)
  if (sum(isLm) > 1L) {
    stop(paste("\nrestriktor ERROR: multiple objects of class lm found, only 1 is allowed."), call. = FALSE)
  }
  
  if (missing == "two.stage") {
    if (family(object)$family != "gaussian") {
      stop(paste("\nrestriktor ERROR: \"two.stage\" is not available in the categorical setting"), call. = FALSE)
    }
    #
    # Check on type and possibly change
    if (type == "goric") {
      message("\nrestriktor Message: missing = \"two.stage\" is only (for now)", 
              "available for type = 'gorica(c)'. The GORICA will be used, not the the GORIC.")
      type = "gorica"
    } else if (type == "goricc") {
      message("\nrestriktor Message: missing = \"two.stage\" is only (for now)", 
              "available for type = 'gorica(c)'. The GORICAC will be used, not the the GORICC.")
      type = "goricac"
    } else if (!c(type %in% c("gorica", "goricac"))) {
      message("\nrestriktor Message: missing = \"two.stage\" is only (for now)", 
              "available for type = 'gorica(c)'. The GORICA will be used.")
      type = "gorica"
    }
    #
    tsm  <- two_stage_matrices(object, auxiliary = auxiliary, emControl = emControl)
    vcov <- two_stage_sandwich(tsm)
    est  <- coef.named.vector(tsm$fitTarget, VCOV = vcov)
    #N <- tsm$N
    objectList$sample_nobs <- tsm$N
    objectList$object <- est
    objectList$VCOV   <- vcov
  } else {
    objectList$object <- object
    objectList$VCOV   <- NULL
  } 
    objectList$sample_nobs <- length(residuals(object))
    objectList$hypotheses  <- hypotheses
    objectList$comparison  <- comparison
    objectList$type        <- type
    objectList$debug       <- debug
    objectList$missing     <- NULL
    objectList$auxiliary   <- NULL
    objectList$emControl   <- NULL
    
    if (missing == "two.stage") {
      res <- do.call(goric.numeric, objectList)
    } else {    
      res <- do.call(goric.default, objectList)
    }
  res
}


# object of class numeric -------------------------------------------------
goric.numeric <- function(object, ..., hypotheses = NULL,
                          comparison = NULL,
                          type = "gorica", 
                          VCOV = NULL,
                          sample_nobs = NULL,
                          debug = FALSE) {
  
  # Check on type and possibly change
  check.type <- check.type(type, class = "numeric") 
  
  if (is.null(VCOV)) {
    stop(paste("\nrestriktor ERROR: the argument VCOV is not found."), call. = FALSE)
  } 
  
  #  if (!is.null(dim(object)) && dim(object)[1] == 1) {
  #  object <- c(object)
  #}
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = c(object),
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
  
  # Check on type and possibly change
  check.type <- check.type(type, class = "lavaan") 

  est <- con_gorica_est_lav(object, standardized)
  # Function is defined in 'gorica_est'.
  
  Amat <- con_constraints(model = est$estimate, VCOV = est$VCOV, 
                          constraints = hypotheses)$Amat
  involved_params <- colSums(abs(Amat)) != 0
  # Check whether all parameters are either all defined or all labelled (undefined) ones:
  if (any(est$defined[involved_params] == 1) & any(est$defined[involved_params] == 0)) {
    # Then mix of defined and regular (labelled) parameters
    message(
"restriktor Message: The hypotheses contain labelled and defined parameters. ",
"To account for their covariances, you should also define the labelled ones: ",
"\n", names(est$estimate[involved_params][est$defined[involved_params] == 0]), ".\n",
"For example, if 'c' is the label for the direct effect in a mediaton model, ",
"you could add 'direct := c' to the model syntax.",
# TO DO: "Example code can be found on \url[https://github.com/rebeccakuiper/Tutorials/tree/main/GORICA%20for%20mediation]."
"\n",
"Press [Enter] to continue without covariances ",  
"between the labelled and defined parameters. ",
"Press [Esc] to stop and define the labelled parameters of interest as well.")
    # NOTE: This message needs to be printed
    # so also when other messages are gathered in a list or so
   readline(prompt = "Press [Enter] to continue without covariances. ")
  } 
  
  # message VCOV
  message.VCOV()
  
  objectList <- list(
    object = est$estimate[involved_params],                 
    VCOV = est$VCOV[involved_params, involved_params],
    sample_nobs = lavInspect(object, what = "ntotal"),
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

# object of class CTmeta --------------------------------------------------
goric.CTmeta <- function(object, ..., hypotheses = NULL,
                         comparison = NULL,
                         type = "gorica", 
                         sample_nobs = NULL,
                         debug = FALSE) {
  
  # Check on type and possibly change
  check.type <- check.type(type, class = "CTmeta") 
  
  # message VCOV
  message.VCOVvb()
  
  # sample_nobs (needed for type = "goricac")
  if (type == "goricac") {
    N <- sum(object$N) # or object$NrStudies # TO DO check:
    # Sample sizes for each study - of juist aantal studies!!!!
    sample_nobs <- check_N_with_sample_nobs(N, sample_nobs)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = coef(object), 
    VCOV = vcov(object), 
    sample_nobs = sample_nobs,
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
                      comparison = NULL,
                      type = "gorica", 
                      VCOV = NULL,
                      sample_nobs = NULL,
                      debug = FALSE) {
  
  if (!inherits(object, c("rma.uni", "rma.mv"))) {
    stop(paste("\nrestriktor ERROR: the object must be of class 'rma.uni' or 'rma.mv'."), call. = FALSE)
  }
  
  # Check on type and possibly change
  check.type <- check.type(type, class = "rma (i.e., 'rma.uni' or 'rma.mv')") 
  
  # message VCOV
  message.VCOV()
  
  # sample_nobs (needed for type = "goricac")
  if (type == "goricac") {
    N <- N <- sum(object$ni) # or length(object$ni) # TO DO check:
    # Sample sizes for each study - of juist aantal studies!!!!
    sample_nobs <- check_N_with_sample_nobs(N, sample_nobs)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = coef(object), 
    VCOV = vcov(object), 
    sample_nobs = sample_nobs,
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
# TO DO dit kan ws een functie worden toch? 
#       Maar wellicht op termijn wel verschillen in de 3 (of nu ook al)?


# object of class nlmerMod -----------------------------------------------------
goric.nlmerMod <- function(object, ..., hypotheses = NULL,
                           comparison = NULL,
                           type = "gorica", 
                           VCOV = NULL, sample_nobs = NULL,
                           debug = FALSE) {
  
  # Check on type and possibly change
  check.type <- check.type(type, class = "nlmerMod") 
  
  # message VCOV
  message.VCOV()
  # Note We can get the Variance-Covariance matrix of the random effects (ranef()) 
  # with VarCorr; see https://lmiratrix.github.io/MLM/lmer_extract.html
  
  # sample_nobs (needed for type = "goricac")
  if (type == "goricac") {
    N <- nobs(object)
    sample_nobs <- check_N_with_sample_nobs(N, sample_nobs)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = object@beta,  
    # object@beta equals fixef(object); latter contains labels.
    VCOV = suppressWarnings(vcov(object)), # dpoMatrix, will be made a matrix in goric.default
    sample_nobs = sample_nobs,
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    debug = debug
  )
  
  # Label 'object = object@beta' (not needed when: 'object = fixef(object)')
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
                           comparison = NULL,
                           type = "gorica", 
                           VCOV = NULL, sample_nobs = NULL,
                           debug = FALSE) {
  
  # Check on type and possibly change
  check.type <- check.type(type, class = "glmerMod") 
  
  # message VCOV
  message.VCOV()
  
  # sample_nobs (needed for type = "goricac")
  if (type == "goricac") {
    N <- nobs(object)
    sample_nobs <- check_N_with_sample_nobs(N, sample_nobs)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = object@beta, 
    # object@beta equals fixef(object); latter contains labels.
    VCOV = suppressWarnings(vcov(object)), 
    sample_nobs = sample_nobs,
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    debug = debug
  )
  
  # Label 'object = object@beta' (not needed when: 'object = fixef(object)')
  names(objectList$object) <- colnames(vcov(object))
  
  # Voeg extra argumenten toe aan de objectList
  extraArgs <- list(...)
  objectList <- c(objectList, extraArgs)
  
  # Roep de goric.default functie aan met de samengestelde lijst
  res <- do.call(goric.default, objectList)
  
  res
}


# object of class lmerMod -----------------------------------------------------
goric.lmerMod <- function(object, ..., hypotheses = NULL,
                          comparison = NULL, type = "gorica", 
                          VCOV = NULL, sample_nobs = NULL,
                          debug = FALSE) {
  
  # Check on type and possibly change
  check.type <- check.type(type, class = "lmerMod") 
  
  # message VCOV
  message.VCOV()
  
  # sample_nobs (needed for type = "goricac")
  if (type == "goricac") {
    N <- nobs(object)
    sample_nobs <- check_N_with_sample_nobs(N, sample_nobs)
  }
  
  # Maak de objectList aan en voeg de vereiste elementen toe
  objectList <- list(
    object = object@beta, 
    # object@beta equals fixef(object); latter contains labels.
    VCOV = suppressWarnings(vcov(object)), 
    sample_nobs = sample_nobs,
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    debug = debug
  )
  
  # Label 'object = object@beta' (not needed when: 'object = fixef(object)')
  names(objectList$object) <- colnames(vcov(object))
  
  # Voeg extra argumenten toe aan de objectList
  extraArgs <- list(...)
  objectList <- c(objectList, extraArgs)
  
  # Roep de goric.default functie aan met de samengestelde lijst
  res <- do.call(goric.default, objectList)
  
  res
}
