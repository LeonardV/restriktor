con_constraints <- function(model, VCOV, est, constraints, bvec = NULL, meq = 0L, 
                            debug = FALSE, ...) {
  
  # VCOV en est worden momenteel niet gebruikt
  
  if (!is.character(constraints) && !is.matrix(constraints) && !is.vector(constraints)) {
    stop("restriktor ERROR: 'constraints' must be character or a constraint matrix.", call. = FALSE)
  }
  
  ## build a bare-bones parameter table for this model
  # if model is a numeric vector
  if (is.numeric(model)) {
    parTable <- con_partable_est(model, est = TRUE, label = TRUE)
    parTable_org <- parTable
  } else {
    # if model is a fitted unrestricted object
    parTable <- con_partable(model, est = TRUE, label = TRUE)  
    parTable_org <- parTable
  }
  
  # unlist constraints
  constraints <- unlist(constraints, use.names = FALSE)
  
  if (is.character(constraints)) {
    
    # 'sd' in the hypotheses (e.g., 'g2 - g1 > 0.2*sd'): replaced by the 
    # (residual) standard deviation of the fitted model; see con_sd_substitute()
    constraints <- con_sd_substitute(constraints, model, parTable)
    
    constraints_usr <- constraints
    
    # these operators are typical lavaan operators and are not allowed in the
    # restriktor syntax. If a lavaan model is fitted, the constraints are already
    # in the fitted object of class lavaan. 
    operators <- c("=~", "<~", "~*~", "~~", "~", "\\|", "%")
    
    if (any(grepl(paste(operators, collapse = "|"), constraints)) || any(grepl("[><]{2,}", constraints))) {
        stop(
        paste(
          "restriktor ERROR: error in constraint syntax.",
           "Only the operators \'<, >, ==, =, :=\' are allowed.",
           "See ?restriktor for details or check the website:",
           "https://restriktor.org/tutorial/syntax.html."
        ), call. = FALSE) 
    }
    
    # Initialize output
    OUT <- vector("list", length(constraints))
      
    # Process constraints
    for (i in seq_along(constraints)) {
      # Clean up the constraints
      constraint.syntax <- clean_constraints(constraints[i])
      
      # Process the constraint syntax
      syntax <- process_constraint_syntax(constraint.syntax)
      
      # Expand parentheses and compound constraints
      #expanded_list <- lapply(syntax, function(x) expand_parentheses(x))
      expanded_list <- lapply(syntax, process_abs_and_expand)
      expanded_list <- lapply(expanded_list, expand_compound_constraints)
      
      # Flatten and remove duplicates
      unLIST <- unique(unlist(expanded_list))
      
      # Handle definitions and reformat equality
      def.idx <- grepl(":=", unLIST)
      unLIST[!def.idx] <- gsub("(?<![=:!])=(?!=)", "==",  unLIST[!def.idx], perl = TRUE) #gsub("=", "==", unLIST[!def.idx])
      
      # Store in output
      OUT[[i]] <- paste(unLIST, collapse = "\n")
    }
    
    constraints <- unlist(OUT)
    
    
    # parse the constraints 
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable    = parTable,
                                 debug       = debug,
                                 theta       = parTable$est)
    CON$constraints <- constraints_usr
    FLAT <- lavParseModelString(constraints)
    CON_FLAT <- attr(FLAT, "constraints")
    LIST <- list()
    lhs <- unlist(lapply(CON_FLAT, "[[", "lhs"))
    op  <- unlist(lapply(CON_FLAT, "[[", "op"))
    rhs <- unlist(lapply(CON_FLAT, "[[", "rhs"))
    LIST$lhs <- lhs # names parameters used in hypotheses
    LIST$op  <- op
    LIST$rhs <- rhs #c(LIST$rhs, rhs) 
    
    # The constraints are added to the parameter table, such that (1) the check
    # below on unsupported constraint syntax (e.g., remaining '<', '>', or '='
    # in lhs or rhs) can be performed, and (2) the parameter table that is 
    # stored in the output (e.g., in restriktor objects) also contains them.
    parTable$lhs   <- c(parTable$lhs, LIST$lhs)
    parTable$op    <- c(parTable$op, LIST$op)
    parTable$rhs   <- c(parTable$rhs, LIST$rhs)
    parTable$label <- c(parTable$label, rep("", length(lhs)))

    # nr of equality constraints
    meq  <- nrow(con_constraints_ceq_amat(model, constraints = constraints))
    # right-hand-side
    bvec <- con_constraints_rhs_bvec(model, constraints = constraints)
    # Matrix with all (in)equality constraints,
    # with first meq rows are the equality contraints
    Amat <- con_constraints_con_amat(model, constraints = constraints)
    
    # check for not supported constraint syntax or an empty matrix
    nsc_lhs.idx <- sum(grepl("<|>|=", parTable$lhs))
    nsc_rhs.idx <- sum(grepl("<|>|=", parTable$rhs))
    
    if (all(Amat == 0) || nsc_lhs.idx > 0 || nsc_rhs.idx > 0) {
      stop(paste("restriktor ERROR: I have no idea how to deal with this constraint syntax.",
           "See ?restriktor for details on how to specify the constraint syntax or check the website", 
            "https://restriktor.org/tutorial/syntax.html."), sep = "", call. = FALSE
      )
    }
    
  
    ## In case of abs() the constraints may incorrectly be considered as non-linear. 
    ## Here, we remove the abs() from the constraint function which is redundant 
    ## for determining if the constraints are linear. 
    
    # check if any abs() function exists in string. 
    if (any(grepl("abs\\(.*\\)", c(LIST$lhs, LIST$rhs)))) {
      LIST2 <- LIST
      
      # remove abs( and ) from string
      LIST2$lhs <- gsub("abs\\(|\\)", "", LIST2$lhs)
      LIST2$rhs <- gsub("abs\\(|\\)", "", LIST2$rhs)
      
      parTable_org$free <- seq_len(length(parTable_org$lhs))
      cin.function <- lav_partable_constraints_ciq(partable = parTable_org, con = LIST2)
      ceq.function <- lav_partable_constraints_ceq(partable = parTable_org, con = LIST2)
      
      CON$cin.nonlinear.idx <- lav_constraints_nonlinear_idx(func = cin.function, 
                                                             npar = length(parTable_org$est))
      CON$ceq.nonlinear.idx <- lav_constraints_nonlinear_idx(func = ceq.function, 
                                                             npar = length(parTable_org$est))
    }
  } else if (!is.character(constraints) && !is.null(constraints)) {
    if (is.vector(constraints) ) {
      constraints <- rbind(constraints)
    }
    CON <- NULL
    Amat <- constraints
    bvec <- if (is.null(bvec)) { rep(0L, nrow(Amat)) } else { bvec }
    meq  <- if (is.null(meq)) { 0L } else { meq }
  } else { 
    stop("restriktor ERROR: No restrictions were specified.", call. = FALSE) 
  }
  
  # Check before deleting redundant constraints
  if (meq > nrow(Amat)) { 
    stop(sprintf("restriktor ERROR: The number of equality restrictions (meq = %d) cannot exceed the number of constraints (nrow(Amat) = %d).", meq, nrow(Amat)),
         "\n There might be conflicting and perhaps redundant constraints. \n",
         "When deleting those, goric() should be able to evaluated the adjusted set of hypotheses.", call. = FALSE)
  }
  
  # correct user errors, like x1 < 2 & x1 < 1, x1 < 2 is removed to get a 
  # full row-rank matrix.
  # Note: this only deletes duplicated restrictions; 
  # so, not all redundant restrictions.
  rrc <- remove_redundant_constraints(Amat, bvec, meq)
  Amat <- rrc$constraints
  bvec <- rrc$rhs
  meq <- rrc$meq
  
  if (!(nrow(Amat) == length(bvec))) {
    warning(paste("restriktor WARNING: The number of constraints does not match", 
                  "the \'rhs\' (nrow(Amat) != length(rhs))."))
  }
  if (meq > nrow(Amat)) { 
    stop(sprintf("restriktor ERROR: The number of equality restrictions (meq = %d) cannot exceed the number of constraints (nrow(Amat) = %d).", meq, nrow(Amat)),
         "\n There might be conflicting and perhaps redundant constraints. \n",
         "When deleting those, goric() should be able to evaluated the adjusted set of hypotheses.", call. = FALSE)
  }
  
  if (length(CON$ceq.nonlinear.idx) > 0L || length(CON$cin.nonlinear.idx) > 0L) {
    stop("restriktor ERROR: nonlinear (in)equality restrictions are not supported.", call. = FALSE)
  }
  
  if (debug && is.character(constraints)) {
    print(as.data.frame(parTable, stringsAsFactors = FALSE))
    print(CON)
  }
  
  OUT <- list(CON      = CON, 
              parTable = parTable,
              Amat     = Amat,
              bvec     = bvec, 
              meq      = meq)
  
  OUT
}


con_constraints_ceq_amat <- function(object, constraints = NULL) {

  # build a bare-bones parameter table for this object
  if ("numeric" %in% class(object)) {
    lavpartable <- con_partable_est(object, est = TRUE, label = TRUE)
  } else {
  # if object is a fitted unrestricted object
    lavpartable <- con_partable(object, est = TRUE, label = TRUE)  
  }

  #lavpartable <- con_partable(object, est = TRUE, label = TRUE)
  

  # parse the constraints
  CON <- lav_constraints_parse(constraints = constraints,
                               partable    = lavpartable,
                               theta       = lavpartable$est)

  CON$ceq.JAC
}


con_constraints_con_amat <- function(object, constraints = NULL) {
  
  # build a bare-bones parameter table for this object
  if ("numeric" %in% class(object)) {
    lavpartable <- con_partable_est(object, est = TRUE, label = TRUE)
  } else {
    # if object is a fitted unrestricted object
    lavpartable <- con_partable(object, est = TRUE, label = TRUE)  
  }
  
  #lavpartable <- con_partable(object, est = TRUE, label = TRUE)

  # parse the constraints
  CON <- lav_constraints_parse(constraints = constraints,
                               partable    = lavpartable, 
                               theta       = lavpartable$est)

  rbind(CON$ceq.JAC, CON$cin.JAC)
}



con_constraints_rhs_bvec <- function(object, constraints = NULL) {

  # build a bare-bones parameter table for this object
  #lavpartable <- con_partable(object, est = TRUE, label = TRUE)

  # build a bare-bones parameter table for this object
  if ("numeric" %in% class(object)) {
    lavpartable <- con_partable_est(object, est = TRUE, label = TRUE)
  } else {
    # if object is a fitted unrestricted object
    lavpartable <- con_partable(object, est = TRUE, label = TRUE)  
  }
  
  
  # parse the constraints
  CON <- lav_constraints_parse(constraints = constraints,
                               partable    = lavpartable,
                               theta       = lavpartable$est)

  c(CON$ceq.rhs, CON$cin.rhs)
}


# 'sd' in hypotheses ---------------------------------------------------------
# Hypotheses can be specified in terms of an effect size, like Cohen's d, 
# using 'sd': for example, 'g2 - g1 > 0.2*sd' denotes that the difference 
# between the means of groups 2 and 1 is at least small (d = (mu2 - mu1)/sd
# of at least 0.2). Here, 'sd' is replaced by the residual standard deviation
# of the fitted model, which (in an ANOVA model) is the pooled within-group 
# standard deviation used in Cohen's d. With covariates, it is the residual
# standard deviation given the covariates. For a multivariate model (mlm), 
# there is one per outcome, denoted by 'sd.<outcome>' (e.g., 'sd.y1').
# This is only available for fitted model objects (lm, aov, mlm, rlm, and 
# glm with family gaussian); with estimates as input, one can plug in the 
# value oneself (e.g., 'g2 - g1 > 0.2*1.5').
con_sd_substitute <- function(constraints, model, parTable = NULL) {
  pattern <- "(?<![A-Za-z0-9._])sd(\\.[A-Za-z0-9_]+)?(?![A-Za-z0-9._(])"
  if (!any(grepl(pattern, constraints, perl = TRUE))) return(constraints)
  tokens <- unique(unlist(regmatches(constraints, gregexpr(pattern, constraints, perl = TRUE))))
  # a parameter (label) called 'sd' (or 'sd.<...>') is not replaced
  labels <- if (!is.null(parTable)) parTable$label else character(0)
  tokens <- setdiff(tokens, labels)
  if (length(tokens) == 0L) return(constraints)
  
  # Check that 'sd' is used as an effect-size cut-off (as in Cohen's d): each 
  # (part of a) constraint with 'sd' should contain at least one parameter 
  # (so, not, e.g., 'sd > 2'), and 'sd' should not multiply a parameter (so,
  # not, e.g., 'g2 - g1 > sd*g3'); allowed are, e.g., 'g2 - g1 > 0.2*sd' and 
  # '(g2 - g1)/sd > 0.2'.
  params <- labels[!is.na(labels) & labels != ""]
  esc <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  par_core <- if (length(params) > 0) {
    paste0("(", paste(esc(params[order(-nchar(params))]), collapse = "|"), ")")
  } else NULL
  par_re <- if (!is.null(par_core)) paste0("(?<![A-Za-z0-9._])", par_core, "(?![A-Za-z0-9._])") else NULL
  tok_re <- paste0("(", paste(gsub(".", "\\\\.", tokens, fixed = TRUE), collapse = "|"), ")")
  parts <- trimws(unlist(strsplit(constraints, "[;\n]")))
  parts <- parts[grepl(paste0("(?<![A-Za-z0-9._])", tok_re, "(?![A-Za-z0-9._(])"), parts, perl = TRUE)]
  sd_use <- paste0("'sd' should be used as an effect-size cut-off (as in Cohen's d), that is, a ",
                   "number times sd (e.g., 'g2 - g1 > 0.2*sd') or a difference divided by sd ",
                   "(e.g., '(g2 - g1)/sd > 0.2').")
  for (p in parts) {
    if (!is.null(par_re) && !grepl(par_re, p, perl = TRUE)) {
      stop(paste0("restriktor ERROR: The constraint '", p, "' does not contain any parameter. ", sd_use),
           call. = FALSE)
    }
    if (!is.null(par_re)) {
      # sd * parameter (or '('), parameter (or ')') * sd, and sd / parameter (or '(');
      # a number times or divided by sd (e.g., 0.2*sd or sd/2) is fine
      mult <- c(paste0(tok_re, "\\s*\\*\\s*(\\(|", par_core, "(?![A-Za-z0-9._]))"),
                paste0("(\\)|(?<![A-Za-z0-9._])", par_core, ")\\s*\\*\\s*", tok_re),
                paste0(tok_re, "\\s*/\\s*(\\(|", par_core, "(?![A-Za-z0-9._]))"))
      if (any(vapply(mult, function(re) grepl(re, p, perl = TRUE), logical(1)))) {
        stop(paste0("restriktor ERROR: In the constraint '", p, "', 'sd' multiplies (or is divided by) ",
                    "a parameter. ", sd_use), call. = FALSE)
      }
    }
  }
  
  if (is.numeric(model)) {
    stop(paste0("restriktor ERROR: 'sd' in the hypotheses (", paste(tokens, collapse = ", "), 
                ") can only be used with a fitted model object (e.g., lm), from which the ",
                "(residual) standard deviation is taken. With estimates as input, please ",
                "plug in the value of the standard deviation yourself (e.g., 'g2 - g1 > 0.2*1.5')."),
         call. = FALSE)
  }
  # (residual) standard deviation(s) of the fitted model
  if (inherits(model, "rlm")) {
    sds <- model$s
  } else if (inherits(model, "glm") && !identical(family(model)$family, "gaussian")) {
    stop(paste0("restriktor ERROR: 'sd' in the hypotheses is only available for models with a ",
                "residual standard deviation, which is not the case for a glm with family '",
                family(model)$family, "'."), call. = FALSE)
  } else {
    sds <- tryCatch(stats::sigma(model), error = function(e) NULL)
  }
  if (is.null(sds) || !is.numeric(sds) || anyNA(sds)) {
    stop("restriktor ERROR: 'sd' in the hypotheses cannot be used for this model, since its ",
         "residual standard deviation cannot be obtained.", call. = FALSE)
  }
  multi <- length(sds) > 1L
  if (multi && is.null(names(sds))) names(sds) <- paste0("y", seq_along(sds))
  
  values <- vapply(tokens, function(tok) {
    if (tok == "sd") {
      if (multi) {
        stop(paste0("restriktor ERROR: The model has more than one outcome; please use 'sd.<outcome>' ",
                    "in the hypotheses, that is, one of: ", 
                    paste0("sd.", names(sds), collapse = ", "), "."), call. = FALSE)
      }
      return(unname(sds[1]))
    }
    outc <- sub("^sd\\.", "", tok)
    if (!multi || !(outc %in% names(sds))) {
      stop(paste0("restriktor ERROR: '", tok, "' in the hypotheses is unknown. ",
                  if (multi) paste0("Please use one of: ", paste0("sd.", names(sds), collapse = ", "), ".")
                  else "Please use 'sd'."), call. = FALSE)
    }
    unname(sds[outc])
  }, numeric(1))
  
  for (i in seq_along(tokens)) {
    pat_i <- paste0("(?<![A-Za-z0-9._])", gsub(".", "\\.", tokens[i], fixed = TRUE), 
                    "(?![A-Za-z0-9._(])")
    constraints <- gsub(pat_i, format(values[i], digits = 15), constraints, perl = TRUE)
  }
  message(paste0("\nrestriktor Message: In the hypotheses, ",
                 paste0("'", tokens, "' is replaced by ", format(values, digits = 4), collapse = "; "),
                 ", that is, the residual standard deviation of the model (as in Cohen's d)."))
  constraints
}
