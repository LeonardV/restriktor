con_constraints <- function(model, VCOV, est, constraints, bvec = NULL, meq = 0L, 
                            debug = FALSE, ...) {
  
  ## build a bare-bones parameter table for this model
  # if model is a numeric vector
  if ("numeric" %in% class(model)) {
    parTable <- con_partable_est(model, est = TRUE, label = TRUE)
    parTable_org <- parTable
  } else {
    # if model is a fitted unrestricted object
    parTable <- con_partable(model, est = TRUE, label = TRUE)  
    parTable_org <- parTable
  }
  
  # unlist constraints
  constraints <- unlist(constraints)
  # for summary function
  constraints_usr <- constraints
  
  if (is.character(constraints)) {
    # these operators are typical lavaan operators and are not allowed in the
    # restriktor syntax. If a lavaan model is fitted, the contraints are already
    # in the fitted object of class lavaan. 
    operators <- c("=~", "<~", "~*~", "~~", "~", "\\|", "%")
    
    # check for user input error
    if (grepl(paste(operators, collapse = "|"), constraints) || all(grepl("[><]{2,}", constraints))) {
      stop(paste("Restriktor ERROR: error in constraint syntax. Only the operators \'<, >, ==, =, :=\' are allowed.",
           "See ?restriktor for details on how to specify the constraint syntax or check the website:",
           "https://restriktor.org/tutorial/syntax.html."), call. = FALSE) 
    }
    
    # deal with constraints of format x1 < x2 < x3
    OUT <- list()
    for (i in seq_len(length(constraints))) {
      # some constraint cleanup
      constraint.syntax <- gsub("[#!].*(?=\n)", "", constraints  , perl = TRUE)
      constraint.syntax <- gsub(";", "\n", constraint.syntax     , perl = TRUE)
      #constraint.syntax <- gsub(",", "\n", constraint.syntax     , perl = TRUE)
      constraint.syntax <- gsub("&", "\n", constraint.syntax     , perl = TRUE)
      constraint.syntax <- gsub("[ \t]+", "", constraint.syntax  , perl = TRUE)
      constraint.syntax <- gsub("\n{2,}", "\n", constraint.syntax, perl = TRUE)
      # the regular expression "[\\(\\)]" is a pattern that matches either an opening 
      # parenthesis ( or a closing parenthesis ).
      if (length(gregexpr("[\\(\\)]", constraints)[[1]]) %% 2 == 0) {
        constraint.syntax <- gsub("\\),", "\\)\n", constraint.syntax, perl = TRUE)
      } else {
        constraint.syntax <- gsub(",", "\n", constraint.syntax, perl = TRUE)
      }
      
      constraint.syntax[[i]] <- strsplit(constraint.syntax[[i]], split = "\n", perl = TRUE)
      constraint.syntax[[i]] <- unlist(constraint.syntax[[i]])
      constraint.syntax[[i]] <- gsub("==","=", constraint.syntax[[i]], perl = TRUE)
      constraint.syntax[[i]] <- gsub("<=","<", constraint.syntax[[i]], perl = TRUE)
      constraint.syntax[[i]] <- gsub(">=",">", constraint.syntax[[i]], perl = TRUE)
      constraint.syntax <- unlist(constraint.syntax[i])
      
      LIST <- lapply(constraint.syntax, function(x) { sapply(x, expand_parentheses) })
      LIST <- lapply(LIST, function(x) { sapply(x, expand_compound_constraints) })
      
      unLIST <- unique(unlist(LIST))
      def.idx  <- grepl(":=", unLIST)
      unLIST[!def.idx] <- gsub("=", "==", unLIST[!def.idx])
      OUT[[i]] <- paste(unLIST, collapse = '\n')
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
    LIST$lhs <- lhs
    LIST$op  <- op
    LIST$rhs <- c(LIST$rhs, rhs) 
    
    parTable$lhs   <- c(parTable$lhs, LIST$lhs)
    parTable$op    <- c(parTable$op, LIST$op)
    parTable$rhs   <- c(parTable$rhs, LIST$rhs)
    parTable$label <- c(parTable$label, rep("", length(lhs)))
    
    # equality constraints
    meq  <- nrow(con_constraints_ceq_amat(model, constraints = constraints))
    # right-hand-side
    bvec <- con_constraints_rhs_bvec(model, constraints = constraints)
    # inequality constraints
    Amat <- con_constraints_con_amat(model, constraints = constraints)
    
    # check for not supported constraint syntax or an empty matrix
    nsc_lhs.idx <- sum(grepl("<|>|=", parTable$lhs))
    nsc_rhs.idx <- sum(grepl("<|>|=", parTable$rhs))
    
    if (all(Amat == 0) || nsc_lhs.idx > 0 || nsc_rhs.idx > 0) {
      stop(paste("Restriktor ERROR: I have no idea how to deal with this constraint syntax.",
           "See ?restriktor for details on how to specify the constraint syntax or check the website", 
            "https://restriktor.org/tutorial/syntax.html."), sep = "", call. = FALSE
      )
    }
    
    ## In case of abs() the constraints may incorrectly be considered as non-linear. 
    ## Here, we remove the abs() from the constraint function which is redundant 
    ## for determining if the constraints are linear. 
    
    # check if any abs() functie exists in string. 
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
    stop("no restrictions were specified.") 
  }
  
  # correct user errors, like x1 < 2 & x1 < 1, x1 < 2 is removed to get a 
  # full row-rank matrix.
  rrc <- remove_redundant_constraints(Amat, bvec)
  Amat <- rrc$constraints
  bvec <- rrc$rhs
  
  if (!(nrow(Amat) == length(bvec))) {
    warning(paste("restriktor WARNING: The number of constraints does not match", 
                  "the \'rhs\' (nrow(Amat) != length(rhs))."))
  }
  
  if (meq > nrow(Amat)) { 
    stop("restriktor ERROR: The maximum number of equality constraints = ", nrow(Amat), "\n")
  }
  
  if (length(CON$ceq.nonlinear.idx) > 0L || length(CON$cin.nonlinear.idx) > 0L) {
    stop(paste("restriktor ERROR: can not handle (yet) nonlinear (in)equality restriktions"))
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
