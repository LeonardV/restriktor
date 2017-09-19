con_constraints <- function(model, constraints, bvec = NULL, meq = 0L, 
                            debug = FALSE, ...) {
  
  # build a bare-bones parameter table for this model
  parTable <- con_partable(model, est = FALSE, label = TRUE)
  
  if (is.character(constraints)) {
    # parse the constraints
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable = parTable,
                                 debug = debug)
    
    FLAT <- lavParseModelString(constraints)
    CON_FLAT <- attr(FLAT, "constraints")
    LIST <- list()
    lhs <- unlist(lapply(CON_FLAT, "[[", "lhs"))
    op  <- unlist(lapply(CON_FLAT, "[[", "op"))
    rhs <- unlist(lapply(CON_FLAT, "[[", "rhs"))
    LIST$lhs <- lhs
    LIST$op  <- op
    LIST$rhs <- c(LIST$rhs, rhs)
    
    parTable$lhs <- c(parTable$lhs, LIST$lhs)
    parTable$op <- c(parTable$op, LIST$op)
    parTable$rhs <- c(parTable$rhs, LIST$rhs)
    parTable$label <- c(parTable$label, rep("", length(lhs)))
    
    # equality constraints
    meq  <- nrow(con_constraints_ceq_amat(model, constraints = constraints))
    # right-hand-side
    bvec <- con_constraints_rhs_bvec(model, constraints = constraints)
    # inequality constraints
    Amat <- con_constraints_con_amat(model, constraints = constraints)
    
    CON$constraints <- constraints
  } else if (!is.character(constraints) && !is.null(constraints)) {
    if (is.vector(constraints)) {
      constraints <- rbind(constraints)
    }
    CON <- NULL
    Amat <- constraints
    bvec <- if (is.null(bvec)) { rep(0L, nrow(Amat)) } else { bvec }
    meq  <- if (is.null(meq)) { 0L } else { meq }
  } else { 
    stop("no restriktions were specified.") 
  }
  
  if (length(CON$ceq.nonlinear.idx) > 0L || length(CON$cin.nonlinear.idx) > 0L) {
    stop("restriktor ERROR: can not handle (yet) nonlinear (in)equality restriktions")
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
    lavpartable <- con_partable(object, est = TRUE, label = TRUE)

    # parse the constraints
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable    = lavpartable)

    CON$ceq.JAC
}


con_constraints_con_amat <- function(object, constraints = NULL) {

    # build a bare-bones parameter table for this object
    lavpartable <- con_partable(object, est = TRUE, label = TRUE)

    # parse the constraints
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable = lavpartable)

    rbind(CON$ceq.JAC, CON$cin.JAC)
}



con_constraints_rhs_bvec <- function(object, constraints = NULL) {

  # build a bare-bones parameter table for this object
  lavpartable <- con_partable(object, est = TRUE, label = TRUE)

  # parse the constraints
  CON <- lav_constraints_parse(constraints = constraints,
                               partable    = lavpartable)

  c(CON$ceq.rhs, CON$cin.rhs)
}
