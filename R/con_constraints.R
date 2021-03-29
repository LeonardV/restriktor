con_constraints <- function(model, VCOV, est, constraints, bvec = NULL, meq = 0L, 
                            debug = FALSE, ...) {
  
  ## build a bare-bones parameter table for this model
  # if model is a numeric vecter
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
  
  if (is.character(constraints)) {
    # parse the constraints 
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable    = parTable,
                                 debug       = debug,
                                 theta       = parTable$est)
     
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
    
    if (all(Amat == 0)) {
      stop("Restriktor ERROR: constraints are not correctly specified. 
                    See ?restriktor for details.")
    }
    
    # In case of abs() the contraints may incorretly be considered as non-linear. 
    # Here, we remove the abs() from the constraint function which is redundant 
    # for determining if the constraints are linear. 
    
    # check if any abs() functie exists in string. 
    if (any(grepl("abs\\(.*\\)", c(LIST$lhs, LIST$rhs)))) {
      LIST2 <- LIST
      
      # reomve abs( and ) from string
      LIST2$lhs <- gsub("abs\\(|\\)", "", LIST2$lhs)
      LIST2$rhs <- gsub("abs\\(|\\)", "", LIST2$rhs)
      
      parTable_org$free <- seq_len(length(parTable_org$lhs))
      cin.function <- lav_partable_constraints_ciq(partable = parTable_org, con = LIST2)
      ceq.function <- lav_partable_constraints_ceq(partable = parTable_org, con = LIST2)
      
      CON$cin.nonlinear.idx <- con_constraints_nonlinear_idx(func = cin.function, 
                                                             npar = length(parTable_org$est))
      CON$ceq.nonlinear.idx <- con_constraints_nonlinear_idx(func = ceq.function, 
                                                             npar = length(parTable_org$est))
    }
    
    
    CON$constraints <- constraints
  } else if (!is.character(constraints) && !is.null(constraints)) {
    if (is.vector(constraints) ) {
      constraints <- rbind(constraints)
    }
    CON <- NULL
    Amat <- constraints
    bvec <- if (is.null(bvec)) { rep(0L, nrow(Amat)) } else { bvec }
    meq  <- if (is.null(meq)) { 0L } else { meq }
  } else { 
    stop("no restriktions were specified.") 
  }
  
  if (!(nrow(Amat) == length(bvec))) {
    warning("restriktor WARNING: The number of constraints does not match 
                    the \'rhs\' (nrow(Amat) != length(rhs)).")
  }
  
  if (meq > nrow(Amat)) { 
    stop("restriktor ERROR: The maximum number of equality constraints = ", nrow(Amat), "\n")
  }
  
  if (length(CON$ceq.nonlinear.idx) > 0L || length(CON$cin.nonlinear.idx) > 0L) {
    stop("restriktor ERROR: can not handle (yet) nonlinear (in)equality restriktions")
  }
  
  if (debug && is.character(constraints)) {
    print(as.data.frame(parTable, stringsAsFactors = FALSE))
    print(CON)
  }
  
  
  # rAmat <- GaussianElimination(t(Amat))

  ## still to catch 
  #H1 <- 'x1 < 4; x1 > 1' # range restrictie
  #H1 <- 'x1 < 1; x1 > 1' # equality
  #H1 <- 'x1 > 3; x1 > 4' # 
  #H1 <- 'x1 > -1; x1 > 4'#
  
  # if (mix.weights == "pmvnorm") {
  #   if (rAmat$rank < nrow(Amat) && rAmat$rank != 0L) {
  #     ## check for inconsistent constraints: quadprog gives an error if constraints
  #     ## are inconsistent
  #     # consistent.check <- con_solver_gorica(est  = est, 
  #     #                                       VCOV = VCOV, 
  #     #                                       Amat = Amat, 
  #     #                                       bvec = bvec, 
  #     #                                       meq  = meq)
  #     
  #     ## remove any linear dependent rows from the constraint matrix. Amat
  #     ## must be of full row rank.
  #     # remove any zero vectors
  #     allZero.idx <- rowSums(abs(Amat)) == 0
  #     Amat <- Amat[!allZero.idx, , drop = FALSE]
  #     bvec <- bvec[!allZero.idx]
  #     # rank Amat
  #     rank <- qr(Amat)$rank 
  #     # singular value decomposition
  #     s <- svd(Amat)
  #     # continue untill Amat is of full-row rank
  #     while (rank != length(s$d)) {
  #       # check which singular values are zero
  #       zero.idx <- which(zapsmall(s$d) <= 1e-16)
  #       # remove linear dependent rows and reconstruct the constraint matrix
  #       Amat <- s$u[-zero.idx, ] %*% diag(s$d) %*% t(s$v)
  #       # zapping small ones to zero
  #       Amat <- zapsmall(Amat)
  #       bvec <- bvec[-zero.idx]
  #       s <- svd(Amat)
  #       if (debug) {
  #         cat("rank = ", rank, " ... non-zero det. = ", length(s$d), "\n")
  #       }
  #     }
  #   }
  # } else if (rAmat$rank < nrow(Amat) &&
  #            !(se %in% c("none", "boot.model.based", "boot.standard")) &&
  #            rAmat$rank != 0L) {
  #   warning(paste("Restriktor Warning: No standard errors could be computed.
  #                     The constraint matrix must be full row-rank.
  #                     Try to set se = \"none\", \"boot.model.based\" or \"boot.standard\".")) 
  # }
  
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
