################################### TO DO ######################################
# - add weights  + down weights to RLM.
# - implement E-bar test statistic for lm.
# - add options for no intercept models.
# - improve code description.

###############################################################################
restriktor <- function(model, constraints, se = "default",
                       rhs = NULL, neq = NULL, control = NULL,
                       debug = FALSE, ...) {
  #if (any(c("glm", "mlm", "gls") %in% class(model))) {
  if (any(c("glm", "mlm", "gls") %in% class(model))) {
    stop("Restriktor does not work on classes glm, mlm and gls (yet).")
  }
  
  # rename for internal use
  bvec <- rhs 
  meq <- neq
  
  # build a bare-bones parameter table for this model
  parTable <- lav_partable(model, est = FALSE, label = TRUE)
  
  if (is.character(constraints)) {#} && class(model) %in% c("lm", "rlm")) {
    # parse the constraints
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable = parTable,
                                 debug = debug)
    
    FLAT <- lavaan:::lavParseModelString(constraints)
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
    meqw  <- nrow(con_constraints_ceq_amat(model, constraints = constraints))
    # right-hand-side
    bvecw <- con_constraints_rhs_bvec(model, constraints=constraints)
    # inequality constraints
    Amatw <- con_constraints_con_amat(model, constraints = constraints)
  }  else if (!is.character(constraints)) {
      if (is.vector(constraints)) {
        constraints <- rbind(constraints)
      }
      Amatw <- constraints
      bvecw <- if (is.null(bvec)) { rep(0L, nrow(Amatw)) } else { bvec }
      meqw  <- if (is.null(meq)) { 0L } else { meq }
  } else { 
    stop("no constraints were specified.") 
  }

  if (debug && is.character(constraints)) {
    print(as.data.frame(lavparTable, stringsAsFactors = FALSE))
    print(CON)
  }

  if (class(model)[1] %in% c("lm","mlm")) {
    UseMethod("conLM")
  } else if (class(model)[1] %in% "rlm") {
    UseMethod("conRLM")
  }
  
}
