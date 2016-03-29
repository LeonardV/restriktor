restriktor <- function(model, constraints, se = "default",
                       rhs = NULL, neq = NULL, control = NULL,
                       debug = FALSE, ...) {

  if (any(c("glm", "mlm") %in% class(model))) 
    stop("Restriktor does not work (yet) on classes glm or mlm.")
  
  # rename for internal use
  bvec <- rhs 
  meq <- neq
  
  if (is.character(constraints) && class(model) %in% c("lm", "rlm")) {
    # build a bare-bones parameter table for this model
    partable <- lav_partable(model, est = FALSE, label = TRUE)
    
    # parse the constraints
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable = partable,
                                 debug = debug)
    
    FLAT <- lavaan:::lavParseModelString(constraints)
    CON_FLAT <- attr(FLAT, "constraints")
    LIST <- list()
    lhs = unlist(lapply(CON_FLAT, "[[", "lhs"))
    op = unlist(lapply(CON_FLAT, "[[", "op"))
    rhs = unlist(lapply(CON_FLAT, "[[", "rhs"))
    LIST$lhs <- lhs
    LIST$op <- op
    LIST$rhs <- c(LIST$rhs, rhs)
    
    partable$lhs <- c(partable$lhs, LIST$lhs)
    partable$op <- c(partable$op, LIST$op)
    partable$rhs <- c(partable$rhs, LIST$rhs)
    #def.idx <- which(LIST$op == ":=")
    partable$label <- c(partable$label, rep("", length(lhs)))
    #partable$label[def.idx] <- LIST$lhs[def.idx]
    
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
    stop("no constraints specified.") 
  }

  if (debug && is.character(constraints)) {
    print(as.data.frame(lavpartable, stringsAsFactors = FALSE))
    print(CON)
  }

  if ("lm" %in% class(model)[1]) {
    UseMethod("conLM")
  } else if ("rlm" %in% class(model)[1]) {
      UseMethod("conRLM")
  }
  
}
