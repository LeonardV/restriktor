restriktor <- function(model, constraints, se = "default",
                       rhs = NULL, neq = NULL, control = NULL,
                       debug = FALSE, ...) {

  if (!(is.null(weights(model)))) {
    stop("Restriktor ERROR: weights are not possible (yet).")
  }
  bvec <- rhs 
  meq <- neq
  
  if (is.character(constraints) && class(model) %in% c("lm", "rlm", "glm", "mlm")) {
    # build a bare-bones parameter table for this model
    partable <- lav_partable(model, est = TRUE, label = TRUE)
    # parse the constraints
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable = partable,
                                 debug = debug)
    # equality constraints
    meqw  <- nrow(con_constraints_ceq_amat(model, constraints = constraints))
    # right-hand-side
    bvecw <- con_constraints_rhs_bvec(model, constraints=constraints)
    # inequality constraints
    Amatw <- con_constraints_con_amat(model, constraints = constraints)
  }  else if (is.vector(constraints) | is.matrix(constraints)) {
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
