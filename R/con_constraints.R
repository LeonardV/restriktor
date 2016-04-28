con_constraints_ceq_amat <- function(object, constraints = NULL) {

    # build a bare-bones parameter table for this object
    lavpartable <- con_partable(object, est = TRUE, label = TRUE)

    # parse the constraints
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable = lavpartable)

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
                               partable = lavpartable)

  c(CON$ceq.rhs, CON$cin.rhs)

}
