# this is originally a function from the lavaan package.
# slightly adapted by LV (17 may 2024)
lav_constraints_parse <- function(partable = NULL, constraints = NULL, theta = NULL, 
                                  debug = FALSE) {
  if (is.null(partable$free)) {
    partable$free <- seq_len(length(partable$lhs))
  }
  if (!is.null(theta)) {
  }
  else if (!is.null(partable$est)) {
    theta <- partable$est[partable$free > 0L]
  }
  else if (!is.null(partable$start)) {
    theta <- partable$start[partable$free > 0L]
  }
  else {
    theta <- rep(0, length(partable$lhs))
  }
  npar <- length(theta)
  if (is.null(constraints)) {
    LIST <- NULL
  }
  else if (!is.character(constraints)) {
    stop("lavaan ERROR: constraints should be a string")
  }
  else {
    FLAT <- lavParseModelString(constraints)
    CON <- attr(FLAT, "constraints")
    LIST <- list()
    if (length(CON) > 0L) {
      lhs = unlist(lapply(CON, "[[", "lhs"))
      op = unlist(lapply(CON, "[[", "op"))
      rhs = unlist(lapply(CON, "[[", "rhs"))
      LIST$lhs <- c(LIST$lhs, lhs)
      LIST$op <- c(LIST$op, op)
      LIST$rhs <- c(LIST$rhs, rhs)
    }
    else {
      stop("lavaan ERROR: no constraints found in constraints argument")
    }
  }
  ceq.simple <- FALSE
  if (!is.null(partable$unco)) {
    ceq.simple <- TRUE
  }
  def.function <- lav_partable_constraints_def(partable, con = LIST, 
                                               debug = debug)
  ceq.function <- lav_partable_constraints_ceq(partable, con = LIST, 
                                               debug = debug)
  ceq.linear.idx <- lav_constraints_linear_idx(func = ceq.function, 
                                               npar = npar)
  ceq.nonlinear.idx <- lav_constraints_nonlinear_idx(func = ceq.function, 
                                                     npar = npar)
  cin.function <- lav_partable_constraints_ciq(partable, con = LIST, 
                                               debug = debug)
  cin.linear.idx <- lav_constraints_linear_idx(func = cin.function, 
                                               npar = npar)
  cin.nonlinear.idx <- lav_constraints_nonlinear_idx(func = cin.function, 
                                                     npar = npar)
  if (!is.null(body(ceq.function))) {
    # ceq.JAC <- try(lav_func_jacobian_complex(func = ceq.function, 
    #                                          x = theta), silent = TRUE)
#    if (inherits(ceq.JAC, "try-error")) {
      ceq.JAC <- lav_func_jacobian_simple(func = ceq.function, 
                                          x = theta)
#    }
    ceq.rhs <- -1 * ceq.function(numeric(npar))
    ceq.theta <- ceq.function(theta)
  }
  else {
    ceq.JAC <- matrix(0, nrow = 0L, ncol = npar)
    ceq.rhs <- numeric(0L)
    ceq.theta <- numeric(0L)
  }
  if (!is.null(body(cin.function))) {
    # cin.JAC <- try(lav_func_jacobian_complex(func = cin.function, 
    #                                          x = theta), silent = TRUE)
#    if (inherits(cin.JAC, "try-error")) {
      cin.JAC <- lav_func_jacobian_simple(func = cin.function, 
                                          x = theta)
#    }
    cin.rhs <- -1 * cin.function(numeric(npar))
    cin.theta <- cin.function(theta)
  }
  else {
    cin.JAC <- matrix(0, nrow = 0L, ncol = npar)
    cin.rhs <- numeric(0L)
    cin.theta <- numeric(0L)
  }
  ceq.linear.flag <- length(ceq.linear.idx) > 0L
  ceq.nonlinear.flag <- length(ceq.nonlinear.idx) > 0L
  ceq.flag <- ceq.linear.flag || ceq.nonlinear.flag
  cin.linear.flag <- length(cin.linear.idx) > 0L
  cin.nonlinear.flag <- length(cin.nonlinear.idx) > 0L
  cin.flag <- cin.linear.flag || cin.nonlinear.flag
  #ceq.only.flag <- ceq.flag && !cin.flag
  cin.only.flag <- cin.flag && !ceq.flag
  ceq.linear.only.flag <- (ceq.linear.flag && !ceq.nonlinear.flag && 
                             !cin.flag)
  ceq.simple.only <- ceq.simple && !ceq.flag && !cin.flag
  if (ceq.linear.flag) {
    QR <- qr(t(ceq.JAC))
    ranK <- QR$rank
    Q <- qr.Q(QR, complete = TRUE)
    ceq.JAC.NULL <- Q[, -seq_len(ranK), drop = FALSE]
    if (all(ceq.rhs == 0)) {
      ceq.rhs.NULL <- numeric(npar)
    }
    else {
      tmp <- qr.coef(QR, diag(npar))
      NA.idx <- which(is.na(rowSums(tmp)))
      if (length(NA.idx) > 0L) {
        tmp[NA.idx, ] <- 0
      }
      ceq.rhs.NULL <- as.numeric(ceq.rhs %*% tmp)
    }
  }
  else {
    ceq.JAC.NULL <- matrix(0, 0L, 0L)
    ceq.rhs.NULL <- numeric(0L)
  }
  ceq.simple.K <- matrix(0, 0, 0)
  if (ceq.simple.only) {
    n.unco <- max(partable$unco)
    n.free <- max(partable$free)
    ceq.simple.K <- matrix(0, nrow = n.unco, ncol = n.free)
    idx.free <- partable$free[partable$free > 0]
    for (k in 1:n.unco) {
      c <- idx.free[k]
      ceq.simple.K[k, c] <- 1
    }
  }
  ceq.jacobian <- function() NULL
  cin.jacobian <- function() NULL
  OUT <- list(def.function = def.function, ceq.function = ceq.function, 
              ceq.JAC = ceq.JAC, ceq.jacobian = ceq.jacobian, ceq.rhs = ceq.rhs, 
              ceq.theta = ceq.theta, ceq.linear.idx = ceq.linear.idx, 
              ceq.nonlinear.idx = ceq.nonlinear.idx, ceq.linear.flag = ceq.linear.flag, 
              ceq.nonlinear.flag = ceq.nonlinear.flag, ceq.flag = ceq.flag, 
              ceq.linear.only.flag = ceq.linear.only.flag, ceq.JAC.NULL = ceq.JAC.NULL, 
              ceq.rhs.NULL = ceq.rhs.NULL, ceq.simple.only = ceq.simple.only, 
              ceq.simple.K = ceq.simple.K, cin.function = cin.function, 
              cin.JAC = cin.JAC, cin.jacobian = cin.jacobian, cin.rhs = cin.rhs, 
              cin.theta = cin.theta, cin.linear.idx = cin.linear.idx, 
              cin.nonlinear.idx = cin.nonlinear.idx, cin.linear.flag = cin.linear.flag, 
              cin.nonlinear.flag = cin.nonlinear.flag, cin.flag = cin.flag, 
              cin.only.flag = cin.only.flag)
  OUT
}
