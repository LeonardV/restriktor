compute_complement_likelihood <- function(model.org, VCOV, 
                                          Amat, Amat.ciq, Amat.ceq, 
                                          bvec, bvec.ciq, bvec.ceq, 
                                          meq, b.unrestr, type, ldots,
                                          debug = FALSE) {
  
  # construeer zelf Amat.ceq en Amat.ciq obv Amat en meq.
  
  # check if any equality constraint is violated
  check.ceq <- !(all(c(Amat.ceq %*% c(b.unrestr)) - bvec.ceq == 0))
  if (nrow(Amat) > meq) {
    # check if any inequality constraint is violated
    check.ciq <- !(all(c(Amat.ciq %*% c(b.unrestr)) - bvec.ciq >= 0))
  } else {
    check.ciq <- FALSE
  }
  
  # number of parameters
  p <- length(b.unrestr)
  
  if (check.ciq || check.ceq) {    
    if (type %in% c("goric", "goricc")) { 
      llc <- logLik(model.org)
      betasc <- b.unrestr
    } else if (type %in% c("gorica", "goricac")) {
      llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
      betasc <- b.unrestr
    }
    if (debug) {
      cat("log-likelihood_c value =", llc, "\n")
    }
    # if any constraints are violated LL_c = LL_u
  } else if (nrow(Amat) > meq && !(all(c(Amat) == 0L))) {
    nr <- seq_len(nrow(Amat))
    ll <- vector("list", length(nr))
    betas <- vector("list", length(nr))
    if (meq > 0L) { nr <- nr[-c(0:meq)] }
      for (l in seq_len(length(nr))) {
        idx <- c(nr[l], nr[-l])
        Amatx <- Amat[idx, , drop = FALSE]
        if (type %in% c("goric", "goricc")) {          
          Hc.restr <- restriktor(model.org, constraints = Amatx, 
                                 neq = 1, rhs = bvec[idx], 
                                 mix_weights = "none", se = "none")
          betas[[l]] <- coef(Hc.restr)
          ll[[l]]    <- logLik(Hc.restr)
        } else if (type %in% c("gorica", "goricac")) {
          ldots$mix_weights <- "none"
          CALL.restr <- append(list(object      = b.unrestr,
                                    constraints = Amatx,
                                    rhs         = bvec[idx],
                                    neq         = 1,
                                    VCOV        = VCOV),
                               ldots)
          Hc.restr   <- do.call("con_gorica_est", CALL.restr) 
          betas[[l]] <- Hc.restr$b.restr
          ll[[l]]    <- dmvnorm(c(b.unrestr - Hc.restr$b.restr), 
                                sigma = VCOV, log = TRUE)            
        }
      }
    if (debug) {
      cat("log-likelihood value =", ll[[l]], "\n")
    }
    ll.unlist <- unlist(ll)
    ll.idx <- which.max(ll.unlist)
    llc <- max(ll.unlist)
    betasc <- betas[[ll.idx]]
  } else if (nrow(Amat) == meq) { 
    if (type %in% c("goric", "goricc")) {
      llc <- logLik(model.org)
      betasc <- b.unrestr
    } else if (type %in% c("gorica", "goricac")) {
      llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
      betasc <- b.unrestr
    }
    if (debug) {
      cat("log-likelihood_c value =", llc, "\n")
    }
  } else if (all(c(Amat) == 0L)) {
    stop("restriktor ERROR: no complement exists for an unconstrained hypothesis.")
  } else {
    stop("restriktor ERROR: you might have found a bug, please contact me at: info@restriktor.org!")
  }
  
  return(list(llc = llc, betasc = betasc))
}
