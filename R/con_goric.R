goric <- function(object, ..., complement = FALSE, lower = NULL, upper = NULL,
                  digits = max(3, getOption("digits") - 2), debug = FALSE) {
  
  mc <- match.call()
  if (inherits(object, "restriktor")) {
    objectlist <- list(object, ...)
  } else {
    objectlist <- object
  }
  
  isrestr <- sapply(objectlist, function(x) inherits(x, "restriktor"))
  conlist <- objectlist[isrestr]  
  isSummary <- lapply(conlist, function(x) summary(x))
  
  CALL <- as.list(mc)
  CALL[[1]] <- NULL
  
  idx <- which(isrestr)
  objectnames <- vector("character", length(idx))
  for (i in idx) { 
    objectnames[i] <- as.character(CALL[[i]])
  }
  
  wt.check <- unlist(lapply(isSummary, function(x) is.null(x$goric)))
  wt.check.idx <- which(wt.check)
  if (length(wt.check.idx) > 0L) {
    stop("Restriktor ERROR: for model ", wt.check.idx, " no GORIC value is available.")
  }
  
  if (complement && length(conlist) > 1L) {
    complement <- FALSE
    warning("Restriktor WARNING: if complement = TRUE, only one order-constrained hypothesis\n",
            "                      is allowed. Therefore, the complement is set to FALSE.")
  } 
  
  df.c <- NULL
  if (complement) {
    # original unconstrained model
    model.org <- object$model.org
    # unrestricted parameters
    b.unrestr <- object$b.unrestr
    # level probabilities
    wt.bar <- object$wt.bar
    # constraints matrix
    Amat <- object$constraints
    # number of equalities
    meq <- object$neq
    # rhs
    bvec <- object$rhs
    # extract equalities
    Amat.ceq <- Amat[0:meq, ,drop = FALSE]
    bvec.ceq <- bvec[0:meq]
    # extract inequalities
    Amat.ciq <- Amat[-c(0:meq), , drop = FALSE]
    bvec.ciq <- bvec[-c(0:meq)]
    
    ub <- upper
    lb <- lower
    
    if (!is.null(ub) && !is.null(lb) && meq == 0L) {
      warning("restriktor WARNING: bounds are only available for equality constraints \n",
              "                      and are therefore ignored.")
      ub <- lb <- NULL
    } else if (!is.null(ub) && is.null(lb)) {
      stop("restriktor ERROR: no lower-bound is specified.")
    } else if (is.null(ub) && !is.null(lb)) {
      stop("restriktor ERROR: no upper-bound is specified.")
    }
    
    if (!is.null(ub) && !is.null(lb)) {
      ## check if any equality constraints
      if (meq > 0L) {
        ## checks
        # upper-bound
        if (length(ub) == 1L && meq >= 1L)  {
          ub <- rep(ub, meq)
        } else if (meq != length(ub)) {
          stop("restriktor ERROR: the number of upper-bounds is not equal to neq.")
        } else {
          stop("restriktor ERROR: you probably found a bug, please contact me!")
        }
        # lower-bound
        if (length(lb) == 1L && meq >= 1L)  {
          lb <- rep(lb, meq)
        } else if (meq != length(lb)) {
          stop("restriktor ERROR: the number of lower-bounds is not equal to neq.")
        } else {
          stop("restriktor ERROR: you might have found a bug, please contact me!")
        }
        
        Amat.ciq <- Amat[-c(0:meq), , drop = FALSE]
        Amat.ceq <- Amat[0:meq, , drop = FALSE]
        bvec.ciq <- bvec[-c(0:meq)]
        bvec.ceq <- bvec[0:meq]
        
        # check if unconstrained mle are violated
        check.ciq <- all(Amat.ciq %*% c(b.unrestr) - bvec.ciq >= 0) 
        ## check if unconstrained mle lay between the boundaries
        # upper
        check.ub <- all(Amat.ceq %*% c(b.unrestr) <= ub)         # check if correct!
        # lower
        check.lb <- all(Amat.ceq %*% c(b.unrestr) >= lb)         # check if correct!
        
        ## check if unrestricted mle lay in boundary area
        if (check.ciq && check.ub && check.lb) {
          # log-likelihood model
          ll.Hm <- logLik(model.org)
          ## determine log-likelihood_c
          ## 2q2 + 1 combinations
          
          # upper-bound
          ll.ub <- list()
          nr <- 1:meq
          for (l in nr) {
            Amat.ub <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
            bvec.ub <- c(ub[l], bvec.ciq)
            Hc.ub <- restriktor(model.org, constraints = Amat.ub, 
                                neq = 1, rhs = bvec.ub, 
                                mix.weights = "none", se = "none")
            ll.ub[[l]] <- logLik(Hc.ub)
          }
          
          # lower-bound
          ll.lb <- list()
          for (l in nr) {
            Amat.lb <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
            bvec.lb <- c(lb[l], bvec.ciq)
            Hc.lb <- restriktor(model.org, constraints = Amat.lb, 
                                neq = 1, rhs = bvec.lb, 
                                mix.weights = "none", se = "none")
            ll.lb[[l]] <- logLik(Hc.lb)
          }
          
          # ciq --> ceq; ceq --> free                            # all ciq at once or one-by-one?
          # e.g., b1 > 0, b2 > 0, b3 = 0; b4 = 0
          # ll.free <- list()
          # nr.free <- 1:nrow(Amat.ciq)
          # for (l in nr.free) {
          #   Hc.free <- restriktor(model.org, constraints = Amat.ciq[l, , drop = FALSE],
          #                         neq = 1, rhs = bvec.ciq[l],
          #                         mix.weights = "none", se = "none")
             
            Hc.free <- restriktor(model.org, constraints = Amat.ciq,
                                  neq = nrow(Amat.ciq), rhs = bvec.ciq,
                                  mix.weights = "none", se = "none")
            ll.free <- logLik(Hc.free)
          #   ll.free[[l]] <- logLik(Hc.free)
          # }
          
          llc <- unlist(c(ll.ub, ll.lb, ll.free))
          ll.Hc <- max(llc)
        } else {
          # determine the log-likelihood_m in case the unconstrained mle 
          # lay outside the range restriktions. 
          # 2^q2 combinations.
          ll.Hc <- logLik(model.org)
          
          # upper-bound
          ll.ub <- list()
          nr <- 1:meq
          for (l in nr) {
            Amat.ub <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
            bvec.ub <- c(ub[l], bvec.ciq)
            Hc.ub <- restriktor(model.org, constraints = Amat.ub, 
                                neq = 1, rhs = bvec.ub, 
                                mix.weights = "none", se = "none")
            ll.ub[[l]] <- logLik(Hc.ub)
          }
          
          # lower-bound
          ll.lb <- list()
          for (l in nr) {
            Amat.lb <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
            bvec.lb <- c(lb[l], bvec.ciq)
            Hc.lb <- restriktor(model.org, constraints = Amat.lb, 
                                neq = 1, rhs = bvec.lb, 
                                mix.weights = "none", se = "none")
            ll.lb[[l]] <- logLik(Hc.lb)
          }
          
          # ciq --> ceq; ceq --> ciq                             # all ciq at once or one-by-one?
          # llx.ub <- list()
          # for (l in nr) {
          #   Amat.ub <- rbind(Amat.ciq, -1*Amat.ceq[l, , drop = FALSE])
          #   bvec.ub <- c(bvec.ciq, -1*ub[l])
          #   Hcx.ub <- restriktor(model.org, constraints = Amat.ub, 
          #                        neq = nrow(Amat.ciq), rhs = bvec.ub, 
          #                        mix.weights = "none", se = "none")
          #   idx.ub <- which(colSums(abs(Amat.ceq)) > 0L)
          #   ub.restr <- coef(Hcx.ub)[idx.ub]
          #   if (all(ub.restr >= lb)) {
          #     llx.ub[[l]] <- logLik(Hcx.ub)
          #   }
          # }

          Amat.ub <- rbind(Amat.ciq, -1*Amat.ceq)
          bvec.ub <- c(bvec.ciq, -1*ub)
          Hc.ub <- restriktor(model.org, constraints = Amat.ub,
                              neq = nrow(Amat.ciq), rhs = bvec.ub,
                              mix.weights = "none", se = "none")
          idx.ub <- which(colSums(abs(Amat.ceq)) > 0L)
          ub.restr <- coef(Hc.ub)[idx.ub]
          if (all(ub.restr >= lb)) {
            llx.ub <- logLik(Hc.ub)
          } else {
            llx.ub <- NULL
          }
          
          # llx.lb <- list()
          # for (l in nr) {
          #   Amat.lb <- rbind(Amat.ciq, 1*Amat.ceq[l, , drop = FALSE])
          #   bvec.lb <- c(bvec.ciq, 1*lb[l])
          #   Hcx.lb <- restriktor(model.org, constraints = Amat.lb, 
          #                        neq = nrow(Amat.ciq), rhs = bvec.lb, 
          #                        mix.weights = "none", se = "none")
          #   
          #   idx.lb <- which(colSums(abs(Amat.ceq)) > 0L)
          #   lb.restr <- coef(Hcx.lb)[idx.lb]
          #   if (all(lb.restr <= ub)) {
          #     llx.lb[[l]] <- logLik(Hcx.lb)
          #   }
          # }

          Amat.lb <- rbind(Amat.ciq, 1*Amat.ceq)
          bvec.lb <- c(bvec.ciq, 1*lb)
          Hc.lb <- restriktor(model.org, constraints = Amat.ub,
                              neq = nrow(Amat.ciq), rhs = bvec.ub,
                              mix.weights = "none", se = "none")
          idx.lb <- which(colSums(abs(Amat.ceq)) > 0L)
          lb.restr <- coef(Hc.lb)[idx.lb]
          if (all(lb.restr <= ub)) {
            llx.lb <- logLik(Hc.lb)
          } else {
            llx.lb <- NULL
          }
          
          llm <- unlist(c(ll.ub, ll.lb, llx.ub, llx.lb))
          ll.Hm <- max(llm)
        }
      } else {
        # if meq = 0. This should be catched earlier
        stop("restriktor ERROR: you might have found a bug, please contact me!")
      }
    } else if (is.null(ub) && is.null(lb)) {
      # check if any equality constraint is violated
      check.ceq <- !(all(Amat.ceq %*% c(b.unrestr) - bvec.ceq == 0))
      if (nrow(Amat) > meq) {
        # check if any inequality constraint is violated
        check.ciq <- !(all(Amat.ciq %*% c(b.unrestr) - bvec.ciq >= 0))
      } else {
        check.ciq <- FALSE
      }
      ## compute log-likelihood for complement
      # check if any constraint is violated
      if (check.ciq || check.ceq) { 
        ll.Hc <- con_loglik_lm(model.org)
        if (debug) {
          cat("log-likelihood value =", ll.Hc, "\n")
        }
        # if any constraints is violated LL_c = LL_u
      } else if (nrow(Amat) > meq && !(all(c(Amat) == 0L))) {
        ll <- list()
        # number of rows
        nr <- 1:nrow(Amat)
        # remove rows corresponding to equality constraints
        if (meq > 0L) {
          nr <- nr[-c(0:meq)]
        }
        # treat each row of Amat as an equality constraint
        for (l in 1:length(nr)) {
          idx <- c(nr[l], nr[-l])
          Amatx <- Amat[idx,,drop = FALSE]
          Hc.restr <- restriktor(model.org, constraints = Amatx, 
                                 neq = 1, rhs = bvec[idx], 
                                 mix.weights = "none", se = "none")
          ll[[l]] <- Hc.restr$loglik
        }
        if (debug) {
          cat("log-likelihood value =", ll[[l]], "\n")
        }
        # take the highest log-likelihood value as a substitute for 
        # the complement
        ll.Hc <- max(unlist(ll))
      } else if (nrow(Amat) == meq) { 
        # redundant, this will be catched by the first statement.
        # in case of equality constraints only, the complement is 
        # equal to the unconstrained log-likelihood
        ll.Hc <- con_loglik_lm(model.org)
        if (debug) {
          cat("log-likelihood value =", ll.Hc, "\n")
        }
      } else if (all(c(Amat) == 0L)) {
        # unconstrained setting
        stop("Restriktor ERROR: no complement exists for the unconstrained hypothesis.")
      } else {
        stop("Restriktor ERROR: you might have found a bug, please contact me!")
      }
    } 
    
    # compute free parameters f in the complement
    p <- length(b.unrestr)
    # nrow q1
    lq1 <- nrow(Amat.ciq)
    # nrow q2
    lq2 <- nrow(Amat.ceq)
    # free parameters. Note that Amat includes q1 and q2 constraints
    f <- p - nrow(Amat) # p - q1 - q2
    if (debug) { cat("number of free parameters =", (f + lq2), "\n") }
    idx <- length(wt.bar)
    # compute penalty term value PTc
    if (attr(wt.bar, "method") == "boot") {
      PTc <- as.numeric(1 + p - wt.bar[idx-meq] * lq1)
    } else if (attr(wt.bar, "method") == "pmvnorm") {
      # the q2 equalities are not included in wt.bar. Hence, do not have to be subtracted.
      PTc <- as.numeric(1 + p - wt.bar[idx] * lq1)
    } else {
      stop("Restriktor ERROR: no level probabilities (chi-bar-square weights) found.")
    }
    
    if (debug) {
      cat("penalty term value =", PTc, "\n")
    }
    # compute goric.c
    goric.Hc <- -2*(ll.Hc - PTc)
    df.c  <- data.frame(model = "complement", loglik = ll.Hc, penalty = PTc, 
                        goric = goric.Hc)
  }
  
  ll    <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
  PT    <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
  goric <- unlist(lapply(isSummary, function(x) x$goric[1]))
  df    <- data.frame(model = objectnames, loglik = ll, penalty = PT, goric)
  df$model <- as.character(df$model)
  df <- rbind(df, df.c)
  
  # compute evidence ratios
  delta <- df$goric - min(df$goric)
  goric.weights <- exp(-delta / 2) / sum(exp(-delta / 2))
  df$goric.weights <- goric.weights
  
  attr(df, "complement") <- complement
  
  class(df) <- "goric"
  
  df
}


print.goric <- function(x, digits = max(3, getOption("digits") - 2), ...) {
  
  dig <- paste0("%6.", digits, "f")
  x2 <- x[-1]
  df <- as.data.frame(lapply(x2, sprintf, fmt = dig))
  df <- cbind(model = x$model, df)
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  
  complement <- attr(x, "complement")
  if (complement) {
    objectnames <- df$model[1]
    goric.weights <- as.numeric(as.character(df$goric.weights[1:2]))
    Gm <- c(goric.weights[1])
    Gc <- goric.weights[2]
    cat("The order-restricted hypothesis", sQuote(objectnames[1]), "is", 
        sprintf("%1.3f", Gm/Gc), "times more likely than its complement.")
  }
}
