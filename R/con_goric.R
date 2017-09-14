goric <- function(object, ..., complement = FALSE, bound = NULL, 
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
    
    if (!is.null(bound) && meq == 0L) {
      warning("restriktor WARNING: bounds are only available for equality constraints \n",
              "                      and are therefore ignored.")
       bound <- NULL
    } 
    
    if (!is.null(bound)) {
      Amat.ciq <- Amat[-c(0:meq), , drop = FALSE]
      Amat.ceq <- Amat[0:meq, , drop = FALSE]
      bvec.ciq <- bvec[-c(0:meq)]
      bvec.ceq <- bvec[0:meq]
      
      ub <- bvec.ceq + bound
      lb <- bvec.ceq - bound 
      
      ## check if any equality constraints
      if (meq > 0L) {
        ## checks
        # upper-bound
        if (length(ub) == 1L && meq >= 1L)  {
          ub <- rep(ub, meq)
        } else if (meq != length(ub)) {
          stop("restriktor ERROR: the number of bounds is not equal to neq.")
        } 
        
        # lower-bound
        if (length(lb) == 1L && meq >= 1L)  {
          lb <- rep(lb, meq)
        } else if (meq != length(lb)) {
          stop("restriktor ERROR: the number of bounds is not equal to neq.")
        } 
        
        # check if unconstrained mle are violated
        check.ciq <- all(Amat.ciq %*% c(b.unrestr) - bvec.ciq >= 0) 
        ## check if unconstrained mle lay between the boundaries
        check.ub <- all(Amat.ceq %*% c(b.unrestr) <= ub)         
        check.lb <- all(Amat.ceq %*% c(b.unrestr) >= lb)         
        
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
          
          ## only if ciq involved                                               # add check!
          if (length(Amat.ciq) > 0L) {
            # ciq --> ceq; ceq --> free                            
            ll.free <- list()
            nr.free <- 1:nrow(Amat.ciq)
            for (l in nr.free) {
             Hc.free <- restriktor(model.org, constraints = Amat.ciq[l, , drop = FALSE],
                                   neq = 1, rhs = bvec.ciq[l],
                                   mix.weights = "none", se = "none")
             ll.free[[l]] <- logLik(Hc.free)
            }
          } else {
            ll.free <- NULL
          }
          
          llc <- unlist(c(ll.ub, ll.lb, ll.free))
          ll.Hc <- max(llc)
          
          if (debug) {
            print(llc)
          }  
        } else {
          # determine the log-likelihood_m in case the unconstrained mle 
          # lay outside the range restriktions. 
          # 2^q2 combinations.
          ll.Hc <- logLik(model.org)
          
          ## ll.Hm
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
          
          
          # only needed if ciq > 0
          # add extra rows to contraint matrix for checks
          
          if (length(Amat.ciq) > 0L) {
            Amatx <- rbind(Amat.ciq, Amat.ceq, Amat.ceq)
            bvecx <- c(bvec.ciq, ub, lb)
            nr.ciq <- 1:nrow(Amat.ciq)
            
            
            rbind( c(-.05, -.10, -.15), # lb lb lb 
                   c( .05,  .10,  .15), # ub ub ub
                   c(-.05, -.10,  .15), # lb lb ub
                   c(-.05,  .10,  .15), # lb ub ub
                   c( .05, -.10, -.15), # ub lb lb
                   c( .05,  .10, -.15)) # ub ub lb
           
            test <- rbind(lb, ub) 
            
            
            
            nr.perm <- 1:nrow(perm)
            ll.outer <- ll.inner <- list() 
            for (l in nr.ciq) {
              nr.idx <- c(nr.ciq[l], nr.ciq[-l])
              Amat.ciq.idx <- Amat.ciq[nr.idx, , drop = FALSE]
              Amat.nr <- rbind(Amat.ciq.idx, Amat.ceq)
              perm.idx <- ifelse(perm > 0, -1, 1)
              perm2 <- perm * perm.idx
              idx.ceq <- which(colSums(abs(Amat.ceq)) > 0L)
              for (m in nr.perm) {
                Amat.nr[idx.ceq, idx.ceq] <- Amat.ceq[ , idx.ceq] * perm.idx[m, ]
                bvec.nr <- c(bvec.ciq[nr.idx], perm2[m, ])
                Hm <- restriktor(model.org, constraints = Amat.nr,
                                 neq = 1, rhs = bvec.nr,
                                 mix.weights = "none", se = "none")
                Hm.b <- coef(Hm)
                Hm.b[abs(Hm.b) < sqrt(.Machine$double.eps)] <- 0L
                if (all(Amatx %*% c(Hm.b) - bvecx >= 0L)) {
                  ll.inner[[m]] <- logLik(Hm)
                }
              }
              ll.outer[[l]] <- ll.inner
            }
          } else {
            ll.outer <- NULL
          }
          llm <- unlist(c(ll.ub, ll.lb, ll.outer))
          ll.Hm <- max(llm)
          
          if (debug) {
            print(llm)
          }  
        }
      } else {
        # if meq = 0. This should be catched earlier
        stop("restriktor ERROR: you might have found a bug, please contact me!")
      }
    } else if (is.null(bound)) {
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
          cat("log-likelihood_c value =", ll.Hc, "\n")
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
          Amatx <- Amat[idx, , drop = FALSE]
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
          cat("log-likelihood_c value =", ll.Hc, "\n")
        }
      } else if (all(c(Amat) == 0L)) {
        # unconstrained setting
        stop("Restriktor ERROR: no complement exists for the unconstrained hypothesis.")
      } else {
        stop("Restriktor ERROR: you might have found a bug, please contact me!")
      }
      ll.Hm <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
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
  }
  
  if (!complement) {
    ll    <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
    PT    <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    goric <- unlist(lapply(isSummary, function(x) x$goric[1]))
    df    <- data.frame(model = objectnames, loglik = ll, penalty = PT, goric)
    df$model <- as.character(df$model)
  } else {
    # model
    PTm <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    goric.Hm <- -2*(ll.Hm - PTm)
    df.Hm  <- data.frame(model = objectnames, loglik = ll.Hm, penalty = PTm, 
                         goric = goric.Hm)
    df.Hm$model <- as.character(df.Hm$model)
    # complement
    goric.Hc <- -2*(ll.Hc - PTc)
    df.c  <- data.frame(model = "complement", loglik = ll.Hc, penalty = PTc, 
                        goric = goric.Hc)
    df <- rbind(df.Hm, df.c)
  } 
  
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
