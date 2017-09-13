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
    
    if (!is.null(lower) && !is.null(upper) && meq == 0L) {
      warning("restriktor WARNING: bounds are only available for equality constraints \n",
              "                      and are therefore ignored.")
      upper <- lower <- NULL
    } else if (!is.null(upper) && is.null(lower)) {
      stop("restriktor ERROR: no lower-bound is specified.")
    } else if (is.null(upper) && !is.null(lower)) {
      stop("restriktor ERROR: no upper-bound is specified.")
    }
    
    if (!is.null(upper) && !is.null(lower)) {
      
      Amat.ciq <- Amat[-c(0:meq), , drop = FALSE]
      Amat.ceq <- Amat[0:meq, , drop = FALSE]
      bvec.ciq <- bvec[-c(0:meq)]
      bvec.ceq <- bvec[0:meq]
      
      ub <- bvec.ceq + upper
      lb <- bvec.ceq + lower 
      
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
          
          ## only if ciq involved  
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
          
          
          Amatx <- rbind(Amat.ciq, Amat.ceq, Amat.ceq)
          bvecx <- c(bvec.ciq, ub, lb)
          
          #### 
          # ciq --> ceq; ceq --> ciq       
          nrx <- 1:nrow(Amat.ciq)
          perm <- gtools:::permutations(n = length(ub), r = length(ub), v = c(lb, ub), repeats.allowed = TRUE)
          nrm <- 1:nrow(perm)
          llx2.ub <- llx.ub <- list() 
          for (l in nrx) {
            idx <- c(nrx[l], nrx[-l])
            Amat.ciq.idx <- Amat.ciq[idx, , drop = FALSE]
            Amat.ub <- rbind(Amat.ciq.idx, Amat.ceq)
            
            idx2 <- ifelse(perm > 0, -1, 1)
            perm2 <- perm * idx2
            idx.ceq <- which(colSums(Amat.ceq) > 0L)
            for (m in nrm) {
              Amat.ub[idx.ceq, idx.ceq] <- Amat.ceq[, idx.ceq] * idx2[m,]
              bvec.ub <- c(bvec.ciq[idx], perm2[m, ])
              Hcx.ub <- restriktor(model.org, constraints = Amat.ub,
                                   neq = 1, rhs = bvec.ub,
                                   mix.weights = "none", se = "none")
              
              b.ub <- coef(Hcx.ub)
              b.ub[abs(b.ub) < sqrt(.Machine$double.eps)] <- 0L
              if (all(Amatx %*% c(b.ub) - bvecx >= 0L)) {
                llx.ub[[m]] <- logLik(Hcx.ub)
              }
            }
            llx2.ub[[l]] <- llx.ub
          }
          llx <- unlist(llx2.ub)
      
          llm <- unlist(c(ll.ub, ll.lb, llx))
          ll.Hm <- max(llm)
          
          if (debug) {
            print(llm)
          }  
        }
      } else {
        # if meq = 0. This should be catched earlier
        stop("restriktor ERROR: you might have found a bug, please contact me!")
      }
    } else if (is.null(upper) && is.null(lower)) {
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
