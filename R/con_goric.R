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
    # extract equalities and inequalities
    if (meq > 0) {
      Amat.ceq <- Amat[1:meq, ,drop = FALSE]
      bvec.ceq <- bvec[1:meq]
      Amat.ciq <- Amat[-c(1:meq), , drop = FALSE]
      bvec.ciq <- bvec[-c(1:meq)]
    } else {
      Amat.ceq <- matrix( , nrow = 0, ncol = ncol(Amat))
      bvec.ceq <- rep(0, 0)
      Amat.ciq <- Amat[ , , drop = FALSE]
      bvec.ciq <- bvec
    }
    
    if (!is.null(bound) && meq == 0L) {
      warning("restriktor WARNING: bounds are only available for equality constraints \n",
              "                      and are therefore ignored.")
       bound <- NULL
    } 
    
    if (!is.null(bound)) {
      if (meq > 0) {
        Amat.ceq <- Amat[1:meq, ,drop = FALSE]
        bvec.ceq <- bvec[1:meq]
        Amat.ciq <- Amat[-c(1:meq), , drop = FALSE]
        bvec.ciq <- bvec[-c(1:meq)]
      } else {
        Amat.ceq <- matrix( , nrow = 0, ncol = ncol(Amat))
        bvec.ceq <- rep(0, 0)
        Amat.ciq <- Amat[ , , drop = FALSE]
        bvec.ciq <- bvec
      }
      
      # upper-bound = rhs + bound
      ub <- bvec.ceq + bound
      # lower-bound = rhs - bound
      lb <- bvec.ceq - bound 
      
      ## check if any equality constraints
      if (meq > 0L) {
        # correct user error
        if ( ((length(ub) == 1L) || (length(lb) == 1L)) && meq >= 1L )  {
          ub <- rep(ub, meq)
          lb <- rep(lb, meq)
        } 
        # check 
        if ( (length(ub) != meq) || (length(lb) != meq) ) {
          stop("restriktor ERROR: the number of bounds is not equal to neq.")
        } 
        
        # check if unconstrained mle are violated
        check.ciq <- all(Amat.ciq %*% c(b.unrestr) - bvec.ciq >= 0) 
        ## check if unconstrained mle lay between the bounds
        check.ub <- all(Amat.ceq %*% c(b.unrestr) <= (ub + .Machine$double.eps))
        check.lb <- all(Amat.ceq %*% c(b.unrestr) >= (lb - .Machine$double.eps))
  
        ## check if unrestricted mle lay in boundary area
        if (check.ciq && check.ub && check.lb) {
          # log-likelihood model
          llm <- logLik(model.org)
          
          ## determine log-likelihood_c
          nr.meq <- 1:meq
          # upper-bound
          llc.ub <- list()
          for (l in nr.meq) {
            Amat.ub <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
            bvec.ub <- c(ub[l], bvec.ciq)
            Hc.ub <- restriktor(model.org, constraints = Amat.ub, 
                                neq = 1, rhs = bvec.ub, 
                                mix.weights = "none", se = "none")
            llc.ub[[l]] <- logLik(Hc.ub)
          }
          # lower-bound
          llc.lb <- list()
          for (l in nr.meq) {
            Amat.lb <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
            bvec.lb <- c(lb[l], bvec.ciq)
            Hc.lb <- restriktor(model.org, constraints = Amat.lb, 
                                neq = 1, rhs = bvec.lb, 
                                mix.weights = "none", se = "none")
            llc.lb[[l]] <- logLik(Hc.lb)
          }
          
          ## only if any ciq                                     
          if (length(Amat.ciq) > 0L) { 
            Amatx <- rbind(Amat.ciq, Amat.ceq, -Amat.ceq)        # correct?
            bvecx <- c(bvec.ciq, lb, lb)                         # correct?   
            
            llc.s <- list()
            nr.ciq <- 1:nrow(Amat.ciq)
            for (l in nr.ciq) {
             Hc.s <- restriktor(model.org, constraints = Amat.ciq[l, , drop = FALSE],
                                neq = 1, rhs = bvec.ciq[l],
                                mix.weights = "none", se = "none")
             b.s <- coef(Hc.s)
             b.s[abs(b.s) < sqrt(.Machine$double.eps)] <- 0L
             if (all(Amatx %*% c(b.s) >= bvecx)) {               # correct check?
               llc.s[[l]] <- logLik(Hc.s)
             }
            }
          } else {
            llc.s <- NULL
          }
          
          llc <- unlist(c(llc.ub, llc.lb, llc.s))
          llc <- max(llc)
          
          if (debug) {
            print(llc)
          }  
        } else {
          # determine the log-likelihood_m in case the unconstrained mle 
          # lay outside the range restriktions. 
          # 2^q2 combinations.
          llc <- logLik(model.org) 
          
          # if all bounds are zero, then llm = logLik(object)
          if (!all(bound == 0L)) {
            # get all combinations for ub and lb
            
            len.bvec.ceq <- length(ub)
            perm <- rep(list(rep(c(-1,1), (2^len.bvec.ceq)/2)), len.bvec.ceq)
            perm.grid <- unique(as.matrix(expand.grid(perm)))
            nr.perm <- 1:(2^len.bvec.ceq)
            # if ub = lb, then the bounds must be zero.
            # Otherwise, +ub, -lb
            bound.zero.idx <- which(ub == lb)
            
            llm <- list()
            for (m in nr.perm) {
              #perm.vec <- perm.grid[m, ]
              perm.vec <- apply(Amat.ceq, 2, function(x) perm.grid[m, ] * x)
              ub.idx <- which(perm.grid[m, ] == -1)
              lb.idx <- which(perm.grid[m, ] ==  1)
              order.idx <- c(ub.idx, lb.idx)
              bounds.new <- c(-ub[ub.idx], lb[lb.idx])
              bounds.new <- bounds.new[order(order.idx, decreasing = FALSE)]
              
              Amatx <- rbind(Amat.ciq, Amat.ceq, -Amat.ceq)
              bvecx <- c(bvec.ciq, bounds.new, bounds.new)
              
              #Amat.ceq.perm <- t( t(Amat.ceq) * perm.vec ) 
              Amat.ceq.perm <- perm.vec
              Amat.new <- rbind(Amat.ceq.perm, Amat.ciq)
              bvec.new <- c(bounds.new, bvec.ciq)
            
              if (length(bound.zero.idx) == 0L) {
                Amat.new.sort <- Amat.new
                bvec.new.sort <- bvec.new
              } else {
                # constraint rows for which bound == 0 are placed on top
                Amat.new.sort <- rbind(Amat.new[bound.zero.idx, , drop = FALSE],
                                       Amat.new[-bound.zero.idx, , drop = FALSE])
                bvec.new.sort <- c(bvec.new[bound.zero.idx],
                                   bvec.new[-bound.zero.idx])    
              }
              
              Hm <- restriktor(model.org, constraints = Amat.new.sort,
                               neq = length(bound.zero.idx), rhs = bvec.new.sort,
                               mix.weights = "none", se = "none")
              beta.Hm <- coef(Hm)
              beta.Hm[abs(beta.Hm) < sqrt(.Machine$double.eps)] <- 0L
              if (all( Amatx %*% c(beta.Hm) - bvecx + .Machine$double.eps >= 0 )) {
                llm[[m]] <- logLik(Hm)
              } 
            }
            llm <- unlist(llm)
            if (debug) {
              cat("log-likelihood_m =", llm, "\n")
            }
  
            llm <- max(llm)
          } else {
            llm <- logLik(object)
          }
          llm
        }
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
        llc <- con_loglik_lm(model.org)
        if (debug) {
          cat("log-likelihood_c value =", llc, "\n")
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
        llc <- max(unlist(ll))
      } else if (nrow(Amat) == meq) { 
        # redundant, this will be catched by the first statement.
        # in case of equality constraints only, the complement is 
        # equal to the unconstrained log-likelihood
        llc <- con_loglik_lm(model.org)
        if (debug) {
          cat("log-likelihood_c value =", llc, "\n")
        }
      } else if (all(c(Amat) == 0L)) {
        # unconstrained setting
        stop("Restriktor ERROR: no complement exists for the unconstrained hypothesis.")
      } else {
        stop("Restriktor ERROR: you might have found a bug, please contact me!")
      }
      llm <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
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
    goric.Hm <- -2*(llm - PTm)
    df.Hm  <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                         goric = goric.Hm)
    df.Hm$model <- as.character(df.Hm$model)
    # complement
    goric.Hc <- -2*(llc - PTc)
    df.c  <- data.frame(model = "complement", loglik = llc, penalty = PTc, 
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
