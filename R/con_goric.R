goric <- function(object, ..., comparison = c("unrestricted", "complement", "none"), 
                  type = "standard", bound = NULL, 
                  digits = max(3, getOption("digits") - 2), debug = FALSE) {
  
  mc <- match.call()
  CALL <- as.list(mc)
  CALL[[1]] <- NULL
  
  comparison <- match.arg(comparison)
  stopifnot(comparison %in% c("unrestricted", "complement", "none"))
  
  # restriktor objects
  if (inherits(object, "restriktor")) {
    objectList <- list(object, ...)
    isRestr <- sapply(objectList, function(x) inherits(x, "restriktor"))
    if (!any(isRestr)) {
      stop("Restriktor ERROR: all models must be of class restriktor.") 
    }
    conList <- objectList[isRestr]  
    isSummary <- lapply(conList, function(x) summary(x))
    
    idx <- which(isRestr)
    objectnames <- vector("character", length(idx))
    for (i in idx) { 
      objectnames[i] <- as.character(CALL[[i]])
    }
    # check if goric value exists in restriktor object
    wt.check <- unlist(lapply(isSummary, function(x) is.null(x$goric)))
    wt.check.idx <- which(wt.check)
    if (length(wt.check.idx) > 0L) {
      stop("Restriktor ERROR: for model ", wt.check.idx, " no GORIC value is available.")
    }
  } else {
    stop("Restriktor ERROR: object must be of class restriktor.") 
  }

  stopifnot(type %in% c("standard", "goric", "gorica"))
  
  if (type == "goric") {
    type <- "standard"
  }

  if (comparison == "complement" && length(conList) > 1L) {
    comparison <- "none"
    warning("Restriktor WARNING: if comparison = 'complement', only one order-restricted hypothesis\n",
            "                      is allowed (for now). Therefore, comparison is set to 'none'.")
  } 
  
  ans <- list()
  # original unconstrained model
  ans$model.org <- object$model.org
  # unrestricted VCOV
  VCOV <- vcov(object$model.org)

  df.c <- NULL
  if (comparison == "complement") {
      # unrestricted estimates
      b.unrestr <- coef(ans$model.org)
      # restricted estimates
      b.restr <- object$b.restr
      # number of parameters
      p <- length(b.unrestr)
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
      Amat.ceq <- Amat[1:meq, , drop = FALSE]
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
      warning("restriktor WARNING: bounds are only available for equality restrictions \n",
              "                      and are therefore ignored.")
      bound <- NULL
    } 
    
    # for now
    bound <- NULL
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
        if ( ((length(ub) == 1L) || (length(lb) == 1L)) && meq >= 1L ) {
          ub <- rep(ub, meq)
          lb <- rep(lb, meq)
        }
        # check
        if ( (length(ub) != meq) || (length(lb) != meq) ) {
          stop("restriktor ERROR: the number of bounds is not equal to number of equality constraints (neq).")
        }

        # check if unconstrained mle are violated
        check.ciq <- all(Amat.ciq %*% c(b.unrestr) - bvec.ciq >= 0)
        ## check if unconstrained mle lay between the bounds
        check.ub <- all(Amat.ceq %*% c(b.unrestr) <= (ub + .Machine$double.eps))
        check.lb <- all(Amat.ceq %*% c(b.unrestr) >= (lb - .Machine$double.eps))

        ## check if unrestricted mle lay in boundary area
        if (check.ciq && check.ub && check.lb) {
          # log-likelihood model
          if (type == "standard") {
            llm <- logLik(ans$model.org)
          } else if (type == "gorica") {
            llm <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
          }
          ## determine log-likelihood_c
          nr.meq <- 1:meq
          # upper-bound
          llc.ub <- list()
          for (l in nr.meq) {
            Amat.ub <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
            bvec.ub <- c(ub[l], bvec.ciq)
            Hc.ub <- restriktor(ans$model.org, constraints = Amat.ub,
                                neq = 1, rhs = bvec.ub,
                                mix.weights = "none", se = "none")
            if (type == "standard") {  
              llc.ub[[l]] <- logLik(Hc.ub)
            } else if (type == "gorica") {
              llc.ub[[l]] <- dmvnorm(c(b.unrestr - Hc.ub$b.restr), 
                                     sigma = VCOV, log = TRUE)
            }
          }
          # lower-bound
          llc.lb <- list()
          for (l in nr.meq) {
            Amat.lb <- rbind(Amat.ceq[l, , drop = FALSE], Amat.ciq)
            bvec.lb <- c(lb[l], bvec.ciq)
              Hc.lb <- restriktor(ans$model.org, constraints = Amat.lb,
                                  neq = 1, rhs = bvec.lb,
                                  mix.weights = "none", se = "none")
            if (type == "standard") {
              llc.lb[[l]] <- logLik(Hc.lb)
            } else if (type == "gorica") {
              llc.lb[[l]] <- dmvnorm(c(b.unrestr - Hc.lb$b.restr), 
                                     sigma = VCOV, log = TRUE)
            }
          }

          ## only if any ciq
          if (length(Amat.ciq) > 0L) {
            Amatx <- rbind(Amat.ciq, Amat.ceq, -Amat.ceq)        # correct?
            bvecx <- c(bvec.ciq, lb, lb)                         # correct?

            llc.s <- list()
            nr.ciq <- 1:nrow(Amat.ciq)
            for (l in nr.ciq) {
              Hc.s <- restriktor(ans$model.org, 
                                 constraints = Amat.ciq[l, , drop = FALSE],
                                 neq = 1, rhs = bvec.ciq[l],
                                 mix.weights = "none", se = "none")
              b.s <- coef(Hc.s)
              b.s[abs(b.s) < sqrt(.Machine$double.eps)] <- 0L
              if (all(Amatx %*% c(b.s) >= bvecx)) {               # correct check?
                # log-likelihood model
                if (type == "standard") {
                  llc.s[[l]] <- logLik(Hc.s)
                } else if (type == "gorica") {
                  llc.s[[l]] <- dmvnorm(c(b.unrestr - Hc.s$b.restr), 
                                        sigma = VCOV, log = TRUE)
                }
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
          # lay outside the range restrictions.
          # 2^q2 combinations.
          if (type == "standard") {
            llc <- logLik(ans$model.org)
          } else if (type == "gorica") {
            llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
          }

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

              Hm <- restriktor(ans$model.org, constraints = Amat.new.sort,
                               neq = length(bound.zero.idx), rhs = bvec.new.sort,
                               mix.weights = "none", se = "none")
              beta.Hm <- coef(Hm)
              beta.Hm[abs(beta.Hm) < sqrt(.Machine$double.eps)] <- 0L
              if (all( Amatx %*% c(beta.Hm) - bvecx + sqrt(.Machine$double.eps) >= 0 )) {
                if (type == "standard") {
                  llm[[m]] <- logLik(Hm)
                } else if (type == "gorica") {
                  llm[[m]] <- dmvnorm(c(Hm$b.unrestr - Hm$b.restr),
                                      sigma = VCOV, log = TRUE)
                }
              }
            }
            llm <- unlist(llm)
            if (debug) {
              cat("log-likelihood_m =", llm, "\n")
            }

            llm <- max(llm)
          } else {
            if (type == "standard") {
              llm <- logLik(object)
            } else if (type == "gorica") {
              llm <- dmvnorm(c(b.unrestr - b.restr),
                             sigma = VCOV, log = TRUE)
            }
          }
          llm
        }
      }
#################################### no bounds ################################    
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
        #llc <- con_loglik_lm(ans$model.org)
        if (type == "standard") {
          llc <- logLik(ans$model.org)
          betasc <- b.unrestr
        } else if (type == "gorica") {
          llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
          betasc <- b.unrestr
        }

        if (debug) {
          cat("log-likelihood_c value =", llc, "\n")
        }
        # if any constraints is violated LL_c = LL_u
      } else if (nrow(Amat) > meq && !(all(c(Amat) == 0L))) {
        ll <- list()
        betas <- list()
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
          Hc.restr <- restriktor(ans$model.org, constraints = Amatx, 
                                 neq = 1, rhs = bvec[idx], 
                                 mix.weights = "none", se = "none")
          betas[[l]] <- coef(Hc.restr)
          if (type == "standard") {
            ll[[l]] <- logLik(Hc.restr)
          } else if (type == "gorica") {
            ll[[l]] <- dmvnorm(c(b.unrestr - Hc.restr$b.restr), sigma = VCOV, log = TRUE)
          }
        }
        if (debug) {
          cat("log-likelihood value =", ll[[l]], "\n")
        }
        # take the highest log-likelihood value as a substitute for 
        # the complement
        ll.unlist <- unlist(ll)
        ll.idx <- which(ll.unlist == max(ll.unlist))
        llc <- max(ll.unlist)
        betasc <- betas[[ll.idx]]#data.frame(matrix(unlist(betas), nrow=length(betas), byrow=TRUE))
      } else if (nrow(Amat) == meq) { 
        # redundant, this will be catched by the first statement.
        # in case of equality constraints only, the complement is 
        # equal to the unconstrained log-likelihood
        if (type == "standard") {
          llc <- logLik(ans$model.org)
          betasc <- b.unrestr
        } else if (type == "gorica") {
          llc <- dmvnorm(rep(0, p), sigma = VCOV, log = TRUE)
          betasc <- b.unrestr
        }
        if (debug) {
          cat("log-likelihood_c value =", llc, "\n")
        }
      } else if (all(c(Amat) == 0L)) {
        # unconstrained setting
        stop("Restriktor ERROR: no complement exists for the unrestricted hypothesis.")
      } else {
        stop("Restriktor ERROR: you might have found a bug, please contact me!")
      }
      #llm <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
      if (type == "standard") {
        llm <- logLik(object)
      } else if (type == "gorica") {
        llm <- dmvnorm(c(b.unrestr - b.restr), sigma = VCOV, log = TRUE)
      }
    } 
    
    # compute the number of free parameters f in the complement
    p <- ncol(VCOV)#length(b.unrestr)
    # rank q1
    lq1 <- qr(Amat.ciq)$rank
    # rank q2
    lq2 <- qr(Amat.ceq)$rank
    # free parameters. Note that Amat includes q1 and q2 constraints
    f <- p - qr(Amat)$rank # p - q1 - q2
    if (debug) { cat("number of free parameters =", (f + lq2), "\n") }
    # compute penalty term value PTc
    if (type == "standard") {
      idx <- length(wt.bar)
      if (attr(wt.bar, "method") == "boot") {
        PTc <- as.numeric(1 + p - wt.bar[idx-meq] * lq1)
      } else if (attr(wt.bar, "method") == "pmvnorm") {
        # the q2 equalities are not included in wt.bar. Hence, do not have to be subtracted.
        PTc <- as.numeric(1 + p - wt.bar[idx] * lq1)
      } else {
        stop("Restriktor ERROR: no level probabilities (chi-bar-square weights) found.")
      }
    } else if (type == "gorica") {
      wt.bar <- con_weights(cov = Amat %*% VCOV %*% t(Amat), meq = meq)
      idx <- length(wt.bar)
      PTc <- as.numeric(1 + p - wt.bar[idx] * lq1)
    }
    
    if (debug) {
      cat("penalty term value =", PTc, "\n")
    }
  } 
    
  if (comparison == "unrestricted") {
    PTm <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    PTu <- 1 + ncol(VCOV)
    if (type == "standard") {
      # model
      llm <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
      # unrestricted
      llu <- logLik(ans$model.org)
    } else if (type == "gorica") {
      # model
      llm <- unlist(lapply(objectList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
                                                    sigma = VCOV, log = TRUE) )) 
      # unrestricted
      llu <- dmvnorm(rep(0, ncol(VCOV)), sigma = VCOV, log = TRUE) 
    }
    goric.Hm <- -2*(llm - PTm)
    goric.Hu <- -2*(llu - PTu)
    df.Hm <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                        goric = goric.Hm)
    df.Hm$model <- as.character(df.Hm$model)
    df.u <- data.frame(model = "unrestricted", loglik = llu, penalty = PTu, 
                       goric = goric.Hu)
    df <- rbind(df.Hm, df.u)
  } else if (comparison == "none") {
    PT <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    if (type == "standard") {
      ll    <- unlist(lapply(isSummary, function(x) attr(x$goric, "loglik"))) 
      goric.Hm <- unlist(lapply(isSummary, function(x) x$goric[1]))
      df <- data.frame(model = objectnames, loglik = ll, penalty = PT, 
                       goric = goric.Hm)
      df$model <- as.character(df$model)
    } else if (type == "gorica") {
      ll <- unlist(lapply(objectList, function(x) dmvnorm(c(x$b.unrestr - x$b.restr), 
                                                          sigma = VCOV, 
                                                          log = TRUE) )) 
      goric.Hm <- -2*(ll - PT)
      df <- data.frame(model = objectnames, loglik = ll, penalty = PT, 
                       gorica = goric.Hm)
      df$model <- as.character(df$model)
    }
  } else if (comparison == "complement") {
    PTm <- unlist(lapply(isSummary, function(x) attr(x$goric, "penalty")))
    if (type == "standard") {
      # model
      goric.Hm <- -2*(llm - PTm)
      df.Hm  <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                           goric = goric.Hm)
      df.Hm$model <- as.character(df.Hm$model)
      # complement
      goric.Hc <- -2*(llc - PTc)
      df.c  <- data.frame(model = "complement", loglik = llc, penalty = PTc, 
                          goric = goric.Hc)
      df <- rbind(df.Hm, df.c)
    } else if (type == "gorica") {
      # model
      gorica.Hm <- -2*(llm - PTm)
      df.Hm  <- data.frame(model = objectnames, loglik = llm, penalty = PTm, 
                           gorica = gorica.Hm)
      df.Hm$model <- as.character(df.Hm$model)
      # complement
      gorica.Hc <- -2*(llc - PTc)
      df.c  <- data.frame(model = "complement", loglik = llc, penalty = PTc, 
                          gorica = gorica.Hc)
      df <- rbind(df.Hm, df.c)
    }
  } 

  ans$objectList  <- objectList
  ans$objectNames <- objectnames
  # compute goric weights and relative weights
  delta <- df$goric - min(df$goric)
  goric.weights <- exp(-delta / 2) / sum(exp(-delta / 2))
  df$goric.weights <- goric.weights
  ans$df <- df
  ans$rel.gw <- combn(goric.weights, 2, FUN = function(x) x[1]/x[2], simplify = TRUE)  
  
  
  if (comparison == "complement" && is.null(bound)) {
    ans$betasc <- betasc
  }
  ans$comparison <- comparison
  
  class(ans) <- "goric"
  
  ans
}



print.goric <- function(x, digits = max(3, getOption("digits") - 4), 
                        brief = TRUE, ...) {
  
  dig <- paste0("%6.", digits, "f")
  x2 <- lapply(x$df[-1], sprintf, fmt = dig)
  x2$model <- x$df$model
  df <- as.data.frame(x2)
  df <- df[, c(5,1,2,3,4)]
  
  objectnames <- as.character(df$model)
  goric.weights <- as.numeric(as.character(df$goric.weights))
  rel.gw <- x$rel.gw
  comparison <- x$comparison
  rw <- as.data.frame(diag(length(objectnames)))
  rownames(rw) <- objectnames
  colnames(rw) <- objectnames
  rw[upper.tri(rw, diag = FALSE)] <- rel.gw
  rw[lower.tri(rw, diag = FALSE)] <- 1/rel.gw
  
  Amat <- lapply(x$objectList, FUN = function(x) { x$constraints } )
  meq  <- lapply(x$objectList, FUN = function(x) { x$neq } )
  bvec <- lapply(x$objectList, FUN = function(x) { x$rhs } )
  iact <- lapply(x$objectList, FUN = function(x) { x$iact } )
  
  cat("\nRestriktor: generalized order-restriced information criterion:\n")
  if (!brief) {
    vnames <- names(coef(x$model.org))
    fn <- function(Amat, bvec, meq, iact, vnames) {
      colnames(Amat) <- vnames
      out.rest <- cbind(round(Amat, 4), c(rep("   ==", meq), rep("   >=", nrow(Amat) - 
                                                                   meq)), bvec, " ")
      rownames(out.rest) <- paste(1:nrow(out.rest), ":", sep = "")
      colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs", "active")
      idx <- ncol(out.rest)
      out.rest[, idx] <- "no"
      out.rest[iact, idx] <- "yes"
      if (nrow(Amat) == meq) {
        out.rest[1:nrow(Amat), idx] <- "yes"
      }  
      out.rest <- as.data.frame(out.rest)
      
      out.rest
    }
  
    conMat <- list()
    for (i in 1:length(x$objectList)) {
      conMat[[i]] <- fn(Amat = Amat[[i]], bvec = bvec[[i]], meq = meq[[i]], 
                        iact = iact[[i]], vnames = vnames)  
    }
    names(conMat) <- x$objectNames
    cat("\nConstraint matrices:")
    print(conMat, quote = FALSE, scientific = FALSE)
  }
  
  cat("\n\nGoric-weights:\n")
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  
  cat("\n\nRelative goric-weights:\n")
  print(format(rw, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  cat("---\n")
  cat("Note that in case of equal log-likhood values, the relative 
weights only reflect the difference in penalty term values.\n")
  
  if (comparison == "complement") {
    cat("\nThe order-restricted hypothesis", sQuote(objectnames[1]), "is", 
        sprintf("%1.3f", rel.gw), "times more likely \nto be the best model than its complement.\n")
  } 
  
  cat("\n\nRestricted coefficients:\n")
  coefs <- lapply(x$objectList, FUN = function(x) { coef(x) } )
  max.length <- max(sapply(coefs, length))
  coefs <- lapply(coefs, function(v) { c(v, rep("", max.length-length(v)))})
  
  coefs <- as.data.frame(do.call("rbind", coefs))
  rownames(coefs) <- x$objectNames
  print(coefs, digits = digits, scientific = FALSE, print.gap = 2L,
        quote = FALSE)
  
}
