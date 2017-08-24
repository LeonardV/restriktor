goric <- function(object, ..., complement = FALSE, 
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
    #b.restr <- object$b.restr
    # restricted parameters
    b.unrestr <- object$b.unrestr
    # level probabilities
    wt.bar <- object$wt.bar
    
    # constraints matrix
    Amat <- object$constraints
    # number of equalities
    meq  <- object$neq
    # rhs
    bvec <- object$rhs
    # extract equalities
    Amat.ceq <- Amat[0:meq, ,drop = FALSE]
    bvec.ceq <- bvec[0:meq]
    # check if any equality constraint is violated
    check.ceq <- !(all(Amat.ceq %*% c(b.unrestr) - bvec.ceq == 0))
    # extract inequalities
    
    if (nrow(Amat) > meq) {
      # extract inequalities
      Amat.cin <- Amat[(meq + 1):nrow(Amat), ,drop = FALSE]
      bvec.cin <- bvec[(meq + 1)]
      # check if any inequality constraint is violated
      check.cin <- !(all(Amat.cin %*% c(b.unrestr) - bvec.cin >= 0))
    } else {
      check.cin <- FALSE
    }
    
    ## compute log-likelihood for complement
    # check if any constraint is violated
    if (check.cin || check.ceq) { 
      ll.Hc <- con_loglik_lm(object$model.org)
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
        Hc.restr <- restriktor(object$model.org, constraints = Amatx, 
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
      ll.Hc <- con_loglik_lm(object$model.org)
      if (debug) {
        cat("log-likelihood value =", ll.Hc, "\n")
      }
    } else if (all(c(Amat) == 0L)) {
      # unconstrained setting
      stop("Restriktor ERROR: no complement exists for the unconstrained hypothesis.")
    } else {
      stop("Restriktor ERROR: you might have found a bug, please contact me!")
    }
    
    # compute free parameters f in the complement
    p <- length(b.unrestr)
    # nrow q1
    lq1 <- nrow(Amat.cin)
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
    goric.weights <- as.numeric(levels(df$goric.weights))
    Gm <- c(goric.weights[1])
    Gc <- goric.weights[2]
    cat("The order-restricted hypothesis", sQuote(objectnames[1]), "is", 
        sprintf("%1.3f", Gm/Gc), "times more likely than its complement.")
  }
}
