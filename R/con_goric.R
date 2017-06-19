goric <- function(object, ..., complement = FALSE, 
                  digits = max(3, getOption("digits") - 2)) {
  
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
    warning("Restriktor WARNING: if complement = TRUE, only one inequality-constrained hypothesis\n",
            "                      is allowed. Therefore, the complement is set to FALSE.")
  } 
  
  df.c <- NULL
  if (complement) {
    Hm <- conlist[[1]]
    Amat <- Hm$constraints
    bvec <- Hm$rhs
    meq  <- Hm$neq
    wt.bar <- Hm$wt.bar
    b.restr <- Hm$b.restr
    b.unrestr <- Hm$b.unrestr
    
    # compute log-likelihood for complement
    if (!(all(Amat %*% b.unrestr - bvec >= 0 * bvec)) & meq == 0) {
      ll.Hc <- logLik(Hm$model.org)
    } else if (nrow(Amat) > meq && !(all(c(Amat) == 0L))) {
      ll <- list()
      # number of rows
      nr <- 1:nrow(Amat)
      # remove rows corresponding to equality constraints
      if (meq > 0L) {
        nr <- nr[-c(0:meq)]
      }
      # treat each row of Amat as an equality constraint.
      for (l in 1:length(nr)) {
        idx <- c(nr[l], nr[-l])
        Amatx <- Amat[idx,,drop = FALSE]
        Hc.restr <- restriktor(Hm$model.org, constraints = Amatx, 
                               neq = 1, rhs = bvec[idx], 
                               mix.weights = "none", se = "none")
        ll[[l]] <- Hc.restr$loglik
      }
      # take the highest log-likelihood value as a substitute for 
      # the complement.
      ll.Hc <- max(unlist(ll))
    } else if (nrow(Amat) == meq){
      # in case of equality constraints only, the complement is 
      # equal to the unconstrained log-likelihood.
      ll.Hc <- logLik(Hm$model.org)
    } else if (all(c(Amat) == 0L)) {
      # unconstrained setting
      stop("Restriktor ERROR: the complement cannot be computed for the unconstrained hypothesis.")
    }
    
    # compute 'true' free parameters t
    p <- length(b.restr)
    f <- p - nrow(Amat) + meq
    # true free parameters
    t <- p - f
    idx <- length(wt.bar)
    # compute penalty term value PTc
    if (attr(wt.bar, "method") == "boot") {
      PTc <- as.numeric(1 + (p - wt.bar[idx-meq] * t))
    } else if (attr(wt.bar, "method") == "pmvnorm") {
      PTc <- as.numeric(1 + (p - wt.bar[idx] * t))
    } else {
      stop("Restriktor ERROR: no chi-bar-square weights found.")
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
  df <- rbind(df, df.c)
  
  delta <- df$goric - min(df$goric)
  goric.weights <- exp(-delta / 2) / sum(exp(-delta / 2))
  df$goric.weights <- goric.weights
  
  dig <- paste0("%6.", digits, "f")
  df2 <- sapply(df[,-1], sprintf, fmt = dig)
  mnames <- as.character(df[,1])
  df2 <- as.data.frame(cbind(model = mnames, df2))
  print(format(df2, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE)
  
  if (complement) {
    Gw.m <- goric.weights[1]
    Gw.c <- goric.weights[2]
    cat("The order-restricted hypothesis", objectnames[1], "is", 
        sprintf("%1.3f", Gw.m/Gw.c), "times more likely than its complement.")
  }
  
  invisible(df)
}
