# mixture of F distributions.
con_pvalue_Fbar <- function(wt.bar, Ts.org, df.residual, type = "A",
                            Amat, bvec, meq = 0L, meq.alt = 0L) {
  
  # wt.method <- attr(wt.bar, "method")
  # if (wt.method == "boot") {
  #   idx.start <- (ncol(Amat) - nrow(Amat)) + 1
  #   idx.end   <- (ncol(Amat) - meq) + 1
  # }
  
  if (type == "global") {
    # compute df
    bvecG <- attr(bvec, "bvec.global")
    df.bar <- (length(bvecG) - nrow(Amat)):(length(bvecG) - meq)
    
    ## for testing purposes
    # r <- qr(attr(Amat, "Amat_global"))$rank
    # q <- qr(Amat)$rank
    # i <- 0:q
    # 1 - pfbar(Ts.org, df1 = r-q+i, df2 = df.residual, wt.bar = wt.bar)
    
    ## p value based on the f-distribution
    pvalue <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual, wt.bar = wt.bar)
    
    # 
    # if (wt.method == "boot") {
    #   pvalue <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
    #                       wt.bar = wt.bar[idx.start:idx.end])
    # } else if (wt.method == "pmvnorm") {
    #   pvalue <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
    #                       wt.bar = wt.bar)
    # } else {
    #   stop("Restriktor ERROR: mix.weights method ", sQuote(wt.method), " unknown.")
    # }
  } else if(type == "A") {
    # compute df
    df.bar <- 0:(nrow(Amat) - meq)
    ## p value based on F-distribution or chi-square distribution
    pvalue <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
                        wt.bar = wt.bar)
    
    # if (wt.method == "boot") {
    #   pvalue <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
    #                       wt.bar = wt.bar[idx.start:idx.end])
    # } else if (wt.method == "pmvnorm") {
    #   pvalue <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
    #                       wt.bar = wt.bar)
    # } else {
    #   stop("Restriktor ERROR: mix.weights method ", sQuote(wt.method), " unknown.")
    # }
  } else if (type == "B") {
    # compute df
    df.bar <- (meq - meq.alt):(nrow(Amat) - meq.alt)#meq:nrow(Amat)
    ## p value based on F-distribution or chi-square distribution
    pvalue <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
                        wt.bar = rev(wt.bar))
    
    # if (wt.method == "boot") {
    #   pvalue <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
    #                       wt.bar = rev(wt.bar[idx.start:idx.end]))
    # } else if (wt.method == "pmvnorm") {
    #   pvalue <- 1 - pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
    #                       wt.bar = rev(wt.bar))
    # } else {
    #   stop("Restriktor ERROR: mix.weights method ", sQuote(wt.method), " unknown.")
    # }
  } else  {
    stop("hypothesis test type ", sQuote(type), " unknown.")
  }

  out <- pvalue
    attr(out, "wt.bar") <- wt.bar
    attr(out, "df.bar") <- df.bar
    attr(out, "df.residual") <- df.residual
  
  out    
}


# mixture of chi-square distributions
# con_pvalue_Chibar <- function(wt.bar, Ts.org, type = "A",
#                               Amat, bvec, meq = 0L, meq.alt = 0L) {
#   #check
#   #if ((qr(Amat)$rank < nrow(Amat))) {
#   #  stop("Restriktor ERROR: restriktions matrix must have full row-rank")
#   #}
#   
#   if (type == "global") {
#     # compute df
#     bvecG <- attr(bvec, "bvec.global")
#     #    df.bar <- ((ncol(Amat) - 1) - nrow(Amat)):((ncol(Amat) - 1) - meq)    
#     df.bar <- (length(bvecG) - nrow(Amat)):(length(bvecG) - meq)
#     # p value based on the chi-square distribution
#     pvalue <- 1 - pchibar(Ts.org, df1 = df.bar, wt.bar = rev(wt.bar))
#   }  else if (type == "A") {
#     # compute df
#     df.bar <- 0:(nrow(Amat) - meq)
#     # p value based on th chi-square distribution
#     pvalue <- 1 - pchibar(Ts.org, df1 = df.bar, wt.bar = rev(wt.bar))
#   } else if (type == "B") {
#     # compute df
#     df.bar <- (meq - meq.alt):(nrow(Amat) - meq.alt)#meq:nrow(Amat)
#     # p value based on th chi-square distribution
#     pvalue <- 1 - pchibar(Ts.org, df1 = df.bar, wt.bar = wt.bar)
#   } else  {
#     stop("hypothesis test type ", sQuote(type), " unknown.")
#   }
#   
#   out <- pvalue
#     attr(out, "wt.bar") <- wt.bar
#     attr(out, "df.bar") <- df.bar
#     attr(out, "df.residual") <- df.residual
#   
#   out
# }


############################ parametric bootstrap ##############################
con_pvalue_boot_parametric <- function(model, Ts.org, 
                                       type     = "A",
                                       test     = "F", 
                                       neq.alt  = 0L, 
                                       R        = 9999, 
                                       p.distr  = rnorm, 
                                       parallel = "no", 
                                       ncpus    = 1L, 
                                       cl       = NULL,
                                       seed     = NULL, 
                                       warn     = -1L,
                                       control  = NULL,
                                       verbose  = FALSE, ...) {

  old.options <- options(); options(warn = warn)
  
  model.org <- model$model.org
  fam <- family(model.org)
  
  X <- model.matrix(model)[,,drop=FALSE]
  n <- dim(X)[1]

  # constraints   
  Amat <- model$constraints
  bvec <- model$rhs
  meq  <- model$neq
  
  ## shift y by q to relocate the vertex to its origin.
  if (!all(bvec == 0)) {
    if (!(fam$family %in% c("gaussian"))) {
      stop("Restriktor ERROR: for the ", sQuote(fam$family), " family, relocated cones are not supported (yet).")
    }
    #q <- t(R) %*% solve(R%*%t(R)) %*% bvec
    shift.q <- MASS::ginv(Amat) %*% bvec
    shift.q[abs(shift.q) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), control$tol)] <- 0L
    
    # <FIXME>
    if (fam$family == "gaussian") {
      ystar.shift <- (X %*% -shift.q)
    } else if (fam$family == "binomial") {
      stop("Restriktor ERROR: for the ", sQuote(fam$family), " family, relocated cones are not supported (yet).")
    } else if (fam$family == "poission") {
      stop("Restriktor ERROR: for the ", sQuote(fam$family), " family, relocated cones are not supported (yet).")
    }
    # </FIXME>
  } else {
    ystar.shift <- 0L
  }
  
  #parallel housekeeping
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") {
      have_snow <- TRUE
    }  
    if (!have_mc && !have_snow) {
      ncpus <- 1L
    }  
  }

  Ts.boot <- vector("numeric", R)
   fn <- function(b) {
    if (!is.null(seed)) {
	    set.seed(seed + b) 
    }  
    if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1)
      RNGstate <- .Random.seed

    ystar <- p.distr(n)
    # shift ystar
    ystar.shifted <- ystar - ystar.shift
    xcol <- which(rowSums(attr(model.org$terms, "factors")) > 0)
    terms <- attr(model.org$terms, "term.labels")[attr(model.org$terms, "order") == 1]
    DATA <- data.frame(ystar.shifted, model.org$model[ ,xcol])
      colnames(DATA) <- c(as.character("ystar"), terms)  
    form <- formula(model.org)
    form[[2]] <- as.name("ystar")
    boot.model <- update(model.org, formula = form, data = DATA)
    
    CALL <- list(object = boot.model, constraints = Amat, 
                 rhs = bvec, neq = meq, se = "none", 
                 mix.weights = "none", control = control)
    
    restr.boot <- do.call("restriktor", CALL)
    
    boot.conTest <- try(conTest(restr.boot, 
                                type    = type, 
                                test    = test,
                                boot    = "no",
                                neq.alt = neq.alt, 
                                control = control))
    if (inherits(boot.conTest, "try-error")) {
      if (verbose) cat("FAILED: creating test statistic\n")
      options(old.options)
      return(NULL)
    }
    OUT <- boot.conTest$Ts
  
    if (verbose) {
      cat("bootstrap draw =", b, "...Ts =", OUT, "\n")
    }
      OUT
   }

   options(old.options)
   RR <- sum(R)
     res <- if (ncpus > 1L && (have_mc || have_snow)) {
       if (have_mc) {
         parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
       }
       else if (have_snow) {
         if (is.null(cl)) {
           cl <- parallel::makePSOCKcluster(rep("localhost",
                                                ncpus))
           if (RNGkind()[1L] == "L'Ecuyer-CMRG")
             parallel::clusterSetRNGStream(cl)
           res <- parallel::parLapply(cl, seq_len(RR),
                                      fn)
           parallel::stopCluster(cl)
           res
         } else parallel::parLapply(cl, seq_len(RR), fn)
       }
     } else lapply(seq_len(RR), fn)
     error.idx <- integer(0)
     for (b in seq_len(RR)) {
       if (!is.null(res[[b]])) {
         Ts.boot[b] <- res[[b]]
       } else {
         error.idx <- c(error.idx, b)
       }
     }
     na.boot.idx <- which(is.na(Ts.boot), arr.ind = TRUE)
     inf.boot.idx <- which(Ts.boot == Inf, arr.ind = TRUE)
     idx <- c(na.boot.idx, inf.boot.idx, error.idx)
     idx.unique <- unique(idx)
     Rboot.tot <- (R - length(idx.unique))
     if (length(idx.unique) > 0) {
       Ts.boot <- Ts.boot[-idx.unique]
     }
     if (length(idx.unique) > 0L) {
       warning("restriktor WARNING: only ", (R - length(idx.unique)), 
               " bootstrap draws were successful")
     }
    # > or >= ??? 
    pvalue <- sum(Ts.boot >= as.numeric(Ts.org)) / Rboot.tot
      attr(pvalue, "boot.type") <- "parametric"
      attr(pvalue, "R")         <- Rboot.tot
      attr(pvalue, "Ts.boot")   <- Ts.boot
    
    OUT <- pvalue

    OUT
}


########################### model based bootstrap ##############################
con_pvalue_boot_model_based <- function(model, Ts.org, 
                                        type     = "A",
                                        test     = "F", 
                                        neq.alt  = 0L,
                                        R        = 9999, 
                                        parallel = "no", 
                                        ncpus    = 1L,
                                        cl       = NULL, 
                                        seed     = NULL, 
                                        warn     = -1L,
                                        control  = NULL,
                                        verbose  = FALSE, ...) {

  old.options <- options(); options(warn = warn)
  
  model.org <- model$model.org
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  X <- model.matrix(model.org)[,,drop=FALSE]
  n <- dim(X)[1]
  
  # constraints 
  Amat <- model$constraints
  bvec <- model$rhs
  meq  <- model$neq
  
  # parallel housekeeping
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") {
      have_snow <- TRUE
    }  
    if (!have_mc && !have_snow) {
      ncpus <- 1L
    }  
  }
  
  if (!is.null(attr(type, "type.org"))) {
    type <- "global"
  }
  
  # estimate null model under different hypothesis tests
  if (type == "A") {
    CALL <- list(object = model.org, constraints = Amat, 
                 rhs = bvec, neq = nrow(Amat), 
                 control = control, se = "none", mix.weights = "none")
    
    fit <- do.call("restriktor", CALL)
  } else if (type == "B") {
    CALL <- list(object = model.org, constraints = Amat, 
                 rhs = bvec, neq = meq, 
                 control = control, se = "none", mix.weights = "none")
                    
    fit <- do.call("restriktor", CALL)
  } 

  # compute residuals under the null-hypothesis
  if (type != "global") {
    r    <- residuals(fit)
    yhat <- fitted(fit)
  } else { 
    # if type == global, we can skip the restriktor() function
    w <- weights(model.org)
    W <- diag(w)
    if (!is.null(w)) { 
      yhat <- 1/sum(w) * sum(W %*% y) 
    } else {
      yhat <- mean(y)
    }
    yhat <- cbind(rep(yhat, n))
    r <- y - as.numeric(yhat)
  }
  
  Ts.boot <- vector("numeric", R)
  fn <- function(b) {
    if (!is.null(seed))
      set.seed(seed + b)
    if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1) 
    RNGstate <- .Random.seed

    idx   <- sample(dim(X)[1], replace = TRUE)
    ystar <- as.matrix(c(yhat + r[idx]))
    xcol  <- which(rowSums(attr(model.org$terms, "factors")) > 0)
    terms <- attr(model.org$terms, "term.labels")[attr(model.org$terms, "order") == 1]
    DATA  <- data.frame(ystar, model.org$model[,xcol])
      colnames(DATA) <- c(as.character("ystar"), terms)
    form <- formula(model.org)
    form[[2]] <- as.name("ystar")
    
    boot.model <- update(model.org, formula = form, data = DATA)
    CALL <- list(object = boot.model, constraints = Amat, rhs = bvec, 
                 neq = meq, control = control, se = "none",
                 mix.weights = "none")
    
    restr.boot <- do.call("restriktor", CALL)  
    
    boot.conTest <- try(conTest(restr.boot, 
                                type    = type, 
                                test    = test, 
                                boot    = "no",
                                neq.alt = neq.alt, 
                                control = control))
    if (inherits(boot.conTest, "try-error")) {
      if (verbose) cat("FAILED: creating test statistic\n")
      options(old.options)
      return(NULL)
    }
    OUT <- boot.conTest$Ts
  
    if (verbose) {
      cat("iteration =", b, "...Ts =", OUT, "\n")
    }
      OUT
  }
    
  options(old.options)
  RR <- sum(R)
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
    } else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost",
                                             ncpus))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, seq_len(RR),
                                   fn)
        parallel::stopCluster(cl)
        res
      } else parallel::parLapply(cl, seq_len(RR), fn)
    }
  } else { lapply(seq_len(RR), fn) }
  error.idx <- integer(0)
  for (b in seq_len(RR)) {
    if (!is.null(res[[b]])) {
      Ts.boot[b] <- res[[b]]
    }
    else {
      error.idx <- c(error.idx, b)
    }
  }
  na.boot.idx <- which(is.na(Ts.boot), arr.ind = TRUE)
  inf.boot.idx <- which(Ts.boot == Inf, arr.ind = TRUE)
  idx <- c(na.boot.idx, inf.boot.idx, error.idx)
  idx.unique <- unique(idx)
  Rboot.tot <- (R - length(idx.unique))
  if (length(idx.unique) > 0) {
    Ts.boot <- Ts.boot[-idx.unique]
  }
  if (length(idx.unique) > 0L) {
    warning("restriktor WARNING: only ", (R - length(idx.unique)), 
            " bootstrap draws were successful")
  }
  # > or >= ???
  pvalue <- sum(Ts.boot >= as.numeric(Ts.org)) / Rboot.tot
    attr(pvalue, "boot.type") <- "model.based"
    attr(pvalue, "R")         <- Rboot.tot
    attr(pvalue, "Ts.boot")   <- Ts.boot
    
  
  OUT <- pvalue
  
  OUT
}
