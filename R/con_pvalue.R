# mixture of F distributions.
con_pvalue_Fbar <- function(wt, Ts.org, df.residual, type = "A",
                            Amat, bvec, meq = 0L, meq.alt = 0L) {
  #check
  if ((qr(Amat)$rank < nrow(Amat))) {
    stop("restriktions matrix must have full row-rank")
  }
  if (type == "global") {
    # compute df
    df.bar <- ((ncol(Amat) - 1) - nrow(Amat)):((ncol(Amat) - 1) - meq)    
    # p value based on the chi-square distribution
    pvalue <- 1-pfbar(Ts.org, df1 = df.bar, df2 = df.residual, wt = rev(wt))
  } else if(type == "A") {
    # compute df
    df.bar <- 0:(nrow(Amat) - meq)
    # p value based on F-distribution or chi-square distribution
    pvalue <- 1-pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
                       wt = rev(wt))
  } else if (type == "B") {
    # compute df
    df.bar <- (meq - meq.alt):(nrow(Amat) - meq.alt)#meq:nrow(Amat)
    # p value based on F-distribution or chi-square distribution
    pvalue <- 1-pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
                       wt = wt)
  } else  {
    stop("hypothesis test type ", sQuote(type), " unknown.")
  }

  out <- pvalue
    attr(out, "wt") <- wt
    attr(out, "df.bar") <- df.bar
  
  out    
}


# mixture of chi-square distributions
con_pvalue_Chibar <- function(wt, Ts.org, type = "A",
                              Amat, bvec, meq = 0L, meq.alt = 0L) {
  #check
  if ((qr(Amat)$rank < nrow(Amat))) {
    stop("Restriktor ERROR: restriktions matrix must have full row-rank")
  }
  
  if (type == "global") {
    # compute df
    df.bar <- ((ncol(Amat) - 1) - nrow(Amat)):((ncol(Amat) - 1) - meq)    
    # p value based on the chi-square distribution
    pvalue <- 1-pchibar(Ts.org, df1 = df.bar, wt = rev(wt))
  }  else if (type == "A") {
    # compute df
    df.bar <- 0:(nrow(Amat) - meq)
    # p value based on th chi-square distribution
    pvalue <- 1-pchibar(Ts.org, df1 = df.bar, wt = rev(wt))
  } else if (type == "B") {
    # compute df
    df.bar <- (meq - meq.alt):(nrow(Amat) - meq.alt)#meq:nrow(Amat)
    # p value based on th chi-square distribution
    pvalue <- 1-pchibar(Ts.org, df1 = df.bar, wt = wt)
  } else  {
    stop("hypothesis test type ", sQuote(type), " unknown.")
  }
  
  out <- pvalue
    attr(out, "wt") <- wt
    attr(out, "df.bar") <- df.bar
  
  out
}


################################################################################
con_pvalue_boot_parametric <- function(model, Ts.org = NULL, 
                                       meq.alt = meq.alt,
                                       test = "F", R = 9999, 
                                       p.distr = "N", df = 7, warn = -1L,
                                       parallel = "no", ncpus = 1L, cl = NULL,
                                       seed = NULL, control = NULL,
                                       verbose = FALSE, ...) {

  p.distr <- tolower(p.distr)
  old_options <- options(); options(warn = warn)
  
  model.org <- model$model.org
  X <- model.matrix(model)[,,drop=FALSE]
  n <- dim(X)[1]

  # constraints   
  Amat <- model$Amat
  bvec <- model$bvec
  meq  <- model$meq
  
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

    if (p.distr == "n") {
      ystar <- rnorm(n = n, 0, 1)
    } else if (p.distr == "t") {
      ystar <- rt(n = n, df = df)
    } else if (p.distr == "chi") {
      ystar <- rchisq(n = n, df = df)
    }
  
    xcol <- which(rowSums(attr(model.org$terms, "factors")) > 0)
    terms <- attr(model.org$terms , "term.labels")
    DATA <- data.frame(ystar, model.org$model[,xcol])
    colnames(DATA) <- c(as.character("ystar"), terms)  
    form <- formula(model.org)
    form[[2]] <- as.name("ystar")
    boot_model <- update(model.org, formula = form, data = DATA)
    
    CALL <- list(model = boot_model, constraints = Amat, rhs = bvec, neq = meq, control = control, se = "none")
    boot_conLM <- do.call("restriktor", CALL)
    boot_conTest <- try(conTest(boot_conLM, type = type, test = test, meq.alt = meq.alt, control = control))
    if (inherits(boot_conTest, "try-error")) {
      if (verbose) cat("FAILED: creating test statistic\n")
      options(old_options)
      return(NULL)
    }
    OUT <- boot_conTest$Ts
  
    if (verbose) {
      cat("iteration =", b, "...Ts =", OUT, "\n")
    }
      OUT
   }

   options(old_options)
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
    pvalue <- sum(Ts.boot >= Ts.org) / Rboot.tot
      attr(pvalue, "B") <- Rboot.tot
    
    OUT <- pvalue

    OUT
}


###################################################################################
con_pvalue_boot_model_based <- function(model, Ts.org = NULL, 
                                        meq.alt = meq.alt,
                                        test = "F", R = 9999, warn = -1L,
                                        parallel = "no", ncpus = 1L,
                                        cl = NULL, seed = NULL, control = NULL,
                                        verbose = FALSE, ...) {

  old_options <- options(); options(warn = warn)
  
  model.org <- model$model.org
  X <- model.matrix(model)[,,drop=FALSE]
  
  # constraints 
  Amat <- model$Amat
  bvec <- model$bvec
  meq  <- model$meq

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

  if (type == "global") {
    intercept <- model.org$assign[1] == 0L
    g <- length(model$b.constr)
    if (intercept) {
      Amatg <- cbind(rep(0, (g - 1)), diag(rep(1, g - 1))) 
      bvecg <- rep(0, g - 1) 
    } else {
      stop("Restriktor ERROR: test not applicable for models without intercept.")      
    }
    call.my <- list(constraints = Amatg, rhs = bvecg, neq = nrow(Amatg), control = control, se = "none")
    call.lm <- list(model = model.org)
    CALL <- c(call.lm, call.my)
    if (any(duplicated(CALL))) {
      stop("duplicated elements in CALL.list")
    }
    fit <- do.call("restriktor", CALL)
  } else if (type == "A") {
    call.my <- list(constraints = Amat, rhs = bvec, neq = nrow(Amat), control = control, se = "none")
    call.lm <- list(model = model.org)
    CALL <- c(call.lm, call.my)
    if (any(duplicated(CALL))) {
      stop("duplicated elements in CALL.list")
    }
    fit <- do.call("restriktor", CALL)
  } else if (type == "B") {
      call.my <- list(constraints = Amat, rhs = bvec, neq = meq, control = control, se = "none")
      call.lm <- list(model = model.org)
      CALL <- c(call.lm, call.my)
      if (any(duplicated(CALL))) {
        stop("duplicated elements in CALL.list")
      }
      fit <- do.call("restriktor", CALL)
    }

    # compute residuals under H0
    r <- residuals(fit)
    yhat <- fitted(fit)

    Ts.boot <- vector("numeric", R)
    fn <- function(b) {
      if (!is.null(seed))
        set.seed(seed + b)
      if (!exists(".Random.seed", envir = .GlobalEnv))
        runif(1)
      RNGstate <- .Random.seed

      idx <- sample(dim(X)[1], replace=TRUE)
      ystar <- as.matrix(c(yhat + r[idx]))
      xcol <- which(rowSums(attr(model.org$terms, "factors")) > 0)
      terms <- attr(model.org$terms , "term.labels")
      DATA <- data.frame(ystar, model.org$model[,xcol])
      colnames(DATA) <- c(as.character("ystar"), terms)
      DATA <- as.data.frame(DATA)
      form <- formula(model.org)
      form[[2]] <- as.name("ystar")
      
      boot_model <- update(model.org, formula = form, data = DATA)
      CALL <- list(boot_model, constraints = Amat, rhs = bvec, neq = meq, control = control, se = "none")
      boot_conLM <- do.call("restriktor", CALL)  
      boot_conTest <- try(conTest(boot_conLM, type = type, test = test, meq.alt = meq.alt, control = control))
      if (inherits(boot_conTest, "try-error")) {
        if (verbose) cat("FAILED: creating test statistic\n")
        options(old_options)
        return(NULL)
      }
      OUT <- boot_conTest$Ts
    
      if (verbose) {
        cat("iteration =", b, "...Ts =", OUT, "\n")
      }
        OUT
    }
    
    options(old_options)
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
    pvalue <- sum(Ts.boot >= Ts.org) / Rboot.tot
      attr(pvalue, "B") <- Rboot.tot
    OUT <- pvalue
    
    OUT
  }



#REF: Silvapulle and Sen (2005, p. 79). Constrained Statistical Inference: Order, 
# Inequality, and Shape Constraints. Hoboken, {NJ}: Wiley
mix.boot <- function(object, 
                     R = 99999, parallel = c("no", "multicore", "snow"),
                     ncpus = 1L, cl = NULL, seed = 1234, verbose = FALSE, ...) {

  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore")
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow")
      have_snow <- TRUE
    if (!have_mc && !have_snow)
      ncpus <- 1L
  }

  # constraints 
  Amat <- object$Amat
  bvec <- rep(0L, nrow(Amat))
  meq  <- object$meq
  
  s2unc <- object$s2.unc
  X <- model.matrix(object)[,,drop = FALSE]
  invW <- kronecker(solve(s2unc), t(X) %*% X)
  W <- solve(invW)
  Dmat <- 2*invW
  
  iact <- vector("numeric", ncol(Amat))
  fn <- function(b) {
    if (verbose)
      cat("iteration =", b, "\n")
    if (!is.null(seed))
      set.seed(seed + b)
    if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1)
    RNGstate <- .Random.seed
    
    Z <- rmvnorm(n = 1, mean = rep(0, ncol(W)), sigma = W)
    dvec <- 2*(Z %*% invW)
    QP <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = t(Amat),
                   bvec = bvec, meq = meq)
    
    if (QP$iact[1] == 0L) {
      return(0L) 
    } else { 
      return(length(QP$iact))
    }    
  }

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
        res <- parallel::parLapply(cl, seq_len(RR), fn)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, seq_len(RR), fn)
    }
  }
  else lapply(seq_len(RR), fn)
  error.idx <- integer(0)
  for (b in seq_len(R)) {
    if (!is.null(res[[b]])) {
      iact[b] <- res[[b]]
    }
    else {
      error.idx <- c(error.idx, b)
    }
  }
  # compute the number of positive components of W.
  dimL <- ncol(W) - iact
  wt <- sapply(1:(ncol(W) + 1), function(x) sum(x == (dimL + 1))) / R
  
  OUT <- wt

  OUT
}


