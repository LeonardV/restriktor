con_pvalue_Fbar <- function(cov, Ts.org, df.residual, type = "type A",
                            Amat, bvec, meq = 0L, meq.alt = 0L) {
  #check
  if ((qr(Amat)$rank < nrow(Amat))) {
    stop("constraint matrix must have full row-rank")
  }
  #compute weights
  wt.bar <- con_wt(Amat%*%cov%*%t(Amat), meq=meq)

  if (type == "global") {
    # compute df
    df.bar <- ((ncol(Amat) - 1) - nrow(Amat)):((ncol(Amat) - 1) - meq)    
    # p value based on the chi-square distribution
    pvalue <- 1-pfbar(Ts.org, df1 = df.bar, df2 = df.residual, wt = rev(wt.bar))
  } else if(type == "A") {
    # compute df
    df.bar <- 0:(nrow(Amat) - meq)
    # p value based on F-distribution or chi-square distribution
    pvalue <- 1-pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
                       wt = rev(wt.bar))
  } else if (type == "B") {
    # compute df
    df.bar <- (meq - meq.alt):(nrow(Amat) - meq.alt)#meq:nrow(Amat)
    # p value based on F-distribution or chi-square distribution
    pvalue <- 1-pfbar(Ts.org, df1 = df.bar, df2 = df.residual,
                       wt = wt.bar)
  } else if (type == "C") {
    # t-distribution
    pvalue <- 1-pt(Ts.org, df.residual)
    names(pvalue) <- "pt.value"
  }

  pvalue
}



con_pvalue_Chibar <- function(cov, Ts.org, df.residual, type = "A",
                              Amat, bvec, meq = 0L, meq.alt = 0L) {
  #check
  if ((qr(Amat)$rank < nrow(Amat))) {
    stop("Restriktor ERROR: constraint matrix must have full row-rank")
  }
  #compute weights
  wt.bar <- con_wt(Amat%*%cov%*%t(Amat), meq=meq)
  
  if (type == "global") {
    # compute df
    df.bar <- ((ncol(Amat) - 1) - nrow(Amat)):((ncol(Amat) - 1) - meq)    
    # p value based on the chi-square distribution
    pvalue <- 1-pchibar(Ts.org, df1 = df.bar, wt = rev(wt.bar))
  }  else if (type == "A") {
    # compute df
    df.bar <- 0:(nrow(Amat) - meq)
    # p value based on th chi-square distribution
    pvalue <- 1-pchibar(Ts.org, df1 = df.bar, wt = rev(wt.bar))
  } else if (type == "B") {
    # compute df
    df.bar <- (meq - meq.alt):(nrow(Amat) - meq.alt)#meq:nrow(Amat)
    # p value based on th chi-square distribution
    pvalue <- 1-pchibar(Ts.org, df1 = df.bar, wt = wt.bar)
  } else if (type == "C") {
    # t-distribution
    pvalue <- 1-pt(Ts.org, df.residual)
    names(pvalue) <- "pt.value"
  }
  
  pvalue
}


################################################################################
con_pvalue_boot_parametric <- function(model, Ts.org = NULL, type = "A",
                                       test = "F", bvec = NULL,
                                       meq = NULL, meq.alt = 0,
                                       R = 9999, p.distr = "N", df = 7,
                                       parallel = "no", ncpus = 1L, cl = NULL,
                                       seed = NULL, control = NULL,
                                       verbose = FALSE, ...) {

  p.distr <- tolower(p.distr)
  if (type == "C") { 
    stop("Restriktor ERROR: type C is based on a t-distribution. Set boot = \"no\" ") 
  }
  
  model.org <- model$model.org
  X <- model.matrix(model.org)[,,drop=FALSE]
  n <- dim(X)[1]
  
#  constraints <- model$constraints
  Amat <- model$Amat
  bvec <- model$bvec
  meq <- model$meq
  
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
    boot_model <- update(model.org, formula = form, data = DATA, maxit = 5000)
    OUT <- NA
    if (boot_model$converged) {
      CALL <- list(model = boot_model, constraints = Amat, rhs = bvec, neq = meq, control = control, se = "no")
      boot_conLM <- do.call("restriktor", CALL)  
      boot_conTest <- conTest(boot_conLM, type = type, test = test, meq.alt = meq.alt, control = control)
      OUT <- boot_conTest$Ts
    }
    if (verbose) {
      cat("iteration =", b, "...Ts =", OUT, "\n")
    }
      OUT
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
     idx <- c(na.boot.idx, inf.boot.idx)
     idx.unique <- unique(idx)
     Rboot.tot <- (R - length(idx.unique))
     if (length(idx.unique) > 0) {
       Ts.boot <- Ts.boot[-idx.unique, ]
     }

    # > or >= ??? 
    pvalue <- sum(Ts.boot >= Ts.org) / Rboot.tot
      attr(pvalue, "B") <- Rboot.tot
    
    OUT <- pvalue

    OUT
}


###################################################################################
con_pvalue_boot_model_based <- function(model, Ts.org = NULL, type = "A",
                                        test = "F", meq.alt = 0,
                                        R = 9999, parallel = "no", ncpus = 1L,
                                        cl = NULL, seed = NULL, control = NULL,
                                        verbose = FALSE, ...) {

  if (type == "C") { 
    stop("type C is based on a t-distribution. Set boot = \"no\" ") 
  }
  model.org <- model$model.org
  X <- model.matrix(model.org)[,,drop=FALSE]
  
  #constraints <- model$constraints
  Amat <- model$Amat
  bvec <- model$bvec
  meq <- model$meq

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
    call.my <- list(constraints = Amatg, rhs = bvecg, neq = nrow(Amatg), control = control, se ="no")
    call.lm <- list(model = model.org)
    CALL <- c(call.lm, call.my)
    if (any(duplicated(CALL))) {
      stop("duplicated elements in CALL.list")
    }
    fit <- do.call("restriktor", CALL)
  } else if (type == "A") {
    call.my <- list(constraints = Amat, rhs = bvec, neq = nrow(Amat), control = control, se ="no")
    call.lm <- list(model = model.org)
    CALL <- c(call.lm, call.my)
    if (any(duplicated(CALL))) {
      stop("duplicated elements in CALL.list")
    }
    fit <- do.call("restriktor", CALL)
  } else if (type == "B") {
      call.my <- list(constraints = Amat, rhs = bvec, neq = meq, control = control, se = "no")
      call.lm <- list(model = model.org)
      CALL <- c(call.lm, call.my)
      if (any(duplicated(CALL))) {
        stop("duplicated elements in CALL.list")
      }
      fit <- do.call("restriktor", CALL)
    }

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
      
      OUT <- NA
      boot_model <- update(model.org, formula = form, data = DATA, maxit = 5000)
      if (boot_model$converged) {
        CALL <- list(boot_model, constraints = Amat, rhs = bvec, neq = meq, control = control, se = "no")
        boot_conLM <- do.call("restriktor", CALL)  
        boot_conTest <- conTest(boot_conLM, type = type, test = test, meq.alt = meq.alt, control = control)$Ts
        OUT <- boot_conTest
      }
      
      if (verbose) {
        cat("iteration =", b, "...Ts =", OUT, "\n")
      }
        OUT
    }

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
    idx <- c(na.boot.idx, inf.boot.idx)
    idx.unique <- unique(idx)
    Rboot.tot <- (R - length(idx.unique))
    if (length(idx.unique) > 0) {
      Ts.boot <- Ts.boot[-idx.unique, ]
    }
    # > or >= ???
    pvalue <- sum(Ts.boot >= Ts.org) / Rboot.tot
      attr(pvalue, "B") <- Rboot.tot
    OUT <- pvalue
    
    OUT
  }


############## NEED FIX - NOT IMPLEMENTED YET ################
mix.boot <- function(object, R = 999, p.distr = c("N", "t", "Chi"), df = 7,
                     parallel = c("no", "multicore", "snow"),
                     ncpus = 1L, cl = NULL, seed = NULL, verbose = FALSE, ...) {

  if(class(object) != "test.lm") {
    stop("object must be of class test.lm().")
  }
  if (missing(parallel)) {
    parallel <- getOption("boot.parallel", "no")
  }
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

  if (missing(p.distr))
    p.distr <- getOption("p.distr", "N")
  p.distr <- match.arg(p.distr)

  Amat <- object$Amat
  meq <- object$meq

  start.idx <- meq + 2
  idx <- start.idx:ncol(Amat)# + (start.idx - 1))
  p <- length(coef(object$object.lm))
  Amatw <- rbind(diag(p)[idx,])
  idx <- sapply(1:nrow(Amatw), function(x) which(Amatw[x,] == 1))
  bvec <- object$bvec
  if(meq > 0L) {
    bvec <- object$bvec[-meq]
  }

  X <- model.matrix(object$object.lm)[,,drop=FALSE]
  n <- dim(X)[1]
  beta.unconstr <- coef(object$object.lm)

  posPar <- matrix(as.numeric(NA), R, nrow(Amat)-meq)
  fn <- function(b) {
    if (verbose)
      cat("iteration =", b, "\n")
    if (!is.null(seed))
      set.seed(seed + b)
    if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1)
    RNGstate <- .Random.seed
    if (p.distr == "N") {
      ystar <- rnorm(n, 0, 1)
    }
    else if (p.distr == "t") {
      ystar <- rt(n, df = df)
    }
    else if (p.distr == "Chi") {
      ystar <- rchisq(n, df = df)
    }

    beta.constr <- constr.solve(beta.unconstr = beta.unconstr, x = X,
                                y = cbind(ystar),
                                Amat = Amatw, bvec = bvec, meq = 0L)$solution
    beta.constr[abs(beta.constr) < sqrt(.Machine$double.eps)] <- 0L
    OUT <- beta.constr[idx]

    OUT
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
      posPar[b, 1:nrow(Amatw)] <- res[[b]]
    }
    else {
      error.idx <- c(error.idx, b)
    }
  }
  posPar <- sapply(1:R, function(x) sum(posPar[x,] > 0L))
  wt.bar <- sapply(nrow(Amatw):0, function(x) sum(posPar == x)/R)
    names(wt.bar) <- nrow(Amatw):0

  OUT <- wt.bar

  OUT
}
