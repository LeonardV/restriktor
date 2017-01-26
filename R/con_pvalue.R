# mixture of F distributions.
con_pvalue_Fbar <- function(wt, Ts_org, df.residual, type = "A",
                            Amat, bvec, meq = 0L, meq_alt = 0L) {
  if (type == "global") {
    # compute df
    bvecG <- attr(bvec, "bvec_global")
    df_bar <- (length(bvecG) - nrow(Amat)):(length(bvecG) - meq)
    
    # for testing purposes
    # r <- qr(attr(Amat, "Amat_global"))$rank
    # q <- qr(Amat)$rank
    # i <- 0:q
    # 1 - pfbar(Ts_org, df1 = r-q+i, df2 = df.residual, wt = rev(wt))
    
    # p value based on the chi-square distribution
    pvalue <- 1 - pfbar(Ts_org, df1 = df_bar, df2 = df.residual, 
                        wt = rev(wt))
  } else if(type == "A") {
    # compute df
    df_bar <- 0:(nrow(Amat) - meq)
    # p value based on F-distribution or chi-square distribution
    pvalue <- 1 - pfbar(Ts_org, df1 = df_bar, df2 = df.residual,
                        wt = rev(wt))
  } else if (type == "B") {
    # compute df
    df_bar <- (meq - meq_alt):(nrow(Amat) - meq_alt)#meq:nrow(Amat)
    # p value based on F-distribution or chi-square distribution
    pvalue <- 1 - pfbar(Ts_org, df1 = df_bar, df2 = df.residual,
                        wt = wt)
  } else  {
    stop("hypothesis test type ", sQuote(type), " unknown.")
  }

  out <- pvalue
    attr(out, "wt") <- wt
    attr(out, "df_bar") <- df_bar
    attr(out, "df.residual") <- df.residual
  
  out    
}


# mixture of chi-square distributions
con_pvalue_Chibar <- function(wt, Ts_org, type = "A",
                              Amat, bvec, meq = 0L, meq_alt = 0L) {
  #check
  #if ((qr(Amat)$rank < nrow(Amat))) {
  #  stop("Restriktor ERROR: restriktions matrix must have full row-rank")
  #}
  
  if (type == "global") {
    # compute df
    bvecG <- attr(bvec, "bvec_global")
    #    df_bar <- ((ncol(Amat) - 1) - nrow(Amat)):((ncol(Amat) - 1) - meq)    
    df_bar <- (length(bvecG) - nrow(Amat)):(length(bvecG) - meq)
    # p value based on the chi-square distribution
    pvalue <- 1 - pchibar(Ts_org, df1 = df_bar, wt = rev(wt))
  }  else if (type == "A") {
    # compute df
    df_bar <- 0:(nrow(Amat) - meq)
    # p value based on th chi-square distribution
    pvalue <- 1 - pchibar(Ts_org, df1 = df_bar, wt = rev(wt))
  } else if (type == "B") {
    # compute df
    df_bar <- (meq - meq_alt):(nrow(Amat) - meq_alt)#meq:nrow(Amat)
    # p value based on th chi-square distribution
    pvalue <- 1 - pchibar(Ts_org, df1 = df_bar, wt = wt)
  } else  {
    stop("hypothesis test type ", sQuote(type), " unknown.")
  }
  
  out <- pvalue
    attr(out, "wt") <- wt
    attr(out, "df_bar") <- df_bar
    attr(out, "df.residual") <- df.residual
  
  out
}


############################ parametric bootstrap ##############################
con_pvalue_boot_parametric <- function(model, 
                                       Ts_org   = NULL, 
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

  #p.distr <- tolower(p.distr)
  #stopifnot(p.distr %in% c("n","t","chi","binom","pois","gamma"))
  old_options <- options(); options(warn = warn)
  
  model_org <- model$model_org
  X <- model.matrix(model)[,,drop=FALSE]
  n <- dim(X)[1]

  # constraints   
  Amat <- model$constraints
  bvec <- model$rhs
  meq  <- model$neq
  
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

  Ts_boot <- vector("numeric", R)
   fn <- function(b) {
    if (!is.null(seed)) {
	    set.seed(seed + b) 
    }  
    if (!exists(".Random.seed", envir = .GlobalEnv))
      runif(1)
      RNGstate <- .Random.seed

    ystar <- p.distr(n)
    xcol <- which(rowSums(attr(model_org$terms, "factors")) > 0)
    terms <- attr(model_org$terms, "term.labels")[attr(model_org$terms, "order") == 1]
    DATA <- data.frame(ystar, model_org$model[ ,xcol])
      colnames(DATA) <- c(as.character("ystar"), terms)  
    form <- formula(model_org)
    form[[2]] <- as.name("ystar")
    boot_model <- update(model_org, formula = form, data = DATA)
    
    CALL <- list(object = boot_model, constraints = Amat, 
                 rhs = bvec, neq = meq, se = "none", 
                 Wt = "none", control = control)
    
    boot_conLM <- do.call("restriktor", CALL)
    
    boot_conTest <- try(conTest(boot_conLM, 
                                type    = type, 
                                test    = test,
                                boot    = "no",
                                neq.alt = neq.alt, 
                                control = control))
    if (inherits(boot_conTest, "try-error")) {
      if (verbose) cat("FAILED: creating test statistic\n")
      options(old_options)
      return(NULL)
    }
    OUT <- boot_conTest[[1]]$Ts
  
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
     error_idx <- integer(0)
     for (b in seq_len(RR)) {
       if (!is.null(res[[b]])) {
         Ts_boot[b] <- res[[b]]
       } else {
         error_idx <- c(error_idx, b)
       }
     }
     na_boot_idx <- which(is.na(Ts_boot), arr.ind = TRUE)
     inf_boot_idx <- which(Ts_boot == Inf, arr.ind = TRUE)
     idx <- c(na_boot_idx, inf_boot_idx, error_idx)
     idx_unique <- unique(idx)
     Rboot_tot <- (R - length(idx_unique))
     if (length(idx_unique) > 0) {
       Ts_boot <- Ts_boot[-idx_unique]
     }
     if (length(idx_unique) > 0L) {
       warning("restriktor WARNING: only ", (R - length(idx_unique)), 
               " bootstrap draws were successful")
     }
    # > or >= ??? 
    pvalue <- sum(Ts_boot >= as.numeric(Ts_org)) / Rboot_tot
      attr(pvalue, "boot_type") <- "parametric"
      attr(pvalue, "R")         <- Rboot_tot
      attr(pvalue, "Ts_boot")   <- Ts_boot
    
    OUT <- pvalue

    OUT
}


########################### model based bootstrap ##############################
con_pvalue_boot_model_based <- function(model, 
                                        Ts_org   = NULL, 
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

  old_options <- options(); options(warn = warn)
  
  model_org <- model$model_org
  y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  X <- model.matrix(model_org)[,,drop=FALSE]
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
  
  if (!is.null(attr(type, "type_org"))) {
    type <- "global"
  }
  
  # estimate null model under different hypothesis tests
  if (type == "A") {
    CALL <- list(object = model_org, constraints = Amat, 
                 rhs = bvec, neq = nrow(Amat), 
                 control = control, se = "none", Wt = "none")
    
    fit <- do.call("restriktor", CALL)
  } else if (type == "B") {
    CALL <- list(object = model_org, constraints = Amat, 
                 rhs = bvec, neq = meq, 
                 control = control, se = "none", Wt = "none")
                    
    fit <- do.call("restriktor", CALL)
  } 

    # compute residuals under the null-hypothesis
    if (type != "global") {
      r    <- residuals(fit)
      yhat <- fitted(fit)
    } else { # if type == global, we can skip the restriktor() function
      w <- weights(model_org)
      W <- diag(w)
      if (!is.null(w)) { 
        yhat <- 1/sum(w) * sum(W %*% y) 
      } 
      else {
        yhat <- mean(y)
      }
      yhat <- cbind(rep(yhat, n))
      r    <- y - as.numeric(yhat)
    }
  
    Ts_boot <- vector("numeric", R)
    fn <- function(b) {
      if (!is.null(seed))
        set.seed(seed + b)
      if (!exists(".Random.seed", envir = .GlobalEnv))
        runif(1) 
      RNGstate <- .Random.seed

      idx   <- sample(dim(X)[1], replace = TRUE)
      ystar <- as.matrix(c(yhat + r[idx]))
      xcol  <- which(rowSums(attr(model_org$terms, "factors")) > 0)
      terms <- attr(model_org$terms , "term.labels")[attr(model_org$terms, "order") == 1]
      DATA  <- data.frame(ystar, model_org$model[,xcol])
        colnames(DATA) <- c(as.character("ystar"), terms)
      form <- formula(model_org)
      form[[2]] <- as.name("ystar")
      
      boot_model <- update(model_org, formula = form, data = DATA)
      CALL <- list(object = boot_model, constraints = Amat, rhs = bvec, 
                   neq = meq, control = control, se = "none",
                   Wt = "none")
      
      boot_conLM <- do.call("restriktor", CALL)  
      boot_conTest <- try(conTest(boot_conLM, 
                                  type    = type, 
                                  test    = test, 
                                  boot    = "no",
                                  neq.alt = neq.alt, 
                                  control = control))
      if (inherits(boot_conTest, "try-error")) {
        if (verbose) cat("FAILED: creating test statistic\n")
        options(old_options)
        return(NULL)
      }
      OUT <- boot_conTest[[1]]$Ts
    
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
    error_idx <- integer(0)
    for (b in seq_len(RR)) {
      if (!is.null(res[[b]])) {
        Ts_boot[b] <- res[[b]]
      }
      else {
        error_idx <- c(error_idx, b)
      }
    }
    na_boot_idx <- which(is.na(Ts_boot), arr.ind = TRUE)
    inf_boot_idx <- which(Ts_boot == Inf, arr.ind = TRUE)
    idx <- c(na_boot_idx, inf_boot_idx, error_idx)
    idx_unique <- unique(idx)
    Rboot_tot <- (R - length(idx_unique))
    if (length(idx_unique) > 0) {
      Ts_boot <- Ts_boot[-idx_unique]
    }
    if (length(idx_unique) > 0L) {
      warning("restriktor WARNING: only ", (R - length(idx_unique)), 
              " bootstrap draws were successful")
    }
    # > or >= ???
    pvalue <- sum(Ts_boot >= as.numeric(Ts_org)) / Rboot_tot
      attr(pvalue, "boot_type") <- "model.based"
      attr(pvalue, "R")         <- Rboot_tot
      attr(pvalue, "Ts_boot")   <- Ts_boot
      
    
    OUT <- pvalue
    
    OUT
  }
