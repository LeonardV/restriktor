bootstrapD <- function(h0 = NULL, h1 = NULL, constraints, type = "A", 
                       bootstrap.type = "bollen.stine", R = 1000L,  
                       return.D = FALSE, double.bootstrap = "no", 
                       double.bootstrap.R = 500L, double.bootstrap.alpha = 0.05, 
                       verbose = FALSE, warn = -1L, 
                       parallel = c("no", "multicore", "snow"), ncpus = 1L, cl = NULL, 
                       seed = NULL) {
  
  bootstrap.type <- tolower(bootstrap.type)
  stopifnot(inherits(h1, "lavaan"), 
            bootstrap.type %in% c("bollen.stine", "parametric", 
                                  "yuan", "nonparametric", "ordinary"), 
            double.bootstrap %in% c("no", "FDB", "standard"),
            type %in% c("A", "B"))
  
  if (bootstrap.type == "nonparametric") {
    bootstrap.type <- "ordinary"
  }
  if (h1@Model@conditional.x) {
    stop("lavaan ERROR: this function is not (yet) available if conditional.x = TRUE")
  }
  old_options <- options()
  options(warn = warn)
  D <- rep(as.numeric(NA), R)
  
  if (double.bootstrap == "FDB") {
    D.2 <- numeric(R)
  }
  else if (double.bootstrap == "standard") {
    plugin.pvalues <- numeric(R)
  }
  if (missing(parallel)) 
    parallel <- "no"
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  }
  
  data <- h0@Data
  
  if (bootstrap.type == "bollen.stine" || bootstrap.type == "parametric" || 
      bootstrap.type == "yuan") {
    model.implied <- lav_model_implied(h0@Model)
    Sigma.hat <- model.implied$cov
    Mu.hat    <- model.implied$mean
  }
  if (bootstrap.type == "bollen.stine" || bootstrap.type == "yuan") {
    if (h0@Options$missing != "listwise") 
      stop("lavaan ERROR: bollen.stine/yuan bootstrap not available for missing data")
    dataX <- vector("list", length = data@ngroups)
  }
  else {
    dataX <- data@X
  }
  
  lavdata <- h0@Data
  lavoptions <- h0@Options
  dataeXo <- lavdata@eXo
  dataWT <- lavdata@weights
  dataMp <- lavdata@Mp
  dataRp <- lavdata@Rp
  dataLp <- lavdata@Lp
  
  if (bootstrap.type == "bollen.stine") {
    for (g in 1:h0@Data@ngroups) {
      sigma.sqrt <- lav_matrix_symmetric_sqrt(Sigma.hat[[g]])
      S.inv.sqrt <- lav_matrix_symmetric_sqrt(h0@SampleStats@icov[[g]])
      X <- scale(data@X[[g]], center = TRUE, scale = FALSE)
      X <- X %*% S.inv.sqrt %*% sigma.sqrt
      if (h0@Model@meanstructure) 
        X <- scale(X, center = (-1 * Mu.hat[[g]]), scale = FALSE)
      dataX[[g]] <- X
    }
  }
  if (bootstrap.type == "yuan") {
    g.a <- function(a, Sigmahat, Sigmahat.inv, S, tau.hat, 
                    p) {
      S.a <- a * S + (1 - a) * Sigmahat
      tmp.term <- S.a %*% Sigmahat.inv
      res1 <- (sum(diag(tmp.term)) - log(det(tmp.term)) - 
                 p) - tau.hat
      res <- res1 * res1
      attr(res, "gradient") <- sum(diag((S - Sigmahat) %*% 
                                          (Sigmahat.inv - chol2inv(chol(S.a)))))
      res
    }
    for (g in 1:h0@Data@ngroups) {
      S <- h0@SampleStats@cov[[g]]
      ghat <- h0@test[[1]]$stat.group[[g]]
      df <- h0@test[[1]]$df
      Sigmahat <- Sigma.hat[[g]]
      Sigmahat.inv <- solve(Sigmahat)
      nmv <- nrow(Sigmahat)
      n <- data@nobs[[g]]
      tau.hat <- (ghat - df)/(n - 1)
      if (tau.hat >= 0) {
        a <- optimize(g.a, c(0, 1), Sigmahat, Sigmahat.inv, 
                      S, tau.hat, nmv)$minimum
        S.a <- a * S + (1 - a) * Sigmahat
      }
      else {
        S.a <- Sigmahat
      }
      S.a.sqrt   <- lav_matrix_symmetric_sqrt(S.a)
      S.inv.sqrt <- lav_matrix_symmetric_sqrt(h0@SampleStats@icov[[g]])
      X <- data@X[[g]]
      X <- X %*% S.inv.sqrt %*% S.a.sqrt
      dataX[[g]] <- X
    }
  }
  
 
  ## compute original D value
  # create constraint matrix based on the character string containing the constraints.
  parTable <- as.list(parTable(h1)) 
  conInfo  <- lav_constraints_parse(parTable, constraints = constraints, 
                                            theta = lavaan::coef(h1))
  # equality constraints
  Amat.ceq <- conInfo$ceq.JAC
  # inequality constraints
  Amat.ciq <- conInfo$cin.JAC
  # constraints matrix, equality constraints go first!
  Amat <- rbind(Amat.ceq, Amat.ciq)

  # right-hand side 
  rhs.ceq <- conInfo$ceq.rhs
  rhs.ciq <- conInfo$cin.rhs
  rhs <- c(rhs.ceq, rhs.ciq)
  
  theta.hat <- lavaan::coef(h1)
  I.hat <- lavInspect(h1, "information")
  N <- lavaan::nobs(h1) 
  
  if ("A" %in% type) {
  ## hypothesis test Type A
  # fit equality constraint model.
    out.eq <- solve.QP(Dmat = I.hat, dvec = theta.hat %*% I.hat, Amat = t(Amat),
                       bvec = rhs, meq = nrow(Amat))
    # fit inequality constraint model, some equality constraints may be preserved.
    out.ineq <- solve.QP(Dmat = I.hat, dvec = theta.hat %*% I.hat, Amat = t(Amat),
                         bvec = rhs, meq = nrow(Amat.ceq)) 
    
    D.original <- N * 2 * (out.eq$value - out.ineq$value)
  }
  
  if ("B" %in% type) {
  ## hypothesis test Type B
    # fit inequality constraint model, some equality constraints may be preserved.
    out.ineq <- solve.QP(Dmat = I.hat, dvec = theta.hat %*% I.hat, Amat = t(Amat),
                         bvec = rhs, meq = nrow(Amat.ceq)) 
    
    # fit unconstraint model, some equality constraints may be preserved.         FIXME
    out.un <- solve.QP(Dmat = I.hat, dvec = theta.hat %*% I.hat, 
                       Amat = t(matrix(0, nrow = dim(Amat)[1], ncol = dim(Amat)[2])),
                       bvec = rep(0, length(rhs)), meq = 0)
    
  
    D.original <- N * 2 * (out.ineq$value - out.un$value)
  }
  
  
  fn <- function(b) {
    if (bootstrap.type == "bollen.stine" || bootstrap.type == "ordinary" || 
        bootstrap.type == "yuan") {
      BOOT.idx <- vector("list", length = lavdata@ngroups)
      for (g in 1:lavdata@ngroups) {
        stopifnot(lavdata@nobs[[g]] > 1L)
        boot.idx <- sample(x = lavdata@nobs[[g]], size = lavdata@nobs[[g]], 
                           replace = TRUE)
        BOOT.idx[[g]] <- boot.idx
        dataX[[g]] <- dataX[[g]][boot.idx, , drop = FALSE]
      }
      newData <- lav_data_update(lavdata = lavdata, newX = dataX, 
                                 BOOT.idx = BOOT.idx, lavoptions = lavoptions)
    }
    else {
      for (g in 1:lavdata@ngroups) {
        dataX[[g]] <- MASS::mvrnorm(n = lavdata@nobs[[g]], 
                                    Sigma = Sigma.hat[[g]], mu = Mu.hat[[g]])
      }
      newData <- lav_data_update(lavdata = lavdata, newX = dataX, 
                                 lavoptions = lavoptions)
    }
    if (verbose) 
      cat("  ... bootstrap draw number: ", b, "\n")
    bootSampleStats <- try(lav_samplestats_from_data(lavdata = newData, 
                                                     lavoptions = lavoptions), silent = TRUE)
    if (inherits(bootSampleStats, "try-error")) {
      if (verbose) 
        cat("     FAILED: creating h0@SampleStats statistics\n")
      options(old_options)
      return(NULL)
    }

    # h0
    if (verbose) 
      cat("  ... ... model h0: ")
    h0@Options$verbose <- FALSE
    h0@Options$se <- "none"
    h0@Options$test <- "none"
    h0@Options$baseline <- FALSE
    h0@Options$h1 <- FALSE
    
    
    fit.boot.h0 <- lavaan(slotOptions = h0@Options, slotParTable = h0@ParTable, 
                                  slotSampleStats = bootSampleStats, slotData = data)
    
    if (!fit.boot.h0@optim$converged) {
      if (verbose) 
        cat("     FAILED: no convergence\n")
      options(old_options)
      return(NULL)
    }
    if (verbose) 
      cat("     ok -- niter = ", fit.boot.h0@optim$iterations, 
          " fx = ", fit.boot.h0@optim$fx, "\n")
    
    # h1
    if (verbose) 
      cat("  ... ... model h1: ")
    h1@Options$verbose <- FALSE
    h1@Options$se <- "none"
    h1@Options$test <- "none"
    h1@Options$baseline <- FALSE
    h1@Options$h1 <- FALSE
    
    
    fit.boot.h1 <- lavaan(slotOptions = h1@Options, slotParTable = h1@ParTable, 
                                  slotSampleStats = bootSampleStats, slotData = data)
    
    if (!fit.boot.h1@optim$converged) {
      if (verbose) 
        cat("     FAILED: no convergence\n")
      options(old_options)
      return(NULL)
    }
    if (verbose) 
      cat("     ok -- niter = ", fit.boot.h1@optim$iterations, 
          " fx = ", fit.boot.h1@optim$fx, "\n")
    
    # data are sampled under the null-model.
    theta.hat.boot <- lavaan::coef(fit.boot.h1)
    I.hat.boot <- lavInspect(fit.boot.h1, "information")
    

    if ("A" %in% type) {
    # hypothesis test Type A
      out.eq.boot   <- solve.QP(Dmat = I.hat.boot, dvec = theta.hat.boot %*% I.hat.boot, Amat = t(Amat),
                                bvec = rhs, meq = nrow(Amat))
      out.ineq.boot <- solve.QP(Dmat = I.hat.boot, dvec = theta.hat.boot %*% I.hat.boot, Amat = t(Amat),
                                bvec = rhs, meq = nrow(Amat.ceq))
      
      D.boot <- N * 2 * (out.eq.boot$value - out.ineq.boot$value)
    }
    
    if ("B" %in% type) {
      # hypothesis test Type B
      out.ineq.boot <- solve.QP(Dmat = I.hat.boot, dvec = theta.hat.boot %*% I.hat.boot, Amat = t(Amat),
                                bvec = rhs, meq = nrow(Amat.ceq))
      out.un.boot <- solve.QP(Dmat = I.hat.boot, dvec = theta.hat.boot %*% I.hat.boot, 
                              Amat = t(matrix(0, nrow = dim(Amat)[1], ncol = dim(Amat)[2])),
                              bvec = rep(0, length(rhs)), meq = 0)
      
      D.boot <- N * 2 * (out.ineq.boot$value - out.un.boot$value)
    }
    
    
    if (double.bootstrap == "standard") {
      if (verbose)
        cat("  ... ... calibrating p.value - ")
      plugin.pvalue <- bootstrapD(h0 = fit.boot.h0, h1 = fit.boot.h1,
                                  constraints = constraints,
                                  R = double.bootstrap.R,
                                  bootstrap.type = bootstrap.type, verbose = FALSE,
                                  return.D = FALSE,
                                  double.bootstrap = "no", warn = warn,
                                  parallel = parallel, ncpus = ncpus, cl = cl)
      if (verbose)
        cat(sprintf("%5.3f", plugin.pvalue), "\n")
      attr(D.boot, "plugin.pvalue") <- plugin.pvalue
    }
    else if (double.bootstrap == "FDB") {
      plugin.pvalue <- bootstrapD(h0 = fit.boot.h0, h1 = fit.boot.h1,
                                  constraints = constraints,
                                  R = 1L, bootstrap.type = bootstrap.type, 
                                  verbose = FALSE, warn = warn,
                                  return.D = TRUE, parallel = parallel, ncpus = ncpus,
                                  cl = cl, double.bootstrap = "no")
      D.2 <- attr(plugin.pvalue, "D")
      if (verbose)
        cat("  ... ... D2 = ", D.2, "\n")
      attr(D.boot, "D.2") <- D.2
    }
    D.boot
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
          parallel::clusterSetRNGStream(cl, seed = seed)
        
        res <- parallel::parLapply(cl, seq_len(RR), fn)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, seq_len(RR), fn)
    }
  }
  else lapply(seq_len(RR), fn)
  error.idx <- integer(0)
  for (b in seq_len(RR)) {
    if (!is.null(res[[b]])) {
      D[b] <- res[[b]]
      if (double.bootstrap == "standard") {
        plugin.pvalues[b] <- attr(res[[b]], "plugin.pvalue")
      }
      else if (double.bootstrap == "FDB") {
        D.2[b] <- attr(res[[b]], "D.2")
      }
    }
    else {
      error.idx <- c(error.idx, b)
    }
  }
  if (length(error.idx) > 0L) {
    warning("lavaan WARNING: only ", (R - length(error.idx)), 
            " bootstrap draws were successful")
    D <- D[-error.idx]
    if (length(D) == 0) 
      D <- as.numeric(NA)
    if (double.bootstrap == "standard") {
      plugin.pvalues <- plugin.pvalues[-error.idx]
      attr(D, "error.idx") <- error.idx
    }
    if (double.bootstrap == "FDB") {
      D.2 <- D.2[-error.idx]
      attr(D.2, "error.idx") <- error.idx
    }
  }
  else {
    if (verbose) 
      cat("Number of successful bootstrap draws:", 
          (R - length(error.idx)), "\n")
  }
  pvalue <- sum(D > D.original) / length(D)
  if (return.D) {
    attr(pvalue, "D.original") <- D.original
    attr(pvalue, "D") <- D
  }
  if (double.bootstrap == "FDB") {
    Q <- (1 - pvalue)
    d.q <- quantile(D.2, Q, na.rm = TRUE)
    adj.pvalue <- sum(D > d.q)/length(D)
    attr(pvalue, "D.q") <- d.q
    attr(pvalue, "adj.pvalue") <- adj.pvalue
    if (return.D) {
      attr(pvalue, "D.original") <- D.original
      attr(pvalue, "D") <- D
      attr(pvalue, "D2") <- D.2
    }
  }
  else if (double.bootstrap == "standard") {
    adj.alpha <- quantile(plugin.pvalues, double.bootstrap.alpha, 
                          na.rm = TRUE)
    attr(pvalue, "adj.alpha") <- adj.alpha
    adj.pvalue <- sum(plugin.pvalues < pvalue)/length(plugin.pvalues)
    attr(pvalue, "plugin.pvalues") <- plugin.pvalues
    attr(pvalue, "adj.pvalue") <- adj.pvalue
  }
  options(old_options)
  pvalue
}
