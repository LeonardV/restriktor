conLM.lm <- function(model, constraints, se = "default", B = 999,
                     rhs = NULL, neq = 0L, bootWt = FALSE, R = 99999,
                     parallel = "no", ncpus = 1L, cl = NULL, seed = NULL, 
                     control = NULL, verbose = FALSE, debug = FALSE, ...) {
  
  # check class
  if (!(class(model)[1] %in% c("lm", "mlm"))) {
    stop("model must be of class lm.")
  }
  
  cl <- match.call()
  # rename for internal use
  bvec <- rhs; meq <- neq
  # construct constraint matrix/vector.
  restr_OUT <- con_constraints(model, constraints = constraints, 
                             bvec = bvec, meq = meq, debug = debug)  
  # a list with useful information about the restriktions.}
  CON <- restr_OUT$CON
  # a parameter table with information about the observed variables in the model 
  # and the imposed restriktions.}
  parTable <- restr_OUT$parTable
  # constraints matrix
  Amat <- restr_OUT$Amat
  # rhs 
  bvec <- restr_OUT$bvec
  # neq
  meq <- restr_OUT$meq
  
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
  
  # standard error methods
  if (se == "default") {
    se <- "standard"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (!(se %in% c("none","standard","const","boot.model.based","boot.standard",
                  "HC","HC0","HC1","HC2","HC3","HC4","HC4m","HC5"))) {
    stop("standard error method ", sQuote(se), " unknown.")
  }
  if (se == "boot.model.based" & any(Amat[,1] == 1)) { 
    stop("no restriktions on intercept possible for 'se = boot.model.based' bootstrap method.")
  }
  
  # model summary
  so <- summary(model)
  # response variable
  y <- as.matrix(model$model[, attr(model$terms, "response")])
  # model matrix
  X <- model.matrix(model)[,,drop = FALSE]
  # residual variance
  s2.unc <- so$sigma^2
  # weigths
  weights <- weights(model)
  # sample size
  n <- dim(X)[1]
  # unconstrained estimates
  b.unrestr <- coef(model)
  # number of parameters
  p <- length(b.unrestr)
  
  if(ncol(Amat) != length(b.unrestr)) {
    stop("length coefficients and ncol(constraints) must be identical")
  }

  # compute (weighted) log-likelihood
  ll.unc <- con_loglik_lm(X = X, y = y, b = b.unrestr, w = weights)
  LL.unc <- ll.unc$loglik
  
  # ML unconstrained MSE
  s2ml.unc <- (s2.unc * model$df.residual) / n
  invW <- kronecker(solve(s2ml.unc), t(X) %*% X)
  W <- solve(invW)
  
  # compute mixing weights
  if (bootWt) { # compute mixing weights based on simulation
    wt <- mix.boot(VCOV = W,
                   Amat = Amat, 
                   meq = meq, 
                   R = R,
                   parallel = parallel,
                   ncpus = ncpus,
                   cl = cl,
                   seed = seed,
                   verbose = verbose)
  } else if (!bootWt & (meq < nrow(Amat))) { # compute mixing weights based on mvnorm
    wt <- rev(con_wt(Amat %*% W %*% t(Amat), meq = meq))
  } else if (!bootWt & (meq == nrow(Amat))) { # only equality constraints
    wt <- rep(0L, ncol(W) + 1)
    wt.idx <- ncol(W) - meq + 1
    wt[wt.idx] <- 1
  }
  attr(wt, "bootWt") <- bootWt
  
  # check if the constraints are in line with the data
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) & meq == 0) {
    b.restr <- b.unrestr
    s2.restr <- s2.unc
    
    OUT <- list(CON = CON,
                parTable = parTable,
                wt = wt,
                b.unrestr = b.unrestr,
                b.restr = b.unrestr,
                residuals = model$residuals, # unweighted residuals
                fitted = model$fitted,
                weights = model$weights,
                df.residual = model$df.residual,
                R2.org = so$r.squared, R2.reduced = so$r.squared,
                s2.unc = s2.unc, s2.restr = s2.unc, 
                loglik = LL.unc, Sigma = vcov(model),
                constraints = Amat, rhs = bvec, neq = meq, iact = NULL, 
                bootout = NULL, call = cl)  
  } else {
    # compute constrained estimates for lm() and mlm() 
    out.QP <- con_solver(b.unrestr, X = X, y = y, w = weights, Amat = Amat,
                         bvec = bvec, meq = meq, 
                         absval = ifelse(is.null(control$absval), 
                                         sqrt(.Machine$double.eps), 
                                         control$absval),
                         maxit = ifelse(is.null(control$maxit), 1e04, 
                                        control$maxit))
    b.restr <- matrix(out.QP$solution, ncol = ncol(y))
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    
    ll.restr <- con_loglik_lm(X = X, y = y, b = b.restr, w = weights)
    LL.restr <- ll.restr$loglik
    # lm
    if (ncol(y) == 1L) {
      b.restr <- as.vector(b.restr)
        names(b.restr) <- names(b.unrestr)
      fitted <- X %*% b.restr
      residuals <- y - fitted
      
      # compute R^2
      if (is.null(weights)) {
        mss <- if (attr(model$terms, "intercept")) {
          sum((fitted - mean(fitted))^2)
        } else {
          sum(fitted^2)
        }
        rss <- sum(residuals^2)
      } else {
        mss <- if (attr(model$terms, "intercept")) {
          m <- sum(weights * fitted / sum(weights))
          sum(weights * (fitted - m)^2)
        } else {
          sum(weights * fitted^2)
        }
        rss <- sum(weights * residuals^2)
      }
      R2.reduced <- mss / (mss + rss)
      
      # compute residual degreees of freedom, corrected for equality constraints. 
      df.residual <- n - (p - qr(Amat[0:meq,])$rank)
      # compute weighted residuals
      if (is.null(weights)) {
        s2.restr <- sum(residuals^2) / df.residual  
      } else {
        s2.restr <- sum(weights * residuals^2) / df.residual
      }
    } else { #mlm <FIXME>
      residuals <- y - fitted
        rownames(b.restr) <- rownames(b.unrestr)
        colnames(b.restr) <- colnames(b.unrestr)
      fitted <- X %*% b.restr
      df <- model$df.residual
      df.residual <- df + qr(Amat[0:meq,])$rank
      se <- "none"
    }

    OUT <- list(CON = CON,
                parTable = parTable,
                wt = wt,
                b.restr = b.restr, 
                b.unrestr = b.unrestr, 
                residuals = residuals, 
                fitted = fitted, 
                weights = weights,
                df.residual = df.residual, 
                R2.org = so$r.squared, R2.reduced = R2.reduced,
                s2.unc = s2.unc,  
                s2.restr = s2.restr,  
                loglik = LL.restr, 
                Sigma = vcov(model), 
                constraints = Amat, rhs = bvec, neq = meq, iact = out.QP$iact,
                bootout = NULL, call = cl)
  }
  
  # original model
  OUT$model.org <- model
  # type standard error
  OUT$se <- se
  # compute standard errors based on the augmented inverted information matrix or
  # based on the standard bootstrap or model.based bootstrap
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      OUT$information <- 1/s2.restr * crossprod(X)
      information <- con_augmented_information(information = OUT$information,
                                               X = X, 
                                               b.unrestr = b.unrestr, 
                                               b.restr = b.restr,
                                               Amat = Amat, 
                                               bvec = bvec, meq = meq) 
      attr(OUT$information, "inverted.information") <- information$inverted.information
      attr(OUT$information, "augmented.information") <- information$augmented.information
    } else if (se == "boot.model.based") {
      OUT$bootout <- con_boot_lm(model, B = B, fixed = TRUE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none",
                                 parallel = parallel, ncpus = ncpus, cl = cl)
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(model, B = B, fixed = FALSE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none",
                                 parallel = parallel, ncpus = ncpus, cl = cl)
    }
  }
  
  if (ncol(y) == 1L) {
    class(OUT) <- c("conLM", "lm")
  } else if (ncol(y) > 1L) {
    class(OUT) <- c("conMLM", "mlm")
  }
    
    OUT
}
