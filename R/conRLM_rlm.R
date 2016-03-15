#compute constraint estimates
conRLM.rlm <- function(model, constraints, debug = FALSE,
                       se = "default", control = NULL,
                       bvec = NULL, meq = 0L,
                       tol = sqrt(.Machine$double.eps), ...) { 
  
  Amat <- Amatw; bvec <- bvecw; meq <- meqw
  
  # check class
  if (!("rlm" %in% class(model))) {
    stop("ERROR: model must be of class rlm.")
  }
  if (se == "default" | se == "standard") {
    se <- "const"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (missing(constraints) && is.null(bvec)) { 
    bvec <- rep(0L, nrow(Amat)) 
  }
  
  # acknowledgement: code taken from ic.infer package (Ulrike Groemping)
  if (nrow(Amat) - meq - 2 > 2) {
    if (!is.numeric(try(matrix(0, floor((nrow(Amat) - meq -
                                         2)/2), choose(nrow(Amat) - meq, floor((nrow(Amat) - meq -
                                                                                2)/2))), silent = TRUE)))
      stop(paste("test does not work, too many inequality constraints, \n",
                 "interim matrix with ", floor((nrow(Amat) - meq)/2) *
                   choose(nrow(Amat) - meq, floor((nrow(Amat) - meq)/2)),
                 " elements cannot be created", sep = ""))
  }
  
  
  if (is.null(model$model)) {
    Y <- model$residuals + model$fitted
  } else {
#    Data <- model$model
#    mfit <- model$model
    Y <- as.matrix(model$model[, attr(model$terms, "response")])
  }  
  X  <- model.matrix(model)[,,drop = FALSE]
#  n <- dim(X)[1]
  b.unconstr <- coef(model)
  so <- MASS:::summary.rlm(model)
  df.residual <- so$df[2]
  # number of parameters
#  p <- length(b.unconstr)
#  df.residual <- nrow(X)-ncol(X)
  
  if (ncol(Amat) != length(b.unconstr)) {
    stop("length coefficients and ncol(Amat) must be identical")
  }
  
  if (all(Amat %*% c(b.unconstr) - bvec >= 0 * bvec) & meq == 0) {  
  #R2 and adjusted R2, acknowledgement: code taken from lmrob() function.
  resid <- resid(model)
  pred <- model$fitted.values
  resp <- pred + resid 
  wgt <- model$w
  if (is.null(model$model)) {
    df.int <- if (any(colnames(X) == "(Intercept)")) 1L else 0L
  } else {
    df.int <- if (attr(model$terms, "intercept")) 1L else 0L
  } 
  resp.mean <- if (df.int == 1L) sum(wgt * resp)/sum(wgt) else 0
  yMy <- sum(wgt * (resp - resp.mean)^2)
  rMr <- sum(wgt * resid^2)
  # bi-square correction
  correc <- 1.207617 
  R2.org <- r2correc <- (yMy - rMr) / (yMy + rMr * (correc - 1))
  #R2.adjusted <- 1 - (1 - r2correc) * ((n - df.int) / df.residual)
  
  
  call.rlm <- as.list(model$call)
  call.rlm <- call.rlm[-1]
  # the default is M-estimation with psi = psi.huber loss function.
  # only psi.bisquare is applicable for now.
  if (is.null(call.rlm$method) && is.null(call.rlm$psi)) {
    stop("Restiktor ERROR: test only applicable with \"psi=psi.bisquare\".")
  } else if (call.rlm$method == "M" && is.null(call.rlm$psi)) {
    stop("Restiktor ERROR: test only applicable with \"psi=psi.bisquare\".")
  } else if (call.rlm$method == "M" && !(call.rlm$psi == "psi.bisquare")) {
    stop("Restiktor ERROR: test only applicable with \"psi=psi.bisquare\".")
  } else if (call.rlm$method == "MM" && !(call.rlm$psi == "psi.bisquare")) {
    stop("Restiktor ERROR: test only applicable with \"psi=psi.bisquare\".")
  }
  
  #fit inequality constrained robust model
  if (is.null(call.rlm[["formula"]])) {
    call.rlm[["data"]] <- NULL
    call.rlm[["x"]] <- NULL
    call.rlm[["y"]] <- NULL
    call.my <- list(x = X, y = Y)      
    CALL <- c(call.rlm, call.my)
    rfit <- do.call("rlm_fit", CALL)
  }
  else {
    call.my <- list()  
    call.rlm[["data"]] <- as.name("Data")
    CALL <- c(call.rlm, call.my)
    rfit <- do.call("rlm.formula", CALL)
  }
  
  b.unconstr <- coef(rfit)
  # compute log-likelihood
  ll.out <- con_loglik_lm(X = X, y = Y, b = b.unconstr, detU = 1)
  ll <- ll.out$loglik
  b.constr <- b.unconstr
  s2 <- so$stddev^2
  
  OUT <- list(CON = NULL,
              partable = NULL,
              constraints = constraints,
              b.unconstr = b.unconstr,
              b.constr = b.unconstr,
              residuals = resid(model),
              resid0 = rfit$resid0,
              fitted.values = fitted(model),
              weights = model$weights,
              R2.org = R2.org,
              R2.reduced = R2.org,
              df.residual = df.residual,
              scale = model$s, 
              s2 = s2,  
              loglik = ll, 
              Sigma = vcov(model),                                      #probably not robust!
              Amat = Amat, bvec = bvec, meq = meq, iact = 0L,
              converged = model$converged,
              bootout = NULL)
  } else {
    call.rlm <- as.list(model$call)
    call.rlm <- call.rlm[-1]
    # the default is M-estimation with psi = psi.huber loss function.
    # only psi.bisquare is applicable for now.
    if (is.null(call.rlm$method) && is.null(call.rlm$psi)) {
      stop("Restiktor ERROR: test only applicable with \"psi=psi.bisquare\".")
    } else if (call.rlm$method == "M" && is.null(call.rlm$psi)) {
      stop("Restiktor ERROR: test only applicable with \"psi=psi.bisquare\".")
    } else if (call.rlm$method == "M" && !(call.rlm$psi == "psi.bisquare")) {
      stop("Restiktor ERROR: test only applicable with \"psi=psi.bisquare\".")
    } else if (call.rlm$method == "MM" && !(call.rlm$psi == "psi.bisquare")) {
      stop("Restiktor ERROR: test only applicable with \"psi=psi.bisquare\".")
    }
    
    #fit inequality constrained robust model
    if (is.null(call.rlm[["formula"]])) {
      call.rlm[["data"]] <- NULL
      call.rlm[["x"]] <- NULL
      call.rlm[["y"]] <- NULL
      call.my <- list(x = X, y = Y, Amat = Amat, meq = meq, bvec = bvec)      
      CALL <- c(call.rlm, call.my)
      rfit <- do.call("conRLM_fit", CALL)
    }
    else {
      call.my <- list(Amat = Amat, meq = meq, bvec = bvec)  
      call.rlm[["data"]] <- as.name("Data")
      CALL <- c(call.rlm, call.my)
      rfit <- do.call("conRLM.formula", CALL)
    }
    
    b.constr <- coef(rfit)
      b.constr[abs(b.constr) < sqrt(.Machine$double.eps)] <- 0L
    ll.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
    ll <- ll.out$loglik
    
    # vector with the indices of the active constraints
    iact <- rfit$iact
    #scale parameter
  #  scale <- model$s
    #standard deviation
    tau.hat <- MASS:::summary.rlm(rfit)$stddev  
    s2 <- tau.hat^2
    #R2 and adjusted R2, code taken from lmrob() function.
    resid <- rfit$residuals
    pred <- rfit$fitted.values
    resp <- pred + resid 
    wgt <- rfit$w
    if (is.null(model$model)) {
      df.int <- if (any(colnames(X) == "(Intercept)")) 1L else 0L
    } else {
      df.int <- if (attr(model$terms, "intercept")) 1L else 0L
    }
    resp.mean <- if (df.int == 1L) sum(wgt * resp)/sum(wgt) else 0
    yMy <- sum(wgt * (resp - resp.mean)^2)
    rMr <- sum(wgt * resid^2)
    # bi-square correction
    correc <- 1.207617 
    R2.reduced <- r2correc <- (yMy - rMr) / (yMy + rMr * (correc - 1))
    #R2.adjusted <- 1 - (1 - r2correc) * ((n - df.int) / df.residual)
  
    OUT <- list(CON = NULL,
                partable = NULL,
                constraints = constraints,
                b.unconstr = b.unconstr,
                b.constr = b.constr,
                residuals = resid,
                resid0 = rfit$resid0,
                wresid = rfit$wresid,
                fitted.values = pred,
                weights = rfit$weights,
                R2.org = R2.org,
                R2.reduced = R2.reduced,
                df.residual = df.residual,
                scale = rfit$s,                                                #unconstrained or constrained scale estimate???
                s2 = s2, #standard deviation tau.hat                                         
                loglik = ll, 
                Sigma = vcov(rfit),                                             #probably not correct???
                Amat = Amat, bvec = bvec, meq = meq, iact = iact,
                converged = rfit$converged,
                bootout = NULL)
  }
  
  OUT$model.org <- model
  OUT$CON <- if (is.character(constraints)) { CON }
  OUT$partable <- if (is.character(constraints)) { partable }
  OUT$se <- se
  
  if (se != "no") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      information <- solve(vcov.MM(X, rfit$resid0, rfit$residuals, rfit$s))     #fix sandwich for rlm!
      information.inv <- con_augmented_information(information = information,
                                                   X = X, b.unconstr = b.unconstr, 
                                                   b.constr = b.constr,
                                                   s2 = s2, constraints = Amat, 
                                                   bvec = bvec, meq = meq)
      OUT$information.inverted <- information.inv
    } else if (se == "boot.model.based") {
      OUT$bootout <- con_boot_lm(model, B = ifelse(is.null(control$B),
                                                   999, control$B), 
                                 fixed = TRUE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "no")
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(model, B = ifelse(is.null(control$B),
                                                   999, control$B),
                                 fixed = FALSE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "no")
    }
  }
  
  class(OUT) <- c("conRLM","conLM")
  
  return(OUT)
}
