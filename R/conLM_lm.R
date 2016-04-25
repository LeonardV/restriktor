conLM.lm <- function(model, constraints, se = "default",
                     bvec = NULL, meq = 0L, control = NULL,
                     tol = sqrt(.Machine$double.eps), debug = FALSE, ...) {
  
  cl <- match.call()
  Amat <- Amatw; bvec <- bvecw; meq <- meqw
  # check class
  if (!(class(model)[1] %in% c("lm", "mlm"))) {
    stop("model must be of class lm.")
  }
  if (se == "default" | se == "standard") {
    se <- "const"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (!(se %in% c("none","const","boot.model.based","boot.standard","HC","HC0",
                  "HC1","HC2","HC3","HC4","HC4m","HC5"))) {
    stop("standard error method ", sQuote(se), " unknown.")
  }
  if (se == "boot.model.based" & any(Amat[,1] == 1)) { 
    stop("no restriktions on intercept possible for model based bootstrap.")
  }
  if (missing(constraints) && is.null(bvec)) { 
    bvec <- rep(0L, nrow(Amat)) 
  }
  # acknowledgement: this check is taken from the ic.infer package 
  if (nrow(Amat) - meq - 2 > 2) {
    if (!is.numeric(try(matrix(0, floor((nrow(Amat) - meq -
                                           2)/2), choose(nrow(Amat) - meq, floor((nrow(Amat) - meq -
                                                                                  2)/2))), silent = TRUE)))
      stop(paste("test does not work, too many inequality restriktions, \n",
                 "interim matrix with ", floor((nrow(Amat) - meq)/2) *
                   choose(nrow(Amat) - meq, floor((nrow(Amat) - meq)/2)),
                 " elements cannot be created", sep = ""))
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
  b.unconstr <- coef(model)
  # number of parameters
  p <- length(b.unconstr)
  
  if(ncol(Amat) != length(b.unconstr)) {
    stop("length coefficients and ncol(Amat) must be identical")
  }

  # compute (weighted) log-likelihood
  ll.unc <- con_loglik_lm(X = X, y = y, b = b.unconstr, w = weights)
  LL.unc <- ll.unc$loglik
  # ml residual variance
  s2ml.unc <- ll.unc$Sigma / n
  
  # check if the constraints are in line with the data
  if (all(Amat %*% c(b.unconstr) - bvec >= 0 * bvec) & meq == 0) {
    b.constr <- b.unconstr
    s2.con <- s2.unc
    
    OUT <- list(CON = NULL,
                parTable = parTable,
                constraints = constraints,
                b.unconstr = b.unconstr,
                b.constr = b.unconstr,
                residuals = model$residuals, # unweighted residuals
                fitted = model$fitted,
                weights = model$weights,
                df.residual = model$df.residual,
                R2.org = so$r.squared, R2.reduced = so$r.squared,
                s2.unc = s2.unc, s2ml.unc = s2ml.unc,
                s2.con = s2.unc, s2ml.con = s2ml.unc, 
                loglik = LL.unc, Sigma = vcov(model),
                Amat = Amat, bvec = bvec, meq = meq, iact = NULL, 
                bootout = NULL, call = cl)  
  } else {
    # compute constrained estimates for lm() and mlm() 
    out.QP <- con_solver(b.unconstr, X = X, y = y, w = weights, Amat = Amat,
                         bvec = bvec, meq = meq, tol = tol,
                         maxit = ifelse(is.null(control$maxit), 1e04, 
                                        control$maxit))
    b.constr <- matrix(out.QP$solution, ncol = ncol(y))
    b.constr[abs(b.constr) < tol] <- 0L
    
    ll.con <- con_loglik_lm(X = X, y = y, b = b.constr, w = weights)
    LL.con <- ll.con$loglik
    s2ml.con <- ll.con$Sigma / n
    # lm
    if (ncol(y) == 1L) {
      b.constr <- as.vector(b.constr)
        names(b.constr) <- names(b.unconstr)
      fitted <- X %*% b.constr
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
      R2.reduced <- mss/(mss + rss)
      
      # compute residual df corrected for equality constraints. 
      df.residual <- n - (p - qr(Amat[0:meq,])$rank)
      # compute weighted residuals
      if (is.null(weights)) {
        s2.con <- sum(residuals^2) / df.residual  
      } else {
        s2.con <- sum(weights * residuals^2) / df.residual
      }
    } else { #mlm <FIXME>
      residuals <- y - fitted
        rownames(b.constr) <- rownames(b.unconstr)
        colnames(b.constr) <- colnames(b.unconstr)
      fitted <- X %*% b.constr
      df <- model$df.residual
      df.residual <- df + qr(Amat[0:meq,])$rank
      se <- "none"
    }

    OUT <- list(CON = NULL,
                parTable = parTable,
                constraints = constraints,
                b.constr = b.constr, # constrained
                b.unconstr = b.unconstr, # unconstrained
                residuals = residuals, # constrained 
                fitted = fitted, # constrained 
                weights = weights,
                df.residual = model$df.residual, # unconstrained
                #df.residual = df.residual, # constrained
                R2.org = so$r.squared, R2.reduced = R2.reduced,
                s2.unc = s2.unc, s2ml.unc = s2ml.unc, # unconstrained
                s2.con = s2.con, s2ml.con = s2ml.con, # constrained
                loglik = LL.con, # constrained
                Sigma = vcov(model), # unconstrained
                Amat = Amat, bvec = bvec, meq = meq, iact = out.QP$iact,
                bootout = NULL, call = cl)
  }
  
  # original model
  OUT$model.org <- model
  # CON only exists if constraints input is not a matrix
  OUT$CON <- if (is.character(constraints)) { CON }
  # type standard error
  OUT$se <- se
  # compute standard errors based on the augmented inverted information matrix or
  # based on the standard bootstrap or model.based bootstrap
  if (se != "none") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      OUT$information <- 1/s2.con * crossprod(X) 
      inverted.information <- con_augmented_information(information = OUT$information,
                                                        X = X, 
                                                        b.unconstr = b.unconstr, 
                                                        b.constr = b.constr,
                                                        Amat = Amat, 
                                                        bvec = bvec, meq = meq) 
      attr(OUT$information, "inverted.information") <- inverted.information        
    } else if (se == "boot.model.based") {
      OUT$bootout <- con_boot_lm(model, B = ifelse(is.null(control$B),
                                                   999, control$B), 
                                 fixed = TRUE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(model, B = ifelse(is.null(control$B),
                                                   999, control$B),
                                 fixed = FALSE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    }
  }
  
  if (ncol(y) == 1L) {
    class(OUT) <- c("conLM","lm")
  } else if (ncol(y) > 1L) {
    class(OUT) <- c("conMLM","mlm")
  }
    
    return(OUT)
}
