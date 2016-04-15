conLM.lm <- function(model, constraints, se = "default",
                     bvec = NULL, meq = 0L, control = NULL,
                     tol = sqrt(.Machine$double.eps), debug = FALSE, ...) {
  
  Amat <- Amatw; bvec <- bvecw; meq <- meqw
  # check class
  if (!(class(model)[1] %in% c("lm", "mlm"))) {
    stop("ERROR: model must be of class lm.")
  }
  if (se == "default" | se == "standard") {
    se <- "const"
  } else if (se == "boot.residual") {
    se <- "boot.model.based"
  }
  if (missing(constraints) && is.null(bvec)) { 
    bvec <- rep(0L, nrow(Amat)) 
  }
  
  # acknowledgement: this check is taken from the ic.infer package 
  if (nrow(Amat) - meq - 2 > 2) {
    if (!is.numeric(try(matrix(0, floor((nrow(Amat) - meq -
                                           2)/2), choose(nrow(Amat) - meq, floor((nrow(Amat) - meq -
                                                                                  2)/2))), silent = TRUE)))
      stop(paste("test does not work, too many inequality constraints, \n",
                 "interim matrix with ", floor((nrow(Amat) - meq)/2) *
                   choose(nrow(Amat) - meq, floor((nrow(Amat) - meq)/2)),
                 " elements cannot be created", sep = ""))
  }

  # model summary
  so <- summary(model)
  # response variable
  Y <- as.matrix(model$model[, attr(model$terms, "response")])
  # model matrix
  X <- model.matrix(model)[,,drop = FALSE]
  # residual variance
  s2.unc <- so$sigma^2
  # weigths
  w <- weights(model)
  # original R^2
  R2.org <- so$r.squared
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
  ll.unc <- con_loglik_lm(X = X, y = Y, b = b.unconstr, w = w)
  LL.unc <- ll.unc$loglik
  s2.unc.ml <- ll.unc$Sigma / n
  
  # check if the constraints are in line with the data!
  if (all(Amat %*% c(b.unconstr) - bvec >= 0 * bvec) & meq == 0) {
    b.constr <- b.unconstr
    s2.con <- s2.unc
    
    OUT <- list(CON = NULL,
                parTable = parTable,
                constraints = constraints,
                b.unconstr = b.unconstr,
                b.constr = b.unconstr,
                residuals = resid(model),
                fitted.values = fitted(model),
                weights = weights(model),
                R2.org = R2.org,
                R2.reduced = R2.org,
                df.residual = model$df.residual,
                s2.unc = s2.unc, s2.unc.ml = s2.unc.ml,
                s2.con = s2.unc, s2.con.ml = s2.unc.ml, 
                loglik = LL.unc, Sigma = vcov(model),
                Amat = Amat, bvec = bvec, meq = meq, iact = NULL, bootout = NULL)  
  } else {
    # compute constrained estimates for lm() and mlm() 
    out.QP <- con_solver(b.unconstr, X = X, y = Y, w = w, Amat = Amat,
                         bvec = bvec, meq = meq, tol = tol,
                         maxit = ifelse(is.null(control$maxit), 1e04, 
                                        control$maxit))
    b.constr <- matrix(out.QP$solution, ncol = ncol(Y))
    b.constr[abs(b.constr) < tol] <- 0L
    
    ll.con <- con_loglik_lm(X = X, y = Y, b = b.constr, w = w)
    LL.con <- ll.con$loglik
    s2.con.ml <- ll.con$Sigma / n
    
    # lm
    if (ncol(Y) == 1L) {
      b.constr <- as.vector(b.constr)
        names(b.constr) <- names(b.unconstr)
      fitted <- X %*% b.constr
      residuals <- Y - fitted
      # compute residual df corrected for equality constraints. 
      df.residual <- n - (p - qr(Amat[1:meq,])$rank)
      # compute weighted residuals
      if (is.null(w)) {
        s2.con <- sum(residuals^2) / df.residual  
      } else {
        s2.con <- sum(w*residuals^2) / df.residual
      }
      
      # compute reduced R^2.
      R2.reduced <- 1 - sum(residuals^2) / sum((Y - mean(Y))^2)
      # no weights and no intercept
      if (is.null(w) & !attr(model$terms, "intercept")) {
        R2.reduced <- 1 - sum(residuals^2) / sum(Y^2)
      # weights and intercept
      } else if (attr(model$terms, "intercept") & !is.null(w)) {
        R2.reduced <- 1 - sum(w * residuals^2) / sum(w * (Y - weighted.mean(Y, w = w))^2)
      # weights and no intercept 
      } else if (!(attr(model$terms, "intercept") | is.null(w))) {
        R2.reduced <- 1 - sum(w * residuals^2) / sum(w * Y^2)
      }

    } else { #mlm <FIXME>
      residuals <- Y - fitted
        rownames(b.constr) <- rownames(b.unconstr)
        colnames(b.constr) <- colnames(b.unconstr)
      fitted <- X %*% b.constr
      df <- model$df.residual
      df.residual <- df + qr(Amat[1:meq,])$rank
      R2.org <- R2.reduced <- NULL
      se <- "no"
    }

    OUT <- list(CON = NULL,
                parTable = parTable,
                constraints = constraints,
                b.constr = b.constr, # constrained
                b.unconstr = b.unconstr, # unconstrained
                residuals = residuals, # constrained 
                fitted.values = fitted, # constrained 
                weights = w,
                R2.org = R2.org, # unconstrained
                R2.reduced = R2.reduced, # constrained
                df.residual = model$df.residual, # unconstrained
                s2.unc = s2.unc, s2.unc.ml = s2.unc.ml, # unconstrained
                s2.con = s2.con, s2.con.ml = s2.con.ml, # constrained
                loglik = LL.con, # constrained
                Sigma = vcov(model), # unconstrained
                Amat = Amat, bvec = bvec, meq = meq, iact = out.QP$iact,
                bootout = NULL)
  }
  
  # original model
  OUT$model.org <- model
  # CON only exists if constraints input is not a matrix
  OUT$CON <- if (is.character(constraints)) { CON }
  # type standard error
  OUT$se <- se
  # compute standard errors based on the augmented inverted information matrix or
  # based on the standard bootstrap or model.based bootstrap
  if (se != "no") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      information <- 1/s2.con * crossprod(X) 
      OUT$information <- information
      inverted.information <- con_augmented_information(information = information,
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
                                 rhs = bvec, neq = meq, se = "no")
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(model, B = ifelse(is.null(control$B),
                                                   999, control$B),
                                 fixed = FALSE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "no")
    }
  }
  
  if (ncol(Y) == 1L) {
    class(OUT) <- c("conLM","lm")
  } else if (ncol(Y) > 1L) {
    class(OUT) <- c("conLM","mlm")
  }
    
    return(OUT)
}
