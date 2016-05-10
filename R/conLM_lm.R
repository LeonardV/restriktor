conLM.lm <- function(model, constraints, se = "default", B = 999,
                     bvec = NULL, meq = 0L, control = NULL, ...) {
  
  cl <- match.call()
  Amat <- Amatw; bvec <- bvecw; meq <- meqw
  # check class
  if (!(class(model)[1] %in% c("lm", "mlm"))) {
    stop("model must be of class lm.")
  }
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
  if (missing(constraints) && is.null(bvec)) { 
    bvec <- rep(0L, nrow(Amat)) 
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
  
  # check if the constraints are in line with the data
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) & meq == 0) {
    b.restr <- b.unrestr
    s2.restr <- s2.unc
    
    OUT <- list(CON = NULL,
                parTable = parTable,
                #constraints = constraints,
                b.unrestr = b.unrestr,
                b.restr = b.unrestr,
                residuals = model$residuals, # unweighted residuals
                fitted = model$fitted,
                weights = model$weights,
                df.residual = model$df.residual,
                R2.org = so$r.squared, R2.reduced = so$r.squared,
                s2.unc = s2.unc, s2.restr = s2.unc, 
                loglik = LL.unc, Sigma = vcov(model),
                Amat = Amat, bvec = bvec, meq = meq, iact = NULL, 
                bootout = NULL, call = cl)  
  } else {
    # compute constrained estimates for lm() and mlm() 
    out.QP <- con_solver(b.unrestr, X = X, y = y, w = weights, Amat = Amat,
                         bvec = bvec, meq = meq, 
                         tol = ifelse(is.null(control$tol), 
                                      sqrt(.Machine$double.eps), 
                                      control$tol),
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
      R2.reduced <- mss/(mss + rss)
      
      # compute residual df corrected for equality constraints. 
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

    OUT <- list(CON = NULL,
                parTable = parTable,
                #constraints = constraints,
                b.restr = b.restr, #constrained
                b.unrestr = b.unrestr, #unconstrained
                residuals = residuals, #constrained 
                fitted = fitted, #constrained 
                weights = weights,
                #df.residual = model$df.residual, # unconstrained
                df.residual = df.residual, # constrained
                R2.org = so$r.squared, R2.reduced = R2.reduced,
                s2.unc = s2.unc,  #unconstrained
                s2.restr = s2.restr,  #constrained
                loglik = LL.restr, #constrained
                Sigma = vcov(model), #unconstrained
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
      OUT$information <- 1/s2.restr * crossprod(X)
      information <- con_augmented_information(information = OUT$information,
                                               X = X, 
                                               b.unrestr = b.unrestr, 
                                               b.restr = b.restr,
                                               Amat = Amat, 
                                               bvec = bvec, meq = meq) 
      attr(OUT$information, "inverted.information")  <- information$inverted.information
      attr(OUT$information, "augmented.information") <- information$augmented.information
    } else if (se == "boot.model.based") {
      OUT$bootout <- con_boot_lm(model, B = B, 
                                 fixed = TRUE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    } else if (se == "boot.standard") {
      OUT$bootout <- con_boot_lm(model, B = B,
                                 fixed = FALSE, constraints = Amat,
                                 rhs = bvec, neq = meq, se = "none")
    }
  }
  
  if (ncol(y) == 1L) {
    class(OUT) <- c("conLM", "lm")
  } else if (ncol(y) > 1L) {
    class(OUT) <- c("conMLM", "mlm")
  }
    
    return(OUT)
}
