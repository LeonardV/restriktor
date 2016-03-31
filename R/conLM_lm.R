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
  # acknowledgement: code taken from ic.infer package 
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
  # original R^2
  R2.org <- so$r.squared
  # sample size
  n <- dim(X)[1]
  # unconstrained estimates
  b.unconstr <- coef(model)
  # number of parameters
  p <- length(b.unconstr)
  # variable names
  col.names <- names(b.unconstr)
  row.names <- rownames(b.unconstr)

  if(ncol(Amat) != length(b.unconstr)) {
    stop("length coefficients and ncol(Amat) must be identical")
  }

  # compute log-likelihood
  ll.out <- con_loglik_lm(X = X, y = Y, b = b.unconstr, detU = 1)
  ll <- ll.out$loglik
  s2 <- so$sigma^2
  # check if the constraints are in line with the data!
  if (all(Amat %*% c(b.unconstr) - bvec >= 0 * bvec) & meq == 0) {
    b.constr <- b.unconstr

    OUT <- list(CON = NULL,
                partable = NULL,
                constraints = constraints,
                b.unconstr = b.unconstr,
                b.constr = b.unconstr,
                residuals = resid(model),
                fitted.values = fitted(model),
                weights = weights(model),
                R2.org = R2.org,
                R2.reduced = R2.org,
                df.residual = model$df.residual,
                s2.unc = so$sigma^2, s2.unc.ml = ll.out$Sigma,
                s2 = s2, s2.ml = ll.out$Sigma, 
                loglik = ll, Sigma = vcov(model),
                Amat = Amat, bvec = bvec, meq = meq, iact = NULL, bootout = NULL)  
  } else {
    # compute constrained estimates for lm() and mlm() 
    out.qp <- con_solver(b.unconstr, X = X, y = Y, Amat = Amat,
                         bvec = bvec, meq = meq, tol = tol,
                         maxit = ifelse(is.null(control$maxit), 1e04, 
                                        control$maxit))
    b.constr <- matrix(out.qp$solution, ncol = ncol(Y))
    b.constr[abs(b.constr) < tol] <- 0L

    conll.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
    conll <- conll.out$loglik
    fitted <- model.matrix(model) %*% b.constr
    residuals <- Y - fitted
    
    p2 <- 0L                                                                    # is this correct?                             
    if (meq > 0L) { p2 <- length((1:meq)[bvec == 0]) }
    s2 <- sum(residuals^2) / (n-(p-p2))
    
    # lm
    if (ncol(Y) == 1L) {
      b.constr <- c(b.constr)
      names(b.constr) <- col.names

      # compute constrained R^2, acknowledgement: code taken from ic.infer package (Ulrike Groemping)
      R2.reduced <- 1 - sum(residuals^2)/sum((Y - mean(Y))^2)
      if (is.null(weights(model)) & !attr(model$terms, "intercept"))
        R2.reduced <- 1 - sum(residuals^2)/sum(Y^2)
      if (attr(model$terms, "intercept") & !is.null(weights(model)))
        R2.reduced <- 1 - sum(weights(model) * residuals^2)/sum(weights(model) *
                                                           (Y - weighted.mean(Y, w = weights(model)))^2)
      if (!(attr(model$terms, "intercept") | is.null(weights(model))))
        R2.reduced <- 1 - sum(weights(model) * residuals^2)/sum(weights(model) * Y^2)

    } else { #mlm
      #fitted <- model.matrix(model) %*% b.constr
      #residuals <- Y - fitted
      rownames(b.constr) <- row.names
      R2.org <- R2.reduced <- NULL
      se <- "no"
    }

    OUT <- list(CON = NULL,
                partable = NULL,
                constraints = constraints,
                b.constr = b.constr, b.unconstr = b.unconstr,
                residuals = residuals,
                fitted.values = fitted,
                weights = weights(model),
                R2.org = R2.org,
                R2.reduced = R2.reduced,
                df.residual = model$df.residual,
                s2.unc = so$sigma^2, s2.unc.ml = ll.out$Sigma,
                s2 = s2, s2.ml = conll.out$Sigma, 
                loglik = conll, Sigma = vcov(model),
                Amat = Amat, bvec = bvec, meq = meq, iact = out.qp$iact,
                bootout = NULL)
  }

  OUT$model.org <- model
  OUT$CON <- if (is.character(constraints)) { CON }
  OUT$partable <- if (is.character(constraints)) { partable }
  OUT$se <- se
  if (se != "no") {
    if (!(se %in% c("boot.model.based","boot.standard"))) {
      information <- 1/s2 * crossprod(X) 
      information.inv <- con_augmented_information(information = information,
                                                   X = X, b.unconstr = b.unconstr, 
                                                   b.constr = b.constr,
                                                   Amat = Amat, 
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
  if (ncol(Y) == 1L) {
    class(OUT) <- c("conLM","lm")
  } else if (ncol(Y) > 1L) {
    class(OUT) <- c("conMLM","mlm")
  }
    
    return(OUT)
}
