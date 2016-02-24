# computes the F- or LRT-statistic
# To do: implement E-bar test statistic.

conTestF.lm <- function(object, type = "A", boot = "none", meq.alt = 0,
                        control = NULL, tol = sqrt(.Machine$double.eps)) {

  # housekeeping
  type <- toupper(type)
  #test <- tolower(test)

  if (!("conLM" %in% class(object))) {
    stop("object must be of class conLM().")
  }
  if(!(type %in% c("A","B","C"))) {
    stop("type must be \"A\", \"B\" or \"C\"")
  }

  if(!(boot %in% c("none", "residual", "model.based", "parametric"))) {
    stop("ERROR: boot method unknown.")
  }

  if (boot == "residual") {
    boot <- "model.based"
  }

  constraints <- object$constraints
  model.org <- object$model.org
  X <- model.matrix(model.org)[,,drop=FALSE]
  Y <- model.org$model[, attr(model.org$terms, "response")]
  #n <- dim(X)[1]
  cov <- object$Sigma
  b.unconstr <- object$b.unconstr
  vnames <- names(b.unconstr)
  b.constr <- object$b.constr
  b.eqconstr <- NULL
  b.constr.alt <- NULL
  #iact <- NULL
  #s2 <- NULL
  Ts <- as.numeric(NA)
    names(Ts) <- "Fbar"
  Amat <- object$Amat
  bvec <- object$bvec
  meq <- object$meq

  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality constraints only.")
  }
  if (!(length(b.unconstr) == nrow(cov))) {
    stop("length b.unconstr and nrow(cov) must be identical.")
  }

  if (type == "A") {
    # optimizer
    b.eqconstr <- con_solver_lm(b.unconstr, X = X, y = Y, Amat = Amat,
                                bvec = bvec, meq = nrow(Amat),
                                tol = ifelse(is.null(control$tol), 1e-09, control$tol),
                                maxit = ifelse(is.null(control$maxit), 1e04, control$maxit))$solution
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
      names(b.eqconstr) <- vnames
    Ts <- c(t(b.constr - b.eqconstr) %*% solve(cov, b.constr - b.eqconstr))
  }
  else if (type == "B") {
    if (meq.alt == 0L) {
      # Fbar test statistic
      Ts <- as.vector(t(b.unconstr - b.constr) %*% solve(cov, b.unconstr - b.constr))
    }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if(meq.alt != 0L && meq.alt <= meq) {
        b.constr.alt <- con_solver_lm(b.unconstr, X = X, y = Y,
                                      Amat = Amat[1:meq.alt,,drop=FALSE],
                                      bvec = bvec[1:meq.alt], meq = meq.alt,
                                      tol = ifelse(is.null(control$tol), 1e-09, control$tol),
                                      maxit = ifelse(is.null(control$maxit), 1e04, control$maxit))$solution

        Ts <- as.vector(t(b.constr - b.constr.alt) %*% solve(cov, b.constr - b.constr.alt))
      }
      else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  } else if (type == "C") { # intersection-union test (Sasabuchi, 1980)
    if (meq == 0L) {
      Ts <- as.vector(min((Amat %*% b.unconstr - bvec) /
                            sqrt(diag(Amat %*%cov%*% t(Amat)))))
      names(Ts) <- "Tbar"
    }
    else {
      stop("test not applicable with equality constraints.")
    }
  }

  #Ts[abs(Ts) < sqrt(.Machine$double.eps)] <- 0L


  if (boot == "none") {
    pvalue <- con_pvalue_default_lm(cov, Ts.org = Ts, object$df.residual, type = type,
                                    Amat, bvec, meq, meq.alt)
  }
  else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric_lm(X = X, Ts.org = Ts, type = type, test = "Fbar",
                                            constraints = constraints, bvec = bvec, meq = meq, meq.alt = meq.alt,
                                            R = ifelse(is.null(control$B), 9999, control$B),
                                            p.distr = ifelse(is.null(control$p.distr), "N", control$p.distr),
                                            df = ifelse(is.null(control$df), 7, control$df),
                                            parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                            ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                            cl = ifelse(is.null(control$cl), NULL, control$cl),
                                            seed = ifelse(is.null(control$seed), 1234, control$seed),
                                            verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  }
  else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based_lm(object, Ts.org = Ts, type = type, test = "Fbar",
                                             meq.alt = meq.alt,
                                             R = ifelse(is.null(control$B), 9999, control$B),
                                             parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                             ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                             cl = ifelse(is.null(control$cl), NULL, control$cl),
                                             seed = ifelse(is.null(control$seed), 1234, control$seed),
                                             verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  }

  OUT <- list(CON = object$CON,
              type = type,
              boot = boot,
              b.eqconstr = NULL,
              b.unconstr = b.unconstr,
              b.constr = b.constr,
              b.constr.alt = b.constr.alt,
              Amat = Amat,
              bvec = bvec,
              meq = meq,
              meq.alt = meq.alt,
              iact = object$iact,
              s2 = object$s2,
              df.residual = object$df.residual,
              cov = cov,
              Ts = Ts,
              pvalue = pvalue,
              model.org = model.org)


  if (type == "A") { 
    OUT$b.eqconstr <- b.eqconstr 
  }

  class(OUT) <- "conTest"

  OUT

}






conTestLRT.lm <- function(object, type = "A", boot = "none", meq.alt = 0,
                           control = NULL, tol = sqrt(.Machine$double.eps)) {

  # housekeeping
  type <- toupper(type)
  #test <- tolower(test)

  if (!("conLM" %in% class(object))) {
    stop("object must be of class conLM().")
  }
  if(!(type %in% c("A","B","C"))) {
    stop("type must be \"A\", \"B\" or \"C\"")
  }

  if(!(boot %in% c("none", "residual", "model.based", "parametric"))) {
    stop("ERROR: boot method unknown.")
  }

  if (boot == "residual") {
    boot <- "model.based"
  }

  constraints <- object$constraints
  model.org <- object$model.org
  X <- model.matrix(model.org)[,,drop=FALSE]
#  n <- dim(X)[1]
  Y <- cbind(model.org$model[, attr(model.org$terms, "response")])
  cov <- object$Sigma
  b.unconstr <- object$b.unconstr
    var.names <- names(b.unconstr)
  b.constr <- object$b.constr
  b.eqconstr <- NULL
  b.constr.alt <- NULL
 # iact <- NULL
  s2 <- NULL
  Ts <- as.numeric(NA)
    names(Ts) <- "LRT"

  Amat <- object$Amat
  bvec <- object$bvec
  meq <- object$meq

  if (meq == nrow(Amat)) {
    stop("test not applicable for object with equality constraints only.")
  }
  if (!(length(b.unconstr) == nrow(cov))) {
    stop("length b.unconstr and nrow(cov) must be identical.")
  }

  if (type == "A") {
    # optimizer
    b.eqconstr <- con_solver_lm(b.unconstr, X = X, y = Y, Amat = Amat,
                                bvec = bvec, meq = nrow(Amat),
                                tol = ifelse(is.null(control$tol), 1e-09, control$tol),
                                maxit = ifelse(is.null(control$maxit), 1e04, control$maxit))$solution
    b.eqconstr[abs(b.eqconstr) < tol] <- 0L
      names(b.eqconstr) <- var.names
    # test statistics
    ll0.out <- con_loglik_lm(X = X, y = Y, b = b.eqconstr, detU = 1)
    ll0 <- ll0.out$loglik
    s2 <- ll0.out$Sigma

    ll1.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
    ll1 <- ll1.out$loglik
    Ts <- -2*(ll0 - ll1)
  }
  else if (type == "B") {
    if (meq.alt == 0L) {
      # test statistic
      ll0.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
      ll0 <- ll0.out$loglik
      s2 <- ll0.out$Sigma

      ll1.out <- con_loglik_lm(X = X, y = Y, b = b.unconstr, detU = 1)
      ll1 <- ll1.out$loglik
      Ts <- -2*(ll0 - ll1)
    }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if(meq.alt != 0L && meq.alt <= meq) {
        b.constr.alt <- con_solver_lm(b.unconstr, X = X, y = Y,
                                      Amat = Amat[1:meq.alt,,drop=FALSE],
                                      bvec=bvec[1:meq.alt], meq = meq.alt,
                                      tol = ifelse(is.null(control$tol), 1e-09, control$tol),
                                      maxit = ifelse(is.null(control$maxit), 1e04, control$maxit))$solution
          names(b.constr.alt) <- var.names

        ll0.out <- con_loglik_lm(X = X, y = Y, b = b.constr, detU = 1)
        ll0 <- ll0.out$loglik
        s2 <- ll0.out$Sigma

        ll1 <- con_loglik_lm(X = X, y = Y, b = b.constr.alt, detU = 1)
        ll1 <- ll1.out$loglik
        Ts <- -2*(ll0 - ll1)
      }
      else {
        stop("meq.alt must not be larger than meq.")
      }
    }
  }
  # intersection-union test (Sasabuchi, 1980)
  else if (type == "C") {
    if (meq == 0L) {
      Ts <- as.vector(min((Amat %*% b.unconstr - bvec) /
                            sqrt(diag(Amat %*%cov%*% t(Amat)))))
      names(Ts) <- "Tbar"
    }
    else {
      stop("test not applicable with equality constraints.")
    }
  }

  if (boot == "none") {
    pvalue <- con_pvalue_default_lm(cov, Ts.org = Ts, object$df.residual, type = type,
                                    Amat, bvec, meq, meq.alt)
  }
  else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric_lm(X = X, Ts.org = Ts, type = type, test = "LRT",
                                            constraints = constraints, meq.alt = meq.alt,
                                            R = ifelse(is.null(control$B), 9999, control$B),
                                            p.distr = ifelse(is.null(control$p.distr), "N", control$p.distr),
                                            df = ifelse(is.null(control$df), 7, control$df),
                                            parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                            ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                            cl = ifelse(is.null(control$cl), NULL, control$cl),
                                            seed = ifelse(is.null(control$seed), 1234, control$seed),
                                            verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  }
  else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based_lm(object, Ts.org = Ts, type = type, test = "LRT",
                                             meq.alt = meq.alt,
                                             R = ifelse(is.null(control$B), 9999, control$B),
                                             parallel = ifelse(is.null(control$parallel), "no", control$parallel),
                                             ncpus = ifelse(is.null(control$ncpus), 1, control$ncpus),
                                             cl = ifelse(is.null(control$cl), NULL, control$cl),
                                             seed = ifelse(is.null(control$seed), 1234, control$seed),
                                             verbose = ifelse(is.null(control$verbose), FALSE, control$verbose))
  }

  OUT <- list(CON = object$CON,
              type = type,
              boot = boot,
              b.eqconstr = NULL,
              b.unconstr = b.unconstr,
              b.constr = b.constr,
              b.constr.alt = b.constr.alt,
              Amat = Amat,
              bvec = bvec,
              meq = meq,
              meq.alt = meq.alt,
              iact = object$iact,
              s2 = s2,
              df.residual = object$df.residual,
              cov = cov,
              Ts = Ts,
              pvalue = pvalue,
              model.org = object$model.org)

  if(type == "A") { OUT$b.eqconstr <- b.eqconstr }

  class(OUT) <- "conTest"

  OUT

}

