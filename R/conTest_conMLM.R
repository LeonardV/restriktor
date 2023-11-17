conTestLRT.conMLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                              p.distr = rnorm, parallel = "no", ncpus = 1L,
                              cl = NULL, seed = 1234, verbose = FALSE,
                              control = NULL, ...) {
  # rename for internal use
  meq.alt <- neq.alt
  
  # checks
  if (!inherits(object, "conMLM")) {
    stop("Restriktor ERROR: object must be of class conMLM.")
  }
  if (type != "global") {
    type <- toupper(type)
  }  
  if(!(type %in% c("A","B","global"))) {
    stop("Restriktor ERROR: type must be \"A\", \"B\", or \"global\"")
  }
  if(!(boot %in% c("no", "residual", "model.based", "parametric"))) {
    stop("Restriktor ERROR: boot method unknown.")
  }
  if (boot == "residual") {
    boot <- "model.based"
  }
  
  # original model
  model.org <- object$model.org
  # model matrix
  X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model.org$model[, attr(model.org$terms, "response")])
  # unconstrained df
  df.residual <- object$df.residual
  # parameter estimates
  b.unrestr <- object$b.unrestr
  b.restr <- object$b.restr
  b.eqrestr <- NULL
  b.restr.alt <- NULL
  # length parameter vector
  p <- length(b.unrestr)
  # variable names
  cnames <- colnames(b.unrestr)
  rnames <- rownames(b.unrestr)
  # constraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  control <- c(object$control, control)
  # remove duplicated elements from control list
  control <- control[!duplicated(control)]
  # get tolerance for control if exists
  tol <- ifelse(is.null(control$tol), sqrt(.Machine$double.eps), control$tol)
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  }
  
  # check for intercept                                          
  intercept <- any(attr(terms(model.org), "intercept"))
  if (type == "global") {
    if (intercept) { 
      AmatG <- cbind(rep(0, (p - 1)), diag(rep(1, p - 1))) 
    } else {
      AmatG <- diag(1, p)
      for (i in 1:p) {
        AmatG[i,i-1] <- -1
      }
      AmatG <- AmatG[-1,]
    }
    AmatX <- AmatG %*% (diag(rep(1, p)) - t(Amat) %*%            
                          MASS::ginv(Amat %*% t(Amat)) %*% Amat)
    
    if (all(abs(AmatX) < tol)) { 
      type <- "A"
      attr(type, "type.org") <- "global"
    } else {
      # remove all rows with only zeros
      AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
      rAmatX <- GaussianElimination(t(AmatX), tol = tol)
      AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
    }
    AmatG <- rbind(AmatX, Amat)
    bvecG <- c(rep(0, nrow(AmatX)), bvec)
    attr(Amat, "Amat.global") <- AmatG
    attr(bvec, "bvec.global") <- bvecG
  }
  
  if (type == "global") {  
    b.eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               Amat      = AmatG,
                               bvec      = bvecG, 
                               meq       = nrow(AmatG),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$qp$solution
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    b.eqrestr <- matrix(b.eqrestr, ncol = ncol(y))
      colnames(b.eqrestr) <- cnames
      rownames(b.eqrestr) <- rnames
    fitted0 <- X %*% b.eqrestr
    residuals0 <- y - fitted0
    object.eqrestr <- list(residuals = residuals0)
    
    ll0 <- con_loglik_lm(object.eqrestr)
    ll1 <- object$loglik
    
    Ts <- -2*(ll0 - ll1)
  } else if (type == "A") {
    b.eqrestr <- con_solver_lm(X         = X, 
                               y         = y, 
                               Amat      = Amat,
                               bvec      = bvec, 
                               meq       = nrow(Amat),
                               absval    = ifelse(is.null(control$absval), 1e-09, 
                                                  control$absval),
                               maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                  control$maxit))$qp$solution
    b.eqrestr[abs(b.eqrestr) < tol] <- 0L
    b.eqrestr <- matrix(b.eqrestr, ncol = ncol(y))
      colnames(b.eqrestr) <- cnames
      rownames(b.eqrestr) <- rnames
      
    fitted0 <- X %*% b.eqrestr
    residuals0 <- y - fitted0
    object.eqrestr <- list(residuals = residuals0)
    object.eqrestr$weights <- object$weights
    
    ll0 <- con_loglik_lm(object.eqrestr)
    ll1 <- object$loglik
    
    Ts <- -2*(ll0 - ll1)
  } else if (type == "B") {
    if (meq.alt == 0L) {
      
      ll0 <- object$loglik
      
      fitted1 <- X %*% object$b.unrestr
      residuals1 <- y - fitted1
      object.unrestr <- list(residuals = residuals1)
      object.unrestr$weights   <- object$weights
      ll1 <- con_loglik_lm(object.unrestr)
      
      Ts <- -2*(ll0 - ll1)
    }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq.alt > 0L && meq.alt <= meq) {
        b.restr.alt <- con_solver_lm(X         = X, 
                                     y         = y, 
                                     Amat      = Amat[1:meq.alt,,drop=FALSE],
                                     bvec      = bvec[1:meq.alt], 
                                     meq       = meq.alt,
                                     absval    = ifelse(is.null(control$absval), 1e-09, 
                                                        control$absval),
                                     maxit     = ifelse(is.null(control$maxit), 1e04, 
                                                        control$maxit))$qp$solution
        b.restr.alt[abs(b.restr.alt) < tol] <- 0L
        b.restr.alt <- matrix(b.restr.alt, ncol = ncol(y))
          colnames(b.restr.alt) <- cnames
          rownames(b.restr.alt) <- rnames
          
        ll0 <- object$loglik
        
        fitted1 <- X %*% b.restr.alt
        residuals1 <- y - fitted1
        object.restr.alt <- list(residuals = residuals1)
        object.restr.alt$weights <- object$weights
        ll1 <- con_loglik_lm(object.restr.alt)
        
        Ts <- -2*(ll0 - ll1)
      }
      else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  if (!(attr(object$wt.bar, "method") == "none") && boot == "no") { 
    wt.bar <- object$wt.bar
    pvalue <- con_pvalue_Fbar(wt.bar      = wt.bar, 
                              Ts.org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq.alt     = meq.alt)
    
    attr(pvalue, "wt.bar"                     ) <- as.numeric(wt.bar)
    attr(pvalue, "wt.bar.method"              ) <- attr(wt.bar, "method")
    attr(pvalue, "converged"                  ) <- attr(wt.bar, "converged")
    attr(pvalue, "convergence_crit"           ) <- attr(wt.bar, "convergence_crit")
    attr(pvalue, "total_bootstrap_draws"      ) <- attr(wt.bar, "total_bootstrap_draws")
    attr(pvalue, "error.idx"                  ) <- attr(wt.bar, "error.idx")
    attr(pvalue, "mix_weights_bootstrap_limit") <- attr(wt.bar, "mix_weights_bootstrap_limit")
    
    attr(pvalue, "wt_bar_chunk") <- attr(wt.bar, "wt_bar_chunk")
    attr(pvalue, "chunk_size"  ) <- attr(wt.bar, "chunk_size_org")
    attr(pvalue, "total_chunks") <- attr(wt.bar, "total_chunks")
    attr(pvalue, "chunk_iter"  ) <- attr(wt.bar, "chunk_iter")
    
  } else if (boot == "parametric") {
    if (!is.function(p.distr)) {
      p.distr <- get(p.distr, mode = "function")
    }
    arguments <- list(...)
    pnames <- names(formals(p.distr))
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    pm <- names(arguments)[pm > 0L]
    formals(p.distr)[pm] <- unlist(arguments[pm])
    
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts.org   = Ts, 
                                         type     = type, 
                                         test     = "LRT", 
                                         meq.alt  = meq.alt,
                                         R        = R, 
                                         p.distr  = p.distr,
                                         parallel = parallel,
                                         ncpus    = ncpus, 
                                         cl       = cl,
                                         seed     = seed, 
                                         verbose  = verbose)
  } else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based(object, 
                                          Ts.org   = Ts, 
                                          type     = type, 
                                          test     = "LRT",
                                          meq.alt  = meq.alt,
                                          R        = R, 
                                          parallel = parallel,
                                          ncpus    = ncpus, 
                                          cl       = cl,
                                          seed     = seed, 
                                          verbose  = verbose)
  } else {
    pvalue <- as.numeric(NA)
  }  
  
  # necessary for the print function
  if (!is.null(attr(type, "type.org"))) {
    type <- "global"
  }
  
  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq.alt     = meq.alt,
              iact        = object$iact,
              type        = type,
              test        = "LRT",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b.eqrestr   = b.eqrestr,
              b.unrestr   = b.unrestr,
              b.restr     = b.restr,
              b.restr.alt = b.restr.alt,
              R2.org      = object$R2.org,
              R2.reduced  = object$R2.reduced,
              boot        = boot,
              model.org   = model.org)
  
  class(OUT) <- "conTest"
  
  OUT
}