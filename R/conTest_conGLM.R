### computes the F-bar, LRT-bar and score-bar test statistic ####
# Silvapulle and Sen (2005). Constrained statistical inference. Chapter 4.
conTestF.conGLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                            p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                            cl = NULL, seed = 1234, verbose = FALSE,
                            control = NULL, ...) {

  # rename for internal use
  meq_alt <- neq.alt
  
  # checks
  if (!inherits(object, "conGLM")) {
    stop("Restriktor ERROR: object must be of class conGLM.")
  }
  if (type != "global") {
    type <- toupper(type)
  }    
  if(!(type %in% c("A","B","global"))) {
    stop("Restriktor ERROR: type must be \"A\", \"B\", or \"global\"")
  }
  if(!(boot %in% c("no","residual","model.based","parametric","mix.weights"))) {
    stop("Restriktor ERROR: boot method unknown.")
  }
  if (boot == "residual") {
    boot <- "model.based"
  }
  
  # original model
  model_org <- object$model_org
  # model matrix
  X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # sample size
#  n <- dim(X)[1]
  # weights
  w <- weights(model_org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- vcov(model_org) 
  # parameter estimates
  b_unrestr <- object$b_unrestr
  b_restr <- object$b_restr
  b_eqrestr <- NULL
  b_restr_alt <- NULL
  # length parameter vector
  p <- length(b_unrestr)
  # variable names
  vnames <- names(b_unrestr)
  # constraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  #control
  control <- object$control
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  }
  
  if (is.null(control$tol)) {
    tol <- sqrt(.Machine$double.eps)
  } else {
    tol <- control$tol
  }
  
  # check for intercept                                          
  intercept <- any(attr(terms(model_org), "intercept"))
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
      attr(type, "type_org") <- "global"
    } else {
      # remove all rows with only zeros
      AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
      rAmatX <- GaussianElimination(t(AmatX), tol = tol)
      AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
    }
    AmatG <- rbind(AmatX, Amat)
    bvecG <- c(rep(0, nrow(AmatX)), bvec)
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
  }
  
  if (type == "global") {
    CALL <- list(object = model_org, constraints = AmatG, rhs = bvecG, 
                 neq = nrow(AmatG), se = "none", Wt = FALSE)
    glm_fit <- do.call("restriktor", CALL)
    
    b_eqrestr <- coef(glm_fit)
    # compute global test statistic
    Ts <- c(t(b_restr - b_eqrestr) %*% solve(Sigma, b_restr - b_eqrestr))
  } else if (type == "A") {
    CALL <- list(object = model_org, constraints = Amat, rhs = bvec, 
                 neq = nrow(Amat), se = "none", Wt = FALSE)
    glm_fit <- do.call("restriktor", CALL)
    
    b_eqrestr <- coef(glm_fit)
    # compute test statistic for hypothesis test type A
    Ts <- c(t(b_restr - b_eqrestr) %*% solve(Sigma, b_restr - b_eqrestr))
  } else if (type == "B") {
    if (meq_alt == 0L) {
      # compute test statistic for hypothesis test type B when no equalities are
      # preserved in the alternative hypothesis.
      Ts <- c(t(b_unrestr - b_restr) %*% solve(Sigma, b_unrestr - b_restr))
    } else {
      if (meq_alt > 0L && meq_alt <= meq) {
        # compute test statistic for hypothesis test type B when some equalities may 
        # be preserved in the alternative hypothesis.
        CALL <- list(object = model_org, constraints = Amat[1:meq_alt,,drop = FALSE], 
                     rhs = bvec[1:meq_alt], neq = meq_alt, 
                     se = "none", Wt = FALSE)
        glm_fit <- do.call("restriktor", CALL)
        
        b_restr_alt <- coef(glm_fit)
        
        Ts <- c(t(b_restr - b_restr_alt) %*% solve(Sigma, b_restr - b_restr_alt))
      } else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  # The test statistics based on inequality constraints are often
  # distributed as mixtures of chi-squares. These mixing weights can be computed
  # using the multivariate normal distribution with additional Monte Carlo steps
  # or via bootstrapping. The pvalue can also be computed directly via 
  # the parametric bootstrap or model based bootstrap, without fist computing 
  # the mixing weights.
  
  if (!is.null(object$wt) && boot == "no") {
    wt <- object$wt
    # is this fool proof? 
    # The number of bootstrap samples must be large enough to avoid spurious results.
    wt <- rev(wt)
    if (attr(object$wt, "bootWt")) {
      if (attr(object$wt, "bootWt.R") < 999) {
        stop("Restriktor ERROR: increase the number of bootstrap draws. Preferably to a large number e.g., bootWt.R = 99999")
      }
      wt.idx <- which(wt == 0)
      wt <- wt[-wt.idx]
    }
    
    pvalue <- con_pvalue_Fbar(wt          = wt, 
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts_org   = Ts, 
                                         type     = type, 
                                         test     = "F", 
                                         meq_alt  = meq_alt, 
                                         R        = R, 
                                         p.distr  = p.distr,
                                         df       = df, 
                                         parallel = parallel,
                                         ncpus    = ncpus, 
                                         cl       = cl,
                                         seed     = seed, 
                                         verbose  = verbose)
  } else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based(object, 
                                          Ts_org   = Ts, 
                                          type     = type, 
                                          test     = "F", 
                                          meq_alt  = meq_alt,
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
  if (!is.null(attr(type, "type_org"))) {
    type <- "global"
  }
  
  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq_alt     = meq_alt,
              iact        = object$iact,
              type        = type,
              test        = "F",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b_eqrestr   = b_eqrestr,
              b_unrestr   = b_unrestr,
              b_restr     = b_restr,
              b_restr_alt = b_restr_alt,
              Sigma       = Sigma,
              R2_org      = object$R2_org,
              R2_reduced  = object$R2_reduced,
              boot        = boot,
              model_org   = model_org)
  
  OUT <- list(OUT)
  names(OUT) <- type
  
  class(OUT) <- "conTest"
  
  OUT
}


# REF: Silvapulle and Sen (2005). Constrained statistical inference. Chapter 4.
conTestLRT.conGLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                              p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                              cl = NULL, seed = 1234, verbose = FALSE,
                              control = NULL, ...) {
  
  # rename for internal use
  meq_alt <- neq.alt
  
  # checks
  if (!inherits(object, "conGLM")) {
    stop("Restriktor ERROR: object must be of class conGLM.")
  }
  if (type != "global") {
    type <- toupper(type)
  }  
  if(!(type %in% c("A","B","global"))) {
    stop("Restriktor ERROR: type must be \"A\", \"B\", or \"global\"")
  }
  if(!(boot %in% c("no", "residual", "model.based", "parametric", "mix.weights"))) {
    stop("Restriktor ERROR: boot method unknown.")
  }
  if (boot == "residual") {
    boot <- "model.based"
  }
  
  # original model
  model_org <- object$model_org
  # model matrix
  #X <- model.matrix(object)[,,drop=FALSE]
  # response variable
  #y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # weights
  #w <- weights(model_org)
  # unconstrained df
  df.residual <- object$df.residual
  # unconstrained covariance matrix
  Sigma <- vcov(model_org) 
  # parameter estimates
  b_unrestr <- object$b_unrestr
  b_restr <- object$b_restr
  b_eqrestr <- NULL
  b_restr_alt <- NULL
  # length parameter vector
  p <- length(b_unrestr)
  # variable names
  #vnames <- names(b_unrestr)
  # constraints stuff
  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq
  #control
  control <- object$control
  
  # check for equalities only
  if (meq == nrow(Amat)) {
    stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  }
  
  if (is.null(control$tol)) {
    tol <- sqrt(.Machine$double.eps)
  } else {
    tol <- control$tol
  }
  
  # check for intercept                                          
  intercept <- any(attr(terms(model_org), "intercept"))
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
      attr(type, "type_org") <- "global"
    } else {
      # remove all rows with only zeros
      AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
      rAmatX <- GaussianElimination(t(AmatX), tol = tol)
      AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
    }
    AmatG <- rbind(AmatX, Amat)
    bvecG <- c(rep(0, nrow(AmatX)), bvec)
    attr(Amat, "Amat_global") <- AmatG
    attr(bvec, "bvec_global") <- bvecG
  }
  
  if (type == "global") {  
    CALL <- list(object = model_org, constraints = AmatG, 
                 rhs = bvecG, neq = nrow(AmatG), se = "none", Wt = FALSE)
    glm_fit <- do.call("restriktor", CALL)
    
    ll0 <- glm_fit$loglik
    ll1 <- object$loglik
    
    Ts <- -2*(ll0 - ll1)
  } else if (type == "A") {
    CALL <- list(object = model_org, constraints = Amat, 
                 rhs = bvec, neq = nrow(Amat), se = "none", Wt = FALSE)
    glm_fit <- do.call("restriktor", CALL)
    
    ll0 <- glm_fit$loglik
    ll1 <- object$loglik
    
    Ts <- -2*(ll0 - ll1)
  } else if (type == "B") {
    if (meq_alt == 0L) {
      
      ll0 <- object$loglik
      ll1 <- logLik(model_org)
      
      Ts <- -2*(ll0 - ll1)
    }
    else {
      # some equality may be preserved in the alternative hypothesis.
      if (meq_alt > 0L && meq_alt <= meq) {
        CALL <- list(object = model_org, constraints = Amat[1:meq_alt,,drop=FALSE], 
                     rhs = bvec[1:meq_alt], neq = meq_alt, se = "none", Wt = FALSE)
        glm_fit <- do.call("restriktor", CALL)
        
        ll0 <- glm_fit$loglik
        ll1 <- object$loglik
        
        Ts <- -2*(ll0 - ll1)
      }
      else {
        stop("Restriktor ERROR: neq.alt must not be larger than neq.")
      }
    }
  } 
  
  if (!is.null(object$wt) && boot == "no") { 
    wt <- object$wt
    # is this fool proof? 
    # The number of bootstrap samples must be large enough to avoid spurious results.
    wt <- rev(wt)
    if (attr(object$wt, "bootWt")) {
      if (attr(object$wt, "bootWt.R") < 999) {
        stop("Restriktor ERROR: increase the number of bootstrap draws. Preferably to a large number e.g., bootWt.R = 99999")
      }
      wt.idx <- which(wt == 0)
      wt <- wt[-wt.idx]
    }
    
    pvalue <- con_pvalue_Fbar(wt          = wt, 
                              Ts_org      = Ts, 
                              df.residual = df.residual, 
                              type        = type,
                              Amat        = Amat, 
                              bvec        = bvec, 
                              meq         = meq, 
                              meq_alt     = meq_alt)
  } else if (boot == "parametric") {
    pvalue <- con_pvalue_boot_parametric(object, 
                                         Ts_org   = Ts, 
                                         type     = type, 
                                         test     = "LRT", 
                                         meq_alt  = meq_alt,
                                         R        = R, 
                                         p.distr  = p.distr,
                                         df       = df, 
                                         parallel = parallel,
                                         ncpus    = ncpus, 
                                         cl       = cl,
                                         seed     = seed, 
                                         verbose  = verbose)
  } else if (boot == "model.based") {
    pvalue <- con_pvalue_boot_model_based(object, 
                                          Ts_org   = Ts, 
                                          type     = type, 
                                          test     = "LRT",
                                          meq_alt  = meq_alt,
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
  if (!is.null(attr(type, "type_org"))) {
    type <- "global"
  }
  
  OUT <- list(CON         = object$CON,
              Amat        = Amat,
              bvec        = bvec,
              meq         = meq,
              meq_alt     = meq_alt,
              iact        = object$iact,
              type        = type,
              test        = "LRT",
              Ts          = Ts,
              df.residual = df.residual,
              pvalue      = pvalue,
              b_eqrestr   = b_eqrestr,
              b_unrestr   = b_unrestr,
              b_restr     = b_restr,
              b_restr_alt = b_restr_alt,
              Sigma       = Sigma,
              R2_org      = object$R2_org,
              R2_reduced  = object$R2_reduced,
              boot        = boot,
              model_org   = model_org)
  
  OUT <- list(OUT)
  names(OUT) <- type
  
  class(OUT) <- "conTest"
  
  OUT
}



# REF: Robertson, Silvapulle and Sen (2005, chapter 4)
conTestScore.conGLM <- function(object, type = "A", neq.alt = 0, boot = "no", R = 9999, 
                               p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                               cl = NULL, seed = 1234, verbose = FALSE,
                               control = NULL, ...) {
  
  class(object) <- c("restriktor", "conGLM", "conLM")
  conTestScore.conLM(object, ...)
  # rename for internal use
  # meq_alt <- neq.alt
  # 
  # # checks
  # if (!("conGLM" %in% class(object))) {
  #   stop("Restriktor ERROR: object must be of class conGLM.")
  # }
  # if (type != "global") {
  #   type <- toupper(type)
  # }
  # if(!(type %in% c("A","B","global"))) {
  #   stop("Restriktor ERROR: type must be \"A\", \"B\", or \"global\"")
  # }
  # if(!(boot %in% c("no", "residual", "model.based", "parametric", "mix.weights"))) {
  #   stop("Restriktor ERROR: boot method unknown.")
  # }
  # if (boot == "residual") {
  #   boot <- "model.based"
  # }
  # 
  # # original model
  # model_org <- object$model_org
  # # model matrix
  # X <- model.matrix(object)[,,drop=FALSE]
  # # response variable
  # y <- as.matrix(model_org$model[, attr(model_org$terms, "response")])
  # # unconstrained df
  # df.residual <- object$df.residual
  # # unconstrained covariance matrix
  # Sigma <- vcov(model_org)
  # # sample size
  # n <- dim(X)[1]
  # # number of parameters
  # p <- dim(X)[2]
  # # weights
  # w <- weights(model_org)
  # if (is.null(w)) {
  #   w <- rep(1, n)
  # }
  # #W <- diag(w)
  # # parameter estimates
  # b_unrestr <- object$b_unrestr
  # b_restr <- object$b_restr
  # b_eqrestr <- NULL
  # b_restr_alt <- NULL
  # # length parameter vector
  # p <- length(b_unrestr)
  # # variable names
  # vnames <- names(b_unrestr)
  # # restraints stuff
  # Amat <- object$constraints
  # bvec <- object$rhs
  # meq  <- object$neq
  # #control
  # control <- object$control
  # 
  # # check for equalities only
  # if (meq == nrow(Amat)) {
  #   stop("Restriktor ERROR: test not applicable for object with equality restrictions only.")
  # }
  # 
  # if (is.null(control$tol)) {
  #   tol <- sqrt(.Machine$double.eps)
  # } else {
  #   tol <- control$tol
  # }
  # 
  # # check for intercept
  # intercept <- any(attr(terms(model_org), "intercept"))
  # if (type == "global") {
  #   if (intercept) {
  #     AmatG <- cbind(rep(0, (p - 1)), diag(rep(1, p - 1)))
  #   } else {
  #     AmatG <- diag(1, p)
  #     for (i in 1:p) {
  #       AmatG[i,i-1] <- -1
  #     }
  #     AmatG <- AmatG[-1,]
  #   }
  #   AmatX <- AmatG %*% (diag(rep(1, p)) - t(Amat) %*%
  #                         MASS::ginv(Amat %*% t(Amat)) %*% Amat)
  # 
  #   if (all(abs(AmatX) < tol)) {
  #     type <- "A"
  #     attr(type, "type_org") <- "global"
  #   } else {
  #     # remove all rows with only zeros
  #     AmatX  <- AmatX[!rowSums(abs(AmatX) < tol) == p,, drop = FALSE]
  #     rAmatX <- GaussianElimination(t(AmatX), tol = tol)
  #     AmatX  <- AmatX[rAmatX$pivot,, drop = FALSE]
  #   }
  #   AmatG <- rbind(AmatX, Amat)
  #   bvecG <- c(rep(0, nrow(AmatX)), bvec)
  #   attr(Amat, "Amat_global") <- AmatG
  #   attr(bvec, "bvec_global") <- bvecG
  # }
  # 
  # if (type == "global") {
  #   b_eqrestr <- con_solver_lm(X         = X,
  #                           y         = y,
  #                           b_unrestr = b_unrestr,
  #                           w         = w,
  #                           Amat      = AmatG,
  #                           bvec      = bvecG,
  #                           meq       = nrow(AmatG),
  #                           absval    = ifelse(is.null(control$absval), 1e-09,
  #                                              control$absval),
  #                           maxit     = ifelse(is.null(control$maxit), 1e04,
  #                                              control$maxit))$solution
  #   b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol),
  #                                     sqrt(.Machine$double.eps),
  #                                     control$tol)] <- 0L
  #   names(b_eqrestr) <- vnames
  # 
  #   res0 <- y - X %*% b_eqrestr
  #   res1 <- residuals(object)
  #   wres0 <- as.vector(res0) * weights(object, "working")
  #   wres1 <- as.vector(res1) * weights(object, "working")
  # 
  #   dispersion0 <- if (substr(model_org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  #   else sum(wres0^2, na.rm = TRUE) / sum(weights(object, "working"), na.rm = TRUE)
  # 
  #   dispersion1 <- if (substr(model_org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  #   else sum(wres1^2, na.rm = TRUE) / sum(weights(object, "working"), na.rm = TRUE)
  # 
  #   G0 <- colSums(wres0 * X / dispersion0)
  #   G1 <- colSums(wres1 * X / dispersion1)
  #   s20 <- sum((res0)^2) / (n - (p - nrow(Amat)))
  #   I0 <- 1/s20 * (t(X) %*% X)
  #   Ts <- (G0 - G1) %*% solve(I0 , (G0 - G1))
  #   ###############################################
  # } else if (type == "A") {
  #   b_eqrestr <- con_solver_lm(X         = X,
  #                           y         = y,
  #                           b_unrestr = b_unrestr,
  #                           w         = w,
  #                           Amat      = Amat,
  #                           bvec      = bvec,
  #                           meq       = nrow(Amat),
  #                           absval    = ifelse(is.null(control$absval), 1e-09,
  #                                              control$absval),
  #                           maxit     = ifelse(is.null(control$maxit), 1e04,
  #                                              control$maxit))$solution
  #   b_eqrestr[abs(b_eqrestr) < ifelse(is.null(control$tol),
  #                                     sqrt(.Machine$double.eps),
  #                                     control$tol)] <- 0L
  #   names(b_eqrestr) <- vnames
  # 
  #   res0 <- y - X %*% b_eqrestr
  #   res1 <- residuals(object)
  #   wres0 <- as.vector(res0) * weights(object, "working")
  #   wres1 <- as.vector(res1) * weights(object, "working")
  # 
  #   dispersion0 <- if (substr(model_org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  #   else sum(wres0^2, na.rm = TRUE) / sum(weights(object, "working"), na.rm = TRUE)
  # 
  #   dispersion1 <- if (substr(model_org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  #   else sum(wres1^2, na.rm = TRUE) / sum(weights(object, "working"), na.rm = TRUE)
  # 
  #   G0 <- colSums(wres0 * X / dispersion0)
  #   G1 <- colSums(wres1 * X / dispersion1)
  #   s20 <- sum((res0)^2) / (n - (p - nrow(Amat)))
  #   I0 <- 1/s20 * (t(X) %*% X)
  #   Ts <- (G0 - G1) %*% solve(I0 , (G0 - G1))
  #   ###############################################
  # } else if (type == "B") {
  #   if (meq_alt == 0L) {
  #     res0 <- residuals(object)
  #     res1 <- residuals(model_org)
  #     wres0 <- as.vector(res0) * weights(object, "working")
  #     wres1 <- as.vector(res1) * weights(object, "working")
  # 
  #     dispersion0 <- if (substr(model_org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  #     else sum(wres0^2, na.rm = TRUE) / sum(weights(object, "working"), na.rm = TRUE)
  # 
  #     dispersion1 <- if (substr(model_org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  #     else sum(wres1^2, na.rm = TRUE) / sum(weights(object, "working"), na.rm = TRUE)
  # 
  #     G0 <- colSums(wres0 * X / dispersion0)
  #     G1 <- colSums(wres1 * X / dispersion1)
  #     s20 <- sum((res0)^2) / n - (p - qr(Amat[0:meq,])$rank)
  #     I0 <- 1/dispersion0 * (t(X) %*% X)
  #     Ts <- (G0 - G1) %*% solve(I0 , (G0 - G1))
  #     ###############################################
  #   }
  #   else {
  #     # some equality may be preserved in the alternative hypothesis.
  #     if (meq_alt > 0L && meq_alt <= meq) {
  #       b_restr_alt <- con_solver_lm(X         = X,
  #                                 y         = y,
  #                                 b_unrestr = b_unrestr,
  #                                 w         = w,
  #                                 Amat      = Amat[1:meq_alt,,drop=FALSE],
  #                                 bvec      = bvec[1:meq_alt],
  #                                 meq       = meq_alt,
  #                                 absval    = ifelse(is.null(control$absval), 1e-09,
  #                                                    control$absval),
  #                                 maxit     = ifelse(is.null(control$maxit), 1e04,
  #                                                    control$maxit))$solution
  #       b_restr_alt[abs(b_restr_alt) < ifelse(is.null(control$tol),
  #                                             sqrt(.Machine$double.eps),
  #                                             control$tol)] <- 0L
  #       names(b_restr_alt) <- vnames
  # 
  #       res0 <- residuals(object)
  #       res1 <- y - X %*% b_restr_alt
  #       wres0 <- as.vector(res0) * weights(object, "working")
  #       wres1 <- as.vector(res1) * weights(object, "working")
  # 
  #       dispersion0 <- if (substr(model_org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  #       else sum(wres0^2, na.rm = TRUE) / sum(weights(object, "working"), na.rm = TRUE)
  # 
  #       dispersion1 <- if (substr(model_org$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  #       else sum(wres1^2, na.rm = TRUE) / sum(weights(object, "working"), na.rm = TRUE)
  # 
  #       G0 <- colSums(wres0 * X / dispersion0)
  #       G1 <- colSums(wres1 * X / dispersion1)
  #       s20 <- sum((res0)^2) / n - (p - qr(Amat[0:meq,])$rank)
  #       I0 <- 1/s20 * (t(X) %*% X)
  #       Ts <- (G0 - G1) %*% solve(I0 , (G0 - G1))
  #     }
  #     else {
  #       stop("Restriktor ERROR: neq.alt must not be larger than neq.")
  #     }
  #   }
  # }
  # 
  # if (!is.null(object$wt) && boot == "no") {
  #   wt <- object$wt
  #   wt <- rev(wt)
  #   # is this fool proof?
  #   # The number of bootstrap samples must be large enough to avoid spurious results.
  #   if (attr(object$wt, "bootWt")) {
  #     if (attr(object$wt, "bootWt.R") < 999) {
  #       stop("Restriktor ERROR: increase the number of bootstrap draws. Preferably to a large number e.g., bootWt.R = 99999")
  #     }
  #     wt.idx <- which(wt == 0)
  #     wt <- wt[-wt.idx]
  #   }
  # 
  #   pvalue <- con_pvalue_Fbar(wt          = wt,
  #                             Ts_org      = Ts,
  #                             df.residual = df.residual,
  #                             type        = type,
  #                             Amat        = Amat,
  #                             bvec        = bvec,
  #                             meq         = meq,
  #                             meq_alt     = meq_alt)
  # } else if (boot == "parametric") {
  #   pvalue <- con_pvalue_boot_parametric(object,
  #                                        Ts_org   = Ts,
  #                                        type     = type,
  #                                        test     = "score",
  #                                        meq_alt  = meq_alt,
  #                                        R        = R,
  #                                        p.distr  = p.distr,
  #                                        df       = df,
  #                                        parallel = parallel,
  #                                        ncpus    = ncpus,
  #                                        cl       = cl,
  #                                        seed     = seed,
  #                                        verbose  = verbose)
  # } else if (boot == "model.based") {
  #   pvalue <- con_pvalue_boot_model_based(object,
  #                                         Ts_org   = Ts,
  #                                         type     = type,
  #                                         test     = "score",
  #                                         meq_alt  = meq_alt,
  #                                         R        = R,
  #                                         parallel = parallel,
  #                                         ncpus    = ncpus,
  #                                         cl       = cl,
  #                                         seed     = seed,
  #                                         verbose  = verbose)
  # } else {
  #   pvalue <- as.numeric(NA)
  # }
  # 
  # # necessary for the print function
  # if (!is.null(attr(type, "type_org"))) {
  #   type <- "global"
  # }
  # 
  # OUT <- list(CON         = object$CON,
  #             Amat        = Amat,
  #             bvec        = bvec,
  #             meq         = meq,
  #             meq_alt     = meq_alt,
  #             iact        = object$iact,
  #             type        = type,
  #             test        = "Score",
  #             Ts          = Ts,
  #             df.residual = df.residual,
  #             pvalue      = pvalue,
  #             b_eqrestr   = b_eqrestr,
  #             b_unrestr   = b_unrestr,
  #             b_restr     = b_restr,
  #             b_restr_alt = b_restr_alt,
  #             Sigma       = Sigma,
  #             R2_org      = object$R2_org,
  #             R2_reduced  = object$R2_reduced,
  #             boot        = boot,
  #             model_org   = model_org)
  # 
  # OUT <- list(OUT)
  # names(OUT) <- type
  # 
  # class(OUT) <- "conTest"
  # 
  # OUT
}
