## Author: Leonard Vanbrabant
## Last updated: 31 Janurary 2023
## function to implement 2-stage EM


two_stage_sandwich <- function(object, ...) {
  H       <- object$H
  delta   <- object$delta
  satACOV <- object$satACOV
  
  acov <- satACOV * object$N
  
  meat  <- H %*% delta
  bread <- solve(t(delta) %*% meat) 
  VCOV  <- bread %*% t(meat) %*% satACOV %*% meat %*% bread
  
  idx <- grepl("~ 1", row.names(VCOV))
  
  VCOV <- VCOV[idx, idx]
  VCOV[upper.tri(VCOV)] <- t(VCOV)[upper.tri(VCOV)]
  
  VCOV
}



## (ordered) factors need to be specified in the data, not in the model
EM <- function(object, emControl = list(), auxiliary = c(), ...) {

    #warn <- getOption("warn")
    #options(warn = -1)
    
    # arguments for emNorm()
    emControl <- as.list(emControl)
    
    #print(emControl)
    #print(auxiliary)
    
    if (!is.list(emControl)) {
      stop("Restriktor ERROR: emControl must be a named list.") 
    }
    
    # get original data set with missings
    df_amp <- eval(object$call$data)
    
    # model names: 
    modelNames <- all.vars(terms(object)) 
    
    # add auxiliary variables to data
    df_amp <- df_amp[, c(modelNames, auxiliary)]
    
    ## (ordered) factors are lost during emNorm, we restore them below
    # which are (ordered) factors
    isFactor  <- sapply(df_amp, is.factor)
    isOrdered <- sapply(df_amp, is.ordered)
    # (ordered) factor names
    fnames  <- names(isFactor[isFactor])
    ofnames <- names(isOrdered[isOrdered])
    # order levels
    oflevels <- lapply(df_amp[, ofnames], levels)
    # reference levels 
    
    # if a factor is coded e.g., 0/1 then emNorm recodes it to 1/2. To avoid this
    # unwanted recoding, we make all columns numeric.
    df_amp <- apply(df_amp, 2, as.numeric)
    
    # EM algorithm for incomplete multivariate normal data 
    param <- norm2::emNorm(df_amp, 
                           iter.max        = ifelse(is.null(emControl$iter.max), 1000, emControl$iter.max),
                           estimate.worst  = ifelse(is.null(emControl$estimate.worst), TRUE, emControl$estimate.worst),
                           prior           = ifelse(is.null(emControl$prior), "uniform", emControl$prior),
                           criterion       = emControl$criterion,
                           prior.df        = emControl$prior.df,
                           prior.sscp      = emControl$prior.sscp,
                           starting.values = emControl$starting.values)
    
    # extract imputed data set
    df_imp <- param$y.mean.imp
    df_imp <- as.data.frame(df_imp)
    
    # round (ordered) factors. NB. order of variables changes here.
    df_imp <- cbind(df_imp[, !isFactor, drop = FALSE], round(df_imp[, isFactor, drop = FALSE]))
    
    # restore (ordered) factors
    df_imp[, fnames] <- lapply(df_imp[, fnames, drop = FALSE], as.factor) 
    
    # restore ordered levels
    for (vnames in ofnames) {
      df_imp[[vnames]] <- factor(df_imp[[vnames]], ordered = TRUE, 
                                      levels = oflevels[[vnames]])   
    }
    
    #options(warn = warn)
    
    return(df_imp)
}



# -------------------------------------------------------------------------
# Victoria Savalei & Peter M. Bentler (2009) A Two-Stage Approach to
# Missing Data: Theory and Application to Auxiliary Variables, Structural Equation Modeling: A
# Multidisciplinary Journal, 16:3, 477-497, DOI: 10.1080/10705510903008238

computeDelta <- function(object = NULL) {
  
  obj_terms <- terms(object)
  # check if model has an intercept
  intrcpt <- attr(obj_terms, "intercept")
  # check for interaction effects
  #label_terms <- attr(obj_terms, "term.labels")
  
  ## get model-matrix 
  mm <- model.matrix(object)[, -1]
  ## get dependent variable name
  dep.name <- all.vars(formula(object))[1]
  ## get model-frame 
  mf <- model.frame(object)[, dep.name]  
  
  ## add outcome variable
  df_imp <- data.frame(mf, mm)
    names(df_imp)[1] <- dep.name 
    
  # deal with interactions
  #interactions <- label_terms[grep(":", label_terms)]
  #interactions_df <- lapply(strsplit(interactions, ":"), FUN = function(x) Reduce("*", df_imp[, x]))
  #  names(interactions_df) <- interactions
  
  # if (length(interactions) > 0L) {
  #   df_imp <- data.frame(object$model, interactions_df)
  # }
  
  # replace colon by dot  
  #label_terms <- gsub(":", ".", label_terms)
  
  # subset model variables
  #df_imp <- df_imp[, c(all.vars(obj_terms)[1], label_terms)]
  
  p <- length(coef(object))
  # no intercept
  if (intrcpt == 0L) {
    p <- p + 1
  }
  
  
  GLIST <- list()
  
  GLIST$lambda <- diag(p)
  # intercept only model
  if (p == 1) {
    GLIST$lambda <- matrix(, 1,0)  
  }
  GLIST$theta <- matrix(0, p, p)
  if (p == 1) {
    GLIST$theta <- sigma(object)^2
  }
  # include interactions
  GLIST$psi <- lavaan::lav_matrix_bdiag(
    sigma(object)^2,
    cov(as.matrix(df_imp)[, -1])#[, label_terms]))
  )
  if (p == 1) {
    GLIST$psi <- matrix(, 0, 0)
  }
  
  if (p > 1) {
    beta <- matrix(0, p, p)
    if (intrcpt == 0L) {
      beta[1, 2:p] <- coef(object)
    } else {
      beta[1, 2:p] <- coef(object)[-1]
    }
    GLIST$beta <- beta
  }  
  
  GLIST$nu <- matrix(0, p, 1)
  if (p == 1) {
    GLIST$nu <- coef(object)
  }
  
  if (p > 1) {
    alpha <- matrix(0, p, 1)
    if (intrcpt == 0L) {
      alpha[2:p, 1] <- colMeans(df_imp)[-1]
    } else {
      alpha[1, 1  ] <- coef(object)[1]
      alpha[2:p, 1] <- colMeans(df_imp)[-1]
    }
  } else {
    alpha <- matrix(, 0, 1)
  }  
  GLIST$alpha <- alpha
  
  
  ###
  #lavmodel <- fitTarget@Model
  #lavmodel@GLIST
  
  nvar <- p
  num.idx <- seq_len(p)
  
  pstar <- as.integer(nvar * (nvar + 1)/2)
  pstar <- nvar + pstar
  
  NCOL <- 1.5*p + 0.5*p^2#lavmodel@nx.free
  
  m.el.idx <- x.el.idx <- vector("list", length = length(GLIST))
  
  ## check intercept only model
  m.el.idx <- list()
  for (gl in names(GLIST)) {
    if (gl %in% c("lambda", "theta", "nu")) {
      m.el.idx[[gl]] <- integer(0)
    } else {
      m.el.idx[[gl]] <- which(abs(GLIST[[gl]]) > 0)
    }
  }
  
  ## intercept model
  # predict mean y (alpha[1])
  # regressies (beta)
  # (co)varianties (psi)
  # means (alpha[-1])
  

  if (!is.na(GLIST$alpha[1])) {
    x.el.idx <- vector("list", 6)
    x.el.idx[[1]] <- integer(0)
    x.el.idx[[2]] <- integer(0)
    
    x.el.idx_ul <- unlist(m.el.idx)
    alpha_len <- length(x.el.idx_ul[grepl("alpha" , names(x.el.idx_ul))])
    beta_len  <- length(x.el.idx_ul[grepl("beta"  , names(x.el.idx_ul))])
    
    if (GLIST$alpha[1] != 0) {
      ## intercept model
      psi_start <- beta_len + 2
    } else {
      ## no intercept model
      psi_start <- beta_len + 1
    }  
    
    # Remove zero elements from the matrix
    matrix_without_zero <- GLIST$psi[GLIST$psi != 0]
    # Create a vector of unique elements
    unique_elements <- unique(matrix_without_zero)
    
    # Create a vector of integer values starting from psi_start
    integer_values <- psi_start:(psi_start + length(unique_elements - 1))
    
    # Use the match function to match each non-zero element in the matrix with its corresponding integer value
    matrix_with_integer_values <- GLIST$psi
    matrix_with_integer_values[GLIST$psi != 0] <- integer_values[match(matrix_without_zero, unique_elements)]
    
    x.el.idx[[3]] <- lavaan::lav_matrix_vec(matrix_with_integer_values[matrix_with_integer_values != 0])
    
    
    if (GLIST$alpha[1] != 0) {
      x.el.idx[[4]] <- 2:(1+beta_len)
    } else {
      x.el.idx[[4]] <- 1:beta_len
    }
    
    
    psi_beta_max <- max(c(x.el.idx[[3]], x.el.idx[[4]]))
    if (GLIST$alpha[1] != 0) {
      x.el.idx[[6]] <- c(1, seq(from = psi_beta_max+1, length.out = alpha_len-1))
    } else {
      x.el.idx[[6]] <- seq(from = psi_beta_max+1, length.out = alpha_len)
    }
    
    
    for (mm in 1:length(GLIST)) {
      #lavmodel@m.free.idx
      #lavmodel@x.free.idx
      dix <- duplicated(x.el.idx[[mm]])
      if (any(dix)) {
        m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
        x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
      }
    }
    
    nmat <- 6
    Delta.group <- matrix(0, nrow = pstar, ncol = NCOL)
    mm.in.group <- 1:nmat + cumsum(c(0, nmat))[1]
    
    for (mm in mm.in.group) {
      mname <- names(GLIST)[mm]
      if (!length(m.el.idx[[mm]])) 
        next
      DELTA <- derivative.sigma.LISREL(m = mname, 
                                       idx = m.el.idx[[mm]], 
                                       MLIST = GLIST[mm.in.group], 
                                       delta = TRUE)
      
      DELTA.mu <- derivative.mu.LISREL(m = mname, 
                                       idx = m.el.idx[[mm]], 
                                       MLIST = GLIST[mm.in.group])
      DELTA <- rbind(DELTA.mu, DELTA)
      
      Delta.group[, x.el.idx[[mm]]] <- DELTA
    }
    
    Delta <- Delta.group
    
  } else {
    Delta <- diag(2)
  }
  
  # remove zero columns
  Delta <- Delta[, apply(Delta, 2, function(x) !all(x==0))]
  
  
  Delta
}

# -------------------------------------------------------------------------


## This is an adapted version of the function twostageMatrices() written by 
##  Terrence D. Jorgensen from the semTools package.

two_stage_matrices <- function(object, auxiliary = c(), emControl = list(), ...) {
  
  
  obj_terms <- terms(object)
  # check for interaction effects
  #label_terms <- attr(obj_terms, "term.labels")
  
  aux <- auxiliary
  #arguments <- list(...)
  #arguments$auxiliary <- aux
  
  # EM algorithm for incomplete multivariate normal data 
  data_imp <- EM(object     = object,
                 emControl  = emControl, 
                 auxiliary  = aux)
  
  fitTarget <- update(object, data = data_imp)
  
  # # add interactions to the data inherited from the lm model
  # interactions <- label_terms[grep(":", label_terms)]
  # interactions_df <- lapply(strsplit(interactions, ":"), FUN = function(x) Reduce("*", data_imp[, x]))
  #   names(interactions_df) <- interactions
  # 
  # if (length(interactions) > 0L) {
  #   data_imp <- data.frame(data_imp, interactions_df)
  # }
  # 
  # ## here we only select the variables that are part of the model
  # label_terms <- gsub(":", ".", label_terms)
  # data_imp <- data_imp[, c(all.vars(obj_terms)[1], label_terms, aux)]

  ## get model-matrix. This also includes new columns in case of interactions (but no dep.var)
  mm <- model.matrix(fitTarget)[, -1]
  ## get dependent variable name
  dep.name <- all.vars(formula(fitTarget))[1]
  ## get model-frame 
  mf <- model.frame(fitTarget)[, dep.name, drop = FALSE] # drop = FALSE is needed to keep column name  
  ## add outcome variable
  data_imp <- data.frame(mf, mm, data_imp[, aux, drop = FALSE])
  
  # create saturated model
  sat_vnames <- colnames(data_imp) #c(varnames, aux)
  
  # check for duplicate object and aux names
  sat_vnames_dup <- sat_vnames[duplicated(sat_vnames)]
  
  if (length(sat_vnames_dup) > 0L) {
    stop("restriktor ERROR: duplicated variable name(s) found: ", sQuote(sat_vnames_dup))
  }
  
  covstruc <- outer(sat_vnames, sat_vnames, function(x, y) paste(x, "~~", y))
  Satmodel_vnames <- c(covstruc[upper.tri(covstruc, diag = TRUE)], paste(sat_vnames, "~ 1"))
  
  # ## fit saturated model
  # fitSat <- lavaan::lavaan(Satmodel_vnames, data = data_imp, fixed.x = FALSE, meanstructure = TRUE, 
  #                       conditional.x = FALSE, missing = "fiml")
  # # asymptotic information and covariance matrices of saturated model
  # satACOV <- lavaan::vcov(fitSat)
  
  # cov.ov
  Sigma <- cov(data_imp[, sat_vnames])
  
  p <- ncol(Sigma)
  n <- nrow(data_imp)
  #p. <- p*(p+1)/2
  
  V <- kronecker(Sigma, Sigma)
  
  # select the correct columns from V
  idx <- c()
  for(i in 0:(p-1)){
    idx <- c(idx, (i*p)+(i+1):p)
  }
  
  my_fun <- function(i, j, k, l) { (1 + (Sigma[i,l] * Sigma[j,k]) / ( Sigma[i,k] * Sigma[j,l] )) }
  eg  <- expand.grid(i = 1:p, j = 1:p, k = 1:p, l = 1:p)
  out <- do.call(mapply, c(my_fun, eg))
  
  R <- matrix(out, ncol(V), ncol(V))
  m1 <- R * V
  m1 <- m1[idx, idx]
  
  satACOV <- lav_matrix_bdiag(m1, Sigma) / n
    row.names(satACOV) <- Satmodel_vnames
    colnames(satACOV ) <- Satmodel_vnames
  
  satInfo <- solve(satACOV * n)
  
  # # model derivatives
  # delta <- lavaan::lavInspect(fitTarget_lav, "delta")
  delta <- computeDelta(fitTarget)
  
  # -------------------------------------------------------------------------
  #acov <- satACOV * lavaan::nobs(fitSat)
  #2*acov[4:5, 4:5]^2
  
  # -------------------------------------------------------------------------
  ## extract parameter table satuarated model
  #PTsat <- lavaan::parTable(fitSat)
  
  # fit target model
  # form <- formula(object)
  # targetModel <- Reduce(paste, deparse(form))
  # # unconstrained target model: replace with lm()
  # fitTarget_lav <- lavaan::sem(targetModel, data = data_imp, meanstructure = TRUE, 
  #                              fixed.x = FALSE)
  
  target_vnames <- sat_vnames[!sat_vnames %in% aux]
  target_covstruc <- outer(target_vnames, target_vnames, function(x, y) paste(x, "~~", y))
  targetModel_vnames <- c( paste(target_vnames, "~ 1"), target_covstruc[upper.tri(target_covstruc, diag = TRUE)])
  
  colnames(delta)  <- targetModel_vnames
  row.names(delta) <- targetModel_vnames
  # -------------------------------------------------------------------------
  #undebug(lavaan:::computeDelta)
  #undebug(lavaan:::derivative.sigma.LISREL)
  # categorical?
  # -------------------------------------------------------------------------
  
  covparams  <- grep(pattern = "~~", x = rownames(delta))
  meanparams <- grep(pattern = "~ 1", x = rownames(delta))
  delta <- delta[c(covparams, meanparams), ]
  
  ## extract estimated moments from saturated model, and number of moments
  satSigma <- Sigma #lavaan::lavInspect(fitSat, "cov.ov") 
  satMu    <- colMeans(data_imp[, sat_vnames]) #lavaan::lavInspect(fitSat, "mean.ov")
  
  ## in case of auxiliaries
  if (!is.null(aux)) {
    auxNames <- aux
    satSigma <- satSigma[target_vnames, target_vnames]
    satMu <- satMu[target_vnames]
  }
  
  p <- length(satMu)
  pStar <- p*(p + 1) / 2
  
  ## extract model-implied moments
  muHat    <- satMu #lavaan::lavInspect(fitTarget, "mean.ov")
  sigmaHat <- Sigma #lavaan::lavInspect(fitTarget, "cov.ov")
  
  muHat    <- muHat[names(satMu)]
  sigmaHat <- sigmaHat[rownames(satSigma), colnames(satSigma)]
  shinv    <- solve(sigmaHat)
  
  ## assemble complete-data information matrix
  H <- matrix(0, (pStar + p), (pStar + p))
  
  ## observed (non-zero in case of aux-variables)
  dMu <- satMu - muHat
  
  H[1:pStar, 1:pStar                ] <- lav_matrix_duplication_pre_post(shinv %x% (shinv %*% (satSigma + dMu %*% t(dMu)) %*% shinv - .5*shinv))
  H[(pStar + 1):(pStar + p), 1:pStar] <- lav_matrix_duplication_post(shinv %x% (t(dMu) %*% shinv))
  H[1:pStar, (pStar + 1):(pStar + p)] <- t(H[(pStar + 1):(pStar + p), 1:pStar])
  H[(pStar + 1):(pStar + p), (pStar + 1):(pStar + p)] <- shinv
  
  ## in case of auxiliaries
  if (!is.null(aux)) {
    #dimTar  <- !(PTsat$lhs %in% auxNames | PTsat$rhs %in% auxNames)
    #dimAux  <- PTsat$lhs %in% auxNames | PTsat$rhs %in% auxNames
    dimTar <- !grepl(paste(as.list(auxNames), collapse = "|"), Satmodel_vnames, perl = TRUE)
    dimAux <-  grepl(paste(as.list(auxNames), collapse = "|"), Satmodel_vnames, perl = TRUE)
    
    infoTar <- satInfo[dimTar, dimTar]
    infoAux <- satInfo[dimAux, dimAux]
    infoAT  <- satInfo[dimAux, dimTar]
    satInfo <- infoTar - t(infoAT) %*% solve(infoAux) %*% infoAT
    satACOV <- solve(satInfo) / n
  }
  
  OUT <- list(delta = delta, H = H, satACOV = satACOV, satInfo = satInfo, 
              fitTarget = fitTarget, N = n)
  
  OUT
}




# -------------------------------------------------------------------------
# acknowledgement: 
# the functions below are taken from the lavaan package (in agreement with Yves Rosseel)
derivative.sigma.LISREL <- function (m = "lambda", idx = seq_len(length(MLIST[[m]])), 
                                     MLIST = NULL, vech = TRUE, delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  PSI <- MLIST$psi
  v.idx <- lav_matrix_vech_idx(nvar)
  pstar <- nvar * (nvar + 1)/2
  if (m == "nu" || m == "alpha" || m == "tau" || m == "gamma" || 
      m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0, nrow = pstar, ncol = length(idx)))
  }
  delta.flag <- FALSE
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- MLIST$delta
    delta.flag <- TRUE
  }
  else if (m == "delta") {
    return(matrix(0, nrow = pstar, ncol = length(idx)))
  }
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  }
  else {
    IB.inv <- internal_get_IB.inv(MLIST = MLIST)
  }
  if (m == "lambda" || m == "beta") {
    L1 <- LAMBDA %*% IB.inv %*% PSI %*% t(IB.inv)
  }
  if (m == "beta" || m == "psi") {
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }
  if (m == "lambda") {
    KOL.idx <- matrix(1:(nvar * nfac), nvar, nfac, byrow = TRUE)[idx]
    DX <- (L1 %x% diag(nvar))[, idx, drop = FALSE] + (diag(nvar) %x% 
                                                        L1)[, KOL.idx, drop = FALSE]
  }
  else if (m == "beta") {
    KOL.idx <- matrix(1:(nfac * nfac), nfac, nfac, byrow = TRUE)[idx]
    DX <- (L1 %x% LAMBDA..IB.inv)[, idx, drop = FALSE] + 
      (LAMBDA..IB.inv %x% L1)[, KOL.idx, drop = FALSE]
    DX[, which(idx %in% lav_matrix_diag_idx(nfac))] <- 0
  }
  else if (m == "psi") {
    DX <- (LAMBDA..IB.inv %x% LAMBDA..IB.inv)
    lower.idx <- lav_matrix_vech_idx(nfac, diagonal = FALSE)
    upper.idx <- lav_matrix_vechru_idx(nfac, diagonal = FALSE)
    offdiagSum <- DX[, lower.idx] + DX[, upper.idx]
    DX[, c(lower.idx, upper.idx)] <- cbind(offdiagSum, offdiagSum)
    DX <- DX[, idx, drop = FALSE]
  }
  else if (m == "theta") {
    DX <- matrix(0, nvar * nvar, length(idx))
    DX[cbind(idx, seq_along(idx))] <- 1
  }
  else if (m == "delta") {
    Omega <- computeSigmaHat.LISREL(MLIST, delta = FALSE)
    DD <- diag(DELTA[, 1], nvar, nvar)
    DD.Omega <- (DD %*% Omega)
    A <- DD.Omega %x% diag(nvar)
    B <- diag(nvar) %x% DD.Omega
    DX <- A[, lav_matrix_diag_idx(nvar), drop = FALSE] + 
      B[, lav_matrix_diag_idx(nvar), drop = FALSE]
    DX <- DX[, idx, drop = FALSE]
  }
  else {
    stop("wrong model matrix names: ", m, "\n")
  }
  if (delta.flag && !m == "delta") {
    DX <- DX * as.vector(DELTA %x% DELTA)
  }
  if (vech) {
    DX <- DX[v.idx, , drop = FALSE]
  }
  DX
}


derivative.mu.LISREL <- function(m = "alpha", idx = seq_len(length(MLIST[[m]])), 
                                 MLIST = NULL) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  if (m == "gamma" || m == "psi" || m == "theta" || m == "tau" || 
      m == "delta" || m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0, nrow = nvar, ncol = length(idx)))
  }
  if (is.null(MLIST$alpha)) 
    ALPHA <- matrix(0, nfac, 1L)
  else ALPHA <- MLIST$alpha
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  }
  else {
    IB.inv <- internal_get_IB.inv(MLIST = MLIST)
  }
  if (m == "nu") {
    DX <- diag(nvar)
  }
  else if (m == "lambda") {
    DX <- t(IB.inv %*% ALPHA) %x% diag(nvar)
  }
  else if (m == "beta") {
    DX <- t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
    DX[, lav_matrix_diag_idx(nfac)] <- 0
  }
  else if (m == "alpha") {
    DX <- LAMBDA %*% IB.inv
  }
  else {
    stop("wrong model matrix names: ", m, "\n")
  }
  DX <- DX[, idx, drop = FALSE]
  DX
}



internal_get_IB.inv <- function(MLIST = NULL) {
  BETA <- MLIST$beta
  nr <- nrow(MLIST$psi)
  if (!is.null(BETA)) {
    tmp <- -BETA
    tmp[lav_matrix_diag_idx(nr)] <- 1
    IB.inv <- solve(tmp)
  }
  else {
    IB.inv <- diag(nr)
  }
  IB.inv
}

lav_matrix_diag_idx <- function(n = 1L) {
  1L + (seq_len(n) - 1L) * (n + 1L)
}

computeSigmaHat.LISREL <- function(MLIST = NULL, delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PSI <- MLIST$psi
  THETA <- MLIST$theta
  BETA <- MLIST$beta
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  }
  else {
    IB.inv <- internal_get_IB.inv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }
  VYx <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + 
    THETA
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- diag(MLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    VYx <- DELTA %*% VYx %*% DELTA
  }
  VYx
}
