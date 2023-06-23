## input options
# 1. est + vcov
# 2. LL + PT
# 3. IC values (AIC, ORIC, GORIC, GORICA)
# 4. IC weights or (Bayesian) posterior model probs. 

# In case of an equal-evidence approach, aggregating evidence from, say, 5 studies 
# with n=100 observations is the same as obtaining evidence from 1 study 
# (as if it was possible) with n=500 observations (like meta-analysis does).
# In the added-evidence approach, the aggregated evidence from, says, 5 studies 
# is stronger than as if the data were combined (as if that was possible).

# TODO

# can gorica function do more than restriktor, see mail RK

# update VIGNETTE (ik krijg toegang tot .Rmd files)
# ik krijg alle files mbt evSyn om aan te passen ()
# nieuwe versie naar CRAN
# update restriktor website

evSyn <- function(object, ...) { UseMethod("evSyn") }

# -------------------------------------------------------------------------
## est (list + vec) + cov (list + mat)
## LL (list + vec) + PT (list + vec)
## IC weights (list + vec) rowSums = 1
## IC values (list + vec)

## object = est, VCOV = cov
## object = LL, PT = PT
## object = IC weights + rowsum check
## object = IC values
# -------------------------------------------------------------------------


evSyn <- function(object, ...) {
  
  arguments <- list(...)
  if (length(arguments)) {
    pnames <- c("VCOV", "PT", "hypotheses", "type", "comparison", "hypo_names")
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    if (any(pm == 0L)) { 
      pm.idx <- which(pm == 0L)
      stop("Restriktor Error: ", names(arguments[pm.idx]), " invalid argument(s).")
    }
  }

  # object must be a list
  if (!is.list(object)) {
    stop("Restriktor ERROR: object must be a list of numeric vectors.", call. = FALSE)
  }
  # object must be a list with numeric vectors
  obj_isnum <- sapply(object, is.numeric)
  if ( !any(obj_isnum) ) {
    stop("Restriktor ERROR: object must be a list of numeric vectors.", call. = FALSE)
  }
  
  VCOV <- arguments$VCOV
  PT   <- arguments$PT
  
  if (!is.null(VCOV) && !is.list(VCOV)) {
    stop("Restriktor ERROR: VCOV must be a list of matrices.", call. = FALSE)
  } else if (!is.null(VCOV) ) {
    # check if it is a list of matrices
    VCOV_isMat <- sapply(VCOV, is.matrix)
    if (!any(VCOV_isMat)) {
      stop("Restriktor ERROR: VCOV must be a list of matrices.", call. = FALSE)  
    }
  }

  if (!is.null(PT) && !is.list(PT)) {
    stop("Restriktor ERROR: PT must be a list of numeric vectors.", call. = FALSE)
  } else if ( !is.null(PT) ) {
    # check if it is a list of numeric vectors
    PT_isNum <- sapply(PT, is.numeric)
    if (!any(PT_isNum)) {
      stop("Restriktor ERROR: PT must be a list of numeric vectors.", call. = FALSE)  
    }
  }

  if (!is.null(VCOV) & !is.null(PT)) {
    stop("Restriktor ERROR: both VCOV and PT are found, which confuses me.")
  }
  
  # if all vectors sum to 1, object contains IC weights
  obj_isICweights <- FALSE
  if ( all(abs(sapply(object, sum) - 1) <= sqrt(.Machine$double.eps)) ) {
   obj_isICweights <- TRUE  
  } 
  
  
  ## est (vector) + VCOV (matrix)
  ## check if list contains a vector/matrix with the the parameter estimates
  if (!is.null(VCOV)) {
    evSyn.est(object, ...)
  } else if (!is.null(PT)) {
    # LL + PT
    evSyn.LL(object, ...)
  } else if (is.null(VCOV) & is.null(PT) & obj_isICweights) { 
    # IC weights
    if (!is.null(arguments$type) && arguments$type == "equal") {
      message("The added-evidence method is the only available approach when weights are included in the input.")
    }
    evSyn.ICweights(object, ...)
  } else if (is.null(VCOV) & is.null(PT) & !obj_isICweights) { 
    # IC values
    if (!is.null(arguments$type) && arguments$type == "equal") { 
      message("The added-evidence method is the only available approach when the input consists of GORIC(A) values.")
    }
    evSyn.ICvalues(object, ...)
  } else {
    stop("restriktor Error: I don't know how to handle the input.")
  }
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on the (standard) parameter estimates and the covariance matrix
evSyn.est <- function(object, VCOV = list(), hypotheses = list(),
                      type = c("equal", "added"), 
                      comparison = c("unconstrained", "complement", "none")) {
  
  if (missing(comparison)) 
    comparison <- "none"
  comparison <- match.arg(comparison)

  if (missing(type)) 
    type <- "added"
  type <- match.arg(type)

  # check if VCOV and hypotheses are both a non-empty list
  if ( !is.list(VCOV) & length(VCOV) == 0 ) {
    stop("Restriktor ERROR: VCOV must be a list of covariance matrices of the (standardized) parameter estimates of interest.", call. = FALSE)  
  } 

  if ( !is.list(hypotheses) & length(hypotheses) == 0 ) {
    stop("Restriktor ERROR: hypotheses must be a list.", call. = FALSE)  
  } 
  
  # number of primary studies
  S <- length(object)
  V <- length(VCOV)
  
  if (S != V) {
    stop("Restriktor ERROR: the number of items in the object list (i.e., number of (standardized) estimates) must equal the number of items in the VCOV list.", call. = FALSE)
  }

  # check if the matrices are all symmetrical
  VCOV_isSym <- sapply(VCOV, isSymmetric)
  if (!all(VCOV_isSym)) {
    sprintf("Restriktor ERROR: the %sth covariance matrix in VCOV is not symmetric.", which(!VCOV_isSym))
    stop(call. = FALSE)  
  }

  # if only 1 hypothesis is provided, that hypothesis will be applied to all studies
  # This also means that the vector with parameter estimates must have equal names between the studies.
  # check if the user specified the same hypothesis in each study. 
  # number of hypotheses (this should be equal in each study, see check below)
  
  # if hypotheses is a unnamed list (or any list element), I assume that all 
  # hypotheses must be applied to each study. Thus, len_H = 1
  # list_names <- names(hypotheses)
  # 
  # if (any(list_names == "") | length(list_names) == 0) {
  #   stop("Restriktor ERROR: The 'hypotheses' argument must be a named list. Please provide hypotheses in the following format: 
  #        'list(H1 = H1)' or 'list(S1 = list(H11, H12), S2 = list(H21, H22))'", call. = FALSE)
  #   
  # } 

  # determine if same set of hypotheses in each study or different set of hypotheses in each study.
  diffSet_hypotheses <- sapply(hypotheses, is.list)
  if (all(!diffSet_hypotheses)) {
    # same hypotheses in each study
    sameHypo <- TRUE  
  } else if (all(diffSet_hypotheses)) {
    # different hypotheses in each study
    sameHypo <- FALSE  
  } else {
    stop("Restriktor ERROR: To apply hypotheses, use the following format:\n",
         "If you want to apply the same set of hypotheses for each study, use: hypotheses = list(H1, H2, ...).\n",
         "If you want to apply a different set of hypotheses for each study, use: hypotheses = list(S1 = list(H1, H2), S2 = list(H3, H4)).",
         call. = FALSE)
    
  }
  
  ## for testing purposes only
  #hypo = list(H1 = list(hypothesis1), H2 = list(hypothesis1, hypothesis2))
  #hypo = list(H1 = list(hypothesis1, hypothesis2), H2 = list(hypothesis1, hypothesis2))  
  #hypo = list(H1 = list(hypothesis1), H2 = list(hypothesis2))  
  #hypo = list(H1 = hypothesis1)  
  
  # number of hypotheses must be equal for each studie. In each study a set of 
  # shared theories (i.e., hypotheses) are compared.
  len_H <- sapply(hypotheses, length)
  if (!sameHypo) {
    if (length(unique(len_H)) > 1) {
      stop("Restriktor ERROR: The number of hypotheses must be consistent across all studies.", call. = FALSE)
    }
    
    if (length(object) != length(len_H)) {
      stop("Restriktor ERROR: The number of hypotheses does not match the number of studies.", call. = FALSE)
    }
  }

  # if same hypo for each study, then only one hypo is allowed for comparison = complement
  comp_check_same <- sameHypo & length(len_H) == 1
  # if diff hypo for each study, then only one hypo is allowed for comparison = complement
  comp_check_diff <- !sameHypo & all(len_H == 1)
  
  
  if (comparison == "complement") {
    if ((sameHypo && !comp_check_same) | (!sameHypo && !comp_check_diff)) {
      warning("Restriktor Warning: Only one order-restricted hypothesis is allowed (for now) when comparison = 'complement'.",
              " Setting comparison to 'unconstrained' instead.", call. = FALSE)
      comparison <- "unconstrained"
    }
  }

  
  NrHypos <- length(len_H)
  NrHypos_incl <- NrHypos + 1
  if (comparison == "none"){
    NrHypos_incl <- NrHypos
  }

  LL <- PT <- weight_m <- GORICA_m <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
    rownames(LL) <- rownames(PT) <- paste0("Study ", 1:S)

  CumulativeGoricaWeights <- CumulativeGorica <- matrix(NA, nrow = S+1, ncol = NrHypos_incl)
  #  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(paste0("study_", 1:S), "final")
  
  sequence <- vector(mode = "character", length = S)
  for (i in 1:S) {
    if (i == 1) {
      sequence[i] <- "1"
    } else {
      sequence[i] <- paste0("1-", i)
    }
  }
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(paste0("Studies ", sequence), "Final")
     
  if (NrHypos == 1 & comparison == "complement") {
    hnames <- c("H1", "Hc")
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = 1)
  } else if (comparison == "none"){
    hnames <- c(paste0("H", 1:NrHypos))
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  } else {
    hnames <- c(paste0("H", 1:NrHypos), "Hu")
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  }
    
  colnames(GORICA_m) <- colnames(weight_m) <- colnames(LL) <- colnames(PT) <- colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- hnames
  rownames(GORICA_m) <- rownames(ratio.weight_mu) <- rownames(weight_m) <- paste0("Study ", 1:S)
  

  for(s in 1:S) {
    if (sameHypo) {
      res_goric <- goric(object[[s]], VCOV = VCOV[[s]], 
                         hypotheses = as.list(unlist(hypotheses)),
                         type = 'gorica', comparison = comparison)
    } else {
      res_goric <- goric(object[[s]], VCOV = VCOV[[s]],
                         hypotheses = hypotheses[[s]],
                         type = 'gorica', comparison = comparison)
    }
    
    if (comparison == "unconstrained") {
      ratio.weight_mu[s, ] <- res_goric$ratio.gw[, NrHypos_incl]
    } else if (comparison == "complement") {
      ratio.weight_mu[s, ] <- res_goric$ratio.gw[1, NrHypos_incl]
    } 
    
  
    LL[s, ] <- res_goric$result$loglik
    PT[s, ] <- res_goric$result$penalty
    GORICA_m[s, ] <- res_goric$result$gorica
    weight_m[s, ] <- res_goric$result$gorica.weights
  }
  
  sumPT <- sumLL <- 0
  if (type == "added") { 
    # added-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else { 
    # equal-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s,] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  }
  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  
  rownames(Final.ratio.GORICA.weights) <- hnames
  hypotheses <- res_goric$hypotheses_usr
  #names(hypotheses) <- hnames
  
  # Output
  if (NrHypos == 1 & comparison == "complement") {
    colnames(ratio.weight_mu) <- c("H1 vs. Hc1")
    colnames(Final.ratio.GORICA.weights) <- c("vs. H1", "vs. Hc")
    
    out <- list(type = type, 
                GORICA_m = GORICA_m, GORICA.weight_m = weight_m, 
                ratio.GORICA.weight_mc = ratio.weight_mu, LL_m = LL, PT_m = PT,
                Cumulative.GORICA = CumulativeGorica, 
                Cumulative.GORICA.weights = CumulativeGoricaWeights,
                Final.ratio.GORICA.weights = Final.ratio.GORICA.weights,
                hypotheses = hypotheses)
  } else if (comparison == "none") {
    colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. H", 1:NrHypos))
    
    out <- list(type = type,
                GORICA_m = GORICA_m, GORICA.weight_m = weight_m, LL_m = LL, PT_m = PT,
                Cumulative.GORICA = CumulativeGorica, 
                Cumulative.GORICA.weights = CumulativeGoricaWeights,
                Final.ratio.GORICA.weights = Final.ratio.GORICA.weights,
                hypotheses = hypotheses)
  } else { 
    # unconstrained
    colnames(ratio.weight_mu) <- c(paste0("H", 1:NrHypos, " vs. Unc."), "Unc. vs. Unc.")
    colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. H", 1:NrHypos), "vs. Hu")
    
    out <- list(type = type,
                GORICA_m        = GORICA_m, 
                GORICA.weight_m = weight_m, 
                ratio.GORICA.weight_mu = ratio.weight_mu, 
                LL_m = LL, 
                PT_m = PT,
                Cumulative.GORICA = CumulativeGorica, 
                Cumulative.GORICA.weights  = CumulativeGoricaWeights,
                Final.ratio.GORICA.weights = Final.ratio.GORICA.weights)
  }
  
  class(out) <- c("evSyn.est", "evSyn")
  
  return(out)
}




# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on log likelihood and penalty values
evSyn.LL <- function(object, PT = list(), type = c("added", "equal"),
                     hypo_names = NULL) {
  
  if (missing(type)) 
    type <- "added"
  type <- match.arg(type)
  
  # check if PT is a non-empty list
  if ( !is.list(PT) & length(PT) == 0 ) {
    stop("Restriktor ERROR: PT must be a list of penalty weights.", call. = FALSE)  
  } 
  
  
  LL <- object
  S <- length(LL)
  NrHypos <- length(LL[[1]]) - 1
  if (is.null(hypo_names)) {
    hnames <- paste0("H", 1:(NrHypos + 1))
  } else {
    hnames <- hypo_names
  }
  
  weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  
  LL <- do.call(rbind, LL)
  PT <- do.call(rbind, PT)
  #LL <- lapply(LL, function(x) setNames(x, hnames))
  #PT <- lapply(PT, function(x) setNames(x, hnames))
  
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(LL) <- colnames(PT) <- colnames(weight_m) <- hnames
  
  sequence <- vector(mode = "character", length = S)
  for (i in 1:S) {
    if (i == 1) {
      sequence[i] <- "1"
    } else {
      sequence[i] <- paste0("1-", i)
    }
  }
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(paste0("Studies ", sequence), "Final")
  rownames(LL) <- rownames(PT) <- rownames(weight_m) <- paste0("Study ", 1:S)
  
  sumLL <- 0
  sumPT <- 0
  IC <- -2 * LL + 2 * PT
  if (type == "added") { 
    # added-ev approach
    for(s in 1:S) {
      minIC <- min(IC[s, ])
      weight_m[s, ] <- exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC)))
      
      sumLL <- sumLL + LL[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- -2 * sumLL + 2 * sumPT
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- as.matrix(exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric))))
    }
  } else { 
    # equal-ev approach
    for (s in 1:S){
      sumLL <- sumLL + LL[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- -2 * sumLL + 2 * sumPT/s
      
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  }
  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S,]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S,]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  rownames(Final.ratio.GORICA.weights) <- hnames
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hnames)
  
  out <- list(type = type,
              LL_m = LL, 
              PT_m = PT, 
              GORICA_m = IC, 
              GORICA.weight_m   = weight_m,
              Cumulative.GORICA = CumulativeGorica, 
              Cumulative.GORICA.weights = CumulativeGoricaWeights,
              Final.ratio.GORICA.weights  = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn.LL", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA values
evSyn.ICvalues <- function(object, hypo_names = NULL) {
  
  IC <- object
  S  <- length(IC)
  NrHypos <- length(IC[[1]]) - 1
  weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  
  if (is.null(hypo_names)) {
    hnames <- paste0("H", 1:(NrHypos+1))
  } else {
    hnames <- hypo_names
  }
  
  IC <- do.call(rbind, IC)
  
  weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(IC) <- colnames(weight_m) <- hnames
  
  sequence <- vector(mode = "character", length = S)
  for (i in 1:S) {
    if (i == 1) {
      sequence[i] <- "1"
    } else {
      sequence[i] <- paste0("1-", i)
    }
  }
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(paste0("Studies ", sequence), "Final")
  #rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(paste0("Study_", 1:S), "Final")
  rownames(IC) <- rownames(weight_m) <- paste0("Study ", 1:S)
  
  sumIC <- 0
  for(s in 1:S){
    minIC <- min(IC[s, ])
    weight_m[s, ] <- exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC)))
  
    sumIC <- sumIC + IC[s, ]
    CumulativeGorica[s, ] <- sumIC
    minGoric <- min(CumulativeGorica[s, ])
    CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
  }

  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  
  rownames(Final.ratio.GORICA.weights) <- hnames
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hnames)
  
  out <- list(type              = "added",
              GORICA_m          = IC, 
              GORICA.weight_m   = weight_m,
              Cumulative.GORICA = CumulativeGorica, 
              Cumulative.GORICA.weights  = CumulativeGoricaWeights,
              Final.ratio.GORICA.weights = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn.ICvalues", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA weights or (Bayesian) posterior model probabilities
evSyn.ICweights <- function(object, priorWeights = NULL, hypo_names = NULL) {
  
  Weights <- object
  S <- length(object[[1]])
  NrHypos <- S - 1
  
  Weights <- do.call(rbind, Weights)
  
  if (is.null(priorWeights)) {
    priorWeights <- rep(1/(NrHypos + 1), (NrHypos + 1))
  }
  # To make it sum to 1 (if it not already did)
  priorWeights <- priorWeights / sum(priorWeights) 
  
  if (is.null(hypo_names)) {
    hypo_names <- paste0("H", 1:(NrHypos+1))
  }
  
  CumulativeWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  colnames(CumulativeWeights) <- hypo_names
  
  sequence <- vector(mode = "character", length = S)
  for (i in 1:S) {
    if (i == 1) {
      sequence[i] <- "1"
    } else {
      sequence[i] <- paste0("1-", i)
    }
  }
  rownames(CumulativeGorica) <- c(paste0("Studies ", sequence), "Final")
  #rownames(CumulativeWeights) <- c(paste0("Study_", 1:S), "Final")
  
  CumulativeWeights[1, ] <- priorWeights * Weights[1, ] / sum( priorWeights * Weights[1, ] )
  
  for(s in 2:S) {
    CumulativeWeights[s, ] <- CumulativeWeights[(s-1), ] * Weights[s, ] / sum( CumulativeWeights[(s-1), ] * Weights[s, ] )
  }
  
  
  CumulativeWeights[(S+1), ] <- CumulativeWeights[S, ]
  
  Final.weights <- CumulativeWeights[S, ]
  Final.ratio.GORICA.weights <- Final.weights %*% t(1/Final.weights)
  rownames(Final.ratio.GORICA.weights) <- hypo_names
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hypo_names)
  
  out <- list(type               = "added",
              weight_m           = Weights,
              Cumulative.weights = CumulativeWeights,
              Final.ratio.GORICA.weights  = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn.ICweights", "evSyn")
  
  return(out)
  
}



## print function
print.evSyn <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  # added or equal approach
  type <- x$type
  # make the first letter upper-case
  #Type <- paste(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type)), sep="")
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  
  #cat("\n", "restriktor: ", paste(type, "Evidence Synthesis results:\n"), sep = "")
  cat(paste(type, "Evidence Synthesis results:\n"), sep = "")

  cat("\nFinal ratio GORICA weights:\n")  
  print(apply(x$Final.ratio.GORICA.weights, c(1,2), function(x) 
    format(x, scientific = (abs(x) >= 1e3 | (abs(x) <= 1e-3)), digits = digits, nsmall = 3)), 
    print.gap = 2, quote = FALSE, right = TRUE)
  
  cat("\n")
  message(x$messages$mix_weights)
}


## summary function
summary.evSyn <- function(object, digits = max(3, getOption("digits") - 4), ...) {
  
  x <- object
  # added or equal approach
  type <- x$type
  
  # make the first letter upper-case
  #Type <- paste(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type)), sep="")
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  cat(paste(type, "Evidence Synthesis results:\n"), sep = "")
  
  cat("\nGORICA values:\n")  
  
  print(apply(x$GORICA_m, c(1,2), function(x) 
    format(x, scientific = (abs(x) >= 1e3 | (abs(x) <= 1e-3)), digits = digits, nsmall = 3)), 
    print.gap = 2, quote = FALSE, right = TRUE)
  cat("---\n")
  
  if (!all(is.na(x$GORICA.weight_m))) {
    cat("\nGORICA weights:\n")  
    print(apply(x$GORICA.weight_m, c(1,2), function(x) sprintf("%.3f", x)), 
          print.gap = 2, quote = FALSE, right = TRUE)
    cat("---\n")
  }
  
  if (!is.null(x$LL_m)) {
    cat("\nLog-likelihood values:\n")  
    print(apply(x$LL_m, c(1,2), function(x) 
      format(x, scientific = (abs(x) >= 1e3 | (abs(x) <= 1e-3)), digits = digits, nsmall = 3)), 
      print.gap = 2, quote = FALSE, right = TRUE)
    cat("---\n")
  }
  
  if (!is.null(x$LL_m)) {
    cat("\nCumulative Log-likelihood values:\n")  
    Cumulative.LL <- apply(x$LL_m, 2, cumsum)  
    
    S <- nrow(Cumulative.LL)
    sequence <- vector(mode = "character", length = S)
    for (i in 1:S) {
      if (i == 1) {
        sequence[i] <- "1"
      } else {
        sequence[i] <- paste0("1-", i)
      }
    }
    rownames(Cumulative.LL) <- paste0("Studies ", sequence)
    Cumulative.LL <- rbind(Cumulative.LL, Final = Cumulative.LL[S, ])

    print(apply(Cumulative.LL, c(1,2), function(x) 
      format(x, scientific = (abs(x) >= 1e3 | (abs(x) <= 1e-3)), digits = digits, nsmall = 3)), 
      print.gap = 2, quote = FALSE, right = TRUE)
    cat("---\n")
    
    cat("\nFinal ratio Log-likelihood values:\n")  
    Final.ratio.Cumulative.LL <- Cumulative.LL[S, ] %*% t(1/Cumulative.LL[S, ])
    rownames(Final.ratio.Cumulative.LL) <- colnames(Final.ratio.Cumulative.LL)
    colnames(Final.ratio.Cumulative.LL) <- paste0("vs. ", rownames(Final.ratio.Cumulative.LL))
    
      
    print(apply(Final.ratio.Cumulative.LL, c(1,2), function(x) 
      format(x, scientific = (abs(x) >= 1e3 | (abs(x) <= 1e-3)), digits = digits, nsmall = 3)), 
      print.gap = 2, quote = FALSE, right = TRUE)
    cat("---\n")
  }
  
  if (!is.null(x$PT_m)) {
    cat("\nPenalty term values:\n")  
    print(apply(x$PT_m, c(1,2), function(x) sprintf("%.3f", x)), 
          print.gap = 2, quote = FALSE, right = TRUE)
    cat("---\n")
  }
  
  cat("\nCumulative GORICA values:\n")  
  print(apply(x$Cumulative.GORICA, c(1,2), function(x) 
    format(x, scientific = (abs(x) >= 1e3 | (abs(x) <= 1e-3)), digits = digits, nsmall = 3)), 
    print.gap = 2, quote = FALSE, right = TRUE)
  cat("---\n")

  cat("\nCumulative GORICA weights:\n")  
  print(apply(x$Cumulative.GORICA.weights, c(1,2), function(x) sprintf("%.3f", x)), 
        print.gap = 2, quote = FALSE, right = TRUE)
  cat("---\n")

  if (!is.null(x$Final.ratio.GORICA.weights)) {
    cat("\nFinal ratio GORICA weights:\n")  
    print(apply(x$Final.ratio.GORICA.weights, c(1,2), function(x) 
      format(x, scientific = (abs(x) >= 1e3 | (abs(x) <= 1e-3)), digits = digits, nsmall = 3)), 
          print.gap = 2, quote = FALSE, right = TRUE)
    cat("---\n")
  }
  
  cat("\n")
  message(x$messages$mix_weights)
}





## plot function
plot.evSyn <- function(x, ...) {
  namesH <- colnames(x$GORICA_m)
  NrHypos_incl <- ncol(x$GORICA_m)
  S <- nrow(x$GORICA_m)
  Name_studies <- as.factor(1:S)
  
  weight_m <- x$GORICA.weight_m
  CumulativeGoricaWeights <- x$Cumulative.GORICA.weights
  
  # Create data frame for per study weights
  per_study_df <- data.frame(study = rep(Name_studies, NrHypos_incl),
                             weight = c(weight_m))
  per_study_df$weight_type <- "per study"
  
  # Create data frame for cumulative weights
  cumulative_df <- data.frame(study = rep(Name_studies, NrHypos_incl),
                              weight = c(CumulativeGoricaWeights[1:S, ]))
  cumulative_df$weight_type <- "cumulative"
  
  # Combine data frames
  plot_data <- rbind(per_study_df, cumulative_df)
  plot_data$variable <- rep(rep(namesH, each = S), times = 2)
  
  # Create plot
  ggplot(plot_data, aes(x = .data[['study']], 
                        y = .data[['weight']], 
                        shape    = .data[['weight_type']], 
                        linetype = .data[['weight_type']], 
                        color    = .data[['variable']])) +
    geom_point(size = 3) +
    geom_line(data = plot_data[plot_data[['weight_type']] == "cumulative", ], 
              aes(group = .data[['variable']]), linewidth = 1) +
    #scale_shape_manual(values = c(8, 16)) +
    #scale_linetype_manual(values = c("solid", "dashed")) +
    #scale_color_manual(values = c("#D81B60", "#1E88E5")) +
    theme(
      plot.margin = unit(c(1,1,1,1), "cm"),
      legend.position = "bottom",
      legend.margin = margin(t = 10, r = 0, b = 3, l = 0),
      #legend.title = element_text(size = 18),
      legend.key = element_blank(),
      legend.text = element_text(size = 12),
      axis.text.x  = element_text(size = 12), axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, vjust = -3), axis.title.y = element_text(size = 14, vjust = 5),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    scale_y_continuous(limits = c(0, 1), n.breaks = 10) +
    scale_x_discrete(expand = c(0, 0.05)) +
    labs(x = "Studies", y = "GORIC(A) weights",
         title = "GORIC(A) weights per study and cumulative",
         shape = "", color = "", linetype = "") +
    guides(shape = guide_legend(override.aes = list(linetype = 0)))
}


