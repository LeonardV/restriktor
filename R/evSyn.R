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

#evSyn <- function(object, ...) { UseMethod("evSyn") }

evSyn_est       <- function(object, ...) UseMethod("evSyn_est")
evSyn_LL        <- function(object, ...) UseMethod("evSyn_LL")
evSyn_ICweights <- function(object, ...) UseMethod("evSyn_ICweights")
evSyn_ICvalues  <- function(object, ...) UseMethod("evSyn_ICvalues")
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
    evSyn_est(object, ...)
  } else if (!is.null(PT)) {
    # LL + PT
    evSyn_LL(object, ...)
  } else if (is.null(VCOV) & is.null(PT) & obj_isICweights) { 
    # IC weights
    if (!is.null(arguments$type) && arguments$type == "equal") {
      message("The added-evidence method is the only available approach when weights are included in the input.")
    }
    evSyn_ICweights(object, ...)
  } else if (is.null(VCOV) & is.null(PT) & !obj_isICweights) { 
    # IC values
    if (!is.null(arguments$type) && arguments$type == "equal") { 
      message("The added-evidence method is the only available approach when the input consists of GORIC(A) values.")
    }
    evSyn_ICvalues(object, ...)
  } else {
    stop("restriktor Error: I don't know how to handle the input.")
  }
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on the (standard) parameter estimates and the covariance matrix
evSyn_est.list <- function(object, ..., VCOV = list(), hypotheses = list(),
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
    stop(sprintf("Restriktor ERROR: the %sth covariance matrix in VCOV is not symmetric.", which(!VCOV_isSym)), call. = FALSE)  
  }

  # if only 1 hypothesis is provided, that hypothesis will be applied to all studies
  # This also means that the vector with parameter estimates must have equal names between the studies.
  # check if the user specified the same hypothesis in each study. 
  # number of hypotheses (this should be equal in each study, see check below)
  
  # determine if same set of hypotheses in each study or different set of hypotheses in each study.
  diffSet_hypotheses <- sapply(hypotheses, is.list)
  if (all(!diffSet_hypotheses)) {
    # same hypotheses in each study
    sameHypo <- TRUE  
  } else if (all(diffSet_hypotheses)) {
    # different hypotheses in each study
    sameHypo <- FALSE  
  } else {
    stop("Restriktor ERROR: use the following hypotheses format:\n",
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

  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Studies 1"
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
    
  if (NrHypos == 1 & comparison == "complement") {
    exist_hnames <- names(hypotheses)
    if (!is.null(exist_hnames)) {
      exist_hnames <- c(exist_hnames, "Complement")
      hnames_idx <- exist_hnames != ""
    } else {
      exist_hnames <- vector("character", 2L)
      hnames_idx <- exist_hnames != ""
    }
    hnames <- c("H1", "Complement")
    hnames_idx <- exist_hnames != ""
    exist_hnames[!hnames_idx] <- hnames[!hnames_idx]
    hnames <- exist_hnames
    names(hypotheses) <- hnames[-max(length(hnames))] # remove Hc
    
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = 1)
  } else if (comparison == "none"){
    # existing hypo names
    exist_hnames <- names(hypotheses)
    if (!is.null(exist_hnames)) {
      hnames_idx <- exist_hnames != ""
    } else {
      exist_hnames <- vector("character", length(hypotheses))
      hnames_idx <- exist_hnames != ""
    }
    hnames <- c(paste0("H", 1:NrHypos))
    exist_hnames[!hnames_idx] <- hnames[!hnames_idx]
    hnames <- exist_hnames
    names(hypotheses) <- hnames
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  } else {
    exist_hnames <- names(hypotheses)
    if (!is.null(exist_hnames)) {
      exist_hnames <- c(exist_hnames, "Unconstrained")
      hnames_idx <- exist_hnames != ""
    } else {
      exist_hnames <- vector("character", length(hypotheses) + 1L)
      hnames_idx <- exist_hnames != ""
    }
    hnames <- c(paste0("H", 1:NrHypos), "Unconstrained")
    hnames_idx <- exist_hnames != ""
    exist_hnames[!hnames_idx] <- hnames[!hnames_idx]
    hnames <- exist_hnames
    names(hypotheses) <- hnames[-max(length(hnames))] # remove Hu
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
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  }
  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  
  rownames(Final.ratio.GORICA.weights) <- hnames
  hypotheses <- res_goric$hypotheses_usr
  
  # Output
  if (NrHypos == 1 & comparison == "complement") {
    colnames(ratio.weight_mu) <- c(paste0(names(hypotheses), " vs. ", "Complement"))
    colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
    
    out <- list(type = type, 
                hypotheses = hypotheses,
                GORICA_m = GORICA_m, 
                GORICA_weight_m = weight_m, 
                ratio_GORICA_weight_mc = ratio.weight_mu, 
                LL_m = LL, PT_m = PT,
                Cumulative_GORICA = CumulativeGorica, 
                Cumulative_GORICA_weights = CumulativeGoricaWeights,
                Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  } else if (comparison == "none") {
    colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
    
    out <- list(type = type,
                hypotheses = hypotheses,
                GORICA_m = GORICA_m, 
                GORICA_weight_m = weight_m, 
                LL_m = LL, PT_m = PT,
                Cumulative_GORICA = CumulativeGorica, 
                Cumulative_GORICA_weights = CumulativeGoricaWeights,
                Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  } else { 
    # unconstrained
    colnames(ratio.weight_mu) <- c(paste0(colnames(CumulativeGorica), " vs. ", "Unconstrained"))
    colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
    
    out <- list(type = type,
                hypotheses = hypotheses,
                GORICA_m = GORICA_m, 
                GORICA_weight_m = weight_m, 
                ratio_GORICA_weight_mu = ratio.weight_mu, 
                LL_m = LL, 
                PT_m = PT,
                Cumulative_GORICA = CumulativeGorica, 
                Cumulative_GORICA_weights = CumulativeGoricaWeights,
                Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  }
  
  class(out) <- c("evSyn.est", "evSyn")
  
  return(out)
}




# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on log likelihood and penalty values
evSyn_LL.list <- function(object, ..., PT = list(), type = c("equal", "added"),
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
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(LL) <- colnames(PT) <- colnames(weight_m) <- hnames
  
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Studies 1"
  
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
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
    for (s in 1:S) {
      minIC <- min(IC[s, ])
      weight_m[s, ] <- exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC)))
      
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
              GORICA_weight_m   = weight_m,
              Cumulative_GORICA = CumulativeGorica, 
              Cumulative_GORICA_weights = CumulativeGoricaWeights,
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn.LL", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA values
evSyn_ICvalues.list <- function(object, ..., hypo_names = NULL) {
  
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
  
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Studies 1"
  
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
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
              GORICA_weight_m   = weight_m,
              Cumulative_GORICA = CumulativeGorica, 
              Cumulative_GORICA_weights  = CumulativeGoricaWeights,
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn.ICvalues", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA weights or (Bayesian) posterior model probabilities
evSyn_ICweights.list <- function(object, ..., priorWeights = NULL, hypo_names = NULL) {
  
  Weights <- object
  S <- length(Weights)
  Weights <- do.call(rbind, Weights)
  NrHypos <- ncol(Weights)
  
  if (is.null(priorWeights)) {
    priorWeights <- rep(1/(NrHypos), (NrHypos))
  }
  # To make it sum to 1 (if it not already did)
  priorWeights <- priorWeights / sum(priorWeights) 
  
  if (is.null(hypo_names)) {
    hypo_names <- paste0("H", 1:(NrHypos)) 
  }
  
  CumulativeWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos))
  colnames(Weights) <- colnames(CumulativeWeights) <- hypo_names
  
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Studies 1"
  
  rownames(CumulativeWeights) <- c(sequence, "Final")
  CumulativeWeights[1, ] <- priorWeights * Weights[1, ] / sum( priorWeights * Weights[1, ] )
  
  for(s in 2:S) {
    CumulativeWeights[s, ] <- CumulativeWeights[(s-1), ] * Weights[s, ] / sum( CumulativeWeights[(s-1), ] * Weights[s, ] )
  }
  
  CumulativeWeights[(S+1), ] <- CumulativeWeights[S, ]
  
  Final.weights <- CumulativeWeights[S, ]
  Final.ratio.GORICA.weights <- Final.weights %*% t(1/Final.weights)
  rownames(Final.ratio.GORICA.weights) <- hypo_names
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hypo_names)
  
  out <- list(type                       = "added",
              GORICA_weight_m            = Weights,
              Cumulative_GORICA_weights  = CumulativeWeights,
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn.ICweights", "evSyn")
  
  return(out)
  
}


print.evSyn <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  # added or equal approach
  type <- x$type
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  cat(paste(type, "Evidence Synthesis results:\n"), sep = "")
  
  if (!is.null(x$Cumulative_GORICA_weights)) {
    cat("\nFinal GORICA weights:\n") 
    cgw <- sapply(x$Cumulative_GORICA_weights["Final", , drop = FALSE], 
                  FUN = function(x) format_numeric(x, digits = digits))
    names(cgw) <- colnames(x$Cumulative_GORICA_weights)
    print(cgw, print.gap = 2, quote = FALSE, right = TRUE)
    cat("---\n")
  }
  
  cat("\nRatio final GORICA weights:\n")  
  print(apply(x$Final_ratio_GORICA_weights, c(1,2), function(x) format_numeric(x, digits = digits)), 
        print.gap = 2, quote = FALSE, right = TRUE)
  cat("\n")
  
  message(x$messages$mix_weights)
}



summary.evSyn <- function(object, ...) {
  x <- object
  class(x) <- NULL
  
  ans <- list(
    type = x$type,
    n_studies = nrow(x$GORICA_weight_m[,, drop = FALSE]),
    hypotheses = x$hypotheses,
    GORICA_weight_m = x$GORICA_weight_m,
    GORICA_m = x$GORICA_m,
    LL_m = x$LL_m,
    PT_m = x$PT_m,
    Cumulative_GORICA_weights = x$Cumulative_GORICA_weights,
    Cumulative_GORICA = x$Cumulative_GORICA
  )
  
  if (!is.null(x$LL_m)) {
    sequence    <- paste0("Studies 1-", 1:ans$n_studies)
    sequence[1] <- "Studies 1"
    Cumulative_LL <- apply(x$LL_m, 2, cumsum)
    Cumulative_LL <- matrix(Cumulative_LL, nrow = nrow(x$LL_m), 
                            dimnames = list(sequence, colnames(x$LL_m)))
    ans$Cumulative_LogLik <- Cumulative_LL
  }
  
  if (!is.null(x$PT_m)) {
    sequence    <- paste0("Studies 1-", 1:ans$n_studies)
    sequence[1] <- "Studies 1"
    Cumulative_PT <- apply(x$PT_m[,, drop = FALSE], 2, cumsum)  
    Cumulative_PT <- matrix(Cumulative_PT, nrow = nrow(x$PT_m), 
                            dimnames = list(sequence, colnames(x$PT_m)))
    if (x$type == "equal") {
      Cumulative_PT <- Cumulative_PT / seq_len(ans$n_studies)
    } 
    ans$Cumulative_PT <- Cumulative_PT
  }
  
  final <- c()
  
  if (!is.null(x[["Cumulative_GORICA_weights"]])) {
    fcgw <- t(x[["Cumulative_GORICA_weights"]]["Final", ])
    rownames(fcgw) <- "GORICA weights"
    colnames(fcgw)  <- colnames(x[["Cumulative_GORICA_weights"]])
    final <- rbind(final, fcgw)
  }
  
  if (!is.null(x[["Cumulative_GORICA"]])) {
    fcgv <- t(x[["Cumulative_GORICA"]]["Final", ])
    rownames(fcgv) <- "GORICA values"
    final <- rbind(final, fcgv)
  }
  
  if (!is.null(ans$Cumulative_LogLik)) {
    fcllv <- t(ans$Cumulative_LogLik[ans$n_studies, ])
    rownames(fcllv) <- "Log-likelihood values"
    final <- rbind(final, fcllv)
  }
  
  if (!is.null(ans$Cumulative_PT)) { 
    fcptv <- t(ans$Cumulative_PT[ans$n_studies, ])
    rownames(fcptv) <- "Penalty term values"
    final <- rbind(final, fcptv)
  }
  
  ans$Final_Cumulative_results <- final
  
  if (!is.null(x$LL_m)) {
    Final_ratio_Cumulative_LL <- ans$Cumulative_LogLik[ans$n_studies, ] %*% t(1/ans$Cumulative_LogLik[ans$n_studies, ])
    rownames(Final_ratio_Cumulative_LL) <- colnames(ans$Cumulative_LogLik)
    colnames(Final_ratio_Cumulative_LL) <- paste0("vs. ", colnames(ans$Cumulative_LogLik))
    ans$Final_ratio_Cumulative_LL <- Final_ratio_Cumulative_LL  
  }
  
  ans$Ratio_GORICA_weight_mu <- x$ratio_GORICA_weight_mu
  ans$Ratio_GORICA_weight_mc <- x$ratio_GORICA_weight_mc
  ans$Final_ratio_GORICA_weights <- x$Final_ratio_GORICA_weights
  
  ans$messages <- list(mix_weights = x$messages$mix_weights)
  class(ans) <- "summary.evSyn"
  
  ans
}



## print.summary function
print.summary.evSyn <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  # added or equal approach
  type <- x$type
  # number of studies
  S <- x$n_studies
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  cat(paste(type, "Evidence Synthesis results:\n"), sep = "")
  
  indentation <- "    "  # Four spaces for indentation
  
  cat("\nStudy-specific results:\n")
  
  if (!is.null(x$GORICA_weight_m)) {
    cat("\n    GORICA weights:\n")  
    formatted_gw <- apply(x$GORICA_weight_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_gw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$GORICA_m)) {
    cat("\n    GORICA values:\n")  
    formatted_gv <- apply(x$GORICA_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_gv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$LL_m)) {
    cat("\n    Log-likelihood values:\n")  
    formatted_llv <- apply(x$LL_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_llv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$PT_m)) {
    cat("\n    Penalty term values:\n")  
    formatted_ptv <- apply(x$PT_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_ptv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  
  cat("\nCumulative results:\n")
  
  if (!is.null(x[["Cumulative_GORICA_weights"]])) {
    cat("\n    GORICA weights:\n")  
    formatted_cgw <- apply(x$Cumulative_GORICA_weights[1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cgw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["Cumulative_GORICA"]])) {
    cat("\n    GORICA values:\n")  
    formatted_cgv <- apply(x$Cumulative_GORICA[1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cgv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$LL_m)) {
    cat("\n    Log-likelihood values:\n")  
    formatted_cllv <- apply(x$Cumulative_LogLik[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cllv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$PT_m)) {
    cat("\n    Penalty term values:\n")  
    formatted_cptv <- apply(x$Cumulative_PT[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cptv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  
  cat("\nFinal results:\n")
  formatted_final <- apply(x$Final_Cumulative_results[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
  captured_output <- capture.output(print(formatted_final, row.names = TRUE, right = TRUE, quote = "FALSE"))
  adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
  cat(paste0(adjusted_output, "\n"), sep = "")
  
  
  cat("\nFinal ratios:\n")
  
  if (!is.null(x$Final_ratio_GORICA_weights)) {
    cat("\n    GORICA weights:\n")  
    formatted_frgw <- apply(x$Final_ratio_GORICA_weights, c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_frgw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$LL_m)) {
    cat("\n    Log-likelihood values:\n")  
    formatted_frllv <- apply(x$Final_ratio_Cumulative_LL, c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_frllv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  cat("\n")
  message(x$messages$mix_weights)
}



## plot function
plot.evSyn <- function(x, ...) {
  namesH <- colnames(x$GORICA_weight_m)
  NrHypos_incl <- ncol(x$GORICA_weight_m[,,drop = FALSE])
  S <- nrow(x$GORICA_weight_m[,,drop = FALSE])
  Name_studies <- as.factor(1:S)
  
  weight_m <- x$GORICA_weight_m
  CumulativeGoricaWeights <- x$Cumulative_GORICA_weights
  
  # Create data frame for per study weights
  if (all(is.na(weight_m))) {
    per_study_df <- NULL
    times <- 1
  } else {
    per_study_df <- data.frame(study = rep(Name_studies, NrHypos_incl),
                               weight = c(weight_m))
    per_study_df$weight_type <- "per study"
    times <- 2
  }
  # Create data frame for cumulative weights
  cumulative_df <- data.frame(study = rep(Name_studies, NrHypos_incl),
                              weight = c(CumulativeGoricaWeights[1:S, ]))
  cumulative_df$weight_type <- "cumulative"
  
  # Combine data frames
  plot_data <- rbind(per_study_df, cumulative_df)
  plot_data$variable <- rep(rep(namesH, each = S), times = times)
  
  # Create plot
  ggplot(plot_data, aes(x = .data[['study']],  
                        y = .data[['weight']], 
                        shape    = .data[['weight_type']], 
                        linetype = .data[['weight_type']], 
                        color    = .data[['variable']])) +
    geom_point(size = 3) +
    geom_line(data = plot_data[plot_data[['weight_type']] == "cumulative", ], 
              aes(group = .data[['variable']]), linewidth = 1) +
    scale_color_brewer(palette = "Dark2") +
    theme(
      plot.margin = unit(c(1,1,1,1), "cm"),
      legend.position = "bottom",
      legend.margin = margin(t = 10, r = 0, b = 3, l = 0),
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


