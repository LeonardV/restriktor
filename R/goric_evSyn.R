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

evSyn_est       <- function(object, ...) UseMethod("evSyn_est")
evSyn_LL        <- function(object, ...) UseMethod("evSyn_LL")
evSyn_ICweights <- function(object, ...) UseMethod("evSyn_ICweights")
evSyn_ICvalues  <- function(object, ...) UseMethod("evSyn_ICvalues")
evSyn_ICratios  <- function(object, ...) UseMethod("evSyn_ICratios")

# -------------------------------------------------------------------------
## est (list + vec) + cov (list + mat)
## LL (list + vec) + PT (list + vec)
## IC weights (list + vec) rowSums = 1
## IC values (list + vec)

## object = est, VCOV = cov
## object = LL, PT = PT
## object = IC weights + rowsum check
## object = IC values
## object = Ratio IC weights
# -------------------------------------------------------------------------

evSyn <- function(object, ...) {
  
  arguments <- list(...)
  
  VCOV <- arguments$VCOV
  PT   <- arguments$PT
  
  isGoric <- vapply(object, function(x) inherits(x, "con_goric"), logical(1))
  
  if (all(isGoric)) {
    return(evSyn_gorica.list(object, ...))
  }
  
  if (!is.list(object) || !any(vapply(object, is.numeric, logical(1)))) {
    stop("Restriktor ERROR: object must be a list of numeric vectors.")
  }
  
  checkListContent <- function(lst, fun, msg) {
    if (!is.null(lst) && (!is.list(lst) || !any(vapply(lst, fun, logical(1))))) {
      stop("Restriktor ERROR: ", msg)
    }
  }
  
  checkListContent(lst = VCOV, fun = is.matrix, msg = "VCOV must be a list of matrices.")
  checkListContent(lst = PT, fun = is.numeric, msg = "PT must be a list of numeric vectors.")
  
  if (!is.null(VCOV) && !is.null(PT)) {
    stop("Restriktor ERROR: both VCOV and PT are found, which confuses me.")
  }
  
  # if they are weights, the sum of all vectors must be 1.
  obj_isICweights <- all(abs(vapply(object, sum, numeric(1)) - 1) <= sqrt(.Machine$double.eps))
  # Check if they are IC ratios: each vector should end with 1.
  obj_isICratios <- all(vapply(object, function(x) tail(x, n = 1) == 1, logical(1)))
  
  if (!is.null(VCOV)) {
    return(evSyn_est(object, ...))
  } else if (!is.null(PT)) {
    return(evSyn_LL(object, ...))
  } else if (obj_isICweights) {
    if (!is.null(arguments$type) && arguments$type == "equal") {
      message("The added-evidence method is the only available approach when weights are included in the input.")
    }
    return(evSyn_ICweights(object, ...))
  } else if (obj_isICratios) {
    if (!is.null(arguments$type) && arguments$type == "equal") {
      message("The added-evidence method is the only available approach when weights are included in the input.")
    }
    return(evSyn_ICratios(object, ...))
  } else {
    if (!is.null(arguments$type) && arguments$type == "equal") {
      message("The added-evidence method is the only available approach when the input consists of GORIC(A) values.")
    }
    return(evSyn_ICvalues(object, ...))
  }
  
  stop("Restriktor Error: I don't know how to handle the input.")
}


# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on the (standardized) parameter estimates and 
# the covariance matrix
evSyn_est.list <- function(object, ..., VCOV = list(), hypotheses = list(),
                           type = c("added", "equal", "average"), 
                           comparison = c("unconstrained", "complement", "none"),
                           hypo_names = c()) {
  
  if (missing(comparison)) 
    comparison <- "unconstrained"
  comparison <- match.arg(comparison)

  if (missing(type)) 
    type <- "added"
  type <- match.arg(type)

  # number of primary studies
  S <- length(object)
  V <- length(VCOV)
  
  if (S != V) {
    stop("Restriktor ERROR: the number of items in the object list (i.e., number of (standardized) estimates) must equal the number of items in the VCOV list.", call. = FALSE)
  }
  
  # Ensure hypotheses are nested
  if (!all(vapply(hypotheses, is.list, logical(1)))) {
    hypotheses <- rep(list(hypotheses), S)
  }
  
  # check if VCOV and hypotheses are both a non-empty list
  if ( !is.list(VCOV) && length(VCOV) == 0 ) {
    stop("Restriktor ERROR: VCOV must be a list of covariance matrices of the (standardized) parameter estimates of interest.", call. = FALSE)  
  } 

  if ( !is.list(hypotheses) && length(hypotheses) == 0 ) {
    stop("Restriktor ERROR: hypotheses must be a list.", call. = FALSE)  
  } 
  
  
  # check if the matrices are all symmetrical
  VCOV_isSym <- vapply(VCOV, isSymmetric, logical(1), check.attributes = FALSE)
  if (!all(VCOV_isSym)) {
    stop(sprintf("Restriktor ERROR: the %sth covariance matrix in VCOV is not symmetric.", which(!VCOV_isSym)), call. = FALSE)  
  }

  # number of hypotheses must be equal for each studie. In each study a set of 
  # shared theories (i.e., hypotheses) are compared.
  len_H <- vapply(hypotheses, length, integer(1))
  NrHypos <- unique(len_H)
  
  if (length(unique(len_H)) > 1) {
    stop("Restriktor ERROR: The number of hypotheses must be consistent across all studies.", call. = FALSE)
  }
  
  if (length(object) != length(len_H)) {
    stop("Restriktor ERROR: The number of hypothesis sets ", "(",length(len_H), ")", " does not match the number of studies (", S, ")",".", call. = FALSE)
  }

  complement_check <- all(len_H == 1)
  if (comparison == "complement") {
    #if ((sameHypo && !comp_check_same) | (!sameHypo && !comp_check_diff)) {
     if (!complement_check) {  
      warning("Restriktor WARNING: Only one order-restricted hypothesis is allowed (for now) when comparison = 'complement'.",
              " Setting comparison to 'unconstrained' instead.", call. = FALSE)
      comparison <- "unconstrained"
    }
  }

  if (is.null(hypo_names)) {
    list_hypo_names <- lapply(hypotheses, names)
    # each study must have the same hypotheses namen
    element_hypo_names <- list_hypo_names[[1]]
    # check if it is similar is the other lists, but also if the same order is used. 
    # else try to fix the order if the same names are used, but in a different order.
    check_name_in_list <- vapply(list_hypo_names, function(x) all(element_hypo_names == x), logical(1))
    
    # this should trigger the WARNING: list(list(Ha = H11), list(Hb = H21))
    if (!is.null(element_hypo_names) && !all(unlist(check_name_in_list))) {
      l_hnames <- vapply(list_hypo_names, function(x) paste0("(", x, ")", collapse = ', '), character(1))
      warning(sprintf("Restriktor WARNING: The hypothesis names %s in each hypothesis set must be identical and in the same order. Renaming hypotheses to H1, H2, ... instead.", paste(l_hnames, collapse = ' and ')))
      element_hypo_names <- NULL
      # just to be sure that all names are removed
      #hypotheses <- lapply(hypotheses, function(x) { names(x) <- NULL; return(x) })
    }
  } else {
    element_hypo_names <- hypo_names
  }
 
  NrHypos_incl <- NrHypos + 1
  if (comparison == "none") {
    NrHypos_incl <- NrHypos
  }

  LL_weights_m <- LL_m <- PT <- GORICA_weight_m <- GORICA_m <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  rownames(LL_weights_m) <- rownames(GORICA_weight_m) <- rownames(GORICA_m) <- rownames(LL_m) <- rownames(PT) <- paste0("Study ", 1:S)
  
  CumulativeLLWeights <- CumulativeGoricaWeights <- CumulativeGorica <- matrix(NA, nrow = S+1, ncol = NrHypos_incl)

  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Studies 1"
  rownames(CumulativeLLWeights) <- rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
  
  if (NrHypos == 1 && comparison == "complement") {
    if (!is.null(element_hypo_names)) {
      element_hypo_names <- c(element_hypo_names, "Complement")
      hnames_idx <- element_hypo_names != ""
    } else {
      element_hypo_names <- vector("character", 2L)
      hnames_idx <- element_hypo_names != ""
    }
    hnames <- c("H1", "Complement")
    hnames_idx <- element_hypo_names != ""
    element_hypo_names[!hnames_idx] <- hnames[!hnames_idx]
    hnames <- element_hypo_names
    
    hypotheses <- lapply(hypotheses, function(h) {
      names(h)[1:(length(hnames) - 1L)] <- hnames[-max(length(hnames))]  
      return(h)
    })
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = 1)
  } else if (comparison == "none") {
    if (!is.null(element_hypo_names)) {
      hnames_idx <- element_hypo_names != ""
    } else {
      element_hypo_names <- vector("character", unique(len_H))
      hnames_idx <- element_hypo_names != ""
    }
    hnames <- c(paste0("H", 1:NrHypos))
    element_hypo_names[!hnames_idx] <- hnames[!hnames_idx]
    hnames <- element_hypo_names
    
    hypotheses <- lapply(hypotheses, function(h) {
      names(h)[seq_len(length(hnames))] <- hnames 
      return(h)
    })
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  } else {
    if (!is.null(element_hypo_names)) {
      element_hypo_names <- c(element_hypo_names, "Unconstrained")
      hnames_idx <- element_hypo_names != ""
    } else {
      element_hypo_names <- vector("character", unique(len_H) + 1L)
      hnames_idx <- element_hypo_names != ""
    }
    hnames <- c(paste0("H", 1:NrHypos), "Unconstrained")
    hnames_idx <- element_hypo_names != ""
    element_hypo_names[!hnames_idx] <- hnames[!hnames_idx]
    hnames <- element_hypo_names

    hypotheses <- lapply(hypotheses, function(h) {
      names(h)[1:(length(hnames) - 1L)] <- hnames[-max(length(hnames))]
      return(h)
    })
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  }
  
  rownames(ratio.weight_mu) <- paste0("Study ", 1:S)
  colnames(GORICA_m) <- colnames(GORICA_weight_m) <- colnames(LL_m) <- colnames(LL_weights_m) <- colnames(PT) <- hnames
  colnames(CumulativeLLWeights) <- colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- hnames

  for (s in 1:S) {
    res_goric <- goric(object[[s]], VCOV = VCOV[[s]],
                       hypotheses = hypotheses[[s]],
                       type = 'gorica', comparison = comparison,
                       ...)

    if (comparison == "unconstrained") {
      ratio.weight_mu[s, ] <- res_goric$ratio.gw[, NrHypos_incl]
    } else if (comparison == "complement") {
      ratio.weight_mu[s, ] <- res_goric$ratio.gw[1, NrHypos_incl]
    } 
    
    LL_m[s, ] <- res_goric$result$loglik
    LL_weights_m[s, ] <- res_goric$result$loglik.weights
    GORICA_m[s, ] <- res_goric$result$gorica
    GORICA_weight_m[s, ] <- res_goric$result$gorica.weights
    PT[s, ] <- res_goric$result$penalty
  }
  
  sumPT <- sumLL <- 0
  if (type == "added") { 
    # added-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else if (type == "equal") { 
    # equal-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else {
    # average-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL/s + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } 
  
  # cumulative log_likelihood weights
  sumLL <- 0
  for (l in 1:S) {
    sumLL <- sumLL + LL_m[l, ]
    CumulativeLL <- -2 * sumLL
    minLL <- min(CumulativeLL)
    CumulativeLLWeights[l, ] <- exp(-0.5*(CumulativeLL-minLL)) / sum(exp(-0.5*(CumulativeLL-minLL)))
  }
  
  # cumulative log-likilihood values
  Cumulative_LL <- apply(LL_m, 2, cumsum)
  Cumulative_LL <- matrix(Cumulative_LL, nrow = nrow(LL_m), 
                          dimnames = list(sequence, colnames(LL_m)))
  # final cumulative log-likilihood value
  Cumulative_LL_final <- -2*Cumulative_LL[S, , drop = FALSE]
  rownames(Cumulative_LL_final) <- "Final"
  minLL <- min(Cumulative_LL_final)
  # cumulative log-likihood weights
  Final.LL.weights <- exp(-0.5*(Cumulative_LL_final-minLL)) / sum(exp(-0.5*(Cumulative_LL_final-minLL)))
  Final.LL.weights <- Final.LL.weights[,, drop = TRUE]
  Final.ratio.LL.weights <- Final.LL.weights %*% t(1/Final.LL.weights)
  rownames(Final.ratio.LL.weights) <- hnames
  
  # add final row
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  CumulativeLLWeights[(S+1), ] <- CumulativeLLWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  rownames(Final.ratio.GORICA.weights) <- hnames
  
  # Output
  if (NrHypos == 1 && comparison == "complement") {
    colnames(ratio.weight_mu) <- c(paste0(hnames[1], " vs. ", "Complement"))
    colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
  } else if (comparison == "none") {
    #colnames(ratio.weight_mu) <- c(paste0(hnames[1], " vs. ", "Complement"))
    colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
  } else { 
    # unconstrained
    colnames(ratio.weight_mu) <- c(paste0(colnames(CumulativeGorica), " vs. ", "Unconstrained"))
    colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
  }
  
  out <- list(type = type,
              hypotheses = hypotheses,
              PT_m = PT,
              GORICA_weight_m = GORICA_weight_m, 
              LL_weights_m = LL_weights_m,
              GORICA_m = GORICA_m, 
              LL_m = LL_m, 
              Cumulative_GORICA = CumulativeGorica, 
              Cumulative_LL = Cumulative_LL,
              Cumulative_GORICA_weights = CumulativeGoricaWeights,
              Cumulative_LL_weights = CumulativeLLWeights,                  
              ratio_GORICA_weight_mu = ratio.weight_mu, 
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights,
              Final_ratio_LL_weights = Final.ratio.LL.weights
              )
  
  class(out) <- c("evSyn.est", "evSyn")
  
  return(out)
}


# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on log likelihood and penalty values
evSyn_LL.list <- function(object, ..., PT = list(), type = c("added", "equal", "average"),
                          hypo_names = c()) {
  
  if (missing(type)) 
    type <- "added"
  type <- match.arg(type)
  
  # check if PT is a non-empty list
  if ( !is.list(PT) && length(PT) == 0 ) {
    stop("Restriktor ERROR: PT must be a list of penalty weights.", call. = FALSE)  
  } 
  
  
  LL_m <- object
  S <- length(LL_m)
  NrHypos <- length(LL_m[[1]]) - 1
  if (is.null(hypo_names)) {
    hnames <- paste0("H", 1:(NrHypos + 1))
  } else {
    hnames <- hypo_names
  }
  
  LL_weights_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  CumulativeLL <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeLLWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  
  GORICA_weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  
  LL_m <- do.call(rbind, LL_m)
  PT <- do.call(rbind, PT)
  
  colnames(CumulativeLL) <- colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(CumulativeLLWeights) <- hnames
  colnames(LL_m) <- colnames(LL_weights_m) <- colnames(PT) <- colnames(GORICA_weight_m) <- hnames
  
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Studies 1"
  
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
  rownames(CumulativeLL) <- rownames(CumulativeLLWeights) <- c(sequence, "Final")
  rownames(LL_weights_m) <- rownames(LL_m) <- rownames(PT) <- rownames(GORICA_weight_m) <- paste0("Study ", 1:S)
  
  # cumulative log-likilihood values
  Cumulative_LL <- apply(LL_m, 2, cumsum)
  Cumulative_LL <- matrix(Cumulative_LL, nrow = nrow(LL_m), 
                          dimnames = list(sequence, colnames(LL_m)))
  # final cumulative log-likilihood value
  Cumulative_LL_final <- -2*Cumulative_LL[S, , drop = FALSE]
  minLL <- min(Cumulative_LL_final)
  # cumulative log-likihood weights
  Final.LL.weights <- exp(-0.5*(Cumulative_LL_final-minLL)) / sum(exp(-0.5*(Cumulative_LL_final-minLL)))
  Final.LL.weights <- Final.LL.weights[,, drop = TRUE]
  Final.ratio.LL.weights <- Final.LL.weights %*% t(1/Final.LL.weights) 
  
  sumLL <- 0
  for (l in 1:S) {
    LL <- -2*LL_m[l, ]
    delta_LL <- LL - min(LL)
    LL_weights_m[l, ] <- exp(-0.5 * delta_LL) / sum(exp(-0.5 * delta_LL))
    
    sumLL <- sumLL + LL_m[l, ]
    CumulativeLL <- -2 * sumLL
    minLL <- min(CumulativeLL)
    CumulativeLLWeights[l, ] <- exp(-0.5*(CumulativeLL-minLL)) / sum(exp(-0.5*(CumulativeLL-minLL)))
  }

  sumPT <- sumLL <- 0
  IC <- -2 * LL_m + 2 * PT
  if (type == "added") { 
    # added-ev approach
    for(s in 1:S) {
      minIC <- min(IC[s, ])
      GORICA_weight_m[s, ] <- exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC)))
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- -2 * sumLL + 2 * sumPT
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else if (type == "equal") { 
    # equal-ev approach
    for (s in 1:S) {
      minIC <- min(IC[s, ])
      GORICA_weight_m[s, ] <- exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC)))
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- -2 * sumLL + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else {
    # average-ev approach
    for (s in 1:S) {
      minIC <- min(IC[s, ])
      GORICA_weight_m[s, ] <- exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC)))
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- -2 * sumLL/s + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
    
  }

  # fill in the final row  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]

  CumulativeLLWeights[(S+1), ] <- CumulativeLLWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)

  rownames(Final.ratio.LL.weights) <- rownames(Final.ratio.GORICA.weights) <- hnames
  colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hnames)
  
  out <- list(type = type,
              PT_m = PT, 
              GORICA_weight_m = GORICA_weight_m,
              LL_weights_m = LL_weights_m,
              GORICA_m = IC, 
              LL_m = LL_m, 
              Cumulative_GORICA_weights = CumulativeGoricaWeights,
              Cumulative_LL_weights = CumulativeLLWeights,
              Cumulative_GORICA = CumulativeGorica, 
              Cumulative_LL = Cumulative_LL,
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights,
              Final_ratio_LL_weights = Final.ratio.LL.weights
              )
  
  class(out) <- c("evSyn.LL", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA values
evSyn_ICvalues.list <- function(object, ..., hypo_names = c()) {
  
  IC <- object
  S  <- length(IC)
  NrHypos <- length(IC[[1]]) - 1
  GORICA_weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  
  if (is.null(hypo_names)) {
    hnames <- paste0("H", 1:(NrHypos+1))
  } else {
    hnames <- hypo_names
  }
  
  IC <- do.call(rbind, IC)
  
  GORICA_weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(IC) <- colnames(GORICA_weight_m) <- hnames
  
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Studies 1"
  
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
  rownames(IC) <- rownames(GORICA_weight_m) <- paste0("Study ", 1:S)
  
  sumIC <- 0
  for (s in 1:S) {
    minIC <- min(IC[s, ])
    GORICA_weight_m[s, ] <- exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC)))
  
    sumIC <- sumIC + IC[s, ]
    CumulativeGorica[s, ] <- sumIC
    minGoric <- min(CumulativeGorica[s, ])
    CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
      sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
  }

  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  
  rownames(Final.ratio.GORICA.weights) <- hnames
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hnames)
  
  out <- list(type              = "added",
              GORICA_m          = IC, 
              GORICA_weight_m   = GORICA_weight_m,
              Cumulative_GORICA = CumulativeGorica, 
              Cumulative_GORICA_weights  = CumulativeGoricaWeights,
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn.ICvalues", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA weights or 
# (Bayesian) posterior model probabilities
evSyn_ICweights.list <- function(object, ..., priorWeights = NULL, hypo_names = c()) {
  
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
  
  for (s in 2:S) {
    CumulativeWeights[s, ] <- CumulativeWeights[(s-1), ] * Weights[s, ] / 
      sum( CumulativeWeights[(s-1), ] * Weights[s, ] )
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


# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on the ratio of AIC or ORIC or GORIC or GORICA 
# weights or (Bayesian) posterior model probabilities
evSyn_ICratios.list <- function(object, ..., priorWeights = NULL, hypo_names = c()) {
  
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
  
  for (s in 2:S) {
    CumulativeWeights[s, ] <- CumulativeWeights[(s-1), ] * Weights[s, ] / 
      sum( CumulativeWeights[(s-1), ] * Weights[s, ] )
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
  
  class(out) <- c("evSyn.ICratios", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# list with goric objects
evSyn_gorica.list <- function(object, ..., type = c("added", "equal", "average"), 
                              hypo_names = c()) {
  
  # Check if all objects are of type "con_goric"
  if (!all(vapply(object, function(x) inherits(x, "con_goric"), logical(1)))) {
    stop("Restriktor ERROR: the object must be a list with fitted objects from the goric() function", 
         call. = FALSE)
  }
  
  # Create a list for the evSyn_LL.list function
  conList <- list(
    object = lapply(object, function(x) x$result$loglik),
    PT = lapply(object, function(x) x$result$penalty),
    type = type,
    hypo_names = hypo_names
  )
  
  # Call the evSyn_LL.list function and return the result
  result <- do.call(evSyn_LL.list, conList)
  
  return(result)
}
