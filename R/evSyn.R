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
#1. de ratio matrix met de LL-weights toevoegen. 

#2. hypo namen willen meegeven, kan nl ook lastig zijn voor gebruiker (adhv vb in Rmd hierboven genoemd):
  
  # results_Set1 <- evSyn(object = estimates, VCOV = covmats,
  #                       hypotheses = list(H1.1, H1.2),
  #                       comparison = "unconstrained",
  #                       hypo_names = c("H1.1", "H1.2", "Hu")) 


#3. meerdere fit objecten in evSyn toevoegen? Of wellicht meerdere goric() objecten?

  
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
  
  # list(list(H11 = H11), list(H21 = H21)), 
  nested_hypo <- all(sapply(hypotheses, is.list))
  # if not nested, list(H0, Hpos, Hneg = Hneg), make it nested
  if (!nested_hypo) {
    hypotheses <- list(hypotheses)
    hypotheses <- rep(hypotheses, S)
  }
  
  # check if VCOV and hypotheses are both a non-empty list
  if ( !is.list(VCOV) & length(VCOV) == 0 ) {
    stop("Restriktor ERROR: VCOV must be a list of covariance matrices of the (standardized) parameter estimates of interest.", call. = FALSE)  
  } 

  if ( !is.list(hypotheses) & length(hypotheses) == 0 ) {
    stop("Restriktor ERROR: hypotheses must be a list.", call. = FALSE)  
  } 
  
  VCOV[[1]]
  
  # check if the matrices are all symmetrical
  VCOV_isSym <- sapply(VCOV, isSymmetric, check.attributes = FALSE)
  if (!all(VCOV_isSym)) {
    stop(sprintf("Restriktor ERROR: the %sth covariance matrix in VCOV is not symmetric.", which(!VCOV_isSym)), call. = FALSE)  
  }

  # number of hypotheses must be equal for each studie. In each study a set of 
  # shared theories (i.e., hypotheses) are compared.
  len_H <- sapply(hypotheses, length) 
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

  list_hypo_names <- lapply(hypotheses, names)
  # each study must have the same hypotheses namen
  #element_hypo_names <- l_hnames[[1]]
  element_hypo_names <- list_hypo_names[[1]]
  # check if it is similar is the other lists, but also if the same order is used. 
  # else try to fix the order if the same names are used, but in a different order.
  check_name_in_list <- lapply(list_hypo_names, function(x) all(element_hypo_names == x))
  
  # this should trigger the WARNING: list(list(Ha = H11), list(Hb = H21))
  if (!is.null(element_hypo_names) && !all(unlist(check_name_in_list))) {
    l_hnames <- lapply(list_hypo_names, function(x) paste0("(",x,")", collapse = ', '))
    l_hnames <- paste(l_hnames, collapse = ' and ')
    warning(paste0("Restriktor WARNING: The hypothesis names ", l_hnames, " in each hypothesis set must be identical and in the same order. ",
                   "For example, hypotheses = list(list(Ha = H11, Hb = H12), list(Ha = H21, Hb = H22)). ",
                   "Hence, the hypotheses have been renamed to ", paste0("H", 1:NrHypos, collapse = " and "), " instead."),
            call. = FALSE)
    
    element_hypo_names <- NULL
    # just to be sure that all names are removed
    #hypotheses <- lapply(hypotheses, function(x) { names(x) <- NULL; return(x) })
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
  
  if (NrHypos == 1 & comparison == "complement") {
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
      names(h)[1:length(hnames)] <- hnames 
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
    # if (sameHypo) {
    #   res_goric <- goric(object[[s]], VCOV = VCOV[[s]], 
    #                      hypotheses = as.list(unlist(hypotheses)),
    #                      type = 'gorica', comparison = comparison)

    res_goric <- goric(object[[s]], VCOV = VCOV[[s]],
                       hypotheses = hypotheses[[s]],
                       type = 'gorica', comparison = comparison)

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
  } else { 
    # equal-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT/s
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
  if (NrHypos == 1 & comparison == "complement") {
    colnames(ratio.weight_mu) <- c(paste0(hnames[1], " vs. ", "Complement"))
    colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
    
    # out <- list(type = type, 
    #             hypotheses = hypotheses,
    #             PT_m = PT,
    #             GORICA_m = GORICA_m, 
    #             GORICA_weight_m = GORICA_weight_m, 
    #             ratio_GORICA_weight_mc = ratio.weight_mu, 
    #             LL_m = LL_m, 
    #             LL_weights_m = LL_weights_m,
    #             Cumulative_LL = Cumulative_LL,
    #             Final_ratio_LL_weights = Final.ratio.LL.weights,
    #             Cumulative_GORICA = CumulativeGorica, 
    #             Cumulative_GORICA_weights = CumulativeGoricaWeights,
    #             Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  } else if (comparison == "none") {
    colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
    
    # out <- list(type = type,
    #             hypotheses = hypotheses,
    #             PT_m = PT,
    #             GORICA_weight_m = GORICA_weight_m, 
    #             LL_weights_m = LL_weights_m,
    #             GORICA_m = GORICA_m, 
    #             LL_m = LL_m, 
    #             Cumulative_GORICA = CumulativeGorica, 
    #             Cumulative_LL = Cumulative_LL,
    #             
    #             Cumulative_GORICA_weights = CumulativeGoricaWeights,
    #             Cumulative_LL_weights = CumulativeLLWeights,
    #             
    #             Final_ratio_GORICA_weights = Final.ratio.GORICA.weights,
    #             Final_ratio_LL_weights = Final.ratio.LL.weights
    #             )
  } else { 
    # unconstrained
    colnames(ratio.weight_mu) <- c(paste0(colnames(CumulativeGorica), " vs. ", "Unconstrained"))
    colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
    
    # out <- list(type = type,
    #             hypotheses = hypotheses,
    #             PT_m = PT,
    #             GORICA_weight_m = GORICA_weight_m, 
    #             LL_weights_m = LL_weights_m,
    #             GORICA_m = GORICA_m, 
    #             LL_m = LL_m, 
    #             Cumulative_GORICA = CumulativeGorica, 
    #             Cumulative_LL = Cumulative_LL,
    #             
    #             Cumulative_GORICA_weights = CumulativeGoricaWeights,
    #             Cumulative_LL_weights = CumulativeLLWeights,                    
    #             
    #             ratio_GORICA_weight_mu = ratio.weight_mu, 
    #             Final_ratio_GORICA_weights = Final.ratio.GORICA.weights,
    #             Final_ratio_LL_weights = Final.ratio.LL.weights
    #             )
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
evSyn_LL.list <- function(object, ..., PT = list(), type = c("equal", "added"),
                          hypo_names = NULL) {
  
  if (missing(type)) 
    type <- "added"
  type <- match.arg(type)
  
  # check if PT is a non-empty list
  if ( !is.list(PT) & length(PT) == 0 ) {
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
  } else { 
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
  }

  # fill in the final row  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]

  #CumulativeLL[(S+1), ] <- CumulativeLL[S, ]
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
evSyn_ICvalues.list <- function(object, ..., hypo_names = NULL) {
  
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
    LL_weights_m = x$LL_weights_m,
    Cumulative_LL_weights = x$Cumulative_LL_weights,
    Cumulative_LL = x$Cumulative_LL,
    LL_m = x$LL_m,
    PT_m = x$PT_m,
    Cumulative_GORICA_weights = x$Cumulative_GORICA_weights,
    Cumulative_GORICA = x$Cumulative_GORICA
  )

  # if (!is.null(x$LL_m)) {
  #   sequence    <- paste0("Studies 1-", 1:ans$n_studies)
  #   sequence[1] <- "Studies 1"
  #   Cumulative_LL <- apply(x$LL_m, 2, cumsum)
  #   Cumulative_LL <- matrix(Cumulative_LL, nrow = nrow(x$LL_m),
  #                           dimnames = list(sequence, colnames(x$LL_m)))
  #   ans$Cumulative_LL <- Cumulative_LL
  # }
  
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

  if (!is.null(x[["Cumulative_LL_weights"]])) {
    fcgw <- t(x[["Cumulative_LL_weights"]]["Final", ])
    rownames(fcgw) <- "Log-likelihood weights"
    colnames(fcgw)  <- colnames(x[["Cumulative_LL_weights"]])
    final <- rbind(final, fcgw)
  }
  
  if (!is.null(x[["Cumulative_LL"]])) {
    fcllv <- t(x[["Cumulative_LL"]][ans$n_studies, ])
    rownames(fcllv) <- "Log-likelihood values"
    final <- rbind(final, fcllv)
  }
  
  if (!is.null(ans$Cumulative_PT)) { 
    fcptv <- t(ans$Cumulative_PT[ans$n_studies, ])
    rownames(fcptv) <- "Penalty term values"
    final <- rbind(final, fcptv)
  }
  
  ans$Final_Cumulative_results <- final
  
  # if (!is.null(x$LL_m)) {
  #   Final_ratio_Cumulative_LL <- ans$Cumulative_LL[ans$n_studies, ] %*% t(1/ans$Cumulative_LL[ans$n_studies, ])
  #   rownames(Final_ratio_Cumulative_LL) <- colnames(ans$Cumulative_LL)
  #   colnames(Final_ratio_Cumulative_LL) <- paste0("vs. ", colnames(ans$Cumulative_LL))
  #   ans$Final_ratio_Cumulative_LL <- Final_ratio_Cumulative_LL  
  # }
  
  ans$Ratio_GORICA_weight_mu <- x$ratio_GORICA_weight_mu
  ans$Ratio_GORICA_weight_mc <- x$ratio_GORICA_weight_mc
  ans$Final_ratio_GORICA_weights <- x$Final_ratio_GORICA_weights
  ans$Final_ratio_LL_weights <- x$Final_ratio_LL_weights
  
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

  if (!is.null(x$LL_weights_m)) {
    cat("\n    Log-likelihood weights:\n")  
    formatted_llw <- apply(x$LL_weights_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_llw, row.names = TRUE, right = TRUE, quote = "FALSE"))
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
  
  if (!is.null(x[["Cumulative_LL_weights"]])) {
    cat("\n    Log-likelihood weights:\n")  
    formatted_cllw <- apply(x$Cumulative_LL_weights[1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cllw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$LL_m)) {
    cat("\n    Log-likelihood values:\n")  
    formatted_cllv <- apply(x$Cumulative_LL[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
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
  
  if (!is.null(x$Final_ratio_LL_weights)) {
    cat("\n    Log-likelihood weights:\n")  
    formatted_frllw <- apply(x$Final_ratio_LL_weights, c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_frllw, row.names = TRUE, right = TRUE, quote = "FALSE"))
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
  
  GORICA_weight_m <- x$GORICA_weight_m
  CumulativeGoricaWeights <- x$Cumulative_GORICA_weights
  
  # Create data frame for per study weights
  if (all(is.na(GORICA_weight_m))) {
    per_study_df <- NULL
    times <- 1
  } else {
    per_study_df <- data.frame(study = rep(Name_studies, NrHypos_incl),
                               weight = c(GORICA_weight_m))
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


