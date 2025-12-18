## input options
# 1. est + vcov
# 2. LL + PT
# 3. IC values (AIC, ORIC, GORIC, GORICA)
# 4. IC weights or (Bayesian) posterior model probs. 
# 5. Output from gorica function
# 6. Output from escalc (metafor package)

# In case of an equal-evidence approach, aggregating evidence from, say, 5 studies 
# with n=100 observations is the same as obtaining evidence from 1 study 
# (as if it was possible) with n=500 observations (like meta-analysis does).
# In the added-evidence approach, the aggregated evidence from, says, 5 studies 
# is stronger than as if the data were combined (as if that was possible).

# evSyn_est       <- function(object, ...) UseMethod("evSyn_est")
# evSyn_LL        <- function(object, ...) UseMethod("evSyn_LL")
# evSyn_ICweights <- function(object, ...) UseMethod("evSyn_ICweights")
# evSyn_ICvalues  <- function(object, ...) UseMethod("evSyn_ICvalues")
# evSyn_ICratios  <- function(object, ...) UseMethod("evSyn_ICratios")
# evSyn_escalc    <- function(object, ...) UseMethod("evSyn_escalc")
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

evSyn <- function(object, input_type = NULL, ...) {
  
  arguments <- list(...)
  
  VCOV <- arguments$VCOV
  PT   <- arguments$PT
  
  isGoric <- vapply(object, function(x) inherits(x, "con_goric"), logical(1))
  
  # Handle input_type explicitly if provided
  if (!is.null(input_type)) {
    input_type <- tolower(input_type)
    if (input_type == "est_vcov") {
      return(evSyn_est(object, ...))
    } else if (input_type == "ll_pt") {
      return(evSyn_LL(object, ...))
    } else if (input_type == "icweights") {
      return(evSyn_ICweights(object, ...))
    } else if (input_type == "icratios") {
      return(evSyn_ICratios(object, ...))
    } else if (input_type == "icvalues") {
      return(evSyn_ICvalues(object, ...))
    } else if (input_type %in% c("goric", "gorica", "goricc", "goricac")) {
      return(evSyn_gorica(object, ...))
    } else if (input_type == "escalc") {
      return(evSyn_escalc(object, ...))
    } else {
      stop(paste0("\nrestriktor ERROR: Unknown input_type ", sQuote(input_type), "."))
    }
  }
  
  if (all(isGoric)) {
    return(evSyn_gorica(object, ...))
  }
  
  if (inherits(object, "escalc")) {
    return(evSyn_escalc(object, ...))
  } 
  
  if (!is.list(object) || !any(vapply(object, is.numeric, logical(1)))) {
    stop("\nrestriktor ERROR: object must be a list of numeric vectors.", call. = FALSE)
    # TO DO
    # Laat evt volgende toe zodat namen hypos direct meegenomen kan worden:
    #df <- data.frame(myGORICs)
    #goric.values <- lapply(as.list(1:dim(df)[1]), function(x) df[x[1],])
    # Met volgende nl geen namen maar wel numeric vector: goric.values <- as.list(data.frame(t(myGORICs)))
    # Dito als LL & PT input zijn
  }
  
  checkListContent <- function(lst, fun, msg) {
    if (!is.null(lst) && (!is.list(lst) || !any(vapply(lst, fun, logical(1))))) {
      stop("restriktor ERROR: ", msg, call. = FALSE)
    }
  }
  
  #checkListContent(lst = VCOV, fun = is.matrix, msg = "VCOV must be a list of matrices.")
  checkListContent(lst = PT, fun = is.numeric, msg = "PT must be a list of numeric vectors.")
  
  if (!is.null(VCOV) && !is.null(PT)) {
    stop("\nrestriktor ERROR: both VCOV and PT are found, which confuses me.", call. = FALSE)
    # TO DO zeg wat je wel graag wilt, dus iets als 'een van de volgende twee opties: ...'
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
    if (!is.null(arguments$type_ev) && arguments$type_ev == "equal") {
      messageAdded <- "\nrestriktor Message: When the input consists of weights, the equal-evidence approach is not applicable. The added-evidence approach is used instead."
      message(messageAdded)
    }
    return(evSyn_ICweights(object, ...))
  } else if (obj_isICratios) {
    if (!is.null(arguments$type_ev) && arguments$type_ev == "equal") {
      messageAdded <- "\nrestriktor Message: When the input consists of ratios of weights, the equal-evidence approach is not applicable. The added-evidence approach is used instead."
      message(messageAdded)
    }
    return(evSyn_ICratios(object, ...))
  } else { # ICvalues
    if (!is.null(arguments$type_ev) && arguments$type_ev == "equal") {
      messageAdded <- "\nrestriktor Message: When the input consists of IC values, the equal-evidence approach is not applicable. The added-evidence approach is used instead."
      message(messageAdded)
    }
    return(evSyn_ICvalues(object, ...))
    # TO DO werkt niet, ook niet als messageAdded = messageAdded
    # TO DO doe dit dan ook voor bovenstaande!
  }
  stop("\nrestriktor ERROR: I don't know how to handle the input.", call. = FALSE)
}


# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on the (standardized) parameter estimates and 
# the covariance matrix
evSyn_est <- function(object, ..., VCOV = list(), hypotheses = list(),
                      type_ev = c("added", "equal", "average"), 
                      comparison = c("unconstrained", "complement", "none"),
                      hypo_names = c(),
                      type = c("gorica", "goricac"),
                      order_studies = c("input_order", "ascending", "descending"),
                      study_names = c(),
                      study_sample_nobs = NULL
                      ) {
  
  if (missing(comparison)) {
    if (length(hypotheses) == 1) {
      comparison <- "complement"
    } else {
      comparison <- "unconstrained"
    }
  }
  comparison <- match.arg(comparison)

  if (missing(type)) { type <- "gorica" }
  #
  if (type == "goric") {
    message("\nrestriktor Message: Since the input is a list of estimates, the GORICA will be used, not the the GORIC.")
    type = "gorica"
  } else if (type == "goricc") {
    message("\nrestriktor Message: Since the input is a list of estimates, the GORICAC will be used, not the the GORICC.")
    type = "goricac"
  } 
  #
  if (type == "goricac" && is.null(study_sample_nobs)) {
    stop("\nrestriktor ERROR: To compute the GORICAC, the argument 'study_sample_nobs' is required. ",
         "Please provide a numeric vector with the sample sizes of all primary studies.",
         call. = FALSE)
  }
  type <- match.arg(type)
  
  if (missing(type_ev)) 
    type_ev <- "added"
  type_ev <- match.arg(type_ev)
  
  if (missing(order_studies)) 
    order_studies <- "input_order"
  order_studies <- match.arg(order_studies)

  # number of primary studies
  S <- length(object)
  V <- length(VCOV)
  
  if (S != V) {
    stop("\nrestriktor ERROR: The number of elements in the 'object' list (i.e., the number of (standardized) estimates) must match the number of elements in the 'VCOV' list.",
         call. = FALSE)
  }
  
  if (type == "goricac" && length(study_sample_nobs) == 1) {
    study_sample_nobs <- rep(study_sample_nobs, S)
    message("\nrestriktor Message: The argument 'study_sample_nobs' contains a single value; all primary studies are assumed to have the same sample size.")
  } else if (type == "goricac" && length(study_sample_nobs) != S) {
    stop("\nrestriktor ERROR: The argument 'study_sample_nobs' must be a numeric vector containing S = ", S, " values (one for each study).",
         call. = FALSE)
  }
  
  # Ensure hypotheses are nested
  if (!all(vapply(hypotheses, is.list, logical(1)))) {
    hypotheses <- rep(list(hypotheses), S)
  }
  
  # check if VCOV and hypotheses are both a non-empty list
  if ( !is.list(VCOV) && length(VCOV) == 0 ) {
    stop("\nrestriktor ERROR: The argument 'VCOV' must be a list of covariance matrices corresponding to the (standardized) parameter estimates of interest.",
         call. = FALSE)
  } 
  # TO DO check of square matrices in list, kan evt door nu fout in goric() te laten gebeuren, maar
  #       dan is het meegeven van study nr wel fijn!

  if ( !is.list(hypotheses) && length(hypotheses) == 0 ) {
    stop("\nrestriktor ERROR: hypotheses must be a list.", call. = FALSE)  
  } 
  
  
  # check if the matrices are all symmetrical
  #VCOV_isSym <- vapply(VCOV, isSymmetric, logical(1), check.attributes = FALSE)
  #if (!all(VCOV_isSym)) {
  #  stop(sprintf("\nrestriktor ERROR: the %sth covariance matrix in VCOV is not symmetric.", which(!VCOV_isSym)), call. = FALSE)  
  #}

  # number of hypotheses must be equal for each study. In each study a set of 
  # shared theories (i.e., hypotheses) are compared.
  len_H <- vapply(hypotheses, length, integer(1))
  NrHypos <- unique(len_H)
  
  if (length(unique(len_H)) > 1) {
    stop("\nrestriktor ERROR: The number of hypotheses must be identical across all studies.",
         call. = FALSE)
  }
  
  if (length(object) != length(len_H)) {
    stop("\nrestriktor ERROR: The number of hypothesis sets (", length(len_H), 
         ") does not match the number of studies (", S, ").",
         call. = FALSE)
  }

  complement_check <- all(len_H == 1)
  if (comparison == "complement") {
    #if ((sameHypo && !comp_check_same) | (!sameHypo && !comp_check_diff)) {
     if (!complement_check) {  
       warning("\nrestriktor WARNING: Only one order-restricted hypothesis is currently supported when comparison = 'complement'. ",
               "The comparison type has been set to 'unconstrained' instead.",
               call. = FALSE)
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
      warning(sprintf(
        "\nrestriktor WARNING: The hypothesis names (%s) within each hypothesis set must be identical and appear in the same order. ",
        paste(l_hnames, collapse = " and ")
      ),
      "Renaming hypotheses to 'H1', 'H2', ... instead.",
      call. = FALSE)
      element_hypo_names <- NULL
      # just to be sure that all names are removed
      #hypotheses <- lapply(hypotheses, function(x) { names(x) <- NULL; return(x) })
    }
  } else {
    # TO DO check length of hypo_names
    #       doe ws op meer plekken dan, nl in andere functies ook
    element_hypo_names <- hypo_names
  }
 
  NrHypos_incl <- NrHypos + 1
  if (comparison == "none") {
    NrHypos_incl <- NrHypos
  }

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
  
  
  LL_m <- LL_weights_m <- GORICA_m <- GORICA_weight_m <- PT <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  colnames(LL_m) <- colnames(LL_weights_m) <- colnames(GORICA_m) <- colnames(GORICA_weight_m) <- colnames(PT) <- hnames
  # rownames are set after determining the order of the studies
  #
  for (s in 1:S) {
    res_goric <- goric(object[[s]], VCOV = VCOV[[s]],
                       hypotheses = hypotheses[[s]],
                       type = type, comparison = comparison,
                       sample_nobs = study_sample_nobs[s],
                       ...)

    if (comparison == "unconstrained") {
      ratio.weight_mu[s, ] <- res_goric$ratio.gw[, NrHypos_incl]
    } else if (comparison == "complement") {
      ratio.weight_mu[s, ] <- res_goric$ratio.gw[1, NrHypos_incl]
    } 
    
    LL_m[s, ] <- res_goric$result$loglik
    LL_weights_m[s, ] <- res_goric$result$loglik.weights
    GORICA_m[s, ] <- res_goric$result[[type]] #res_goric$result$gorica
    GORICA_weight_m[s, ] <- res_goric$result[[paste0(type, ".weights")]]
    PT[s, ] <- res_goric$result$penalty
  }
  
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  #if (order_studies != "input_order"){
  # Order needs to be changed.
  if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    if (type_ev == "added") { 
      # added-evidence approach
      OverallGoric <- colSums(-LL_m) + colSums(PT)
      OverallPrefHypo <- which(OverallGoric == min(OverallGoric))
    } else if (type_ev == "equal") { 
      # equal-evidence approach
      OverallGoric <- colSums(-LL_m) + colMeans(PT)
      OverallPrefHypo <- which(OverallGoric == min(OverallGoric))
    } else if (type_ev == "average") { 
      # average-evidence approach
      OverallGoric <- colMeans(-LL_m) + colMeans(PT)
      OverallPrefHypo <- which(OverallGoric == min(OverallGoric))
    } #else {}
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    orderStudies <- order(GORICA_weight_m[,OverallPrefHypo], decreasing = decreasing)
    #
    LL_m <- LL_m[orderStudies,]
    LL_weights_m <- LL_weights_m[orderStudies,]
    GORICA_m <- GORICA_m[orderStudies,]
    GORICA_weight_m <- GORICA_weight_m[orderStudies,]
    PT <- PT[orderStudies,]
  }
  
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    if (length(study_names) != S) {
      stop("restriktor ERROR: Length of study_names must match the number of studies.", call. = FALSE)
    }
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(LL_m) <- rownames(LL_weights_m) <- rownames(GORICA_m) <- rownames(GORICA_weight_m) <- rownames(PT) <- study_names
  rownames(ratio.weight_mu) <- study_names
  
  CumulativeLLWeights <- CumulativeGoricaWeights <- CumulativeGorica <- matrix(NA, nrow = S+1, ncol = NrHypos_incl)
  colnames(CumulativeLLWeights) <- colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- hnames
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Study 1"
  rownames(CumulativeLLWeights) <- rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
  #
  sumPT <- sumLL <- 0
  if (type_ev == "added") { 
    # added-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else if (type_ev == "equal") { 
    # equal-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else if (type_ev == "average") { 
    # average-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL/s + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } # else {}
  
  # cumulative log_likelihood weights
  sumLL <- 0
  for (l in 1:S) {
    sumLL <- sumLL + LL_m[l, ]
    CumulativeLL <- -2 * sumLL
    minLL <- min(CumulativeLL)
    CumulativeLLWeights[l, ] <- exp(-0.5*(CumulativeLL-minLL)) / sum(exp(-0.5*(CumulativeLL-minLL)))
  }
  
  # cumulative log-likelihood values
  Cumulative_LL <- apply(LL_m, 2, cumsum)
  Cumulative_LL <- matrix(Cumulative_LL, nrow = nrow(LL_m), 
                          dimnames = list(sequence, colnames(LL_m)))
  # final cumulative log-likelihood value
  Cumulative_LL_final <- -2*Cumulative_LL[S, , drop = FALSE]
  rownames(Cumulative_LL_final) <- "Final"
  minLL <- min(Cumulative_LL_final)
  # cumulative log-likelihood weights
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
              type_ev = type_ev,
              hypotheses = hypotheses,
              order_studies = orderStudies,
              study_names = study_names,
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
  # TO DO welke volgorde en dan in alle functies zo gelijk mogelijk maken ook
  
  class(out) <- c("evSyn_est", "evSyn")
  
  return(out)
}


# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on log likelihood and penalty values
evSyn_LL <- function(object, ..., PT = list(), 
                     type_ev = c("added", "equal", "average"),
                     hypo_names = c(),
                     order_studies = c("input_order", "ascending", "descending"),
                     study_names = c()) {
  
  if (missing(type_ev)) 
    type_ev <- "added"
  type_ev <- match.arg(type_ev)
  
  if (missing(order_studies)) 
    order_studies <- "input_order"
  order_studies <- match.arg(order_studies)
  
  # check if PT is a non-empty list
  if ( !is.list(PT) && length(PT) == 0 ) {
    stop("\nrestriktor ERROR: PT must be a list of penalty weights.", call. = FALSE)  
  } 
  
  LL_m <- object
  S <- length(LL_m)
  NrHypos <- length(LL_m[[1]]) - 1
  if (is.null(hypo_names)) {
    hnames <- paste0("H", 1:(NrHypos + 1))
  } else {
    hnames <- hypo_names
  }
  
  LL_m <- do.call(rbind, LL_m)
  PT <- do.call(rbind, PT)
  IC <- -2 * LL_m + 2 * PT
  #
  GORICA_weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  for(s in 1:S) {
    minIC <- min(IC[s, ])
    GORICA_weight_m[s, ] <- exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC)))
  }
    
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  #if (order_studies != "input_order"){
  # Order needs to be changed.
  if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    if (type_ev == "added") { 
      # added-evidence approach
      OverallGoric <- colSums(-LL_m) + colSums(PT)
      OverallPrefHypo <- which(OverallGoric == min(OverallGoric))
    } else if (type_ev == "equal") { 
      # equal-evidence approach
      OverallGoric <- colSums(-LL_m) + colMeans(PT)
      OverallPrefHypo <- which(OverallGoric == min(OverallGoric))
    } else if (type_ev == "average") { 
      # average-evidence approach
      OverallGoric <- colMeans(-LL_m) + colMeans(PT)
      OverallPrefHypo <- which(OverallGoric == min(OverallGoric))
    } # else {}
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    orderStudies <- order(GORICA_weight_m[,OverallPrefHypo], decreasing = decreasing)
    #
    LL_m <- LL_m[orderStudies,]
    PT <- PT[orderStudies,]
    IC <- IC[orderStudies,]
    GORICA_weight_m <- GORICA_weight_m[orderStudies,]
  }
  
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    # TO DO check
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(LL_m) <- rownames(PT) <- rownames(IC) <- rownames(GORICA_weight_m) <- study_names
  # Set colnames
  colnames(LL_m) <- colnames(PT) <- colnames(IC) <- colnames(GORICA_weight_m) <- hnames
  
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Study 1"
  #
  # cumulative log-likelihood values
  Cumulative_LL <- apply(LL_m, 2, cumsum)
  Cumulative_LL <- matrix(Cumulative_LL, nrow = nrow(LL_m), 
                          dimnames = list(sequence, colnames(LL_m)))
  # final cumulative log-likelihood value
  Cumulative_LL_final <- -2*Cumulative_LL[S, , drop = FALSE]
  minLL <- min(Cumulative_LL_final)
  # cumulative log-likelihood weights
  Final.LL.weights <- exp(-0.5*(Cumulative_LL_final-minLL)) / sum(exp(-0.5*(Cumulative_LL_final-minLL)))
  Final.LL.weights <- Final.LL.weights[,, drop = TRUE]
  Final.ratio.LL.weights <- Final.LL.weights %*% t(1/Final.LL.weights) 
  
  sumLL <- 0
  LL_weights_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  CumulativeLLWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  rownames(LL_weights_m) <- study_names
  rownames(CumulativeLLWeights) <- c(sequence, "Final")
  colnames(LL_weights_m) <- colnames(CumulativeLLWeights) <- hnames
  for (l in 1:S) {
    LL <- -2*LL_m[l, ]
    delta_LL <- LL - min(LL)
    LL_weights_m[l, ] <- exp(-0.5 * delta_LL) / sum(exp(-0.5 * delta_LL))
    #
    sumLL <- sumLL + LL_m[l, ]
    CumulativeLL <- -2 * sumLL
    minLL <- min(CumulativeLL)
    CumulativeLLWeights[l, ] <- exp(-0.5*(CumulativeLL-minLL)) / sum(exp(-0.5*(CumulativeLL-minLL)))
  }
  
  sumPT <- sumLL <- 0
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- hnames
  if (type_ev == "added") { 
    # added-ev approach
    for(s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- -2 * sumLL + 2 * sumPT
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else if (type_ev == "equal") { 
    # equal-ev approach
    for (s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- -2 * sumLL + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else if (type_ev == "average") { 
    # average-ev approach
    for (s in 1:S) {
      sumLL <- sumLL + LL_m[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- -2 * sumLL/s + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
    
  } # else {}

  # fill in the final row  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]

  CumulativeLLWeights[(S+1), ] <- CumulativeLLWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)

  rownames(Final.ratio.LL.weights) <- rownames(Final.ratio.GORICA.weights) <- hnames
  colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hnames)
  
  out <- list(type_ev = type_ev,
              order_studies = orderStudies,
              study_names = study_names,
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

  class(out) <- c("evSyn_LL", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA values
evSyn_ICvalues <- function(object, ..., type_ev = c("added", "average"), 
                           hypo_names = c(),
                           order_studies = c("input_order", "ascending", "descending"),
                           study_names = c()) {
  
  if (missing(type_ev)) 
    type_ev <- "added"
  type_ev <- match.arg(type_ev)
  
  IC <- object
  S  <- length(IC)
  NrHypos <- length(IC[[1]]) - 1
  GORICA_weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  
  if (is.null(hypo_names)) {
    hnames <- paste0("H", 1:(NrHypos+1))
  } else {
    hnames <- hypo_names
  }
  
  if (missing(order_studies)) 
    order_studies <- "input_order"
  order_studies <- match.arg(order_studies)
  
  IC <- do.call(rbind, IC)
  #
  GORICA_weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  for(s in 1:S) {
    minIC <- min(IC[s, ])
    GORICA_weight_m[s, ] <- exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC)))
  }
  
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  #if (order_studies != "input_order"){
  # Order needs to be changed.
  if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    if (type_ev == "average") { 
      # average-evidence approach
      OverallGoric <- colMeans(IC)
      OverallPrefHypo <- which(OverallGoric == min(OverallGoric))
    } else {
      # type_ev == "added" (or when "equal", because then it is overruled to be "added")
      type_ev = "added"
      OverallGoric <- colSums(IC)
      OverallPrefHypo <- which(OverallGoric == min(OverallGoric))
    }
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    orderStudies <- order(GORICA_weight_m[,OverallPrefHypo], decreasing = decreasing)
    #
    IC <- IC[orderStudies,]
    GORICA_weight_m <- GORICA_weight_m[orderStudies,]
  }
  
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  #
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(IC) <- colnames(GORICA_weight_m) <- hnames
  #
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    # TO DO check
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(IC) <- rownames(GORICA_weight_m) <- study_names
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Study 1"
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
  #
  if (type_ev == "average") { 
    # average-ev approach
    sumIC <- 0
    for (s in 1:S) {
      sumIC <- sumIC + IC[s, ]
      CumulativeGorica[s, ] <- sumIC/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  } else {
    # type_ev == "added" (or when "equal", because then it is overruled to be "added")
    sumIC <- 0
    for (s in 1:S) {
      sumIC <- sumIC + IC[s, ]
      CumulativeGorica[s, ] <- sumIC
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / 
        sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
  }

  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  
  rownames(Final.ratio.GORICA.weights) <- hnames
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hnames)
  
  out <- list(type_ev           = type_ev,
              order_studies     = orderStudies,
              study_names       = study_names,
              GORICA_m          = IC, 
              GORICA_weight_m   = GORICA_weight_m,
              Cumulative_GORICA = CumulativeGorica, 
              Cumulative_GORICA_weights  = CumulativeGoricaWeights,
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  
  # if (!is.null(messageAdded)) {
  #   out <- append(out, messageAdded)
  #   # TO DO dit is ws nodig als message meegeven, doe ook bij andere varianten dan nog
  # }
  
  class(out) <- c("evSyn_ICvalues", "evSyn")
  
  return(out)
  
  }



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA weights or 
# (Bayesian) posterior model probabilities
evSyn_ICweights <- function(object, ..., priorWeights = NULL, hypo_names = c(),
                            order_studies = c("input_order", "ascending", "descending"),
                            study_names = c()) {
  
  Weights <- object
  S <- length(Weights)
  Weights <- do.call(rbind, Weights)
  NrHypos <- ncol(Weights)
  
  if (is.null(priorWeights)) {
    priorWeights <- rep(1/(NrHypos), (NrHypos))
  }
  # To make it sum to 1 (if it not already did)
  # TO DO give message that this is done
  priorWeights <- priorWeights / sum(priorWeights) 
  
  if (is.null(hypo_names)) {
    hypo_names <- paste0("H", 1:(NrHypos)) 
  }
  
  if (missing(order_studies)) 
    order_studies <- "input_order"
  order_studies <- match.arg(order_studies)
  
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  #if (order_studies != "input_order"){
  # Order needs to be changed.
  if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    if (type_ev == "average") { 
      # average-evidence approach
      #OverallGoric <- 
      OverallWeight <- apply(Weights, 2, prod)^(1/S)
      OverallPrefHypo <- which(OverallGoric == max(OverallWeight))
    } else {
      # type_ev == "added" (or when "equal", because then it is overruled to be "added")
      type_ev = "added"
      #OverallGoric <- 
      OverallWeight <- apply(Weights, 2, prod)
      OverallPrefHypo <- which(OverallGoric == max(OverallWeight))
    }
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    orderStudies <- order(GORICA_weight_m[,OverallPrefHypo], decreasing = decreasing)
    #
    Weights <- Weights[orderStudies,]
  }
  #
  # Set colnames
  colnames(Weights) <- hypo_names
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    # TO DO check
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(Weights) <- study_names
  
  CumulativeWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos))
  colnames(CumulativeWeights) <- hypo_names
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Study 1"
  rownames(CumulativeWeights) <- c(sequence, "Final")
  #
  if (type_ev == "average") { 
    # average-ev approach
    CumulativeWeights[1, ] <- priorWeights * Weights[1, ] / sum( priorWeights * Weights[1, ] )
    for (s in 2:S) {
      CumulativeWeights[s, ] <- CumulativeWeights[(s-1), ] * Weights[s, ]^(1/s) / 
        sum( CumulativeWeights[(s-1), ] * Weights[s, ]^(1/s) )
    }
    # TO DO check of dit klopt en goed gaat
  } else {
    # type_ev == "added" (or when "equal", because then it is overruled to be "added")
    CumulativeWeights[1, ] <- priorWeights * Weights[1, ] / sum( priorWeights * Weights[1, ] )
    for (s in 2:S) {
      CumulativeWeights[s, ] <- CumulativeWeights[(s-1), ] * Weights[s, ] / 
        sum( CumulativeWeights[(s-1), ] * Weights[s, ] )
    }
  }
  CumulativeWeights[(S+1), ] <- CumulativeWeights[S, ]
  
  Final.weights <- CumulativeWeights[S, ]
  Final.ratio.GORICA.weights <- Final.weights %*% t(1/Final.weights)
  rownames(Final.ratio.GORICA.weights) <- hypo_names
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hypo_names)
  
  out <- list(type_ev                    = type_ev,
              GORICA_weight_m            = Weights,
              Cumulative_GORICA_weights  = CumulativeWeights,
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn_ICweights", "evSyn")
  
  return(out)
}


# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on the ratio of AIC or ORIC or GORIC or GORICA 
# weights or (Bayesian) posterior model probabilities
evSyn_ICratios <- function(object, ..., priorWeights = NULL, hypo_names = c(),
                           order_studies = c("input_order", "ascending", "descending"),
                           study_names = c()) {
  
  Weights <- object # Now, ratio of weights # TO DO check of onderstaande dan wel goed gaat
  S <- length(Weights)
  Weights <- do.call(rbind, Weights)
  NrHypos <- ncol(Weights)
  
  if (is.null(priorWeights)) {
    priorWeights <- rep(1/(NrHypos), (NrHypos))
  }
  # To make it sum to 1 (if it not already did)
  priorWeights <- priorWeights / sum(priorWeights) 
  # TO DO wel message meegeven dat dit gedaan is
  
  if (is.null(hypo_names)) {
    hypo_names <- paste0("H", 1:(NrHypos)) 
  }
  
  if (missing(order_studies)) 
    order_studies <- "input_order"
  order_studies <- match.arg(order_studies)
  
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  #if (order_studies != "input_order"){
  # Order needs to be changed.
  if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    if (type_ev == "average") { 
      # average-evidence approach
      #OverallGoric <- 
      OverallWeight <- apply(Weights, 2, prod)^(1/S)
      OverallPrefHypo <- which(OverallGoric == max(OverallWeight))
    } else {
      # type_ev == "added" (or when "equal", because then it is overruled to be "added")
      type_ev = "added"
      #OverallGoric <- 
      OverallWeight <- apply(Weights, 2, prod)
      OverallPrefHypo <- which(OverallGoric == max(OverallWeight))
    }
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    orderStudies <- order(GORICA_weight_m[,OverallPrefHypo], decreasing = decreasing)
    #
    Weights <- Weights[orderStudies,]
  }
  #
  # Set colnames
  colnames(Weights) <- hypo_names
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    # TO DO check
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(Weights) <- study_names
  
  CumulativeWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos))
  colnames(CumulativeWeights) <- hypo_names
  sequence <- paste0("Studies 1-", 1:S)
  sequence[1] <- "Study 1"
  rownames(CumulativeWeights) <- c(sequence, "Final")
  #
  if (type_ev == "average") { 
    # average-ev approach
    CumulativeWeights[1, ] <- priorWeights * Weights[1, ] / sum( priorWeights * Weights[1, ] )
    for (s in 2:S) {
      CumulativeWeights[s, ] <- CumulativeWeights[(s-1), ] * Weights[s, ]^(1/s) / 
        sum( CumulativeWeights[(s-1), ] * Weights[s, ]^(1/s) )
    }
    # TO DO check of dit klopt en goed gaat
  } else {
    # type_ev == "added" (or when "equal", because then it is overruled to be "added")
    CumulativeWeights[1, ] <- priorWeights * Weights[1, ] / sum( priorWeights * Weights[1, ] )
    for (s in 2:S) {
      CumulativeWeights[s, ] <- CumulativeWeights[(s-1), ] * Weights[s, ] / 
        sum( CumulativeWeights[(s-1), ] * Weights[s, ] )
    }
  }
  CumulativeWeights[(S+1), ] <- CumulativeWeights[S, ]
  
  Final.weights <- CumulativeWeights[S, ]
  Final.ratio.GORICA.weights <- Final.weights %*% t(1/Final.weights)
  rownames(Final.ratio.GORICA.weights) <- hypo_names
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hypo_names)
  
  out <- list(type_ev                    = type_ev,
              order_studies              = orderStudies,
              study_names                = study_names,
              GORICA_weight_m            = Weights,
              Cumulative_GORICA_weights  = CumulativeWeights,
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  # TO DO dit klopt dan toch niet qua naamgeving!
  #       Ms direct al WeightRatios noemen, dat helpt ms.
  
  class(out) <- c("evSyn_ICratios", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# list with goric objects
evSyn_gorica <- function(object, ..., type_ev = c("added", "equal", "average"), 
                         hypo_names = c(),
                         order_studies = c("input_order", "ascending", "descending"),
                         study_names = c()) {

  # Check if all objects are of type "con_goric"
  if (!all(vapply(object, function(x) inherits(x, "con_goric"), logical(1)))) {
    stop("\nrestriktor ERROR: the object must be a list with fitted objects from the goric() function", 
         call. = FALSE)
  }
  # TO DO check of allemaal gelijke type, alleen dan samennemen.
  # TO DO als small sample, dan ook sample_nobs nodig of kan het zonder?
  # TO DO check of in elke zelfde aantal hypotheses (met evt controle of failsafe ook - dit als message dan), anders werkt het ook niet
  
  # Create a list for the evSyn_LL.list function
  conList <- list(
    object = lapply(object, function(x) x$result$loglik),
    PT = lapply(object, function(x) x$result$penalty),
    type_ev = type_ev,
    hypo_names = hypo_names,
    study_names = study_names
  )
  
  # Call the evSyn_LL.list function and return the result
  result <- do.call(evSyn_LL, append(conList, list(...)))
  class(result) <- c(class(result), "evSyn_gorica")
  
  return(result)
}


evSyn_escalc <- function(data, yi_col = "yi", vi_cols = "vi", 
                         cluster_col = c("trial", "study", "author", "authors", "Trial", "Study", "Author", "Authors"),
                         outcome_col = NULL, ...) {
  
  results <- extract_est_vcov_outcomes(
    data = data,
    outcome_col = outcome_col,
    yi_col = yi_col,
    vi_cols = vi_cols,
    cluster_col = cluster_col
  )
  
  # Access the parameter estimates and vcov blocks
  yi_list <- results$yi_list
  vcov_blocks <- results$vcov_blocks
  
  
  conList <- list(
    object = yi_list,
    VCOV = vcov_blocks
  )
  
  # Call the evSyn_LL.list function and return the result
  result <- do.call(evSyn_est, append(conList, list(...)))
  class(result) <- c(class(result), "evSyn_escalc")
  
  return(result)
  
}
