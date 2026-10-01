leave1studyout <- function(object, ...) {
  UseMethod("leave1studyout")
}

leave1studyout.default <- function(object, ...) {
  stop(
    "restriktor ERROR: leave1studyout() takes an 'evSyn' object as input.",
    call. = FALSE
  )
}

leave1studyout.evSyn <- function(object, ...) {
  
  type <- object$type
  
  if (!is.null(type) && type %in% c("goric", "goricc", "gorica", "goricac")) {
    type_missing <- FALSE
  } else {
    type <- "gorica"
    type_missing <- TRUE
  }
  
  type_ev <- object$type_ev
  S <- object$n_studies
  
  if (is.null(object$GORICA_m)) {
    IC_m <- switch(
      type,
      goric   = object$GORIC_m,
      goricc  = object$GORICC_m,
      gorica  = object$GORICA_m,
      goricac = object$GORICAC_m
    )
  } else {
    IC_m <- object$GORICA_m
  }
  
  if (is.null(object$GORICA_weight_m)) {
    ICw_m <- switch(
      type,
      goric   = object$GORIC_weight_m,
      goricc  = object$GORICC_weight_m,
      gorica  = object$GORICA_weight_m,
      goricac = object$GORICAC_weight_m
    )
  } else {
    ICw_m <- object$GORICA_weight_m
  }
  
  # Prior IC weights (used in the same way as in evSyn() when determining the
  # IC weights)
  priorICweights <- object$priorICweights

  # Input consisting of IC weights or ratios of IC weights: no IC values, but
  # the differences in IC values (vs a reference hypothesis) can be used,
  # since these lead to the same IC weights.
  IC_is_diff <- FALSE
  if (is.null(IC_m) && !is.null(object$ICdiff_m)) {
    # ratios of IC weights: differences in IC values vs the reference hypothesis
    IC_m <- object$ICdiff_m
    IC_is_diff <- TRUE
  } else if (is.null(IC_m) && !is.null(ICw_m)) {
    # IC weights: the study-specific weights in the evSyn object include the
    # prior IC weights, which are removed here (and included again below).
    W <- ICw_m[, , drop = FALSE]
    if (!is.null(priorICweights) && length(priorICweights) == ncol(W)) {
      W <- sweep(W, 2, ifelse(priorICweights > 0, priorICweights, NA), "/")
      W[, priorICweights == 0] <- 0
    }
    W <- W / rowSums(W)
    IC_m <- -2 * log(W)
    IC_is_diff <- TRUE
  }

  if (is.null(IC_m)) {
    stop("restriktor ERROR: IC matrix is missing from the evSyn object.")
  }
  if (is.null(ICw_m)) {
    ICw_m <- IC_m
  }
  
  if (!type_ev %in% c("added", "equal", "average")) {
    stop("restriktor ERROR: unknown evidence-synthesis type in 'type_ev'.")
  }
  
  if (S < 2L) {
    stop("restriktor ERROR: leave1studyout() requires at least two studies.")
  }
  
  if (type_ev == "equal") {
    LL_m <- object$LL_m
    PT_m <- object$PT_m
    
    if (is.null(LL_m) || is.null(PT_m)) {
      stop("restriktor ERROR: LL_m and PT_m are required for type_ev = 'equal'.")
    }
  }
  # penalty factor (IC = -2 * LL + penalty_factor * PT), as used in evSyn()
  penalty_factor <- object$penalty_factor
  if (is.null(penalty_factor)) {
    penalty_factor <- 2
  }
  
   
  if(all(object$study_names == 1:S)) {
    rownames <- paste0("Leave Study ", object$study_names, " out:")
  } else {
    #rownames <- paste0("Leave '", object$study_names, "' (i.e., study nr. ", 1:S, ") out:")
    rownames <- paste0("Leave study nr. ", 1:S, " called '", object$study_names, "' out:")
  }
  
  OverallGoric <- matrix(
    NA_real_,
    nrow = S,
    ncol = ncol(IC_m),
    dimnames = list(
      ##paste0("Leave Study ", object$order_studies, " out:"),
      #paste0("Leave '", object$study_names, "' out:"),
      rownames,
      colnames(IC_m)
    )
  )
  
  OverallGoricWeights <- matrix(
    NA_real_,
    nrow = S,
    ncol = ncol(ICw_m),
    dimnames = list(
      #paste0("Leave Study ", object$order_studies, " out:"),
      #paste0("Leave '", object$study_names, "' out:"),
      rownames,
      colnames(ICw_m)
    )
    # TO DO Als weight bijna een, dan geeft het 1 en bijv niet 1.000 of .999.
    #       Ik denk dat je met 1 denkt dat het exact 1 is... 
    #       Graag aanpassen, maar hoe&Waar?
    # TO DO Sowieso overal kijken naar hoeveel decimalen printen....
  )
  
  OverallPrefHypo <- matrix(
    NA_character_,
    nrow = S,
    ncol = 1L,
    dimnames = list(rownames(OverallGoric), "")
  )
  
  # Study weights and prior IC weights, used in the same way as in evSyn():
  # - the study weights of the remaining studies are rescaled such that they 
  #   sum to the number of remaining studies (S - 1), 
  # - the prior IC weights are used when determining the IC weights.
  study_weights <- object$study_weights
  if (is.null(study_weights)) {
    study_weights <- rep(1/S, S)
  }
  if (is.null(priorICweights) || length(priorICweights) != ncol(IC_m)) {
    priorICweights <- rep(1/ncol(IC_m), ncol(IC_m))
  }
  
  for (s in seq_len(S)) {
    
    keep <- seq_len(S) != s
    # number of remaining studies with a positive study weight (a study with 
    # weight 0 contributes nothing, as in evSyn())
    S_keep <- sum(study_weights[keep] > 0)
    # rescaled study weights (sum to S_keep); if the remaining studies all 
    # have weight zero, there is no evidence (IC values of 0, so the IC 
    # weights equal the prior IC weights), as in evSyn().
    if (S_keep > 0) {
      w_keep <- S_keep * study_weights[keep] / sum(study_weights[keep])
    } else {
      S_keep <- 1
      w_keep <- rep(0, sum(keep))
    }
    
    OverallGoric[s, ] <- switch(
      type_ev,
      added   = colSums(IC_m[keep, , drop = FALSE] * w_keep),
      equal   = -2 * colSums(LL_m[keep, , drop = FALSE] * w_keep) +
        penalty_factor * colSums(PT_m[keep, , drop = FALSE] * w_keep) / S_keep,
      average = colSums(IC_m[keep, , drop = FALSE] * w_keep) / S_keep
    )
    
    # IC weights (incl. prior IC weights; computed on the log scale) and the
    # preferred hypothesis, i.e., the one with the highest (prior-weighted)
    # IC weight (and not the one with the lowest IC value).
    OverallGoricWeights[s, ] <- ic_weights_log(OverallGoric[s, ], priorICweights)
    best <- which(OverallGoricWeights[s, ] == max(OverallGoricWeights[s, ], na.rm = TRUE))
    OverallPrefHypo[s, 1L] <- paste(colnames(IC_m)[best], collapse = ", ")
  }
  
  resultIC <- switch(
    type,
    goric   = list(OverallGoric = OverallGoric),
    goricc  = list(OverallGoricc = OverallGoric),
    gorica  = list(OverallGorica = OverallGoric),
    goricac = list(OverallGoricac = OverallGoric)
  )
  
  resultICw <- switch(
    type,
    goric   = list(OverallGoricWeights = OverallGoricWeights),
    goricc  = list(OverallGoriccWeights = OverallGoricWeights),
    gorica  = list(OverallGoricaWeights = OverallGoricWeights),
    goricac = list(OverallGoricacWeights = OverallGoricWeights)
  )
  
  result <- c(
    resultIC,
    resultICw,
    list(
      OverallPrefHypo = OverallPrefHypo,
      type_ev = type_ev
    )
  )
  if (IC_is_diff) {
    # Based on IC weights or ratios of IC weights: the 'IC values' are
    # differences in IC values (vs a reference hypothesis), not IC values.
    result$IC_is_diff <- TRUE
  }
  
  if (!type_missing) {
    result$type <- type
  }
  
  class(result) <- "leave1studyout.evSyn"
  
  result
}

# TO DO:
# NB som gaat sowieso omlaag en het hangt van N af hoeveel!!!
# Dus eigenlijk zou je ook voorspelling moeten geven
# en/of verdisconteren!

