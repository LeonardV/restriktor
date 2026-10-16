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
  
  if (is.null(IC_m)) {
    stop("restriktor ERROR: IC matrix is missing from the evSyn object.")
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
  
   
  # 'Study' or, when type_ev = "average" (multiple analyses of one study), 'Analysis'
  unit_cap <- .evSyn_unit_label(type_ev)
  unit_low <- .evSyn_unit_label(type_ev, capital = FALSE)
  if(all(object$study_names == 1:S)) {
    rownames <- paste0("Leave ", unit_cap, " ", object$study_names, " out:")
  } else {
    #rownames <- paste0("Leave '", object$study_names, "' (i.e., study nr. ", 1:S, ") out:")
    rownames <- paste0("Leave ", unit_low, " nr. ", 1:S, " called '", object$study_names, "' out:")
  }

  # Penalty factor (default 2), prior model probabilities (priorICweights),
  # and study weights, as used in the evidence synthesis itself.
  penalty_factor <- if (is.null(object$penalty_factor)) 2 else object$penalty_factor
  priorICweights <- object$priorICweights
  if (is.null(priorICweights)) {
    priorICweights <- rep(1 / ncol(IC_m), ncol(IC_m))
  }
  study_weights <- object$study_weights
  if (is.null(study_weights)) {
    study_weights <- rep(1 / S, S)
  }
  # The study weights are stored in the original (input) order of the studies;
  # re-order them such that they follow the (possibly re-ordered) rows of IC_m.
  if (!is.null(object$order_studies)) {
    study_weights <- study_weights[object$order_studies]
  }
  # Study weights summing to S (as in the evSyn functions)
  study_weights_S <- S * study_weights / sum(study_weights)
  
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
  
  for (s in seq_len(S)) {
    
    keep <- seq_len(S) != s
    nr_keep <- sum(keep)

    # Study weights of the remaining studies, rescaled such that they average
    # to 1 (as in the evSyn functions, where the weighted sums are divided by
    # the mean of the study weights of the studies included so far).
    w_keep <- study_weights_S[keep] / mean(study_weights_S[keep])

    # Same IC as in the evSyn functions, but without study s:
    # added:   weighted sum of IC values
    # equal:   -2 x weighted sum of LL values + penalty_factor x weighted average of PT values
    # average: weighted average of IC values
    OverallGoric[s, ] <- switch(
      type_ev,
      added   = colSums(IC_m[keep, , drop = FALSE] * w_keep),
      equal   = -2 * colSums(LL_m[keep, , drop = FALSE] * w_keep) +
        penalty_factor * colSums(PT_m[keep, , drop = FALSE] * w_keep) / nr_keep,
      average = colSums(IC_m[keep, , drop = FALSE] * w_keep) / nr_keep
    )

    best <- which(OverallGoric[s, ] == min(OverallGoric[s, ], na.rm = TRUE))
    OverallPrefHypo[s, 1L] <- paste(colnames(IC_m)[best], collapse = ", ")
    # IC weights, taking into account the prior model probabilities (priorICweights)
    expGW <- priorICweights * exp(-0.5 * (OverallGoric[s, ] - min(OverallGoric[s, ], na.rm = TRUE)))
    OverallGoricWeights[s, ] <- expGW / sum(expGW)
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

