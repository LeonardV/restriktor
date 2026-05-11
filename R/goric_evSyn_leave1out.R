leave1out <- function(object, ...) {
  UseMethod("leave1out")
}


leave1out.evSyn <- function(object, ...) {
  
  type <- object$type
  
  if (!is.null(type) && type %in% c("goric", "goricc", "gorica", "goricac")) {
    type_missing <- FALSE
  } else {
    type <- "gorica"
    type_missing <- TRUE
  }
  
  type_ev <- object$type_ev
  S <- object$n_studies
  
  IC_m <- switch(
    type,
    goric   = object$GORIC_m,
    goricc  = object$GORICC_m,
    gorica  = object$GORICA_m,
    goricac = object$GORICAC_m
  )
  
  if (is.null(IC_m)) {
    stop("restriktor ERROR: IC matrix is missing from the evSyn object.")
  }
  
  if (!type_ev %in% c("added", "equal", "average")) {
    stop("restriktor ERROR: unknown evidence-synthesis type in 'type_ev'.")
  }
  
  if (S < 2L) {
    stop("restriktor ERROR: leave1out() requires at least two studies.")
  }
  
  if (type_ev == "equal") {
    LL_m <- object$LL_m
    PT_m <- object$PT_m
    
    if (is.null(LL_m) || is.null(PT_m)) {
      stop("restriktor ERROR: LL_m and PT_m are required for type_ev = 'equal'.")
    }
  }
  
  OverallGoric <- matrix(
    NA_real_,
    nrow = S,
    ncol = ncol(IC_m),
    dimnames = list(
      paste0("Leave Study ", object$order_studies, " out:"),
      colnames(IC_m)
    )
  )
  
  OverallPrefHypo <- matrix(
    NA_character_,
    nrow = S,
    ncol = 1L,
    dimnames = list(rownames(OverallGoric), "")
  )
  
  for (s in seq_len(S)) {
    
    keep <- seq_len(S) != s
    
    OverallGoric[s, ] <- switch(
      type_ev,
      added   = colSums(IC_m[keep, , drop = FALSE]),
      equal   = colSums(-LL_m[keep, , drop = FALSE]) +
        colMeans(PT_m[keep, , drop = FALSE]),
      average = colMeans(IC_m[keep, , drop = FALSE])
    )
    
    best <- which(OverallGoric[s, ] == min(OverallGoric[s, ], na.rm = TRUE))
    OverallPrefHypo[s, 1L] <- paste(colnames(IC_m)[best], collapse = ", ")
  }
  
  resultIC <- switch(
    type,
    goric   = list(OverallGoric = OverallGoric),
    goricc  = list(OverallGoricc = OverallGoric),
    gorica  = list(OverallGorica = OverallGoric),
    goricac = list(OverallGoricac = OverallGoric)
  )
  
  result <- c(
    resultIC,
    list(
      OverallPrefHypo = OverallPrefHypo,
      type_ev = type_ev
    )
  )
  
  if (!type_missing) {
    result$type <- type
  }
  
  class(result) <- "leave1out.evSyn"
  
  result
}


leave1out.evSyn_est <- function(object, ...) {
  leave1out.evSyn(object, ...)
}


leave1out.default <- function(object, ...) {
  stop(
    "restriktor ERROR: leave1out() takes an 'evSyn' or 'evSyn_est' object as input.",
    call. = FALSE
  )
}