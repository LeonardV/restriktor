leave1out.evSyn <- function(object) {
  # TO DO 
  # een fie leave1out() voor een evSyn object
  
  if (!any(class(object) %in% c("evSyn_est", "evSyn"))) {
    stop("restriktor ERROR: The leave1out function takes an 'evSyn' (or 'evSyn_est') object as input.")
  }
  
  type <- object$type
  type_ev <- object$type_ev
  S <- object$n_studies
  
  # IC values needed 
  if(object$type == "goric") {
    IC_m <- object$GORIC_m
  } else if(object$type == "goricc") {
    IC_m <- object$GORICC_m
  } else if(object$type == "gorica") {
    IC_m <- object$GORICA_m
  } else if(object$type == "goricac") {
    IC_m <- object$GORICAC_m
  } 
  if (type_ev == "equal") { 
    # equal-evidence approach
    # Then LL and PT needed and available
    LL_m <- object$LL_m
    PT_m <- object$PT_m
  }
  
  OverallGoric <- matrix(NA, nrow = S, ncol = dim(IC_m)[2])
  OverallPrefHypo <- matrix(NA, nrow = S, ncol = 1) # rep(NA, S)
  for (s in 1:S) {
    if (type_ev == "added") { 
      # added-evidence approach
      OverallGoric[s,] <- colSums(IC_m[-s,])
    } else if (type_ev == "equal") { 
      # equal-evidence approach
      OverallGoric[s,] <- colSums(-LL_m[-s,]) + colMeans(PT[-s,])
    } else if (type_ev == "average") { 
      # average-evidence approach
      OverallGoric[s,] <- colMeans(IC_m[-s,])
    }
    which <- which(OverallGoric[s,] == min(OverallGoric[s,]))
    OverallPrefHypo[s] <- colnames(IC_m)[which]
  }

  colnames(OverallGoric) <- colnames(IC_m)
  names <- paste0("Leave Study ", object$order_studies, " out:")
  rownames(OverallGoric) <- names
  #names(OverallPrefHypo) <- names
  rownames(OverallPrefHypo) <- names
  colnames(OverallPrefHypo) <- ""
  
  # Output
  if(object$type == "goric") {
    resultIC <- list(OverallGoric = OverallGoric)
  } else if(object$type == "goricc") {
    result <- list(OverallGoricc = OverallGoric)
  } else if(object$type == "gorica") {
    result <- list(OverallGorica = OverallGoric)
  } else if(object$type == "goricac") {
    result <- list(OverallGoricac = OverallGoric)
  } 
  result <- append(resultIC, 
                   list(OverallPrefHypo = OverallPrefHypo,
                        type = type, type_ev = type_ev)
                   )
  
  class(result) <- c("leave1out.evSyn")
  
  return(result)
  
}
