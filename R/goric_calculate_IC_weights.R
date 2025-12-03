# Calculating IC weights based on IC values (AIC, ORIC, GORIC(A), BIC, SIC, ...)
# This function transforms IC values into IC weights: IC values denote the ordering 
# of hypotheses/models, while IC weights quantify the relative strength of 
# hypotheses/models.
calc_ICweights <- calculate_IC_weights <- function(IC, hypo_names = NULL) {
  
  # Check if IC is a vector or a matrix with one column
  if (!is.vector(IC) && !(is.matrix(IC) && ncol(IC) == 1)) {
    stop("The argument IC should either be a vector or a matrix with one column.")
  }
  
  NrHypos <- length(IC)
  
  # Set default hypothesis names if not provided
  if (is.null(hypo_names)) {
    hypo_names <- paste0("H", 1:NrHypos)
  }
  
  # Check if hypo_names has the correct number of elements and is all characters
  if (length(hypo_names) != NrHypos || any(!is.character(hypo_names))) {
    stop(paste("The argument 'hypo_names' should consist of NrHypos =", NrHypos, " elements (all characters)."))
  }
  
  # Set names for IC vector/matrix
  names(IC) <- hypo_names
  
  # Calculate minimum IC value
  minIC <- min(IC)
  
  # Calculate IC weights
  weight_m <- exp(-0.5 * (IC - minIC)) / sum(exp(-0.5 * (IC - minIC)))
    names(weight_m) <- hypo_names

  # Calculate ratio of IC weights
  ratio_IC_weights <- weight_m %*% t(1/weight_m)
    rownames(ratio_IC_weights) <- hypo_names
    colnames(ratio_IC_weights) <- paste0("vs_", hypo_names)
  
    
#  if (use_scientific) {
#    weight_m <- format(weight_m, scientific = use_scientific)
#    ratio_IC_weights <- format(ratio_IC_weights, scientific = use_scientific)
#  }
    
  out <- list(IC = IC, IC_weights = weight_m, ratio_IC_weights = ratio_IC_weights)
  class(out) <- c("goric_ICw", "list")
  
  return(out)
}


print.goric_ICw <- function(x, digits = max(3, getOption("digits") - 4), 
                            use_scientific = TRUE, ...) {
  df <- as.data.frame(x)
 
  # Format and print the dataframe
  formatted_df <- format(df, digits = digits, nsmall = digits, scientific = use_scientific)
  print(formatted_df, print.gap = 2, quote = FALSE)
}
# TO DO werkt niet naar behoren
#use_scientific = FALSE
#digits = 3
# vooral digits dacht ik.
# Sommig eblijven die 7 (= getOption("digits"))
# Komt omdat dan sommige getallen nul lijken te zijn ws.
# Ms dan alleen die wel scientific maken en dan met gevraagde aantal digits?
# Dus check ws of in de column een getal zit dat kleiner is dan 0,0...1 met aantal 0en afh van digits,
# als, dan scientific maken met die digits (evt +1?)

# Voorbeeld:
# library(lavaan)
# # specify the model
# HS.model_mgcfa <- "
# visual =~ x1 + x2 + x3
# textual =~ x4 + x5 + x6
# speed =~ x7 + x8 + x9
# "
# # configural invariance
# fit_ci <- cfa(HS.model_mgcfa, data = HolzingerSwineford1939, group = "school")
# # weak invariance
# fit_wi <- cfa(HS.model_mgcfa, data = HolzingerSwineford1939, group = "school", group.equal = "loadings")
# # strong invariance
# fit_si <- cfa(HS.model_mgcfa, data = HolzingerSwineford1939, group = "school", group.equal = c("intercepts",
#                                                                                                "loadings"))
# # strict invariance
# fit_strict <- cfa(HS.model_mgcfa, data = HolzingerSwineford1939, group = "school",
#                   group.equal = c("intercepts", "loadings", "residuals"))
# # model comparison with AIC
# AIC_meas.invar <- lavTestLRT(fit_ci, fit_wi, fit_si, fit_strict)$AIC
# hypo.names <- c("configural", "weak", "strong", "strict")
# AIC_weights <- calc_ICweights(AIC_meas.invar, hypo.names)
# AIC_weights$ratio_IC_weights
# AIC_weights$ratio_IC_weights["weak", ]
# #
# print(AIC_weights, use_scientific = FALSE, digits = 1)