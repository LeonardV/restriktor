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
    colnames(ratio_IC_weights) <- paste0("vs. ", hypo_names)
  
  out <- list(IC = IC, IC_weights = weight_m, ratio_IC_weights = ratio_IC_weights)
  class(out) <- c("goric_ICw", "list")
  
  return(out)
}


print.goric_ICw <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  x <- as.list(x)
  model_names <- attr(x$IC, "names")
  rownames(x$ratio_IC_weights) <- NULL
  
  # Create the dataframe
  df <- data.frame(model = model_names,
                   IC.values = x$IC,
                   IC.weights = x$IC_weights,
                   ratio.IC.weights = x$ratio_IC_weights)
  rownames(df) <- NULL
  
  # Correcting the column names
  names(df)[4:ncol(df)] <- gsub("vs..", "vs.", names(df)[4:ncol(df)]) 
  
  # Format and print the dataframe
  formatted_df <- format(df, digits = digits, nsmall = digits, scientific = FALSE)
  print(formatted_df, print.gap = 2, quote = FALSE)
}



