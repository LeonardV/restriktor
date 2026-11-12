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
 
  # Format each column separately, with 'digits' decimals. If use_scientific 
  # = TRUE, a column is printed in scientific notation (with 'digits' 
  # decimals) only when it contains a non-zero value that would otherwise be
  # printed as zero (i.e., smaller than 0.5 x 10^-digits in absolute value) 
  # or a value of 1e5 or larger. (Formatting the data frame as a whole, as 
  # done before, gave more digits than requested for some columns.)
  formatted_df <- as.data.frame(lapply(df, function(v) {
    if (!is.numeric(v)) return(v)
    nz <- v[!is.na(v) & v != 0]
    use_sci <- isTRUE(use_scientific) && length(nz) > 0 &&
      (any(abs(nz) < 0.5 * 10^(-digits)) || any(abs(nz) >= 1e5))
    out <- formatC(v, format = if (use_sci) "e" else "f", digits = digits)
    out[is.na(v)] <- ""
    out
  }), row.names = rownames(df), check.names = FALSE, stringsAsFactors = FALSE)
  print(formatted_df, print.gap = 2, quote = FALSE, right = TRUE)
  invisible(x)
}
