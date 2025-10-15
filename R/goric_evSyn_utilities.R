extract_est_vcov_outcomes <- function(data, outcome_col = NULL, 
                                      yi_col = "yi", 
                                      vi_cols = "vi", 
                                      cluster_col = c("trial", "study", "author", "authors", "Trial", "Study", "Author", "Authors")
                                      ) {

  cluster_col_sum <- sum(cluster_col %in% names(data))
  cluster_col_length <- length(cluster_col)
  if (cluster_col_sum == 0 && cluster_col_length > 1) {
    # The data set does not contain a variable with one of the possible cluster column names
    stop(paste("\nrestriktor ERROR: The data does not contain a variable with one of the following names", sQuote(cluster_col), 
               ". Make sure to either match that or the specify the name using the cluster_col argument."))
  } else if (cluster_col_sum == 0 && cluster_col_length == 1) {
    # One specified column, which is not available in the data set
    stop(paste("\nrestriktor ERROR: The specified cluster_col", sQuote(cluster_col), 
               "is not found in the data."))
  } else {
    # Thus, the data contains at least on of the cluster_col names;
    # either one of the pre-specified ones or one(s) user-specified
    cluster_col <- cluster_col[cluster_col %in% names(data)][1]
    if (cluster_col_sum > 1) {
      # If there are more matches, then use the first one
      message(paste("\nrestriktor Message: The data contains multiple variables with a 'cluster_col' name, the following one is used:", sQuote(cluster_col)))
    }
  }
  
  if (!yi_col %in% names(data)) {
    stop(paste("\nrestriktor ERROR: The yi_col", sQuote(yi_col), 
               "is not found in the data."))
  }

  if (all(!vi_cols %in% names(data))) {
    stop(paste("\nrestriktor ERROR: The vi_cols", sQuote(vi_cols), 
               "are not found in the data."))
  }
  
  if (is.null(outcome_col) && "outcome" %in% names(data)) {
    # If no names specified, but there is a column with names, use that
    outcome_col <- "outcome" #unique(data$outcome)
  }
    
  yi_list <- list()
  vcov_blocks <- list()
  
  for (cluster_id in unique(data[[cluster_col]])) {
    cluster_data <- data[data[[cluster_col]] == cluster_id, ]
    
    # Extract outcomes and yi values for each trial
    if (is.null(outcome_col)) {
      outcomes <- "theta" # "yi"
      yi_vals <- setNames(cluster_data[[yi_col]], outcomes)
      yi_list[[paste(cluster_col, cluster_id)]] <- yi_vals
    } else {
      outcomes <- cluster_data[[outcome_col]]
      yi_vals <- setNames(cluster_data[[yi_col]], outcomes)
      yi_list[[paste(cluster_col, cluster_id)]] <- yi_vals
    }
    
    # Determine the number of outcomes for this trial
    num_outcomes <- nrow(cluster_data)
    if (num_outcomes > 1) {
      vi_cols <- paste0("v", 1:num_outcomes, "i")
    }
    vcov_matrix <- matrix(0, nrow = num_outcomes, ncol = num_outcomes)
    rownames(vcov_matrix) <- colnames(vcov_matrix) <- outcomes
    
    # Fill in the variance-covariance matrix
    if (num_outcomes == 1) {
      # Only one outcome: single variance, no covariances
      vcov_matrix[1, 1] <- cluster_data[[vi_cols[1]]][1]
    } else {
      # Multiple outcomes: fill diagonal with variances and off-diagonals with covariances
      vcov_matrix <- as.matrix(cluster_data[vi_cols])
    }
    
    # Store the matrix in vcov_blocks list
    vcov_blocks[[paste(cluster_col, cluster_id)]] <- vcov_matrix
  }
  
  # Return both lists as a named list
  return(list(yi_list = yi_list, vcov_blocks = vcov_blocks))
}
