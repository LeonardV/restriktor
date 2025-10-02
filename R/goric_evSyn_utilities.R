extract_est_vcov_outcomes <- function(data, outcome_col = NULL, yi_col = "yi", 
                                      vi_cols = "vi", cluster_col = "trial") {
  # TO DO
  #cluster_col = c("trial", "study", "author", "authors",
  #"Trial", "Study", "Author", "Authors")
  
  if (!cluster_col %in% names(data)) {
    stop(paste("Restriktor ERROR: The cluster_col", sQuote(cluster_col), 
               "is not found in the data."))
  }
  
  if (!yi_col %in% names(data)) {
    stop(paste("Restriktor ERROR: The yi_col", sQuote(yi_col), 
               "is not found in the data."))
  }

  if (all(!vi_cols %in% names(data))) {
    stop(paste("Restriktor ERROR: The vi_cols", sQuote(vi_cols), 
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
      outcomes <- "yi"
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
