extract_est_vcov_outcomes <- function(data, outcome_col = NULL, 
                                      yi_col = "yi", 
                                      vi_cols = "vi", 
                                      type = NULL,
                                      #type = c("gorica", "goricac"),
                                      study_sample_nobs = "ni",
                                      cluster_col = c("trial", "study", "author", "authors", "Trial", "Study", "Author", "Authors")) {
  
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
    stop(paste("\nrestriktor ERROR: The yi_col column called", sQuote(yi_col), 
               "is not found in the data."))
  }
  
  if (!is.null(type) && type %in% c("goric", "goricc", "gorica", "goricac")) {
    type_missing <- FALSE
  } else {
    type <- "gorica"
    type_missing <- TRUE
  }
  #
  if (type %in% c("goricc", "goricac") && !(study_sample_nobs %in% names(data))) {
    stop(paste("\nrestriktor ERROR: The study_sample_nobs column called", sQuote(study_sample_nobs), 
               "is not found in the data. \n",
               "This is required to compute the GORICAC."))
  }
  # Note that later on 'goricc' will be changed to 'goricac'.
  
  if (all(!vi_cols %in% names(data))) {
    stop(paste("\nrestriktor ERROR: The vi_cols column called", sQuote(vi_cols), 
               "are not found in the data."))
  }
  
  if (is.null(outcome_col) && "outcome" %in% names(data)) {
    # If no names specified, but there is a column with names, use that
    outcome_col <- "outcome" #unique(data$outcome)
  }
  
  yi_list <- list()
  if (type %in% c("goricc", "goricac")) {
    ni_list <- list()
  } else {
    ni_list <- NULL
  }
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
    
    # If type == 'goricac' or 'goricc': Extract ni values for each trial
    if (type %in% c("goricc", "goricac")) {
      ni_list[[paste(cluster_col, cluster_id)]] <- cluster_data[[study_sample_nobs]][1]
    }
      
    # Determine the number of outcomes for this trial
    num_outcomes <- nrow(cluster_data)
    # Note: the columns used for this cluster are stored in 'vi_cols_s', such
    # that the user-specified 'vi_cols' is not overwritten (and is thus also
    # used for the next clusters).
    if (num_outcomes > 1) {
      # Multiple outcomes: use the user-specified vi_cols if these consist of
      # one column per outcome; otherwise, use the escalc-type names v1i, v2i, ...
      if (length(vi_cols) == num_outcomes && all(vi_cols %in% names(data))) {
        vi_cols_s <- vi_cols
      } else {
        vi_cols_s <- paste0("v", 1:num_outcomes, "i")
      }
      if (!all(vi_cols_s %in% names(data))) {
        stop(paste0("\nrestriktor ERROR: For ", cluster_col, " ", sQuote(cluster_id), ", there are ",
                    num_outcomes, " outcomes, but the (co)variance columns ",
                    paste(sQuote(vi_cols_s[!vi_cols_s %in% names(data)]), collapse = ", "),
                    " are not found in the data. Specify one (co)variance column per outcome ",
                    "in 'vi_cols' or use the names ", paste(sQuote(paste0("v", 1:num_outcomes, "i")), collapse = ", "), "."),
             call. = FALSE)
      }
    } else {
      # One outcome: use the first user-specified vi_col (by default, 'vi'); if
      # that one is not available (or missing) for this cluster, use 'v1i' (if available).
      if (vi_cols[1] %in% names(data) && !is.na(cluster_data[[vi_cols[1]]][1])) {
        vi_cols_s <- vi_cols[1]
      } else if ("v1i" %in% names(data) && !is.na(cluster_data[["v1i"]][1])) {
        vi_cols_s <- "v1i"
      } else {
        stop(paste0("\nrestriktor ERROR: For ", cluster_col, " ", sQuote(cluster_id),
                    ", no (non-missing) variance is found in column ", sQuote(vi_cols[1]), "."),
             call. = FALSE)
      }
    }

    # Fill in the variance-covariance matrix
    if (num_outcomes == 1) {
      # Only one outcome: single variance, no covariances
      vcov_matrix <- matrix(cluster_data[[vi_cols_s]][1], nrow = 1, ncol = 1)
    } else {
      # Multiple outcomes: fill diagonal with variances and off-diagonals with covariances
      vcov_matrix <- as.matrix(cluster_data[vi_cols_s])
    }
    # Label the rows and columns by the outcomes (i.e., the parameter names)
    rownames(vcov_matrix) <- colnames(vcov_matrix) <- outcomes
    
    # Store the matrix in vcov_blocks list
    vcov_blocks[[paste(cluster_col, cluster_id)]] <- vcov_matrix
  }
  
  # Return both lists as a named list
  return(list(yi_list = yi_list, vcov_blocks = vcov_blocks, ni_list = ni_list))
}