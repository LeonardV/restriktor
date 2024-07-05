plot.benchmark <- function(x, output_type = c("rgw", "rlw", "gw", "ld"), 
                           percentiles = c(0.05, 0.95), x_lim = c(), 
                           alpha = 0.50, nrow_grid = NULL, ncol_grid = 1, 
                           distr_grid = FALSE, ...) {
  
  # Ensure the object is of class 'benchmark_means'
  if (!inherits(x, "benchmark")) {
    stop("Invalid object. The object should be of class 'benchmark'.", call. = FALSE)
  }
  
  # Check if the output_type is valid
  output_type <- match.arg(output_type, c("rgw", "rlw", "gw", "ld"))
  
  legend_lab <- names(x$benchmarks_goric_weights)
  # first letter to upper-case
  #paste0(toupper(substring(x$type, 1, 1)), substring(x$type, 2))
  goric_type <- toupper(x$type) 
  comparison <- x$comparison
  pref_hypo_name <- x$pref_hypo_name
  
  # Define the labels based on the output_type
  if (output_type == "gw") {
    DATA <- x$combined_values$gw_combined
    sample_value <- x$benchmarks_goric_weights[[1]][1]
    xlabel <- paste(goric_type, "weights")
    title <- paste0("Benchmark: ", goric_type, "-weight Distribution for Preferred Hypothesis ", pref_hypo_name)
  } else if (output_type == "rgw") {
    DATA <- x$combined_values$rgw_combined
    sample_value <- x$benchmarks_ratio_goric_weights[[1]][, 1, drop = FALSE]
    xlabel <- paste0("Ratio ", goric_type, " weights")
    title <- paste0("Benchmark: Ratio ", goric_type, "-weight Distribution for Preferred Hypothesis ", pref_hypo_name)
  } else if (output_type == "rlw") {
    DATA <- x$combined_values$rlw_combined  
    sample_value <- x$benchmarks_ratio_ll_weights[[1]][, 1, drop = FALSE]
    xlabel <- "Ratio likelihood weights"
    title <- paste0("Benchmark: Ratio Likelihood Weight Distribution for Preferred Hypothesis ", pref_hypo_name)
  } else if (output_type == "ld") {
    DATA <- x$combined_values$ld_combined
    sample_value <- x$benchmarks_difLL[[1]][, 1, drop = FALSE]
    xlabel <- "Likelihood difference"
    title <- paste0("Benchmark: Likelihood Difference Distribution for Preferred Hypothesis ", pref_hypo_name)
  }
  
  # -------------------------------------------------------------------------
  new_combined_values <- lapply(names(DATA), function(list_name) { 
    colnames <- colnames(DATA[[list_name]])
    new_colnames <- construct_colnames(list_name, colnames, pref_hypo_name)
    colnames(DATA[[list_name]]) <- new_colnames
    DATA[[list_name]]
  })
  
  DATA <- new_combined_values
  
  combine_matrices_cbind <- function(lst) {
    df_list <- lapply(lst, as.data.frame)
    combined_df <- do.call(cbind, df_list)
    return(combined_df)
  }
  
  df <- combine_matrices_cbind(DATA)
  # -------------------------------------------------------------------------
  
  # Reshape the dataframe 
  df_long <- reshape(df, 
                     varying = list(1:length(df)),
                     v.names = "Value", 
                     timevar = "Group", 
                     times = names(df), 
                     direction = "long")
  row.names(df_long) <- NULL
  
  if (inherits(x, "benchmark_asymp")) {
    df_long$Group <- gsub("pop_est", "Population estimates", df_long$Group)
  } else {
    df_long$Group <- gsub("pop_es", "Effect-size", df_long$Group)
  } 
  
  # Rename the Group column to replace triple dots with equals sign
  df_long$Group <- factor(df_long$Group)
  
  percentile_df <- aggregate(Value ~ Group, data = df_long, function(x) {
    c(lower = quantile(x, probs = percentiles[1]), 
      upper = quantile(x, probs = percentiles[2]), mean = mean(x))
  })
  
  percentile_df <- data.frame(Group = percentile_df$Group, percentile_df$Value)
  names(percentile_df) <- c("Group", "lower", "upper", "mean")
  
  first_group <- levels(factor(df_long$Group))[1]
  first_group_data <- subset(df_long, Group == first_group)
  percentile_first_group <- percentile_df[percentile_df$Group == first_group, ]
  
  group_color <- scales::brewer_pal(palette = "Set3")(length(unique(df_long$Group)))
  first_group_color <- group_color[1]
  
  percentile_first_group_lower <- round(percentile_first_group$lower, 3)
  percentile_first_group_upper <- round(percentile_first_group$upper, 3)
  
  if (!output_type == "gw") {
    # sub plots per hypo, dus h1 vs h2, maar wel alle ES in dezelfde plot
    df_long$Group_pop_values <- sub("\\s*\\(.*\\)", "", df_long$Group)
    df_long$Group_hypo_comparison <- trimws(gsub("\\(|\\)", "",  extract_in_parentheses(df_long$Group)))
    sample_value <- setNames(as.vector(t(sample_value)), rownames(sample_value))
    
    sample_value_df <- data.frame(
      Group_hypo_comparison = names(sample_value),
      sample_value = as.numeric(sample_value),
      lb_first_group = percentile_first_group_lower,
      ub_first_group = percentile_first_group_upper,
      first_group_color = first_group_color,
      stringsAsFactors = FALSE
    )
    df_long <- merge(df_long, sample_value_df, by = "Group_hypo_comparison", all.x = TRUE)
  } else {
    df_long$Group <- paste(df_long$Group, "()")
    df_long$Group_pop_values <- sub("\\s*\\(.*\\)", "", df_long$Group)
    df_long$Group_hypo_comparison <- gsub("\\(|\\)", "",  extract_in_parentheses(df_long$Group))
    sample_value <- as.vector(sample_value)
    df_long <- suppressWarnings(cbind(df_long, sample_value, 
                                      lb_first_group = percentile_first_group_lower,
                                      ub_first_group = percentile_first_group_upper,
                                      first_group_color = first_group_color))
  }
  
  group <- unique(df_long$Group_hypo_comparison)
  # here we start to create the plot ----------------------------------------
  plots <- plot_all_groups(plot_df = df_long, 
                           groups = group, 
                           title = title,
                           xlabel = xlabel,
                           x_lim = x_lim,
                           alpha = alpha,
                           distr_grid = distr_grid,
                           percentiles = percentiles)
  
  plots <- do.call(grid.arrange, c(plots, ncol = ncol_grid, nrow = nrow_grid))
  
  return(plots)
}
