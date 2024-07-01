plot.benchmark <- function(x, output_type = c("rgw", "rlw", "gw", "ld"), 
                           hypothesis_comparison = NULL,
                           percentiles = c(0.025, 0.975), x_lim = c(), ...) {
  
  # Ensure the object is of class 'benchmark_means'
  if (!inherits(x, "benchmark")) {
    stop("Invalid object. The object should be of class 'benchmark'.", call. = FALSE)
  }
  
  # Check if the output_type is valid
  output_type <- match.arg(output_type, c("rgw", "rlw", "gw", "ld"))
  # if (!output_type %in% valid_types) {
  #   stop("Restriktor ERROR: Invalid output_type. Valid types are 'gw', 'rgw', 'rlw', 'ld'.", call. = FALSE)
  # }
  
  legend_lab <- names(x$benchmarks_goric_weights)
  # first letter to upper-case
  goric_type <- toupper(x$type) #paste0(toupper(substring(x$type, 1, 1)), substring(x$type, 2))
  
  comparison <- x$comparison
  pref_hypo_name <- x$pref_hypo_name
  
  # Define the labels based on the output_type
  if (output_type == "gw") {
    df <- as.data.frame(x$combined_values$gw_combined)  
    sample_value <- x$benchmarks_goric_weights[[1]][1]
    xlabel <- paste(goric_type, "weights")
    title <- paste0("Benchmark: ", goric_type, "-weight distribution for preferred hypothesis ", pref_hypo_name)
  } else if (output_type == "rgw") {
    df <- as.data.frame(x$combined_values$rgw_combined)  
    sample_value <- x$benchmarks_ratio_goric_weights[[1]][, 1]
    xlabel <- paste0("Ratio ", goric_type, " weights")
    title <- paste0("Benchmark: Ratio ", goric_type, "-weight distribution for preferred hypothesis ", pref_hypo_name)
  } else if (output_type == "rlw") {
    df <- as.data.frame(x$combined_values$rlw_combined)  
    sample_value <- x$benchmarks_ratio_ll_weights[[1]][, 1]
    xlabel <- "Ratio likelihood weights"
    title <- paste0("Benchmark: Ratio likelihood weight distribution for preferred hypothesis ", pref_hypo_name)
  } else if (output_type == "ld") {
    df <- as.data.frame(x$combined_values$ld_combined)
    sample_value <- x$benchmarks_difLL[[1]][, 1]
    xlabel <- "Likelihood difference"
    title <- paste0("Benchmark: Likelihood difference distribution for preferred hypothesis ", pref_hypo_name)
  }
  
  if (output_type != "gw") {
    ratio_names <- rownames(x$benchmarks_ratio_goric_weights[[1]])
    result_vector <- outer(legend_lab, ratio_names, function(x, y) paste(x, " (", y, ")", sep = ""))
    result_vector <- as.vector(result_vector)
    names(df) <- result_vector
    if (!is.null(hypothesis_comparison) && output_type != "gw") {
      sample_value <- filter_vector(sample_value, hypothesis_comparison)
      df <- filter_columns(df, hypothesis_comparison)
      result_vector <- names(df)
    }
  } else {
    result_vector <- legend_lab
    names(df) <- result_vector
  }

  # Ensure percentiles is valid
  if (any(percentiles <= 0 | percentiles >= 1)) {
    stop("Restriktor ERROR: a percentile should be between 0 and 1.", call. = FALSE)
  }
  
  # Reshape the dataframe 
  df_long <- reshape(df, 
                     varying = list(1:length(df)),#list(names(df)[grepl("^pop", names(df))]), 
                     v.names = "Value", 
                     timevar = "Group", 
                     times = result_vector, #names(df)[grepl("^pop", names(df))], 
                     direction = "long")
  
  row.names(df_long) <- NULL
  
  # Rename the Group column to replace triple dots with equals sign
  #df_long$Group <- gsub("\\.\\.\\.", " = ", df_long$Group)
  df_long$Group <- factor(df_long$Group, levels = result_vector)
  
  # Calculate the custom confidence intervals
  # ci_df <- aggregate(Value ~ Group, data = df_long, function(x) {
  #   c(lower = quantile(x, probs = (1 - CI) / 2), 
  #     upper = quantile(x, probs = 1 - (1 - CI) / 2), mean = mean(x))
  # })
  # Calculate the custom percentile intervals
  percentile_df <- aggregate(Value ~ Group, data = df_long, function(x) {
    c(lower = quantile(x, probs = percentiles[1]), 
      upper = quantile(x, probs = percentiles[2]), mean = mean(x))
  })
  
  percentile_df <- data.frame(Group = percentile_df$Group, percentile_df$Value)
  names(percentile_df) <- c("Group", "lower", "upper", "mean")
  
  # Calculate maximum density values for each group 
  max_density <- max(sapply(split(df_long$Value, df_long$Group), function(x) max(density(x)$y)))
  
  # Prepare text for the annotations with fixed digit formatting
  # percentile_df$label <- paste(percentile_df$Group, ": ", sprintf("%.0f", quant * 100), "th percentile = [", 
  #                      sprintf("%.3f", percentile_df$lower), ", ", 
  #                      sprintf("%.3f", percentile_df$upper), "]", sep = "")

  first_group <- levels(factor(df_long$Group))[1]
  first_group_data <- subset(df_long, Group == first_group)
  percentile_first_group <- percentile_df[percentile_df$Group == first_group, ]
  
  group_color <- scales::brewer_pal(palette = "Set3")(length(unique(df_long$Group)))
  first_group_color <- group_color[1]
  
  percentile_first_group_lower <- round(percentile_first_group$lower, 3)
  percentile_first_group_upper <- round(percentile_first_group$upper, 3)
  

  #############################################
  sample_segments <- lapply(seq_along(sample_value), function(i) {
    label <- paste0(format_value(sample_value[i]), " (", names(sample_value)[i], ")")
    geom_segment(aes(x = sample_value[i], xend = sample_value[i], y = 0, yend = max_density,
                     color = label, linetype = label), linewidth = 1, show.legend = TRUE)
  })
  
  # Generate percentile labels dynamically
  percentile_label <- paste0("[", percentiles[1]*100, ", ", percentiles[2]*100, "]th Percentiles")
  observed_label <- paste0(sapply(sample_value, format_value), " (", names(sample_value), ")")
  
  linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  observed_linetypes <- linetypes[1:length(sample_value)]
  
  # Create the plot
  p <- ggplot(df_long, aes(x = Value, fill = Group)) +
    geom_density(alpha = 0.5) +  # Density plots for multiple groups
    geom_segment(aes(x = percentile_first_group_lower, 
                     xend = percentile_first_group_lower, y = 0, yend = max_density,
                     color = percentile_label, linetype = percentile_label), linewidth = 1) + 
    geom_segment(aes(x = percentile_first_group_upper, 
                     xend = percentile_first_group_upper, y = 0, yend = max_density,
                     color = percentile_label, linetype = percentile_label), linewidth = 1) + 
    ggtitle(title) + 
    xlab(xlabel) + ylab("Density") + 
    theme(axis.text = element_text(size = 11),
          axis.title.x = element_text(size = 12, margin = margin(t = 10)),
          axis.title.y = element_text(size = 12, margin = margin(r = 10)),
          plot.title = element_text(size = 12)) +
    scale_fill_brewer(palette = "Set3") +
    scale_color_manual(name = "", 
                       values = c(setNames(rep("red", length(sample_value)), observed_label),
                                  setNames(first_group_color, percentile_label)),
                       breaks = c(observed_label, percentile_label),
                       guide = "none") +
    scale_linetype_manual(name = "Observed Values", 
                          values = c(setNames(observed_linetypes, observed_label),
                                     setNames("solid", percentile_label))) + 
    theme(legend.key = element_rect(fill = "white")) +
    guides(fill = guide_legend(override.aes = list(color = NA)),  # Maintain filled boxes in legend
           linetype = guide_legend(override.aes = list(color = c(first_group_color, rep("red", length(sample_value))))))  # Set legend colors
  
  # Add the segments for sample values
  for (segment in sample_segments) {
    p <- p + segment
  }
  
  if (inherits(x, "benchmark_asymp")) {
     p <- p + labs(fill = "Effect-Size Estimates")
  } else {
     p <- p + labs(fill = "Effect-Size Cohens f")
  } 

  
  if (!is.null(x_lim) && length(x_lim) == 2) {
    p <- p + coord_cartesian(xlim = x_lim)
  } 

  # Render the plot to get the y-axis limits
  #plot_lims <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  
  # Calculate the maximum y position for the annotations
  #y_max <- plot_lims[2]
  
  # Calculate the spacing for the annotations
  #num_annotations <- nrow(percentile_df) + 1
  #spacing <- y_max / 20#(num_annotations + 1)
  
  # Determine the position for the annotation to ensure it stays within the plot area
  # if (sample_value > mean(range(df_long$Value))) {
  #   annotate_position <- sample_value - (0.1 * diff(range(df_long$Value)))
  #   hjust_value <- 1  # Align text to the left
  # } else {
  #   annotate_position <- sample_value + (0.1 * diff(range(df_long$Value)))
  #   hjust_value <- 0  # Align text to the right
  # }
  
  # Add sample value annotation with adjusted position and alignment
  #sample_value <- sprintf("%.3f", sample_value)
  
  # 
  # p <- p + annotate("text", x = annotate_position, y = y_max, 
  #                   label = paste("Observed value = ", sprintf("%.3f", sample_value), sep = ""), 
  #                   hjust = hjust_value, vjust = 1, size = 3, color = "black", family = "mono")
  
  return(p)
}
