plot.benchmark <- function(x, output_type = "gw", CI = 0.95, x_lim = c(), ...) {
  
  # Ensure the object is of class 'benchmark_means'
  if (!inherits(x, "benchmark")) {
    stop("Invalid object. The object should be of class 'benchmark'.", call. = FALSE)
  }
  
  # Check if the output_type is valid
  valid_types <- c("gw", "rgw", "rlw", "ld")
  if (!output_type %in% valid_types) {
    stop("Restriktor ERROR: Invalid output_type. Valid types are 'gw', 'rgw', 'rlw', 'ld'.", call. = FALSE)
  }
  
  quant <- CI
  legend_lab <- names(x$benchmarks_goric_weights)
  # first letter to upper-case
  goric_type <- paste0(toupper(substring(x$type, 1, 1)), substring(x$type, 2))
  
  comparison <- x$comparison
  pref_hypo_name <- x$pref_hypo_name
  
  # Define the labels based on the output_type
  if (output_type == "gw") {
    df <- as.data.frame(x$combined_values$gw_combined)  
    sample_value <- x$benchmarks_goric_weights[[1]][1]
    xlabel <- paste0(goric_type, "-weights")
    title <- paste0("Benchmark: ", goric_type, "-weight distribution for preferred hypothesis ", pref_hypo_name)
  } else if (output_type == "rgw") {
    df <- as.data.frame(x$combined_values$rgw_combined)  
    df <- df[, !colMeans(df) == 1, drop = FALSE]
    sample_value <- x$benchmarks_ratio_goric_weights[[1]][1, 1]
    xlabel <- paste0("Ratio ", goric_type, "-weights")
    title <- paste0("Benchmark: Ratio ", goric_type, "-weight distribution for preferred hypothesis ", pref_hypo_name)
  } else if (output_type == "rlw") {
    df <- as.data.frame(x$combined_values$rlw_combined)  
    df <- df[, !colMeans(df) == 1, drop = FALSE]
    sample_value <- x$benchmarks_ratio_ll_weights[[1]][1, 1]
    xlabel <- "Ratio likelihood-weights"
    title <- paste0("Benchmark: Ratio likelihood weight distribution for preferred hypothesis ", pref_hypo_name)
  } else if (output_type == "ld") {
    df <- as.data.frame(x$combined_values$ld_combined)
    df <- df[, !colMeans(df) == 0, drop = FALSE]
    sample_value <- x$benchmarks_difLL[[1]][1, 1]
    xlabel <- "Likelihood difference"
    title <- paste0("Benchmark: Likelihood difference distribution for preferred hypothesis ", pref_hypo_name)
  }
  
  names(df) <- legend_lab
  
  # Ensure quant is valid
  if (quant <= 0 || quant >= 1) {
    stop("Restriktor ERROR: a confidence interval should be between 0 and 1.", call. = FALSE)
  }
  
  # Reshape the dataframe 
  df_long <- reshape(df, 
                     varying = list(1:length(df)),#list(names(df)[grepl("^pop", names(df))]), 
                     v.names = "Value", 
                     timevar = "Group", 
                     times = legend_lab, #names(df)[grepl("^pop", names(df))], 
                     direction = "long")
  
  row.names(df_long) <- NULL
  
  # Rename the Group column to replace triple dots with equals sign
  df_long$Group <- gsub("\\.\\.\\.", " = ", df_long$Group)
  df_long$Group <- factor(df_long$Group, levels = legend_lab)
  
  # Calculate the custom confidence intervals
  ci_df <- aggregate(Value ~ Group, data = df_long, function(x) {
    c(lower = quantile(x, probs = (1 - quant) / 2), 
      upper = quantile(x, probs = 1 - (1 - quant) / 2), mean = mean(x))
  })
  ci_df <- data.frame(Group = ci_df$Group, ci_df$Value)
  names(ci_df) <- c("Group", "lower", "upper", "mean")
  
  # Calculate maximum density values for each group 
  max_density <- max(sapply(split(df_long$Value, df_long$Group), function(x) max(density(x)$y)))
  
  # Prepare text for the annotations with fixed digit formatting
  ci_df$label <- paste(ci_df$Group, ": ", sprintf("%.0f", quant * 100), "% CI = [", 
                       sprintf("%.3f", ci_df$lower), ", ", 
                       sprintf("%.3f", ci_df$upper), "]", sep = "")

  first_group <- levels(factor(df_long$Group))[1]
  first_group_data <- subset(df_long, Group == first_group)
  ci_first_group <- ci_df[ci_df$Group == first_group, ]
  
  first_group_color <- scales::brewer_pal(palette = "Set3")(length(unique(df_long$Group)))[1]
  
  sample_value <- round(sample_value, 3)
  ci_first_group_lower <- round(ci_first_group$lower, 3)
  ci_first_group_upper <- round(ci_first_group$upper, 3)
  
  # Plotting with ggplot
  p <- ggplot(df_long, aes(x = Value, fill = Group)) +
    geom_density(alpha = 0.5) +  # Density plots for multiple groups
    geom_segment(aes(x = sample_value, xend = sample_value, y = 0, yend = max_density),
                 linetype = "solid", color = "red", linewidth = 1,
                 show.legend = FALSE) +
    geom_segment(aes(x = ci_first_group_lower, xend = ci_first_group_lower, y = 0, yend = max_density),
                 linetype = "solid", color = first_group_color, linewidth = 1) + 
    geom_segment(aes(x = ci_first_group_upper, xend = ci_first_group_upper, y = 0, yend = max_density),
                 linetype = "solid", color = first_group_color, linewidth = 1) + 
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 0, color = paste0(CI * 100, "% CI")),
                 linetype = "solid", linewidth = 0) +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 0, color = paste0("Observed value ", sample_value)),
                 linetype = "solid", linewidth = 0) +
    ggtitle(title) + 
    xlab(xlabel) + ylab("Density") + 
    {if (inherits(x, "benchmark_asymp")) {
      labs(fill = "Effect-Size")
    } else {
      labs(fill = "Effect-Size f")
    }}  +
    theme(axis.text = element_text(size = 11),
          axis.title.x = element_text(size = 12, margin = margin(t = 10)),
          axis.title.y = element_text(size = 12, margin = margin(r = 10)),
          plot.title = element_text(size = 12)) +
    scale_fill_brewer(palette = "Set3") +
    scale_color_manual(name = c("", ""), 
                       values = c(setNames(first_group_color, paste0(CI * 100, "% CI")),
                                  setNames("red", paste0("Observed value ", sample_value)))) +
    theme(legend.key = element_rect(fill = "white"))  
  
   if (!is.null(x_lim) && length(x_lim) == 2) {
     p <- p + coord_cartesian(xlim = x_lim)
   }
  
  # Render the plot to get the y-axis limits
  #plot_lims <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  
  # Calculate the maximum y position for the annotations
  #y_max <- plot_lims[2]
  
  # Calculate the spacing for the annotations
  #num_annotations <- nrow(ci_df) + 1
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
