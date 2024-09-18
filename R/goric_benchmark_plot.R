plot.benchmark <- function(x, output_type = c("rgw", "rlw", "gw", "ld"), 
                           percentiles = NULL, x_lim = c(), log_scale = FALSE,
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
  
  # values are extracted from benchmark
  if (is.null(percentiles)) {
    # exclude sample
    percentile_names <- colnames(x$benchmarks_goric_weights[[1]])[-1]
    percentiles <- as.numeric(sub("%", "", percentile_names)) / 100
  }
  
  # Define the labels based on the output_type
  if (output_type == "gw") {
    DATA <- x$combined_values$gw_combined
    sample_value <- x$benchmarks_goric_weights[[1]][1]
    xlabel <- paste(goric_type, "Weights")
    title <- paste0("Benchmark: ", goric_type, "-Weights Distribution for Preferred Hypothesis ", pref_hypo_name)
  } else if (output_type == "rgw") {
    DATA <- x$combined_values$rgw_combined
    sample_value <- x$benchmarks_ratio_goric_weights[[1]][, 1, drop = FALSE]
    xlabel <- paste("Ratio", goric_type, "Weights")
    title <- paste0("Benchmark: Ratio-", goric_type, "-Weights Distribution for Preferred Hypothesis ", pref_hypo_name)
  } else if (output_type == "rlw") {
    DATA <- x$combined_values$rlw_combined  
    sample_value <- x$benchmarks_ratio_ll_weights[[1]][, 1, drop = FALSE]
    xlabel <- "Ratio Log-likelihood Weights"
    title <- paste0("Benchmark: Ratio-Likelihood-Weights Distribution for Preferred Hypothesis ", pref_hypo_name)
  } else if (output_type == "ld") {
    DATA <- x$combined_values$ld_combined
    sample_value <- x$benchmarks_difLL[[1]][, 1, drop = FALSE]
    xlabel <- "Log-likelihood Difference"
    title <- paste0("Benchmark: Log-likelihood-Difference Distribution for Preferred Hypothesis ", pref_hypo_name)
  }
  
  # -------------------------------------------------------------------------
  new_combined_values <- lapply(names(DATA), function(list_name) { 
    colnames <- colnames(DATA[[list_name]])
    new_colnames <- construct_colnames(list_name, colnames, pref_hypo_name)
    colnames(DATA[[list_name]]) <- new_colnames
    DATA[[list_name]]
  })
  
  DATA <- new_combined_values
  

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
    quantile(x, probs = percentiles, names = TRUE, na.rm = TRUE)
  })
  
  percentile_df <- data.frame(Group = percentile_df$Group, percentile_df$Value, check.names = FALSE)
  
  n_plots <- nrow(x$benchmarks_ratio_goric_weights[[1]])
  first_group <- levels(factor(df_long$Group))[1:n_plots]
  first_group_data <- subset(df_long, Group %in% first_group)
  # first groups are the ones with no-effect or ES = 0
  percentile_first_group <- percentile_df[percentile_df$Group %in% first_group, ]
  rownames(percentile_first_group) <- NULL
  
  group_color <- scales::brewer_pal(palette = "Set3")(length(unique(df_long$Group)))
  first_group_color <- group_color[1]

  if (!output_type == "gw") {
    # sub plots per hypo, dus h1 vs h2, maar wel alle ES in dezelfde plot
    df_long$Group_pop_values <- sub("\\s*\\(.*\\)", "", df_long$Group)
    df_long$Group_hypo_comparison <- trimws(gsub("\\(|\\)", "",  extract_in_parentheses(df_long$Group)))
    sample_value <- setNames(as.vector(t(sample_value)), rownames(sample_value))
    
    sample_value_df <- data.frame(
      Group_hypo_comparison = names(sample_value),
      sample_value = as.numeric(sample_value),
      first_group_color = first_group_color,
      stringsAsFactors = FALSE
    )
    sample_value_df <- cbind(sample_value_df, percentile_first_group[-1])
    df_long <- merge(df_long, sample_value_df, by = "Group_hypo_comparison", all.x = TRUE)
  } else {
    df_long$Group <- paste(df_long$Group, "()")
    df_long$Group_pop_values <- sub("\\s*\\(.*\\)", "", df_long$Group)
    df_long$Group_hypo_comparison <- gsub("\\(|\\)", "",  extract_in_parentheses(df_long$Group))
    sample_value <- as.vector(sample_value)
    df_long <- suppressWarnings(cbind(df_long, sample_value, 
                                      percentile_first_group,
                                      first_group_color = first_group_color)
    )
  }
  
  attr(df_long, "output_type") <- output_type
  
  group <- unique(df_long$Group_hypo_comparison)
  # here we start to create the plot ----------------------------------------
  plot_list <- plot_all_groups(plot_df = df_long, 
                               groups = group, 
                               title = title,
                               xlabel = xlabel,
                               x_lim = x_lim,
                               alpha = alpha,
                               distr_grid = distr_grid,
                               percentiles = percentiles,
                               log_scale = log_scale)
  
  combined_plots <- list(
    plots = plot_list,
    ncol = ncol_grid,
    nrow = nrow_grid
  )
  
  #combined_plots <- do.call(gridExtra::arrangeGrob, c(plot_list, ncol = ncol_grid, nrow = nrow_grid))

  class(combined_plots) <- c("benchmark_plot", class(combined_plots))
  
  return(combined_plots)
}


# benchmark plots
create_density_plot <- function(plot_df, group_comparison, title, xlabel,
                                x_lim = NULL, alpha = 0.5, distr_grid = FALSE,
                                percentiles = NULL, log_scale = FALSE) {
  
  output_type <- attr(plot_df, "output_type")
  df_subset <- subset(plot_df, Group_hypo_comparison == group_comparison)
  
  if (!is.null(df_subset$Group_hypo_comparison)) {
    title <- paste(title, "vs.", sub(".*vs\\. ", "", unique(df_subset$Group_hypo_comparison)))
  }

  # Calculate densities for each group
  # density_data <- df_subset %>%
  #   group_by(Group_pop_values) %>%
  #   do(calculate_density(., "Value", unique(.$sample_value))) %>%
  #   ungroup()
  
  group_levels <- unique(df_subset$Group_pop_values)
  density_data <- do.call(rbind, lapply(group_levels, function(group) {
    data <- subset(df_subset, Group_pop_values == group)
    data$Group_pop_values <- factor(group)
    data$sample_value <- unique(data$sample_value)
    data
    # dens <- calculate_density(data, var = "Value", sample_value = data$sample_value[1])
    # dens$Group_pop_values <- factor(group)
    # dens$sample_value <- unique(data$sample_value)
    # dens
  }))
  
  # 
  percentile_values <- as.numeric(df_subset[, paste0(percentiles*100, "%")][1, ])
  #percentile_labels <- paste0("[", percentiles * 100, "]th Percentile = ", sprintf("%.3f", percentile_values))
  percentile_labels <- paste0(percentiles * 100, "th Percentile = ", sprintf("%.3f", percentile_values))
  formatted_sample_value <- sprintf("Sample Value = %.3f", unique(df_subset$sample_value)[1])
  
  
  percentile_df <- data.frame(
    percentile_value = percentile_values,
    percentile_label = percentile_labels, 
    color = df_subset$first_group_color[1],
    width = 2
  )
  # add sample value
  percentile_df <- rbind(data.frame(percentile_value = unique(df_subset$sample_value)[1], 
                                    percentile_label = formatted_sample_value, 
                                    color = "red", width = 2),  
                         percentile_df)
  
  # add segment: if gw or ld, then x = 0, else x = 1
  vert_line_width <- 1
  if (output_type == "ld") {
    vert_line_label <- "Equal Fit"
    percentile_df <- rbind(percentile_df, 
                           data.frame(percentile_value = 0, 
                                      percentile_label = vert_line_label, 
                                      color = "darkgrey", width = vert_line_width))
  } else if (output_type %in% c("rgw", "rlw")) {
    if (output_type == "rgw") {
      vert_line_label <- "Equal Support"  
    } else if (output_type == "rlw") {
      vert_line_label <- "Equal Fit"  
    }
    percentile_df <- rbind(percentile_df, 
                           data.frame(percentile_value = 1, 
                                      percentile_label = vert_line_label, 
                                      color = "darkgrey", width = vert_line_width))
  }
  
  percentile_df$Group_pop_values <- unique(df_subset$Group_pop_values)[1]
  # make sure the right order is preserved in the legend
  percentile_df$percentile_label <- factor(percentile_df$percentile_label, levels = percentile_df$percentile_label) 
  
  #sample_value_df <- percentile_df[percentile_df$percentile_label == "sample_value", ]
  
  # Plot maken
  p <- ggplot(density_data, aes(x = Value, fill = Group_pop_values)) +
      stat_density(
        aes(ymin = 0, ymax = after_stat(density)),
        geom = "ribbon",
        position = "identity",
        alpha = alpha,
        adjust = 0.5,        # Pas aan als je de gladheid wilt veranderen
        trim = TRUE,
        bw = "nrd0",         # Bandwidth selector
        kernel = "gaussian", # Kernel voor dichtheidschatting
        na.rm = TRUE
      ) +
    #geom_ribbon(aes(ymin = 0, ymax = y), alpha = alpha) +
    geom_segment(data = percentile_df[nrow(percentile_df):1, ], aes(x = percentile_value, xend = percentile_value,
                                           y = 0, yend = Inf, linetype = percentile_label, 
                                           color = percentile_label),
                                           size = 1) +
    ggtitle(title) +
    xlab(xlabel) + 
    ylab("Density") + 
    theme(axis.text = element_text(size = 11),
          axis.title.x = element_text(size = 12, margin = margin(t = 10)),
          axis.title.y = element_text(size = 12, margin = margin(r = 10)),
          plot.title = element_text(size = 12)) +
    scale_fill_brewer(palette = "Set2", name = "Distribution under:") +
    scale_linetype_manual(values = setNames(c("dotted", rep(c("dotted", "dotdash", 
                                                              "dashed", "twodash", "longdash"), 
                                                            length.out = length(percentiles)),
                                              "solid"), 
                                            percentile_df$percentile_label), 
                          name = "") +
    scale_color_manual(values = setNames(percentile_df$color, percentile_df$percentile_label),
                       name = "") +
    theme(legend.key = element_rect(fill = "white")) +
    labs(fill = "Distribution under:", linetype = "Legend", color = "Legend") +
    guides(fill = guide_legend(order = 1),
           linetype = guide_legend(order = 2),
           color = guide_legend(order = 2),
           size = "none") +
    theme(
      panel.background = element_rect(fill = "gray97", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),   
      panel.grid.major = element_line(color = "gray90"),            
      panel.grid.minor = element_line(color = "gray90")             
    ) 
  # + geom_segment(data = sample_value_df,
  #                    aes(x = percentile_value, xend = percentile_value,
  #                        y = 0, yend = Inf, linetype = percentile_label, 
  #                        color = percentile_label),
  #                    size = 1)
  
  if (log_scale) {
    p <- p + scale_x_log10() + 
      labs(caption = "Note: The x-axis is on a log10 scale. The data values themselves are not transformed.")

  }
  
  if (distr_grid) {
    p <- p + facet_grid(. ~ Group_pop_values, scales = "free_x")
  }
  
  if (!is.null(x_lim) && length(x_lim) == 2) {
    p <- p + coord_cartesian(xlim = x_lim)
   } 
  #  else {
  #   if (log_scale) {
  #     value <- log10(df_subset$Value)
  #   } else {
  #     value <- df_subset$Value
  #   }
  #   iqr <- IQR(value, na.rm = TRUE)
  #   q1 <- quantile(value, 0.05, na.rm = TRUE)
  #   q3 <- quantile(value, 0.95, na.rm = TRUE)
  #   lower_limit <- q1 - 1.5 * iqr
  #   upper_limit <- q3 + 1.5 * iqr
  #   if (log_scale) {
  #     lower_limit <- max(lower_limit, min(value[value > 0], na.rm = TRUE)) 
  #   }
  #   p <- p + coord_cartesian(xlim = c(lower_limit, upper_limit))
  # }
  
  return(p)
}


# 
plot_all_groups <- function(plot_df, groups, title, xlabel, x_lim = NULL, 
                            alpha = 0.5, distr_grid = FALSE, 
                            percentiles = NULL, log_scale = FALSE) {
  plot_list <- list()
  for (group in groups) {
    plot <- create_density_plot(plot_df, group, title, xlabel, x_lim, alpha, 
                                distr_grid, percentiles, log_scale)
    plot_list[[group]] <- plot
  }
  return(plot_list)
}

# avoid printing the plot 
print.benchmark_plot <- function(x, ...) {
  plot_list <- x$plots
  ncol <- x$ncol
  nrow <- x$nrow
  
  # Arrange the plots and draw them
  grid::grid.newpage()
  do.call(gridExtra::grid.arrange, c(plot_list, ncol = ncol, nrow = nrow))
  #grid::grid.draw(x)
}


