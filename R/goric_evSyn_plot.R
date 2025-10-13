plot.evSyn <- function(x, output_type = "gorica_weights", 
                       xlab = NULL, ...) {
  if (!output_type %in% c("gorica_weights", "ll_weights")) {
    stop("Restriktor ERROR: output_type must be gorica_weights or ll_weights", call. = FALSE)
  }
  
  if (inherits(x, "evSyn_est")) {
    data_list <- switch(output_type,
                        "gorica_weights" = list(weight_m = x$GORICA_weight_m, 
                                                cumulative_weights = x$Cumulative_GORICA_weights,
                                                y_label = "GORIC(A) weights"),
                        "ll_weights" = list(weight_m = x$LL_weights_m, 
                                            cumulative_weights = x$Cumulative_LL_weights,
                                            y_label = "Log-Likelihood weights")
    )
  } else if (inherits(x, "evSyn_LL")) {
    data_list <- switch(output_type,
                        "gorica_weights" = list(weight_m = x$GORICA_weight_m, 
                                                cumulative_weights = x$Cumulative_GORICA_weights,
                                                y_label = "GORIC(A) weights"),
                        "ll_weights" = list(weight_m = x$LL_weights_m, 
                                            cumulative_weights = x$Cumulative_LL_weights,
                                            y_label = "Log-Likelihood weights")
    )
  } else if (inherits(x, "evSyn_ICvalues")) {
    data_list <- switch(output_type,
                        "gorica_weights" = list(weight_m = x$GORICA_weight_m, 
                                                cumulative_weights = x$Cumulative_GORICA_weights,
                                                y_label = "GORIC(A) weights")
    )
  } else if (inherits(x, "evSyn_ICweights")) {
    data_list <- switch(output_type,
                        "gorica_weights" = list(weight_m = x$GORICA_weight_m, 
                                                cumulative_weights = x$Cumulative_GORICA_weights,
                                                y_label = "GORIC(A) weights"))
  } 
  
  # Base data and labels for the plot on selected data_list
  weight_m <- data_list$weight_m
  cumulative_weights <- data_list$cumulative_weights
  y_label <- data_list$y_label
  
  # Determine names of hypotheses and studies
  namesH <- colnames(weight_m)
  NrHypos_incl <- ncol(weight_m[,,drop = FALSE])
  S <- nrow(weight_m[,,drop = FALSE])
  if (!is.null(xlab)) {
    if (length(xlab) != S) {
      stop("Restriktor ERROR: Length of xlab must match the number of studies", call. = FALSE)
    }
    Name_studies <- as.factor(xlab) 
    angle_x = 30
    vjust_x = 1 # 0: bottom-aligned to 1: is top-aligned
    hjust_x = 1 # right-aligned
  } else {
    #Name_studies <- as.factor(1:S)
    Name_studies <- factor(x$orderStudies, levels = x$orderStudies)
    angle_x = 0
    vjust_x = 0.5 # centered vertically
    hjust_x = 0.5 # centered horizontally
  }

  # Create data frames for individual and cumulative weights
  if (all(is.na(weight_m))) {
    per_study_df <- NULL
    times <- 1
  } else {
    per_study_df <- data.frame(study = rep(Name_studies, NrHypos_incl),
                               weight = c(weight_m))
    per_study_df$weight_type <- "per study"
    times <- 2
  }
  
  cumulative_df <- data.frame(study = rep(Name_studies, NrHypos_incl),
                              weight = c(cumulative_weights[1:S, ]))
  cumulative_df$weight_type <- "cumulative"
  
  # Combine the data frames
  plot_data <- rbind(per_study_df, cumulative_df)
  plot_data$variable <- rep(rep(namesH, each = S), times = times)
  
  # Plot
  ggplot(plot_data, aes(x = .data[['study']],
                        y = .data[['weight']])) +
    geom_point(size = 3, aes(color = factor(.data[['variable']], levels = namesH), shape = .data[['weight_type']])) + 
    theme(legend.position = "bottom", 
          legend.box = "vertical",
          legend.margin=margin(unit(0, "cm"))) +
    geom_line(data = plot_data[plot_data[['weight_type']] == "cumulative", ],
              aes(group = .data[['variable']], color = .data[['variable']]), linewidth = 1) +
    scale_color_brewer(palette = "Dark2") +
    theme(
      plot.margin = unit(c(1,1,1,1), "cm"),
      legend.text = element_text(size = 12),
      axis.text.x  = element_text(size = 12, 
                                  angle = angle_x, vjust = vjust_x, hjust = hjust_x), 
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, vjust = -3), 
      axis.title.y = element_text(size = 14, vjust = 5),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    scale_y_continuous(limits = c(0, 1), n.breaks = 10) +
    scale_x_discrete(expand = c(0, 0.05)) +
    labs(x = "Studies", y = y_label,
         title = paste("Cumulative", y_label, "and", y_label, "per study"),
         subtitle = paste(x$type_ev, "Evidence Synthesis results"),
         shape = "", color = "") + 
    guides(color = guide_legend(order = 2),
           shape = guide_legend(order = 1)) 
  
}
