## plot function
plot.evSyn <- function(x, output_type = "gorica_weights", ...) {
  
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
  
  # Zet de data en labels voor de plot op basis van geselecteerde data_list
  weight_m <- data_list$weight_m
  cumulative_weights <- data_list$cumulative_weights
  y_label <- data_list$y_label
  
  # Bepaal de namen van hypothesen en studies
  namesH <- colnames(weight_m)
  NrHypos_incl <- ncol(weight_m[,,drop = FALSE])
  S <- nrow(weight_m[,,drop = FALSE])
  Name_studies <- as.factor(1:S)
  
  # Creëer data frames voor individuele en cumulatieve gewichten
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
  
  # Combineer de data frames
  plot_data <- rbind(per_study_df, cumulative_df)
  plot_data$variable <- rep(rep(namesH, each = S), times = times)
  
  # Plot aanmaken
  ggplot(plot_data, aes(x = .data[['study']],  
                        y = .data[['weight']], 
                        shape = .data[['weight_type']], 
                        linetype = .data[['weight_type']], 
                        color = .data[['variable']])) +
    geom_point(size = 3) +
    geom_line(data = plot_data[plot_data[['weight_type']] == "cumulative", ], 
              aes(group = .data[['variable']]), linewidth = 1) +
    scale_color_brewer(palette = "Dark2") +
    theme(
      plot.margin = unit(c(1,1,1,1), "cm"),
      legend.position = "bottom",
      legend.margin = margin(t = 10, r = 0, b = 3, l = 0),
      legend.key = element_blank(),
      legend.text = element_text(size = 12),
      axis.text.x  = element_text(size = 12), axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, vjust = -3), 
      axis.title.y = element_text(size = 14, vjust = 5),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    scale_y_continuous(limits = c(0, 1), n.breaks = 10) +
    scale_x_discrete(expand = c(0, 0.05)) +
    labs(x = "Studies", y = y_label,
         title = paste("Cumulative", y_label, "and", y_label, "per study"),
         shape = "", color = "", linetype = "") +
    guides(shape = guide_legend(override.aes = list(linetype = 0)))
}
