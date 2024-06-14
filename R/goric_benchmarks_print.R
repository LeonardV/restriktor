print.benchmark <- function(x, 
                            output_type = c("all", "gw", "rgw", "rlw", "ld"), 
                            ...) {

  # Ensure the object is of class 'benchmark_means'
  if (!inherits(x, "benchmark")) {
    stop("Invalid object. The object should be of class 'benchmark'.", call. = FALSE)
  }
  
  goric_type <- x$type
  pref_hypo <- x$pref_hypo_name
  ratio_pop_means <- x$ratio_pop_means
  group_size <- x$group_size
  ngroups<- x$ngroups
  
  # ANSI escape codes for colors
  blue <- "\033[34m"
  green <- "\033[32m"
  reset <- "\033[0m"
  # Helper function to print rounded data frames
  print_rounded <- function(df, pop_es) {
    pop_es_value <- gsub("pop_es = ", "", pop_es)
    cat(sprintf("Population effect size = %s%s%s\n", green, pop_es_value, reset))
    
    rounded_column <- sprintf("%.3f", df)
    rounded_df <- `dim<-`(rounded_column, dim(df))
    
    rownames(rounded_df) <- rownames(df)
    colnames(rounded_df) <- colnames(df)
    
    print(rounded_df, row.names = TRUE, quote = FALSE)
    cat("\n")
  }
  
  text_gw <- paste0("Benchmark Analysis: Goric Weights Percentiles for Preferred Hypothesis '", pref_hypo, "'")
  text_rgw <- paste0("Benchmark Analysis: Ratio of Goric Weights Percentiles for Preferred Hypothesis '", pref_hypo, "'")
  text_rlw <- paste0("Benchmark Analysis: Ratio of Likelihood Weights Percentiles for Preferred Hypothesis '", pref_hypo, "'")
  text_ld <- paste0("Benchmark Analysis: Difference in Likelihood Values Percentiles for Preferred Hypothesis '", pref_hypo, "'")
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat(paste0(blue, "Benchmark Results", reset), "\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("Preferred Hypothesis: %s%s%s\n", green, pref_hypo, reset))
  cat(sprintf("Number of Groups: %s%s%s\n", green, ngroups, reset))
  cat(sprintf("Group Size: %s%s%s\n", green, paste(group_size, collapse = ", "), reset))
  cat(sprintf("Ratio of Population Means: %s%s%s\n", green, paste(sprintf("%.3f", ratio_pop_means), collapse = ", "), reset))
  #cat(strrep("-", 70), "\n")
  cat("\n")
  
  # Function to print section headers and footers
  print_section <- function(header, content_printer, nchar) {
    cat("\n")
    cat(strrep("=", nchar), "\n")
    cat(paste0(blue, header, reset), "\n")
    cat(strrep("-", nchar), "\n")
    content_printer()
    #cat(strrep("-", nchar), "\n")
    cat("\n")
  }
  
  # Normalize output_type to lowercase
  output_type <- tolower(output_type)
  
  # Print benchmarks for Goric-weights percentiles
  if ("all" %in% output_type || "gw" %in% output_type) {
    print_section(
      text_gw,
      function() {
        for (pop_es in names(x$benchmarks_goric_weights)) {
          print_rounded(x$benchmarks_goric_weights[[pop_es]], pop_es)
        }
      }, nchar(text_gw)
    )
  }
  
  # Print benchmarks for ratio Goric-weights percentiles
  if ("all" %in% output_type || "rgw" %in% output_type) {
    print_section(
      text_rgw,
      function() {
        for (pop_es in names(x$benchmarks_ratio_goric_weights)) {
          print_rounded(x$benchmarks_ratio_goric_weights[[pop_es]], pop_es)
        }
      }, nchar(text_rgw)
    )
  }
  
  # Print benchmarks for ratio likelihood-weights percentiles
  if ("all" %in% output_type || "rlw" %in% output_type) {
    print_section(
      text_rlw,
      function() {
        for (pop_es in names(x$benchmarks_ratio_ll_weights)) {
          print_rounded(x$benchmarks_ratio_ll_weights[[pop_es]], pop_es)
        }
      }, nchar(text_rlw)
    )
  }
  
  # Print benchmarks for difference likelihood-values percentiles
  if ("all" %in% output_type || "ld" %in% output_type) {
    print_section(
      text_ld,
      function() {
        for (pop_es in names(x$benchmarks_difLL)) {
          print_rounded(x$benchmarks_difLL[[pop_es]], pop_es)
        }
      }, nchar(text_ld)
    )
  }
}

