summary.evSyn <- function(object, ...) {
  x <- object
  class(x) <- NULL
  
  ans <- list(
    type = x$type,
    type_ev = x$type_ev,
    #n_studies = nrow(x$GORICA_weight_m[,, drop = FALSE]),
    n_studies = x$n_studies, 
    hypotheses = x$hypotheses,
    GORICA_weight_m = x$GORICA_weight_m,
    GORICA_m = x$GORICA_m,
    LL_weights_m = x$LL_weights_m,
    Cumulative_LL_weights = x$Cumulative_LL_weights,
    Cumulative_LL = x$Cumulative_LL,
    LL_m = x$LL_m,
    PT_m = x$PT_m,
    Cumulative_GORICA_weights = x$Cumulative_GORICA_weights,
    Cumulative_GORICA = x$Cumulative_GORICA
  )
  
  if (!is.null(x$PT_m)) {
    sequence    <- paste0("Studies 1-", 1:ans$n_studies)
    sequence[1] <- "Studies 1"
    Cumulative_PT <- apply(x$PT_m[,, drop = FALSE], 2, cumsum)  
    Cumulative_PT <- matrix(Cumulative_PT, nrow = nrow(x$PT_m), 
                            dimnames = list(sequence, colnames(x$PT_m)))
    if (x$type_ev %in% c("equal", "average")) {
      Cumulative_PT <- Cumulative_PT / seq_len(ans$n_studies)
    } 
    ans$Cumulative_PT <- Cumulative_PT
  }
  
  final <- c()
  
  if (!is.null(x[["Cumulative_GORICA_weights"]])) {
    fcgw <- t(x[["Cumulative_GORICA_weights"]]["Final", ])
    if (x$type == "gorica") {
      rownames(fcgw) <- "GORICA weights"
    } else if (x$type == "goricac") {
      rownames(fcgw) <- "GORICAC weights"
    }
    colnames(fcgw)  <- colnames(x[["Cumulative_GORICA_weights"]])
    final <- rbind(final, fcgw)
  }
  
  if (!is.null(x[["Cumulative_GORICA"]])) {
    fcgv <- t(x[["Cumulative_GORICA"]]["Final", ])
    if (x$type == "gorica") {
      rownames(fcgv) <- "GORICA values"
    } else if (x$type == "goricac") {
      rownames(fcgv) <- "GORICAC values"
    }
    final <- rbind(final, fcgv)
  }
  
  if (!is.null(x[["Cumulative_LL_weights"]])) {
    fcgw <- t(x[["Cumulative_LL_weights"]]["Final", ])
    rownames(fcgw) <- "Log-likelihood weights"
    colnames(fcgw)  <- colnames(x[["Cumulative_LL_weights"]])
    final <- rbind(final, fcgw)
  }
  
  if (!is.null(x[["Cumulative_LL"]])) {
    fcllv <- t(x[["Cumulative_LL"]][ans$n_studies, ])
    rownames(fcllv) <- "Log-likelihood values"
    final <- rbind(final, fcllv)
  }
  
  if (!is.null(ans$Cumulative_PT)) { 
    fcptv <- t(ans$Cumulative_PT[ans$n_studies, ])
    rownames(fcptv) <- "Penalty term values"
    final <- rbind(final, fcptv)
  }
  
  ans$Final_Cumulative_results <- final
  
  ans$Ratio_GORICA_weight_mu <- x$ratio_GORICA_weight_mu
  ans$Ratio_GORICA_weight_mc <- x$ratio_GORICA_weight_mc
  ans$Final_ratio_GORICA_weights <- x$Final_ratio_GORICA_weights
  ans$Final_ratio_LL_weights <- x$Final_ratio_LL_weights
  
  ans$messages <- list(mix_weights = x$messages$mix_weights)
  class(ans) <- "summary.evSyn"
  
  ans
  
  # TO DO x$messageAdded

}



## print.summary function
print.summary.evSyn <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  # added or equal approach
  type_ev <- x$type_ev
  # number of studies
  S <- x$n_studies
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  cat(paste(type_ev, "Evidence Synthesis results:\n"), sep = "")
  
  indentation <- "    "  # Four spaces for indentation
  
  cat("\nStudy-specific results:\n")
  
  if (!is.null(x$GORICA_weight_m)) {
    if (x$type == "gorica") {
      cat("\n    GORICA weights:\n")
    } else if (x$type == "goricac") {
      cat("\n    GORICAC weights:\n") 
    }
    formatted_gw <- apply(x$GORICA_weight_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_gw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$GORICA_m)) {
    if (x$type == "gorica") {
      cat("\n    GORICA values:\n")
    } else if (x$type == "goricac") {
      cat("\n    GORICAC values:\n")
    }
    formatted_gv <- apply(x$GORICA_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_gv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$LL_weights_m)) {
    cat("\n    Log-likelihood weights:\n")  
    formatted_llw <- apply(x$LL_weights_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_llw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$LL_m)) {
    cat("\n    Log-likelihood values:\n")  
    formatted_llv <- apply(x$LL_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_llv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$PT_m)) {
    cat("\n    Penalty term values:\n")  
    formatted_ptv <- apply(x$PT_m[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_ptv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  cat("\nCumulative results:\n")
  
  if (!is.null(x[["Cumulative_GORICA_weights"]])) {
    if (x$type == "gorica") {
      cat("\n    GORICA weights:\n")
    } else if (x$type == "goricac") {
      cat("\n    GORICAC weights:\n") 
    } 
    formatted_cgw <- apply(x$Cumulative_GORICA_weights[1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cgw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["Cumulative_GORICA"]])) {
    if (x$type == "gorica") {
      cat("\n    GORICA values:\n")
    } else if (x$type == "goricac") {
      cat("\n    GORICAC values:\n")
    }  
    formatted_cgv <- apply(x$Cumulative_GORICA[1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cgv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["Cumulative_LL_weights"]])) {
    cat("\n    Log-likelihood weights:\n")  
    formatted_cllw <- apply(x$Cumulative_LL_weights[1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cllw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$LL_m)) {
    cat("\n    Log-likelihood values:\n")  
    formatted_cllv <- apply(x$Cumulative_LL[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cllv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$PT_m)) {
    cat("\n    Penalty term values:\n")  
    formatted_cptv <- apply(x$Cumulative_PT[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cptv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  
  cat("\nFinal results:\n")
  formatted_final <- apply(x$Final_Cumulative_results[,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
  captured_output <- capture.output(print(formatted_final, row.names = TRUE, right = TRUE, quote = "FALSE"))
  adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
  cat(paste0(adjusted_output, "\n"), sep = "")
  
  
  cat("\nFinal ratios:\n")
  
  if (!is.null(x$Final_ratio_GORICA_weights)) {
    if (x$type == "gorica") {
      cat("\n    GORICA weights:\n")
    } else if (x$type == "goricac") {
      cat("\n    GORICAC weights:\n") 
    } 
    formatted_frgw <- apply(x$Final_ratio_GORICA_weights, c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_frgw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x$Final_ratio_LL_weights)) {
    cat("\n    Log-likelihood weights:\n")  
    formatted_frllw <- apply(x$Final_ratio_LL_weights, c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_frllw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  cat("\n")
  message(x$messages$mix_weights)
  
  # # TO DO
  # if (!is.null(x$messageAdded)) {
  #   cat("\n")
  #   message(x$messageAdded)
  # }
  
}