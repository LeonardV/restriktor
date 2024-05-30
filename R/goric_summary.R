summary.con_goric <- function(object, brief = TRUE, 
                              digits = max(3, getOption("digits") - 4), ...) {
  
  x <- object
  type <- x$type
  comparison <- x$comparison
  dig <- paste0("%6.", digits, "f")
  
  ratio.gw <- x$ratio.gw
  rownames(ratio.gw) <- rownames(x$ratio.gw)
  
  x2 <- lapply(x$result[-1], sprintf, fmt = dig)
  df <- data.frame(model = x$result$model, x2)
  
  objectnames <- as.character(df$model)
  
  Amat <- x$constraints
  meq  <- x$neq
  bvec <- x$rhs
  iact <- lapply(x$objectList, FUN = function(x) { x$iact } )
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  
  if (type == "goric") {
    cat("generalized order-restricted information criterion: \n\n")
  } else if (type == "gorica") {
    cat("generalized order-restricted information criterion approximation:\n\n")
  } else if (type == "goricc") {
    cat("small sample generalized order-restricted information criterion:\n\n")
  } else if (type == "goricac") {
    cat("small sample generalized order-restricted information criterion approximation:\n\n")
  }
  
  # wt_bar <- sapply(x$objectList, function(x) attr(x$wt.bar, "method") == "boot")
  # # we need a check if the hypothesis is equalities only
  # ceq_only <- sapply(x$objectList, function(x) nrow(x$constraints) == x$neq)
  # wt_bar <- as.logical(wt_bar * !ceq_only)
  # 
  # if (sum(wt_bar) > 0) {
  #   wt_method_boot <- x$objectList[wt_bar]
  #   wt_bootstrap_draws  <- sapply(wt_method_boot, function(x) attr(x$wt.bar, "total_bootstrap_draws"))
  #   wt_bootstrap_errors <- lapply(wt_method_boot, function(x) attr(x$wt.bar, "error.idx"))
  #   max_nchar <- max(nchar(names(wt_method_boot)))
  #   
  #   len <- length(wt_method_boot)
  #   if (len > 0) { 
  #     cat("\n")
  #     cat("Bootstrap-based penalty term calculation:\n")
  #     cat("  Number of bootstrap draws", wt_bootstrap_draws[1], "\n")
  #     for (i in 1:len) {
  #       cat(paste0("  Number of successful bootstrap draws for ", sprintf(paste0("%", max_nchar, "s"), names(wt_method_boot)[i]), ": ", (wt_bootstrap_draws[1] - length(wt_bootstrap_errors[[i]])), "\n"))
  #     }
  #   }
  # }
  
  wt_bar_attributes <- lapply(x$objectList, function(obj) {
    list(
      method = attr(obj$wt.bar, "method"),
      converged = attr(obj$wt.bar, "converged"),
      total_draws = attr(obj$wt.bar, "total_bootstrap_draws"),
      errors = attr(obj$wt.bar, "error.idx")
    )
  })
  
  # Compute indicators
  wt_bar <- vapply(wt_bar_attributes, function(attr) attr$method == "boot", logical(1))
  ceq_only <- vapply(x$objectList, function(obj) nrow(obj$constraints) == obj$neq, logical(1))
  wt_bar <- wt_bar & !ceq_only
  
  if (any(wt_bar)) {
    wt_bar_attributes <- wt_bar_attributes[wt_bar]
    wt_method_boot <- x$objectList[wt_bar]
    max_nchar <- max(nchar(names(wt_method_boot)))
    
    # Summarize bootstrap information
    bootstrap_summary <- vapply(wt_bar_attributes, function(attr) {
      successful_draws <- attr$total_draws - length(attr$errors)
      paste0(successful_draws, ifelse(attr$converged, " (Converged)", " (Not converged)"))
    }, character(1))
    
    converged <- vapply(wt_bar_attributes, function(attr) attr$converged, logical(1))
    total_bootstrap_draws <- vapply(wt_bar_attributes, function(attr) attr$total_draws, integer(1))
    wt_bootstrap_errors <- sapply(wt_bar_attributes, function(attr) attr$errors)
    
    if (length(wt_method_boot) > 0) {
      #cat("\n")
      successful_draws <- total_bootstrap_draws - sapply(wt_bootstrap_errors, length)
      has_errors <- vapply(wt_bootstrap_errors, function(errors) length(errors) > 0, logical(1))
      not_all_converged <- !all(converged)
      not_all_draws_successful <- !all(successful_draws == total_bootstrap_draws)
      
      if (any(has_errors) || not_all_converged) { 
        if (not_all_draws_successful || not_all_converged) {
          cat("Bootstrap-based penalty term calculation:\n")
          cat("  Number of bootstrap draws:", sapply(wt_bar_attributes, `[[`, "total_draws"), "\n")
          for (i in seq_along(bootstrap_summary)) {
            cat(sprintf("  Number of successful bootstrap draws for %*s: %s\n", 
                        max_nchar, names(wt_method_boot)[i], bootstrap_summary[i]))
          }
          cat("\n")
        } 
      }
    }
  }
  
  cat("Results:\n")  
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE, right = TRUE)
  cat("---\n")
  
  if (comparison == "complement") {
    cat("The order-restricted hypothesis", sQuote(objectnames[1]), "has", 
        sprintf("%.2f", as.numeric(ratio.gw[1,2])), "times more support than its complement.\n\n")
  } 
  
  
  if (!is.null(x$ratio.gw)) {
    if (type == "goric") {
      cat("\nRatio GORIC-weights:\n")
    } else if (type == "gorica") {
      cat("\nRatio GORICA-weights:\n")
    } else if (type == "goricc") {
      cat("\nRatio GORICC-weights:\n")
    } else if (type == "goricac") {
      cat("\nRatio GORICAC-weights:\n") 
    }
    
    ratio.gw <- apply(x$ratio.gw, 2, sprintf, fmt = dig)
    rownames(ratio.gw) <- rownames(x$ratio.gw)
    class(ratio.gw) <- "numeric"
    
    if (max(ratio.gw, na.rm = TRUE) >= 1e4) {
      print(format(ratio.gw, digits = digits, scientific = TRUE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = TRUE) 
    } else {
      print(format(ratio.gw, digits = digits, scientific = FALSE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = TRUE)
    }
    cat("---\n")
  }
  
  if (!is.null(x$ratio.lw)) {
    cat("\nRatio loglik-weights:\n")
    ratio.lw <- apply(x$ratio.lw, 2, sprintf, fmt = dig)
    rownames(ratio.lw) <- rownames(x$ratio.lw) 
    class(ratio.lw) <- "numeric"
    
    if (max(ratio.lw, na.rm = TRUE) >= 1e4) {
      print(format(ratio.lw, digits = digits, scientific = TRUE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = TRUE)
    } else {
      print(format(ratio.lw, digits = digits, scientific = FALSE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = TRUE)
    }
    cat("---\n")
  }
  
  if (!is.null(x$ratio.pw)) {
    cat("\nRatio penalty-weights:\n")
    ratio.pw <- apply(x$ratio.pw, 2, sprintf, fmt = dig)
    rownames(ratio.pw) <- rownames(x$ratio.pw)
    class(ratio.pw) <- "numeric"
    
    if (max(ratio.pw, na.rm = TRUE) >= 1e4) {
      print(format(x$ratio.pw, digits = digits, scientific = TRUE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = TRUE)
    } else {
      print(format(ratio.pw, digits = digits, scientific = FALSE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = TRUE)
    }
    cat("---\n")
  }
  
  if (!brief) {
    cat("\n\nOrder-restricted coefficients:\n")
    coefs <- trimws(apply(x$ormle$b.restr, 2, sprintf, fmt = dig))
    coefs[coefs == "NA"] <- ""
    rownames(coefs) <- rownames(x$ormle$b.restr)
    # print(format(coefs, digits = digits, scientific = TRUE, trim = TRUE), 
    #       print.gap = 2, quote = FALSE, right = TRUE) 
    # 
    print(coefs, scientific = TRUE, right = TRUE, quote = FALSE, print.gap = 2)
    cat("---\n")
    
    vnames <- names(x$ormle$b.restr)
    vnames_len <- length(x$objectList)
    first_na <- apply(x$ormle$b.restr[1:vnames_len, , drop = FALSE], 1, function(x) { which(is.na(x))[1] }) -1 
    first_na[is.na(first_na)] <- 0
    
    selected_names <- list()
    for (i in seq_len(length(first_na))) {
      if (first_na[i] == 0) {
        selected_names[[i]] <- vnames
      } else {
        selected_names[[i]] <- vnames[1:first_na[i]]
      }
    }
    
    fn <- function(Amat, bvec, meq, iact, vnames) {
      colnames(Amat) <- vnames
      out.rest <- cbind(round(Amat, 4), c(rep("   ==", meq), rep("   >=", nrow(Amat) - 
                                                                   meq)), bvec, " ")
      rownames(out.rest) <- paste(seq_len(nrow(out.rest)), ":", sep = "")
      colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs", "active")
      idx <- ncol(out.rest)
      out.rest[, idx] <- "no"
      out.rest[iact, idx] <- "yes"
      if (nrow(Amat) == meq) {
        out.rest[seq_len(nrow(Amat)), idx] <- "yes"
      }  
      out.rest <- as.data.frame(out.rest)
      
      out.rest
    }
    
    conMat <- list()
    for (i in 1:vnames_len) {
      conMat[[i]] <- fn(Amat = Amat[[i]], bvec = bvec[[i]], meq = meq[[i]], 
                        iact = iact[[i]], vnames = selected_names[[i]])  
    }
    names(conMat) <- x$objectNames
    
    if (comparison == "complement") {
      conMat$complement <- paste("not", x$objectNames) 
    }
    
    cat("\nRestriction matrices:\n")
    print(conMat, quote = FALSE, scientific = FALSE)
    
    #invisible(x)
  } else {
    if (!is.null(object$hypotheses_usr)) {
      cat("\norder-restricted hypotheses:\n\n")
      hypotheses_usr <- object$hypotheses_usr
      for (i in seq_len(length(x$objectList))) {
        text <- gsub("(\\n\\s+)+", "\n", hypotheses_usr[[i]])
        cat(paste0(objectnames[i],":\n", trimws(gsub("\\h+", " ", text, perl = TRUE))), "\n\n")
      }
    }
  }
  #cat("\n")
  #message(x$messages$mix_weights)
}
