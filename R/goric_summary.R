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
    cat("generalized order-restricted information criterion: \n")
  } else if (type == "gorica") {
    cat("generalized order-restricted information criterion approximation:\n")
  } else if (type == "goricc") {
    cat("small sample generalized order-restricted information criterion:\n")
  } else if (type == "goricac") {
    cat("small sample generalized order-restricted information criterion approximation:\n")
  }
  
  wt_bar <- sapply(x$objectList, function(x) attr(x$wt.bar, "method") == "boot")
  # we need a check if the hypothesis is equalities only
  ceq_only <- sapply(x$objectList, function(x) nrow(x$constraints) == x$neq)
  wt_bar <- as.logical(wt_bar * !ceq_only)
  
  if (sum(wt_bar) > 0) {
    wt_method_boot <- x$objectList[wt_bar]
    wt_bootstrap_draws  <- sapply(wt_method_boot, function(x) attr(x$wt.bar, "mix.bootstrap"))
    wt_bootstrap_errors <- lapply(wt_method_boot, function(x) attr(x$wt.bar, "error.idx"))
    max_nchar <- max(nchar(names(wt_method_boot)))
    
    len <- length(wt_method_boot)
    if (len > 0) { 
      cat("\n")
      cat("Level probabilities:\n")
      cat("  Number of requested bootstrap draws", wt_bootstrap_draws[1], "\n")
      for (i in 1:len) {
        #cat("Number of successful bootstrap draws for", names(wt_method_boot)[i], ":", (wt_bootstrap_draws[1] - length(wt_bootstrap_errors[[i]])), "\n")
        cat(paste0("  Number of successful bootstrap draws for ", sprintf(paste0("%", max_nchar, "s"), names(wt_method_boot)[i]), ": ", (wt_bootstrap_draws[1] - length(wt_bootstrap_errors[[i]])), "\n"))
      }
    }
  }
  
  cat("\nResults:\n")  
  print(format(df, digits = digits, scientific = FALSE), 
        print.gap = 2, quote = FALSE, right = TRUE)
  cat("---\n")
  
  # Hoe groot mag het verschil zijn. Gebruik ratio en refeer naar par. estimates
  if (any(combn(as.numeric(df$loglik), 2, FUN = function(x) abs(diff(x))) <= 1e-5)) { 
    cat("Note: If log-likelihood values are equal or close to each other, the goric ratio weights are determined", 
        "only by the difference in penalty values. Please check the ratio penalty-weights.\n\n")
  }
  if (comparison == "complement") {
    cat("The order-restricted hypothesis", sQuote(objectnames[1]), "has", 
        sprintf("%.3f", as.numeric(ratio.gw[1,2])), "times more support than its complement.\n\n")
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
            print.gap = 2, quote = FALSE, right = FALSE) 
    } else {
      print(format(ratio.gw, digits = digits, scientific = FALSE, trim = TRUE), 
            print.gap = 5, quote = FALSE, right = FALSE)
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
            print.gap = 2, quote = FALSE, right = FALSE)
    } else {
      print(format(ratio.lw, digits = digits, scientific = FALSE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = FALSE)
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
            print.gap = 2, quote = FALSE, right = FALSE)
    } else {
      print(format(ratio.pw, digits = digits, scientific = FALSE, trim = TRUE), 
            print.gap = 2, quote = FALSE, right = FALSE)
    }
    cat("---\n")
  }
  
  if (!brief) {
    cat("\n\nOrder-restricted coefficients:\n")
    coefs <- trimws(apply(x$ormle$b.restr, 2, sprintf, fmt = dig))
    coefs[coefs == "NA"] <- ""
    rownames(coefs) <- rownames(x$ormle$b.restr)
    print(format(coefs, digits = digits, scientific = FALSE), print.gap = 2,
          quote = FALSE)
    cat("---\n")
    
    vnames <- names(x$ormle$b.restr)
    vnames_len <- length(x$objectList)
    first_na <- apply(x$ormle$b.restr[1:vnames_len, ], 1, function(x) { which(is.na(x))[1] }) -1
    first_na[is.na(first_na)] <- 0
    
    selected_names <- list()
    for (i in 1:length(first_na)) {
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
      rownames(out.rest) <- paste(1:nrow(out.rest), ":", sep = "")
      colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs", "active")
      idx <- ncol(out.rest)
      out.rest[, idx] <- "no"
      out.rest[iact, idx] <- "yes"
      if (nrow(Amat) == meq) {
        out.rest[1:nrow(Amat), idx] <- "yes"
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
    
    invisible(x)
  } else {
    if (!is.null(object$hypotheses_usr)) {
      cat("\norder-restricted hypotheses:\n\n")
      hypotheses_usr <- object$hypotheses_usr
      #hypotheses_usr <- lapply(hypotheses_usr, function(x) unlist(strsplit(x, "\n")))
      # remove user specified paramters
      #hypotheses_usr <- lapply(hypotheses_usr, function(x) x[!grepl(":=", x)])
      
      for (i in 1:length(x$objectList)) {
        text <- gsub("(\\n\\s+)+", "\n", hypotheses_usr[[i]])
        cat(paste0(objectnames[i],":\n", trimws(gsub("\\h+", " ", text, perl = TRUE))), "\n\n")
      }
      
      #calculate max length of vectors
      #max_length <- max(sapply(hypotheses_usr, function(x) length(x)))
      # for (j in 1:length(hypotheses_usr)) {
      #   length(hypotheses_usr[[j]]) <- max_length
      # }
      # hypotheses_usr <- do.call(rbind, hypotheses_usr)
      # row.names(hypotheses_usr) <- paste0(objectnames[1:length(x$objectList)], ":")
      # hypotheses_usr[is.na(hypotheses_usr)] <- ""
      # name.width <- max(sapply(hypotheses_usr, nchar))
      # hypotheses_usr <- format(hypotheses_usr, width = name.width, justify = "left")
      # hypotheses_usr <- as.data.frame(hypotheses_usr)
      # names(hypotheses_usr) <- NULL
      # print(hypotheses_usr)
    }
  }
  cat("\n")
  message(x$messages$mix_weights)
}