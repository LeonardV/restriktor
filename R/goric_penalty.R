## compute asymptotic and small sample penalty term value for the goric(a)
penalty_goric <- function(Amat, meq, LP, correction = FALSE, 
                          sample.nobs = NULL, ...) {

  num_cols <- ncol(Amat)
  method <- attr(LP, "method")
  
  compute_PT <- function(lPT_values) {
    if (correction) {
      N <- sample.nobs
      return(sum(((N * (lPT_values + 1) / (N - lPT_values - 2))) * LP))
    } else {
      return(1 + sum(lPT_values * LP))
    }
  }
  
  if (all(Amat == 0)) { 
    lPT_values <- ifelse(correction, num_cols, 0:num_cols)
    return(compute_PT(lPT_values))
  }
  
  switch(method,
         boot = {
           lPT_values <- 0:num_cols
           return(compute_PT(lPT_values))
         },
         pmvnorm = {
           min_col <- num_cols - nrow(Amat)
           max_col <- num_cols - meq
           lPT_values <- min_col:max_col
           return(compute_PT(lPT_values))
         },
         stop("Unknown method specified in LP attribute.")
  )
  
  # if (correction) {
  #   N <- sample.nobs
  #   # unconstrained case
  #   if (all(c(Amat) == 0)) {
  #     lPT <- ncol(Amat)
  #     PT  <- ( (N * (lPT + 1) / (N - lPT - 2)) )
  #   } else {
  #     if (attr(LP, "method") == "boot") {
  #       lPT <- 0 : ncol(Amat)
  #       PT  <- sum( ( (N * (lPT + 1) / (N - lPT - 2) ) ) * LP)
  #     } else if (attr(LP, "method") == "pmvnorm") {
  #       min.col <- ncol(Amat) - nrow(Amat) # p - q1 - q2
  #       max.col <- ncol(Amat) - meq        # p - q2
  #       lPT     <- min.col : max.col
  #       PT      <- sum( ( (N * (lPT + 1) / (N - lPT - 2) ) ) * LP)
  #     }
  #   }
  # } else {
  #   if (all(c(Amat) == 0)) {
  #     PT <- 1 + ncol(Amat)
  #   } else {
  #     if (attr(LP, "method") == "boot") {
  #       PT <- 1 + sum(0 : ncol(Amat) * LP)
  #     } else if (attr(LP, "method") == "pmvnorm") {
  #       min.col <- ncol(Amat) - nrow(Amat) # p - q1 - q2
  #       max.col <- ncol(Amat) - meq        # p - q2
  #       PT <- 1 + sum(min.col : max.col * LP)
  #     }
  #   }
  # }

  #return(PT)
}


# [CHANGE 2026-10 | audit] E13: TO DO on the mlm penalty for the residual (co)variance; goricc/goricac blocked for mlm
# TO DO mlm: the 1 in the penalty terms below is the penalty for the residual
#       variance. For mlm objects ny*(ny+1)/2 (co)variances are estimated. A 
#       choice still has to be made whether to use 1 or ny*(ny+1)/2.
#       The small-sample correction (goricc/goricac, correction = TRUE) for
#       mlm objects is blocked in goric.lm() until it has been derived for a
#       multivariate residual covariance matrix.
# [/CHANGE 2026-10]
penalty_complement_goric <- function(Amat, meq, type, wt.bar, 
                                     sample.nobs = NULL, debug = FALSE) {
  # compute the number of free parameters f in the complement
  p <- ncol(Amat)
  # rank q1
  #lq1 <- qr(Amat.ciq)$rank
  lq1 <- qr(Amat)$rank - meq
  
  # compute penalty term value PTc
  if (type %in% c("goric", "gorica")) {
    idx <- length(wt.bar)
    if (attr(wt.bar, "method") == "boot") {
      PTc <- as.numeric(1 + p - wt.bar[idx-meq] * lq1)  
    } else if (attr(wt.bar, "method") == "pmvnorm") {
      PTc <- as.numeric(1 + p - wt.bar[idx] * lq1) 
    } else {
      stop("restriktor ERROR: no level probabilities (chi-bar-square weights) found.")
    }
  } else if (type %in% c("goricc", "goricac")) {
    idx <- length(wt.bar) 
    if (is.null(sample.nobs)) {
      stop("restriktor ERROR: the argument sample.nobs is not found.")
    }
    N <- sample.nobs
    # small sample correction
    if (attr(wt.bar, "method") == "boot") {
      PTc <- 1 + wt.bar[idx-meq] * (N * (p - lq1) + (p - lq1) + 2) / (N - (p - lq1) - 2) + 
        (1 - wt.bar[idx-meq]) * (N * p + p + 2) / (N - p - 2)
    } else if (attr(wt.bar, "method") == "pmvnorm") {
      PTc <- 1 + wt.bar[idx] * (N * (p - lq1) + (p - lq1) + 2) / (N - (p - lq1) - 2) + 
        (1 - wt.bar[idx]) * (N * p + p + 2) / (N - p - 2) 
    } 
  }
  
  if (debug) {
    cat("penalty term value =", PTc, "\n")
  }
  
  # correction for gorica(c): -1 because no 
  if (type %in% c("gorica", "goricac")) {
    PTc <- PTc - 1
  }
  
  return(PTc)
}
