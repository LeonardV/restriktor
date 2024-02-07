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