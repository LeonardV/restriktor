## compute asymptotic and small sample penalty term value for the goric(a)
penalty_goric <- function(Amat, meq, LP, correction = FALSE, 
                          sample.nobs = NULL, ...) {

  if (correction) { 
    N <- sample.nobs  
    # unconstrained case
    if (all(c(Amat) == 0)) {
      lPT <- ncol(Amat)
      PT  <- ( (N * (lPT + 1) / (N - lPT - 2)) ) 
    } else {
      if (attr(LP, "method") == "boot") { 
        lPT <- 0 : ncol(Amat)
        PT  <- sum( ( (N * (lPT + 1) / (N - lPT - 2) ) ) * LP)
      } else if (attr(LP, "method") == "pmvnorm") {
        min.col <- ncol(Amat) - nrow(Amat) # p - q1 - q2
        max.col <- ncol(Amat) - meq        # p - q2
        lPT     <- min.col : max.col
        PT      <- sum( ( (N * (lPT + 1) / (N - lPT - 2) ) ) * LP) 
      }
    }
  } else {
    if (all(c(Amat) == 0)) {
      PT <- 1 + ncol(Amat)
    } else {
      if (attr(LP, "method") == "boot") {  
        PT <- 1 + sum(0 : ncol(Amat) * LP)  
      } else if (attr(LP, "method") == "pmvnorm") {
        min.col <- ncol(Amat) - nrow(Amat) # p - q1 - q2
        max.col <- ncol(Amat) - meq        # p - q2
        PT <- 1 + sum(min.col : max.col * LP) 
      }
    }
  }
  
  return(PT)
}