## Penalty terms of the GORIC(A)(C) ------------------------------------------
##
## Notation: N = sample size (number of subjects), l = number of free
## (regression) parameters in the estimate vector (averaged over with the
## level probabilities LP_l), d = number of nuisance parameters (parameters of
## the model that are not in the estimate vector, e.g., (co)variances, or
## coefficients that are left out of 'est'), K_l = l + d, and, for a
## multivariate regression (mlm), q = number of response variables (then
## d = q(q+1)/2).
##
##   GORIC  : sum_l LP_l * K_l                         (= d + sum_l LP_l * l)
##   GORICC : sum_l LP_l * N K_l / (N - K_l - 1)
##            (mlm: sum_l LP_l * N K_l / (N - l/q - q - 1); Bedrick and Tsai,
##            1994, Biometrics, 50, 226-231)
##   GORICA : sum_l LP_l * l                           (= GORIC - d)
##   GORICAC: sum_l LP_l * [N K_l / (N - K_l - 1) - d] (= GORICC - d)
##
## For a univariate regression (d = 1), the GORICC term is the usual
## N (l + 1) / (N - l - 2). The GORICA(C) penalties are obtained by
## subtracting d from the GORIC(C) penalties (in the callers): then the
## nuisance parameters do not count asymptotically, but they do affect the
## small sample correction (via K_l).

# small sample (AICc-type) penalty for f free parameters, d nuisance
# parameters, and q responses (q > 1: multivariate regression, mlm)
small_sample_PT <- function(f, N, n_nuis = 1, n_resp = 1L) {
  if (n_resp > 1L) {
    N * (f + n_nuis) / (N - f / n_resp - n_resp - 1)
  } else if (n_nuis == 1) {
    N * (f + 1) / (N - f - 2) # the usual expression (univariate regression)
  } else {
    N * (f + n_nuis) / (N - f - n_nuis - 1)
  }
}

# number of nuisance parameters (d) and responses (q) of a fitted restriktor
# object (lm, mlm, rlm, glm): the (co)variance/dispersion parameters, since
# all regression coefficients are in the estimate vector.
nuisance_info_restriktor <- function(z) {
  if (inherits(z, "conMLM")) {
    q <- NCOL(residuals(z$model.org))
    list(n_nuis = q * (q + 1) / 2, n_resp = q)
  } else if (inherits(z, "conGLM")) {
    has_disp <- z$model.org$family$family %in% c("gaussian", "Gamma", "inverse.gaussian")
    list(n_nuis = as.numeric(has_disp), n_resp = 1L)
  } else {
    list(n_nuis = 1, n_resp = 1L)
  }
}

## compute the penalty term value of the GORIC(C) (see above); the callers
## subtract n_nuis for the GORICA(C).
penalty_goric <- function(Amat, meq, LP, correction = FALSE, 
                          sample.nobs = NULL, n_resp = 1L, n_nuis = NULL, ...) {

  num_cols <- ncol(Amat)
  method <- attr(LP, "method")
  if (is.null(n_nuis)) {
    n_nuis <- n_resp * (n_resp + 1) / 2 # (co)variance parameters
  }
  
  compute_PT <- function(lPT_values) {
    if (correction) {
      N <- sample.nobs
      return(sum(small_sample_PT(lPT_values, N, n_nuis, n_resp) * LP))
    } else {
      return(n_nuis + sum(lPT_values * LP))
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


penalty_complement_goric <- function(Amat, meq, type, wt.bar, 
                                     sample.nobs = NULL, debug = FALSE,
                                     n_resp = 1L, n_nuis = NULL) {
  # n_nuis: number of nuisance parameters (d); n_resp: number of response
  # variables (q > 1 for an mlm); see penalty_goric()
  if (is.null(n_nuis)) {
    n_nuis <- n_resp * (n_resp + 1) / 2
  }
  # compute the number of free parameters f in the complement
  p <- ncol(Amat)
  # rank q1
  #lq1 <- qr(Amat.ciq)$rank
  lq1 <- qr(Amat)$rank - meq
  
  idx <- length(wt.bar)
  if (attr(wt.bar, "method") == "boot") {
    w_c <- wt.bar[idx-meq]
  } else if (attr(wt.bar, "method") == "pmvnorm") {
    w_c <- wt.bar[idx]
  } else {
    stop("restriktor ERROR: no level probabilities (chi-bar-square weights) found.")
  }
  
  # compute penalty term value PTc (of the GORIC(C))
  if (type %in% c("goric", "gorica")) {
    PTc <- as.numeric(n_nuis + p - w_c * lq1)
  } else if (type %in% c("goricc", "goricac")) {
    if (is.null(sample.nobs)) {
      stop("restriktor ERROR: the argument sample.nobs is not found.")
    }
    N <- sample.nobs
    # small sample correction: mixture of the penalty with p - lq1 and with
    # p free parameters
    PTc <- as.numeric(w_c * small_sample_PT(p - lq1, N, n_nuis, n_resp) +
                        (1 - w_c) * small_sample_PT(p, N, n_nuis, n_resp))
  }
  
  if (debug) {
    cat("penalty term value =", PTc, "\n")
  }
  
  # gorica(c): the nuisance parameters are not counted (asymptotically)
  if (type %in% c("gorica", "goricac")) {
    PTc <- PTc - n_nuis
  }
  
  return(PTc)
}
