# =============================================================================
# Confidence intervals (CIs) for the aggregated (i.e., final) results of a GORICA
# evidence synthesis based on (standardized) estimates and their covariance
# matrices (i.e., input type 'est_vcov' or 'escalc').
# =============================================================================
#
# The studies are independent, so the sampling (co)variance of the aggregated
# log-likelihood values is the sum of the study-specific ones. For each study
# s, the covariance matrix Sigma_s of its log-likelihood values (one per
# hypothesis) and/or the parametric-bootstrap draws of its IC values are
# obtained as in goric_CI() (see ?goric_CI); the penalties and priors are fixed.
#
# Aggregated IC values (as reported by evSyn; without study weights):
#   added:  IC_m = -2 sum_s LL_ms + c sum_s PT_ms
#   equal:  IC_m = -2 sum_s LL_ms + c mean_s PT_ms
# Aggregated IC weights (as reported by evSyn; with study weights a_s, which
# are rescaled to sum to S):
#   eta_m = log(p_m) - IC^w_m / 2, with
#   added:  IC^w_m = -2 sum_s a_s LL_ms + c sum_s a_s PT_ms
#   equal:  IC^w_m = -2 sum_s a_s LL_ms + c sum_s a_s PT_ms / S
# Hence (only the LL values vary):
#   Cov(IC values)   = 4 sum_s Sigma_s
#   Cov(eta)         =   sum_s a_s^2 Sigma_s
# - method = "delta": normal-based CIs for the IC values and their
#   differences, and for the ratios of IC weights (via the CI for log(ratio)),
#   and delta-method CIs (logit scale) for the IC weights.
# - method = "mc": as "delta", but the CIs for the IC weights are obtained by
#   Monte Carlo draws eta* ~ N(eta, Cov(eta)).
# - method = "boot": for each study, B parametric-bootstrap draws
#   theta*_s ~ N(theta^_s, V_s), independently over the studies; the aggregated
#   quantities are computed per draw. This is recommended when hypotheses 
#   overlap (i.e., have the same maximum order-restricted log likelihood, LL in 
#   evSyn()) in (some of) the studies, since the differences in IC values are 
#   then not normally distributed. Two types of bootstrap CIs (boot_type):
#   - "percentile": the percentiles of the aggregated bootstrap draws;
#   - "centered": the bootstrap draws are centered at the observed value 
#     (i.e., their mean shift is removed) and the basic bootstrap CI is used: 
#     [obs - q_hi(Z), obs - q_lo(Z)], with Z = draw - mean(draws), on the 
#     scale of the IC values / differences, log(ratio), and logit(weight).
#   The bootstrap distribution of a difference in IC values is shifted relative 
#   to its sampling distribution (by roughly the number of active constraints 
#   per study). In one study, this shift is small relative to the spread, but 
#   in aggregated results the shifts add up over the studies (growing with S), 
#   whereas the spread grows with sqrt(S). In a simulation study (see 
#   experiments/coverage_CI), the percentile CIs had a coverage of about .95 
#   (or higher) with S = 3 studies, but undercovered with S = 10 studies 
#   (down to about .75), which the centered CIs (largely) corrected. Hence, by 
#   default (boot_type = "auto"), the centered CIs are used from S = 5 studies 
#   onwards and the percentile CIs otherwise.
# The inversion method of goric_CI() (method = "invert") is not available for
# aggregated results, since its one-dimensional reduction is based on one vector
# of estimates; then, the bootstrap is used.
#
# Not available for type_ev = "average": the analyses are then based on the
# same data and are thus dependent, so that their sampling (co)variances
# cannot be added (and they cannot be bootstrapped independently).
# =============================================================================

evSyn_CI <- function(object, level = 0.95, method = c("boot", "delta", "mc"),
                     boot_type = c("auto", "percentile", "centered"),
                     nsim = 10000L, B = 1000L, seed = NULL, tol = 1e-8, rho = 1,
                     MI = isTRUE(object$MI)) {

  if (!inherits(object, "evSyn")) {
    stop("\nrestriktor ERROR: The object should be of class 'evSyn' ",
         "(i.e., output from the evSyn() function).", call. = FALSE)
  }
  goric_list <- object$goric_studies
  # Dependent studies (evSyn(..., cluster = )): two-level evidence synthesis 
  # ("cluster") or robust variance estimation ("rve"); see .evSyn_cluster().
  cl_info   <- object$cluster_info
  rve       <- !is.null(cl_info) && identical(cl_info$dependence, "rve")
  two_level <- !is.null(cl_info) && identical(cl_info$dependence, "cluster")
  if (rve && !is.null(cl_info$study_LL)) {
    # RVE only requires the study-specific log-likelihood and penalty values
  } else if (!inherits(object, c("evSyn_est", "evSyn_gorica")) || is.null(goric_list)) {
    stop("\nrestriktor ERROR: Confidence intervals for the aggregated results are only available when ",
         "the input of evSyn() consists of (standardized) estimates and their covariance matrices ",
         "(input type 'est_vcov' or 'escalc') or of goric objects (input type 'gorica'); for the other ",
         "input types, there is no information about the sampling variability.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("\nrestriktor ERROR: The argument 'level' should be a single number in between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || rho < 0 || rho > 1) {
    stop("\nrestriktor ERROR: The argument 'rho' should be a single number in between 0 and 1.",
         call. = FALSE)
  }
  if (identical(method, "invert")) {
    message("\nrestriktor Message: method = 'invert' is not available for aggregated results; ",
            "the parametric bootstrap (method = 'boot') is used instead.")
    method <- "boot"
  }
  method <- match.arg(method)
  type_ev <- object$type_ev
  # Average approach: the analyses are based on the same data and are thus 
  # dependent. Their covariances are unknown (and the estimates of the analyses 
  # need not be comparable), so the analyses are treated as fully dependent: the 
  # (co)variance of the averaged log-likelihood values is the mean of the 
  # analysis-specific (co)variances, W-bar (an upper bound, by Cauchy-Schwarz, 
  # so the CIs are conservative). They cannot be bootstrapped independently, so 
  # the delta method (with Monte Carlo draws for the weights) is used. The 
  # between-analysis variability (B) is reported separately; it reflects the 
  # analysis choices, not sampling variability, and is not part of the CIs.
  # TO DO: joint case bootstrap (re-running all analyses on the same bootstrap 
  # samples) when the data are available.
  if (!is.null(cl_info) && method == "boot") {
    message("\nrestriktor Message: With dependent studies (argument 'cluster' in evSyn()), the studies ",
            "cannot be bootstrapped independently; method = 'mc' is used",
            if (rve) " (with the cluster-robust covariance matrix)." else
              " (with the studies within a cluster treated as fully dependent).")
    method <- "mc"
  }
  # Multiple imputation (MI = TRUE): the analyses are the analyses of multiple 
  # imputed data sets (NOT different analyses of the same data, as in a 
  # multiverse analysis). The aggregated values and weights are those of the 
  # average approach; their (co)variance is obtained with Rubin's rules: 
  # T = W-bar + (1 + 1/m) B, with W-bar the mean within-imputation covariance 
  # matrix and B the between-imputation covariance matrix (of the log-likelihood 
  # values minus the penalty part), and t quantiles with Rubin's degrees of freedom.
  MI <- isTRUE(MI)
  if (MI && type_ev != "average") {
    stop("\nrestriktor ERROR: MI = TRUE (multiple imputation) requires type_ev = 'average', with one ",
         "analysis per imputed data set.", call. = FALSE)
  }
  if (MI && method == "boot") {
    message("\nrestriktor Message: With multiple imputation (MI = TRUE), Rubin's rules are used ",
            "(method = 'mc').")
    method <- "mc"
  }
  if (type_ev == "average" && method == "boot") {
    message("\nrestriktor Message: For type_ev = 'average', the analyses are dependent and cannot be ",
            "bootstrapped independently; method = 'mc' is used (with a correlation of rho = ", format(rho),
            " between the sampling errors of the analyses).")
    method <- "mc"
  }
  if (!is.null(seed)) set.seed(seed)

  S      <- if (rve) cl_info$n_studies else if (two_level) cl_info$n_clusters else length(goric_list)
  type   <- object$type
  pf     <- if (is.null(object$penalty_factor)) 2 else object$penalty_factor
  prior  <- object$priorICweights
  models <- colnames(object$LL_m)
  M      <- length(models)
  if (is.null(prior) || length(prior) != M) prior <- rep(1 / M, M)
  prior  <- prior / sum(prior)

  # Study weights (input order; as in evSyn: rescaled to sum to S)
  sw <- object$study_weights
  if (is.null(sw)) sw <- rep(1 / S, S)
  a  <- S * sw / sum(sw)

  # Study-specific log-likelihood and penalty values (input order of the studies)
  if (rve) {
    LL_s <- unname(as.matrix(cl_info$study_LL)); PT_s <- unname(as.matrix(cl_info$study_PT))
  } else if (two_level) {
    # cluster level (input order of the clusters)
    pos_c <- match(seq_len(S), object$order_studies)
    if (anyNA(pos_c)) pos_c <- seq_len(S)
    LL_s <- unname(as.matrix(object$LL_m))[pos_c, , drop = FALSE]
    PT_s <- unname(as.matrix(object$PT_m))[pos_c, , drop = FALSE]
  } else {
    LL_s <- t(vapply(goric_list, function(g) g$result$loglik, numeric(M)))
    PT_s <- t(vapply(goric_list, function(g) g$result$penalty, numeric(M)))
    if (M == 1L) { LL_s <- t(LL_s); PT_s <- t(PT_s) }
  }
  colnames(LL_s) <- colnames(PT_s) <- models

  # Aggregated IC values (without study weights) and the IC values used for the
  # weights (with study weights), as in evSyn
  pool_values <- function(LL) {   # LL: S x M (or a 3-dim array: draws x S x M)
    if (type_ev == "added") {
      -2 * colSums(LL) + pf * colSums(PT_s)
    } else if (type_ev == "equal") {
      -2 * colSums(LL) + pf * colMeans(PT_s)
    } else {   # average
      (-2 * colSums(LL) + pf * colSums(PT_s)) / S
    }
  }
  pool_weighted <- function(LL) {
    if (type_ev == "added") {
      -2 * colSums(LL * a) + pf * colSums(PT_s * a)
    } else if (type_ev == "equal") {
      -2 * colSums(LL * a) + pf * colSums(PT_s * a) / S
    } else {   # average
      (-2 * colSums(LL * a) + pf * colSums(PT_s * a)) / S
    }
  }
  IC   <- pool_values(LL_s);   names(IC) <- models
  ICw  <- pool_weighted(LL_s); names(ICw) <- models
  eta  <- log(prior) - ICw / 2
  w    <- exp(eta - max(eta)); w <- w / sum(w)

  # check: same as the final results of the evSyn object
  chk_IC <- object$Cumulative_GORICA["Final", models]
  chk_w  <- object$Cumulative_GORICA_weights["Final", models]
  if (max(abs(chk_IC - IC)) > 1e-6 * max(1, abs(IC)) || max(abs(chk_w - w)) > 1e-6) {
    warning("\nrestriktor WARNING: The aggregated ", toupper(type), " values and/or weights could not be ",
            "reproduced from the study-specific goric objects; please check the results.", call. = FALSE)
  }

  # Study-specific covariance matrices of the log-likelihood values (delta
  # method based on the estimates, as in goric_CI()), and the number of studies
  # in which a pair of hypotheses has the same maximum order-restricted log likelihood (LL)
  # Note: the hypotheses are matched by position, since the labels in the goric 
  # objects can differ from those in the evSyn object (e.g., 'unconstrained' 
  # versus 'Unconstrained', or when 'hypo_names' is used).
  pairs <- expand.grid(m = seq_len(M), k = seq_len(M))
  pairs <- pairs[pairs$m != pairs$k, , drop = FALSE]
  pairs <- pairs[order(pairs$m, pairs$k), , drop = FALSE]
  same_LL <- function(LL) {   # per pair: number of rows (studies/clusters) with the same LL
    vapply(seq_len(nrow(pairs)), function(j) {
      sum(abs(LL[, pairs$m[j]] - LL[, pairs$k[j]]) < tol * pmax(1, abs(LL[, pairs$m[j]])))
    }, numeric(1))
  }
  df_t <- Inf
  if (rve) {
    # Robust variance estimation: cluster-robust covariance matrix of the sum 
    # of the log-likelihood values (residuals around their (weighted) mean), 
    # with the small-sample factor C / (C - 1) and t(C - 1) quantiles.
    cl <- cl_info$cluster; levs <- cl_info$levels; C <- length(levs)
    rob <- function(wts) {
      mu <- colSums(LL_s * wts) / sum(wts)
      E  <- sweep(LL_s, 2, mu, "-") * wts
      G  <- rowsum(E, cl)                      # C x M
      C / (C - 1) * crossprod(G)
    }
    Sigma_v <- rob(rep(1, S))
    Sigma_w <- rob(a)
    df_t <- C - 1
    n_same <- same_LL(LL_s)
  } else {
  models_g <- lapply(goric_list, function(g) as.character(g$result$model))
  if (any(lengths(models_g) != M)) {
    stop("\nrestriktor ERROR: The number of hypotheses in the study-specific goric objects ",
         "does not match the evSyn object.", call. = FALSE)
  }
  # (approach "auto": for GORIC(C) objects, the casewise log-likelihood 
  # contributions; for GORICA(C) objects, the estimates; as in goric_CI())
  ci_s <- lapply(goric_list, function(g) suppressMessages(suppressWarnings(
    goric_CI(g, level = level, method = "delta", approach = "auto"))))
  if (two_level) {
    # studies within a cluster: dependent, with correlation rho between their 
    # (log-likelihood) sampling errors (rho = 1: fully dependent); the clusters 
    # are independent.
    cl <- cl_info$cluster; levs <- cl_info$levels; sw_st <- cl_info$study_weights
    if (is.null(sw_st)) sw_st <- rep(1, length(cl))
    Sig_c <- lapply(levs, function(lv) {
      idx <- which(cl == lv); ww <- sw_st[idx] / sum(sw_st[idx])
      .evSyn_dep_cov(lapply(ci_s[idx], function(ci) unname(ci$Sigma)), ww, rho)
    })
    Sigma_v <- Reduce(`+`, Sig_c)
    Sigma_w <- Reduce(`+`, Map(function(Sg, as) as^2 * Sg, Sig_c, a))
  } else if (MI) {
    # Rubin's rules (multiple imputation): T = W-bar + (1 + 1/m) B
    W_mi <- Reduce(`+`, lapply(ci_s, function(ci) unname(ci$Sigma))) / S
    Lt <- LL_s - (pf / 2) * PT_s                      # -IC / 2 per imputation
    B_mi <- if (S > 1) stats::cov(Lt) else matrix(0, M, M)
    Sigma_v <- Sigma_w <- W_mi + (1 + 1 / S) * B_mi
    if (length(unique(round(a, 12))) > 1) {
      message("\nrestriktor Message: With multiple imputation, the imputed data sets should have equal ",
              "(study) weights; Rubin's rules are applied without weights.")
    }
  } else if (type_ev == "average") {
    # dependent analyses, with correlation rho between their (log-likelihood) 
    # sampling errors (rho = 1: fully dependent): covariance matrix of the 
    # (weighted) mean of the log-likelihood values
    Sig_list <- lapply(ci_s, function(ci) unname(ci$Sigma))
    Sigma_v <- .evSyn_dep_cov(Sig_list, rep(1 / S, S), rho)
    Sigma_w <- .evSyn_dep_cov(Sig_list, a / S, rho)
  } else {
    Sigma_v <- Reduce(`+`, lapply(ci_s, function(ci) unname(ci$Sigma)))
    Sigma_w <- Reduce(`+`, Map(function(ci, as) as^2 * unname(ci$Sigma), ci_s, a))
  }
  n_same <- if (two_level) same_LL(LL_s) else Reduce(`+`, Map(function(ci, mg) {
    d <- ci$difference
    key <- paste(mg[pairs$m], mg[pairs$k], sep = "\r")
    as.integer(d$same_loglik[match(key, paste(d$model, d$vs, sep = "\r"))])
  }, ci_s, models_g))
  }
  dimnames(Sigma_v) <- dimnames(Sigma_w) <- list(models, models)

  z <- if (is.finite(df_t)) stats::qt(1 - (1 - level) / 2, df = df_t) else qnorm(1 - (1 - level) / 2)
  # Multiple imputation: Rubin's degrees of freedom per contrast c (of the 
  # log-likelihood values): nu = (m - 1) (1 + 1 / r)^2, r = (1 + 1/m) c'Bc / c'Wc
  rubin_df <- function(cvec) {
    if (!MI) return(Inf)
    wv <- drop(t(cvec) %*% W_mi %*% cvec); bv <- drop(t(cvec) %*% B_mi %*% cvec)
    if (!is.finite(bv) || bv <= 0) return(Inf)
    if (!is.finite(wv) || wv <= 0) return(S - 1)
    r <- (1 + 1 / S) * bv / wv
    (S - 1) * (1 + 1 / r)^2
  }
  qz <- function(nu) if (is.finite(nu)) stats::qt(1 - (1 - level) / 2, df = nu) else z
  df_IC <- vapply(seq_len(M), function(m) { cv <- numeric(M); cv[m] <- 1; rubin_df(cv) }, numeric(1))
  df_d  <- vapply(seq_len(nrow(pairs)), function(j) {
    cv <- numeric(M); cv[pairs$m[j]] <- 1; cv[pairs$k[j]] <- -1; rubin_df(cv) }, numeric(1))
  z_IC <- vapply(df_IC, qz, numeric(1)); z_d <- vapply(df_d, qz, numeric(1))
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)

  var_pair <- function(Sig) Sig[cbind(pairs$m, pairs$m)] + Sig[cbind(pairs$k, pairs$k)] -
    2 * Sig[cbind(pairs$m, pairs$k)]
  se_d   <- 2 * sqrt(pmax(var_pair(Sigma_v), 0))   # SE of differences in aggregated IC values
  se_dw  <- 2 * sqrt(pmax(var_pair(Sigma_w), 0))   # idem, for the (weighted) log ratios
  dIC    <- IC[pairs$m] - IC[pairs$k]
  log_ratio <- eta[pairs$m] - eta[pairs$k]
  scale_ref <- max(1, sqrt(max(abs(diag(Sigma_v)))))
  same_all  <- !(se_d > tol * scale_ref)   # same maximum order-restricted LL in all studies

  use_boot <- method == "boot"
  boot_type <- match.arg(boot_type)
  if (boot_type == "auto") boot_type <- if (S >= 5) "centered" else "percentile"
  if (!use_boot) boot_type <- NA_character_
  # Bootstrap CI from draws 'x' and observed value 'obs' (on a scale where the 
  # centering is applied: IC values / differences, log(ratio), logit(weight))
  boot_ci <- function(x, obs) {
    if (identical(boot_type, "centered")) {
      z <- x - mean(x)
      obs - rev(quantile(z, probs = probs, names = FALSE))
    } else {
      quantile(x, probs = probs, names = FALSE)
    }
  }
  if (use_boot) {
    # Bootstrap draws per study (independent over the studies):
    # IC*_s = const_s + D*_s + c * PT_s, so LL*_s = -(IC*_s - c * PT_s) / 2
    IC_vals_b <- matrix(0, B, M); IC_w_b <- matrix(0, B, M)
    # For GORIC(C) objects, the bootstrap uses the GORICA-type log-likelihood 
    # (asymptotically equivalent for differences; as in goric_CI()); its draws 
    # are shifted such that they are centered at the observed GORIC(C) 
    # log-likelihood values (for GORICA(C) objects, the shift is zero).
    for (s in seq_len(S)) {
      ICb <- unname(goric_CI_boot(goric_list[[s]], models_g[[s]], B = B, penalty_factor = pf))  # B x M
      LLb <- -(sweep(ICb, 2, pf * PT_s[s, ], "-")) / 2
      if (type %in% c("goric", "goricc")) {
        IC0 <- unname(goric_CI_boot(goric_list[[s]], models_g[[s]], B = 0L, penalty_factor = pf))
        LL0 <- -(IC0 - pf * PT_s[s, ]) / 2
        LLb <- sweep(LLb, 2, LL_s[s, ] - LL0, "+")
      }
      IC_vals_b <- IC_vals_b - 2 * LLb
      IC_w_b    <- IC_w_b - 2 * a[s] * LLb
    }
    if (type_ev == "added") {
      IC_vals_b <- sweep(IC_vals_b, 2, pf * colSums(PT_s), "+")
      IC_w_b    <- sweep(IC_w_b, 2, pf * colSums(PT_s * a), "+")
    } else {
      IC_vals_b <- sweep(IC_vals_b, 2, pf * colMeans(PT_s), "+")
      IC_w_b    <- sweep(IC_w_b, 2, pf * colSums(PT_s * a) / S, "+")
    }
    eta_b <- sweep(-IC_w_b / 2, 2, log(prior), "+")
  }

  # 1. Aggregated IC values
  se_IC <- 2 * sqrt(pmax(diag(Sigma_v), 0))
  q_IC <- if (use_boot && type %in% c("gorica", "goricac")) {
    t(vapply(seq_len(M), function(m) boot_ci(IC_vals_b[, m], IC[m]), numeric(2)))
  } else {
    cbind(IC - z_IC * se_IC, IC + z_IC * se_IC)
  }
  IC_tab <- data.frame(model = models, IC = IC, SE = se_IC,
                       lower = q_IC[, 1], upper = q_IC[, 2],
                       row.names = NULL, stringsAsFactors = FALSE)
  names(IC_tab)[2] <- type

  # 2. Differences in aggregated IC values and 3. ratios of aggregated IC weights
  if (use_boot) {
    q_d <- t(vapply(seq_len(nrow(pairs)), function(j) {
      boot_ci(IC_vals_b[, pairs$m[j]] - IC_vals_b[, pairs$k[j]], dIC[j])
    }, numeric(2)))
    q_r <- t(vapply(seq_len(nrow(pairs)), function(j) {
      exp(boot_ci(eta_b[, pairs$m[j]] - eta_b[, pairs$k[j]], log_ratio[j]))
    }, numeric(2)))
  } else {
    q_d <- cbind(ifelse(!same_all, dIC - z_d * se_d, NA_real_),
                 ifelse(!same_all, dIC + z_d * se_d, NA_real_))
    q_r <- cbind(ifelse(!same_all, exp(log_ratio - z_d * se_dw / 2), NA_real_),
                 ifelse(!same_all, exp(log_ratio + z_d * se_dw / 2), NA_real_))
  }
  ci_method <- rep(if (use_boot) { if (boot_type == "centered") "boot (centered)" else "boot" } else "analytic", 
                   nrow(pairs))
  diff_tab <- data.frame(model = models[pairs$m], vs = models[pairs$k],
                         difference = dIC, SE = se_d,
                         lower = q_d[, 1], upper = q_d[, 2],
                         n_same_loglik = n_same,
                         CI_method = ci_method,
                         row.names = NULL, stringsAsFactors = FALSE)
  ratio_tab <- data.frame(model = models[pairs$m], vs = models[pairs$k],
                          ratio = exp(log_ratio),
                          lower = q_r[, 1], upper = q_r[, 2],
                          n_same_loglik = n_same,
                          CI_method = ci_method,
                          row.names = NULL, stringsAsFactors = FALSE)

  # 4. Aggregated IC weights
  if (method == "delta") {
    ww <- matrix(NA_real_, M, 2)
    for (m in seq_len(M)) {
      if (M == 1L) break
      g <- -w / (1 - w[m]); g[m] <- 1
      if (!is.finite(sum(g))) { g <- rep(0, M); g[m] <- 1 }
      se_logit <- sqrt(max(drop(t(g) %*% Sigma_w %*% g), 0))
      wm <- min(max(w[m], 1e-300), 1 - 1e-16)
      ww[m, ] <- plogis(qlogis(wm) + c(-1, 1) * qz(rubin_df(g)) * se_logit)
    }
  } else {
    if (method == "mc") {
      ev <- eigen((Sigma_w + t(Sigma_w)) / 2, symmetric = TRUE)
      A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), M)
      Z <- matrix(rnorm(nsim * M), nsim, M)
      if (is.finite(df_t)) Z <- Z / sqrt(stats::rchisq(nsim, df = df_t) / df_t)   # multivariate t (RVE)
      if (MI && is.finite(min(df_d))) {                                             # multivariate t (MI)
        df_mi <- max(min(df_d), 1); Z <- Z / sqrt(stats::rchisq(nsim, df = df_mi) / df_mi)
      }
      eta_sim <- sweep(Z %*% t(A), 2, eta, "+")
    } else {
      eta_sim <- eta_b
    }
    eta_sim <- eta_sim - apply(eta_sim, 1, max)
    w_sim <- exp(eta_sim); w_sim <- w_sim / rowSums(w_sim)
    if (identical(boot_type, "centered")) {
      # centered bootstrap on the logit scale of the weights
      clip <- function(p) pmin(pmax(p, 1e-300), 1 - 1e-16)
      ww <- t(vapply(seq_len(M), function(m) {
        plogis(boot_ci(qlogis(clip(w_sim[, m])), qlogis(clip(w[m]))))
      }, numeric(2)))
    } else {
      ww <- t(apply(w_sim, 2, quantile, probs = probs, names = FALSE))
    }
  }
  weight_tab <- data.frame(model = models, weight = w,
                           lower = ww[, 1], upper = ww[, 2],
                           CI_method = rep(if (identical(boot_type, "centered")) "boot (centered)" else method, M),
                           row.names = NULL, stringsAsFactors = FALSE)
  names(weight_tab)[2] <- paste0(type, ".weights")

  if (MI) {
    IC_tab$df <- df_IC
    diff_tab$df <- df_d
    ratio_tab$df <- df_d
  }
  # Average approach: between-analysis variability (SD over the analyses of the 
  # IC values, differences in IC values, and log ratios of IC weights); not part 
  # of the CIs (see above).
  if (type_ev == "average") {
    IC_s  <- -2 * LL_s + pf * PT_s                               # S x M
    eta_s <- sweep(-IC_s / 2, 2, log(prior), "+")
    IC_tab$SD_between <- apply(IC_s, 2, stats::sd)
    diff_tab$SD_between <- vapply(seq_len(nrow(pairs)), function(j)
      stats::sd(IC_s[, pairs$m[j]] - IC_s[, pairs$k[j]]), numeric(1))
    ratio_tab$SD_between_log <- vapply(seq_len(nrow(pairs)), function(j)
      stats::sd(eta_s[, pairs$m[j]] - eta_s[, pairs$k[j]]), numeric(1))
  }

  out <- list(type = type, type_ev = type_ev, level = level, method = method, boot_type = boot_type,
              n_studies = S, B = if (use_boot) B else NULL, nsim = if (method == "mc") nsim else NULL,
              penalty_factor = pf, priorICweights = prior, study_weights = sw,
              Sigma_values = Sigma_v, Sigma_eta = Sigma_w,
              dependence = if (is.null(cl_info)) NULL else cl_info$dependence,
              rho = if ((type_ev == "average" && !MI) || two_level) rho else NULL,
              MI = MI, df_MI = if (MI) list(IC = df_IC, difference = df_d) else NULL,
              n_clusters = if (is.null(cl_info)) NULL else cl_info$n_clusters, df = df_t,
              IC = IC_tab, difference = diff_tab,
              weights = weight_tab, ratio_weights = ratio_tab,
              best = models[which.max(w)],
              any_same_loglik = any(n_same > 0),
              weighted = length(unique(round(a, 12))) > 1)
  class(out) <- "evSyn_CI"
  out
}


print.evSyn_CI <- function(x, digits = max(3, getOption("digits") - 4),
                           all_pairs = FALSE, ...) {
  type_label <- toupper(x$type)
  pct <- paste0(format(100 * x$level), "%")
  cat(sprintf("restriktor: %s confidence intervals for the aggregated (final) %s values, differences, weights, and ratios of weights\n",
              pct, type_label))
  units <- .evSyn_unit_label(x$type_ev, plural = TRUE, capital = FALSE)
  cat(sprintf("(%s evidence synthesis of %s %s; %s)\n", x$type_ev, x$n_studies, units,
              switch(x$method,
                     boot  = paste0("parametric-bootstrap ", 
                                    if (identical(x$boot_type, "centered")) "CIs centered at the observed values (basic bootstrap)" else "percentile CIs",
                                    ", B = ", x$B,
                                    " draws from N(estimates, VCOV) per study, independently over the studies",
                                    if (x$type %in% c("goric", "goricc")) 
                                      paste0("; CIs of the aggregated ", type_label, " values: casewise approach") else ""),
                     delta = "delta method based on the estimates and their covariance matrices; weights: delta method on the logit scale",
                     mc    = paste0("delta method based on the estimates and their covariance matrices; weights: Monte Carlo (", x$nsim, " draws)"))))
  if (isTRUE(x$MI)) {
    cat(paste0("(Multiple imputation: the ", x$n_studies, " analyses are analyses of multiple imputed data sets; ",
               "the aggregated values and weights are those of the average approach, and their CIs are based on ",
               "Rubin's rules: total variance = within-imputation variance + (1 + 1/m) x between-imputation ",
               "variance, with t quantiles with Rubin's degrees of freedom (df).)\n"))
  } else if (x$type_ev == "average") {
    cat(paste0("(The analyses are based on the same data and are dependent: the CIs assume a correlation of rho = ",
               format(x$rho), " between the sampling errors of the analyses",
               if (isTRUE(x$rho == 1)) " (fully dependent; conservative)" else "", ".)\n"))
  }
  if (identical(x$dependence, "cluster")) {
    cat(paste0("(Two-level evidence synthesis of ", x$n_clusters, " clusters of dependent studies: the CIs assume ",
               "a correlation of rho = ", format(x$rho), " between the sampling errors of the studies within a cluster",
               if (isTRUE(x$rho == 1)) " (fully dependent; conservative)" else "",
               ", and independent clusters.)\n"))
  } else if (identical(x$dependence, "rve")) {
    cat(paste0("(Robust variance estimation (RVE) with ", x$n_clusters, " clusters of dependent studies: ",
               "cluster-robust covariance matrix and t(", x$n_clusters - 1, ") quantiles; these CIs also ",
               "include the heterogeneity of the studies.)\n"))
  }

  fmt <- function(df) .format_CI_table(df, digits = digits)

  cat(sprintf("\nAggregated %s values:\n", type_label))
  print(fmt(x$IC), row.names = FALSE, right = TRUE)

  sel <- if (all_pairs) rep(TRUE, nrow(x$difference)) else x$difference$model == x$best
  cat(sprintf("\nDifferences in aggregated %s values%s:\n", type_label,
              ifelse(all_pairs, "", paste0(" (best hypothesis ", x$best, " versus the others)"))))
  print(fmt(x$difference[sel, , drop = FALSE]), row.names = FALSE, right = TRUE)

  cat(sprintf("\nAggregated %s weights:\n", type_label))
  print(fmt(x$weights), row.names = FALSE, right = TRUE)

  cat(sprintf("\nRatios of aggregated %s weights%s:\n", type_label,
              ifelse(all_pairs, "", paste0(" (best hypothesis ", x$best, " versus the others)"))))
  print(fmt(x$ratio_weights[sel, , drop = FALSE]), row.names = FALSE, right = TRUE)

  cat(paste0("\n(n_same_loglik: number of ", units, " in which the two hypotheses have the same maximum ",
             "order-restricted log likelihood (i.e., LL in evSyn()), and thus the same fit.)\n"))
  if (isTRUE(x$MI)) {
    cat(paste0("(SD_between, SD_between_log: standard deviation over the imputed data sets of the ", type_label,
               " values, of their differences, and of the log ratios of ", type_label, " weights, respectively; ",
               "this between-imputation variability is part of the CIs (Rubin's rules).)\n"))
    message("\nrestriktor Note: MI = TRUE is only appropriate for analyses of multiple imputed data sets, ",
            "not for different analyses of the same data (e.g., a multiverse analysis); for the latter, use ",
            "MI = FALSE (with the argument rho).")
  } else if (x$type_ev == "average") {
    cat(paste0("(SD_between, SD_between_log: standard deviation over the analyses of the ", type_label,
               " values, of their differences, and of the log ratios of ", type_label, " weights, respectively; ",
               "this between-analysis variability reflects the analysis choices, not sampling variability, ",
               "and is not part of the CIs.)\n"))
  }
  if (isTRUE(x$weighted)) {
    message("\nrestriktor Note: The study weights are used in the aggregated ", type_label,
            " weights and ratios of weights (and their CIs), not in the aggregated ", type_label,
            " values and their differences (as in evSyn).")
  }
  if (identical(x$boot_type, "percentile") && x$n_studies >= 5) {
    message("\nrestriktor Note: With ", x$n_studies, " studies, the percentile bootstrap CIs can be too ",
            "narrow (the shift of the bootstrap distribution adds up over the studies); ",
            "boot_type = \"centered\" (the default from 5 studies onwards) corrects this largely.")
  }
  if (!is.null(x$priorICweights) && length(unique(round(x$priorICweights, 12))) > 1) {
    message("\nrestriktor Note: Unequal prior IC weights (priorICweights) are used. These are part of the ",
            "aggregated ", type_label, " weights and the ratios of ", type_label, " weights (and their CIs), ",
            "but not of the aggregated ", type_label, " values and their differences (and their CIs). ",
            "Hence, conclusions based on the (CIs of the) ", type_label, " values may differ from ",
            "those based on the (CIs of the) ", type_label, " weights.")
  }
  if (isTRUE(x$any_same_loglik) && x$method != "boot" && x$type_ev != "average") {
    message("\nrestriktor Note: In some studies, hypotheses have the same maximum order-restricted ",
            "log likelihood (i.e., LL in evSyn()) (n_same_loglik > 0), for example, a hypothesis and the unconstrained hypothesis ",
            "when the estimates are in agreement with that hypothesis. The differences in ", type_label,
            " values are then not normally distributed; the normal-based CIs ignore this (order 1) ",
            "variability for these studies and are too narrow (no CI is given if this holds in all studies). ",
            "Use method = \"boot\" instead.")
  }
  if (x$type %in% c("gorica", "goricac") && x$method != "boot") {
    message("\nrestriktor Note: The CI of a single aggregated ", type_label, " value is based on ",
            "the squared Mahalanobis distances to the hypothesis; its SE is 0 for hypotheses ",
            "in agreement with the estimates in all studies. Only the relative values ",
            "(differences, weights, ratios) are of interest anyway.")
  }
  invisible(x)
}


# Covariance matrix of a weighted mean of dependent (log-likelihood) vectors, 
# with covariance matrices Sig_list[[s]] and weights w (summing to 1), when the 
# sampling errors of two of them have correlation rho (cross-covariance 
# rho * A_s A_t, with A_s the symmetric square root of Sig_list[[s]]):
#   (1 - rho) * sum_s w_s^2 Sigma_s + rho * (sum_s w_s A_s)^2.
# For rho = 1 (fully dependent), this is the maximal variance (for one 
# dimension: (sum_s w_s SD_s)^2); for rho = 0, the analyses are independent.
.evSyn_dep_cov <- function(Sig_list, w, rho = 1) {
  sqrtm <- function(S) {
    e <- eigen((S + t(S)) / 2, symmetric = TRUE)
    e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(S)) %*% t(e$vectors)
  }
  indep <- Reduce(`+`, Map(function(S, w1) w1^2 * S, Sig_list, w))
  if (rho == 0) return(indep)
  Abar <- Reduce(`+`, Map(function(S, w1) w1 * sqrtm(S), Sig_list, w))
  full <- Abar %*% Abar
  (1 - rho) * indep + rho * (full + t(full)) / 2
}
