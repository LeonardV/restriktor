# =============================================================================
# Confidence intervals (CIs) for the pooled (i.e., final) results of a GORICA
# evidence synthesis based on (standardized) estimates and their covariance
# matrices (i.e., input type 'est_vcov' or 'escalc').
# =============================================================================
#
# The studies are independent, so the sampling (co)variance of the pooled
# log-likelihood values is the sum of the study-specific ones. For each study
# s, the covariance matrix Sigma_s of its log-likelihood values (one per
# hypothesis) and/or the parametric-bootstrap draws of its IC values are
# obtained as in goric_CI() (see ?goric_CI); the penalties and priors are fixed.
#
# Pooled IC values (as reported by evSyn; without study weights):
#   added:  IC_m = -2 sum_s LL_ms + c sum_s PT_ms
#   equal:  IC_m = -2 sum_s LL_ms + c mean_s PT_ms
# Pooled IC weights (as reported by evSyn; with study weights a_s, which
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
#   theta*_s ~ N(theta^_s, V_s), independently over the studies; the pooled
#   quantities are computed per draw, and the percentile CIs are reported.
#   This is recommended when hypotheses overlap (i.e., have the same
#   order-restricted MLEs) in (some of) the studies, since the differences in
#   IC values are then not normally distributed.
# The inversion method of goric_CI() (method = "invert") is not available for
# pooled results, since its one-dimensional reduction is based on one vector
# of estimates; then, the bootstrap is used.
#
# Not available for type_ev = "average": the analyses are then based on the
# same data and are thus dependent, so that their sampling (co)variances
# cannot be added (and they cannot be bootstrapped independently).
# =============================================================================

evSyn_CI <- function(object, level = 0.95, method = c("boot", "delta", "mc"),
                     nsim = 10000L, B = 1000L, seed = NULL, tol = 1e-8) {

  if (!inherits(object, "evSyn")) {
    stop("\nrestriktor ERROR: The object should be of class 'evSyn' ",
         "(i.e., output from the evSyn() function).", call. = FALSE)
  }
  goric_list <- object$goric_studies
  if (!inherits(object, c("evSyn_est", "evSyn_gorica")) || is.null(goric_list)) {
    stop("\nrestriktor ERROR: Confidence intervals for the pooled results are only available when ",
         "the input of evSyn() consists of (standardized) estimates and their covariance matrices ",
         "(input type 'est_vcov' or 'escalc') or of goric objects (input type 'gorica'); for the other ",
         "input types, there is no information about the sampling variability.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("\nrestriktor ERROR: The argument 'level' should be a single number in between 0 and 1.",
         call. = FALSE)
  }
  if (identical(method, "invert")) {
    message("\nrestriktor Message: method = 'invert' is not available for pooled results; ",
            "the parametric bootstrap (method = 'boot') is used instead.")
    method <- "boot"
  }
  method <- match.arg(method)
  type_ev <- object$type_ev
  if (type_ev == "average") {
    stop("\nrestriktor ERROR: Confidence intervals for the pooled results are not available for ",
         "type_ev = 'average'. The analyses are then based on the same data and are thus dependent, ",
         "so that their sampling variability cannot be combined from the analysis-specific results. ",
         "(This would require re-sampling the data and re-running all analyses.)", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  S      <- length(goric_list)
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
  LL_s <- t(vapply(goric_list, function(g) g$result$loglik, numeric(M)))
  PT_s <- t(vapply(goric_list, function(g) g$result$penalty, numeric(M)))
  if (M == 1L) { LL_s <- t(LL_s); PT_s <- t(PT_s) }
  colnames(LL_s) <- colnames(PT_s) <- models

  # Pooled IC values (without study weights) and the IC values used for the
  # weights (with study weights), as in evSyn
  pool_values <- function(LL) {   # LL: S x M (or a 3-dim array: draws x S x M)
    if (type_ev == "added") {
      -2 * colSums(LL) + pf * colSums(PT_s)
    } else {
      -2 * colSums(LL) + pf * colMeans(PT_s)
    }
  }
  pool_weighted <- function(LL) {
    if (type_ev == "added") {
      -2 * colSums(LL * a) + pf * colSums(PT_s * a)
    } else {
      -2 * colSums(LL * a) + pf * colSums(PT_s * a) / S
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
    warning("\nrestriktor WARNING: The pooled ", toupper(type), " values and/or weights could not be ",
            "reproduced from the study-specific goric objects; please check the results.", call. = FALSE)
  }

  # Study-specific covariance matrices of the log-likelihood values (delta
  # method based on the estimates, as in goric_CI()), and the number of studies
  # in which a pair of hypotheses has the same or-mle's
  # Note: the hypotheses are matched by position, since the labels in the goric 
  # objects can differ from those in the evSyn object (e.g., 'unconstrained' 
  # versus 'Unconstrained', or when 'hypo_names' is used).
  models_g <- lapply(goric_list, function(g) as.character(g$result$model))
  if (any(lengths(models_g) != M)) {
    stop("\nrestriktor ERROR: The number of hypotheses in the study-specific goric objects ",
         "does not match the evSyn object.", call. = FALSE)
  }
  # (approach "auto": for GORIC(C) objects, the casewise log-likelihood 
  # contributions; for GORICA(C) objects, the estimates; as in goric_CI())
  ci_s <- lapply(goric_list, function(g) suppressMessages(suppressWarnings(
    goric_CI(g, level = level, method = "delta", approach = "auto"))))
  Sigma_v <- Reduce(`+`, lapply(ci_s, function(ci) unname(ci$Sigma)))
  Sigma_w <- Reduce(`+`, Map(function(ci, as) as^2 * unname(ci$Sigma), ci_s, a))
  dimnames(Sigma_v) <- dimnames(Sigma_w) <- list(models, models)

  pairs <- expand.grid(m = seq_len(M), k = seq_len(M))
  pairs <- pairs[pairs$m != pairs$k, , drop = FALSE]
  pairs <- pairs[order(pairs$m, pairs$k), , drop = FALSE]
  n_same <- Reduce(`+`, Map(function(ci, mg) {
    d <- ci$difference
    key <- paste(mg[pairs$m], mg[pairs$k], sep = "\r")
    as.integer(d$same_ormle[match(key, paste(d$model, d$vs, sep = "\r"))])
  }, ci_s, models_g))

  z <- qnorm(1 - (1 - level) / 2)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)

  var_pair <- function(Sig) Sig[cbind(pairs$m, pairs$m)] + Sig[cbind(pairs$k, pairs$k)] -
    2 * Sig[cbind(pairs$m, pairs$k)]
  se_d   <- 2 * sqrt(pmax(var_pair(Sigma_v), 0))   # SE of differences in pooled IC values
  se_dw  <- 2 * sqrt(pmax(var_pair(Sigma_w), 0))   # idem, for the (weighted) log ratios
  dIC    <- IC[pairs$m] - IC[pairs$k]
  log_ratio <- eta[pairs$m] - eta[pairs$k]
  scale_ref <- max(1, sqrt(max(abs(diag(Sigma_v)))))
  same_all  <- !(se_d > tol * scale_ref)   # same or-mle's in all studies

  use_boot <- method == "boot"
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

  # 1. Pooled IC values
  se_IC <- 2 * sqrt(pmax(diag(Sigma_v), 0))
  q_IC <- if (use_boot && type %in% c("gorica", "goricac")) {
    t(apply(IC_vals_b, 2, quantile, probs = probs, names = FALSE))
  } else {
    cbind(IC - z * se_IC, IC + z * se_IC)
  }
  IC_tab <- data.frame(model = models, IC = IC, SE = se_IC,
                       lower = q_IC[, 1], upper = q_IC[, 2],
                       row.names = NULL, stringsAsFactors = FALSE)
  names(IC_tab)[2] <- type

  # 2. Differences in pooled IC values and 3. ratios of pooled IC weights
  if (use_boot) {
    q_d <- t(vapply(seq_len(nrow(pairs)), function(j) {
      quantile(IC_vals_b[, pairs$m[j]] - IC_vals_b[, pairs$k[j]], probs = probs, names = FALSE)
    }, numeric(2)))
    q_r <- t(vapply(seq_len(nrow(pairs)), function(j) {
      exp(quantile(eta_b[, pairs$m[j]] - eta_b[, pairs$k[j]], probs = probs, names = FALSE))
    }, numeric(2)))
  } else {
    q_d <- cbind(ifelse(!same_all, dIC - z * se_d, NA_real_),
                 ifelse(!same_all, dIC + z * se_d, NA_real_))
    q_r <- cbind(ifelse(!same_all, exp(log_ratio - z * se_dw / 2), NA_real_),
                 ifelse(!same_all, exp(log_ratio + z * se_dw / 2), NA_real_))
  }
  ci_method <- rep(if (use_boot) "boot" else "analytic", nrow(pairs))
  diff_tab <- data.frame(model = models[pairs$m], vs = models[pairs$k],
                         difference = dIC, SE = se_d,
                         lower = q_d[, 1], upper = q_d[, 2],
                         n_same_ormle = n_same,
                         CI_method = ci_method,
                         row.names = NULL, stringsAsFactors = FALSE)
  ratio_tab <- data.frame(model = models[pairs$m], vs = models[pairs$k],
                          ratio = exp(log_ratio),
                          lower = q_r[, 1], upper = q_r[, 2],
                          n_same_ormle = n_same,
                          CI_method = ci_method,
                          row.names = NULL, stringsAsFactors = FALSE)

  # 4. Pooled IC weights
  if (method == "delta") {
    ww <- matrix(NA_real_, M, 2)
    for (m in seq_len(M)) {
      if (M == 1L) break
      g <- -w / (1 - w[m]); g[m] <- 1
      if (!is.finite(sum(g))) { g <- rep(0, M); g[m] <- 1 }
      se_logit <- sqrt(max(drop(t(g) %*% Sigma_w %*% g), 0))
      wm <- min(max(w[m], 1e-300), 1 - 1e-16)
      ww[m, ] <- plogis(qlogis(wm) + c(-1, 1) * z * se_logit)
    }
  } else {
    if (method == "mc") {
      ev <- eigen((Sigma_w + t(Sigma_w)) / 2, symmetric = TRUE)
      A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), M)
      Z <- matrix(rnorm(nsim * M), nsim, M)
      eta_sim <- sweep(Z %*% t(A), 2, eta, "+")
    } else {
      eta_sim <- eta_b
    }
    eta_sim <- eta_sim - apply(eta_sim, 1, max)
    w_sim <- exp(eta_sim); w_sim <- w_sim / rowSums(w_sim)
    ww <- t(apply(w_sim, 2, quantile, probs = probs, names = FALSE))
  }
  weight_tab <- data.frame(model = models, weight = w,
                           lower = ww[, 1], upper = ww[, 2],
                           CI_method = rep(method, M),
                           row.names = NULL, stringsAsFactors = FALSE)
  names(weight_tab)[2] <- paste0(type, ".weights")

  out <- list(type = type, type_ev = type_ev, level = level, method = method,
              n_studies = S, B = if (use_boot) B else NULL, nsim = if (method == "mc") nsim else NULL,
              penalty_factor = pf, priorICweights = prior, study_weights = sw,
              Sigma_values = Sigma_v, Sigma_eta = Sigma_w,
              IC = IC_tab, difference = diff_tab,
              weights = weight_tab, ratio_weights = ratio_tab,
              best = models[which.max(w)],
              any_same_ormle = any(n_same > 0),
              weighted = length(unique(round(a, 12))) > 1)
  class(out) <- "evSyn_CI"
  out
}


print.evSyn_CI <- function(x, digits = max(3, getOption("digits") - 4),
                           all_pairs = FALSE, ...) {
  type_label <- toupper(x$type)
  pct <- paste0(format(100 * x$level), "%")
  cat(sprintf("restriktor: %s confidence intervals for the pooled (final) %s values, differences, weights, and ratios of weights\n",
              pct, type_label))
  cat(sprintf("(%s evidence synthesis of %s studies; %s)\n", x$type_ev, x$n_studies,
              switch(x$method,
                     boot  = paste0("parametric-bootstrap percentile CIs, B = ", x$B,
                                    " draws from N(estimates, VCOV) per study, independently over the studies",
                                    if (x$type %in% c("goric", "goricc")) 
                                      paste0("; CIs of the pooled ", type_label, " values: casewise approach") else ""),
                     delta = "delta method based on the estimates and their covariance matrices; weights: delta method on the logit scale",
                     mc    = paste0("delta method based on the estimates and their covariance matrices; weights: Monte Carlo (", x$nsim, " draws)"))))

  fmt <- function(df) {
    num <- vapply(df, function(v) is.numeric(v) && !is.integer(v), logical(1))
    df[num] <- lapply(df[num], function(v) formatC(v, digits = digits, format = "f"))
    df[] <- lapply(df, function(v) { v[v == "NA"] <- "-"; v })
    df
  }

  cat(sprintf("\nPooled %s values:\n", type_label))
  print(fmt(x$IC), row.names = FALSE, right = TRUE)

  sel <- if (all_pairs) rep(TRUE, nrow(x$difference)) else x$difference$model == x$best
  cat(sprintf("\nDifferences in pooled %s values%s:\n", type_label,
              ifelse(all_pairs, "", paste0(" (best hypothesis ", x$best, " versus the others)"))))
  print(fmt(x$difference[sel, , drop = FALSE]), row.names = FALSE, right = TRUE)

  cat(sprintf("\nPooled %s weights:\n", type_label))
  print(fmt(x$weights), row.names = FALSE, right = TRUE)

  cat(sprintf("\nRatios of pooled %s weights%s:\n", type_label,
              ifelse(all_pairs, "", paste0(" (best hypothesis ", x$best, " versus the others)"))))
  print(fmt(x$ratio_weights[sel, , drop = FALSE]), row.names = FALSE, right = TRUE)

  cat(paste0("\n(n_same_ormle: number of studies in which the two hypotheses have the same ",
             "order-restricted MLEs (or-mle's), and thus the same fit.)\n"))
  if (isTRUE(x$weighted)) {
    message("\nrestriktor Note: The study weights are used in the pooled ", type_label,
            " weights and ratios of weights (and their CIs), not in the pooled ", type_label,
            " values and their differences (as in evSyn).")
  }
  if (isTRUE(x$any_same_ormle) && x$method != "boot") {
    message("\nrestriktor Note: In some studies, hypotheses have the same order-restricted MLEs ",
            "(or-mle's) (n_same_ormle > 0), for example, a hypothesis and the unconstrained hypothesis ",
            "when the estimates are in agreement with that hypothesis. The differences in ", type_label,
            " values are then not normally distributed; the normal-based CIs ignore this (order 1) ",
            "variability for these studies and are too narrow (no CI is given if this holds in all studies). ",
            "Use method = \"boot\" instead.")
  }
  if (x$type %in% c("gorica", "goricac") && x$method != "boot") {
    message("\nrestriktor Note: The CI of a single pooled ", type_label, " value is based on ",
            "the squared Mahalanobis distances to the hypothesis; its SE is 0 for hypotheses ",
            "in agreement with the estimates in all studies. Only the relative values ",
            "(differences, weights, ratios) are of interest anyway.")
  }
  invisible(x)
}
