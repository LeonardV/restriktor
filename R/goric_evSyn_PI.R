# =============================================================================
# Prediction interval (PI) for the ratio of IC weights in a new study
# =============================================================================
# The frequentist counterpart of the 'Bayes factor for the next person'
# (Klaassen et al., 2018): the interval in which the ratio of IC weights of the
# overall best hypothesis versus another hypothesis is expected to lie in a
# new (comparable) study. It summarizes the heterogeneity of the studies in an
# interpretable way.
#
# Model (added / equal approach; independent studies):
#   log(w_b / w_k) in study s = d_s + p_s, with
#     d_s = LL_b - LL_k                  (difference in the maximum
#                                          order-restricted log likelihoods)
#     p_s = -(c / 2) * (PT_b - PT_k)     (penalty part; c = penalty factor)
#   The log-likelihood difference grows (roughly linearly) with the sample size
#   n_s, while the penalty part does not. Therefore, the per-observation
#   difference z_s = d_s / n_s is modelled with a random-effects model:
#     z_s ~ N(mu, v_s + tau^2),   v_s = (SE(d_s) / n_s)^2,
#   with SE(d_s) the bootstrap standard deviation of the log-likelihood
#   difference in study s (as in goric_CI()), and tau^2 estimated with REML.
#   For a new study of size n_new (by default the median sample size of the
#   studies), the PI of the (observed) log ratio is
#     n_new * (mu^ +- t_{S-2} * sqrt(SE(mu^)^2 + tau^2 + v_new)) + p_new,
#   with v_new = median(n_s * v_s) / n_new (the sampling variance in the new
#   study; omitted for type = "true") and p_new the median penalty part
#   (cf. Higgins, Thompson, & Spiegelhalter, 2009).
#
# Ceiling: when the best hypothesis and the other one have the same maximum
# log likelihood in a study (d_s = 0; the hypotheses overlap there and, when
# nested, the ratio is at its bound), d_s lies at its bound and its sampling
# distribution is not normal. Then, a two-part model is used: (1) the
# probability that a new study is at the ceiling, estimated by the proportion
# of studies at the ceiling, and (2) the random-effects model above for the
# other studies (with t_{k-2} quantiles, k the number of these studies; at
# least 1 df), truncated at the bound. The PI consists of the quantiles of
# this mixture. (A censored, Tobit-type, random-effects model was tried as
# well, but the bootstrap SEs of studies at the ceiling are very small, which
# makes it unstable and dominated by those studies.)
#
# When the sample sizes are unknown, the model is fitted to d_s itself (i.e.,
# n_s = n_new = 1), which assumes studies of comparable size.
#
# Average approach (analyses of the same data): the PI regards another
# analysis of the same data: mean(y) +- t_{S-1} * sqrt((1 + 1/S) * B), with y
# the analysis-specific log ratios and B their (between-analysis) variance.
# The sampling variance is shared by all analyses (same data) and is not part
# of this interval.
#
# The study-specific ratios do not include the prior IC weights (as the
# study-specific weights in evSyn); unless priorICweights_study = TRUE (which
# is not active for now).
# =============================================================================

evSyn_PI <- function(object, vs = NULL, best = NULL, level = 0.95, n_new = NULL,
                     type = c("observed", "true"), B = 1000L, seed = NULL, tol = 1e-8) {

  if (!inherits(object, "evSyn")) {
    stop("\nrestriktor ERROR: The object should be of class 'evSyn' ",
         "(i.e., output from the evSyn() function).", call. = FALSE)
  }
  type <- match.arg(type)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("\nrestriktor ERROR: The argument 'level' should be a single number in between 0 and 1.",
         call. = FALSE)
  }
  type_ev <- object$type_ev
  models  <- colnames(object$LL_m)
  if (is.null(models)) models <- colnames(object$Cumulative_GORICA_weights)
  M <- length(models)
  if (M < 2L) {
    stop("\nrestriktor ERROR: A prediction interval for a ratio of IC weights requires at least two ",
         "hypotheses.", call. = FALSE)
  }
  final_w <- object$Cumulative_GORICA_weights["Final", ]
  b <- if (is.null(best)) which.max(final_w) else match(best, models)
  if (is.na(b)) stop("\nrestriktor ERROR: 'best' should be one of: ", paste(models, collapse = ", "), ".",
                     call. = FALSE)
  if (is.null(vs)) vs <- models[-b]
  kk_all <- match(vs, models)
  if (anyNA(kk_all) || any(kk_all == b)) {
    stop("\nrestriktor ERROR: 'vs' should contain hypotheses other than ", sQuote(models[b]), ".",
         call. = FALSE)
  }
  pf    <- if (is.null(object$penalty_factor)) 2 else object$penalty_factor
  prior <- object$priorICweights
  if (is.null(prior) || length(prior) != M) prior <- rep(1 / M, M)
  prior <- prior / sum(prior)
  lp_study <- if (isTRUE(object$priorICweights_study)) log(prior) else rep(0, M)
  S     <- object$n_studies
  units <- .evSyn_unit_label(type_ev, plural = TRUE, capital = FALSE)
  unit1 <- .evSyn_unit_label(type_ev, capital = FALSE)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)

  LL_m <- object$LL_m; PT_m <- object$PT_m
  if (is.null(LL_m) || is.null(PT_m)) {
    stop("\nrestriktor ERROR: A prediction interval requires the ", unit1, "-specific log-likelihood ",
         "and penalty values (i.e., input of evSyn() consisting of estimates and their covariance ",
         "matrices, of goric objects, or of log-likelihood and penalty values).", call. = FALSE)
  }
  LL_m <- as.matrix(LL_m); PT_m <- as.matrix(PT_m)

  res <- vector("list", length(kk_all))
  # -------------------------------------------------------------------------
  # Average approach: another analysis of the same data
  # -------------------------------------------------------------------------
  if (type_ev == "average") {
    if (S < 2L) stop("\nrestriktor ERROR: At least two analyses are needed.", call. = FALSE)
    for (i in seq_along(kk_all)) {
      k <- kk_all[i]
      d <- LL_m[, b] - LL_m[, k]
      p <- -(pf / 2) * (PT_m[, b] - PT_m[, k]) + (lp_study[b] - lp_study[k])
      y <- d + p
      same <- abs(d) < tol * max(1, abs(LL_m[, c(b, k)]))
      Bv <- stats::var(y)
      sd_pred <- sqrt((1 + 1 / S) * Bv)
      tq <- stats::qt(probs, df = S - 1)
      lim <- mean(y) + tq * sd_pred
      bound <- .evSyn_PI_bound(same, d, tol)
      if (!is.na(bound$side)) {
        cap <- stats::median(p[same])
        lim <- if (bound$side == "upper") pmin(lim, cap) else pmax(lim, cap)
      }
      res[[i]] <- data.frame(model = models[b], vs = models[k],
                             ratio_typical = exp(mean(y)), lower = exp(lim[1]), upper = exp(lim[2]),
                             P_favours_best = if (sd_pred > 0) 1 - stats::pt(-mean(y) / sd_pred, df = S - 1) else as.numeric(mean(y) > 0),
                             n_ceiling = sum(same), P_ceiling = mean(same),
                             tau = sqrt(Bv), n_new = NA_real_, k = S, df = S - 1,
                             method = "between-analysis variance", stringsAsFactors = FALSE)
    }
    out <- list(PI = do.call(rbind, res), type_ev = type_ev, level = level, type = "another analysis",
                n_studies = S, per_obs = FALSE, best = models[b], ic_type = object$type)
    class(out) <- "evSyn_PI"
    return(out)
  }

  # -------------------------------------------------------------------------
  # Added / equal approach: a new (independent) study
  # -------------------------------------------------------------------------
  goric_list <- object$goric_studies
  if (is.null(goric_list)) {
    stop("\nrestriktor ERROR: A prediction interval requires the ", unit1, "-specific sampling ",
         "variability, which is only available when the input of evSyn() consists of (standardized) ",
         "estimates and their covariance matrices (input type 'est_vcov' or 'escalc') or of goric ",
         "objects (input type 'gorica').", call. = FALSE)
  }
  cl_info <- object$cluster_info
  if (S < 3L) {
    stop("\nrestriktor ERROR: A prediction interval requires at least three ",
         if (!is.null(cl_info) && identical(cl_info$dependence, "cluster")) "clusters" else units, ".",
         call. = FALSE)
  }
  S <- length(goric_list)   # number of studies (also with clusters of studies)
  if (!is.null(seed)) set.seed(seed)
  # the stored LL_m / PT_m follow the (possibly re-)ordered studies; the goric
  # objects are stored in the input order
  ord <- object$order_studies
  if (is.null(ord) || length(ord) != S) ord <- seq_len(S)
  LL_in <- t(vapply(goric_list, function(g) g$result$loglik, numeric(M)))
  PT_in <- t(vapply(goric_list, function(g) g$result$penalty, numeric(M)))

  # sample sizes
  n_s <- vapply(goric_list, function(g) { n <- g$sample_nobs; if (is.null(n)) NA_real_ else sum(n) },
                numeric(1))
  if (anyNA(n_s) && !is.null(object$study_sample_nobs) && length(object$study_sample_nobs) == S) {
    n_s <- as.numeric(unlist(object$study_sample_nobs))
  }
  per_obs <- !anyNA(n_s) && all(n_s > 0)
  if (!per_obs) {
    n_s <- rep(1, S)
    if (!is.null(n_new)) {
      message("\nrestriktor Message: The sample sizes of the ", units, " are unknown; 'n_new' is ignored.")
    }
    n_new <- 1
  } else if (is.null(n_new)) {
    n_new <- stats::median(n_s)
  }

  # bootstrap SDs of the log-likelihood differences (per study)
  models_g <- lapply(goric_list, function(g) as.character(g$result$model))
  IC_b <- lapply(seq_len(S), function(s) {
    unname(goric_CI_boot(goric_list[[s]], models_g[[s]], B = B, penalty_factor = pf))
  })

  tq <- stats::qt(probs, df = S - 2)
  het <- vector("list", length(kk_all))
  n_med <- if (per_obs) stats::median(n_s) else 1   # median sample size (for I^2 and tau on the ratio scale)
  for (i in seq_along(kk_all)) {
    k <- kk_all[i]
    d  <- LL_in[, b] - LL_in[, k]
    p  <- -(pf / 2) * (PT_in[, b] - PT_in[, k]) + (lp_study[b] - lp_study[k])
    # bootstrap draws of the log-likelihood difference (the draws of the IC values include the
    # penalties, which are fixed: -(IC_b - IC_k) / 2 = d + p_pen)
    p_pen <- -(pf / 2) * (PT_in[, b] - PT_in[, k])
    dB <- vapply(seq_len(S), function(s) -(IC_b[[s]][, b] - IC_b[[s]][, k]) / 2 - p_pen[s], numeric(B))   # B x S
    sd_d <- apply(dB, 2, stats::sd)
    same <- abs(d) < tol * pmax(1, abs(LL_in[, b]))
    bound <- .evSyn_PI_bound(same, dB, tol)
    n_u <- n_s
    if (!is.null(cl_info)) {
      # Dependent studies: the units are the clusters (a new cluster, or a new 
      # independent study): (weighted) means of d, p, and n within a cluster; the 
      # studies within a cluster are treated as fully dependent, so the SD of the 
      # mean is the weighted mean of the SDs (comonotone).
      cl <- cl_info$cluster; levs <- cl_info$levels; sw_st <- cl_info$study_weights
      if (is.null(sw_st)) sw_st <- rep(1, S)
      agg <- t(vapply(levs, function(lv) {
        idx <- which(cl == lv); ww <- sw_st[idx] / sum(sw_st[idx])
        c(sum(ww * d[idx]), sum(ww * p[idx]), sum(ww * sd_d[idx]), sum(ww * n_s[idx]), all(same[idx]))
      }, numeric(5)))
      d <- agg[, 1]; p <- agg[, 2]; sd_d <- agg[, 3]; n_u <- agg[, 4]; same <- agg[, 5] == 1
    }
    z <- d / n_u
    v <- (sd_d / n_u)^2
    v[v <= 0] <- min(v[v > 0], 1e-12)
    side <- bound$side
    p_new <- stats::median(p)
    # Two-part model: (1) the probability that a new study is at the ceiling 
    # (estimated by the proportion of studies at the ceiling) and (2) a 
    # random-effects model for the other studies (truncated at the bound, if any).
    nc <- !same
    k_nc <- sum(nc)
    pi_c <- if (is.na(side)) 0 else mean(same)
    if (k_nc < 2L && !(k_nc == 1L && !is.na(side))) {
      at_c <- all(same)
      res[[i]] <- data.frame(model = models[b], vs = models[k],
                             ratio_typical = if (at_c) exp(p_new) else NA_real_,
                             lower = if (at_c) exp(p_new) else NA_real_,
                             upper = if (at_c) exp(p_new) else NA_real_,
                             P_favours_best = if (at_c) as.numeric(p_new > 0) else NA_real_,
                             n_ceiling = sum(same), P_ceiling = if (at_c) 1 else pi_c, tau = NA_real_,
                             n_new = n_new, k = k_nc, df = NA_real_,
                             method = if (at_c) "all at the ceiling" else "not estimable",
                             stringsAsFactors = FALSE)
      next
    }
    fit <- if (k_nc >= 2L) {
      .evSyn_PI_reml(z[nc], v[nc])
    } else {
      # only one study not at the ceiling: tau cannot be estimated (set to 0)
      list(mu = z[nc], se_mu = sqrt(v[nc]), tau = 0, method = "REML")
    }
    df  <- max(k_nc - 2L, 1L)
    v_new <- if (type == "observed") stats::median(n_u[nc] * v[nc]) / n_new else 0
    sd_pred <- sqrt(fit$se_mu^2 + fit$tau^2 + v_new)
    Gt  <- function(x) stats::pt((x / n_new - fit$mu) / sd_pred, df = df)       # CDF of d_new (not truncated)
    Gti <- function(q) n_new * (fit$mu + sd_pred * stats::qt(q, df = df))     # its quantile function
    G0  <- Gt(0)
    # quantile function and CDF of d_new (mixture of the point mass at the ceiling and the truncated part)
    Q <- function(pr) {
      if (is.na(side)) return(Gti(pr))
      if (side == "upper") {   # d <= 0; ceiling (d = 0) is the maximum
        if (pr >= 1 - pi_c) 0 else Gti(pr / (1 - pi_c) * G0)
      } else {                 # d >= 0; ceiling (d = 0) is the minimum
        if (pr <= pi_c) 0 else Gti(G0 + (pr - pi_c) / (1 - pi_c) * (1 - G0))
      }
    }
    Fd <- function(x) {
      if (is.na(side)) return(Gt(x))
      if (side == "upper") {
        if (x >= 0) 1 else (1 - pi_c) * Gt(x) / G0
      } else {
        if (x < 0) 0 else pi_c + (1 - pi_c) * (Gt(x) - G0) / (1 - G0)
      }
    }
    lim_d <- c(Q(probs[1]), Q(probs[2]))
    typ_d <- Q(0.5)
    P_fav <- 1 - Fd(-p_new)   # probability that a new study favours the best hypothesis (log ratio > 0)
    # Heterogeneity (studies not at the ceiling), on the per-observation scale:
    # - tau, and tau translated to the log-ratio scale for a study of median size (n_med * tau);
    # - I^2 for a study of median size: tau^2 / (tau^2 + v_med), with v_med the 
    #   typical within-study variance at the median sample size;
    # - Cochran's Q (fixed-effect weights 1 / v_s), with k - 1 df.
    v_med <- stats::median(n_u[nc] * v[nc]) / n_med
    I2 <- if (k_nc >= 2L) fit$tau^2 / (fit$tau^2 + v_med) else NA_real_
    if (k_nc >= 2L) {
      w0 <- 1 / v[nc]; zbar <- sum(w0 * z[nc]) / sum(w0)
      Qs <- sum(w0 * (z[nc] - zbar)^2); Qdf <- as.integer(k_nc - 1L)
      Qp <- stats::pchisq(Qs, df = Qdf, lower.tail = FALSE)
    } else {
      Qs <- Qp <- NA_real_; Qdf <- NA_integer_
    }
    het[[i]] <- data.frame(model = models[b], vs = models[k], k = k_nc, n_ceiling = sum(same),
                           tau = if (k_nc >= 2L) fit$tau else NA_real_,
                           tau_logratio = if (k_nc >= 2L) n_med * fit$tau else NA_real_,
                           I2 = I2, Q = Qs, Q_df = Qdf, Q_p = Qp, n_med = n_med,
                           stringsAsFactors = FALSE)
    res[[i]] <- data.frame(model = models[b], vs = models[k],
                           ratio_typical = exp(typ_d + p_new),
                           lower = exp(lim_d[1] + p_new), upper = exp(lim_d[2] + p_new),
                           P_favours_best = P_fav,
                           n_ceiling = sum(same), P_ceiling = pi_c,
                           tau = fit$tau, n_new = n_new, k = k_nc, df = df,
                           method = if (k_nc == 1L) "ceiling (+ 1 study; tau = 0)" else
                             if (any(same) && !is.na(side)) "REML + ceiling" else "REML",
                           stringsAsFactors = FALSE)
  }
  out <- list(PI = do.call(rbind, res), type_ev = type_ev, level = level, type = type,
              n_studies = if (!is.null(cl_info)) cl_info$n_clusters else S, per_obs = per_obs,
              clusters = !is.null(cl_info), best = models[b], B = B, ic_type = object$type,
              sample_nobs = if (per_obs) n_s else NULL,
              heterogeneity = do.call(rbind, het[!vapply(het, is.null, logical(1))]),
              heterogeneity_sd = tryCatch(
                .evSyn_signed_distance_het(goric_list, best = models[b], models = models, n_s = n_s,
                                           per_obs = per_obs, cl_info = cl_info, B = B, pf = pf, tol = tol),
                error = function(e) NULL))
  class(out) <- "evSyn_PI"
  out
}


# Side of the bound (ceiling) of the log-likelihood difference d = LL_b - LL_k:
# "upper" when d <= 0 in all (bootstrap) draws (best hypothesis nested in the
# other one; at the ceiling, d = 0 is the maximum), "lower" when d >= 0 in all
# draws, and NA when there is no bound.
.evSyn_PI_bound <- function(same, draws, tol) {
  scale <- max(1, abs(draws), na.rm = TRUE)
  if (all(draws <= tol * scale, na.rm = TRUE)) return(list(side = "upper"))
  if (all(draws >= -tol * scale, na.rm = TRUE)) return(list(side = "lower"))
  list(side = NA_character_)
}

# Heterogeneity based on the signed distance of the estimates to the best 
# hypothesis, which uses all studies (also those at the ceiling).
# The ceiling of the log-likelihood difference arises because, for estimates 
# inside a hypothesis, its log likelihood equals that of the unconstrained 
# hypothesis, whether the estimates lie just inside or far inside. The signed 
# (Mahalanobis) distance keeps this information: positive inside the best 
# hypothesis (the depth, i.e., the distance to its complement) and negative 
# outside (minus the distance to the hypothesis):
#   s = sign(D) * sqrt(2 |D|),  D = LL(H_b) - LL(complement of H_b)
# (for a hypothesis with only equality constraints, the complement is the 
# unconstrained hypothesis, so that s <= 0). It is continuous, without a bound, 
# and roughly normally distributed with SD about 1 (its SE is obtained with the 
# parametric bootstrap); s / sqrt(n) is a standardized margin per observation 
# (effect-size like), on which the sample size cancels out. tau, I^2 (for a 
# study of median sample size), and Cochran's Q are computed on this scale, 
# using all studies.
.evSyn_signed_distance_het <- function(goric_list, best, models, n_s, per_obs, cl_info = NULL,
                                       B = 1000L, pf = 2, tol = 1e-8) {
  S <- length(goric_list)
  hyps <- goric_list[[1]]$hypotheses_usr
  sgn <- 1
  hb <- hyps[[best]]
  if (is.null(hb)) {
    # best hypothesis is a failsafe hypothesis
    if (tolower(best) == "complement" && length(hyps) == 1L) {
      hb <- hyps[[1]]; sgn <- -1          # distance to the complement = minus that to H1
    } else {
      return(NULL)                        # unconstrained: no informative hypothesis
    }
  }
  sd_s <- s_obs <- numeric(S)
  for (s in seq_len(S)) {
    g <- goric_list[[s]]
    gc <- suppressMessages(suppressWarnings(
      goric(g$b.unrestr, VCOV = g$VCOV, hypotheses = list(Hb = hb), comparison = "complement",
            type = "gorica", penalty_factor = pf)))
    LL <- gc$result$loglik; PT <- gc$result$penalty
    D  <- LL[1] - LL[2]
    s_obs[s] <- sgn * sign(D) * sqrt(2 * abs(D))
    ICb <- unname(goric_CI_boot(gc, as.character(gc$result$model), B = B, penalty_factor = pf))
    Db <- -(ICb[, 1] - ICb[, 2]) / 2 + (pf / 2) * (PT[1] - PT[2])
    sd_s[s] <- stats::sd(sgn * sign(Db) * sqrt(2 * abs(Db)))
  }
  n_u <- if (per_obs) n_s else rep(1, S)
  if (!is.null(cl_info)) {
    cl <- cl_info$cluster; levs <- cl_info$levels; sw_st <- cl_info$study_weights
    if (is.null(sw_st)) sw_st <- rep(1, S)
    agg <- t(vapply(levs, function(lv) {
      idx <- which(cl == lv); ww <- sw_st[idx] / sum(sw_st[idx])
      c(sum(ww * s_obs[idx]), sum(ww * sd_s[idx]), sum(ww * n_u[idx]))
    }, numeric(3)))
    s_obs <- agg[, 1]; sd_s <- agg[, 2]; n_u <- agg[, 3]
  }
  k <- length(s_obs)
  u <- s_obs / sqrt(n_u)                     # per-observation (standardized) margin
  v <- sd_s^2 / n_u
  v[v <= 0] <- min(c(v[v > 0], 1e-12))
  fit <- .evSyn_PI_reml(u, v)
  n_med <- if (per_obs) stats::median(n_s) else 1
  v_med <- stats::median(n_u * v) / n_med
  w0 <- 1 / v; ubar <- sum(w0 * u) / sum(w0)
  Qs <- sum(w0 * (u - ubar)^2)
  data.frame(hypothesis = best, k = k, n_inside = sum(s_obs > tol), n_outside = sum(s_obs < -tol),
             mean_margin = fit$mu, tau = fit$tau, tau_distance = sqrt(n_med) * fit$tau,
             I2 = fit$tau^2 / (fit$tau^2 + v_med), Q = Qs, Q_df = as.integer(k - 1L),
             Q_p = stats::pchisq(Qs, df = k - 1, lower.tail = FALSE), n_med = n_med,
             stringsAsFactors = FALSE)
}

# Random-effects model with REML estimate of tau^2
.evSyn_PI_reml <- function(z, v) {
  if (length(z) < 2L) return(NULL)
  rl <- function(t2) {
    w <- 1 / (v + t2); mu <- sum(w * z) / sum(w)
    -0.5 * (sum(log(v + t2)) + log(sum(w)) + sum(w * (z - mu)^2))
  }
  upper <- max(10 * stats::var(z), 10 * max(v), 1e-8)
  t2 <- stats::optimize(rl, c(0, upper), maximum = TRUE)$maximum
  if (rl(0) >= rl(t2)) t2 <- 0
  w <- 1 / (v + t2)
  list(mu = sum(w * z) / sum(w), se_mu = sqrt(1 / sum(w)), tau = sqrt(t2), method = "REML")
}

print.evSyn_PI <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  type_label <- if (!is.null(x$ic_type)) toupper(x$ic_type) else "GORICA"
  pct <- paste0(format(100 * x$level), "%")
  units <- .evSyn_unit_label(x$type_ev, plural = TRUE, capital = FALSE)
  unit1 <- .evSyn_unit_label(x$type_ev, capital = FALSE)
  pi <- x$PI
  if (x$type_ev == "average") {
    cat(sprintf("restriktor: %s prediction interval for the ratio of %s weights in another analysis of the same data\n\n",
                pct, type_label))
  } else {
    cat(sprintf("restriktor: %s prediction interval for the ratio of %s weights in a new %s\n\n",
                pct, type_label, unit1))
  }
  tab <- data.frame(comparison = paste0(pi$model, " vs. ", pi$vs),
                    typical = pi$ratio_typical, lower = pi$lower, upper = pi$upper,
                    "P(favours best)" = pi$P_favours_best,
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (any(pi$n_ceiling > 0)) {
    tab[["at ceiling"]] <- pi$n_ceiling
    tab[["P(ceiling)"]] <- pi$P_ceiling
  }
  tab[["method"]] <- pi$method
  print(.format_CI_table(tab, digits), row.names = FALSE, right = TRUE)
  cat("\n")
  if (x$type_ev == "average") {
    cat(paste0("typical: ratio based on the mean of the analysis-specific log ratios; lower, upper: bounds of ",
               "the prediction interval for another analysis of the same data, based on the between-analysis ",
               "variance (the sampling variance is shared by all analyses and is not part of this interval); ",
               "P(favours best): probability that such an analysis favours ", sQuote(x$best), ".\n"))
  } else {
    cat(paste0("typical: ratio expected in a typical new ", unit1, "; lower, upper: bounds of the prediction ",
               "interval for the ", if (x$type == "observed") "observed " else "true ",
               "ratio in a new ", unit1,
               if (x$per_obs) paste0(" of size n = ", format(pi$n_new[1]), " (median sample size of the ", units,
                                     ", unless specified otherwise)") else "",
               "; P(favours best): probability that a new ", unit1, " favours ", sQuote(x$best),
               " (ratio > 1).\n"))
    cat(paste0("Based on a random-effects model for the per-observation difference in maximum order-restricted ",
               "log likelihoods (i.e., loglik in goric() and LL in evSyn()), with the ", unit1,
               "-specific standard errors from the parametric bootstrap (as in goric_CI()) and ",
               "t(S - 2) quantiles.\n"))
    if (isTRUE(x$clusters)) {
      cat(paste0("With clusters of dependent studies, the (", x$n_studies, ") clusters are the units: the ",
                 "cluster-specific evidence is the (weighted) mean over its studies, which are treated as fully ",
                 "dependent.\n"))
    }
    if (!x$per_obs) {
      message("\nrestriktor Note: The sample sizes of the ", units, " are unknown (specify 'study_sample_nobs' ",
              "in evSyn() or 'sample_nobs' in goric()); the prediction interval assumes ", units,
              " of comparable size.")
    }
    if (x$n_studies < 10) {
      message("\nrestriktor Note: With ", x$n_studies, " ", units, ", the heterogeneity (tau) is estimated ",
              "imprecisely, and so is the prediction interval.")
    }
  }
  # Heterogeneity measures (added / equal approach)
  het <- x$heterogeneity
  if (!is.null(het) && nrow(het) > 0) {
    nm <- het$n_med[1]
    cat(paste0("\nHeterogeneity of the ", unit1, "-specific evidence (per-observation difference in maximum ",
               "order-restricted log likelihoods", if (any(het$n_ceiling > 0)) paste0("; ", units,
               " not at the ceiling") else "", "):\n"))
    ht <- data.frame(comparison = paste0(het$model, " vs. ", het$vs), k = het$k,
                     tau = het$tau, "tau (log ratio)" = het$tau_logratio,
                     "range (factor)" = exp(stats::qnorm(0.975) * het$tau_logratio),
                     "I2 (%)" = 100 * het$I2, Q = het$Q, df = het$Q_df, "p(Q)" = het$Q_p,
                     check.names = FALSE, stringsAsFactors = FALSE)
    print(.format_CI_table(ht, digits), row.names = FALSE, right = TRUE)
    cat("\n")
    cat(paste0("k: number of ", units, " used; tau: between-", unit1, " standard deviation of the per-observation ",
               "evidence; tau (log ratio): the same for the log ratio of weights in a ", unit1,
               if (x$per_obs) paste0(" of median sample size (n = ", format(nm), ")") else "",
               ", i.e., n x tau; range (factor): the true ratios of ", units,
               if (x$per_obs) paste0(" of median sample size (n = ", format(nm), ")") else "",
               " differ by about this factor (up or down) from the typical ratio (95% of them); I2: percentage of ",
               "the variance in the ", unit1, "-specific evidence that is due to true differences between the ",
               units, " rather than sampling error, for a ", unit1,
               if (x$per_obs) paste0(" of median sample size (n = ", format(nm), ")") else " of typical size",
               "; Q: Cochran's test of homogeneity (with few ", units, " it has little power, so a non-significant ",
               "Q does not show that the ", units, " agree).\n"))
    hs <- x$heterogeneity_sd
    if (!is.null(hs) && nrow(hs) > 0) {
      cat(paste0("\nHeterogeneity based on the signed distance of the estimates to ", sQuote(hs$hypothesis[1]),
                 " (all ", units, "):\n"))
      hst <- data.frame(hypothesis = hs$hypothesis, k = hs$k, inside = hs$n_inside, outside = hs$n_outside,
                        tau = hs$tau, "tau (distance)" = hs$tau_distance, "I2 (%)" = 100 * hs$I2,
                        Q = hs$Q, df = hs$Q_df, "p(Q)" = hs$Q_p, check.names = FALSE, stringsAsFactors = FALSE)
      print(.format_CI_table(hst, digits), row.names = FALSE, right = TRUE)
      cat("\n")
      cat(paste0("Signed distance: the (Mahalanobis) distance of the estimates to the boundary of ",
                 sQuote(hs$hypothesis[1]), ", positive inside (depth) and negative outside, divided by sqrt(n) ",
                 "(a standardized margin per observation, on which the sample size cancels out); inside, outside: ",
                 "number of ", units, " with estimates inside or outside ", sQuote(hs$hypothesis[1]),
                 "; tau (distance): the between-", unit1, " SD of the signed distance for a ", unit1,
                 if (x$per_obs) paste0(" of median sample size (n = ", format(hs$n_med[1]), ")") else "",
                 ", in standard-error units; I2: for a ", unit1,
                 if (x$per_obs) paste0(" of median sample size (n = ", format(hs$n_med[1]), ")") else " of typical size",
                 ". Unlike the table above, this uses all ", units, ", also those at the ceiling: there, the ",
                 "log-likelihood difference is 0 whether the estimates lie just inside or far inside the ",
                 "hypothesis, whereas the signed distance distinguishes these. It regards how clearly the ",
                 "estimates lie inside (or outside) ", sQuote(hs$hypothesis[1]),
                 ", not the ratio of weights versus a specific other hypothesis.\n"))
    }
    if (any(het$n_ceiling > 0)) {
      message("\nrestriktor Note: The heterogeneity measures in the first table (and the prediction interval) ",
              "are based on the ", units, " not at the ceiling only (k); with many ", units, " at the ceiling, ",
              "these measures can understate the heterogeneity. The measures based on the signed distance use ",
              "all ", units, ".")
    }
    cat(paste0("Notably, this regards heterogeneity in the evidence (which combines direction and size of the ",
               "effects); for heterogeneity in the parameters themselves, a (multivariate) meta-analysis of the ",
               "estimates (e.g., with the metafor package) is the appropriate tool.\n"))
  }
  if (any(pi$n_ceiling > 0)) {
    cat(paste0("at ceiling: number of ", units, " in which both hypotheses have the same maximum order-restricted ",
               "log likelihood (the ratio is then at its bound, if nested); P(ceiling): probability that a new ",
               unit1, " is at the ceiling. The prediction interval is capped at the ceiling; ",
               "a two-part model is used: the proportion of ", units, " at the ceiling and a random-effects model for the other ", units, ".\n"))
  }
  invisible(x)
}
