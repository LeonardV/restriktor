# =============================================================================
# Confidence intervals (CIs) for GORIC(A) values, differences in GORIC(A)
# values, GORIC(A) weights, and ratios of GORIC(A) weights.
# =============================================================================
#
# Two types of CIs are available:
# - Analytic (normal-based) CIs (method = "delta" or "mc"), based on the
#   within-sample variability of the individual (casewise) log-likelihood
#   contributions, as proposed for AIC/BIC (differences) by Merkle, You, and
#   Preacher (2016; building on Vuong, 1989), and as used for the standard
#   errors of WAIC/LOO (differences) by Vehtari, Gelman, and Gabry (2017).
#   For the GORICA, which is not a sum over cases, the analogous delta-method
#   result in terms of the (asymptotically normal) estimates is used.
# - Parametric-bootstrap percentile CIs (method = "boot"), based on draws
#   theta* ~ N(theta^, V), recommended when hypotheses are close or
#   (nearly) indistinguishable (see below).
#
# -----------------------------------------------------------------------------
# Notation
# -----------------------------------------------------------------------------
# Hypotheses m = 1, ..., M (including the unconstrained hypothesis / the
# complement / Heq, if in the set).
#   IC_m   = -2 * loglik_m + c * PT_m        (c = penalty_factor, default 2)
#   eta_m  = log(p_m) - IC_m / 2             (p_m = priorICweights)
#   w_m    = exp(eta_m) / sum_k exp(eta_k)   (GORIC(A) weights)
#   w_m / w_k = exp(eta_m - eta_k)           (ratio of GORIC(A) weights)
#
# The penalty PT_m depends on the constraints and the level probabilities
# (based on the design / covariance matrix), not on the outcome values; it is
# treated as fixed. The priors are fixed as well. Hence, all sampling
# variability comes from the log-likelihood values, and
#   Cov(eta) = Cov(loglik) =: Sigma   (an M x M matrix),
#   Var(IC_m) = 4 * Sigma[m, m].
#
# (A) GORIC(C): loglik_m = sum_i l_im, with l_im = log f(y_i | theta_m~) the
#     contribution of case i at the order-restricted MLE theta_m~ of H_m.
#     As in Vuong (1989) and Merkle et al. (2016), the difference
#     d_i = l_im - l_ik gives n^(-1/2) * (sum_i d_i - n E[d_i]) -> N(0, omega^2),
#     with omega^2 = Var(d_i) (estimating theta only adds O_p(1) terms).
#     More generally, Cov(loglik_m, loglik_k) is estimated by
#       Sigma[m, k] = sum_i (l_im - mean_m) (l_ik - mean_k) = n * S[m, k],
#     with S the sample covariance matrix (divisor n) of the casewise
#     contributions. Hence
#       SE(IC_m)        = 2 * sqrt(n) * s_m          (cf. Vehtari et al., 2017, Eq. 23)
#       SE(IC_m - IC_k) = 2 * sqrt(n) * omega_mk^    (Merkle et al., 2016, Eq. 16;
#                                                     Vehtari et al., 2017, Eq. 24)
#     (Vehtari et al. give these on the elpd scale, i.e., without the 2, and
#     with divisor n - 1 in the sample variance; Merkle et al., Eq. 9, use
#     divisor n, as here.) Note: the analytic CI for a single BIC in Preacher
#     and Merkle (2012) is based on inverting a noncentral chi-square
#     distribution instead; the casewise SE for a single value is the one of
#     Vehtari et al. (2017).
#     Merkle et al. (2016) first apply Vuong's (1989) formal test of
#     distinguishability (n * omega^2 versus a weighted sum of chi-squares).
#     That test assumes estimates that are not on the boundary of the
#     parameter space, which does not hold for order-restricted estimates;
#     therefore, here, pairs are only flagged when their SE is (numerically)
#     zero, and the bootstrap is recommended for close hypotheses.
#     Casewise contributions are available for:
#     - lm (also with prior weights) and rlm: normal density with the ML
#       variance of the residuals (as in con_loglik_lm()),
#     - mlm: multivariate normal density with the ML covariance matrix of the
#       residuals (note: goric() itself does not (yet) run for mlm objects),
#     - glm: gaussian, poisson, binomial, Gamma, and inverse.gaussian family
#       (with the ML/deviance-based dispersion as used in the glm AIC),
#     - lm with missing = "fiml": the observed-data (multivariate normal)
#       density of each case under the implied mean and covariance matrix of
#       the (restricted) FIML solution (for the complement: of the boundary
#       solution with the highest log-likelihood, as in goric()).
#
# (B) GORICA(C): loglik_m = -1/2 log|2 pi V| - 1/2 D_m, with
#       D_m = (theta^ - theta_m~)' V^{-1} (theta^ - theta_m~),
#     the squared Mahalanobis distance from the estimates theta^ to H_m, and
#     theta_m~ its projection onto H_m. Because D_m is the squared distance to
#     a closed convex set (cone), it is continuously differentiable with
#     gradient 2 V^{-1} (theta^ - theta_m~) (no derivative of theta_m~ itself
#     is needed). With theta^ ~ N(theta, V), the delta method gives, with
#     r_m = theta^ - theta_m~,
#       Sigma[m, k] = r_m' V^{-1} r_k        (i.e., Sigma = R V^{-1} R'),
#     and therefore
#       SE(IC_m)        = 2 * sqrt(D_m),
#       SE(IC_m - IC_k) = 2 * sqrt( (theta_m~ - theta_k~)' V^{-1} (theta_m~ - theta_k~) ).
#     (For the complement, D_c is the minimum over a finite union of convex
#     sets; the same gradient holds almost everywhere.)
#     For a linear model with known error variance, (A) and (B) coincide to
#     first order: the casewise d_i is then dominated by e_i x_i'(theta_m~ - theta_k~)/sigma^2,
#     whose summed variance is (theta_m~ - theta_k~)' (X'X/sigma^2) (theta_m~ - theta_k~),
#     with V = sigma^2 (X'X)^{-1}.
#
# Then, with z the (1 + level)/2 quantile of the standard normal:
#   1. IC value:         IC_m +/- z * 2 * sqrt(Sigma[m,m])
#   2. IC difference:    (IC_m - IC_k) +/- z * SE_mk,
#                        SE_mk = 2 * sqrt(Sigma[m,m] + Sigma[k,k] - 2 Sigma[m,k])
#   3. Ratio of weights: log(w_m / w_k) = log(p_m / p_k) - (IC_m - IC_k) / 2, so
#                        CI = exp( log(w_m / w_k) +/- z * SE_mk / 2 )
#                        (a monotone transformation of the CI in 2; asymmetric,
#                        always positive).
#   4. Weight:           on the logit scale (delta method):
#                        logit(w_m) = eta_m - log(sum_{k != m} exp(eta_k)),
#                        gradient g: g_m = 1, g_k = -w_k / (1 - w_m) (k != m),
#                        CI = plogis( logit(w_m) +/- z * sqrt(g' Sigma g) );
#                        or (method = "mc") by Monte Carlo: eta* ~ N(eta^, Sigma),
#                        w* = softmax(eta*), and taking the percentiles.
#
# (C) Parametric bootstrap (method = "boot"): draw theta*_b ~ N(theta^, V),
#     b = 1, ..., B, with V = object$VCOV; for each draw, project theta*_b onto
#     each hypothesis (quadratic programming; for the complement: the
#     closest boundary solution, as in goric()), giving D*_mb, and compute
#       IC*_mb = log|2 pi V| + D*_mb + c * PT_m   (PT_m fixed, as in the object).
#     The percentile CIs are the (1 - level)/2 and (1 + level)/2 quantiles of
#     IC*_mb - IC*_kb (differences), exp(eta*_mb - eta*_kb) (ratios), and
#     softmax(eta*_b)_m (weights). These draws are the same as those of the
#     'Observed' population in benchmark_asymp(), but here all hypotheses and
#     all pairs are kept (and no re-fitting of the hypotheses' penalties is
#     needed, which makes it fast).
#     For a GORIC(C) object, the bootstrap uses the GORICA log-likelihood
#     (which is asymptotically equivalent for differences; as in benchmark()),
#     so then only the CIs for differences, weights, and ratios are
#     bootstrapped; those for the single GORIC(C) values are the analytic
#     casewise ones.
#
# (D) Inversion (method = "invert"; cf. the noncentral chi-square based CIs
#     for (differences in) expected BIC in Preacher & Merkle, 2012, Eqs. 8-9):
#     for a hypothesis H_m versus (i) the unconstrained hypothesis or (ii) its
#     own complement, the difference T = IC_m - IC_k depends on the true
#     parameters (mainly) via their signed distance s to the boundary of H_m
#     (in the V^{-1} metric; s > 0 outside H_m, s < 0 inside). Let b^ be the
#     closest boundary point of H_m to theta^ (the projection of theta^ onto
#     H_m if theta^ is outside H_m, otherwise the closest point on its
#     boundary, i.e., the complement projection) and u the outward unit
#     normal (u = +/- (theta^ - b^) / ||theta^ - b^||). On the line
#       theta(s) = b^ + s * u   (so theta(s_obs) = theta^),
#     the distribution F_s of T* (with theta* ~ N(theta(s), V)) is
#     stochastically increasing in s. The confidence set for s is
#       {s : (1 - level)/2 <= F_s(T_obs) <= (1 + level)/2}
#     (with mid-p for the point mass of T, e.g., at c * (PT_m - PT_u) when
#     theta* is in H_m), and the CI for the expected difference E_s[T] is
#     [E_{s_L}[T], E_{s_U}[T]] (and, via exp(log(p_m/p_k) - ./2), for the
#     corresponding ratio of weights).
#     This one-dimensional reduction is valid when H_m is (a) equality-only 
#     (then the distribution depends on the true parameters only via the 
#     noncentrality s^2) or (b) a single inequality constraint (a half-space;
#     then only via the signed distance s). For a cone with two or more 
#     inequality constraints (e.g., A < B < C), the distribution also depends
#     on where along the boundary the parameters are: the line through b^ 
#     leaves the cone through another face, F_s is then not monotone in s, 
#     and the inversion fails (in a simulation, coverage dropped to .01-.22).
#     Therefore, the inversion is only used for (a) and (b); all other pairs
#     keep the bootstrap CI.
#     For an equality-only H_m versus the
#     unconstrained hypothesis, T - c * (PT_m - PT_u) ~ noncentral
#     chi-square(q, s^2) exactly, so this is exactly the noncentrality-
#     parameter inversion of Preacher and Merkle (2012); for inequality
#     constraints it is its chi-bar-square analogue. The distributions are
#     obtained by simulation (B draws, common random numbers over s). Unlike
#     the percentile bootstrap (which uses the distribution at s_obs itself,
#     shifted upward by about the number of active constraints), the
#     inversion targets the expected value directly (note: the sample 
#     difference can lie outside this CI, since, e.g., E[D_m] >= the number 
#     of equality constraints, while D_m itself can be close to 0). As for the bootstrap, a GORIC(C) object is evaluated with
#     the GORICA-type log-likelihood.
#
# (E) Weights and ratios of weights: these are functions of eta, so for each
#     method their CIs follow from (the distribution of) the IC values:
#     - ratio w_m / w_k = exp(eta_m - eta_k): monotone in IC_m - IC_k, so its CI
#       is exp(log(p_m/p_k) - ./2) of the CI for that difference (all methods);
#     - weight w_m: depends on all differences jointly; "delta": delta method
#       on the logit scale (see 4. above); "mc": Monte Carlo from N(eta^, Sigma);
#       "boot"/"invert": percentiles of softmax(eta*) over the bootstrap draws.
#       With only two hypotheses, w_1 = plogis(log(p_1/p_2) - (IC_1 - IC_2)/2)
#       is monotone in the difference, so its CI is the transformation of the
#       CI for the difference; for "invert" (eligible pair) the inversion CI.
#     The population counterparts are the weights/ratios based on the
#     expected IC values.
#
# Choosing a method (based on the simulations below):
#   single GORIC(C) value ............ analytic casewise ("delta")
#   single GORICA(C) value ........... "boot" (analytic CI degenerate)
#   difference / ratio, H_m with only equalities or one inequality, versus
#     the unconstrained hypothesis or its complement ... "invert"
#   difference / ratio, other pairs .. "boot"
#   weights, two hypotheses .......... as for their difference ("invert" if
#                                      eligible, otherwise "boot")
#   weights, more hypotheses ......... "boot"
#   => in general: method = "invert" (inversion where valid, bootstrap
#      otherwise); "delta"/"mc" only for clearly distinguishable hypotheses
#      in large samples (fast, no simulation).
#
# -----------------------------------------------------------------------------
# Important caveat: overlapping (indistinguishable) hypotheses
# -----------------------------------------------------------------------------
# The normal approximation requires that the two hypotheses are
# distinguishable (omega^2 > 0; Vuong, 1989; Merkle et al., 2016). When two
# hypotheses have the same maximum order-restricted log likelihood (i.e., loglik in
# goric(); because they lead to the same order-restricted MLEs) -- e.g., H1 versus the
# unconstrained hypothesis when the estimates are in agreement with H1, or
# overlapping hypotheses -- the casewise differences are (asymptotically) zero,
# SE_mk = 0, and the difference in IC values has a non-normal (chi-bar-square
# type) distribution with a spread of order 1. A normal-based CI would then
# wrongly collapse to a point. Such pairs are therefore flagged
# ('same_loglik = TRUE') and get no analytic CI (NA). The analytic CIs of
# the weights ignore this extra (order 1) variability as well, so these are
# too narrow when the set contains overlapping hypotheses (a note is
# given). The bootstrap (method = "boot") does capture this variability and
# gives CIs for all pairs.
#
# Accuracy (small simulation, 95% CIs, three groups, H1: a < b < c vs.
# H2: a < b, true means (0, .5, .2), so H2 true and H1 false):
# - analytic CIs: coverage of the CI for the difference in GORIC(A) values
#   about .83-.86 for n = 60-150 (when the hypotheses were distinguishable in
#   the sample; in 6-17% of the samples they were not), .87 for n = 600, and
#   .92 for n = 2000 (similar for means (0, 1, .4): .89 for n = 150, .92 for
#   n = 600). So, asymptotically correct, but too narrow when the hypotheses
#   are close relative to the sample size (i.e., when the true parameters are
#   near the boundaries of the hypotheses), because the difference in IC
#   values is then skewed (noncentral chi-bar-square like).
# - parametric-bootstrap percentile CIs (method = "boot", all samples, B = 500):
#   GORICA .93 (n = 60 and n = 150); GORIC .91 (n = 60) and .92 (n = 150);
#   in a first prototype also .95 for n = 600 and for means (0, 1, .4),
#   n = 150. The 'basic' and bias-corrected (BC) bootstrap CIs performed
#   worse (.51-.79 and .82-.95, resp.), so the percentile CI is used.
#
# - inversion CIs (method = "invert"; 95%, B = 500, n = 60/150), for
#   H2: a < b (single constraint) versus the unconstrained hypothesis:
#   coverage .98-1.00 in all settings (true means (0, .5, .2), (0, .2, .4), 
#   (0, 0, 0)), where the bootstrap gave .58, .81, .91, and 1.00 (it 
#   undercovers here because its distribution is degenerate when theta^ is 
#   far inside the hypothesis). For an equality-only hypothesis, it 
#   reproduced the exact noncentral chi-square inversion (e.g., [-2.02, 6.65]
#   versus [-2.00, 6.71]). For H1: a < b < c, the bootstrap (.95-.995) is used.
#
# - weights (H1, H2, unconstrained; true means (0, .5, .2), (0, .3, .6),
#   (0, 1, .4); n = 60/150): coverage of the analytic CIs .29-.89 ("delta")
#   and .29-.89 ("mc"), of the bootstrap CIs .93-.97, except for H2 when H1
#   (a subset of H2) is true (.73).
# - weights with two hypotheses (H: a < b versus the unconstrained hypothesis
#   or its complement; n = 60): coverage .05-.82 ("delta"), .61-.95 ("boot"),
#   and .96-.98 ("invert").
#
# Further notes:
# - The CIs concern the sampling variability of the IC (difference / weight /
#   ratio) as an estimate of its population counterpart; they are not
#   posterior intervals.
# - The penalty terms (and level probabilities) are treated as known.
#
# References
# Merkle, E. C., You, D., & Preacher, K. J. (2016). Testing nonnested
#   structural equation models. Psychological Methods, 21(2), 151-163.
#   doi:10.1037/met0000038
# Preacher, K. J., & Merkle, E. C. (2012). The problem of model selection
#   uncertainty in structural equation modeling. Psychological Methods, 17(1),
#   1-14. doi:10.1037/a0026804
# Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#   evaluation using leave-one-out cross-validation and WAIC. Statistics and
#   Computing, 27(5), 1413-1432. doi:10.1007/s11222-016-9696-4
# Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
#   non-nested hypotheses. Econometrica, 57(2), 307-333.
# =============================================================================


goric_CI <- function(object, level = 0.95, method = c("invert", "boot", "delta", "mc"),
                     approach = c("auto", "casewise", "estimates"),
                     nsim = 10000L, B = 1000L, seed = NULL, tol = 1e-8) {

  if (!inherits(object, "con_goric")) {
    stop("\nrestriktor ERROR: The object should be of class 'con_goric' ",
         "(i.e., output from the goric() function).", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("\nrestriktor ERROR: The argument 'level' should be a single number in between 0 and 1.",
         call. = FALSE)
  }
  method   <- match.arg(method)
  approach <- match.arg(approach)
  if (!is.null(seed)) set.seed(seed)

  type   <- object$type
  res    <- object$result
  models <- as.character(res$model)
  M      <- length(models)
  IC     <- res[[type]]
  names(IC) <- models
  pf     <- if (is.null(object$penalty_factor)) 2 else object$penalty_factor

  if (approach == "auto") {
    approach <- if (type %in% c("goric", "goricc")) "casewise" else "estimates"
  }
  if (approach == "casewise" && !(type %in% c("goric", "goricc"))) {
    stop("\nrestriktor ERROR: approach = 'casewise' is only available for type = 'goric' or 'goricc'.",
         call. = FALSE)
  }

  # prior model probabilities (as used in the goric object)
  prior <- object$priorICweights
  if (is.null(prior) || length(prior) != M) {
    prior <- rep(1 / M, M)
  }
  prior <- prior / sum(prior)

  # restricted estimates for each hypothesis in the set (rows), for the
  # (unrestricted) parameters (columns)
  theta_hat <- object$b.unrestr
  B_restr <- as.matrix(object$ormle$b.restr)
  if (!all(models %in% rownames(B_restr))) {
    stop("\nrestriktor ERROR: The restricted estimates of (some of) the hypotheses could not be found.",
         call. = FALSE)
  }
  B_restr <- B_restr[models, , drop = FALSE]
  if (!is.null(names(theta_hat)) && all(names(theta_hat) %in% colnames(B_restr))) {
    B_restr <- B_restr[, names(theta_hat), drop = FALSE] # drop defined (:=) parameters
  } else {
    B_restr <- B_restr[, seq_along(theta_hat), drop = FALSE]
  }

  # Sigma = Cov(loglik_1, ..., loglik_M) (= Cov(eta))
  if (approach == "casewise") {
    L <- goric_casewise_loglik(object, B_restr, models)   # n x M
    # check: the casewise contributions should add up to the loglik values
    # (up to a constant that is the same for all hypotheses, which does not
    # affect any of the results)
    dif <- colSums(L) - res$loglik
    if (any(!is.finite(dif)) ||
        max(abs(dif - mean(dif))) > 1e-6 * max(1, abs(res$loglik))) {
      stop("\nrestriktor ERROR: The casewise log-likelihood contributions do not add up to ",
           "the log-likelihood values (", paste(sQuote(models), collapse = ", "), "). ",
           "Use approach = 'estimates' instead.", call. = FALSE)
    }
    Lc <- sweep(L, 2, colMeans(L))
    Sigma <- crossprod(Lc)                  # = n * S (S with divisor n)
    n_obs <- nrow(L)
  } else {
    V <- as.matrix(object$VCOV)
    if (nrow(V) != length(theta_hat)) {
      stop("\nrestriktor ERROR: The dimensions of VCOV and the estimates do not match.", call. = FALSE)
    }
    Vinv <- solve(V)
    R <- sweep(-B_restr, 2, theta_hat, "+")       # rows: r_m = theta^ - theta_m~
    Sigma <- R %*% Vinv %*% t(R)
    if (type %in% c("gorica", "goricac")) {
      # check: D_m = diag(Sigma) and -2 loglik_m = log|2 pi V| + D_m
      ll_check <- -0.5 * (determinant(2 * pi * V, logarithm = TRUE)$modulus + diag(Sigma))
      check <- abs(ll_check - res$loglik) > 1e-6 * pmax(1, abs(res$loglik))
      if (any(check)) {
        warning("\nrestriktor WARNING: The GORICA log-likelihood values could not be reproduced ",
                "from the estimates for ", paste(sQuote(models[check]), collapse = ", "), ".",
                call. = FALSE)
      }
    }
    n_obs <- NA
  }
  dimnames(Sigma) <- list(models, models)

  z <- qnorm(1 - (1 - level) / 2)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  eta <- log(prior) - IC / 2
  w   <- exp(eta - max(eta)); w <- w / sum(w)

  # parametric bootstrap draws (method = "boot", and also used for "invert")
  use_boot <- method %in% c("boot", "invert")
  if (use_boot) {
    IC_boot <- goric_CI_boot(object, models, B = B, penalty_factor = pf)  # B x M
    # check: the bootstrap ICs at theta^ itself reproduce the (GORICA) IC
    # differences of the object
    IC_obs <- goric_CI_boot(object, models, B = 0L, penalty_factor = pf)
    if (type %in% c("gorica", "goricac")) {
      chk <- IC_obs - IC
      if (max(abs(chk)) > 1e-6 * max(1, abs(IC))) {
        warning("\nrestriktor WARNING: The bootstrap could not reproduce the ", toupper(type),
                " values of the object; please check the results.", call. = FALSE)
      }
    }
    eta_boot <- sweep(-IC_boot / 2, 2, log(prior), "+")
  }

  # 1. IC values
  se_IC <- 2 * sqrt(pmax(diag(Sigma), 0))
  if (use_boot && type %in% c("gorica", "goricac")) {
    q_IC <- t(apply(IC_boot, 2, quantile, probs = probs, names = FALSE))
  } else {
    q_IC <- cbind(IC - z * se_IC, IC + z * se_IC)
  }
  IC_tab <- data.frame(model = models, IC = IC, SE = se_IC,
                       lower = q_IC[, 1], upper = q_IC[, 2],
                       row.names = NULL, stringsAsFactors = FALSE)
  names(IC_tab)[2] <- type

  # 2. differences in IC values and 3. ratios of weights (all ordered pairs)
  pairs <- expand.grid(m = seq_len(M), k = seq_len(M))
  pairs <- pairs[pairs$m != pairs$k, , drop = FALSE]
  pairs <- pairs[order(pairs$m, pairs$k), , drop = FALSE]
  var_d <- Sigma[cbind(pairs$m, pairs$m)] + Sigma[cbind(pairs$k, pairs$k)] -
    2 * Sigma[cbind(pairs$m, pairs$k)]
  se_d  <- 2 * sqrt(pmax(var_d, 0))
  dIC   <- IC[pairs$m] - IC[pairs$k]
  log_ratio <- eta[pairs$m] - eta[pairs$k]
  # same_loglik: the two hypotheses have the same maximum order-restricted log
  # likelihood (i.e., loglik in goric()) in this sample (since they lead to the 
  # same order-restricted MLEs), so that the SE of their difference is (numerically) zero
  # relative to the IC scale ('indistinguishable' in Vuong's terms).
  scale_ref <- max(1, sqrt(max(abs(diag(Sigma)))))
  same_loglik <- !(se_d > tol * scale_ref)

  ci_method <- rep(if (use_boot) "boot" else "analytic", nrow(pairs))
  if (use_boot) {
    q_d <- t(vapply(seq_len(nrow(pairs)), function(j) {
      quantile(IC_boot[, pairs$m[j]] - IC_boot[, pairs$k[j]], probs = probs, names = FALSE)
    }, numeric(2)))
    q_r <- t(vapply(seq_len(nrow(pairs)), function(j) {
      exp(quantile(eta_boot[, pairs$m[j]] - eta_boot[, pairs$k[j]], probs = probs, names = FALSE))
    }, numeric(2)))
  } else {
    q_d <- cbind(ifelse(!same_loglik, dIC - z * se_d, NA_real_),
                 ifelse(!same_loglik, dIC + z * se_d, NA_real_))
    q_r <- cbind(ifelse(!same_loglik, exp(log_ratio - z * se_d / 2), NA_real_),
                 ifelse(!same_loglik, exp(log_ratio + z * se_d / 2), NA_real_))
  }

  # method = "invert": for a hypothesis versus the unconstrained hypothesis or
  # versus its own complement, the CI for the (expected) difference is
  # obtained by inverting its sampling distribution (see (D) above); the
  # other pairs keep the bootstrap CIs.
  if (method == "invert") {
    inv <- goric_CI_invert(object, models, pairs, level = level, B = B,
                           penalty_factor = pf)
    ok <- !is.na(inv[, 1])
    if (any(ok)) {
      q_d[ok, ] <- inv[ok, ]
      lp <- log(prior[pairs$m]) - log(prior[pairs$k])
      q_r[ok, ] <- cbind(exp(lp[ok] - inv[ok, 2] / 2), exp(lp[ok] - inv[ok, 1] / 2))
      ci_method[ok] <- "invert"
    }
  }

  diff_tab <- data.frame(model = models[pairs$m], vs = models[pairs$k],
                         difference = dIC, SE = se_d,
                         lower = q_d[, 1], upper = q_d[, 2],
                         same_loglik = same_loglik,
                         CI_method = ci_method,
                         row.names = NULL, stringsAsFactors = FALSE)
  ratio_tab <- data.frame(model = models[pairs$m], vs = models[pairs$k],
                          ratio = exp(log_ratio),
                          lower = q_r[, 1], upper = q_r[, 2],
                          same_loglik = same_loglik,
                          CI_method = ci_method,
                          row.names = NULL, stringsAsFactors = FALSE)

  # 4. weights
  if (method == "delta") {
    ww <- matrix(NA_real_, M, 2)
    for (m in seq_len(M)) {
      if (M == 1L) break
      g <- -w / (1 - w[m])
      g[m] <- 1
      if (!is.finite(sum(g))) { # w_m (numerically) 1
        g <- rep(0, M); g[m] <- 1
      }
      se_logit <- sqrt(max(drop(t(g) %*% Sigma %*% g), 0))
      wm <- min(max(w[m], 1e-300), 1 - 1e-16)
      ww[m, ] <- plogis(qlogis(wm) + c(-1, 1) * z * se_logit)
    }
  } else {
    if (method == "mc") {
      ev <- eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE)
      A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), M)
      Z <- matrix(rnorm(nsim * M), nsim, M)
      eta_sim <- sweep(Z %*% t(A), 2, eta, "+")
    } else {
      eta_sim <- eta_boot
    }
    eta_sim <- eta_sim - apply(eta_sim, 1, max)
    w_sim <- exp(eta_sim); w_sim <- w_sim / rowSums(w_sim)
    ww <- t(apply(w_sim, 2, quantile, probs = probs, names = FALSE))
  }
  w_method <- rep(switch(method, delta = "delta", mc = "mc", boot = "boot", invert = "boot"), M)
  # With two hypotheses, the weight is a monotone function of the difference
  # in IC values: w_1 = plogis(log(p_1/p_2) - (IC_1 - IC_2)/2) and w_2 = 1 - w_1.
  # Then, the CI for the weights follows directly from the CI for the
  # difference; for method = "invert", this gives the inversion CI for the
  # weights (when the pair is eligible; see (D)).
  if (method == "invert" && M == 2L) {
    j12 <- which(pairs$m == 1 & pairs$k == 2)
    if (identical(ci_method[j12], "invert")) {
      lp12 <- log(prior[1]) - log(prior[2])
      w1_ci <- plogis(lp12 - rev(q_d[j12, ]) / 2)
      ww <- rbind(w1_ci, 1 - rev(w1_ci))
      w_method <- rep("invert", 2)
    }
  }
  weight_tab <- data.frame(model = models, weight = w,
                           lower = ww[, 1], upper = ww[, 2],
                           CI_method = w_method,
                           row.names = NULL, stringsAsFactors = FALSE)
  names(weight_tab)[2] <- paste0(type, ".weights")

  out <- list(type = type, level = level, method = method, approach = approach,
              n = n_obs, B = if (use_boot) B else NULL,
              penalty_factor = pf,
              priorICweights = prior,
              Sigma = Sigma,
              IC = IC_tab, difference = diff_tab,
              weights = weight_tab, ratio_weights = ratio_tab,
              best = models[which.max(w)],
              any_same_loglik = any(same_loglik))
  class(out) <- "goric_CI"
  out
}


# -----------------------------------------------------------------------------
# Constraint information (Amat, rhs, neq) per row of the result table, in the
# parameter space of object$b.unrestr. Returns a list with, per model,
# list(kind = "hypothesis" / "unconstrained" / "complement", Amat, rhs, neq).
goric_CI_row_specs <- function(object, models) {
  hyp_names <- names(object$constraints)
  get_h <- function(h) {
    Amat <- as.matrix(object$constraints[[h]])
    rhs  <- object$rhs[[h]]
    if (is.null(rhs)) rhs <- rep(0, nrow(Amat))
    neq  <- object$neq[[h]]
    if (is.null(neq)) neq <- 0L
    list(Amat = Amat, rhs = rhs, neq = neq)
  }
  lapply(models, function(m) {
    if (m == "unconstrained") {
      list(kind = "unconstrained")
    } else if (m %in% hyp_names) {
      c(list(kind = "hypothesis"), get_h(m))
    } else if (m == "complement") {
      ref <- setdiff(hyp_names, "Heq")[1]
      c(list(kind = "complement"), get_h(ref))
    } else if (!is.null(object$add_Hc) && m == paste0("Complement of ", object$add_Hc)) {
      c(list(kind = "complement"), get_h(object$add_Hc))
    } else {
      stop("\nrestriktor ERROR: Could not determine the constraints of hypothesis ",
           sQuote(m), ".", call. = FALSE)
    }
  })
}

# Pieces (Amat, rhs, neq) of the boundary of a hypothesis, used for its
# complement: as in compute_complement_likelihood() (goric_compute_loglik_complement.R),
# each inequality constraint in turn is set to an equality, with the other
# inequality constraints kept; returns NULL if the complement equals the
# unconstrained hypothesis (only equality constraints).
goric_CI_complement_pieces <- function(spec) {
  Amat <- spec$Amat; rhs <- spec$rhs; meq <- spec$neq
  keep <- apply(Amat, 1, function(x) !all(x == 0))
  Amat <- Amat[keep, , drop = FALSE]; rhs <- rhs[keep]
  if (nrow(Amat) == meq) return(NULL)
  nr <- seq_len(nrow(Amat))
  if (meq > 0L) nr <- nr[-seq_len(meq)]
  lapply(seq_along(nr), function(l) {
    idx <- c(nr[l], nr[-l])
    list(Amat = Amat[idx, , drop = FALSE], rhs = rhs[idx], neq = 1L,
         Amat_full = Amat, rhs_full = rhs, meq_full = meq)
  })
}

# Squared Mahalanobis distance from 'th' to {theta: Amat theta (=/>=) rhs}
goric_CI_proj_D <- function(th, Vinv, Amat, rhs, neq) {
  sol <- quadprog::solve.QP(Dmat = Vinv, dvec = drop(Vinv %*% th),
                            Amat = t(Amat), bvec = rhs, meq = neq)
  d <- th - sol$solution
  drop(t(d) %*% Vinv %*% d)
}

# Parametric bootstrap: B x M matrix of IC* values (GORICA-type log-likelihood
# plus the fixed penalties of the object). With B = 0, the IC values at
# theta^ itself are returned (a vector), used as a check.
goric_CI_boot <- function(object, models, B, penalty_factor = 2) {
  theta_hat <- c(object$b.unrestr)
  V <- as.matrix(object$VCOV)
  if (nrow(V) != length(theta_hat)) {
    stop("\nrestriktor ERROR: The dimensions of VCOV and the estimates do not match.", call. = FALSE)
  }
  Vinv <- solve(V)
  specs <- goric_CI_row_specs(object, models)
  for (j in seq_along(specs)) {
    if (!is.null(specs[[j]]$Amat) && ncol(specs[[j]]$Amat) != length(theta_hat)) {
      stop("\nrestriktor ERROR: The constraint matrix of ", sQuote(models[j]),
           " does not match the number of parameters.", call. = FALSE)
    }
    if (specs[[j]]$kind == "complement") {
      specs[[j]]$pieces <- goric_CI_complement_pieces(specs[[j]])
    }
  }
  PT <- object$result$penalty
  const <- as.numeric(determinant(2 * pi * V, logarithm = TRUE)$modulus)
  tol <- sqrt(.Machine$double.eps)

  D_one <- function(th) {
    vapply(specs, function(sp) {
      if (sp$kind == "unconstrained") return(0)
      if (sp$kind == "hypothesis") {
        return(goric_CI_proj_D(th, Vinv, sp$Amat, sp$rhs, sp$neq))
      }
      # complement
      if (is.null(sp$pieces)) return(0)
      A <- sp$pieces[[1]]$Amat_full; b <- sp$pieces[[1]]$rhs_full; meq <- sp$pieces[[1]]$meq_full
      Ath <- drop(A %*% th) - b
      viol_eq  <- meq > 0L && any(abs(Ath[seq_len(meq)]) > tol)
      viol_ineq <- any(Ath[setdiff(seq_along(Ath), seq_len(meq))] < -tol)
      if (viol_eq || viol_ineq) return(0)   # th itself is in the complement
      min(vapply(sp$pieces, function(pc)
        goric_CI_proj_D(th, Vinv, pc$Amat, pc$rhs, pc$neq), numeric(1)))
    }, numeric(1))
  }

  if (B == 0L) {
    return(const + D_one(theta_hat) + penalty_factor * PT)
  }
  draws <- mvtnorm::rmvnorm(B, mean = theta_hat, sigma = V)
  D <- t(apply(draws, 1, D_one))
  if (length(models) == 1L) D <- t(D)
  IC_b <- sweep(const + D, 2, penalty_factor * PT, "+")
  colnames(IC_b) <- models
  IC_b
}


# -----------------------------------------------------------------------------
# Inversion CIs (method = "invert") for the pairs (H_m, unconstrained) and
# (H_m, complement of H_m), and their reverses; NA for the other pairs.
# Returns a (number of pairs) x 2 matrix with the lower and upper bounds of
# the CI for the expected difference in IC values.
goric_CI_invert <- function(object, models, pairs, level, B, penalty_factor = 2) {
  theta_hat <- c(object$b.unrestr)
  V <- as.matrix(object$VCOV)
  Vinv <- solve(V)
  p <- length(theta_hat)
  specs <- goric_CI_row_specs(object, models)
  PT <- object$result$penalty
  alpha <- 1 - level
  tol <- sqrt(.Machine$double.eps)
  out <- matrix(NA_real_, nrow(pairs), 2)

  same_constraints <- function(a, b) {
    isTRUE(all.equal(a$Amat, b$Amat, check.attributes = FALSE)) &&
      isTRUE(all.equal(a$rhs, b$rhs, check.attributes = FALSE)) &&
      identical(as.integer(a$neq), as.integer(b$neq))
  }
  # one-dimensional family available: equality-only hypothesis, or a single
  # inequality constraint (half-space); see (D) in the header for why not for
  # cones with two or more inequality constraints
  one_dim <- function(sp) {
    nz <- apply(sp$Amat, 1, function(x) !all(x == 0))
    nr <- sum(nz)
    nr > 0 && (nr == sp$neq || (sp$neq == 0 && nr == 1))
  }
  # eligible: m such a hypothesis, k unconstrained or the complement of m
  eligible <- function(m, k) {
    sm <- specs[[m]]; sk <- specs[[k]]
    sm$kind == "hypothesis" && one_dim(sm) &&
      (sk$kind == "unconstrained" || (sk$kind == "complement" && same_constraints(sm, sk)))
  }

  proj <- function(th, Amat, rhs, neq) {
    sol <- quadprog::solve.QP(Dmat = Vinv, dvec = drop(Vinv %*% th),
                              Amat = t(Amat), bvec = rhs, meq = neq)$solution
    d <- th - sol
    list(D = drop(t(d) %*% Vinv %*% d), sol = sol)
  }
  # D for hypothesis m and for its complement (over the rows of TH)
  D_hyp <- function(TH, sp) apply(TH, 1, function(th) proj(th, sp$Amat, sp$rhs, sp$neq)$D)
  D_compl <- function(TH, sp) {
    pieces <- goric_CI_complement_pieces(sp)
    if (is.null(pieces)) return(rep(0, nrow(TH)))
    A <- pieces[[1]]$Amat_full; b <- pieces[[1]]$rhs_full; meq <- pieces[[1]]$meq_full
    apply(TH, 1, function(th) {
      Ath <- drop(A %*% th) - b
      viol <- (meq > 0L && any(abs(Ath[seq_len(meq)]) > tol)) ||
        any(Ath[setdiff(seq_along(Ath), seq_len(meq))] < -tol)
      if (viol) return(0)
      min(vapply(pieces, function(pc) proj(th, pc$Amat, pc$rhs, pc$neq)$D, numeric(1)))
    })
  }

  # common random numbers (the same for all pairs and all s)
  E <- mvtnorm::rmvnorm(B, mean = rep(0, p), sigma = V)

  done <- list()
  for (j in seq_len(nrow(pairs))) {
    m <- pairs$m[j]; k <- pairs$k[j]; sgn <- 1
    if (!eligible(m, k)) {
      if (eligible(k, m)) { tmp <- m; m <- k; k <- tmp; sgn <- -1 } else next
    }
    key <- paste(m, k)
    if (is.null(done[[key]])) {
      sm <- specs[[m]]; sk <- specs[[k]]
      cPT <- penalty_factor * (PT[m] - PT[k])
      Tstat <- function(TH) {
        Dm <- D_hyp(TH, sm)
        Dk <- if (sk$kind == "unconstrained") 0 else D_compl(TH, sk)
        Dm - Dk + cPT
      }
      T_obs <- Tstat(matrix(theta_hat, nrow = 1))

      # closest boundary point b^ of H_m and outward unit normal u
      pr <- proj(theta_hat, sm$Amat, sm$rhs, sm$neq)
      if (pr$D > tol) {                       # theta^ outside H_m
        b_hat <- pr$sol; dist <- sqrt(pr$D); u <- (theta_hat - b_hat) / dist; s_obs <- dist
      } else {                                # theta^ in H_m
        pieces <- goric_CI_complement_pieces(sm)
        bd <- if (is.null(pieces)) NULL else
          lapply(pieces, function(pc) proj(theta_hat, pc$Amat, pc$rhs, pc$neq))
        dist <- if (is.null(bd)) 0 else min(vapply(bd, function(x) x$D, numeric(1)))
        if (dist > tol) {
          b_hat <- bd[[which.min(vapply(bd, function(x) x$D, numeric(1)))]]$sol
          dist <- sqrt(dist); u <- (b_hat - theta_hat) / dist; s_obs <- -dist
        } else {                              # (numerically) on the boundary
          A <- sm$Amat; res <- abs(drop(A %*% theta_hat) - sm$rhs)
          a <- A[which.min(res), ]
          u <- -drop(V %*% a) / sqrt(drop(t(a) %*% V %*% a))
          b_hat <- theta_hat; s_obs <- 0
        }
      }

      T_at <- function(s) Tstat(sweep(E, 2, b_hat + s * u, "+"))
      F_at <- function(s) {                   # mid-p cdf of T at T_obs
        Ts <- T_at(s)
        mean(Ts < T_obs - 1e-10) + 0.5 * mean(abs(Ts - T_obs) <= 1e-10)
      }
      W <- 8 + abs(s_obs)
      lo_end <- s_obs - W; hi_end <- s_obs + W
      # an equality-only hypothesis has no 'inside': the distribution depends
      # on s only via s^2 (noncentrality), so s >= 0
      if (nrow(sm$Amat) == sm$neq) lo_end <- max(lo_end, 0)
      solve_s <- function(target) {
        f_lo <- F_at(lo_end) - target; f_hi <- F_at(hi_end) - target
        if (f_lo < 0) return(lo_end)          # target never reached: limit
        if (f_hi > 0) return(hi_end)
        uniroot(function(s) F_at(s) - target, c(lo_end, hi_end), tol = 1e-3)$root
      }
      s_L <- solve_s(1 - alpha / 2)
      s_U <- solve_s(alpha / 2)
      # For an equality-only hypothesis, s = 0 is the boundary of the parameter 
      # space of s; when a limit equals 0, the expected difference is known exactly: 
      # T - c * (PT_m - PT_k) is then central chi-square with q = rank(Amat) 
      # degrees of freedom, so E_0[T] = q + c * (PT_m - PT_k). This exact value is 
      # used instead of its Monte Carlo estimate (which would otherwise lie above 
      # or below the true value about half of the time when the hypothesis is true).
      exact0 <- if (nrow(sm$Amat) == sm$neq) qr(sm$Amat)$rank + cPT else NA_real_
      E_T <- function(s) if (!is.na(exact0) && s <= 0) exact0 else mean(T_at(s))
      done[[key]] <- c(E_T(s_L), E_T(s_U))
    }
    ci <- done[[key]]
    out[j, ] <- if (sgn == 1) ci else rev(-ci)
  }
  out
}


# -----------------------------------------------------------------------------
# Casewise log-likelihood contributions l_im (n x M matrix), evaluated at the
# restricted estimates in the rows of B (in the order of the result table),
# or -- for missing = "fiml" -- at the implied mean and covariance matrix of
# the (restricted) FIML solutions.
goric_casewise_loglik <- function(object, B, models = rownames(B)) {
  fit <- object$model.org
  if (is.null(fit) || !inherits(fit, "lm")) {
    stop("\nrestriktor ERROR: The casewise approach needs the fitted model (lm, mlm, rlm, or glm) ",
         "in the goric object. Use approach = 'estimates' instead.", call. = FALSE)
  }
  if (!is.null(object$objectList[[1]]$fiml)) {
    return(goric_casewise_loglik_fiml(object, models))
  }

  X <- model.matrix(fit)
  mf <- model.frame(fit)
  # for a glm, fit$y is the response as used by the family (e.g., proportions
  # for a binomial with cbind(successes, failures); then the prior weights
  # are the numbers of trials)
  y <- if (inherits(fit, "glm")) fit$y else model.response(mf)
  off <- model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(X))
  pw <- if (inherits(fit, "glm")) fit$prior.weights else fit$weights
  if (is.null(pw)) pw <- rep(1, nrow(X))
  keep <- pw > 0

  # mlm: multivariate normal with the ML residual covariance matrix
  if (inherits(fit, "mlm")) {
    Y <- as.matrix(y); q <- ncol(Y); p <- ncol(X)
    if (ncol(B) != p * q) {
      stop("\nrestriktor ERROR: The number of parameters does not match the mlm coefficients.", call. = FALSE)
    }
    L <- apply(B, 1, function(b) {
      Bm <- matrix(b, nrow = p, ncol = q)    # vec(coef) (column-major, per response)
      Rs <- Y - X %*% Bm - off
      S  <- crossprod(Rs) / nrow(Rs)
      mvtnorm::dmvnorm(Rs, sigma = S, log = TRUE)
    })
    L <- as.matrix(L); colnames(L) <- rownames(B)
    return(L)
  }

  family <- if (inherits(fit, "glm")) fit$family$family else "gaussian"
  linkinv <- if (inherits(fit, "glm")) fit$family$linkinv else identity
  dev_resids <- if (inherits(fit, "glm")) fit$family$dev.resids else NULL

  if (!all(colnames(B) %in% colnames(X))) {
    stop("\nrestriktor ERROR: The parameter names do not match the model matrix.", call. = FALSE)
  }
  X <- X[, colnames(B), drop = FALSE]

  L <- apply(B, 1, function(b) {
    mu <- linkinv(drop(X %*% b) + off)
    if (family == "gaussian") {
      # lm, rlm, gaussian glm: ML variance (as in con_loglik_lm() and the glm AIC)
      r  <- y - mu
      s2 <- sum(pw[keep] * r[keep]^2) / sum(keep)
      0.5 * log(pw[keep]) - 0.5 * log(2 * pi * s2) - pw[keep] * r[keep]^2 / (2 * s2)
    } else if (family == "poisson") {
      dpois(y[keep], mu[keep], log = TRUE) * pw[keep]
    } else if (family == "binomial") {
      m <- pw[keep]
      dbinom(round(m * y[keep]), round(m), mu[keep], log = TRUE)
    } else if (family == "Gamma") {
      # dispersion as in Gamma()$aic: deviance / sum of the prior weights
      disp <- sum(dev_resids(y, mu, pw)) / sum(pw)
      dgamma(y[keep], shape = 1 / disp, scale = mu[keep] * disp, log = TRUE) * pw[keep]
    } else if (family == "inverse.gaussian") {
      # dispersion as in inverse.gaussian()$aic: deviance / sum of the prior weights
      disp <- sum(dev_resids(y, mu, pw)) / sum(pw)
      yk <- y[keep]; mk <- mu[keep]; wk <- pw[keep]
      wk * (-0.5 * log(2 * pi * disp * yk^3) - (yk - mk)^2 / (2 * disp * mk^2 * yk))
    } else {
      stop("\nrestriktor ERROR: The casewise approach is not available for the ",
           sQuote(family), " family. Use approach = 'estimates' instead.", call. = FALSE)
    }
  })
  L <- as.matrix(L)
  colnames(L) <- rownames(B)
  L
}

# FIML: casewise observed-data log-likelihood of each case, under the implied
# mean and covariance matrix of each (restricted) FIML solution.
goric_casewise_loglik_fiml <- function(object, models) {
  fit  <- object$model.org
  aux  <- object$objectList[[1]]$auxiliary
  ctrl <- object$objectList[[1]]$control$fiml
  if (is.null(aux)) aux <- c()
  if (is.null(ctrl)) ctrl <- list()
  fiml <- con_fiml_lm(fit, auxiliary = aux, control = ctrl)

  Z <- fiml$Z
  obs <- !is.na(Z)
  key <- apply(obs, 1, function(o) paste(as.integer(o), collapse = ""))
  case_ll <- function(mu, Sigma) {
    out <- numeric(nrow(Z))
    for (k in unique(key)) {
      idx <- which(key == k)
      o <- obs[idx[1], ]
      if (!any(o)) next
      out[idx] <- mvtnorm::dmvnorm(Z[idx, o, drop = FALSE], mean = mu[o],
                                   sigma = Sigma[o, o, drop = FALSE], log = TRUE)
    }
    out
  }
  ecm <- function(Amat, rhs, neq) {
    fiml_ecm_restricted(fiml, Amat = Amat, bvec = rhs, meq = neq,
                        max.iter = fiml$control$max.iter, tol = fiml$control$tol)
  }

  specs <- goric_CI_row_specs(object, models)
  L <- vapply(specs, function(sp) {
    if (sp$kind == "unconstrained") {
      return(case_ll(fiml$mu, fiml$Sigma))
    }
    if (sp$kind == "hypothesis") {
      e <- ecm(sp$Amat, sp$rhs, sp$neq)
      return(case_ll(e$mu, e$Sigma))
    }
    # complement (as in compute_complement_likelihood())
    pieces <- goric_CI_complement_pieces(sp)
    if (is.null(pieces)) return(case_ll(fiml$mu, fiml$Sigma))
    A <- pieces[[1]]$Amat_full; b <- pieces[[1]]$rhs_full; meq <- pieces[[1]]$meq_full
    Ath <- drop(A %*% fiml$est) - b
    viol_eq   <- meq > 0L && !all(Ath[seq_len(meq)] == 0)
    viol_ineq <- !all(Ath[setdiff(seq_along(Ath), seq_len(meq))] >= 0)
    if (viol_eq || viol_ineq) return(case_ll(fiml$mu, fiml$Sigma))
    fits <- lapply(pieces, function(pc) ecm(pc$Amat, pc$rhs, pc$neq))
    best <- which.max(vapply(fits, function(e) e$loglik, numeric(1)))
    case_ll(fits[[best]]$mu, fits[[best]]$Sigma)
  }, numeric(nrow(Z)))
  L <- as.matrix(L)
  colnames(L) <- models
  L
}


print.goric_CI <- function(x, digits = max(3, getOption("digits") - 4),
                           all_pairs = FALSE, ...) {
  type_label <- toupper(x$type)
  pct <- paste0(format(100 * x$level), "%")
  cat(sprintf("restriktor: %s confidence intervals for %s values, differences, weights, and ratios of weights\n",
              pct, type_label))
  if (x$method == "invert") {
    cat(sprintf("(inversion CIs for the hypotheses versus the unconstrained hypothesis / their complement; %s)\n",
                paste0("other pairs (and the weights, if more than two hypotheses): parametric-bootstrap percentile CIs, B = ", x$B, "; see the CI_method columns")))
  } else if (x$method == "boot") {
    cat(sprintf("(parametric-bootstrap percentile CIs, B = %s draws from N(estimates, VCOV)%s)\n", x$B,
                ifelse(x$type %in% c("goric", "goricc"),
                       "; CIs of the single GORIC values: casewise approach", "")))
  } else {
    cat(sprintf("(approach: %s; weights: %s)\n",
                ifelse(x$approach == "casewise",
                       "casewise log-likelihood contributions (cf. Merkle, You, & Preacher, 2016)",
                       "delta method based on the estimates and their covariance matrix"),
                ifelse(x$method == "delta", "delta method on the logit scale", "Monte Carlo")))
  }

  fmt <- function(df) .format_CI_table(df, digits = digits)

  cat(sprintf("\n%s values:\n", type_label))
  print(fmt(x$IC), row.names = FALSE, right = TRUE)

  sel <- if (all_pairs) rep(TRUE, nrow(x$difference)) else x$difference$model == x$best
  cat(sprintf("\nDifferences in %s values%s:\n", type_label,
              ifelse(all_pairs, "", paste0(" (best hypothesis ", x$best, " versus the others)"))))
  print(fmt(x$difference[sel, , drop = FALSE]), row.names = FALSE, right = TRUE)

  cat(sprintf("\n%s weights:\n", type_label))
  print(fmt(x$weights), row.names = FALSE, right = TRUE)

  cat(sprintf("\nRatios of %s weights%s:\n", type_label,
              ifelse(all_pairs, "", paste0(" (best hypothesis ", x$best, " versus the others)"))))
  print(fmt(x$ratio_weights[sel, , drop = FALSE]), row.names = FALSE, right = TRUE)

  if (isTRUE(x$any_same_loglik)) {
    if (x$method %in% c("boot", "invert")) {
      message("\nrestriktor Note: Some hypotheses have the same maximum order-restricted log ",
              "likelihood (i.e., loglik in goric()) in this sample (same_loglik = TRUE). For these pairs, ",
              "only the bootstrap and inversion CIs are meaningful (the SE is 0, but the ",
              "difference in ", type_label, " values is not fixed).")
    } else {
      message("\nrestriktor Note: Some hypotheses have the same maximum order-restricted log ",
              "likelihood (i.e., loglik in goric()) (same_loglik = TRUE), for example, a hypothesis and the unconstrained hypothesis ",
              "when the estimates are in agreement with that hypothesis. For these pairs, the ",
              "difference in ", type_label, " values is not normally distributed and no CI is given. ",
              "The CIs of the weights do not reflect this (order 1) variability and are too narrow then. ",
              "Use method = \"boot\" instead.")
    }
  }
  if (!is.null(x$priorICweights) && length(unique(round(x$priorICweights, 12))) > 1) {
    message("\nrestriktor Note: Unequal prior IC weights (priorICweights) are used. These are part of ",
            "the ", type_label, " weights and the ratios of ", type_label, " weights (and their CIs), ",
            "but not of the ", type_label, " values and their differences (and their CIs). ",
            "Hence, conclusions based on the (CIs of the) ", type_label, " values may differ from ",
            "those based on the (CIs of the) ", type_label, " weights.")
  }
  if (x$type %in% c("gorica", "goricac") && !(x$method %in% c("boot", "invert"))) {
    message("\nrestriktor Note: The CI of a single ", type_label, " value is based on ",
            "2 x sqrt(squared Mahalanobis distance to the hypothesis); it is degenerate ",
            "(SE = 0) for hypotheses in agreement with the estimates. Only the relative ",
            "values (differences, weights, ratios) are of interest anyway.")
  }
  invisible(x)
}

# -----------------------------------------------------------------------------
# Text lines describing which CI method was used for the IC values, the IC
# weights, and the ratios of IC weights (used in print.con_goric() and
# summary.con_goric() when goric() computed the CIs).
goric_CI_method_text <- function(ci) {
  IC_lab <- toupper(ci$type)
  B_txt <- if (!is.null(ci$B)) paste0(", B = ", ci$B, " draws") else ""
  B_par <- if (!is.null(ci$B)) paste0(" (B = ", ci$B, " draws)") else ""
  lab <- function(m) {
    switch(m,
      analytic = "analytic (normal-based)",
      delta    = "analytic (delta method on the logit scale)",
      mc       = "analytic (Monte Carlo)",
      boot     = paste0("parametric bootstrap (percentile", B_txt, ")"),
      invert   = paste0("inversion of the sampling distribution of the difference in ",
                        IC_lab, " values", B_par),
      m)
  }
  # IC values
  if (ci$type %in% c("goric", "goricc") && ci$approach == "casewise") {
    ic_m <- "analytic (based on the casewise log-likelihood contributions)"
  } else if (ci$method %in% c("boot", "invert")) {
    ic_m <- lab("boot")
  } else {
    ic_m <- "analytic (delta method based on the estimates)"
  }
  # weights
  w_m <- unique(ci$weights$CI_method)
  w_txt <- paste(vapply(w_m, lab, character(1)), collapse = " / ")
  # ratios of weights: per unordered pair; list the pairs if methods differ
  rd <- ci$ratio_weights
  models <- unique(c(rd$model, rd$vs))
  keep <- match(rd$model, models) < match(rd$vs, models)
  rd <- rd[keep, , drop = FALSE]
  r_m <- unique(rd$CI_method)
  if (length(r_m) == 1L) {
    r_txt <- lab(r_m)
  } else {
    parts <- vapply(r_m, function(m) {
      pr <- rd[rd$CI_method == m, , drop = FALSE]
      paste0(lab(m), " for ", paste(paste0(sQuote(pr$model), " vs. ", sQuote(pr$vs)), collapse = ", "))
    }, character(1))
    r_txt <- paste(parts, collapse = "; ")
  }
  c(paste0("CI methods (see ?goric_CI for details):"),
    paste0("  - ", IC_lab, " values: ", ic_m),
    paste0("  - ", IC_lab, " weights: ", w_txt),
    paste0("  - ratios of ", IC_lab, " weights: ", r_txt))
}

# -----------------------------------------------------------------------------
# Helper: format the numeric columns of a CI table (used by print.goric_CI() and 
# print.evSyn_CI()): fixed notation, but scientific notation for (absolute) 
# values of 1000 or larger and for values smaller than 0.001 (other than 0); 
# integer columns (e.g., counts) are kept as they are; NA is shown as "-".
.format_CI_number <- function(v, digits = 3) {
  vapply(v, function(x) {
    if (is.na(x)) return("NA")
    if (!is.finite(x)) return(as.character(x))
    if (abs(x) >= 1000 || (abs(x) < 0.001 && x != 0)) {
      formatC(x, digits = max(2, digits - 1), format = "e")
    } else {
      formatC(x, digits = digits, format = "f")
    }
  }, character(1))
}
.format_CI_table <- function(df, digits = 3) {
  num <- vapply(df, function(v) is.numeric(v) && !is.integer(v), logical(1))
  df[num] <- lapply(df[num], .format_CI_number, digits = digits)
  df[] <- lapply(df, function(v) { v[v == "NA"] <- "-"; v })
  df
}
