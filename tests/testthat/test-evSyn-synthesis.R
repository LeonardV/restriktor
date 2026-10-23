# Tests for the evidence synthesis (evSyn) computations: consistency across
# input types, study weights, type_ev, leave1studyout(), variability measures
# and conclusion, confidence intervals (evSyn_CI(), CI = TRUE), and evSyn_forest().

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
set.seed(1)
S_ev <- 5
est_ev <- lapply(seq_len(S_ev), function(i)
  c(b1 = rnorm(1, .3, .1), b2 = rnorm(1, .1, .1), b3 = rnorm(1, 0, .1)))
vcov_ev <- lapply(seq_len(S_ev), function(i) diag(c(.010, .012, .011)))
H_ev <- list(H1 = "b1 > b2 > b3", H2 = "b1 > b2 = b3")
sw_ev <- c(.1, .15, .2, .25, .3)

fit_ev <- function(...) suppressMessages(evSyn(est_ev, VCOV = vcov_ev, hypotheses = H_ev, ...))
rows_to_list <- function(m) lapply(seq_len(nrow(m)), function(s) m[s, ])
final_w <- function(x) unname(x$Cumulative_GORICA_weights["Final", ])

# -----------------------------------------------------------------------------
# Consistency across input types
# -----------------------------------------------------------------------------
test_that("all input types give the same final weights (added and average, with study weights)", {
  for (tev in c("added", "average")) {
    r_est <- fit_ev(type_ev = tev, study_weights = sw_ev)
    r_ll  <- suppressMessages(evSyn(rows_to_list(r_est$LL_m), PT = rows_to_list(r_est$PT_m),
                                    type_ev = tev, study_weights = sw_ev))
    r_ic  <- suppressMessages(evSyn(rows_to_list(r_est$GORICA_m), type_ev = tev, study_weights = sw_ev))
    W     <- rows_to_list(r_est$GORICA_weight_m)
    r_w   <- suppressMessages(evSyn(W, type_ev = tev, study_weights = sw_ev))
    r_r   <- suppressMessages(evSyn(lapply(W, function(w) w / w[3]), type_ev = tev, study_weights = sw_ev))
    expect_equal(final_w(r_ll), final_w(r_est), tolerance = 1e-8)
    expect_equal(final_w(r_ic), final_w(r_est), tolerance = 1e-8)
    expect_equal(final_w(r_w),  final_w(r_est), tolerance = 1e-8)
    expect_equal(final_w(r_r),  final_w(r_est), tolerance = 1e-8)
  }
})

test_that("estimates and LL + PT input agree for all type_ev (also with re-ordered studies)", {
  for (tev in c("added", "equal", "average")) {
    r <- fit_ev(type_ev = tev, study_weights = sw_ev, order_studies = "descending")
    r_ll <- suppressMessages(evSyn(rows_to_list(r$LL_m), PT = rows_to_list(r$PT_m), type_ev = tev,
                                   study_weights = sw_ev[r$order_studies]))
    expect_equal(unname(r$Cumulative_GORICA_weights), unname(r_ll$Cumulative_GORICA_weights),
                 tolerance = 1e-8)
  }
})

# -----------------------------------------------------------------------------
# type_ev and study weights
# -----------------------------------------------------------------------------
test_that("type_ev = 'average' averages the LL: with penalty_factor = 0, GORICA weights = LL weights", {
  r <- fit_ev(type_ev = "average", penalty_factor = 0)
  expect_equal(unname(r$Cumulative_LL_weights["Final", ]), final_w(r), tolerance = 1e-8)
})

test_that("IC values equal -2 x LL + penalty_factor x penalty, also with study weights", {
  for (tev in c("added", "equal", "average")) {
    r <- fit_ev(type_ev = tev, study_weights = sw_ev)
    f <- summary(r)$Final_Cumulative_results
    expect_equal(unname(f["GORICA values", ]),
                 unname(-2 * f["Log-likelihood values", ] + 2 * f["Penalty term values", ]),
                 tolerance = 1e-8)
  }
})

test_that("study weights: only their relative sizes matter; invalid weights give an error", {
  r1 <- fit_ev(study_weights = sw_ev)
  r2 <- fit_ev(study_weights = 7 * sw_ev)
  expect_equal(final_w(r1), final_w(r2), tolerance = 1e-10)
  expect_message(evSyn(est_ev, VCOV = vcov_ev, hypotheses = H_ev, study_weights = 7 * sw_ev),
                 "do not add up to 1 or to S")
  expect_error(suppressMessages(evSyn(est_ev, VCOV = vcov_ev, hypotheses = H_ev,
                                      study_weights = c(.5, -.1, .2, .2, .2))),
               "larger than 0")
  # floating-point sums such as 0.1 + 0.2 + 0.7 do not give a message
  expect_no_message(restriktor:::.check_study_weights(c(.1, .2, .7), 3))
})

test_that("one hypothesis: by default, the complement is the failsafe hypothesis", {
  r <- suppressMessages(evSyn(lapply(est_ev, function(e) e[1:2]),
                              VCOV = lapply(vcov_ev, function(v) v[1:2, 1:2]),
                              hypotheses = list(H1 = "b1 > b2")))
  expect_equal(colnames(r$GORICA_weight_m), c("H1", "Complement"))
})

# -----------------------------------------------------------------------------
# Edge cases
# -----------------------------------------------------------------------------
test_that("one study works for all input types", {
  r <- fit_ev()
  ref <- unname(r$GORICA_weight_m[1, ])
  W1 <- r$GORICA_weight_m[1, ]
  expect_equal(final_w(suppressMessages(evSyn(est_ev[1], VCOV = vcov_ev[1], hypotheses = H_ev))),
               ref, tolerance = 1e-8)
  expect_equal(final_w(evSyn(list(W1))), ref, tolerance = 1e-8)
  expect_equal(final_w(evSyn(list(W1 / W1[3]))), ref, tolerance = 1e-8)
  expect_equal(final_w(evSyn(list(r$GORICA_m[1, ]))), ref, tolerance = 1e-8)
})

test_that("ICratios: reference hypothesis is found when the first study has two ratios of 1", {
  R2 <- list(c(1, 1, .5), c(.5, .8, 1), c(.2, .4, 1))
  expect_message(r <- evSyn(R2), "same reference hypothesis")
  expect_equal(unname(r$Href), 1L)
  r5 <- evSyn(list(c(2, .5, 1), c(3, .8, 1)), hypo_names = c("Ha", "Hb", "Hu"))
  expect_equal(names(r5$Href), "Hu")
  expect_equal(rownames(r5$GORICA_weight_m), r5$study_names)
})

test_that("tied overall preferred hypotheses in ordering give a message, not an error", {
  expect_message(r <- evSyn(list(c(.4, .4, .2), c(.3, .3, .4), c(.45, .45, .1)),
                            order_studies = "ascending"),
                 "equally preferred overall")
  expect_equal(sort(r$order_studies), 1:3)
})

test_that("escalc-type data: a single-outcome study after a multi-outcome study uses 'vi'", {
  d <- data.frame(study = c(1, 1, 2), outcome = c("a", "b", "a"), yi = c(.2, .3, .25),
                  vi = c(.01, NA, .02), v1i = c(.01, .002, NA), v2i = c(.002, .012, NA))
  x <- restriktor:::extract_est_vcov_outcomes(d)
  expect_equal(unname(x$vcov_blocks[["study 2"]][1, 1]), .02)
  expect_equal(dimnames(x$vcov_blocks[["study 1"]]), list(c("a", "b"), c("a", "b")))
})

# -----------------------------------------------------------------------------
# leave1studyout()
# -----------------------------------------------------------------------------
test_that("leave1studyout() equals re-running evSyn() without that study", {
  for (tev in c("added", "equal", "average")) {
    r <- fit_ev(type_ev = tev, study_weights = sw_ev, order_studies = "descending")
    l1 <- leave1studyout(r)
    ref <- t(sapply(seq_len(S_ev), function(k) {
      orig <- r$order_studies[-k]
      final_w(suppressMessages(evSyn(est_ev[orig], VCOV = vcov_ev[orig], hypotheses = H_ev,
                                     type_ev = tev, study_weights = sw_ev[orig])))
    }))
    expect_equal(unname(l1$OverallGoricaWeights), unname(ref), tolerance = 1e-8)
  }
})

test_that("leave1studyout() works for IC weights and IC ratios (with priors and study weights)", {
  W <- rows_to_list(fit_ev()$GORICA_weight_m)
  pr <- c(.5, .3, .2)
  for (obj in list(W, lapply(W, function(w) w / w[3]))) {
    r <- suppressMessages(evSyn(obj, study_weights = sw_ev, priorICweights = pr))
    l1 <- leave1studyout(r)
    ref <- t(sapply(seq_len(S_ev), function(k) {
      final_w(suppressMessages(evSyn(obj[-k], study_weights = sw_ev[-k], priorICweights = pr)))
    }))
    expect_equal(unname(l1$OverallGoricaWeights), unname(ref), tolerance = 1e-8)
  }
})

# -----------------------------------------------------------------------------
# Variability measures and conclusion
# -----------------------------------------------------------------------------
test_that("variability: geometric mean, ER, and SR", {
  r <- fit_ev(type_ev = "average")
  v <- evSyn_variability(r)
  W <- r$GORICA_weight_m
  R12 <- W[, 1] / W[, 2]
  gm <- exp(mean(log(R12)))
  expect_equal(unname(v$ratios_rates[1, "Geom. mean"]), unname(r$Final_ratio_GORICA_weights[1, 2]),
               tolerance = 1e-8)
  expect_equal(unname(v$ratios_rates[1, "ER"]), if (gm >= 1) mean(R12 > 1) else mean(R12 < 1))
  expect_equal(unname(v$ratios_rates[1, "SR"]), if (gm >= 1) mean(R12 > gm) else mean(R12 < gm))
  # weighted versions
  rw <- fit_ev(type_ev = "average", study_weights = sw_ev, order_studies = "descending")
  vw <- summary(rw, variability = TRUE)$variability
  w <- sw_ev[rw$order_studies] / sum(sw_ev)
  Rw <- rw$GORICA_weight_m[, 1] / rw$GORICA_weight_m[, 2]
  expect_equal(unname(vw$ratios_rates[1, "Geom. mean"]), unname(rw$Final_ratio_GORICA_weights[1, 2]),
               tolerance = 1e-8)
  expect_equal(unname(vw$ratios_rates[1, "ER"]), sum(w * (Rw > 1)), tolerance = 1e-10)
  # no variability section with one study
  expect_null(summary(suppressMessages(evSyn(est_ev[1], VCOV = vcov_ev[1], hypotheses = H_ev)),
                      variability = TRUE)$variability)
  expect_message(v1 <- evSyn_variability(suppressMessages(evSyn(est_ev[1], VCOV = vcov_ev[1], hypotheses = H_ev))),
                 "at least two")
  expect_null(v1)
})

test_that("conclusion: two steps with a failsafe; assumption text without a failsafe", {
  # by default (variability = FALSE), no variability section in print and summary
  expect_false(any(grepl("^Conclusion \\(incl. variability", capture.output(print(summary(fit_ev()))))))
  expect_false(any(grepl("^Variability over", capture.output(print(fit_ev())))))
  # with variability = TRUE
  out <- capture.output(print(summary(fit_ev(variability = TRUE))))
  expect_true(any(grepl("^Conclusion \\(incl. variability", out)))
  expect_true(any(grepl("Step 1", out)))
  expect_true(any(grepl("ER >= 0.90 and 0.40 <= SR <= 0.60", out)))
  # blank line before the Ceiling note
  i <- grep("\\(Ceiling:", out)
  expect_true(length(i) == 1 && out[i - 1] == "")
  # print() shows the variability section and conclusion as well
  out_p <- capture.output(print(fit_ev(variability = TRUE)))
  # the separate function gives the same section
  out_v <- capture.output(print(evSyn_variability(fit_ev())))
  expect_true(any(grepl("^Conclusion \\(incl. variability", out_v)))
  expect_true(any(grepl("^Conclusion \\(incl. variability", out_p)))
  # no failsafe in the set (IC weights input with default names)
  W <- rows_to_list(fit_ev()$GORICA_weight_m)
  out2 <- paste(capture.output(print(summary(evSyn(W, variability = TRUE)))), collapse = " ")
  expect_match(out2, "cover all possibilities")
})

# -----------------------------------------------------------------------------
# Confidence intervals: evSyn_CI() and evSyn(..., CI = TRUE)
# -----------------------------------------------------------------------------
test_that("evSyn_CI(): point values equal the final results of evSyn", {
  r <- fit_ev(type_ev = "equal", study_weights = sw_ev)
  ci <- evSyn_CI(r, B = 100, seed = 1)
  expect_s3_class(ci, "evSyn_CI")
  expect_equal(unname(ci$IC$gorica), unname(r$Cumulative_GORICA["Final", ]), tolerance = 1e-8)
  expect_equal(unname(ci$weights$gorica.weights), final_w(r), tolerance = 1e-8)
  expect_true(all(ci$weights$lower <= ci$weights$gorica.weights + 1e-12))
})

test_that("evSyn_CI() with one study equals goric_CI()", {
  r1 <- suppressMessages(evSyn(est_ev[1], VCOV = vcov_ev[1], hypotheses = H_ev))
  g1 <- goric(est_ev[[1]], VCOV = vcov_ev[[1]], hypotheses = H_ev)
  a <- evSyn_CI(r1, method = "delta")
  b <- goric_CI(g1, method = "delta")
  expect_equal(a$weights$lower, b$weights$lower, tolerance = 1e-8)
  expect_equal(a$ratio_weights$upper, b$ratio_weights$upper, tolerance = 1e-8)
})

test_that("evSyn_CI() is not available for input without estimates", {
  expect_error(evSyn_CI(evSyn(rows_to_list(fit_ev()$GORICA_weight_m))), "only available")
})

test_that("evSyn(..., CI = TRUE) adds CIs and leaves the random number generator as without CIs", {
  set.seed(9); invisible(fit_ev()); a <- runif(1)
  set.seed(9); r <- fit_ev(CI = TRUE, CI_B = 50); b <- runif(1)
  expect_identical(a, b)
  expect_s3_class(r$CI, "evSyn_CI")
  out <- capture.output(print(r))
  expect_true(any(grepl("95%CI", out)))
  out_s <- capture.output(print(summary(r)))
  expect_true(any(grepl("GORICA weights: 95%CI lower", out_s)))
  # message for input without information about the sampling variability
  expect_message(evSyn(rows_to_list(r$GORICA_weight_m), CI = TRUE), "only available")
})

test_that("CIs for goric objects as input", {
  set.seed(42)
  fits <- lapply(1:3, function(i) {
    n <- 60; g <- factor(rep(1:3, each = n / 3)); y <- c(0, .3, .6)[g] + rnorm(n)
    lm(y ~ -1 + g)
  })
  gor <- lapply(fits, function(f) goric(f, hypotheses = list(H1 = "g1 < g2 < g3")))
  r <- evSyn(gor, CI = TRUE, CI_B = 100)
  expect_s3_class(r$CI, "evSyn_CI")
  expect_equal(unname(r$CI$weights$goric.weights), final_w(r), tolerance = 1e-8)
  # one GORIC object: delta-method CIs equal goric_CI()
  a <- evSyn_CI(evSyn(gor[1]), method = "delta")
  b <- goric_CI(gor[[1]], method = "delta")
  expect_equal(a$ratio_weights$lower, b$ratio_weights$lower, tolerance = 1e-8)
  # GORICA objects give the same as the estimates + VCOV input
  gs <- lapply(seq_len(S_ev), function(s) goric(est_ev[[s]], VCOV = vcov_ev[[s]], hypotheses = H_ev))
  c1 <- evSyn_CI(evSyn(gs), B = 100, seed = 7)
  c2 <- evSyn_CI(fit_ev(), B = 100, seed = 7)
  expect_equal(c1$weights[, 2:4], c2$weights[, 2:4], tolerance = 1e-8, ignore_attr = TRUE)
})

# -----------------------------------------------------------------------------
# evSyn_forest()
# -----------------------------------------------------------------------------
test_that("evSyn_forest() returns a ggplot with the study-specific and final ratios", {
  r <- fit_ev(study_weights = sw_ev, order_studies = "descending")
  p <- evSyn_forest(r, vs = "H2", B = 100, seed = 1)
  expect_s3_class(p, "ggplot")
  d <- attr(p, "forest_data")
  expect_equal(nrow(d), S_ev + 1)
  W <- r$GORICA_weight_m
  expect_equal(unname(d$ratio[d$kind == "study"]), unname(W[, 1] / W[, 2]), tolerance = 1e-8)
  # average: final ratio without CI
  pa <- evSyn_forest(fit_ev(type_ev = "average"), B = 50, seed = 1)
  da <- attr(pa, "forest_data")
  expect_true(all(is.na(da$lower[da$kind == "pooled"])))
  expect_error(evSyn_forest(evSyn(rows_to_list(W))), "only available")
})

# -----------------------------------------------------------------------------
# Centered bootstrap, scientific notation, prior weights
# -----------------------------------------------------------------------------
test_that("evSyn_CI(): boot_type 'auto' uses the centered bootstrap from 5 studies onwards", {
  r5 <- fit_ev()                                   # 5 studies
  ci5 <- evSyn_CI(r5, B = 100, seed = 1)
  expect_identical(ci5$boot_type, "centered")
  expect_true(all(ci5$ratio_weights$CI_method == "boot (centered)"))
  r3 <- suppressMessages(evSyn(est_ev[1:3], VCOV = vcov_ev[1:3], hypotheses = H_ev))
  expect_identical(evSyn_CI(r3, B = 100, seed = 1)$boot_type, "percentile")
  ci_p <- evSyn_CI(r5, B = 100, seed = 1, boot_type = "percentile")
  expect_identical(ci_p$boot_type, "percentile")
  # the point values do not depend on the type of CI
  expect_equal(ci5$weights$gorica.weights, ci_p$weights$gorica.weights)
  expect_true(all(ci5$weights$lower <= ci5$weights$upper))
  expect_true("n_same_loglik" %in% names(ci5$ratio_weights))
})

test_that("CI tables print large values in scientific notation", {
  expect_identical(restriktor:::.format_CI_number(123456.789, digits = 3), "1.23e+05")
  expect_identical(restriktor:::.format_CI_number(12.3456, digits = 3), "12.346")
  r <- fit_ev()
  out <- capture.output(print(evSyn_CI(r, B = 100, seed = 1)))
  expect_false(any(grepl("[0-9]{6,}\\.[0-9]", out)))
})

test_that("goric_CI(): same_loglik and the exact lower bound at s = 0 (equality-only hypothesis)", {
  set.seed(3)
  g <- goric(c(b1 = 0.01, b2 = -0.02, b3 = 0.015), VCOV = diag(0.04, 3),
             hypotheses = list(H1 = "b1 = b2 = b3"), comparison = "unconstrained")
  ci <- goric_CI(g, method = "invert", B = 500)
  expect_true("same_loglik" %in% names(ci$difference))
  j <- which(ci$difference$model == "H1" & ci$difference$vs == "unconstrained")
  # exact expected difference at s = 0: 2 equalities + 2 x (PT_H1 - PT_u)
  exact0 <- 2 + 2 * (g$result$penalty[1] - g$result$penalty[2])
  expect_equal(ci$difference$lower[j], exact0, tolerance = 1e-10)
})

test_that("priorICweights: extra column in goric print/summary and row/column in evSyn output", {
  g <- goric(est_ev[[1]], VCOV = vcov_ev[[1]], hypotheses = H_ev, priorICweights = c(.5, .3, .2))
  out <- capture.output(print(g))
  expect_true(any(grepl("prior.gorica.weights", out, fixed = TRUE)))
  out_s <- capture.output(summary(g))
  expect_true(any(grepl("prior.gorica.weights", out_s, fixed = TRUE)))
  g0 <- goric(est_ev[[1]], VCOV = vcov_ev[[1]], hypotheses = H_ev)
  expect_false(any(grepl("prior.gorica.weights", capture.output(print(g0)), fixed = TRUE)))
  r <- fit_ev(priorICweights = c(.5, .3, .2))
  expect_true(any(grepl("prior GORICA weights", capture.output(print(r)))))
  expect_true("prior GORICA weights" %in% rownames(summary(r)$Final_Cumulative_results))
  expect_false(any(grepl("prior GORICA weights", capture.output(print(fit_ev())))))
})

test_that("prior IC weights are used once (aggregated results), not in the study-specific weights", {
  pr <- c(.5, .3, .2)
  r0 <- fit_ev(); r <- fit_ev(priorICweights = pr)
  # study-specific IC weights: the same as without prior IC weights
  expect_equal(unname(r$GORICA_weight_m), unname(r0$GORICA_weight_m), tolerance = 1e-12)
  # final weights: prior IC weights times the final weights without them (renormalized)
  w <- pr * final_w(r0); expect_equal(final_w(r), w / sum(w), tolerance = 1e-10)
  # added approach, equal study weights: final ratio = prior ratio x product of study-specific ratios
  W <- r$GORICA_weight_m
  expect_equal(final_w(r)[1] / final_w(r)[3], (pr[1] / pr[3]) * prod(W[, 1] / W[, 3]), tolerance = 1e-8)
  # same for IC weights as input
  rw <- suppressMessages(evSyn(rows_to_list(r0$GORICA_weight_m), priorICweights = pr))
  expect_equal(unname(rw$GORICA_weight_m), unname(r0$GORICA_weight_m), tolerance = 1e-10)
  expect_equal(unname(rw$Cumulative_GORICA_weights["Final", ]), final_w(r), tolerance = 1e-8)
})

test_that("priorICweights_study (not active for now) gives the earlier behaviour when switched on", {
  pr <- c(.5, .3, .2)
  r_new <- fit_ev(priorICweights = pr)
  ns <- asNamespace("restriktor")
  old <- get(".evSyn_priorICweights_study", envir = ns)
  assignInNamespace(".evSyn_priorICweights_study", TRUE, ns = "restriktor")
  on.exit(assignInNamespace(".evSyn_priorICweights_study", old, ns = "restriktor"))
  r_old <- fit_ev(priorICweights = pr)
  W <- sweep(r_new$GORICA_weight_m, 2, pr, "*"); W <- W / rowSums(W)
  expect_equal(unname(r_old$GORICA_weight_m), unname(W), tolerance = 1e-10)
  expect_equal(final_w(r_old), final_w(r_new), tolerance = 1e-12)  # aggregated results unchanged
  expect_true(isTRUE(r_old$priorICweights_study)); expect_false(isTRUE(r_new$priorICweights_study))
})

# -----------------------------------------------------------------------------
# Prediction interval, average-approach CIs, dependent studies
# -----------------------------------------------------------------------------
test_that("evSyn_PI() gives a prediction interval (per observation, ceiling, average)", {
  r <- fit_ev(study_sample_nobs = rep(100, S_ev))
  pi <- evSyn_PI(r, B = 200, seed = 1)
  expect_s3_class(pi, "evSyn_PI")
  expect_false(anyNA(pi$PI$lower))
  expect_true(all(pi$PI$lower <= pi$PI$ratio_typical + 1e-12 & pi$PI$ratio_typical <= pi$PI$upper + 1e-12))
  expect_true(all(pi$PI$P_favours_best >= 0 & pi$PI$P_favours_best <= 1))
  expect_true(pi$per_obs)
  # ceiling: the PI of the best (nested) hypothesis versus the unconstrained one is capped at its maximum
  j <- which(tolower(pi$PI$vs) == "unconstrained")
  g1 <- r$goric_studies[[1]]
  ceil <- exp(-(g1$result$penalty[1] - g1$result$penalty[3]))
  expect_true(pi$PI$upper[j] <= ceil * (1 + 1e-8))
  # no sample sizes: raw scale (with a note when printed)
  pi0 <- evSyn_PI(fit_ev(), B = 200, seed = 1)
  expect_false(pi0$per_obs)
  expect_output(suppressMessages(print(pi0)), "prediction interval")
  # wider for 'observed' than for 'true'
  pit <- evSyn_PI(r, B = 200, seed = 1, type = "true")
  k <- which(tolower(pi$PI$vs) != "unconstrained")[1]
  expect_true(pi$PI$upper[k] / pi$PI$lower[k] >= pit$PI$upper[k] / pit$PI$lower[k] - 1e-8)
  # average approach: another analysis of the same data
  pa <- evSyn_PI(fit_ev(type_ev = "average"))
  expect_identical(pa$type, "another analysis")
})

test_that("evSyn_CI() for the average approach: rho, between-analysis SDs", {
  ra <- fit_ev(type_ev = "average")
  c1 <- suppressMessages(evSyn_CI(ra, seed = 1))
  c0 <- suppressMessages(evSyn_CI(ra, rho = 0, seed = 1))
  expect_identical(c1$method, "mc"); expect_equal(c1$rho, 1)
  expect_true(all(c1$IC$SE >= c0$IC$SE - 1e-12))
  expect_true("SD_between" %in% names(c1$IC)); expect_true("SD_between_log" %in% names(c1$ratio_weights))
  # rho = 1, one dimension: SE of the mean = mean of the SEs
  S1 <- lapply(ra$goric_studies, function(g) suppressMessages(goric_CI(g, method = "delta"))$Sigma)
  se_mean <- 2 * mean(sapply(S1, function(S) sqrt(S[1, 1])))
  expect_equal(c1$IC$SE[1], se_mean, tolerance = 1e-8)
  # evSyn(..., CI = TRUE) works for the average approach
  rci <- suppressMessages(fit_ev(type_ev = "average", CI = TRUE))
  expect_false(is.null(rci$CI))
  expect_error(evSyn_CI(ra, rho = 2))
})

test_that("dependent studies: two-level synthesis (cluster) and RVE", {
  cl <- c(1, 1, 2, 3, 3)
  r  <- fit_ev()
  rc <- suppressMessages(fit_ev(cluster = cl))
  expect_s3_class(rc, "evSyn_cluster")
  expect_equal(rc$n_studies, 3)
  # cluster level = average of the studies within a cluster (LL), then added over clusters
  LL <- r$LL_m[match(1:5, r$order_studies), ]
  LLc <- rbind(colMeans(LL[1:2, ]), LL[3, ], colMeans(LL[4:5, ]))
  expect_equal(unname(colSums(LLc)), unname(colSums(rc$LL_m)), tolerance = 1e-10)
  # each study its own cluster: same as without clusters
  r1 <- suppressMessages(fit_ev(cluster = 1:5))
  expect_equal(final_w(r1), final_w(r), tolerance = 1e-10)
  # CIs: mc with rho; RVE keeps the point evidence of the studies as such
  ci <- suppressMessages(evSyn_CI(rc, seed = 1))
  expect_identical(ci$dependence, "cluster")
  rr <- suppressMessages(fit_ev(cluster = cl, dependence = "rve"))
  expect_equal(final_w(rr), final_w(r), tolerance = 1e-10)
  cr <- suppressMessages(evSyn_CI(rr, seed = 1))
  expect_identical(cr$dependence, "rve"); expect_equal(cr$df, 2)
  # PI with clusters: the clusters are the units
  pc <- evSyn_PI(rc, B = 100, seed = 1)
  expect_true(pc$clusters)
  expect_error(suppressMessages(fit_ev(type_ev = "average", cluster = cl)))
  expect_error(suppressMessages(fit_ev(cluster = c(1, 2))))
})

test_that("multiple imputation (MI = TRUE): average results, Rubin's rules for the CIs", {
  ra  <- fit_ev(type_ev = "average")
  rmi <- suppressMessages(fit_ev(type_ev = "average", MI = TRUE, CI = TRUE))
  expect_true(isTRUE(rmi$MI))
  expect_equal(final_w(rmi), final_w(ra), tolerance = 1e-12)       # same aggregated results
  ci <- rmi$CI
  expect_true(ci$MI); expect_true("df" %in% names(ci$difference))
  # Rubin: T = W + (1 + 1/m) B (on the LL scale), for the first IC value
  m <- S_ev
  W <- mean(sapply(ra$goric_studies, function(g) suppressMessages(goric_CI(g, method = "delta"))$Sigma[1, 1]))
  Lt <- sapply(ra$goric_studies, function(g) g$result$loglik[1] - g$result$penalty[1])
  expect_equal(ci$IC$SE[1], 2 * sqrt(W + (1 + 1 / m) * var(Lt)), tolerance = 1e-8)
  expect_error(suppressMessages(fit_ev(MI = TRUE)))                  # only for type_ev = "average"
  expect_error(suppressMessages(evSyn_CI(fit_ev(), MI = TRUE)))
})

test_that("evSyn_PI(): heterogeneity measures (tau, I2, Q) on the per-observation scale", {
  r <- fit_ev(study_sample_nobs = c(50, 100, 100, 200, 400))
  pi <- evSyn_PI(r, B = 200, seed = 1)
  h <- pi$heterogeneity
  expect_true(all(c("tau", "tau_logratio", "I2", "Q", "Q_df", "Q_p", "n_med") %in% names(h)))
  expect_equal(h$n_med[1], 100)
  ok <- h$k >= 2
  expect_true(all(h$I2[ok] >= 0 & h$I2[ok] <= 1))
  expect_equal(h$tau_logratio[ok], h$n_med[ok] * h$tau[ok])
  expect_true(all(h$Q_p[ok] >= 0 & h$Q_p[ok] <= 1))
  expect_output(print(pi), "median sample size")
})

test_that("evSyn_PI(): heterogeneity based on the signed distance uses all studies", {
  r <- fit_ev(study_sample_nobs = rep(100, S_ev))
  pi <- evSyn_PI(r, B = 200, seed = 1)
  hs <- pi$heterogeneity_sd
  expect_false(is.null(hs))
  expect_equal(hs$k, S_ev)                                  # all studies
  expect_equal(hs$n_inside + hs$n_outside <= S_ev, TRUE)
  expect_true(hs$I2 >= 0 && hs$I2 <= 1)
  expect_equal(hs$tau_distance, sqrt(hs$n_med) * hs$tau)
  # signed distance: positive inside, negative outside (study 1)
  g <- r$goric_studies[[1]]
  gc <- suppressMessages(goric(g$b.unrestr, VCOV = g$VCOV, hypotheses = list(H = H_ev[[pi$best]]),
                               comparison = "complement"))
  inside <- gc$result$loglik[1] > gc$result$loglik[2]
  expect_true(is.logical(inside))
  expect_output(suppressMessages(print(pi)), "signed distance")
})
