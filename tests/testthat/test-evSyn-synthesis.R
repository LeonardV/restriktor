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
  v <- summary(r)$variability
  W <- r$GORICA_weight_m
  R12 <- W[, 1] / W[, 2]
  gm <- exp(mean(log(R12)))
  expect_equal(unname(v$ratios_rates[1, "Geom. mean"]), unname(r$Final_ratio_GORICA_weights[1, 2]),
               tolerance = 1e-8)
  expect_equal(unname(v$ratios_rates[1, "ER"]), if (gm >= 1) mean(R12 > 1) else mean(R12 < 1))
  expect_equal(unname(v$ratios_rates[1, "SR"]), if (gm >= 1) mean(R12 > gm) else mean(R12 < gm))
  # weighted versions
  rw <- fit_ev(type_ev = "average", study_weights = sw_ev, order_studies = "descending")
  vw <- summary(rw)$variability
  w <- sw_ev[rw$order_studies] / sum(sw_ev)
  Rw <- rw$GORICA_weight_m[, 1] / rw$GORICA_weight_m[, 2]
  expect_equal(unname(vw$ratios_rates[1, "Geom. mean"]), unname(rw$Final_ratio_GORICA_weights[1, 2]),
               tolerance = 1e-8)
  expect_equal(unname(vw$ratios_rates[1, "ER"]), sum(w * (Rw > 1)), tolerance = 1e-10)
  # no variability section with one study
  expect_null(summary(suppressMessages(evSyn(est_ev[1], VCOV = vcov_ev[1], hypotheses = H_ev)))$variability)
})

test_that("conclusion: two steps with a failsafe; assumption text without a failsafe", {
  out <- capture.output(print(summary(fit_ev())))
  expect_true(any(grepl("^Conclusion \\(incl. variability", out)))
  expect_true(any(grepl("Step 1", out)))
  expect_true(any(grepl("ER >= 0.90 and 0.40 <= SR <= 0.60", out)))
  # blank line before the Ceiling note
  i <- grep("\\(Ceiling:", out)
  expect_true(length(i) == 1 && out[i - 1] == "")
  # print() shows the variability section and conclusion as well
  out_p <- capture.output(print(fit_ev()))
  expect_true(any(grepl("^Conclusion \\(incl. variability", out_p)))
  # no failsafe in the set (IC weights input with default names)
  W <- rows_to_list(fit_ev()$GORICA_weight_m)
  out2 <- paste(capture.output(print(summary(evSyn(W)))), collapse = " ")
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

test_that("evSyn_CI() is not available for type_ev = 'average' or input without estimates", {
  expect_error(evSyn_CI(fit_ev(type_ev = "average")), "not available for type_ev = 'average'")
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
