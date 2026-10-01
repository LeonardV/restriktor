# =============================================================================
# Tests: benchmark() (benchmark_means / benchmark_asymp), print en plot
# =============================================================================
# Smoke- en regressietests; klein aantal iteraties om de looptijd te beperken.

# Onderdrukt de voortgangsoutput (cat/progressbar) en messages van benchmark()
quiet <- function(expr) {
  res <- NULL
  utils::capture.output(res <- suppressMessages(expr))
  res
}

# --- Gedeelde testdata ---
set.seed(42)
n_g <- c(20, 25, 30)
df_bm <- data.frame(group = factor(rep(1:3, times = n_g)), x = rnorm(sum(n_g)))
df_bm$y <- c(0, 0.4, 0.8)[df_bm$group] + 0.5 * df_bm$x + rnorm(sum(n_g))
fit_bm <- lm(y ~ -1 + group, data = df_bm)
h1_bm <- "group1 < group2 < group3"

g_compl <- goric(fit_bm, hypotheses = list(H1 = h1_bm), comparison = "complement")
g_unc   <- goric(fit_bm, hypotheses = list(H1 = h1_bm), comparison = "unconstrained")

# controleer dat tabellen, draws en rates per populatie dezelfde hypothesen bevatten
expect_aligned <- function(b) {
  for (p in names(b$benchmarks$ratio_goric_weights)) {
    rn <- rownames(b$benchmarks$ratio_goric_weights[[p]])
    expect_equal(rn, paste(b$pref_hypo_name, colnames(b$combined_values$rgw_combined[[p]])))
    expect_equal(rownames(b$benchmarks$ratio_ll_weights[[p]]),
                 paste(b$pref_hypo_name, colnames(b$combined_values$rlw_combined[[p]])))
    expect_equal(rownames(b$benchmarks$difLL[[p]]),
                 paste(b$pref_hypo_name, colnames(b$combined_values$ld_combined[[p]])))
    expect_length(b$hypothesis_rate[[p]], length(rn))
    expect_length(b$rate_rlw[[p]], nrow(b$benchmarks$ratio_ll_weights[[p]]))
  }
}


test_that("benchmark_asymp: lm met complement en unconstrained", {
  b1 <- quiet(benchmark(g_compl, model_type = "asymp", iter = 30, seed = 1))
  expect_s3_class(b1, c("benchmark_asymp", "benchmark"))
  expect_equal(b1$pref_hypo_name, "H1")
  expect_equal(b1$iter, 30)
  expect_equal(rownames(b1$benchmarks$ratio_goric_weights[[1]]), "H1 vs. complement")
  expect_aligned(b1)

  b2 <- quiet(benchmark(g_unc, model_type = "asymp", iter = 30, seed = 1))
  expect_equal(rownames(b2$benchmarks$ratio_goric_weights[[1]]), "H1 vs. unconstrained")
  expect_aligned(b2)
})


test_that("benchmark_means: ANCOVA met continue covariaat", {
  fit_ancova <- lm(y ~ -1 + group + x, data = df_bm)
  g_ancova <- goric(fit_ancova, hypotheses = list(H1 = h1_bm), comparison = "complement")
  b <- quiet(benchmark(g_ancova, model_type = "means", iter = 30, seed = 1))
  expect_s3_class(b, "benchmark_means")
  # groepen = alleen factorcellen (geen covariaat)
  expect_equal(b$ngroups, 3)
  expect_equal(unname(b$group_size), n_g)
  # covariaat blijft op de geschatte waarde in elke populatie
  expect_equal(unname(b$pop_group_means[, "x"]),
               rep(unname(coef(fit_ancova)["x"]), 2), tolerance = 1e-6)
  expect_aligned(b)
})


test_that("benchmark_means: schattingen + VCOV vereist group_size", {
  g_est <- goric(coef(fit_bm), VCOV = vcov(fit_bm), hypotheses = list(H1 = h1_bm),
                 comparison = "complement")
  expect_error(quiet(benchmark(g_est, model_type = "means", iter = 20, seed = 1)),
               "group_size")
  b <- quiet(benchmark(g_est, model_type = "means", iter = 20, seed = 1,
                       group_size = n_g))
  expect_s3_class(b, "benchmark_means")
  expect_equal(unname(b$group_size), n_g)
  expect_error(quiet(benchmark(g_est, model_type = "means", iter = 20,
                               group_size = c(20, 25))), "group_size")
})


test_that("benchmark_means: scalaire alt_group_size gelijk aan rep(alt, ngroups)", {
  b_s <- quiet(benchmark(g_compl, model_type = "means", iter = 20, seed = 1,
                         alt_group_size = 50))
  b_v <- quiet(benchmark(g_compl, model_type = "means", iter = 20, seed = 1,
                         alt_group_size = rep(50, 3)))
  expect_equal(b_s$cohens_f_observed, b_v$cohens_f_observed)
  expect_equal(b_s$benchmarks, b_v$benchmarks)
  expect_equal(unname(b_s$group_size), rep(50, 3))
  expect_error(quiet(benchmark(g_compl, model_type = "means", iter = 20,
                               alt_group_size = c(50, 60))), "alt_group_size")
})


test_that("benchmark: priorICweights worden doorgegeven (pref_hypo_name)", {
  hypos <- list(H1 = h1_bm, H2 = "group1 = group2 = group3")
  g_prior <- goric(fit_bm, hypotheses = hypos, comparison = "none",
                   priorICweights = c(0.001, 0.999))
  best <- g_prior$result$model[which.max(g_prior$result[, 7])]
  g_flat <- goric(fit_bm, hypotheses = hypos, comparison = "none")
  # de priors veranderen de voorkeurshypothese
  expect_false(identical(best, g_flat$result$model[which.max(g_flat$result[, 7])]))

  b_m <- quiet(benchmark(g_prior, model_type = "means", iter = 20, seed = 1))
  expect_equal(b_m$pref_hypo_name, best)
  b_a <- quiet(benchmark(g_prior, model_type = "asymp", iter = 20, seed = 1))
  expect_equal(b_a$pref_hypo_name, best)
})


test_that("print.benchmark: output_type = 'all' en percentiles", {
  b <- quiet(benchmark(g_unc, model_type = "asymp", iter = 30, seed = 1))
  expect_no_warning(out <- utils::capture.output(
    print(b, output_type = "all", color = FALSE)))
  expect_true(any(grepl("Hypothesis rate", out)))
  expect_no_warning(out2 <- utils::capture.output(
    print(b, output_type = "all", color = FALSE, percentiles = c(0.1, 0.9))))
  expect_true(any(grepl("10%", out2)))
  expect_no_warning(utils::capture.output(
    print(b, output_type = "rgw", color = FALSE, percentiles = 0.5)))
})


test_that("print/plot: geen uitlijningsfouten bij (vrijwel) constante draws", {
  # 'strong' populatie: H1 is (vrijwel) altijd waar, dus rlw = 1 en ld = 0
  # voor bijna alle draws; deze rijen mogen niet (alleen in de tabel of
  # alleen in de draws) wegvallen
  set.seed(11)
  d2 <- data.frame(value = c(rnorm(30, 0, 1), rnorm(30, .5, 1), rnorm(30, .45, 1)),
                   group = factor(rep(1:3, each = 30)))
  f2 <- lm(value ~ -1 + group, d2)
  g2 <- goric(f2, hypotheses = list(H1 = h1_bm), comparison = "unconstrained",
              type = "gorica")
  pe <- rbind(null = c(0, 0, 0), strong = c(0, 5, 10))
  b <- quiet(benchmark(g2, model_type = "asymp", seed = 1, iter = 30, pop_est = pe))
  expect_aligned(b)
  for (ot in c("rlw", "ld", "all")) {
    expect_no_warning(expect_no_error(utils::capture.output(
      print(b, output_type = ot, color = FALSE, percentiles = c(.1, .9)))))
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (ot in c("rgw", "rlw", "gw", "ld")) {
    p <- plot(b, output_type = ot)
    expect_s3_class(p, "benchmark_plot")
    expect_no_error(print(p))
  }
})


test_that("plot.benchmark: alle output_types", {
  b <- quiet(benchmark(g_compl, model_type = "means", iter = 30, seed = 1))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (ot in c("rgw", "rlw", "gw", "ld")) {
    p <- plot(b, output_type = ot, opacity = 0.3)
    expect_s3_class(p, "benchmark_plot")
    expect_no_error(print(p))
  }
  expect_warning(plot(b, alpha = 0.3), "opacity")
})


test_that("adaptieve iter: stopt pas na twee stabiele rondes", {
  # met een ruime tolerantie is de percentiel direct stabiel: stoppen na
  # iter_min + 2 * iter_step draws (twee stabiele rondes), niet eerder
  b <- quiet(benchmark(g_compl, model_type = "asymp", seed = 1, iter_min = 30,
                       iter_step = 20, iter_max = 200, iter_stability_tol = 100))
  expect_equal(b$iter, 30 + 2 * 20)
  # tolerantie 0: nooit stabiel, dus tot iter_max (laatste batch afgekapt)
  msgs <- character(0)
  utils::capture.output(b0 <- withCallingHandlers(
    benchmark(g_compl, model_type = "asymp", seed = 1, iter_min = 30,
              iter_step = 20, iter_max = 75, iter_stability_tol = 0),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  expect_equal(b0$iter, 75)
  expect_true(any(grepl("iter_stability_tol = 0", msgs)))
})


test_that("adaptieve iter: ongeldige argumenten geven een duidelijke fout", {
  expect_error(benchmark(g_compl, iter_step = 0), "iter_step")
  expect_error(benchmark(g_compl, iter_min = 0), "iter_min")
  expect_error(benchmark(g_compl, iter_max = 0), "iter_max")
  expect_error(benchmark(g_compl, iter = 0), "iter")
  expect_error(benchmark(g_compl, iter_stability_tol = -1), "iter_stability_tol")
  expect_error(benchmark(g_compl, iter_adequacy_band = c(0.6, 0.4)), "iter_adequacy_band")
})


test_that("generate_scaled_means: volgorde behouden en geen NaN bij negatief/nul minimum", {
  gsm <- restriktor:::generate_scaled_means
  cf <- restriktor:::compute_cohens_f
  N <- c(10, 10, 10)
  V <- diag(3) * 0.1
  # negatief minimum: vroeger draaide de volgorde om
  m_neg <- c(a = -2, b = -1, c = 1)
  r_neg <- gsm(m_neg, target_f = 0.3, N, V)
  expect_true(all(is.finite(r_neg)))
  expect_equal(order(r_neg), order(m_neg))
  expect_equal(sign(r_neg - weighted.mean(r_neg, N)), sign(m_neg - weighted.mean(m_neg, N)))
  expect_equal(cf(r_neg, N, V), 0.3)
  # nul minimum: vroeger NaN/Inf
  m_zero <- c(a = 0, b = 1, c = 2)
  r_zero <- gsm(m_zero, target_f = 0.3, N, V)
  expect_true(all(is.finite(r_zero)))
  expect_equal(order(r_zero), order(m_zero))
  expect_equal(cf(r_zero, N, V), 0.3)
  # alleen een verschuiving maakt niet uit (zelfde patroon)
  expect_equal(gsm(c(3, 2, 1), 0.3, N, V), gsm(c(1, 0, -1), 0.3, N, V))
  # f = 0 geeft nullen; gelijke means met f > 0 geeft een duidelijke fout
  expect_equal(unname(gsm(m_neg, 0, N, V)), c(0, 0, 0))
  expect_error(gsm(c(1, 1, 1), 0.3, N, V), "all equal")
})


test_that("benchmark_means: negatieve groepsgemiddelden behouden de volgorde; Observed = schattingen", {
  set.seed(7)
  d_neg <- data.frame(group = factor(rep(1:3, times = n_g)))
  d_neg$y <- c(0.5, -0.6, 0.9)[d_neg$group] + rnorm(sum(n_g))
  f_neg <- lm(y ~ -1 + group, d_neg)
  expect_true(any(coef(f_neg) < 0))
  g_neg <- goric(f_neg, hypotheses = list(H1 = h1_bm), comparison = "complement")
  b <- quiet(benchmark(g_neg, model_type = "means", iter = 10, seed = 1))
  pm <- b$pop_group_means
  expect_true(all(is.finite(pm)))
  expect_equal(unname(pm["pop_es = 0", ]), c(0, 0, 0))
  obs_row <- grep("Observed", names(b$benchmarks$goric_weights), value = TRUE)
  expect_length(obs_row, 1)
  expect_equal(unname(pm[2, ]), unname(coef(f_neg)))
  expect_equal(order(pm[2, ]), order(coef(f_neg)))
  # ook een andere effectgrootte behoudt de volgorde van de data
  b2 <- quiet(benchmark(g_neg, model_type = "means", iter = 10, seed = 1,
                        pop_es = c(0.2, 0.5)))
  for (i in 1:2) {
    expect_equal(order(b2$pop_group_means[i, ]), order(coef(f_neg)))
    expect_equal(restriktor:::compute_cohens_f(b2$pop_group_means[i, ], b2$group_size,
                                               g_neg$VCOV), unname(b2$pop_es[i]))
  }
})


test_that("benchmark: goricac en gorica geven verschillende benchmarks", {
  est <- c(x = 1, y = 2, z = 3)
  V <- diag(3)
  H <- list(H = "x < y < z")
  ga <- goric(est, VCOV = V, hypotheses = H, type = "gorica")
  gc <- goric(est, VCOV = V, hypotheses = H, type = "goricac", sample_nobs = 8)
  ba <- quiet(benchmark(ga, iter = 20, seed = 1))
  bc <- quiet(benchmark(gc, iter = 20, seed = 1))
  expect_equal(ba$type, "gorica")
  expect_equal(bc$type, "goricac")
  expect_false(isTRUE(all.equal(ba$benchmarks$goric_weights, bc$benchmarks$goric_weights)))
  # de 'Sample'-waarde is op hetzelfde criterium gebaseerd als de draws
  expect_equal(bc$benchmarks$goric_weights[[1]][1, "Sample"],
               gc$result$goricac.weights[1])
  # goricac zonder steekproefgrootte: duidelijke fout
  gc_noN <- goric(est, VCOV = V, hypotheses = H, type = "gorica")
  gc_noN$type <- "goricc"
  expect_error(quiet(benchmark(gc_noN, iter = 10, seed = 1)), "sample_size")
  # means: goricc wordt goricac, met sum(N) als steekproefgrootte
  gcc <- goric(fit_bm, hypotheses = list(H1 = h1_bm), comparison = "complement",
               type = "goricc")
  bcc <- quiet(benchmark(gcc, model_type = "means", iter = 10, seed = 1))
  expect_equal(bcc$type, "goricac")
  bgm <- quiet(benchmark(g_compl, model_type = "means", iter = 10, seed = 1))
  expect_false(isTRUE(all.equal(bcc$benchmarks$goric_weights, bgm$benchmarks$goric_weights)))
})


test_that("benchmark: oneindige ratio (gewicht complement 0) draait, print en plot", {
  g_inf <- goric(c(a = 0, b = 20), VCOV = diag(2) * .01, hypotheses = list(H = "a < b"))
  expect_no_error(b <- quiet(benchmark(g_inf, iter = 10, seed = 2)))
  obs <- "pop_est = Observed"
  # de ratio zelf is Inf, maar de log-ratio is eindig
  expect_true(is.infinite(b$benchmarks$ratio_goric_weights[[obs]][1, "Sample"]))
  expect_true(all(is.finite(b$benchmarks$ratio_goric_weights_log[[obs]])))
  expect_true(all(is.finite(b$combined_values$rgw_log_combined[[obs]])))
  expect_equal(b$benchmarks$ratio_goric_weights_log[[obs]][1, "Sample"], 10000)
  # Inf telt mee als 'groter dan elke eindige waarde' in de rates/percentielen
  # (rgw = gewicht voorkeur / gewicht alternatief, dus Inf > 1 in elke draw)
  expect_equal(unname(b$hypothesis_rate[[obs]]), 1)
  expect_equal(unname(b$pctl_Sample$ratio_goric_weights[[obs]][1, 1]), 100)
  # overlap op basis van een dichtheid kan niet met Inf: NA, maar niet op de log-schaal
  expect_true(is.na(b$overlap$ratio_goric_weights[["pop_est = No-effect"]]))
  expect_false(is.na(b$overlap$ratio_goric_weights_log[["pop_est = No-effect"]]))
  expect_no_error(utils::capture.output(print(b, output_type = "all", color = FALSE)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (ot in c("rgw", "rlw", "gw", "ld")) {
    expect_no_error(print(plot(b, output_type = ot)))
  }
})


test_that("benchmark_means: ratio_pop_means bepaalt het patroon van de populatiegemiddelden", {
  b_123 <- quiet(benchmark(g_compl, model_type = "means", iter = 10, seed = 1,
                           ratio_pop_means = c(1, 2, 3)))
  b_321 <- quiet(benchmark(g_compl, model_type = "means", iter = 10, seed = 1,
                           ratio_pop_means = c(3, 2, 1)))
  expect_equal(unname(b_123$ratio_pop_means), c(1, 2, 3))
  expect_false(isTRUE(all.equal(b_123$pop_group_means, b_321$pop_group_means)))
  pm <- b_123$pop_group_means[2, ]
  # gevraagde verhoudingen: opeenvolgende verschillen gelijk, oplopend
  expect_equal(unname(diff(pm)), rep(unname(diff(pm))[1], 2))
  expect_true(all(diff(pm) > 0))
  expect_true(all(diff(b_321$pop_group_means[2, ]) < 0))
  expect_equal(unname(b_321$pop_group_means[2, ]), -unname(pm))
  # effectgrootte klopt
  expect_equal(restriktor:::compute_cohens_f(pm, b_123$group_size, g_compl$VCOV),
               unname(b_123$pop_es[2]))
  # verschuiving maakt niet uit
  b_shift <- quiet(benchmark(g_compl, model_type = "means", iter = 10, seed = 1,
                             ratio_pop_means = c(11, 12, 13)))
  expect_equal(b_shift$pop_group_means, b_123$pop_group_means)
  expect_equal(b_shift$benchmarks, b_123$benchmarks)
  # verkeerde lengte: duidelijke fout
  expect_error(quiet(benchmark(g_compl, model_type = "means", iter = 10,
                               ratio_pop_means = c(1, 2))), "ratio_pop_means")
})
