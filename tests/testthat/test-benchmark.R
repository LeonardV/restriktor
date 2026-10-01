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
