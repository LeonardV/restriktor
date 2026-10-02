# =============================================================================
# Tests: bevindingen uit de audit van de benchmark-functies (ronde 4)
# =============================================================================
# Regressietests met de reproducers uit de audit (A9-A12, B26-B35); klein
# aantal iteraties om de looptijd te beperken.

quiet <- function(expr) {
  res <- NULL
  utils::capture.output(res <- suppressMessages(expr))
  res
}

# vangt de messages van benchmark() op (de voortgangsoutput wordt weggegooid)
with_msgs <- function(expr) {
  msgs <- character(0)
  utils::capture.output(res <- withCallingHandlers(expr,
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  list(value = res, msgs = msgs)
}

# Onafhankelijke referentieformule voor Cohen's f
cohens_f_ref <- function(m, N, s2) {
  mu <- sum(N * m) / sum(N)
  unname(sqrt(sum(N * (m - mu)^2) / sum(N)) / sqrt(s2))
}

set.seed(42)
n_g <- c(20, 25, 30)
df_bm <- data.frame(group = factor(rep(1:3, times = n_g)), x = rnorm(sum(n_g)))
df_bm$y <- c(0, 0.4, 0.8)[df_bm$group] + 0.5 * df_bm$x + rnorm(sum(n_g))
fit_bm <- lm(y ~ -1 + group, data = df_bm)
h1_bm <- "group1 < group2 < group3"
g_compl <- goric(fit_bm, hypotheses = list(H1 = h1_bm), comparison = "complement")

# tracer op goric() die alleen binnen een benchmark-draw (parallel_function_asymp)
# iets doet; 'i' is het draw-nummer
in_draw_tracer <- function(expr) {
  bquote({
    .calls <- sys.calls()
    .k <- which(vapply(.calls, function(cl) {
      identical(cl[[1]], as.name("parallel_function_asymp"))
    }, logical(1)))
    if (length(.k) > 0) {
      i <- get("i", sys.frames()[[.k[length(.k)]]])
      .(expr)
    }
  })
}
trace_goric <- function(tracer) {
  suppressMessages(trace("goric", where = asNamespace("restriktor"), print = FALSE,
                         tracer = tracer))
}
untrace_goric <- function() {
  suppressMessages(untrace("goric", where = asNamespace("restriktor")))
}


# -----------------------------------------------------------------------------
# A9: warnings in een draw worden gedempt (draw blijft), alleen errors
# verwijderen een draw; iter telt de geslaagde draws
# -----------------------------------------------------------------------------

test_that("A9: benchmark draait met ':='-definities in de hypothesen", {
  gd <- goric(fit_bm, hypotheses = list(
    H1 = "d1 := group2 - group1; d2 := group3 - group2; d1 > 0; d2 > d1"),
    comparison = "complement")
  for (mt in c("means", "asymp")) {
    b <- quiet(benchmark(gd, model_type = mt, iter = 5, seed = 3))
    expect_s3_class(b, "benchmark")
    expect_equal(b$iter, 5)
    expect_equal(unname(b$n_failed_draws), c(0, 0))
    expect_equal(nrow(b$combined_values$rgw_combined[[1]]), 5)
    expect_no_error(utils::capture.output(print(b, color = FALSE)))
  }
})


test_that("A9: een warning in een draw dempt, maar de draw blijft behouden en wordt geteld", {
  trace_goric(in_draw_tracer(quote(warning("synthetische warning"))))
  on.exit(untrace_goric(), add = TRUE)
  r <- with_msgs(benchmark(g_compl, model_type = "asymp", iter = 6, seed = 1))
  b <- r$value
  # alle draws zijn behouden
  expect_equal(b$iter, 6)
  expect_equal(b$iter_requested, 6)
  expect_equal(unname(b$n_failed_draws), c(0, 0))
  expect_equal(unname(b$n_warned_draws), c(6, 6))
  expect_equal(nrow(b$combined_values$rgw_combined[[1]]), 6)
  expect_match(unname(b$draw_warnings[1]), "synthetische warning")
  # een samenvattende melding (aantal draws met warning + eerste melding)
  expect_true(any(grepl("12 of the 12 benchmark draws", r$msgs)))
  expect_true(any(grepl("synthetische warning", r$msgs)))
  out <- utils::capture.output(print(b, color = FALSE))
  expect_true(any(grepl("muffled", out)))
  # dezelfde draws als zonder warning (de warning verandert niets aan het resultaat)
  untrace_goric()
  b0 <- quiet(benchmark(g_compl, model_type = "asymp", iter = 6, seed = 1))
  expect_equal(b$combined_values$gw_combined, b0$combined_values$gw_combined)
  expect_equal(b$benchmarks, b0$benchmarks)
})


test_that("A9: een error in een draw verwijdert alleen die draw; iter telt de geslaagde draws", {
  trace_goric(in_draw_tracer(quote(if (i %% 4 == 0) stop("synthetische fout"))))
  on.exit(untrace_goric(), add = TRUE)
  r <- with_msgs(benchmark(g_compl, model_type = "asymp", iter = 8, seed = 1))
  b <- r$value
  expect_equal(b$iter, 6)
  expect_equal(b$iter_requested, 8)
  expect_equal(unname(b$n_failed_draws), c(2, 2))
  expect_match(unname(b$draw_errors[1]), "synthetische fout")
  expect_equal(unname(sapply(b$combined_values$gw_combined, length)), c(6, 6))
  expect_true(any(grepl("4 of the 16 benchmark draws", r$msgs)))
  expect_true(any(grepl("synthetische fout", r$msgs)))
  # print: aantal geslaagde draws en de reden
  out <- utils::capture.output(suppressMessages(print(b, color = FALSE)))
  expect_true(any(grepl("requested bootstrap draws: 8", out)))
  expect_true(any(grepl("successful bootstrap draws for pop_est = No-effect: 6", out)))
  expect_true(any(grepl("synthetische fout", out)))
  # percentielen op basis van de behouden draws
  draws <- b$combined_values$rgw_combined[[2]][, 1]
  expect_equal(unname(b$pctl_Sample$ratio_goric_weights[[2]][1, 1]),
               100 * mean(draws <= b$benchmarks$ratio_goric_weights[[2]][1, "Sample"]))
  # alle draws mislukt: duidelijke fout
  untrace_goric()
  trace_goric(in_draw_tracer(quote(stop("altijd fout"))))
  expect_error(quiet(benchmark(g_compl, model_type = "asymp", iter = 3, seed = 1)),
               "All 3 benchmark draws.*failed.*altijd fout")
})


test_that("A9/B32: adaptieve iter telt geslaagde draws; progressr-handlers blijven ongemoeid", {
  old_handlers <- progressr::handlers()
  on.exit(progressr::handlers(old_handlers), add = TRUE)
  progressr::handlers("void")
  before <- progressr::handlers()
  trace_goric(in_draw_tracer(quote(if (i %% 5 == 0) stop("synthetische fout"))))
  on.exit(untrace_goric(), add = TRUE)
  r <- with_msgs(benchmark(g_compl, model_type = "asymp", seed = 1, iter_min = 10,
                           iter_step = 5, iter_max = 20, iter_stability_tol = 100))
  b <- r$value
  expect_equal(b$iter_requested, 20)
  expect_equal(b$iter, 16)
  expect_true(any(grepl("Using iter = 16 draws", r$msgs)))
  untrace_goric()
  # de globale handlers zijn niet overschreven (B32)
  expect_identical(progressr::handlers(), before)
})


# -----------------------------------------------------------------------------
# A10 / B37: vector sample_size + alt_sample_size (benchmark_asymp)
# -----------------------------------------------------------------------------

test_that("A10: vector sample_size met alt_sample_size schaalt VCOV scalair en symmetrisch", {
  est <- c(a = 1, b = 2, c = 3)
  V <- diag(c(0.3, 0.2, 0.1))
  V[1, 2] <- V[2, 1] <- 0.05
  g <- goric(est, VCOV = V, hypotheses = list(H = "a < b < c"), comparison = "complement",
             type = "gorica")
  b <- quiet(benchmark(g, model_type = "asymp", iter = 5, seed = 1,
                       sample_size = c(20, 25, 30), alt_sample_size = 150))
  expect_equal(unname(b$pop_VCOV), V * 75 / 150)
  expect_true(isSymmetric(unname(b$pop_VCOV)))
  expect_equal(b$sample_size, 150)
  b_s <- quiet(benchmark(g, model_type = "asymp", iter = 5, seed = 1,
                         sample_size = 75, alt_sample_size = 150))
  expect_equal(b$benchmarks, b_s$benchmarks)
  # vector alt_sample_size wordt gesommeerd
  b_v <- quiet(benchmark(g, model_type = "asymp", iter = 5, seed = 1,
                         sample_size = c(20, 25, 30), alt_sample_size = c(50, 50, 50)))
  expect_equal(b_v$benchmarks, b_s$benchmarks)
  expect_error(quiet(benchmark(g, model_type = "asymp", iter = 5, seed = 1,
                               sample_size = c(20, 25, 30), alt_sample_size = -1)),
               "alt_sample_size")
  expect_error(quiet(benchmark(g, model_type = "asymp", iter = 5, seed = 1,
                               alt_sample_size = 150)), "sample_size")
})


test_that("B26: sample_size voor een lm-object wordt gebruikt, met een melding", {
  r <- with_msgs(benchmark(g_compl, model_type = "asymp", iter = 5, seed = 1,
                           sample_size = 30, alt_sample_size = 150))
  expect_true(any(grepl("'sample_size' \\(total: 30\\) differs", r$msgs)))
  expect_equal(unname(r$value$pop_VCOV), unname(g_compl$VCOV) * 30 / 150)
  # zonder afwijking geen melding; zonder sample_size: N uit het object (75)
  r2 <- with_msgs(benchmark(g_compl, model_type = "asymp", iter = 5, seed = 1,
                            sample_size = 75, alt_sample_size = 150))
  expect_false(any(grepl("differs", r2$msgs)))
  b3 <- quiet(benchmark(g_compl, model_type = "asymp", iter = 5, seed = 1,
                        alt_sample_size = 150))
  expect_equal(unname(b3$pop_VCOV), unname(g_compl$VCOV) * 75 / 150)
})


# -----------------------------------------------------------------------------
# A11 / B30: tweede factor en intercept-modellen (benchmark_means)
# -----------------------------------------------------------------------------

test_that("A11: alleen de eerste factor levert groepsgemiddelden; tweede factor is covariaat", {
  set.seed(3)
  df_bm$sex <- factor(sample(c("f", "m"), nrow(df_bm), TRUE))
  f2 <- lm(y ~ -1 + group + sex, data = df_bm)
  g2 <- goric(f2, hypotheses = list(H1 = h1_bm), comparison = "complement")
  r <- with_msgs(benchmark(g2, model_type = "means", iter = 5, seed = 1))
  b <- r$value
  expect_true(any(grepl("term 'group' \\(group1, group2, group3\\) are treated as the group means",
                        r$msgs)))
  expect_equal(b$ngroups, 3)
  expect_equal(unname(b$group_size), n_g)
  # Cohen's f over de 3 (gecorrigeerde) groepsgemiddelden, sigma2 = RSS/N
  s2 <- sum(residuals(f2)^2) / nrow(df_bm)
  expect_equal(b$res_var, s2)
  expect_equal(b$cohens_f_observed, cohens_f_ref(coef(f2)[1:3], n_g, s2))
  # het sekse-contrast blijft in elke populatie op de geschatte waarde
  expect_equal(unname(b$pop_group_means[, "sexm"]), rep(unname(coef(f2)["sexm"]), 2))
  expect_equal(unname(b$pop_group_means["pop_es = 0", 1:3]), c(0, 0, 0))
  b2 <- quiet(benchmark(g2, model_type = "means", iter = 5, seed = 1, pop_es = 0.3))
  expect_equal(unname(b2$pop_group_means[1, "sexm"]), unname(coef(f2)["sexm"]))
  expect_equal(cohens_f_ref(b2$pop_group_means[1, 1:3], n_g, s2), 0.3)
  # verkeerde groepsgroottelengte geeft een duidelijke fout (3 groepen, niet 4)
  expect_error(quiet(benchmark(g2, model_type = "means", iter = 5,
                               group_size = c(20, 25, 30, 40))), "group_size")
})


test_that("B30/A11: intercept-model geeft een duidelijke fout, ook via benchmark_means()", {
  fi <- lm(y ~ group, data = df_bm)
  gi <- goric(fi, hypotheses = list(H1 = "group2 < group3"), comparison = "complement")
  expect_error(quiet(benchmark_means(gi, iter = 5, seed = 1)), "intercept")
  expect_error(quiet(benchmark(gi, model_type = "means", iter = 5, seed = 1)), "intercept")
  # schattingen + VCOV: de foutmelding noemt de reden en 'group_size'
  g_est <- goric(coef(fit_bm), VCOV = vcov(fit_bm), hypotheses = list(H1 = h1_bm),
                 comparison = "complement")
  expect_error(quiet(benchmark(g_est, model_type = "means", iter = 5)),
               "only estimates and their covariance matrix.*group_size")
})


# -----------------------------------------------------------------------------
# A12: f-definitie == simulatievariantie (RSS/N), onafhankelijke referentie
# -----------------------------------------------------------------------------

test_that("A12: Cohen's f en de populatiegemiddelden horen bij de foutvariantie van de draws", {
  # onafhankelijk: RSS/N uit de data zelf
  rss <- sum((df_bm$y - ave(df_bm$y, df_bm$group))^2)
  s2 <- rss / nrow(df_bm)
  b <- quiet(benchmark(g_compl, model_type = "means", iter = 5, seed = 1, pop_es = 0.3))
  expect_equal(b$res_var, s2)
  # == de variantie waarmee de draws worden gegenereerd (N_g * VCOV[g, g])
  expect_equal(b$res_var, mean(n_g * diag(g_compl$VCOV)))
  expect_equal(cohens_f_ref(b$pop_group_means[1, ], n_g, s2), 0.3)
  # dus niet meer f * sqrt(N / (N - k)) t.o.v. de simulatie (vroeger sigma(fit)^2)
  expect_false(isTRUE(all.equal(cohens_f_ref(b$pop_group_means[1, ], n_g, sigma(fit_bm)^2), 0.3)))
  # kleine N: hetzelfde
  set.seed(2)
  d_small <- data.frame(group = factor(rep(1:3, each = 4)))
  d_small$y <- rnorm(12)
  f_small <- lm(y ~ -1 + group, d_small)
  g_small <- goric(f_small, hypotheses = list(H1 = h1_bm), comparison = "complement")
  bs <- quiet(benchmark(g_small, model_type = "means", iter = 5, seed = 1, pop_es = 0.3))
  expect_equal(bs$res_var, sum(residuals(f_small)^2) / 12)
  expect_equal(cohens_f_ref(bs$pop_group_means[1, ], rep(4, 3), bs$res_var), 0.3)
  expect_equal(unname(rep(4, 3) * diag(g_small$VCOV)), rep(bs$res_var, 3))
})


test_that("B27: Observed Effect-Size met de groepsgroottes van de data, ook met alt_group_size", {
  rss <- sum((df_bm$y - ave(df_bm$y, df_bm$group))^2)
  f_data <- cohens_f_ref(coef(fit_bm), n_g, rss / nrow(df_bm))
  b <- quiet(benchmark(g_compl, model_type = "means", iter = 5, seed = 1,
                       alt_group_size = c(100, 10, 10)))
  expect_equal(b$cohens_f_observed, f_data)
  expect_equal(b$cohens_f_alt_group_size, cohens_f_ref(coef(fit_bm), c(100, 10, 10), b$res_var))
  out <- utils::capture.output(print(b, color = FALSE))
  expect_true(any(grepl(sprintf("Observed Effect-Size \\(Cohens f\\): %.3f", f_data), out)))
  expect_true(any(grepl("with alt_group_size", out)))
  # alternatieve groepsgroottes gelijk aan die van de data: zelfde f, geen extra regel
  b2 <- quiet(benchmark(g_compl, model_type = "means", iter = 5, seed = 1, alt_group_size = n_g))
  expect_equal(b2$cohens_f_observed, f_data)
  expect_equal(b2$cohens_f_alt_group_size, f_data)
  expect_false(any(grepl("with alt_group_size", utils::capture.output(print(b2, color = FALSE)))))
})


# -----------------------------------------------------------------------------
# B28: control met 1 element; pop_est als data.frame
# -----------------------------------------------------------------------------

test_that("B28: een control-lijst met een optie wordt doorgegeven; pop_est als data.frame werkt", {
  log_env <- new.env()
  log_env$ctrl <- list()
  assign("bm_audit_log_env", log_env, envir = .GlobalEnv)
  on.exit(rm("bm_audit_log_env", envir = .GlobalEnv), add = TRUE)
  suppressMessages(trace("run_benchmark_simulation", where = asNamespace("restriktor"),
        tracer = quote({
          e <- get("bm_audit_log_env", envir = .GlobalEnv)
          e$ctrl[[length(e$ctrl) + 1]] <- control
        }), print = FALSE))
  on.exit(suppressMessages(untrace("run_benchmark_simulation", where = asNamespace("restriktor"))),
          add = TRUE)
  invisible(quiet(benchmark(g_compl, model_type = "asymp", iter = 2, seed = 1,
                            control = list(mix_weights_bootstrap = 50))))
  invisible(quiet(benchmark(g_compl, model_type = "asymp", iter = 2, seed = 1)))
  expect_equal(log_env$ctrl[[1]], list(mix_weights_bootstrap = 50))
  expect_equal(log_env$ctrl[[2]], g_compl$objectList[[1]]$control)
  # pop_est als data.frame == als matrix
  pe <- rbind(a = c(0, .5, 1), b = c(0, 0, 0))
  b_m <- quiet(benchmark(g_compl, model_type = "asymp", iter = 3, seed = 1, pop_est = pe))
  b_d <- quiet(benchmark(g_compl, model_type = "asymp", iter = 3, seed = 1,
                         pop_est = data.frame(pe)))
  expect_equal(b_d$pop_est, b_m$pop_est)
  expect_equal(b_d$benchmarks, b_m$benchmarks)
  expect_error(quiet(benchmark(g_compl, model_type = "asymp", iter = 3,
                               pop_est = data.frame(a = "x", b = "y", c = "z"))), "numeric")
})


# -----------------------------------------------------------------------------
# B31: NaN-gewichten
# -----------------------------------------------------------------------------

test_that("B31: NaN-gewichten (goricac met oneindige penalty) geven een duidelijke fout", {
  set.seed(5)
  d4 <- data.frame(group = factor(rep(1:4, each = 15)))
  d4$y <- c(0, .1, .2, .3)[d4$group] + rnorm(60)
  f4 <- lm(y ~ -1 + group, d4)
  # goric() weigert zelf al een te kleine N voor de goricac (core-fix)
  expect_error(suppressWarnings(goric(coef(f4), VCOV = vcov(f4),
                                      hypotheses = list(H1 = "group1 < group2 < group3 < group4"),
                                      comparison = "complement", type = "goricac",
                                      sample_nobs = 6)),
               "sample size")
  gc60 <- goric(coef(f4), VCOV = vcov(f4),
                hypotheses = list(H1 = "group1 < group2 < group3 < group4"),
                comparison = "complement", type = "goricac", sample_nobs = 60)
  expect_no_error(quiet(benchmark(gc60, iter = 5, seed = 1)))
  # object met NaN-gewichten (zoals vroeger bij een oneindige penalty): de
  # benchmark weigert dit met een duidelijke fout in plaats van NaN-tabellen
  gc <- gc60
  gc$result$goricac.weights[] <- NaN
  gc$result$goricac[] <- NaN
  expect_error(quiet(benchmark(gc, iter = 5, seed = 1)), "NA/NaN")
  expect_error(quiet(benchmark(gc, iter = 5, seed = 1)), "goricac weights of the GORIC\\(A\\) object")
})


# -----------------------------------------------------------------------------
# B33: 1 hypothese + comparison = "none"
# -----------------------------------------------------------------------------

test_that("B33: plot() bij 1 hypothese met comparison = 'none' crasht niet", {
  g1 <- goric(fit_bm, hypotheses = list(H1 = h1_bm), comparison = "none")
  b1 <- quiet(benchmark(g1, iter = 5, seed = 1))
  expect_equal(nrow(b1$benchmarks$ratio_goric_weights[[1]]), 0)
  expect_no_error(utils::capture.output(print(b1, color = FALSE)))
  for (ot in c("rgw", "rlw", "ld")) {
    expect_error(plot(b1, output_type = ot), "no alternative hypothesis.*'gw'")
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(print(plot(b1, output_type = "gw")))
})


# -----------------------------------------------------------------------------
# B34: dubbele populatienamen
# -----------------------------------------------------------------------------

test_that("B34: dubbele populatienamen worden uniek gemaakt (met melding)", {
  est <- c(a = .1, b = .3, c = .6)
  V <- diag(.04, 3)
  g <- goric(est, VCOV = V, hypotheses = list(H1 = "a < b < c"), comparison = "complement")
  r <- with_msgs(benchmark(g, iter = 3, seed = 1, pop_est = rbind(A = est, A = 2 * est)))
  expect_true(any(grepl("not unique", r$msgs)))
  expect_equal(rownames(r$value$pop_est), c("A", "A.1"))
  expect_equal(names(r$value$benchmarks$goric_weights), c("pop_est = A", "pop_est = A.1"))
  expect_equal(unname(sapply(r$value$combined_values$gw_combined, length)), c(3, 3))
  # de resultaten horen bij de juiste populatie (A.1 heeft 2x zo grote schattingen)
  b_sep <- quiet(benchmark(g, iter = 3, seed = 1, pop_est = rbind(B = 2 * est)))
  expect_equal(unname(r$value$pop_est["A.1", ]), unname(2 * est))
  # means: named duplicate pop_es
  r2 <- with_msgs(benchmark(g_compl, model_type = "means", iter = 3, seed = 1,
                            pop_es = c(A = 0.2, A = 0.3)))
  expect_equal(names(r2$value$pop_es), c("A", "A.1"))
  expect_true(any(grepl("not unique", r2$msgs)))
  expect_equal(unname(r2$value$pop_es), c(0.2, 0.3))
})


# -----------------------------------------------------------------------------
# B29 / B35: plot: percentiellijnen gelabeld; niet-eindige populaties
# -----------------------------------------------------------------------------

test_that("B29: de legenda van de percentiellijnen noemt de populatie", {
  b3 <- quiet(benchmark(g_compl, model_type = "means", iter = 10, seed = 1,
                        pop_es = c(0.1, 0.4)))
  p <- plot(b3, output_type = "rgw")
  labs <- levels(p$plots[[1]]$layers[[2]]$data$percentile_label)
  expect_true(any(grepl("5th Percentile \\(Effect-size = PES_1\\) = ", labs)))
  # de lijnen zijn de kwantielen van de eerste populatie
  q <- quantile(b3$combined_values$rgw_combined[[1]][, 1], c(.05, .95))
  expect_true(any(grepl(sprintf("5th Percentile \\(Effect-size = PES_1\\) = %.3f", q[1]), labs)))
  expect_true(any(grepl(sprintf("95th Percentile \\(Effect-size = PES_1\\) = %.3f", q[2]), labs)))
})


test_that("B35: plot() crasht niet als de eerste populatie alleen Inf-draws heeft", {
  est <- c(a = 0, b = 60)
  V <- diag(1e-4, 2)
  dimnames(V) <- list(names(est), names(est))
  g <- goric(est, VCOV = V, hypotheses = list(H1 = "a < b"), comparison = "complement")
  b <- quiet(benchmark(g, iter = 10, seed = 1, pop_est = rbind(Big = est, Null = c(30, 30))))
  expect_equal(sum(is.finite(b$combined_values$rgw_combined[["pop_est = Big"]])), 0)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (ot in c("rgw", "rlw")) {
    msgs <- character(0)
    p <- withCallingHandlers(plot(b, output_type = ot),
                             message = function(m) {
                               msgs <<- c(msgs, conditionMessage(m))
                               invokeRestart("muffleMessage")
                             })
    expect_true(any(grepl("Population estimates = Big", msgs)))
    expect_true(any(grepl("not shown", msgs)))
    expect_no_error(print(p))
    # percentiellijnen van de eerste populatie MET eindige draws (Null)
    labs <- levels(p$plots[[1]]$layers[[2]]$data$percentile_label)
    expect_true(any(grepl("Percentile \\(Population estimates = Null\\)", labs)))
  }
  # ook als de niet-eindige populatie niet de eerste is: melding
  b2 <- quiet(benchmark(g, iter = 10, seed = 1, pop_est = rbind(Null = c(30, 30), Big = est)))
  expect_message(plot(b2, output_type = "rgw"), "Big")
  # alle populaties niet-eindig: duidelijke fout
  b3 <- quiet(benchmark(g, iter = 5, seed = 1, pop_est = rbind(Big = est, Big2 = 2 * est)))
  expect_error(plot(b3, output_type = "rgw"), "rgw_log")
})


# -----------------------------------------------------------------------------
# B24: label van de referentiepopulatie bij ratio_pop_means
# -----------------------------------------------------------------------------

test_that("B24: print noemt bij ratio_pop_means dat de 'Observed' populatie dat patroon volgt", {
  b <- quiet(benchmark(g_compl, model_type = "means", iter = 5, seed = 1,
                       ratio_pop_means = c(1, 2, 3)))
  out <- utils::capture.output(print(b, color = FALSE))
  expect_true(any(grepl("Observed \\(Reference population; observed effect size, means from ratio_pop_means\\)", out)))
  b0 <- quiet(benchmark(g_compl, model_type = "means", iter = 5, seed = 1))
  out0 <- utils::capture.output(print(b0, color = FALSE))
  expect_true(any(grepl("Observed \\(Reference population\\)", out0)))
  expect_false(any(grepl("ratio_pop_means\\)", out0)))
})


# -----------------------------------------------------------------------------
# B11: dode code verwijderd
# -----------------------------------------------------------------------------

test_that("B11: parallel_function_means bestaat niet meer", {
  expect_false(exists("parallel_function_means", envir = asNamespace("restriktor")))
})


# -----------------------------------------------------------------------------
# W3-N-01: lm-object met afwijkende sample_nobs in goric(): de N van het model
# -----------------------------------------------------------------------------

test_that("W3-N-01: sigma2 en N komen van het model (nobs), niet van een afwijkende sample_nobs", {
  g_n <- suppressMessages(goric(fit_bm, hypotheses = list(H1 = h1_bm),
                                comparison = "complement", sample_nobs = 1000))
  expect_equal(g_n$sample_nobs, 1000)
  rss <- sum(residuals(fit_bm)^2)
  r <- with_msgs(benchmark(g_n, model_type = "means", iter = 5, seed = 1, pop_es = 0.3))
  b <- r$value
  expect_true(any(grepl("'sample_nobs' of the GORIC\\(A\\) object \\(1000\\) differs", r$msgs)))
  # res_var = RSS/N met N = 75 (niet RSS/1000)
  expect_equal(b$res_var, rss / nrow(df_bm))
  expect_equal(b$cohens_f_observed, cohens_f_ref(coef(fit_bm), n_g, rss / nrow(df_bm)))
  # de populatie met pop_es = 0.3 heeft f = 0.3 t.o.v. de foutvariantie van de draws
  s2_draws <- mean(n_g * diag(g_n$VCOV))
  expect_equal(cohens_f_ref(b$pop_group_means[1, ], n_g, s2_draws), 0.3)
  # identiek aan het object zonder sample_nobs
  b_ref <- quiet(benchmark(g_compl, model_type = "means", iter = 5, seed = 1, pop_es = 0.3))
  expect_equal(b$pop_group_means, b_ref$pop_group_means)
  expect_equal(b$benchmarks, b_ref$benchmarks)
  # asymp: alt_sample_size schaalt met de N van het model (75), niet met 1000
  ra <- with_msgs(benchmark(g_n, model_type = "asymp", iter = 5, seed = 1, alt_sample_size = 150))
  expect_true(any(grepl("'sample_nobs' of the GORIC\\(A\\) object \\(1000\\) differs", ra$msgs)))
  expect_equal(ra$value$pop_VCOV, g_n$VCOV * 75 / 150)
  expect_equal(unname(ra$value$sample_size), 150)
  # zonder afwijking: geen melding
  r0 <- with_msgs(benchmark(g_compl, model_type = "asymp", iter = 3, seed = 1))
  expect_false(any(grepl("differs from the sample size of the fitted model", r0$msgs)))
})


# -----------------------------------------------------------------------------
# W3-N-02: de groeperingsterm is de term waar de hypothesen over gaan
# -----------------------------------------------------------------------------

test_that("W3-N-02: de factorterm van de hypothesen wordt de groeperingsterm; anders duidelijke fout", {
  set.seed(11)
  df_bm$sex <- factor(sample(c("f", "m"), nrow(df_bm), TRUE))
  f3 <- lm(y ~ -1 + sex + group, data = df_bm)
  # hypothese over sex (de celgemiddelden): sex is de groeperingsterm
  g_sex <- goric(f3, hypotheses = list(H1 = "sexf < sexm"), comparison = "complement")
  r <- with_msgs(benchmark(g_sex, model_type = "means", iter = 5, seed = 1))
  expect_true(any(grepl("term 'sex' \\(sexf, sexm\\) are treated as the group means", r$msgs)))
  expect_equal(r$value$ngroups, 2)
  expect_equal(unname(r$value$group_size), unname(c(table(df_bm$sex))))
  s2 <- sum(residuals(f3)^2) / nrow(df_bm)
  expect_equal(r$value$cohens_f_observed,
               cohens_f_ref(coef(f3)[1:2], c(table(df_bm$sex)), s2))
  expect_equal(unname(r$value$pop_group_means[, "group2"]), rep(unname(coef(f3)["group2"]), 2))
  # hypothese over group (hier contrasten, group is niet de eerste factorterm): fout
  g_grp <- goric(f3, hypotheses = list(H1 = "group2 < group3"), comparison = "complement")
  expect_error(quiet(benchmark(g_grp, model_type = "means", iter = 5, seed = 1)),
               "single factor term.*group2, group3.*term 'group'.*contrasts.*model_type = 'asymp'")
  # ook met group_size: fout (de coefficienten zijn geen gemiddelden)
  expect_error(quiet(benchmark(g_grp, model_type = "means", iter = 5, seed = 1,
                               group_size = c(20, 25, 30))), "contrasts")
  # group eerst: group is de groeperingsterm; hypothese over sexm (contrast): fout
  f2 <- lm(y ~ -1 + group + sex, data = df_bm)
  g_sexm <- goric(f2, hypotheses = list(H1 = "sexm < 0"), comparison = "complement")
  expect_error(quiet(benchmark(g_sexm, model_type = "means", iter = 5, seed = 1)),
               "single factor term.*sexm.*term 'sex'.*contrasts")
  # hypothesen over meerdere termen: fout
  g_mix <- goric(f2, hypotheses = list(H1 = "group1 < group2; sexm < 0"), comparison = "complement")
  expect_error(quiet(benchmark(g_mix, model_type = "means", iter = 5, seed = 1)),
               "group1, group2, sexm.*do not belong to a single factor")
  # hypothese over een covariaat: fout
  fx <- lm(y ~ -1 + x + group, data = df_bm)
  g_x <- goric(fx, hypotheses = list(H1 = "x > 0"), comparison = "complement")
  expect_error(quiet(benchmark(g_x, model_type = "means", iter = 5, seed = 1)),
               "coefficient\\(s\\) x, which do not belong")
  # hypothese over group bij x eerst: group blijft de groeperingsterm
  g_xg <- goric(fx, hypotheses = list(H1 = h1_bm), comparison = "complement")
  r2 <- with_msgs(benchmark(g_xg, model_type = "means", iter = 5, seed = 1))
  expect_true(any(grepl("term 'group' \\(group1, group2, group3\\) are treated", r2$msgs)))
  expect_equal(unname(r2$value$group_size), n_g)
  # ':='-definities: de verwezen coefficienten worden uit de constraintmatrix gehaald
  g_def <- goric(fit_bm, hypotheses = list(H1 = "d := group2 - group1; d > 0"),
                 comparison = "complement")
  expect_equal(quiet(benchmark(g_def, model_type = "means", iter = 3, seed = 1))$ngroups, 3)
})


# -----------------------------------------------------------------------------
# W3-N-03: plot: niet-eindige draws per populatie x vergelijking
# -----------------------------------------------------------------------------

test_that("W3-N-03: plot() meldt per populatie x vergelijking welke niet getoond wordt", {
  est <- c(a = 0, b = 60)
  V <- diag(1e-4, 2)
  dimnames(V) <- list(names(est), names(est))
  g3 <- goric(est, VCOV = V, hypotheses = list(H1 = "a < b", H2 = "a > b", H3 = "a < b + 100"),
              comparison = "none")
  b3 <- quiet(benchmark(g3, iter = 10, seed = 1, pop_est = rbind(Big = est, Null = c(30, 30))))
  # Big: 0 eindige draws voor vs. H2, 10 voor vs. H3
  expect_equal(sum(is.finite(b3$combined_values$rgw_combined[["pop_est = Big"]][, "vs. H2"])), 0)
  expect_equal(sum(is.finite(b3$combined_values$rgw_combined[["pop_est = Big"]][, "vs. H3"])), 10)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  r <- with_msgs(plot(b3, output_type = "rgw"))
  expect_true(any(grepl("Population estimates = Big \\(H1 vs. H2\\)", r$msgs)))
  expect_false(any(grepl("Big \\(H1 vs. H2; H1 vs. H3\\)|H1 vs. H3\\)", r$msgs)))
  expect_no_error(print(r$value))
  # paneel vs. H2: percentiellijnen van Null; paneel vs. H3: van Big (de eerste populatie)
  labs2 <- levels(r$value$plots[[1]]$layers[[2]]$data$percentile_label)
  labs3 <- levels(r$value$plots[[2]]$layers[[2]]$data$percentile_label)
  expect_true(any(grepl("Percentile \\(Population estimates = Null\\)", labs2)))
  expect_true(any(grepl("Percentile \\(Population estimates = Big\\)", labs3)))
  # populatie zonder enige eindige draw: zonder vergelijking tussen haakjes
  g1 <- goric(est, VCOV = V, hypotheses = list(H1 = "a < b"), comparison = "complement")
  b1 <- quiet(benchmark(g1, iter = 5, seed = 1, pop_est = rbind(Big = est, Null = c(30, 30))))
  r1 <- with_msgs(plot(b1, output_type = "rgw"))
  expect_true(any(grepl("plot \\(if only for some comparisons, these are given in parentheses\\): Population estimates = Big\\. See", r1$msgs)))
})


# -----------------------------------------------------------------------------
# W3-N-04: gewogen lm-fit wordt geweigerd in benchmark_means
# -----------------------------------------------------------------------------

test_that("W3-N-04: een gewogen lm-fit geeft een duidelijke fout bij model_type = 'means'", {
  set.seed(4)
  w <- runif(nrow(df_bm), 0.5, 2)
  fw <- lm(y ~ -1 + group, data = df_bm, weights = w)
  gw <- goric(fw, hypotheses = list(H1 = h1_bm), comparison = "complement")
  expect_error(quiet(benchmark(gw, model_type = "means", iter = 5, seed = 1)),
               "weighted fit.*model_type = 'asymp'")
  # de asymptotische route werkt wel
  bw <- quiet(benchmark(gw, model_type = "asymp", iter = 3, seed = 1))
  expect_s3_class(bw, "benchmark_asymp")
  # gewichten die allemaal 1 zijn: gewoon toegestaan (sigma2 = RSS/N)
  fw1 <- lm(y ~ -1 + group, data = df_bm, weights = rep(1, nrow(df_bm)))
  gw1 <- goric(fw1, hypotheses = list(H1 = h1_bm), comparison = "complement")
  b1 <- quiet(benchmark(gw1, model_type = "means", iter = 3, seed = 1))
  expect_equal(b1$res_var, sum(residuals(fit_bm)^2) / nrow(df_bm))
})

test_that("W5-N-01: benchmark() geeft een duidelijke fout bij hypothesen als constraintmatrix", {
  set.seed(5)
  d <- data.frame(g = factor(rep(1:3, each = 10))); d$y <- rnorm(30) + as.numeric(d$g)
  fit <- lm(y ~ -1 + g, data = d)
  A <- rbind(c(-1, 1, 0), c(0, -1, 1))
  g <- suppressMessages(goric(fit, hypotheses = list(H1 = list(constraints = A, rhs = c(0, 0), neq = 0))))
  expect_error(benchmark(g, iter = 5, seed = 1), "specified as text")
  expect_error(benchmark(g, model_type = "means", iter = 5, seed = 1), "specified as text")
})

# -----------------------------------------------------------------------------
# FX6-B: de foutkans (error_prob_pref_hypo) gebruikt de penalty_factor van het object
# -----------------------------------------------------------------------------

test_that("FX6-B: error_prob_pref_hypo gebruikt de penalty_factor van het goric-object", {
  est <- c(x = 1, y = 2, z = 3)
  hyp <- list(H1 = "x < y < z", H2 = "x > y > z")
  # asymptotische route, penalty_factor = 6: gelijk aan het complement-gewicht
  # van de voorkeurshypothese met dezelfde penalty_factor (onafhankelijk berekend)
  g6 <- goric(est, VCOV = diag(3), hypotheses = hyp, type = "gorica", penalty_factor = 6)
  b6 <- quiet(benchmark(g6, iter = 5, seed = 1))
  gc6 <- goric(est, VCOV = diag(3), hypotheses = hyp[1], type = "gorica",
               penalty_factor = 6, comparison = "complement")
  expect_equal(b6$error_prob_pref_hypo, gc6$result$gorica.weights[2], tolerance = 1e-8)
  expect_equal(b6$error_prob_pref_hypo, 0.06008665, tolerance = 1e-6)
  # standaard penalty_factor: ongewijzigd (zelfde waarde als voorheen)
  g2 <- goric(est, VCOV = diag(3), hypotheses = hyp, type = "gorica")
  b2 <- quiet(benchmark(g2, iter = 5, seed = 1))
  gc2 <- goric(est, VCOV = diag(3), hypotheses = hyp[1], type = "gorica", comparison = "complement")
  expect_equal(b2$error_prob_pref_hypo, gc2$result$gorica.weights[2], tolerance = 1e-8)
  expect_equal(b2$error_prob_pref_hypo, 0.2528757, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(b6$error_prob_pref_hypo, b2$error_prob_pref_hypo)))
  # route met groepsgemiddelden (lm-fit; benchmark herberekent met gorica)
  hyp_m <- list(H1 = h1_bm, H2 = "group1 > group2 > group3")
  gm4 <- goric(fit_bm, hypotheses = hyp_m, penalty_factor = 4)
  bm4 <- quiet(benchmark(gm4, iter = 5, seed = 1))
  gmc4 <- goric(fit_bm, hypotheses = hyp_m[1], type = "gorica", penalty_factor = 4,
                comparison = "complement")
  expect_equal(bm4$error_prob_pref_hypo, gmc4$result$gorica.weights[2], tolerance = 1e-8)
  gm2 <- goric(fit_bm, hypotheses = hyp_m)
  bm2 <- quiet(benchmark(gm2, iter = 5, seed = 1))
  gmc2 <- goric(fit_bm, hypotheses = hyp_m[1], type = "gorica", comparison = "complement")
  expect_equal(bm2$error_prob_pref_hypo, gmc2$result$gorica.weights[2], tolerance = 1e-8)
  # object met comparison = 'complement': 1 - gewicht van de voorkeurshypothese
  gcm <- goric(est, VCOV = diag(3), hypotheses = hyp[1], type = "gorica",
               penalty_factor = 6, comparison = "complement")
  bcm <- quiet(benchmark(gcm, iter = 5, seed = 1))
  expect_equal(bcm$error_prob_pref_hypo, 1 - gcm$result$gorica.weights[1], tolerance = 1e-8)
  expect_equal(bcm$error_prob_pref_hypo, gc6$result$gorica.weights[2], tolerance = 1e-8)
})
