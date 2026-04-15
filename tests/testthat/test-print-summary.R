# =============================================================================
# Tests: print() en summary() methoden voor alle klassen
# =============================================================================
# Zorg dat vereiste imports beschikbaar zijn in het huidige library-pad
# for (pkg in c("lavaan", "mvtnorm", "tmvtnorm", "quadprog", "norm")) {
#   if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# }

# # Laad het restriktor package vanuit broncode
# setwd("/Workspace/Users/l.vanbrabant@ggdwestbrabant.nl/restriktor")
# devtools::load_all(".")

# --- Gedeelde testdata ---

set.seed(123)
n <- 60
group <- factor(rep(1:3, each = n / 3))
y <- c(rnorm(n / 3, mean = 5), rnorm(n / 3, mean = 3), rnorm(n / 3, mean = 1))
df <- data.frame(y = y, group = group)
fit_lm <- lm(y ~ -1 + group, data = df)

est_3 <- c(x1 = 5.0, x2 = 3.0, x3 = 1.0)
VCOV_3 <- diag(c(0.04, 0.04, 0.04))


# =============================================================================
# print.con_goric / summary.con_goric
# =============================================================================

test_that("print.con_goric: werkt voor gorica", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_output(print(result), "restriktor")
})

test_that("print.con_goric: werkt voor goric met lm", {
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = "group1 > group2 > group3"))
  expect_output(print(result), "restriktor")
})

test_that("summary.con_goric: werkt voor gorica", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_output(summary(result), "restriktor")
})

test_that("summary.con_goric: werkt voor goric met meerdere hypothesen", {
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = "group1 > group2 > group3",
                                    H2 = "group1 > group2 = group3"))
  expect_output(summary(result), "restriktor")
})

test_that("summary.con_goric: toont ratio-tabel", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3",
                                    H2 = "x3 > x2 > x1"))
  # Summary moet ratio informatie bevatten
  out <- capture.output(summary(result))
  expect_true(any(grepl("ratio|weight|Ratio", out)))
})


# =============================================================================
# coef.con_goric
# =============================================================================

test_that("coef.con_goric: retourneert matrix met geschatte coëfficiënten", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  coefs <- coef(result)
  expect_true(is.matrix(coefs) || is.data.frame(coefs))
  expect_equal(nrow(coefs), nrow(result$result))
})


# =============================================================================
# print.restriktor / summary.restriktor
# =============================================================================

test_that("print.restriktor: werkt voor conLM", {
  fit_con <- restriktor(fit_lm, constraints = "group1 > group2 > group3")
  expect_output(print(fit_con), "restriktor|Constraints")
})

test_that("summary.restriktor: toont coëfficiënten", {
  fit_con <- restriktor(fit_lm, constraints = "group1 > group2 > group3")
  out <- capture.output(summary(fit_con))
  expect_true(any(grepl("group1|group2|group3|Coefficients|coef", out)))
})


# =============================================================================
# coef.restriktor / logLik.restriktor / model.matrix.restriktor
# =============================================================================

test_that("coef.restriktor: retourneert benoemde vector", {
  fit_con <- restriktor(fit_lm, constraints = "group1 > group2 > group3")
  b <- coef(fit_con)
  expect_true(is.numeric(b))
  expect_true(!is.null(names(b)))
  expect_true("group1" %in% names(b))
})

test_that("logLik.restriktor: retourneert numerieke loglikelihood", {
  fit_con <- restriktor(fit_lm, constraints = "group1 > group2 > group3")
  ll <- logLik(fit_con)
  expect_true(is.numeric(ll))
  expect_true(is.finite(ll))
})

test_that("model.matrix.restriktor: retourneert designmatrix", {
  fit_con <- restriktor(fit_lm, constraints = "group1 > group2 > group3")
  mm <- model.matrix(fit_con)
  expect_true(is.matrix(mm))
  expect_equal(nrow(mm), nrow(df))
})


# =============================================================================
# print.conTest
# =============================================================================

test_that("print.conTest: werkt voor conTestF", {
  fit_con <- restriktor(fit_lm, constraints = "group1 > group2 > group3")
  ct <- conTest(fit_con)
  expect_output(print(ct), "Type A|Type B|F-bar|p-value")
})


# =============================================================================
# print.evSyn / summary.evSyn
# =============================================================================

test_that("print.evSyn: werkt", {
  est1 <- c(x1 = 5, x2 = 3)
  est2 <- c(x1 = 4, x2 = 2)
  VCOV1 <- diag(c(0.1, 0.1))
  VCOV2 <- diag(c(0.15, 0.15))

  g1 <- goric(est1, VCOV = VCOV1, type = "gorica",
               hypotheses = list(H1 = "x1 > x2"))
  g2 <- goric(est2, VCOV = VCOV2, type = "gorica",
               hypotheses = list(H1 = "x1 > x2"))

  es <- evSyn(object = list(g1, g2), type = "added")
  expect_output(print(es), "evidence|Evidence|synthesis|Synthesis")
})


# =============================================================================
# print.benchmark
# =============================================================================

test_that("print.benchmark: werkt voor benchmark object", {
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = "group1 > group2 > group3"))
  bm <- benchmark(result, model_type = "means", iter = 10)
  expect_output(print(bm), "benchmark|Benchmark|ratio|pop")
})


# =============================================================================
# print.goric_ICw
# =============================================================================

test_that("print.goric_ICw: werkt", {
  ic_vals <- c(model1 = 10, model2 = 12)
  result <- calculate_IC_weights(ic_vals)
  expect_output(print(result), "weight|Weight")
})

