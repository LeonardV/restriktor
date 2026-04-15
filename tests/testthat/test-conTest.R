# =============================================================================
# Tests: conTest(), conTestF, conTestLRT, conTestScore, conTest_ceq
#
# Snelheidsoptimalisatie:
#   - restriktor() met mix_weights = "pmvnorm" zodat wt.bar 1x berekend en
#     gecachet wordt in het object. Alle conTest() calls hergebruiken dit.
#   - Maximaal 1 conTest-summary call per (object, test)-combinatie.
#   - Alle assertions op gecachte resultaten.
# =============================================================================

# # Zorg dat vereiste imports beschikbaar zijn in het huidige library-pad
# for (pkg in c("lavaan", "mvtnorm", "tmvtnorm", "quadprog", "norm")) {
#   if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# }

# # Laad het restriktor package vanuit broncode
# setwd("/Workspace/Users/l.vanbrabant@ggdwestbrabant.nl/restriktor")
# devtools::load_all(".")

# --- Gedeelde testdata ---
set.seed(42)
n <- 60
group <- factor(rep(1:3, each = n / 3))
x1 <- rnorm(n)
y <- 2 + 1.5 * (group == "2") + 3 * (group == "3") + rnorm(n)
df_test <- data.frame(y = y, group = group, x1 = x1)

fit_lm  <- lm(y ~ -1 + group, data = df_test)
fit_glm <- glm(y ~ -1 + group, family = gaussian, data = df_test)
fit_rlm <- MASS::rlm(y ~ -1 + group, data = df_test, method = "MM")

# mix_weights = "pmvnorm" berekent wt.bar 1x en cachet het in het object.
# Alle volgende conTest() calls hergebruiken de gecachte gewichten.
con_lm  <- restriktor(fit_lm,  constraints = "group1 < group2 < group3")
con_glm <- restriktor(fit_glm, constraints = "group1 < group2 < group3")
con_rlm <- restriktor(fit_rlm, constraints = "group1 < group2 < group3")
con_lm_ceq <- restriktor(fit_lm, constraints = "group1 = group2; group2 = group3")


# =============================================================================
# conTest op conLM: F-toets (1 summary call, alle typen)
# =============================================================================

ct_lm_F <- conTest(con_lm, type = "summary", test = "F")

test_that("conTest conLM F-toets: summary, type A, B, global", {
  expect_s3_class(ct_lm_F, "conTest")
  # Type A
  expect_true(ct_lm_F$A$Ts >= -1e-10)
  expect_true(ct_lm_F$A$pvalue >= 0 & ct_lm_F$A$pvalue <= 1)
  # Type B
  expect_true(ct_lm_F$B$pvalue >= 0 & ct_lm_F$B$pvalue <= 1)
  # Global
  expect_true(ct_lm_F$global$pvalue >= 0 & ct_lm_F$global$pvalue <= 1)
})


# =============================================================================
# conTest op conLM: LRT (1 summary call)
# =============================================================================

ct_lm_LRT <- conTest(con_lm, type = "summary", test = "LRT")

test_that("conTest conLM LRT: type A, B, global", {
  expect_true(ct_lm_LRT$A$Ts >= -1e-10)
  expect_true(ct_lm_LRT$A$pvalue >= 0 & ct_lm_LRT$A$pvalue <= 1)
  expect_true(ct_lm_LRT$B$pvalue >= 0 & ct_lm_LRT$B$pvalue <= 1)
  expect_true(!is.null(ct_lm_LRT$global$pvalue))
})


# =============================================================================
# conTest op conLM: Score-toets (1 call)
# =============================================================================

ct_lm_score <- conTest(con_lm, type = "A", test = "score")

test_that("conTest conLM Score: type A", {
  expect_true(ct_lm_score$Ts >= -1e-10)
  expect_true(ct_lm_score$pvalue >= 0 & ct_lm_score$pvalue <= 1)
})


# =============================================================================
# conTest op conGLM (1 summary call)
# =============================================================================

ct_glm_F <- conTest(con_glm, type = "summary", test = "F")

test_that("conTest conGLM F-toets: summary", {
  expect_s3_class(ct_glm_F, "conTest")
  expect_true(ct_glm_F$A$pvalue >= 0 & ct_glm_F$A$pvalue <= 1)
})


# =============================================================================
# conTest op conRLM (1 summary call)
# =============================================================================

ct_rlm_F <- conTest(con_rlm, type = "summary", test = "F")

test_that("conTest conRLM F-toets", {
  expect_true(ct_rlm_F$A$pvalue >= 0 & ct_rlm_F$A$pvalue <= 1)
})


# =============================================================================
# conTest_ceq: gelijkheidsbeperkingen
# =============================================================================

test_that("conTest_ceq conLM geeft testresultaat", {
  result <- conTest_ceq(con_lm_ceq)
  expect_true(!is.null(result))
})


# =============================================================================
# iht() alias (hergebruik gecachte gewichten)
# =============================================================================

test_that("iht() is alias voor conTest()", {
  r2 <- iht(con_lm, type = "A", test = "F")
  expect_equal(ct_lm_F$A$Ts, r2$Ts, tolerance = 1e-10)
  expect_equal(ct_lm_F$A$pvalue, r2$pvalue, tolerance = 1e-10)
})


# =============================================================================
# Wiskundige consistentie (hergebruik)
# =============================================================================

test_that("conTest: teststatistiek niet-negatief, p-waarden in [0,1]", {
  expect_true(ct_lm_F$A$Ts >= -1e-10)
  expect_true(ct_lm_LRT$A$Ts >= -1e-10)
  expect_true(ct_lm_score$Ts >= -1e-10)
  expect_true(ct_lm_F$A$pvalue >= 0 & ct_lm_F$A$pvalue <= 1)
  expect_true(ct_lm_F$B$pvalue >= 0 & ct_lm_F$B$pvalue <= 1)
})


# =============================================================================
# Foutafhandeling
# =============================================================================

test_that("conTest fout bij object zonder restricties", {
  expect_error(conTest(fit_lm), "restriktor")
})

test_that("conTest fout bij onbekende test methode", {
  expect_error(conTest(con_lm, type = "A", test = "onbekend"), "unknown")
})

test_that("conTest fout bij ongeldige type", {
  expect_error(conTest(con_lm, type = "xyz", test = "F"))
})

test_that("conTest conRLM fout bij LRT (niet beschikbaar)", {
  expect_error(conTest(con_rlm, type = "A", test = "LRT"), "unknown")
})

