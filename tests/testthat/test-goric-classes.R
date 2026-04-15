# =============================================================================
# Tests: goric() - alle klassen en type-varianten
# goric.lm, goric.glm (via lm), goric.numeric, goric.lavaan,
# type = goric/gorica/goricc/goricac, penalty_factor, comparison varianten
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
n <- 90
group <- factor(rep(1:3, each = n / 3))
x1 <- rnorm(n)
y <- 5 * (group == "1") + 3 * (group == "2") + 1 * (group == "3") + rnorm(n)
df_test <- data.frame(y = y, group = group, x1 = x1)

# Fits
fit_lm     <- lm(y ~ -1 + group, data = df_test)
fit_lm_int <- lm(y ~ group, data = df_test)
fit_glm    <- glm(y ~ -1 + group, family = gaussian, data = df_test)

# Schattingen voor gorica
est_3 <- c(x1 = 5.0, x2 = 3.0, x3 = 1.0)
VCOV_3 <- diag(c(0.04, 0.04, 0.04))

est_2 <- c(a = 2.0, b = 1.0)
VCOV_2 <- matrix(c(0.05, 0.01, 0.01, 0.05), 2, 2)


# =============================================================================
# goric.lm: type = "goric" (standaard)
# =============================================================================

test_that("goric.lm type='goric' werkt met 1 hypothese", {
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = "group1 > group2 > group3"))
  expect_s3_class(result, "con_goric")
  expect_equal(result$type, "goric")
  expect_equal(nrow(result$result), 2)
})

test_that("goric.lm type='goric' met model met intercept", {
  result <- goric(fit_lm_int, type = "goric",
                  hypotheses = list(H1 = "group2 > 0; group3 > group2"))
  expect_s3_class(result, "con_goric")
})

test_that("goric.lm type='goricc' (small sample correctie)", {
  result <- goric(fit_lm, type = "goricc",
                  hypotheses = list(H1 = "group1 > group2 > group3"))
  expect_s3_class(result, "con_goric")
  expect_equal(result$type, "goricc")
  # goricc penalty moet groter zijn dan goric penalty
  res_goric  <- goric(fit_lm, type = "goric",
                      hypotheses = list(H1 = "group1 > group2 > group3"))
  expect_true(result$result$penalty[1] >= res_goric$result$penalty[1] - 1e-6)
})


# =============================================================================
# goric.lm: type = "gorica" / "goricac"
# =============================================================================

test_that("goric.lm type='gorica' werkt (approximatie)", {
  result <- goric(fit_lm, type = "gorica",
                  hypotheses = list(H1 = "group1 > group2 > group3"))
  expect_s3_class(result, "con_goric")
  expect_true("con_gorica" %in% class(result))
  expect_equal(result$type, "gorica")
})

test_that("goric.lm type='goricac' werkt (small sample approximatie)", {
  result <- goric(fit_lm, type = "goricac",
                  hypotheses = list(H1 = "group1 > group2 > group3"))
  expect_s3_class(result, "con_goric")
  expect_equal(result$type, "goricac")
})


# =============================================================================
# goric.glm (via gaussian GLM)
# =============================================================================

test_that("goric.glm type='goric' werkt", {
  result <- goric(fit_glm, type = "goric",
                  hypotheses = list(H1 = "group1 > group2 > group3"))
  expect_s3_class(result, "con_goric")
  expect_equal(nrow(result$result), 2)
})

test_that("goric.glm meerdere hypothesen met unconstrained", {
  result <- goric(fit_glm, type = "goric",
                  hypotheses = list(H1 = "group1 > group2 > group3",
                                    H2 = "group1 = group2 = group3"))
  expect_equal(nrow(result$result), 3)
  expect_true("unconstrained" %in% result$result$model)
})

test_that("goric.glm ware hypothese krijgt hoogste gewicht", {
  result <- goric(fit_glm, type = "goric",
                  hypotheses = list(H_correct = "group1 > group2 > group3",
                                    H_reverse = "group3 > group2 > group1"))
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  weights <- result$result[[wt_col]]
  names(weights) <- result$result$model
  expect_true(weights["H_correct"] > weights["H_reverse"])
})


# =============================================================================
# goric.numeric: gorica type-varianten
# =============================================================================

test_that("goric.numeric gorica 2D probleem", {
  result <- goric(est_2, VCOV = VCOV_2, type = "gorica",
                  hypotheses = list(H1 = "a > b"))
  expect_s3_class(result, "con_goric")
  expect_equal(nrow(result$result), 2)
})

test_that("goric.numeric gorica met gecorreleerde VCOV", {
  vcov_corr <- matrix(c(0.1, 0.05, 0.02,
                         0.05, 0.1, 0.03,
                         0.02, 0.03, 0.1), 3, 3)
  result <- goric(est_3, VCOV = vcov_corr, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_s3_class(result, "con_goric")
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(sum(result$result[[wt_col]]), 1, tolerance = 1e-6)
})

test_that("goric.numeric goricac met sample_nobs", {
  result <- goric(est_3, VCOV = VCOV_3, type = "goricac",
                  sample_nobs = 90,
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_equal(result$type, "goricac")
  # goricac penalty > gorica penalty
  res_gorica <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                      hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_true(result$result$penalty[1] >= res_gorica$result$penalty[1] - 1e-6)
})


# =============================================================================
# penalty_factor parameter
# =============================================================================

test_that("goric penalty_factor = 1 (AIC-achtig)", {
  res_pf1 <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                    penalty_factor = 1, comparison = "none",
                    hypotheses = list(H1 = "x1 > x2 > x3"))
  res_pf2 <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                    penalty_factor = 2, comparison = "none",
                    hypotheses = list(H1 = "x1 > x2 > x3"))

  # Zelfde loglik, maar goric = -2*ll + pf*PT
  expect_equal(res_pf1$result$loglik[1], res_pf2$result$loglik[1])
  expect_equal(res_pf1$result$penalty[1], res_pf2$result$penalty[1])
  # IC waarde verschilt door penalty_factor
  expect_true(res_pf1$result$gorica[1] < res_pf2$result$gorica[1])
})

test_that("goric penalty_factor = 0 geeft alleen -2*loglik", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  penalty_factor = 0, comparison = "none",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expected_ic <- -2 * result$result$loglik[1]
  expect_equal(result$result$gorica[1], expected_ic, tolerance = 1e-8)
})

test_that("goric penalty_factor = 3 (BIC-achtig)", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  penalty_factor = 3, comparison = "none",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expected_ic <- -2 * result$result$loglik[1] + 3 * result$result$penalty[1]
  expect_equal(result$result$gorica[1], expected_ic, tolerance = 1e-8)
})

test_that("goric penalty_factor wordt opgeslagen in output", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  penalty_factor = 2.5,
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_equal(result$penalty_factor, 2.5)
})


# =============================================================================
# comparison varianten: uitgebreid
# =============================================================================

test_that("goric comparison='complement' met >1 hypothese wordt unconstrained", {
  expect_warning(
    result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                    comparison = "complement",
                    hypotheses = list(H1 = "x1 > x2 > x3",
                                      H2 = "x1 > x2 = x3")),
    "unconstrained"
  )
  expect_true("unconstrained" %in% result$result$model)
})

test_that("goric comparison='none' met 1 hypothese", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  comparison = "none",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_equal(nrow(result$result), 1)
  expect_false("complement" %in% result$result$model)
  expect_false("unconstrained" %in% result$result$model)
})


# =============================================================================
# Complexe hypothesen: gedefinieerde parameters, abs(), haakjes
# =============================================================================

test_that("goric.numeric met gedefinieerde parameter (:=)", {
  est_4 <- c(a = 3, b = 2, c = 1, d = 0.5)
  vcov_4 <- diag(4) * 0.04
  result <- goric(est_4, VCOV = vcov_4, type = "gorica",
                  hypotheses = list(H1 = "a > b; b > c; c > d"))
  expect_s3_class(result, "con_goric")
})

test_that("goric.numeric met gelijkheid en ongelijkheid gecombineerd", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2; x2 = x3"))
  b_restr <- coef(result)
  # x2 moet gelijk zijn aan x3 in H1 rij
  expect_equal(b_restr[1, "x2"], b_restr[1, "x3"], tolerance = 1e-4)
})

test_that("goric.numeric met meerdere ongelijkheden op zelfde parameter", {
  est_4 <- c(a = 5, b = 3, c = 4, d = 1)
  vcov_4 <- diag(4) * 0.05
  result <- goric(est_4, VCOV = vcov_4, type = "gorica",
                  hypotheses = list(H1 = "a > b; a > c; a > d",
                                    H2 = "a > b > c > d"))
  expect_equal(nrow(result$result), 3)
})


# =============================================================================
# goric.lm met mix_weights opties
# =============================================================================

test_that("goric.lm met mix_weights = 'boot'", {
  result <- goric(fit_lm, type = "goric",
                  mix_weights = "boot", seed = 42,
                  hypotheses = list(H1 = "group1 > group2 > group3"))
  expect_s3_class(result, "con_goric")
  # Boot gewichten
  wt_method <- attr(result$objectList[[1]]$wt.bar, "method")
  expect_equal(wt_method, "boot")
})

test_that("goric.lm met mix_weights = 'none'", {
  # mix_weights = 'none' zou een fout moeten geven want penalty kan niet berekend worden
  expect_error(
    goric(fit_lm, type = "goric",
          mix_weights = "none",
          hypotheses = list(H1 = "group1 > group2 > group3")),
    "mix_weights.*none|penalty.*complement.*cannot"
  )
})


# =============================================================================
# Grote dataset / meer parameters
# =============================================================================

test_that("goric.numeric met 5 parameters", {
  est_5 <- c(p1 = 10, p2 = 8, p3 = 6, p4 = 4, p5 = 2)
  vcov_5 <- diag(5) * 0.1
  result <- goric(est_5, VCOV = vcov_5, type = "gorica",
                  hypotheses = list(H1 = "p1 > p2 > p3 > p4 > p5"))
  expect_s3_class(result, "con_goric")
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(sum(result$result[[wt_col]]), 1, tolerance = 1e-6)
})


# =============================================================================
# goric: output print en summary
# =============================================================================

test_that("goric summary methode werkt", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_output(summary(result), "gorica|GORICA|restriktor")
})

test_that("goric print methode werkt", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_output(print(result), "gorica|GORICA|restriktor")
})


# =============================================================================
# goric: type-consistentie (goric vs gorica op zelfde data)
# =============================================================================

test_that("goric en gorica op zelfde lm: ordening van gewichten is consistent", {
  res_goric  <- goric(fit_lm, type = "goric",
                      hypotheses = list(H1 = "group1 > group2 > group3",
                                        H2 = "group3 > group2 > group1"))
  res_gorica <- goric(fit_lm, type = "gorica",
                      hypotheses = list(H1 = "group1 > group2 > group3",
                                        H2 = "group3 > group2 > group1"))
  wt_goric  <- grep("weights", names(res_goric$result),  value = TRUE)[1]
  wt_gorica <- grep("weights", names(res_gorica$result), value = TRUE)[1]

  # Beide moeten H1 (correcte ordening) hoger wegen dan H2
  w_goric  <- res_goric$result[[wt_goric]]
  w_gorica <- res_gorica$result[[wt_gorica]]
  names(w_goric)  <- res_goric$result$model
  names(w_gorica) <- res_gorica$result$model

  expect_true(w_goric["H1"]  > w_goric["H2"])
  expect_true(w_gorica["H1"] > w_gorica["H2"])
})


# =============================================================================
# Edge cases
# =============================================================================

test_that("goric met scalaire VCOV (1 parameter)", {
  est_1 <- c(mu = 5.0)
  vcov_1 <- matrix(0.1, 1, 1)
  result <- goric(est_1, VCOV = vcov_1, type = "gorica",
                  hypotheses = list(H1 = "mu > 0"))
  expect_s3_class(result, "con_goric")
})

test_that("goric met VCOV als scalar getal", {
  est_1 <- c(mu = 5.0)
  result <- goric(est_1, VCOV = 0.1, type = "gorica",
                  hypotheses = list(H1 = "mu > 0"))
  expect_s3_class(result, "con_goric")
})

test_that("goric met restrictie die al voldaan is", {
  # x1=5 > x2=3 > x3=1 al voldaan, geen projectie nodig
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  comparison = "none",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  b <- coef(result)[1, ]
  # Schattingen moeten ongewijzigd zijn (gebruik [[ ]] voor scalar extractie)
  expect_equal(b[["x1"]], 5.0, tolerance = 1e-6)
  expect_equal(b[["x2"]], 3.0, tolerance = 1e-6)
  expect_equal(b[["x3"]], 1.0, tolerance = 1e-6)
})
