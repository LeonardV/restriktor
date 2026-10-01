# =============================================================================
# Tests: goric() voor multivariate lm-objecten (mlm)
# =============================================================================

set.seed(1)
df_mlm <- data.frame(Group = factor(rep(c("Active", "No", "Passive"), c(6, 6, 5)),
                                    levels = c("Active", "No", "Passive")))
df_mlm$Age  <- 30 + 2 * (df_mlm$Group == "No") + rnorm(17, sd = 3)
df_mlm$Age2 <- df_mlm$Age^2 / 30 + rnorm(17)
fit_mv <- lm(cbind(Age, Age2) ~ Group, data = df_mlm)

H_mv <- "Age.GroupPassive = 0; Age.GroupPassive < Age.GroupNo"

# handmatig berekende multivariaat-normale loglik van het ongerestricteerde model
E_mv <- residuals(fit_mv)
ll_mv_ref <- sum(mvtnorm::dmvnorm(E_mv, sigma = crossprod(E_mv) / nrow(E_mv), log = TRUE))

run_goric_mv <- function(...) suppressMessages(goric(fit_mv, ...))

test_that("goric mlm: alle combinaties van comparison en type draaien", {
  for (cmp in c("complement", "unconstrained", "none")) {
    for (ty in c("goric", "goricc", "gorica")) {
      res <- run_goric_mv(hypotheses = list(H1 = H_mv), comparison = cmp, type = ty)
      expect_s3_class(res, "con_goric")
      expect_true(all(is.finite(res$result[[ty]])))
    }
  }
})

test_that("goric mlm: loglik van complement en unconstrained is de MVN-loglik", {
  res_c <- run_goric_mv(hypotheses = list(H1 = H_mv), comparison = "complement")
  expect_equal(res_c$result$loglik[2], ll_mv_ref, tolerance = 1e-8)
  res_u <- run_goric_mv(hypotheses = list(H1 = H_mv), comparison = "unconstrained")
  expect_equal(res_u$result$loglik[2], ll_mv_ref, tolerance = 1e-8)
})

test_that("goric mlm: sample size is N en niet N x aantal responsen", {
  res <- run_goric_mv(hypotheses = list(H1 = H_mv), comparison = "unconstrained",
                      type = "goricc")
  N <- nrow(df_mlm); p <- length(coef(fit_mv))
  expect_equal(res$result$penalty[2], N * (p + 1) / (N - p - 2))
  expect_equal(res$sample_nobs, N)
})

test_that("goric mlm: complement van alleen ongelijkheden (in lijn met de data)", {
  b <- coef(fit_mv)
  # kies de richting zodat de ongerestricteerde schattingen aan H voldoen
  H_in <- if (b["GroupNo", "Age"] > b["GroupPassive", "Age"]) {
    "Age.GroupNo > Age.GroupPassive"
  } else {
    "Age.GroupNo < Age.GroupPassive"
  }
  res <- run_goric_mv(hypotheses = list(H1 = H_in), comparison = "complement")
  expect_equal(nrow(res$result), 2L)
  expect_true(res$result$loglik[2] <= ll_mv_ref + 1e-8)
  expect_true(all(is.finite(unlist(res$ormle$b.restr))))
})

test_that("goric mlm: Heq = TRUE werkt", {
  res <- run_goric_mv(hypotheses = list(H1 = "Age.GroupPassive < Age.GroupNo"),
                      comparison = "complement", Heq = TRUE)
  expect_equal(nrow(res$result), 3L)
})

test_that("goric mlm: geschatte coefficienten als vector met namen van vcov()", {
  res <- run_goric_mv(hypotheses = list(H1 = H_mv), comparison = "unconstrained")
  expect_identical(colnames(res$ormle$b.restr), rownames(vcov(fit_mv)))
  expect_equal(unlist(res$ormle$b.restr["unconstrained", ]),
               setNames(as.vector(coef(fit_mv)), rownames(vcov(fit_mv))))
})

test_that("goric mlm: gorica op mlm gelijk aan gorica op coef + VCOV", {
  res_obj <- run_goric_mv(hypotheses = list(H1 = H_mv), type = "gorica")
  est <- setNames(as.vector(coef(fit_mv)),
                  gsub("[:()]", ".", rownames(vcov(fit_mv))))
  N <- nrow(df_mlm)
  VCOV <- vcov(fit_mv) * fit_mv$df.residual / N
  res_est <- suppressMessages(goric(est, VCOV = VCOV, hypotheses = list(H1 = H_mv),
                                    type = "gorica"))
  expect_equal(res_obj$result$gorica, res_est$result$gorica, tolerance = 1e-6)
})
