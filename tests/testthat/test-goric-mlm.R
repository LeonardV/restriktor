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

# onafhankelijke referentiewaarden voor H_mv: de ongerestricteerde schattingen
# als vector (volgorde en namen van vcov(), zoals in ormle$b.restr) en de
# restricties (Age:GroupPassive = 0 en Age:GroupNo - Age:GroupPassive >= 0)
b_mv_unr <- setNames(as.vector(coef(fit_mv)), rownames(vcov(fit_mv)))
check_H_mv <- function(b) {
  isTRUE(all.equal(unname(b["Age:GroupPassive"]), 0, tolerance = 1e-8)) &&
    b["Age:GroupNo"] - b["Age:GroupPassive"] >= -1e-8
}

test_that("goric mlm: alle combinaties van comparison en type draaien", {
  for (cmp in c("complement", "unconstrained", "none")) {
    for (ty in c("goric", "gorica")) {
      res <- run_goric_mv(hypotheses = list(H1 = H_mv), comparison = cmp, type = ty)
      expect_s3_class(res, "con_goric")
      expect_true(all(is.finite(res$result[[ty]])))
      # IC = -2 * loglik + 2 * penalty
      expect_equal(res$result[[ty]], -2 * res$result$loglik + 2 * res$result$penalty,
                   tolerance = 1e-10)
      # de IC-gewichten sommeren tot 1 en volgen uit de IC-waarden
      w <- res$result[[paste0(ty, ".weights")]]
      expect_equal(sum(w), 1, tolerance = 1e-10)
      ic <- res$result[[ty]]
      expect_equal(w, exp(-(ic - min(ic)) / 2) / sum(exp(-(ic - min(ic)) / 2)),
                   tolerance = 1e-8)
      # de gerestricteerde schattingen van H1 voldoen aan de restricties
      b_H1 <- unlist(res$ormle$b.restr["H1", names(b_mv_unr)])
      expect_true(check_H_mv(b_H1))
      # loglik van H1 <= loglik van het ongerestricteerde model (bij goric de
      # MVN-loglik; bij gorica 0 = dmvnorm(0) is het maximum)
      ll_max <- if (ty == "goric") ll_mv_ref else
        mvtnorm::dmvnorm(rep(0, length(b_mv_unr)), sigma = res$VCOV, log = TRUE)
      expect_true(res$result$loglik[1] <= ll_max + 1e-8)
      if (cmp != "none") {
        # complement/unconstrained: loglik <= ongerestricteerd en >= H1
        expect_true(res$result$loglik[2] <= ll_max + 1e-8)
        expect_true(res$result$loglik[2] >= res$result$loglik[1] - 1e-8)
        # de ongerestricteerde schattingen staan in de laatste rij
        if (cmp == "unconstrained") {
          expect_equal(unlist(res$ormle$b.restr["unconstrained", names(b_mv_unr)]),
                       b_mv_unr, tolerance = 1e-8)
        }
        # PT van het ongerestricteerde model = 1 + p (goric) resp. p (gorica)
        if (cmp == "unconstrained") {
          expect_equal(res$result$penalty[2],
                       length(b_mv_unr) + (ty == "goric"), tolerance = 1e-10)
        }
      }
    }
  }
})

test_that("goric mlm: loglik van complement en unconstrained is de MVN-loglik", {
  res_c <- run_goric_mv(hypotheses = list(H1 = H_mv), comparison = "complement")
  expect_equal(res_c$result$loglik[2], ll_mv_ref, tolerance = 1e-8)
  res_u <- run_goric_mv(hypotheses = list(H1 = H_mv), comparison = "unconstrained")
  expect_equal(res_u$result$loglik[2], ll_mv_ref, tolerance = 1e-8)
})

test_that("goric mlm: goricc en goricac geven een duidelijke fout", {
  # de small-sample correctie is (nog) niet afgeleid voor een multivariate
  # residuele covariantiematrix; de univariate correctie mag niet gebruikt worden
  for (ty in c("goricc", "goricac", "GORICC")) {
    expect_error(run_goric_mv(hypotheses = list(H1 = H_mv), type = ty),
                 "not \\(yet\\) available for objects of class mlm")
  }
})

test_that("goric mlm: sample size is N en niet N x aantal responsen", {
  res <- run_goric_mv(hypotheses = list(H1 = H_mv), comparison = "unconstrained")
  expect_equal(res$sample_nobs, nrow(df_mlm))
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
