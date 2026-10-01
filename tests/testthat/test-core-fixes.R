# =============================================================================
# Tests: fixes in goric.default en remove_redundant_constraints
# =============================================================================

set.seed(4)
d_cf <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
d_cf$y <- 0.2 * d_cf$x1 + 0.6 * d_cf$x2 + rnorm(80)
fit_cf <- lm(y ~ x1 + x2, data = d_cf)

est_cf  <- c(x1 = 0.5, x2 = 0.8, x3 = 0.3)
VCOV_cf <- diag(3) * 0.05


# -----------------------------------------------------------------------------
# remove_redundant_constraints: gelijkheden blijven vooraan
# -----------------------------------------------------------------------------

test_that("remove_redundant_constraints: gelijkheden blijven de eerste meq rijen", {
  # x2 = 0 (eq, rhs 0) en x1 >= 1 (ineq, rhs 1): sorteren op rhs zou de
  # ongelijkheid voor de gelijkheid zetten
  Amat <- rbind(c(0, 1), c(1, 0))
  res <- restriktor:::remove_redundant_constraints(Amat, c(0, 1), meq = 1)
  expect_equal(res$meq, 1)
  expect_equal(res$constraints[1, ], c(0, 1))
  expect_equal(res$rhs[1], 0)
  expect_equal(res$constraints[2, ], c(1, 0))
  expect_equal(res$rhs[2], 1)
})

test_that("restriktor: gelijkheidsrestrictie wordt exact gehandhaafd", {
  b <- coef(restriktor(fit_cf, constraints = "x1 > 1; x2 = 0"))
  expect_equal(unname(b["x2"]), 0, tolerance = 1e-10)
  expect_true(b["x1"] >= 1 - 1e-10)
})

test_that("goric: gelijkheidsrestrictie wordt exact gehandhaafd (gorica)", {
  g <- goric(est_cf, VCOV = VCOV_cf,
             hypotheses = list(H1 = "x1 > 1; x2 = 0"))
  expect_equal(unname(coef(g)["H1", "x2"]), 0, tolerance = 1e-10)
  expect_true(coef(g)["H1", "x1"] >= 1 - 1e-10)
})

test_that("remove_redundant_constraints: meq heeft default 0", {
  res <- restriktor:::remove_redundant_constraints(rbind(c(1, 0), c(1, 0)), c(1, 2))
  expect_equal(res$meq, 0)
  expect_equal(res$rhs, 2)
})


# -----------------------------------------------------------------------------
# remove_redundant_constraints: gelijkheid + ongelijkheid op dezelfde rij
# -----------------------------------------------------------------------------

test_that("remove_redundant_constraints: ongelijkheid geimpliceerd door gelijkheid wordt verwijderd", {
  # x1 = 2 en x1 >= 1: x1 >= 1 is overbodig, de gelijkheid blijft
  res <- restriktor:::remove_redundant_constraints(rbind(c(1, 0), c(1, 0)),
                                                   c(2, 1), meq = 1)
  expect_equal(nrow(res$constraints), 1)
  expect_equal(res$meq, 1)
  expect_equal(res$rhs, 2)

  # x1 = 2 en x1 >= 2 (zelfde rhs)
  res <- restriktor:::remove_redundant_constraints(rbind(c(1, 0), c(1, 0)),
                                                   c(2, 2), meq = 1)
  expect_equal(nrow(res$constraints), 1)
  expect_equal(res$meq, 1)
})

test_that("remove_redundant_constraints: gelijkheid gaat niet verloren bij grotere rhs ongelijkheid", {
  # x1 = 0 en x1 >= 1 zijn strijdig -> fout i.p.v. stil verwijderen van de gelijkheid
  expect_error(
    restriktor:::remove_redundant_constraints(rbind(c(1, 0), c(1, 0)),
                                              c(0, 1), meq = 1),
    "conflicting"
  )
})

test_that("remove_redundant_constraints: strijdige gelijkheden geven fout", {
  expect_error(
    restriktor:::remove_redundant_constraints(rbind(c(1, 0), c(1, 0)),
                                              c(1, 2), meq = 2),
    "conflicting"
  )
  # identieke gelijkheden: een blijft over
  res <- restriktor:::remove_redundant_constraints(rbind(c(1, 0), c(1, 0)),
                                                   c(1, 1), meq = 2)
  expect_equal(nrow(res$constraints), 1)
  expect_equal(res$meq, 1)
})

test_that("restriktor: strijdige gelijkheid en ongelijkheid geven fout", {
  expect_error(restriktor(fit_cf, constraints = "x1 = 0; x1 > 1"), "conflicting")
})

test_that("goric: redundante ongelijkheden met Heq = TRUE werken", {
  res <- goric(fit_cf, hypotheses = list(H1 = "x1 > 0.5; x1 > 0.2; x2 > 0"),
               Heq = TRUE)
  expect_s3_class(res, "con_goric")
  expect_equal(res$result$model, c("Heq", "H1", "complement"))
  b_heq <- as.matrix(coef(res))["Heq", ]
  expect_equal(unname(b_heq["x1"]), 0.5, tolerance = 1e-8)
  expect_equal(unname(b_heq["x2"]), 0, tolerance = 1e-8)
})


# -----------------------------------------------------------------------------
# priorICweights validatie
# -----------------------------------------------------------------------------

test_that("goric: priorICweights validatie geeft duidelijke fouten", {
  hyp <- list(H1 = "x1 > x2 > x3")  # complement => 2 modellen
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                     priorICweights = c(-1, 2)), "non-negative")
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                     priorICweights = c(NA, 1)), "missing")
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                     priorICweights = c(0, 0)), "all be zero")
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                     priorICweights = c(0.2, 0.3, 0.5)), "should consist of 2 elements")
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                     priorICweights = c("a", "b")), "numeric")
})

test_that("goric: priorICweights worden herschaald naar som 1", {
  hyp <- list(H1 = "x1 > x2 > x3")
  expect_message(
    res <- goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                 priorICweights = c(3, 1)),
    "rescaled"
  )
  expect_equal(res$priorICweights, c(0.75, 0.25))

  # som is (numeriek) 1: geen herschaling en geen melding
  expect_no_message(
    res2 <- goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                  priorICweights = c(0.1 + 0.2, 0.7))
  )
  expect_equal(res2$priorICweights, c(0.3, 0.7))

  # met Heq: 3 gewichten
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp, Heq = TRUE,
                     priorICweights = c(0.5, 0.5)), "complement")
})


# -----------------------------------------------------------------------------
# ratio-matrices en best_hypo
# -----------------------------------------------------------------------------

test_that("goric: ratio.gw is te indexeren op naam en best_hypo is aanwezig", {
  res <- goric(est_cf, VCOV = VCOV_cf,
               hypotheses = list(H1 = "x1 > x2 > x3", H2 = "x2 > x1 > x3"))
  expect_true(all(c("H1", "H2", "unconstrained") %in% rownames(res$ratio.gw)))
  expect_false(any(grepl("(best)", rownames(res$ratio.gw), fixed = TRUE)))
  expect_false(any(grepl("(best)", rownames(res$ratio.pw), fixed = TRUE)))
  expect_false(any(grepl("(best)", rownames(res$ratio.lw), fixed = TRUE)))
  expect_length(res$ratio.gw["H1", ], 3)
  expect_equal(unname(res$ratio.gw["H1", "vs. H1"]), 1)

  expect_true(!is.null(res$best_hypo))
  expect_equal(res$best_hypo, which.max(res$result$gorica.weights))

  # de markering wordt alleen in de geprinte summary toegevoegd
  out <- capture.output(summary(res))
  best_name <- res$result$model[res$best_hypo]
  expect_true(any(grepl(paste0(best_name, " (best)"), out, fixed = TRUE)))
})


# -----------------------------------------------------------------------------
# verwijderde argumenten add_Hc en posthoc
# -----------------------------------------------------------------------------

test_that("goric: add_Hc en posthoc zijn geen argumenten meer", {
  expect_false("add_Hc" %in% names(formals(restriktor:::goric.default)))
  expect_false("posthoc" %in% names(formals(restriktor:::goric.default)))
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = list(H1 = "x1 > x2"),
                     add_Hc = 1), "Unknown argument")
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = list(H1 = "x1 > x2"),
                     posthoc = TRUE), "Unknown argument")
})
