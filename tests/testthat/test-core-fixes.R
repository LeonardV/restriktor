# [CHANGE 2026-10 | audit] new file: regression tests for the fixes in goric.default() and remove_redundant_constraints() (E4, E7, E13, R6, R7)
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
                     priorICweights = c(NA, 1)), "NA")
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                     priorICweights = c(Inf, 1)), "finite")
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                     priorICweights = c(0, 0)), "at least one positive")
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp,
                     priorICweights = c(0.2, 0.3, 0.5)), "failsafe")
  expect_error(goric(est_cf, VCOV = VCOV_cf, hypotheses = hyp, comparison = "none",
                     priorICweights = c(0.5, 0.5)), "should consist of 1 element,")
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


# -----------------------------------------------------------------------------
# IC-gewichten: nul-prior en grote IC-verschillen (log-sum-exp)
# -----------------------------------------------------------------------------

test_that("goric: priorICweight 0 geeft gewicht 0 en 1 (geen NaN)", {
  res <- goric(c(x = 100), VCOV = matrix(1), hypotheses = list(H = "x > 0"),
               priorICweights = c(0, 1))
  expect_equal(res$result$gorica.weights, c(0, 1))
  expect_false(anyNA(res$result$gorica.weights))
  expect_equal(res$ratio.gw["H", "vs. complement"], 0)
  expect_equal(res$ratio.gw["complement", "vs. H"], Inf)
  expect_equal(diag(res$ratio.gw), c(1, 1))
  expect_output(print(res))
  expect_output(print(summary(res)))
})

test_that("goric: nul-prior bij unconstrained vergelijking", {
  res <- goric(est_cf, VCOV = VCOV_cf,
               hypotheses = list(H1 = "x1 > x2", H2 = "x1 < x2"),
               priorICweights = c(0.5, 0.5, 0))
  expect_equal(res$result$gorica.weights[3], 0)
  expect_equal(sum(res$result$gorica.weights), 1)
  expect_false(anyNA(res$result$gorica.weights))
  expect_output(print(res))
  expect_output(print(summary(res)))
})

test_that("goric: enorm IC-verschil geeft eindige gewichten 1/0 zonder NaN", {
  # x = 100 met variantie 1: loglik-verschil H vs complement is enorm
  res <- goric(c(x = 100), VCOV = matrix(1), hypotheses = list(H = "x < 0"))
  w <- res$result$gorica.weights
  expect_false(anyNA(w))
  expect_equal(w, c(0, 1))
  expect_equal(res$ratio.gw["complement", "vs. H"], Inf)
  expect_equal(res$ratio.gw["H", "vs. complement"], 0)
  expect_false(anyNA(res$result$loglik.weights))
  expect_false(anyNA(res$result$penalty.weights))
  expect_output(print(res))
  expect_output(print(summary(res)))
})

test_that("ic_weights_log: gelijk aan de directe formule bij gewone waarden", {
  IC <- c(10, 12, 15); prior <- c(0.5, 0.3, 0.2)
  w_direct <- prior * exp(-IC / 2) / sum(prior * exp(-IC / 2))
  expect_equal(restriktor:::ic_weights_log(IC, prior), w_direct)
  expect_equal(restriktor:::ic_weights_log(IC), exp(-IC / 2) / sum(exp(-IC / 2)))
  expect_equal(restriktor:::ic_weights_log(c(10, 2000)), c(1, 0))
})


# -----------------------------------------------------------------------------
# gorica voor lavaan: gestandaardiseerde schattingen via merge op parameter
# -----------------------------------------------------------------------------

test_that("goric.lavaan standardized: := en == rijen worden correct gekoppeld", {
  skip_if_not_installed("lavaan")
  set.seed(1)
  d <- data.frame(x = rnorm(100))
  d$y <- 0.5 * d$x + rnorm(100)
  d$z <- 0.4 * d$x + rnorm(100)
  fit <- lavaan::sem("y ~ a*x\nz ~ b*x\na == b\ntotal := a+b", data = d)
  est <- restriktor:::con_gorica_est_lav(fit, standardized = TRUE)
  std <- lavaan::standardizedSolution(fit)
  expect_equal(unname(est$estimate["total"]), std$est.std[std$label == "total"])
  expect_equal(unname(est$estimate["a"]), std$est.std[std$label == "a"])
  # dezelfde structuur als het ongestandaardiseerde pad
  est_u <- restriktor:::con_gorica_est_lav(fit, standardized = FALSE)
  expect_identical(names(est$estimate), names(est_u$estimate))
  expect_identical(dim(est$VCOV), dim(est_u$VCOV))
  res <- suppressMessages(goric(fit, hypotheses = list(H = "total > 0"),
                                standardized = TRUE))
  expect_s3_class(res, "con_gorica")
  expect_equal(unname(res$b.unrestr["total"]), std$est.std[std$label == "total"])
})

test_that("goric.lavaan standardized: multigroep met group.equal", {
  skip_if_not_installed("lavaan")
  set.seed(2)
  d <- data.frame(x = rnorm(120), g = rep(c("A", "B"), 60))
  d$y <- 0.5 * d$x + rnorm(120)
  d$z <- 0.4 * d$x + rnorm(120)
  fit <- lavaan::sem("y ~ c(a1, a2)*x\nz ~ b*x\ntot := a1 + a2", data = d,
                     group = "g", group.equal = "regressions")
  est_s <- restriktor:::con_gorica_est_lav(fit, standardized = TRUE)
  est_u <- restriktor:::con_gorica_est_lav(fit, standardized = FALSE)
  expect_identical(names(est_s$estimate), names(est_u$estimate))
  expect_identical(est_s$rhs, est_u$rhs)
  std <- lavaan::standardizedSolution(fit)
  expect_equal(unname(est_s$estimate["tot"]), std$est.std[std$label == "tot"])
  expect_equal(unname(est_s$estimate["a1"]),
               std$est.std[std$label == "a1" & std$group == 1])
  expect_equal(unname(est_s$estimate["a2"]),
               std$est.std[std$label == "a2" & std$group == 2])
  res <- suppressMessages(goric(fit, hypotheses = list(H = "a1 > a2"),
                                standardized = TRUE))
  expect_s3_class(res, "con_gorica")
})

test_that("goric.lavaan standardized: model zonder restricties of := ", {
  skip_if_not_installed("lavaan")
  set.seed(3)
  d <- data.frame(x = rnorm(100))
  d$y <- 0.5 * d$x + rnorm(100)
  d$z <- 0.4 * d$x + rnorm(100)
  fit <- lavaan::sem("y ~ a*x\nz ~ b*x", data = d)
  est_s <- restriktor:::con_gorica_est_lav(fit, standardized = TRUE)
  est_u <- restriktor:::con_gorica_est_lav(fit, standardized = FALSE)
  expect_identical(names(est_s$estimate), names(est_u$estimate))
  expect_identical(dim(est_s$VCOV), dim(est_u$VCOV))
  std <- lavaan::standardizedSolution(fit)
  expect_equal(unname(est_s$estimate["a"]), std$est.std[std$label == "a"])
  expect_equal(unname(est_s$estimate["b"]), std$est.std[std$label == "b"])
  expect_equal(unname(est_u$estimate["a"]), lavaan::coef(fit)[["a"]])
})


# -----------------------------------------------------------------------------
# mlm: goricc/goricac geblokkeerd
# -----------------------------------------------------------------------------

test_that("goric: goricc/goricac geven een fout voor mlm-objecten", {
  set.seed(5)
  d <- data.frame(x = rnorm(40))
  d$y1 <- d$x + rnorm(40); d$y2 <- 2 * d$x + rnorm(40)
  fit_mv <- lm(cbind(y1, y2) ~ x, data = d)
  expect_error(goric(fit_mv, hypotheses = list(H = "y1.x > 0"), type = "goricc"),
               "not \\(yet\\) available for objects of class mlm")
  expect_error(goric(fit_mv, hypotheses = list(H = "y1.x > 0"), type = "goricac"),
               "not \\(yet\\) available for objects of class mlm")
  # goric en gorica blijven werken
  expect_s3_class(suppressMessages(goric(fit_mv, hypotheses = list(H = "y1.x > 0"),
                                         type = "goric")), "con_goric")
})
