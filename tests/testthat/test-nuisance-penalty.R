# Penalties with nuisance parameters (d) and the multivariate (mlm) case
#   GORIC  : sum_l LP_l K_l,                  K_l = l + d
#   GORICC : sum_l LP_l N K_l / (N - K_l - 1)   (mlm: N - l/q - q - 1)
#   GORICA : GORIC  - d
#   GORICAC: GORICC - d

set.seed(5)
tg <- ToothGrowth
tg$dose <- factor(tg$dose)
tg$len2 <- tg$len + rnorm(nrow(tg), 2)
fit_lm  <- lm(len ~ 0 + dose, data = tg)
fit_mlm <- lm(cbind(len, len2) ~ 0 + dose, data = tg)
H_lm  <- list(H1 = "dose0.5 < dose1 < dose2")
H_mlm <- list(H1 = "len.dose0.5 < len.dose1")
N <- nrow(tg)

test_that("lm: GORICA(C) = GORIC(C) - 1 (d = 1)", {
  pen <- function(ty) suppressMessages(goric(fit_lm, hypotheses = H_lm, type = ty))$result$penalty
  expect_equal(pen("gorica"), pen("goric") - 1)
  expect_equal(pen("goricac"), pen("goricc") - 1)
})

test_that("lm: unconstrained GORICC penalty is N (p + 1) / (N - p - 2)", {
  res <- suppressMessages(goric(fit_lm, hypotheses = list(H1 = "dose0.5 < dose1", H2 = "dose1 < dose2"),
                                type = "goricc"))
  p <- 3
  expect_equal(res$result$penalty[3], N * (p + 1) / (N - p - 2))
})

test_that("mlm: d = q(q+1)/2 and the Bedrick-Tsai correction", {
  pen <- function(ty) suppressMessages(goric(fit_mlm, hypotheses = H_mlm, type = ty))$result$penalty
  q <- 2; m <- q * (q + 1) / 2
  expect_equal(pen("goric"), pen("gorica") + m)
  expect_equal(pen("goricac"), pen("goricc") - m)
  # unconstrained GORICC equals the corrected AIC of Bedrick and Tsai (1994)
  res <- suppressMessages(goric(fit_mlm, hypotheses = list(H1 = "len.dose0.5 < len.dose1",
                                                         H2 = "len2.dose0.5 > len2.dose1"),
                                type = "goricc"))
  k <- 3
  S <- crossprod(resid(fit_mlm)) / N
  aicc_bt <- N * log(det(S)) + N * (N + k) * q / (N - k - q - 1) + N * q * log(2 * pi)
  expect_equal(res$result$goricc[3], aicc_bt)
  expect_equal(as.vector(res$nr_nuis_param), m)
  expect_equal(attr(res$nr_nuis_param, "n_resp"), q)
})

test_that("binomial glm: d = 0, so GORICC = GORICAC and GORIC = GORICA", {
  set.seed(1)
  dd <- data.frame(g = factor(rep(c("a", "b", "c"), each = 40)))
  dd$z <- rbinom(120, 1, plogis(as.numeric(dd$g) - 2))
  gl <- glm(z ~ 0 + g, family = binomial, data = dd)
  pen <- function(ty) suppressMessages(goric(gl, hypotheses = list(H1 = "ga < gb < gc"),
                                             type = ty))$result$penalty
  expect_equal(pen("goric"), pen("gorica"))
  expect_equal(pen("goricc"), pen("goricac"))
})

test_that("est + VCOV: goricac requires nr_nuis_param; d = 1 equals the lm object", {
  est  <- coef(fit_lm)
  VCOV <- vcov(fit_lm) * (N - 3) / N
  expect_error(goric(est, VCOV = VCOV, hypotheses = H_lm, type = "goricac", sample_nobs = N),
               "nr_nuis_param")
  r_est <- suppressMessages(goric(est, VCOV = VCOV, hypotheses = H_lm, type = "goricac",
                                  sample_nobs = N, nr_nuis_param = 1))
  r_lm  <- suppressMessages(goric(fit_lm, hypotheses = H_lm, type = "goricac"))
  expect_equal(r_est$result$goricac, r_lm$result$goricac)
  # more nuisance parameters: larger small-sample correction
  r_est3 <- suppressMessages(goric(est, VCOV = VCOV, hypotheses = H_lm, type = "goricac",
                                   sample_nobs = N, nr_nuis_param = 3))
  expect_true(all(r_est3$result$penalty > r_est$result$penalty))
  # the GORICA does not depend on it
  r_a1 <- suppressMessages(goric(est, VCOV = VCOV, hypotheses = H_lm, type = "gorica"))
  r_a3 <- suppressMessages(goric(est, VCOV = VCOV, hypotheses = H_lm, type = "gorica",
                                 nr_nuis_param = 3))
  expect_equal(r_a1$result$penalty, r_a3$result$penalty)
  expect_error(goric(est, VCOV = VCOV, hypotheses = H_lm, type = "goricac", sample_nobs = N,
                     nr_nuis_param = -1), "nr_nuis_param")
})

test_that("lavaan: d = free parameters not involved in the hypotheses", {
  skip_if_not_installed("lavaan")
  set.seed(1)
  dl <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
  dl$y <- .3 * dl$x1 + .2 * dl$x2 + rnorm(100)
  fit_lav <- lavaan::sem("y ~ b1*x1 + b2*x2 + b3*x3", data = dl)
  res <- suppressMessages(goric(fit_lav, hypotheses = list(H1 = "b1 > b2"), type = "goricac"))
  # 4 free parameters (b1, b2, b3, residual variance), 2 involved
  expect_equal(as.vector(res$nr_nuis_param), 2)
})
