# Tests for goric.Arima (goric_arima.R), arme() and goric.arme (arme.R),
# and the shared machinery (goric_tsme.R)

sim_ar1_me <- function(n, phi, lambda, mu = 0) {
  sx2 <- 1 / (1 - phi^2)
  r   <- sx2 * (1 - lambda) / lambda
  mu + as.numeric(arima.sim(list(ar = phi), n = n)) + rnorm(n, sd = sqrt(r))
}

test_that("back-transformation ARMA(1,1) -> AR(1) + ME reproduces the population values", {
  phi <- 0.6; s2 <- 1; r <- 0.5
  sx2 <- s2 / (1 - phi^2)
  g0 <- sx2 + r
  c0 <- s2 + (1 + phi^2) * r; c1 <- -phi * r
  rho <- c1 / c0
  theta <- (1 - sqrt(1 - 4 * rho^2)) / (2 * rho)      # invertible root
  sa2 <- c1 / theta
  tr <- restriktor:::arima_me_transform(c(phi, theta, sa2))
  expect_equal(unname(tr["r"]), r, tolerance = 1e-8)
  expect_equal(unname(tr["s2"]), s2, tolerance = 1e-8)
  expect_equal(unname(tr["gamma0"]), g0, tolerance = 1e-8)
  expect_equal(unname(tr["pme"]), r / g0, tolerance = 1e-8)
  expect_equal(unname(tr["lambda"] + tr["pme"]), 1, tolerance = 1e-12)
})

test_that("default 'relevance' set runs and returns a goric_arima object", {
  set.seed(1)
  y <- sim_ar1_me(300, 0.6, 0.7, mu = 5)
  fit <- arima(y, order = c(1, 0, 1), method = "ML")
  out <- suppressMessages(goric(fit, type = "gorica"))
  expect_s3_class(out, "goric_arima")
  expect_s3_class(out, "con_goric")
  expect_equal(out$result$model, c("H_negl", "H_ME", "H_neg"))
  expect_equal(sum(out$result$gorica.weights), 1, tolerance = 1e-8)
  expect_true(out$tsme$me)
  expect_true(all(c("r", "gamma0", "lambda", "pme") %in% rownames(out$tsme$estimates)))
  expect_equal(out$tsme$delta, 0.1)
  expect_output(print(out), "measurement error")
})

test_that("'delta' is inserted and the sensitivity analysis is consistent", {
  set.seed(2)
  y <- sim_ar1_me(300, 0.6, 0.7)
  fit <- arima(y, order = c(1, 0, 1), method = "ML")
  out <- suppressMessages(goric(fit, type = "gorica", delta = 0.2, sensitivity = TRUE,
                                delta_values = c(0.1, 0.2, 0.3)))
  expect_true(grepl("0.2", out$tsme$hypotheses_used$H_ME))
  s <- out$tsme$sensitivity
  expect_equal(s$delta, c(0.1, 0.2, 0.3))
  expect_equal(unname(unlist(s[s$delta == 0.2, c("H_negl", "H_ME", "H_neg")])),
               out$result$gorica.weights, tolerance = 1e-6)
  expect_true(all(diff(s$H_negl) >= -1e-8))
})

test_that("other default sets and user hypotheses work", {
  set.seed(3)
  y <- sim_ar1_me(300, 0.7, 0.6)
  fit <- arima(y, order = c(1, 0, 1), method = "ML")
  z <- suppressMessages(goric(fit, default_hypotheses = "zero"))
  expect_equal(z$result$model, c("H0", "Hpos", "Hneg"))
  s <- suppressMessages(goric(fit, default_hypotheses = "sign"))
  expect_equal(s$result$model, c("Hpos", "complement"))
  m <- suppressMessages(goric(fit, default_hypotheses = "minimum"))
  expect_equal(m$result$model, c("H_ME", "complement"))
  u <- suppressMessages(goric(fit, hypotheses = list(H1 = "lambda > 1 - delta; ar1 > 0"),
                              delta = 0.2))
  expect_true(all(c("lambda", "ar1") %in%
                    rownames(u$tsme$estimates)[u$tsme$estimates$involved]))
})

test_that("non-ARMA(1,1) models: coefficient hypotheses work, ME parameters give an error", {
  fit <- arima(lh, order = c(2, 0, 0))
  out <- suppressMessages(goric(fit, hypotheses = list(H1 = "ar1 > 0; ar2 < 0")))
  expect_s3_class(out, "goric_arima")
  expect_false(out$tsme$me)
  expect_error(suppressMessages(goric(fit)), "only available for models with measurement error")
  expect_error(suppressMessages(goric(fit, hypotheses = list(H1 = "r > 0"))),
               "only available for")
})

test_that("combining functionally dependent parameters gives an informative error", {
  set.seed(4)
  y <- sim_ar1_me(300, 0.6, 0.7)
  fit <- arima(y, order = c(1, 0, 1), method = "ML")
  expect_error(suppressMessages(goric(fit, hypotheses = list(H1 = "lambda > 0.5; pme < 0.5"))),
               "singular")
})

test_that("GORIC: log-likelihoods are the (order-restricted) maximum log-likelihoods", {
  set.seed(5)
  y <- sim_ar1_me(250, 0.6, 0.7)
  fit <- arima(y, order = c(1, 0, 1), method = "ML")
  # equality restriction on a coefficient: compare with arima(..., fixed = )
  g <- suppressMessages(goric(fit, hypotheses = list(H1 = "ar1 = 0.3"), type = "goric",
                              posthoc = FALSE))
  f03 <- arima(y, order = c(1, 0, 1), fixed = c(0.3, NA, NA), transform.pars = FALSE,
               method = "ML")
  expect_equal(g$result$loglik[1], f03$loglik, tolerance = 1e-4)
  expect_equal(g$result$loglik[g$result$model == "unconstrained"], fit$loglik, tolerance = 1e-4)
  # default set: the best hypothesis contains the unconstrained estimate
  g2 <- suppressMessages(goric(fit, type = "goric"))
  expect_equal(g2$type, "goric")
  expect_equal(max(g2$result$loglik), fit$loglik, tolerance = 1e-4)
  # restricted estimates satisfy the restrictions
  expect_true(all(abs(g2$ormle$b.restr[c("H_negl", "H_neg"), "pme"]) <= 0.1 + 1e-4))
})

test_that("arme(): AR(1) + ME equals the ARMA(1,1) fit; likelihood equals the Gaussian one", {
  set.seed(6)
  y <- sim_ar1_me(150, 0.6, 0.7, mu = 2)
  fa <- arima(y, order = c(1, 0, 1), method = "ML")
  fm <- arme(y)
  expect_equal(fm$loglik, fa$loglik, tolerance = 1e-4)
  # Kalman likelihood versus multivariate normal (also for a negative r)
  sp <- restriktor:::arme_spec_uni(1L)
  th <- c(mu = 2, phi = 0.6, s2 = 1, r = -0.1)
  g <- (1 / (1 - 0.36)) * 0.6^(0:149); g[1] <- g[1] - 0.1
  S <- stats::toeplitz(g)
  ll_mvn <- -0.5 * (150 * log(2 * pi) + as.numeric(determinant(S)$modulus) +
                      sum((y - 2) * solve(S, y - 2)))
  expect_equal(sp$loglik(th, y), ll_mvn, tolerance = 1e-8)
})

test_that("arme(): VAR(1) + ME fit and goric()", {
  skip_on_cran()
  set.seed(10)
  sp <- restriktor:::arme_spec_var(2)
  th <- c(mu1 = 0, mu2 = 1, Phi11 = .5, Phi12 = .2, Phi21 = -.1, Phi22 = .4,
          Sigma11 = 1, Sigma21 = .3, Sigma22 = 1, r1 = .6, r2 = .4)
  Y <- sp$simulate(th, 300)
  fv <- arme(Y)
  expect_s3_class(fv, "arme")
  expect_equal(length(coef(fv)), 11)
  out <- suppressMessages(goric(fv))
  expect_equal(out$result$model, c("H_negl", "H_ME1", "H_ME2", "H_ME12", "unconstrained"))
})

test_that("bootstrap covariance matrix and benchmark_tsme() run and are reproducible", {
  skip_on_cran()
  set.seed(8)
  y <- sim_ar1_me(200, 0.6, 0.7)
  fit <- arima(y, order = c(1, 0, 1), method = "ML")
  b1 <- suppressMessages(goric(fit, type = "gorica", VCOV_method = "bootstrap", B = 50, seed = 10))
  b2 <- suppressMessages(goric(fit, type = "gorica", VCOV_method = "bootstrap", B = 50, seed = 10))
  expect_equal(b1$result$gorica.weights, b2$result$gorica.weights)
  out <- suppressMessages(goric(fit, type = "gorica", default_hypotheses = "minimum"))
  bm <- benchmark_tsme(out, B = 30, seed = 1)
  expect_s3_class(bm, "benchmark_tsme")
  expect_true(all(bm$percentile >= 0 & bm$percentile <= 100))
})

test_that("type: GORICA by default; GORIC on request (GORICA if the data are not found)", {
  set.seed(9)
  y <- sim_ar1_me(200, 0.6, 0.7)
  fit <- arima(y, order = c(1, 0, 1), method = "ML")
  expect_equal(suppressMessages(goric(fit))$type, "gorica")
  expect_equal(suppressMessages(goric(fit, type = "goric"))$type, "goric")
  expect_equal(suppressMessages(goric(fit, type = NULL))$type, "goric")
  # the data cannot be found from the call -> GORICA, with a message
  make_fit <- function() { yy <- sim_ar1_me(200, 0.6, 0.7); arima(yy, order = c(1, 0, 1), method = "ML") }
  fit2 <- make_fit()
  expect_message(g2 <- goric(fit2, type = NULL), "GORICA is used instead")
  expect_equal(g2$type, "gorica")
  expect_equal(suppressMessages(goric(fit2, type = "goric"))$type, "gorica")
})
