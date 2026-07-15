## tests for missing = "fiml" in restriktor() and goric()

# simulate data with MAR missingness
gen_fiml_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- 0.5 * x1 + rnorm(n); x3 <- rnorm(n)
  y  <- 1 + 0.8 * x1 + 0.4 * x2 - 0.3 * x3 + rnorm(n)
  aux <- 0.6 * y + rnorm(n)
  df <- data.frame(y, x1, x2, x3, aux)
  # MAR: missingness depends on observed variables
  df$y [df$x1 + rnorm(n, sd = 0.8) >  1.0] <- NA
  df$x2[df$aux + rnorm(n, sd = 0.8) < -1.2] <- NA
  df$x3[runif(n) < 0.1] <- NA
  df
}

test_that("unrestricted FIML estimates and SEs match lavaan", {
  skip_if_not_installed("lavaan")
  df <- gen_fiml_data()
  fit <- lm(y ~ x1 + x2 + x3, data = df)

  fiml <- restriktor:::con_fiml_lm(fit)

  lav <- lavaan::sem("y ~ x1 + x2 + x3", data = df, missing = "ml",
                     fixed.x = FALSE, meanstructure = TRUE)
  pe <- lavaan::parameterEstimates(lav)
  lav_coef <- c(pe[pe$op == "~1" & pe$lhs == "y", "est"],
                pe[pe$op == "~", "est"])
  lav_se <- c(pe[pe$op == "~1" & pe$lhs == "y", "se"],
              pe[pe$op == "~", "se"])

  expect_equal(unname(fiml$est), lav_coef, tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(fiml$VCOV))), lav_se, tolerance = 1e-3)
  # the linear model is just-identified, so its FIML log-likelihood equals
  # the saturated log-likelihood
  expect_equal(fiml$loglik, as.numeric(lavaan::fitMeasures(lav, "logl")), tolerance = 1e-5)
})


test_that("FIML with auxiliary variable matches lavaan saturated-correlates", {
  skip_if_not_installed("lavaan")
  df <- gen_fiml_data()
  fit <- lm(y ~ x1 + x2 + x3, data = df)

  fiml <- restriktor:::con_fiml_lm(fit, auxiliary = "aux")

  model <- "y ~ x1 + x2 + x3
            aux ~~ y + x1 + x2 + x3
            x1 ~~ x2 + x3
            x2 ~~ x3"
  lav <- lavaan::sem(model, data = df, missing = "ml", fixed.x = FALSE,
                     meanstructure = TRUE)
  pe <- lavaan::parameterEstimates(lav)
  lav_coef <- c(pe[pe$op == "~1" & pe$lhs == "y", "est"],
                pe[pe$op == "~", "est"])
  lav_se <- c(pe[pe$op == "~1" & pe$lhs == "y", "se"],
              pe[pe$op == "~", "se"])

  expect_equal(unname(fiml$est), lav_coef, tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(fiml$VCOV))), lav_se, tolerance = 1e-3)
})


test_that("restricted FIML (ECM) log-likelihood matches lavaan at active constraint", {
  skip_if_not_installed("lavaan")
  df <- gen_fiml_data()
  fit <- lm(y ~ x1 + x2 + x3, data = df)
  fiml <- restriktor:::con_fiml_lm(fit)

  # active constraint b_x1 <= 0.2
  ecm <- restriktor:::fiml_ecm_restricted(fiml, Amat = rbind(c(0, -1, 0, 0)),
                                          bvec = -0.2, meq = 0)
  expect_lte(ecm$b.restr[2], 0.2 + 1e-8)
  expect_lte(ecm$loglik, fiml$loglik)

  # since the constraint is active, the solution equals the equality-restricted fit
  lav <- lavaan::sem("y ~ b1*x1 + x2 + x3\n b1 == 0.2", data = df,
                     missing = "ml", fixed.x = FALSE, meanstructure = TRUE)
  expect_equal(ecm$loglik, as.numeric(lavaan::fitMeasures(lav, "logl")), tolerance = 1e-4)
})


test_that("complete data: fiml equals listwise/OLS results", {
  df <- gen_fiml_data()
  df_c <- df[complete.cases(df), ]
  fit <- lm(y ~ x1 + x2 + x3, data = df_c)

  r_fiml <- restriktor(fit, constraints = "x1 < 0.2", missing = "fiml")
  r_none <- restriktor(fit, constraints = "x1 < 0.2")

  expect_equal(coef(r_fiml), coef(r_none), tolerance = 1e-4)
  expect_equal(r_fiml$b.unrestr, r_none$b.unrestr, tolerance = 1e-6)

  # GORIC weights are invariant to the constant f(x) part of the joint likelihood
  g_fiml <- goric(fit, hypotheses = list(H1 = "x1 > x2 > 0"),
                  comparison = "complement", missing = "fiml")
  g_none <- goric(fit, hypotheses = list(H1 = "x1 > x2 > 0"),
                  comparison = "complement")
  expect_equal(diff(g_fiml$result$loglik), diff(g_none$result$loglik),
               tolerance = 1e-3)
  expect_equal(g_fiml$result$goric.weights, g_none$result$goric.weights,
               tolerance = 1e-3)
})


test_that("restriktor() with missing = fiml returns a valid object", {
  df <- gen_fiml_data()
  fit <- lm(y ~ x1 + x2 + x3, data = df)

  r1 <- restriktor(fit, constraints = "x1 < 0.2", missing = "fiml")
  expect_s3_class(r1, "restriktor")
  expect_identical(r1$missing, "fiml")
  expect_lte(coef(r1)["x1"], 0.2 + 1e-8)
  # restricted loglik cannot exceed the unrestricted loglik
  expect_lte(logLik(r1), r1$fiml$loglik.unrestr)
  # N = all rows with at least one observed variable
  expect_equal(r1$fiml$N, sum(rowSums(!is.na(df[, 1:4])) > 0))
  # summary and print work
  expect_error(capture.output(print(summary(r1))), NA)
  s1 <- summary(r1)
  expect_false(any(is.na(s1$coefficients[, "Estimate"])))

  # inactive constraint: restricted = unrestricted
  r2 <- restriktor(fit, constraints = "x1 > 0", missing = "fiml")
  expect_equal(coef(r2), r2$b.unrestr, tolerance = 1e-8)

  # auxiliary variable changes (improves) the estimates
  r3 <- restriktor(fit, constraints = "x1 > 0", missing = "fiml",
                   auxiliary = "aux")
  expect_identical(r3$auxiliary, "aux")
  expect_false(isTRUE(all.equal(coef(r2), coef(r3))))
})


test_that("goric() with missing = fiml supports goric, goricc and gorica", {
  df <- gen_fiml_data()
  fit <- lm(y ~ x1 + x2 + x3, data = df)
  N <- sum(rowSums(!is.na(df[, 1:4])) > 0)

  # loglik-based GORIC with complement
  g1 <- goric(fit, hypotheses = list(H1 = "x1 > x2 > 0"), missing = "fiml",
              comparison = "complement")
  expect_false(any(is.na(g1$result$goric)))
  expect_equal(g1$sample_nobs, N)

  # unconstrained comparison: unconstrained loglik is the maximum
  g2 <- goric(fit, hypotheses = list(H1 = "x1 > x2", H2 = "x1 < x2"),
              missing = "fiml")
  llu <- g2$result$loglik[g2$result$model == "unconstrained"]
  expect_true(all(g2$result$loglik <= llu + 1e-6))
  # H1 holds in the data, so it gets the same loglik as the unconstrained model
  expect_equal(g2$result$loglik[g2$result$model == "H1"], llu, tolerance = 1e-6)

  # gorica and goricc
  g3 <- goric(fit, hypotheses = list(H1 = "x1 > x2 > 0"), missing = "fiml",
              type = "gorica", comparison = "complement")
  expect_false(any(is.na(g3$result$gorica)))
  g4 <- goric(fit, hypotheses = list(H1 = "x1 > x2"), missing = "fiml",
              type = "goricc")
  expect_false(any(is.na(g4$result$goricc)))

  # auxiliary variables
  g5 <- goric(fit, hypotheses = list(H1 = "x1 > x2 > 0"), missing = "fiml",
              auxiliary = "aux", comparison = "complement")
  expect_false(any(is.na(g5$result$goric)))

  # deprecated two.stage alias maps onto fiml
  expect_message(
    g6 <- goric(fit, hypotheses = list(H1 = "x1 > x2"), missing = "two.stage"),
    "fiml"
  )
  g7 <- goric(fit, hypotheses = list(H1 = "x1 > x2"), missing = "fiml")
  expect_equal(g6$result$goric, g7$result$goric, tolerance = 1e-6)
})


test_that("fiml works with factors and interactions", {
  set.seed(7)
  n <- 300
  x1 <- rnorm(n)
  g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  y <- 1 + 0.8 * x1 + 0.5 * (g == "b") - 0.5 * (g == "c") + rnorm(n)
  df <- data.frame(y, x1, g)
  df$y[runif(n) < 0.15] <- NA
  df$x1[runif(n) < 0.15] <- NA
  fit <- lm(y ~ x1 + g, data = df)

  r <- restriktor(fit, constraints = "gb > gc", missing = "fiml")
  expect_identical(names(coef(r)), names(coef(fit)))
  expect_gte(coef(r)["gb"], coef(r)["gc"])
})


test_that("fiml error handling", {
  df <- gen_fiml_data()
  fit <- lm(y ~ x1 + x2 + x3, data = df)

  # bootstrapped standard errors are not available
  expect_error(restriktor(fit, constraints = "x1 > 0", missing = "fiml",
                          se = "boot.standard"), "fiml")
  # informative hypothesis tests are not available
  r <- restriktor(fit, constraints = "x1 > 0", missing = "fiml")
  expect_error(iht(r), "fiml")
  # unknown missing method
  expect_error(restriktor(fit, constraints = "x1 > 0", missing = "foo"),
               "unknown")
  # auxiliary without fiml
  expect_error(restriktor(fit, constraints = "x1 > 0", auxiliary = "aux"),
               "auxiliary")
  # no intercept
  fit0 <- lm(y ~ x1 + x2 + x3 - 1, data = df)
  expect_error(restriktor(fit0, constraints = "x1 > 0", missing = "fiml"),
               "intercept")
  # weights
  fitw <- lm(y ~ x1 + x2 + x3, data = df, weights = rep(1, nrow(df)))
  expect_error(restriktor(fitw, constraints = "x1 > 0", missing = "fiml"),
               "weighted")
})
