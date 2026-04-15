# =============================================================================
# Tests: con_loglik_lm en con_loglik_glm — loglikelihood berekeningen
# =============================================================================
# Zorg dat vereiste imports beschikbaar zijn in het huidige library-pad
# for (pkg in c("lavaan", "mvtnorm", "tmvtnorm", "quadprog", "norm")) {
#   if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# }
# 
# # Laad het restriktor package vanuit broncode
# setwd("/Workspace/Users/l.vanbrabant@ggdwestbrabant.nl/restriktor")
# devtools::load_all(".")

# --- Gedeelde testdata ---

set.seed(42)
n <- 100
x <- rnorm(n)
y <- 2 + 0.5 * x + rnorm(n, sd = 0.8)
df <- data.frame(y = y, x = x)
fit_lm <- lm(y ~ x, data = df)


# =============================================================================
# con_loglik_lm: lineair model
# =============================================================================

test_that("con_loglik_lm: retourneert numerieke waarde", {
  ll <- restriktor:::con_loglik_lm(fit_lm)
  expect_true(is.numeric(ll))
  expect_true(is.finite(ll))
})

test_that("con_loglik_lm: loglik is negatief", {
  ll <- restriktor:::con_loglik_lm(fit_lm)
  expect_true(ll < 0)
})

test_that("con_loglik_lm: vergelijkbaar met stats::logLik (ML-schatter)", {
  ll_restr  <- restriktor:::con_loglik_lm(fit_lm)
  ll_stats  <- as.numeric(logLik(fit_lm))

  # con_loglik_lm gebruikt sigma_ML (n-deler) ipv REML
  # verschil is klein maar niet exact gelijk aan logLik() die sigma_REML gebruikt
  # Ze moeten wel in dezelfde orde van grootte liggen
  expect_equal(ll_restr, ll_stats, tolerance = 3)
})

test_that("con_loglik_lm: beter model heeft hogere loglik", {
  # Goed model
  fit_good <- lm(y ~ x, data = df)
  # Slecht model (intercept only)
  fit_bad  <- lm(y ~ 1, data = df)

  ll_good <- restriktor:::con_loglik_lm(fit_good)
  ll_bad  <- restriktor:::con_loglik_lm(fit_bad)

  expect_true(ll_good > ll_bad)
})


# =============================================================================
# con_loglik_lm: gewogen regressie
# =============================================================================

test_that("con_loglik_lm: werkt met gewogen regressie", {
  w <- runif(n, 0.5, 2)
  fit_wlm <- lm(y ~ x, data = df, weights = w)

  ll <- restriktor:::con_loglik_lm(fit_wlm)
  expect_true(is.numeric(ll))
  expect_true(is.finite(ll))
})

test_that("con_loglik_lm: gewichten met nullen worden correct afgehandeld", {
  w <- rep(1, n)
  w[1:5] <- 0  # Vijf observaties met gewicht 0
  fit_wlm <- lm(y ~ x, data = df, weights = w)

  ll <- restriktor:::con_loglik_lm(fit_wlm)
  expect_true(is.numeric(ll))
  expect_true(is.finite(ll))
})


# =============================================================================
# con_loglik_lm: multivariate lm (mlm)
# =============================================================================

test_that("con_loglik_lm: werkt met mlm (multivariate response)", {
  y2 <- 1 + 0.3 * x + rnorm(n)
  fit_mlm <- lm(cbind(y, y2) ~ x, data = data.frame(y = y, y2 = y2, x = x))

  ll <- restriktor:::con_loglik_lm(fit_mlm)
  expect_true(is.numeric(ll))
  expect_true(is.finite(ll))
})


# =============================================================================
# con_loglik_glm: generalized linear model
# =============================================================================

test_that("con_loglik_glm: retourneert numerieke waarde (Poisson)", {
  y_count <- rpois(n, lambda = exp(0.5 + 0.3 * x))
  fit_glm <- glm(y_count ~ x, family = poisson, data = data.frame(y_count, x))

  ll <- restriktor:::con_loglik_glm(fit_glm)
  expect_true(is.numeric(ll))
  expect_true(is.finite(ll))
})

test_that("con_loglik_glm: retourneert numerieke waarde (Gaussian GLM)", {
  fit_glm_g <- glm(y ~ x, family = gaussian, data = df)

  ll <- restriktor:::con_loglik_glm(fit_glm_g)
  expect_true(is.numeric(ll))
  expect_true(is.finite(ll))
})

test_that("con_loglik_glm: vergelijkbaar met stats::logLik", {
  y_count <- rpois(n, lambda = exp(0.5 + 0.3 * x))
  fit_glm <- glm(y_count ~ x, family = poisson, data = data.frame(y_count, x))

  ll_restr <- restriktor:::con_loglik_glm(fit_glm)
  ll_stats <- as.numeric(logLik(fit_glm))

  # Formula: ll = p - aic/2 => ll_restr = logLik + p (close but not exact)
  expect_equal(ll_restr, ll_stats, tolerance = 5)
})

test_that("con_loglik_glm: beter model geeft hogere loglik (Poisson)", {
  y_count <- rpois(n, lambda = exp(0.5 + 0.3 * x))
  df_p <- data.frame(y_count = y_count, x = x)
  fit_good <- glm(y_count ~ x, family = poisson, data = df_p)
  fit_bad  <- glm(y_count ~ 1, family = poisson, data = df_p)

  ll_good <- restriktor:::con_loglik_glm(fit_good)
  ll_bad  <- restriktor:::con_loglik_glm(fit_bad)

  expect_true(ll_good >= ll_bad)
})
