# =============================================================================
# Tests: con_sandwich.R (sandwich, bread, estfun, meatHC)
#
# Twee exacte referenties:
#  1. inactieve restricties: alle HC-types moeten sandwich::vcovHC() exact
#     reproduceren (b.restr == b.unrestr, augmentatie is een no-op).
#  2. equality constraint x1 == x2: identiek aan het geherparametriseerde
#     model y ~ I(x1 + x2) + ..., dus V == J %*% vcovHC(reduced) %*% t(J).
# =============================================================================

skip_if_not_installed("sandwich")

set.seed(1234)
n_sw  <- 100
x1_sw <- rnorm(n_sw); x2_sw <- rnorm(n_sw); x3_sw <- rnorm(n_sw)
# heteroskedastische fouten zodat de HC-correcties er echt toe doen
e_sw  <- rnorm(n_sw, sd = 0.5 + 0.8 * abs(x1_sw))
y_sw  <- 1 + 0.5 * x1_sw + 0.6 * x2_sw - 0.3 * x3_sw + e_sw
df_sw <- data.frame(y = y_sw, x1 = x1_sw, x2 = x2_sw, x3 = x3_sw)

# J-matrix voor herparametrisatie b(x1) == b(x2): b_vol = J %*% b_gereduceerd
J_sw <- rbind(c(1, 0, 0),
              c(0, 1, 0),
              c(0, 1, 0),
              c(0, 0, 1))

hc_types_sw <- c("const", "HC0", "HC1", "HC2", "HC3", "HC4", "HC4m", "HC5")

# -----------------------------------------------------------------------------
# conLM
# -----------------------------------------------------------------------------

test_that("conLM: inactieve restricties reproduceren sandwich::vcovHC exact", {
  fit <- lm(y ~ x1 + x2 + x3, data = df_sw)
  for (tp in hc_types_sw) {
    z <- restriktor(fit, constraints = "x1 > -100; x2 > -100", se = tp,
                    mix_weights = "none")
    expect_equal(unname(summary(z)$V), unname(sandwich::vcovHC(fit, type = tp)),
                 tolerance = 1e-10, label = paste("conLM", tp))
  }
})

test_that("conLM: equality constraint == geherparametriseerd model", {
  fit     <- lm(y ~ x1 + x2 + x3, data = df_sw)
  fit_red <- lm(y ~ I(x1 + x2) + x3, data = df_sw)
  # HC2/HC3 dekken ook de hatvalues in de gereduceerde kolomruimte (X %*% N)
  for (tp in c("const", "HC0", "HC1", "HC2", "HC3")) {
    z <- restriktor(fit, constraints = "x1 == x2", se = tp, mix_weights = "none")
    expect_equal(unname(summary(z)$V),
                 unname(J_sw %*% sandwich::vcovHC(fit_red, type = tp) %*% t(J_sw)),
                 tolerance = 1e-10, label = paste("conLM eq", tp))
  }
  z <- restriktor(fit, constraints = "x1 == x2", se = "standard",
                  mix_weights = "none")
  expect_equal(unname(summary(z)$V),
               unname(J_sw %*% vcov(fit_red) %*% t(J_sw)),
               tolerance = 1e-10)
})

test_that("conLM: gewogen lm (WLS) geeft correcte standard en HC standaardfouten", {
  set.seed(99)
  w_sw <- runif(n_sw, 0.2, 3)
  fit  <- lm(y ~ x1 + x2 + x3, data = df_sw, weights = w_sw)
  z <- restriktor(fit, constraints = "x1 > -100", se = "standard",
                  mix_weights = "none")
  expect_equal(unname(summary(z)$coefficients[, "Std. Error"]),
               unname(summary(fit)$coefficients[, "Std. Error"]),
               tolerance = 1e-10)
  for (tp in c("HC0", "HC3")) {
    z <- restriktor(fit, constraints = "x1 > -100", se = tp,
                    mix_weights = "none")
    expect_equal(unname(summary(z)$V), unname(sandwich::vcovHC(fit, type = tp)),
                 tolerance = 1e-10, label = paste("conLM WLS", tp))
  }
})

test_that("conLM: nulgewicht geeft geen NaN in HC2/HC3 standaardfouten", {
  set.seed(99)
  w0 <- runif(n_sw, 0.2, 3); w0[1] <- 0
  fit <- lm(y ~ x1 + x2 + x3, data = df_sw, weights = w0)
  z <- restriktor(fit, constraints = "x1 > -100", se = "HC3",
                  mix_weights = "none")
  expect_true(all(is.finite(summary(z)$coefficients[, "Std. Error"])))
})

# -----------------------------------------------------------------------------
# conGLM
# -----------------------------------------------------------------------------

test_that("conGLM poisson: inactieve restricties reproduceren vcovHC exact", {
  set.seed(42)
  yp  <- rpois(n_sw, exp(0.3 + 0.4 * x1_sw + 0.5 * x2_sw - 0.2 * x3_sw))
  fit <- glm(yp ~ x1 + x2 + x3, data = df_sw, family = poisson)
  for (tp in hc_types_sw) {
    z <- restriktor(fit, constraints = "x1 > -100", se = tp,
                    mix_weights = "none")
    expect_equal(unname(summary(z)$V), unname(sandwich::vcovHC(fit, type = tp)),
                 tolerance = 1e-10, label = paste("conGLM poisson", tp))
  }
})

test_that("conGLM gaussian (dispersie != 1): inactief reproduceert vcovHC exact", {
  fit <- glm(y ~ x1 + x2 + x3, data = df_sw, family = gaussian)
  for (tp in c("const", "HC0", "HC3")) {
    z <- restriktor(fit, constraints = "x1 > -100", se = tp,
                    mix_weights = "none")
    expect_equal(unname(summary(z)$V), unname(sandwich::vcovHC(fit, type = tp)),
                 tolerance = 1e-10, label = paste("conGLM gaussian", tp))
  }
})

test_that("conGLM poisson: equality constraint == geherparametriseerd model", {
  set.seed(42)
  yp      <- rpois(n_sw, exp(0.3 + 0.4 * x1_sw + 0.5 * x2_sw - 0.2 * x3_sw))
  fit     <- glm(yp ~ x1 + x2 + x3, data = df_sw, family = poisson)
  fit_red <- glm(yp ~ I(x1 + x2) + x3, data = df_sw, family = poisson)
  z <- restriktor(fit, constraints = "x1 == x2", se = "HC0",
                  mix_weights = "none")
  # iteratieve fit: iets ruimere tolerantie
  expect_equal(unname(summary(z)$V),
               unname(J_sw %*% sandwich::vcovHC(fit_red, type = "HC0") %*% t(J_sw)),
               tolerance = 1e-4)
})

# -----------------------------------------------------------------------------
# conRLM
# -----------------------------------------------------------------------------

test_that("conRLM: inactieve restricties reproduceren vcovHC exact (const/HC0/HC1)", {
  set.seed(7)
  yr <- 1 + 0.5 * x1_sw + 0.6 * x2_sw - 0.3 * x3_sw + rnorm(n_sw)
  yr[1:5] <- yr[1:5] + 8  # outliers, door bisquare volledig neergewogen
  df_r <- cbind(df_sw, yr = yr)
  fit <- MASS::rlm(yr ~ x1 + x2 + x3, data = df_r, method = "MM")
  # "const" dekt de terugval op modelresiduen bij estfun-rijen die exact 0 zijn
  for (tp in c("const", "HC0", "HC1")) {
    z <- restriktor(fit, constraints = "x1 > -100", se = tp,
                    mix_weights = "none")
    expect_equal(unname(summary(z)$V), unname(sandwich::vcovHC(fit, type = tp)),
                 tolerance = 1e-8, label = paste("conRLM", tp))
  }
})

test_that("conRLM: HC2 t/m HC5 draaien zonder error (regressie: x$wgt in meatHC)", {
  set.seed(7)
  yr <- 1 + 0.5 * x1_sw + 0.6 * x2_sw - 0.3 * x3_sw + rnorm(n_sw)
  yr[1:5] <- yr[1:5] + 8
  df_r <- cbind(df_sw, yr = yr)
  fit <- MASS::rlm(yr ~ x1 + x2 + x3, data = df_r, method = "MM")
  for (tp in c("HC2", "HC3", "HC4", "HC4m", "HC5")) {
    z <- restriktor(fit, constraints = "x1 > -100", se = tp,
                    mix_weights = "none")
    s <- summary(z)
    expect_true(all(is.finite(s$coefficients[, "Std. Error"])),
                label = paste("conRLM", tp))
  }
})

# -----------------------------------------------------------------------------
# HAC: vcovHAC / kernHAC / NeweyWest (issue #8)
#
# Exacte referenties, analoog aan de HC-tests hierboven:
#  1. inactieve restricties: se = "HAC"/"kernHAC"/"NeweyWest" moet
#     sandwich::vcovHAC/kernHAC/NeweyWest op het originele model exact
#     reproduceren, en de sandwich-generics moeten ook direct op het
#     restriktor-object werken.
#  2. equality constraint x1 == x2: identiek aan het geherparametriseerde
#     model, mits de datagedreven stappen (bandbreedte, prewhitening) worden
#     vastgezet, want die zien 4 vs. 3 estfun-kolommen.
# -----------------------------------------------------------------------------

# AR(1) fouten zodat de autocorrelatie-correctie er echt toe doet
set.seed(321)
e_ar  <- as.numeric(stats::arima.sim(list(ar = 0.6), n = n_sw))
y_ar  <- 1 + 0.5 * x1_sw + 0.6 * x2_sw - 0.3 * x3_sw + e_ar
df_ar <- data.frame(y = y_ar, x1 = x1_sw, x2 = x2_sw, x3 = x3_sw)

test_that("conLM: inactieve restricties reproduceren vcovHAC/kernHAC/NeweyWest exact", {
  fit <- lm(y ~ x1 + x2 + x3, data = df_ar)
  ref <- list(HAC       = sandwich::vcovHAC(fit),
              kernHAC   = sandwich::kernHAC(fit),
              NeweyWest = sandwich::NeweyWest(fit))
  for (tp in names(ref)) {
    z <- restriktor(fit, constraints = "x1 > -100; x2 > -100", se = tp,
                    mix_weights = "none")
    expect_equal(unname(summary(z)$V), unname(ref[[tp]]),
                 tolerance = 1e-10, label = paste("conLM", tp))
  }
})

test_that("sandwich generics werken direct op restriktor-objecten (conLM)", {
  fit <- lm(y ~ x1 + x2 + x3, data = df_ar)
  z <- restriktor(fit, constraints = "x1 > -100", se = "standard",
                  mix_weights = "none")
  expect_equal(unname(sandwich::estfun(z)), unname(sandwich::estfun(fit)),
               tolerance = 1e-10)
  expect_equal(unname(sandwich::bread(z)), unname(sandwich::bread(fit)),
               tolerance = 1e-10)
  expect_equal(unname(sandwich::vcovHAC(z)), unname(sandwich::vcovHAC(fit)),
               tolerance = 1e-10)
  expect_equal(unname(sandwich::kernHAC(z)), unname(sandwich::kernHAC(fit)),
               tolerance = 1e-10)
  expect_equal(unname(sandwich::NeweyWest(z)), unname(sandwich::NeweyWest(fit)),
               tolerance = 1e-10)
})

test_that("conLM: equality constraint == geherparametriseerd model (HAC)", {
  fit     <- lm(y ~ x1 + x2 + x3, data = df_ar)
  fit_red <- lm(y ~ I(x1 + x2) + x3, data = df_ar)

  # kernHAC met vaste bandbreedte en zonder prewhitening
  z <- restriktor(fit, constraints = "x1 == x2", se = "kernHAC",
                  mix_weights = "none")
  expect_equal(unname(summary(z, bw = 2, prewhite = FALSE)$V),
               unname(J_sw %*% sandwich::kernHAC(fit_red, bw = 2, prewhite = FALSE) %*% t(J_sw)),
               tolerance = 1e-10)

  # vcovHAC met vaste gewichten
  w_hac <- c(1, 0.7, 0.4, 0.1)
  z <- restriktor(fit, constraints = "x1 == x2", se = "HAC",
                  mix_weights = "none")
  expect_equal(unname(summary(z, weights = w_hac)$V),
               unname(J_sw %*% sandwich::vcovHAC(fit_red, weights = w_hac) %*% t(J_sw)),
               tolerance = 1e-10)

  # NeweyWest met vaste lag en zonder prewhitening (adjust = FALSE default)
  z <- restriktor(fit, constraints = "x1 == x2", se = "NeweyWest",
                  mix_weights = "none")
  expect_equal(unname(summary(z, lag = 3, prewhite = FALSE)$V),
               unname(J_sw %*% sandwich::NeweyWest(fit_red, lag = 3, prewhite = FALSE) %*% t(J_sw)),
               tolerance = 1e-10)
})

test_that("conGLM poisson: inactieve restricties reproduceren vcovHAC exact", {
  set.seed(42)
  yp  <- rpois(n_sw, exp(0.3 + 0.4 * x1_sw + 0.5 * x2_sw - 0.2 * x3_sw))
  fit <- glm(yp ~ x1 + x2 + x3, data = df_sw, family = poisson)
  z <- restriktor(fit, constraints = "x1 > -100", se = "HAC",
                  mix_weights = "none")
  expect_equal(unname(summary(z)$V), unname(sandwich::vcovHAC(fit)),
               tolerance = 1e-10)
  z <- restriktor(fit, constraints = "x1 > -100", se = "NeweyWest",
                  mix_weights = "none")
  expect_equal(unname(summary(z)$V), unname(sandwich::NeweyWest(fit)),
               tolerance = 1e-10)
})

test_that("conRLM: inactieve restricties reproduceren vcovHAC exact", {
  set.seed(7)
  yr <- 1 + 0.5 * x1_sw + 0.6 * x2_sw - 0.3 * x3_sw + rnorm(n_sw)
  yr[1:5] <- yr[1:5] + 8
  df_r <- cbind(df_sw, yr = yr)
  fit <- MASS::rlm(yr ~ x1 + x2 + x3, data = df_r, method = "MM")
  z <- restriktor(fit, constraints = "x1 > -100", se = "HAC",
                  mix_weights = "none")
  expect_equal(unname(summary(z)$V), unname(sandwich::vcovHAC(fit)),
               tolerance = 1e-8)
})

test_that("conLM: equality constraint geeft gelijke HAC standaardfouten voor x1 en x2", {
  fit <- lm(y ~ x1 + x2 + x3, data = df_ar)
  for (tp in c("HAC", "kernHAC", "NeweyWest")) {
    z <- restriktor(fit, constraints = "x1 == x2", se = tp,
                    mix_weights = "none")
    s <- summary(z)
    expect_true(all(is.finite(s$coefficients[, "Std. Error"])),
                label = paste("conLM eq", tp))
    expect_equal(unname(s$coefficients["x1", "Std. Error"]),
                 unname(s$coefficients["x2", "Std. Error"]),
                 tolerance = 1e-10, label = paste("conLM eq", tp))
  }
})

test_that("conRLM: equality constraint geeft eindige HC standaardfouten", {
  set.seed(7)
  yr <- 1 + 0.5 * x1_sw + 0.6 * x2_sw - 0.3 * x3_sw + rnorm(n_sw)
  yr[1:5] <- yr[1:5] + 8
  df_r <- cbind(df_sw, yr = yr)
  fit <- MASS::rlm(yr ~ x1 + x2 + x3, data = df_r, method = "MM")
  for (tp in c("HC0", "HC3")) {
    z <- restriktor(fit, constraints = "x1 == x2", se = tp,
                    mix_weights = "none")
    s <- summary(z)
    expect_true(all(is.finite(s$coefficients[, "Std. Error"])))
    # onder x1 == x2 zijn de standaardfouten van x1 en x2 gelijk
    expect_equal(unname(s$coefficients["x1", "Std. Error"]),
                 unname(s$coefficients["x2", "Std. Error"]), tolerance = 1e-10)
  }
})
