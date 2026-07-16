# =============================================================================
# Tests: restriktor() / conLM / conGLM / conRLM / conMLM
# =============================================================================
# Zorg dat vereiste imports beschikbaar zijn in het huidige library-pad
# for (pkg in c("lavaan", "mvtnorm", "tmvtnorm", "quadprog", "norm")) {
#   if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# }

# # Laad het restriktor package vanuit broncode
# setwd("/Workspace/Users/l.vanbrabant@ggdwestbrabant.nl/restriktor")
# devtools::load_all(".")

# --- Gedeelde testdata ---

set.seed(42)
n <- 60
group <- factor(rep(1:3, each = n / 3))
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 2 + 1.5 * (group == "2") + 3 * (group == "3") + 0.5 * x1 + rnorm(n)

df_test <- data.frame(y = y, group = group, x1 = x1, x2 = x2)

# --- lm fits ---
fit_lm   <- lm(y ~ -1 + group, data = df_test)
fit_lm2  <- lm(y ~ x1 + x2, data = df_test)

# --- glm fit (Poisson) ---
y_count <- rpois(n, lambda = exp(0.5 + 0.3 * x1))
fit_glm <- glm(y_count ~ x1 + x2, family = poisson, data = df_test)

# --- rlm fit ---
fit_rlm <- MASS::rlm(y ~ -1 + group, data = df_test, method = "MM")

# --- mlm fit ---
y2 <- 1 + 0.8 * (group == "2") + 2 * (group == "3") + rnorm(n)
Y_mat <- cbind(y, y2)
fit_mlm <- lm(Y_mat ~ -1 + group, data = df_test)


# =============================================================================
# conLM.lm: basisgedrag
# =============================================================================

test_that("conLM.lm: string-restrictie geeft restriktor object", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3")
  expect_s3_class(result, "restriktor")
  expect_true("conLM" %in% class(result))
})

test_that("conLM.lm: beperkte schattingen respecteren ongelijkheden", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3")
  b <- coef(result)
  expect_true(b["group1"] <= b["group2"] + 1e-6)
  expect_true(b["group2"] <= b["group3"] + 1e-6)
})

test_that("conLM.lm: gelijkheidsbeperking wordt gerespecteerd", {
  result <- restriktor(fit_lm, constraints = "group1 = group2")
  b <- coef(result)
  expect_equal(unname(b["group1"]), unname(b["group2"]), tolerance = 1e-5)
})

test_that("conLM.lm: matrixrestricties werken", {
  Amat <- rbind(c(-1, 1, 0), c(0, -1, 1))
  result <- restriktor(fit_lm, constraints = Amat)
  b <- coef(result)
  expect_true(b["group1"] <= b["group2"] + 1e-6)
  expect_true(b["group2"] <= b["group3"] + 1e-6)
})

test_that("conLM.lm: matrixrestrictie met rhs en neq", {
  Amat <- rbind(c(1, -1, 0), c(0, 1, -1))
  result <- restriktor(fit_lm, constraints = Amat, rhs = c(0, 0), neq = 1)
  b <- coef(result)
  # Eerste restrictie is gelijkheid: group1 = group2
  expect_equal(unname(b["group1"]), unname(b["group2"]), tolerance = 1e-4)
})


# =============================================================================
# conLM.lm: output structuur
# =============================================================================

test_that("conLM.lm: output bevat alle verwachte velden", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3")

  expect_true(!is.null(result$b.restr))
  expect_true(!is.null(result$b.unrestr))
  expect_true(!is.null(result$constraints))
  expect_true(!is.null(result$neq))
  expect_true(!is.null(result$rhs))
  expect_true(!is.null(result$wt.bar))
  expect_true(!is.null(result$model.org))
  expect_true(!is.null(result$residuals))
  expect_true(!is.null(result$fitted))
})

test_that("conLM.lm: wt.bar gewichten tellen op tot 1", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3")
  expect_equal(sum(result$wt.bar), 1, tolerance = 1e-4)
})

test_that("conLM.lm: logLik methode werkt", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3")
  ll <- logLik(result)
  expect_true(is.numeric(ll))
  expect_true(is.finite(ll))
})

test_that("conLM.lm: coef methode geeft benoemde vector", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3")
  b <- coef(result)
  expect_true(is.numeric(b))
  expect_true(!is.null(names(b)))
  expect_equal(length(b), length(coef(fit_lm)))
})

test_that("conLM.lm: model.matrix methode werkt", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3")
  mm <- model.matrix(result)
  expect_true(is.matrix(mm))
  expect_equal(nrow(mm), nrow(model.matrix(fit_lm)))
})


# =============================================================================
# conLM.lm: summary
# =============================================================================

test_that("conLM.lm: summary geeft verwacht object", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3")
  s <- summary(result)
  expect_true(!is.null(s$coefficients))
  expect_true(is.matrix(s$coefficients))
  expect_equal(nrow(s$coefficients), length(coef(fit_lm)))
})

test_that("conLM.lm: summary met goric = 'none' werkt", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3")
  s <- summary(result, goric = "none")
  expect_true(!is.null(s$coefficients))
})


# =============================================================================
# conLM.lm: se-opties
# =============================================================================

test_that("conLM.lm: se = 'none' geeft geen standaardfouten", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3",
                       se = "none")
  s <- summary(result)
  # Bij se = none: geen standard errors
  expect_true(is.null(s$coefficients) || ncol(s$coefficients) <= 1 ||
              all(is.na(s$coefficients[, "Std. Error"])) ||
              !("Std. Error" %in% colnames(s$coefficients)))
})

test_that("conLM.lm: se = 'standard' geeft standaardfouten", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3",
                       se = "standard")
  s <- summary(result)
  expect_true("Std. Error" %in% colnames(s$coefficients))
  expect_true(all(s$coefficients[, "Std. Error"] >= 0))
})

test_that("conLM.lm: HC standaardfouten werken", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3",
                       se = "HC0")
  s <- summary(result)
  expect_true("Std. Error" %in% colnames(s$coefficients))
})


# =============================================================================
# conLM.lm: mix_weights opties
# =============================================================================

test_that("conLM.lm: mix_weights = 'pmvnorm' werkt", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3",
                       mix_weights = "pmvnorm")
  expect_equal(attr(result$wt.bar, "method"), "pmvnorm")
})

test_that("conLM.lm: mix_weights = 'boot' werkt", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3",
                       mix_weights = "boot", seed = 42)
  expect_equal(attr(result$wt.bar, "method"), "boot")
})

test_that("conLM.lm: mix_weights = 'none' slaat gewichten over", {
  result <- restriktor(fit_lm, constraints = "group1 < group2 < group3",
                       mix_weights = "none")
  # Bij 'none': wt.bar is NA of NULL
  expect_true(is.null(result$wt.bar) || all(is.na(result$wt.bar)))
})


# =============================================================================
# conLM.lm: regressiemodel (met predictors, niet alleen groepsgemiddelden)
# =============================================================================

test_that("conLM.lm: restricties op regressiecoefficienten", {
  result <- restriktor(fit_lm2, constraints = "x1 > 0; x2 > 0")
  b <- coef(result)
  expect_true(b["x1"] >= -1e-6)
  expect_true(b["x2"] >= -1e-6)
})

test_that("conLM.lm: gelijkheid op regressiecoefficienten", {
  result <- restriktor(fit_lm2, constraints = "x1 = x2")
  b <- coef(result)
  expect_equal(unname(b["x1"]), unname(b["x2"]), tolerance = 1e-6)
})


# =============================================================================
# conGLM.glm
# =============================================================================

test_that("conGLM.glm: Poisson met ongelijkheidsrestricties", {
  result <- restriktor(fit_glm, constraints = "x1 > 0")
  expect_s3_class(result, "restriktor")
  expect_true("conGLM" %in% class(result))
  b <- coef(result)
  expect_true(b["x1"] >= -1e-6)
})

test_that("conGLM.glm: gelijkheidsbeperking", {
  result <- restriktor(fit_glm, constraints = "x1 = x2")
  b <- coef(result)
  expect_equal(unname(b["x1"]), unname(b["x2"]), tolerance = 1e-4)
})

test_that("conGLM.glm: output heeft restriktor klasse", {
  result <- restriktor(fit_glm, constraints = "x1 > 0")
  expect_true(inherits(result, "restriktor"))
  expect_true(!is.null(result$b.restr))
  expect_true(!is.null(result$wt.bar))
})

test_that("conGLM.glm: summary werkt", {
  result <- restriktor(fit_glm, constraints = "x1 > 0")
  s <- summary(result)
  expect_true(!is.null(s$coefficients))
})


# =============================================================================
# conRLM.rlm
# =============================================================================

test_that("conRLM.rlm: robuuste lm met ongelijkheidsrestricties", {
  result <- restriktor(fit_rlm, constraints = "group1 < group2 < group3")
  expect_s3_class(result, "restriktor")
  expect_true("conRLM" %in% class(result))
  b <- coef(result)
  expect_true(b["group1"] <= b["group2"] + 1e-6)
  expect_true(b["group2"] <= b["group3"] + 1e-6)
})

test_that("conRLM.rlm: gelijkheidsbeperking", {
  result <- restriktor(fit_rlm, constraints = "group1 = group2")
  b <- coef(result)
  expect_equal(unname(b["group1"]), unname(b["group2"]), tolerance = 1e-4)
})

test_that("conRLM.rlm: output bevat verwachte velden", {
  result <- restriktor(fit_rlm, constraints = "group1 < group2 < group3")
  expect_true(!is.null(result$b.restr))
  expect_true(!is.null(result$b.unrestr))
  expect_true(!is.null(result$constraints))
  expect_true(!is.null(result$wt.bar))
})



# =============================================================================
# conMLM.mlm
# =============================================================================
# Coefficientvolgorde vec(B): eerst de 3 coefficienten van y, dan die van y2.

test_that("conMLM.mlm: matrixrestrictie geeft restriktor object en respecteert restricties", {
  # geschonden restrictie: group1 (y) >= group3 (y)
  Amat <- rbind(c(1, 0, -1, 0, 0, 0))
  result <- restriktor(fit_mlm, constraints = Amat, mix_weights = "none")
  expect_s3_class(result, "restriktor")
  expect_true("conMLM" %in% class(result))
  b <- result$b.restr
  expect_true(b["group1", "y"] >= b["group3", "y"] - 1e-6)
})

test_that("conMLM.mlm: gelijkheidsbeperking wordt gerespecteerd", {
  Amat <- rbind(c(1, -1, 0, 0, 0, 0))
  result <- restriktor(fit_mlm, constraints = Amat, neq = 1, mix_weights = "none")
  b <- result$b.restr
  expect_equal(unname(b["group1", "y"]), unname(b["group2", "y"]), tolerance = 1e-5)
})

test_that("conMLM.mlm: R2.reduced komt overeen met handmatige berekening per respons", {
  Amat <- rbind(c(1, 0, -1, 0, 0, 0))
  result <- restriktor(fit_mlm, constraints = Amat, mix_weights = "none")
  fitted <- result$fitted
  resid  <- result$residuals
  # fit_mlm heeft geen intercept, dus mss = colSums(fitted^2)
  mss <- colSums(fitted^2)
  rss <- colSums(resid^2)
  expect_equal(unname(result$R2.reduced), unname(mss / (mss + rss)), tolerance = 1e-8)
})

test_that("conMLM.mlm: R2.reduced correct bij model met intercept en responsen op verschillende schaal", {
  # regressietest voor de recycling-bug in de gecentreerde kwadratensom
  y2_scaled <- 100 + df_test$y + rnorm(n)
  fit_mlm_int <- lm(cbind(y, y2_scaled) ~ x1, data = cbind(df_test, y2_scaled = y2_scaled))
  # geschonden restrictie op de slope van y: b_x1 (y) >= 1
  Amat <- rbind(c(0, 1, 0, 0))
  result <- restriktor(fit_mlm_int, constraints = Amat, rhs = 1, mix_weights = "none")
  fitted <- result$fitted
  resid  <- result$residuals
  mss <- colSums(sweep(fitted, 2, colMeans(fitted))^2)
  rss <- colSums(resid^2)
  expect_equal(unname(result$R2.reduced), unname(mss / (mss + rss)), tolerance = 1e-8)
  expect_true(all(result$R2.reduced < 1))
})

test_that("conMLM.mlm: loglik neemt niet toe door restricties", {
  Amat <- rbind(c(1, 0, -1, 0, 0, 0))
  result <- restriktor(fit_mlm, constraints = Amat, mix_weights = "none")
  ll.unrestr <- restriktor(fit_mlm, mix_weights = "none")$loglik
  expect_true(result$loglik <= ll.unrestr + 1e-8)
})

test_that("conMLM.mlm: se = 'standard' zonder restricties reproduceert vcov(mlm)", {
  result <- restriktor(fit_mlm, se = "standard", mix_weights = "none")
  s <- summary(result)
  se.restr <- s$coefficients[, "Std. Error"]
  se.vcov  <- sqrt(diag(vcov(fit_mlm)))
  expect_equal(unname(se.restr), unname(se.vcov), tolerance = 1e-6)
})

test_that("conMLM.mlm: se = 'standard' met actieve restrictie geeft bruikbare summary", {
  Amat <- rbind(c(1, 0, -1, 0, 0, 0))
  result <- restriktor(fit_mlm, constraints = Amat, se = "standard",
                       mix_weights = "none")
  s <- summary(result)
  coefs <- s$coefficients
  expect_equal(nrow(coefs), length(coef(fit_mlm)))
  expect_true(all(c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
                  %in% colnames(coefs)))
  expect_true(all(is.finite(coefs[, "Std. Error"])))
})

test_that("conMLM.mlm: summary met se = 'none' behoudt coefficientenmatrix", {
  Amat <- rbind(c(1, 0, -1, 0, 0, 0))
  result <- restriktor(fit_mlm, constraints = Amat, mix_weights = "none")
  s <- summary(result)
  expect_equal(ncol(s$coefficients), ncol(coef(fit_mlm)))
})

test_that("conMLM.mlm: fout bij niet-ondersteunde se methode", {
  expect_error(
    restriktor(fit_mlm, constraints = rbind(c(1, 0, -1, 0, 0, 0)),
               se = "boot.standard"),
    "not \\(yet\\) available"
  )
})

test_that("conMLM.mlm: fout bij gewogen mlm fit", {
  w <- runif(n, 0.5, 1.5)
  fit_mlm_w <- lm(Y_mat ~ -1 + group, data = df_test, weights = w)
  expect_error(
    restriktor(fit_mlm_w, constraints = rbind(c(1, 0, -1, 0, 0, 0))),
    "weighted mlm"
  )
})


# =============================================================================
# Foutafhandeling
# =============================================================================

test_that("conLM fout bij verkeerde object klasse", {
  expect_error(conLM.lm(list()), "class lm")
})

test_that("conGLM fout bij verkeerde object klasse", {
  expect_error(conGLM.glm(fit_lm), "class glm")
})

test_that("conRLM fout bij verkeerde object klasse", {
  expect_error(conRLM.rlm(fit_lm), "class rlm")
})

test_that("conMLM fout bij verkeerde object klasse", {
  expect_error(conMLM.mlm(fit_lm), "class mlm")
})

test_that("restriktor fout bij ongeldige se methode", {
  expect_error(
    restriktor(fit_lm, constraints = "group1 < group2", se = "onbekend"),
    "unknown"
  )
})

test_that("restriktor fout bij ongeldige mix_weights", {
  expect_error(
    restriktor(fit_lm, constraints = "group1 < group2", mix_weights = "xyz"),
    "unknown"
  )
})
