# =============================================================================
# Tests for fixes of former TO DOs in the code (second round)
# =============================================================================

# --- leave1studyout(): order of the output, N_tot, and note ---------------------

test_that("leave1studyout(): IC weights first, N_tot (if available), and a note on the IC values", {
  set.seed(1)
  est <- lapply(1:3, function(i) c(b1 = rnorm(1, .3, .05), b2 = rnorm(1, 0, .05)))
  vc  <- lapply(1:3, function(i) diag(c(.01, .01)))
  r   <- suppressMessages(evSyn(est, VCOV = vc, hypotheses = list(H1 = "b1 > b2"),
                                study_sample_nobs = c(50, 100, 150)))
  l1  <- leave1studyout(r)
  expect_identical(names(l1)[1:3], c("OverallGoricaWeights", "OverallPrefHypo", "N_tot"))
  expect_equal(as.vector(l1$N_tot), c(250, 200, 150))
  out <- capture.output(print(l1))
  expect_true(grep("weights", out)[1] < grep("values", out)[1])
  expect_true(any(grepl("N_tot", out)))
  expect_true(any(grepl("cannot be compared across the", out)))
  # without sample sizes: no N_tot
  r2 <- suppressMessages(evSyn(est, VCOV = vc, hypotheses = list(H1 = "b1 > b2")))
  expect_null(leave1studyout(r2)$N_tot)
})


# --- lavaan: hypotheses with lavaan's own parameter names ------------------------

test_that("goric() for lavaan: lavaan names (with/without spaces) and restriktor names give the same result as labels", {
  set.seed(3)
  n <- 200
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- .4 * d$x1 + .2 * d$x2 + rnorm(n)
  fit  <- lavaan::sem("y ~ x1 + x2 + x3", data = d)
  fitl <- lavaan::sem("y ~ a*x1 + b*x2 + c*x3", data = d)
  r_lab <- suppressMessages(goric(fitl, hypotheses = list(H1 = "a > b > c")))
  r1 <- suppressMessages(goric(fit, hypotheses = list(H1 = "y~x1 > y~x2 > y~x3")))
  r2 <- suppressMessages(goric(fit, hypotheses = list(H1 = "y ~ x1 > y ~ x2 > y ~ x3")))
  r3 <- suppressMessages(goric(fit, hypotheses = list(H1 = "y.on.x1 > y.on.x2 > y.on.x3")))
  for (r in list(r1, r2, r3)) expect_equal(r$result[, 2:4], r_lab$result[, 2:4])
  # goricac: same number of nuisance parameters as with labels
  ra  <- suppressMessages(goric(fit,  hypotheses = list(H1 = "y~x1 > y~x2"), type = "goricac"))
  ral <- suppressMessages(goric(fitl, hypotheses = list(H1 = "a > b"), type = "goricac"))
  expect_equal(ra$result$penalty, ral$result$penalty)
  # unknown name: informative error with the available names
  expect_error(suppressMessages(goric(fit, hypotheses = list(H1 = "y~x1 > y~x9"))),
               "available names are")
})

test_that("lavaan names for multiple groups and latent variables", {
  set.seed(3)
  d <- data.frame(x1 = rnorm(200), x2 = rnorm(200), g = rep(1:2, each = 100))
  d$y <- .3 * d$x1 + rnorm(200)
  fitg <- lavaan::sem("y ~ x1 + x2", data = d, group = "g")
  est <- restriktor:::con_gorica_est_lav(fitg)
  expect_true(all(c("y.on.x1", "y.on.x1.g2") %in% names(est$estimate)))
  expect_s3_class(suppressMessages(goric(fitg, hypotheses = list(H1 = "y~x1 > y~x1.g2"))),
                  "con_goric")
  fcfa <- lavaan::cfa("visual =~ x1 + x2 + x3", data = lavaan::HolzingerSwineford1939)
  expect_s3_class(suppressMessages(goric(fcfa, hypotheses = list(H1 = "visual=~x2 > visual=~x3"))),
                  "con_goric")
})


# --- benchmark: group_size is only for model_type = "means" ------------------------

test_that("benchmark(model_type = 'asymp') gives an error for group_size", {
  fit <- goric(c(a = .3, b = .1), VCOV = diag(2) * .01, hypotheses = list(H1 = "a > b"))
  expect_error(benchmark(fit, model_type = "asymp", group_size = 50, iter = 10),
               "only used for model_type")
})


# --- complement: log-likelihood based on the original constraints ---------------

test_that("range restrictions: the complement's penalty is lower than the unconstrained one", {
  est <- c(b1 = .3, b2 = .1)
  fc <- goric(est, VCOV = diag(2) * .01, hypotheses = list(H1 = "-0.2 < b1 < 0.5; b2 > 0"))
  fu <- goric(est, VCOV = diag(2) * .01, hypotheses = list(H1 = "-0.2 < b1 < 0.5; b2 > 0"),
              comparison = "unconstrained")
  expect_lt(fc$result$penalty[2], fu$result$penalty[2])
  # the complement's log-likelihood is the maximum outside the region of H1:
  # here at b2 = 0 (b1 stays at its estimate), so loglik(H1) - 0.5 * .1^2 / .01
  expect_equal(fc$result$loglik[2], fc$result$loglik[1] - 0.5, tolerance = 1e-6)
})


# --- evSyn(): the same order of the output elements in all evSyn functions -------

test_that("evSyn output: settings, studies, per study, cumulative, final (values before weights)", {
  set.seed(1)
  est <- lapply(1:3, function(i) c(b1 = rnorm(1, .3, .05), b2 = rnorm(1, 0, .05)))
  vc  <- lapply(1:3, function(i) diag(c(.01, .01)))
  r <- suppressMessages(evSyn(est, VCOV = vc, hypotheses = list(H1 = "b1 > b2")))
  nm <- names(r)
  pos <- function(x) match(x, nm)
  expect_true(pos("type") < pos("n_studies"))
  expect_true(pos("n_studies") < pos("LL_m"))
  expect_true(pos("GORICA_m") < pos("LL_weights_m"))
  expect_true(pos("GORICA_weight_m") < pos("Cumulative_LL"))
  expect_true(pos("Cumulative_GORICA") < pos("Cumulative_LL_weights"))
  expect_true(pos("Cumulative_GORICA_weights") < pos("Final_ratio_LL_weights"))
  # the same relative order for other input types (that have fewer elements)
  W <- lapply(asplit(r$GORICA_weight_m, 1), c)
  rw <- suppressMessages(evSyn(W))
  common <- intersect(nm, names(rw))
  expect_identical(names(rw)[names(rw) %in% common], nm[nm %in% common])
})


# --- 'sd' in the hypotheses (effect size, like Cohen's d) ---------------------------

test_that("'sd' in the hypotheses is replaced by the residual SD of the model", {
  set.seed(1)
  d <- data.frame(y = rnorm(90, rep(c(0, .3, .6), each = 30)), g = gl(3, 30))
  m <- lm(y ~ g - 1, data = d)
  s <- sigma(m)
  fit_sd  <- suppressMessages(goric(m, hypotheses = list(H1 = "g2 - g1 > 0.2*sd; g3 - g2 > 0.2*sd")))
  fit_num <- goric(m, hypotheses = list(H1 = paste0("g2 - g1 > ", 0.2 * s, "; g3 - g2 > ", 0.2 * s)))
  expect_equal(fit_sd$result, fit_num$result)
  expect_message(goric(m, hypotheses = list(H1 = "g2 - g1 > 0.2*sd")), "residual standard deviation")
  # stored hypotheses contain the value (such that, e.g., benchmark() can use them)
  expect_false(grepl("sd", fit_sd$hypotheses_usr$H1))
  # mlm: one per outcome
  mm <- lm(cbind(y1 = y, y2 = y + rnorm(90)) ~ g - 1, data = d)
  expect_s3_class(suppressMessages(goric(mm, hypotheses = list(H1 = "y1.g2 - y1.g1 > 0.5*sd.y1"))),
                  "con_goric")
  expect_error(suppressMessages(goric(mm, hypotheses = list(H1 = "y1.g2 - y1.g1 > 0.5*sd"))),
               "sd.y1")
  # estimates as input: not available
  expect_error(goric(c(a = 1, b = .5), VCOV = diag(2), hypotheses = list(H1 = "a - b > 0.2*sd")),
               "fitted model object")
  # glm (non-gaussian): not available
  dg <- data.frame(y = rbinom(90, 1, .5), g = gl(3, 30))
  mg <- glm(y ~ g - 1, family = binomial, data = dg)
  expect_error(suppressMessages(goric(mg, hypotheses = list(H1 = "g2 - g1 > 0.2*sd"))), "binomial")
})


# --- 'sd' must be used as an effect-size cut-off ------------------------------------

test_that("'sd' in the hypotheses: allowed and not allowed forms", {
  set.seed(1)
  d <- data.frame(y = rnorm(90, rep(c(0, .3, .6), each = 30)), g = gl(3, 30))
  m <- lm(y ~ g - 1, data = d)
  ok <- c("g1 > sd", "g2 - g1 > 0.2*sd + g1", "(g2 - g1)/sd > 0.2", "g2 - g1 > sd/2",
          "g1 < g2 < g3; g3 - g1 > 0.8*sd")
  for (h in ok) expect_s3_class(suppressMessages(goric(m, hypotheses = list(H1 = h))), "con_goric")
  expect_error(suppressMessages(goric(m, hypotheses = list(H1 = "sd > 2"))), "does not contain any parameter")
  for (h in c("g2 - g1 > sd*g3", "g2 > g1*sd", "sd*(g2 - g1) > 1", "(g2 - g1)*sd > 1", "g2 - g1 > sd/g3")) {
    expect_error(suppressMessages(goric(m, hypotheses = list(H1 = h))), "multiplies")
  }
  # '(g2 - g1)/sd > 0.2' is the same as 'g2 - g1 > 0.2*sd'
  f1 <- suppressMessages(goric(m, hypotheses = list(H1 = "(g2 - g1)/sd > 0.2")))
  f2 <- suppressMessages(goric(m, hypotheses = list(H1 = "g2 - g1 > 0.2*sd")))
  expect_equal(f1$result, f2$result)
})


# --- evSyn_CI(): joint case bootstrap (type_ev = "average") --------------------------

test_that("evSyn_CI(): joint case bootstrap with data and fit_fun", {
  set.seed(10)
  n <- 120
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  dat$y <- .4 * dat$x1 + .2 * dat$x2 + .1 * dat$x3 + rnorm(n)
  fit_fun <- function(d) {
    m1 <- lm(y ~ x1 + x2, data = d)
    m2 <- lm(y ~ x1 + x2 + x3, data = d)
    lapply(list(m1, m2), function(m)
      list(est = coef(m)[c("x1", "x2")], VCOV = vcov(m)[c("x1", "x2"), c("x1", "x2")]))
  }
  fits <- fit_fun(dat)
  ev <- suppressMessages(evSyn(lapply(fits, `[[`, "est"), VCOV = lapply(fits, `[[`, "VCOV"),
                               hypotheses = list(H1 = "x1 > x2 > 0"), type_ev = "average"))
  ci <- evSyn_CI(ev, data = dat, fit_fun = fit_fun, B = 50, seed = 1)
  expect_true(ci$joint)
  expect_identical(ci$n_failed, 0L)
  expect_true(all(ci$weights$lower <= ci$weights[[2]] + 1e-12))
  expect_true(all(ci$weights$upper >= ci$weights[[2]] - 1e-12))
  expect_output(print(ci), "joint case bootstrap")
  # the final values equal those of evSyn
  expect_equal(unname(ci$IC[[2]]), unname(ev$Cumulative_GORICA["Final", ]), tolerance = 1e-8)
  # checks
  expect_error(evSyn_CI(ev, data = dat), "both 'data'")
  ev_add <- suppressMessages(evSyn(lapply(fits, `[[`, "est"), VCOV = lapply(fits, `[[`, "VCOV"),
                                   hypotheses = list(H1 = "x1 > x2 > 0")))
  expect_error(evSyn_CI(ev_add, data = dat, fit_fun = fit_fun), "type_ev = 'average'")
  bad_fun <- function(d) fit_fun(d)[1]
  expect_error(evSyn_CI(ev, data = dat, fit_fun = bad_fun, B = 5), "one element per analysis")
  # fit_fun that does not reproduce the estimates: warning
  other_fun <- function(d) lapply(fit_fun(d), function(f) { f$est <- f$est + 1; f })
  expect_warning(evSyn_CI(ev, data = dat, fit_fun = other_fun, B = 5, seed = 1), "does not render")
})


# --- evSyn(): priorICweights_study ----------------------------------------------------

test_that("priorICweights_study: study-specific weights include the prior; final weights use it once", {
  set.seed(1)
  est <- lapply(1:4, function(i) c(b1 = rnorm(1, .3, .1), b2 = rnorm(1, .1, .1)))
  vc  <- lapply(1:4, function(i) diag(c(.01, .01)))
  H <- list(H1 = "b1 > b2", H2 = "b1 < b2")
  pr <- c(.6, .2, .2)
  r0 <- suppressMessages(evSyn(est, VCOV = vc, hypotheses = H, priorICweights = pr))
  r1 <- suppressMessages(evSyn(est, VCOV = vc, hypotheses = H, priorICweights = pr,
                               priorICweights_study = TRUE))
  expect_false(r0$priorICweights_study)
  expect_true(r1$priorICweights_study)
  # final (and cumulative) weights: the same
  expect_equal(r0$Cumulative_GORICA_weights, r1$Cumulative_GORICA_weights)
  # study-specific weights: with the prior
  chk <- sweep(r0$GORICA_weight_m, 2, pr, "*")
  expect_equal(unname(r1$GORICA_weight_m), unname(chk / rowSums(chk)))
  # other input types
  W <- lapply(asplit(r0$GORICA_weight_m, 1), c)
  a <- suppressMessages(evSyn(W, priorICweights = pr))
  b <- suppressMessages(evSyn(W, priorICweights = pr, priorICweights_study = TRUE))
  expect_equal(a$Cumulative_GORICA_weights["Final", ], b$Cumulative_GORICA_weights["Final", ])
  # print states the choice; plot and helpers work
  expect_true(any(grepl("priorICweights_study = TRUE", capture.output(print(r1)))))
  expect_s3_class(plot(r1), "ggplot")
  expect_s3_class(leave1studyout(r1), "leave1studyout.evSyn")
  # cluster: the choice is kept
  rc <- suppressMessages(evSyn(est, VCOV = vc, hypotheses = H, priorICweights = pr,
                               priorICweights_study = TRUE, cluster = c(1, 1, 2, 3)))
  expect_true(rc$priorICweights_study)
  expect_error(evSyn(est, VCOV = vc, hypotheses = H, priorICweights_study = "yes"), "TRUE or FALSE")
})

test_that("evSyn() with goric objects uses the names of the hypotheses of the first object", {
  est <- list(c(b1 = .3, b2 = .1), c(b1 = .2, b2 = .15), c(b1 = .25, b2 = .05))
  gl <- lapply(est, function(e) goric(e, VCOV = diag(2) * .01,
                                      hypotheses = list(H1 = "b1 > b2", H2 = "b1 < b2")))
  g <- suppressMessages(evSyn(gl))
  expect_identical(colnames(g$GORICA_weight_m), c("H1", "H2", "unconstrained"))
  g2 <- suppressMessages(evSyn(gl, hypo_names = c("A", "B", "C")))
  expect_identical(colnames(g2$GORICA_weight_m), c("A", "B", "C"))
})
