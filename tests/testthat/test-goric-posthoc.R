# =============================================================================
# Tests: post-hoc analysis (best hypothesis vs. its complement)
#        goric_posthoc() and goric(..., posthoc = TRUE)
# =============================================================================

est_ph  <- c(x1 = 0.5, x2 = 0.3, x3 = 0.1)
VCOV_ph <- diag(c(0.01, 0.01, 0.01))
H_ph    <- list(H1 = "x1 > x2 > x3", H2 = "x1 < x2 < x3")

set.seed(123)
dat_ph <- data.frame(y = rnorm(90, rep(c(0, 0.5, 1), each = 30)), g = gl(3, 30))
fit_lm_ph <- lm(y ~ g - 1, data = dat_ph)


# --- argument posthoc ----------------------------------------------------------

test_that("default (posthoc = TRUE): performed with several hypotheses when one of these is the best", {
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph)
  expect_s3_class(fit$posthoc, "goric_posthoc")
  expect_true(fit$posthoc$applicable)
  expect_true(any(grepl("Post-hoc analysis", capture.output(print(fit)))))
})

test_that("default: skipped without a message for a single hypothesis", {
  # vs. its complement (default comparison) and vs. the unconstrained
  expect_no_message(fit1 <- goric(est_ph, VCOV = VCOV_ph, hypotheses = list(H1 = "x1 > x2 > x3")))
  expect_null(fit1$posthoc)
  fit2 <- goric(est_ph, VCOV = VCOV_ph, hypotheses = list(H1 = "x1 > x2 > x3"),
                comparison = "unconstrained")
  expect_null(fit2$posthoc)
  expect_false(any(grepl("Post-hoc analysis", capture.output(print(fit2)))))
})

test_that("default: skipped without a message when the unconstrained or Heq is the best", {
  est_u <- c(x1 = 0.5, x2 = 0.0, x3 = 0.5)
  expect_no_message(fit <- goric(est_u, VCOV = VCOV_ph / 10, hypotheses = H_ph))
  expect_identical(as.character(fit$result$model[which.max(fit$result[, 7])]), "unconstrained")
  expect_null(fit$posthoc)
  expect_false(any(grepl("Post-hoc analysis", capture.output(print(fit)))))
  fit_heq <- goric(est_ph, VCOV = VCOV_ph, hypotheses = list(H1 = "x1 > x2 > x3"), Heq = TRUE)
  expect_null(fit_heq$posthoc)
})

test_that("posthoc = FALSE: no post-hoc analysis is stored or printed", {
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph, posthoc = FALSE)
  expect_null(fit$posthoc)
  expect_false(any(grepl("Post-hoc analysis", capture.output(print(fit)))))
})

test_that("posthoc = TRUE (specified) with a single hypothesis vs. the unconstrained", {
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = list(H1 = "x1 > x2 > x3"),
               comparison = "unconstrained", posthoc = TRUE)
  expect_true(fit$posthoc$applicable)
  direct <- goric(est_ph, VCOV = VCOV_ph, hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_equal(fit$posthoc$ratio.gw, unname(direct$ratio.gw[1, 2]))
})

test_that("internal goric() calls (evSyn, benchmark, add_Hc) do not run a post-hoc analysis", {
  n_calls <- 0
  local_mocked_bindings(goric_posthoc = function(...) { n_calls <<- n_calls + 1; NULL },
                        .package = "restriktor")
  set.seed(1)
  est_s  <- lapply(1:3, function(i) c(b1 = rnorm(1, .3, .1), b2 = rnorm(1, .1, .1), b3 = rnorm(1, 0, .1)))
  vcov_s <- lapply(1:3, function(i) diag(c(.010, .012, .011)))
  ev <- suppressMessages(evSyn(est_s, VCOV = vcov_s,
                               hypotheses = list(H1 = "b1 > b2 > b3", H2 = "b1 < b2 < b3")))
  expect_equal(n_calls, 0)

  fit <- goric(fit_lm_ph, hypotheses = list(H1 = "g1 < g2 < g3", H2 = "g1 > g2 > g3"),
               posthoc = FALSE)
  bm <- suppressMessages(benchmark(fit, model_type = "means", iter = 10))
  expect_equal(n_calls, 0)

  fit_hc <- suppressWarnings(goric(est_ph, VCOV = VCOV_ph,
                                   hypotheses = list(H1 = "x1 > x2 > x3", H2 = "x1 > x2"),
                                   add_Hc = "H2", posthoc = TRUE))
  expect_equal(n_calls, 1) # only for the final, combined set
})

test_that("posthoc must be a single TRUE or FALSE", {
  expect_error(goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph, posthoc = NA),
               "must be TRUE or FALSE")
  expect_error(goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph, posthoc = "yes"),
               "must be TRUE or FALSE")
  expect_error(goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph, posthoc = c(TRUE, FALSE)),
               "must be TRUE or FALSE")
})

test_that("posthoc = TRUE does not change the main result", {
  fit0 <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph)
  fit1 <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph, posthoc = TRUE)
  expect_equal(fit0$result, fit1$result)
  expect_equal(fit0$ratio.gw, fit1$ratio.gw)
})


# --- results -------------------------------------------------------------------

test_that("post-hoc result equals goric() of the best hypothesis vs. its complement (gorica)", {
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph, posthoc = TRUE)
  ph  <- fit$posthoc
  expect_s3_class(ph, "goric_posthoc")
  expect_true(ph$applicable)
  expect_identical(ph$best, "H1")
  expect_identical(ph$comparison, "complement")

  direct <- goric(est_ph, VCOV = VCOV_ph, hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_equal(ph$ratio.gw, unname(direct$ratio.gw[1, 2]))
  expect_equal(ph$error_prob, unname(direct$result[2, 7]))
  expect_equal(ph$error_prob, 1 / (1 + ph$ratio.gw))
  expect_equal(ph$priorICweights, c(0.5, 0.5))
})

test_that("goric_posthoc() gives the same as goric(..., posthoc = TRUE)", {
  fit0 <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph)
  fit1 <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph, posthoc = TRUE)
  ph0  <- goric_posthoc(fit0)
  expect_equal(ph0$ratio.gw, fit1$posthoc$ratio.gw)
  expect_equal(ph0$error_prob, fit1$posthoc$error_prob)
})

test_that("post-hoc works for an lm object (type = 'goric')", {
  fit <- goric(fit_lm_ph, hypotheses = list(H1 = "g1 < g2 < g3", H2 = "g1 < g2 = g3"),
               posthoc = TRUE)
  ph <- fit$posthoc
  expect_true(ph$applicable)
  expect_identical(ph$best, as.character(fit$result$model[which.max(fit$result[, 7])]))
  direct <- goric(fit_lm_ph, hypotheses = list(H1 = "g1 < g2 < g3"))
  expect_equal(ph$ratio.gw, unname(direct$ratio.gw[1, 2]))
  expect_equal(ph$result$loglik, direct$result$loglik)
})

test_that("post-hoc works with type = 'gorica' for an lm object and with unnamed hypotheses", {
  fit <- goric(fit_lm_ph, type = "gorica",
               hypotheses = list("g1 < g2 < g3", "g1 > g2 > g3"), posthoc = TRUE)
  expect_true(fit$posthoc$applicable)
  expect_identical(fit$posthoc$best, "H1")
})

test_that("best hypothesis with only equality restrictions: complement is the unconstrained", {
  est_eq <- c(x1 = 0.30, x2 = 0.31, x3 = 0.29)
  fit <- goric(est_eq, VCOV = VCOV_ph,
               hypotheses = list(Heq1 = "x1 == x2 == x3", H2 = "x1 > x2 > x3"))
  skip_if(fit$result$model[which.max(fit$result[, 7])] != "Heq1")
  ph <- suppressMessages(goric_posthoc(fit))
  expect_true(ph$applicable)
  expect_identical(ph$comparison, "unconstrained")
  expect_output(print(ph), "unconstrained")
})

test_that("post-hoc works together with add_Hc", {
  H_sub <- list(H1 = "x1 > x2 > x3", H2 = "x1 > x2")
  fit <- suppressWarnings(goric(est_ph, VCOV = VCOV_ph, hypotheses = H_sub,
                                add_Hc = "H2", posthoc = TRUE))
  ph <- fit$posthoc
  best <- as.character(fit$result$model[which.max(fit$result[, 7])])
  expect_identical(ph$best, best)
  if (best %in% c("H1", "H2")) {
    expect_true(ph$applicable)
    direct <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_sub[best])
    expect_equal(ph$ratio.gw, unname(direct$ratio.gw[1, 2]))
  } else {
    expect_false(ph$applicable)
  }
})


# --- prior IC weights ------------------------------------------------------------

test_that("unequal prior weights: best keeps its prior weight, complement gets the rest", {
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph,
               priorICweights = c(0.5, 0.25, 0.25), posthoc = TRUE)
  ph <- fit$posthoc
  expect_true(ph$priorICweights_rescaled)
  expect_equal(ph$priorICweights, c(0.5, 0.5))
  fit2 <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph,
                priorICweights = c(0.4, 0.4, 0.2), posthoc = TRUE)
  expect_equal(fit2$posthoc$priorICweights, c(0.4, 0.6))
  expect_output(print(fit2$posthoc), "prior")
})

test_that("goric_posthoc(): user-specified priorICweights are used and checked", {
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph)
  ph <- goric_posthoc(fit, priorICweights = c(3, 1))
  expect_equal(ph$priorICweights, c(0.75, 0.25))
  expect_false(ph$priorICweights_rescaled)
  expect_error(goric_posthoc(fit, priorICweights = c(1, 1, 1)), "2")
})


# --- not applicable ----------------------------------------------------------------

test_that("not applicable when the unconstrained hypothesis is the best", {
  est_u <- c(x1 = 0.5, x2 = 0.0, x3 = 0.5)
  fit <- goric(est_u, VCOV = VCOV_ph / 10, hypotheses = H_ph, posthoc = TRUE)
  expect_identical(as.character(fit$result$model[which.max(fit$result[, 7])]), "unconstrained")
  expect_false(fit$posthoc$applicable)
  expect_match(fit$posthoc$reason, "unconstrained")
  expect_output(print(fit), "Post-hoc analysis")
})

test_that("not applicable when the output is already the hypothesis vs. its complement", {
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = list(H1 = "x1 > x2 > x3"), posthoc = TRUE)
  expect_false(fit$posthoc$applicable)
  expect_match(fit$posthoc$reason, "complement")
})

test_that("not applicable when the complement is the best (comparison = 'complement')", {
  est_c <- c(x1 = 0.1, x2 = 0.3, x3 = 0.5)
  fit <- goric(est_c, VCOV = VCOV_ph, hypotheses = list(H1 = "x1 > x2 > x3"), posthoc = TRUE)
  expect_false(fit$posthoc$applicable)
})

test_that("not applicable with Heq = TRUE", {
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = list(H1 = "x1 > x2 > x3"),
               Heq = TRUE, posthoc = TRUE)
  expect_false(fit$posthoc$applicable)
})


# --- errors and output ----------------------------------------------------------------

test_that("goric_posthoc(): error for a non-goric object or an object without rerun_args", {
  expect_error(goric_posthoc(list(a = 1)), "con_goric")
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph)
  fit$rerun_args <- NULL
  expect_error(goric_posthoc(fit), "rerun_args")
})

test_that("print() and summary() show the post-hoc analysis", {
  fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph, posthoc = TRUE)
  out_p <- capture.output(print(fit), type = "output")
  out_s <- capture.output(summary(fit))
  expect_true(any(grepl("Post-hoc analysis", out_p)))
  expect_true(any(grepl("conditional error probability", out_p)))
  expect_true(any(grepl("Post-hoc analysis", out_s)))
  expect_output(print(fit$posthoc), "times more support than its complement")
})

test_that("an error in the post-hoc analysis does not stop goric()", {
  local_mocked_bindings(goric_posthoc = function(...) stop("test error"),
                        .package = "restriktor")
  fit <- NULL
  expect_message(
    fit <- goric(est_ph, VCOV = VCOV_ph, hypotheses = H_ph, posthoc = TRUE),
    "post-hoc analysis")
  expect_s3_class(fit, "con_goric")
  expect_null(fit$posthoc)
})
