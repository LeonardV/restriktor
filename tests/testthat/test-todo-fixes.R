# =============================================================================
# Tests for fixes of former TO DOs in the code
# =============================================================================

# --- calc_ICweights(): print with the requested number of digits ----------------

test_that("print.goric_ICw uses the requested number of decimals", {
  x <- calc_ICweights(c(100, 102, 130, 160), c("configural", "weak", "strong", "strict"))
  out <- capture.output(print(x, use_scientific = FALSE, digits = 1))
  expect_true(any(grepl("100.0", out, fixed = TRUE)))
  expect_false(any(grepl("0.0000008", out, fixed = TRUE)))
  out_sci <- capture.output(print(x, digits = 3))
  # columns with values that would be printed as zero are in scientific notation
  expect_true(any(grepl("e-07", out_sci, fixed = TRUE)))
  # a column without such values is not
  expect_true(any(grepl("100.000", out_sci, fixed = TRUE)))
})


# --- evSyn(): checks on VCOV per study and on the goric objects ------------------

test_that("evSyn() reports the study with a non-square or non-matching VCOV", {
  est  <- list(c(a = .3, b = .1), c(a = .2, b = .1))
  vc   <- list(diag(2) * .01, matrix(.01, 2, 3))
  expect_error(evSyn(est, VCOV = vc, hypotheses = list(H1 = "a > b")),
               "study 2 is not a square matrix")
  vc2  <- list(diag(2) * .01, diag(3) * .01)
  expect_error(evSyn(est, VCOV = vc2, hypotheses = list(H1 = "a > b")),
               "study 2 is 3 x 3")
})

test_that("evSyn() with goric objects checks the number of hypotheses", {
  g1 <- goric(c(a = .3, b = .1), VCOV = diag(2) * .01, hypotheses = list(H1 = "a > b"))
  g2 <- goric(c(a = .3, b = .1, c = 0), VCOV = diag(3) * .01,
              hypotheses = list(H1 = "a > b", H2 = "a > c"))
  expect_error(evSyn(list(g1, g2)), "same number of hypotheses")
  g3 <- goric(c(a = .2, b = .1), VCOV = diag(2) * .01, hypotheses = list(Hx = "a > b"))
  expect_message(evSyn(list(g1, g3)), "names of the hypotheses")
})


# --- leave1studyout(): print with fixed decimals ---------------------------------

test_that("print of leave1studyout() shows weights with fixed decimals (e.g., 1.000)", {
  set.seed(1)
  est <- lapply(1:3, function(i) c(b1 = rnorm(1, 1, .05), b2 = rnorm(1, 0, .05)))
  vc  <- lapply(1:3, function(i) diag(c(.001, .001)))
  r   <- suppressMessages(evSyn(est, VCOV = vc, hypotheses = list(H1 = "b1 > b2")))
  l1  <- leave1studyout(r)
  out <- capture.output(print(l1))
  expect_true(any(grepl("1.000", out, fixed = TRUE)))
  expect_true(any(grepl("Preferred hypothesis", out)))
})


# --- goric(): a matrix of estimates is turned into a named vector ----------------

test_that("goric() accepts a numeric matrix of estimates (names from its dimnames or VCOV)", {
  m <- matrix(c(.5, .3, .1, .2), 2, dimnames = list(c("a", "b"), c("c", "d")))
  fit <- NULL
  expect_message(fit <- goric(m, VCOV = diag(4) * .01, hypotheses = list(H1 = "a.c > b.c")),
                 "converted into a vector")
  ref <- goric(c(a.c = .5, b.c = .3, a.d = .1, b.d = .2), VCOV = diag(4) * .01,
               hypotheses = list(H1 = "a.c > b.c"))
  expect_equal(fit$result, ref$result)
  V <- diag(4) * .01; rownames(V) <- colnames(V) <- c("p1", "p2", "p3", "p4")
  fit2 <- suppressMessages(goric(m, VCOV = V, hypotheses = list(H1 = "p1 > p2")))
  expect_equal(fit2$result$loglik, ref$result$loglik)
  expect_error(suppressMessages(goric(matrix(c(.5, .3), 1), VCOV = diag(2) * .01,
                                      hypotheses = list(H1 = "a > b"))),
               "named vector")
})

test_that("goric(): complement with a selection of the parameters in the hypothesis", {
  est <- c(x1 = .5, x2 = .3, x3 = .1)
  V <- diag(3) * .01; V[1, 2] <- V[2, 1] <- .003
  fit_all <- goric(est, VCOV = V, hypotheses = list(H1 = "x1 > x2"))
  fit_sub <- goric(est[1:2], VCOV = V[1:2, 1:2], hypotheses = list(H1 = "x1 > x2"))
  expect_equal(fit_all$result$gorica.weights, fit_sub$result$gorica.weights)
})


# --- goric() for rlm objects --------------------------------------------------------

test_that("goric() works for an rlm object (method = 'MM')", {
  set.seed(2)
  d <- data.frame(y = rnorm(90, rep(c(0, .5, 1), each = 30)), g = gl(3, 30))
  m <- MASS::rlm(y ~ g - 1, data = d, method = "MM")
  fit <- goric(m, hypotheses = list(H1 = "g1 < g2 < g3"))
  expect_s3_class(fit, "con_goric")
  expect_equal(sum(fit$result$goric.weights), 1)
})


# --- lavaan: duplicate labels (equality-constrained parameters) ---------------------

test_that("goric() for lavaan with a label used for more than one parameter", {
  set.seed(3)
  n <- 200
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- .3 * d$x1 + .3 * d$x2 + .1 * d$x3 + rnorm(n)
  fit_lav <- lavaan::sem("y ~ a*x1 + a*x2 + c*x3", data = d)
  est <- restriktor:::con_gorica_est_lav(fit_lav)
  expect_false(any(duplicated(names(est$estimate)[names(est$estimate) != ""])))
  expect_equal(nrow(est$VCOV), length(est$estimate))
  fit <- suppressMessages(goric(fit_lav, hypotheses = list(H1 = "a > c")))
  expect_s3_class(fit, "con_goric")
})


# --- coef_named_vector(): names from VCOV -------------------------------------------

test_that("coef_named_vector() gives an informative error when names are unavailable", {
  fake <- structure(list(coefficients = matrix(1:4, 2)), class = "fakefit")
  coef.fakefit <- function(object, ...) object$coefficients
  vcov.fakefit <- function(object, ...) diag(4)  # no names
  registerS3method("coef", "fakefit", coef.fakefit)
  registerS3method("vcov", "fakefit", vcov.fakefit)
  expect_error(restriktor:::coef_named_vector(fake), "cannot be obtained")
  V <- diag(4); rownames(V) <- colnames(V) <- paste0("p", 1:4)
  expect_equal(names(suppressMessages(restriktor:::coef_named_vector(fake, VCOV = V))),
               paste0("p", 1:4))
})


# --- benchmark print: number of iterations -----------------------------------------

test_that("print of benchmark() shows the number of iterations", {
  set.seed(4)
  d <- data.frame(y = rnorm(60, rep(c(0, .3, .6), each = 20)), g = gl(3, 20))
  fit <- goric(lm(y ~ g - 1, data = d), hypotheses = list(H1 = "g1 < g2 < g3"))
  bm <- suppressMessages(benchmark(fit, model_type = "means", iter = 10))
  expect_output(print(bm), "Number of iterations \\(draws\\) per benchmark")
})
