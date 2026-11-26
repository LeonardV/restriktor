# =============================================================================
# Tests: wrapping of long output text and collecting messages/warnings in goric()
# =============================================================================

test_that(".rk_wrap() wraps running text, but not tables or short lines", {
  wrap <- restriktor:::.rk_wrap
  long <- paste(rep("word", 40), collapse = " ")
  out <- wrap(long, width = 76)
  expect_true(all(nchar(strsplit(out, "\n")[[1]]) <= 76))
  expect_identical(gsub("\n", " ", out), long)
  # table-like line (aligned with several spaces): as is
  tab <- paste(rep("col     ", 15), collapse = "")
  expect_identical(wrap(tab, width = 76), tab)
  # short lines and trailing line breaks are kept
  expect_identical(wrap("short\n", width = 76), "short\n")
  expect_identical(wrap(paste0("\n", long, "\n"), width = 76), paste0("\n", out, "\n"))
  # bullets: continuation lines line up with the text of the bullet
  b <- wrap(paste("- ", long), width = 40)
  expect_true(all(startsWith(strsplit(b, "\n")[[1]][-1], "   ")))
  # Inf: no wrapping
  expect_identical(wrap(long, width = Inf), long)
})

test_that("messages, warnings, and errors of restriktor are wrapped (option restriktor.linewidth)", {
  old <- options(restriktor.linewidth = 50)
  on.exit(options(old))
  m <- tryCatch(restriktor:::message(paste(rep("word", 30), collapse = " ")),
                message = function(m) conditionMessage(m))
  expect_true(all(nchar(strsplit(m, "\n")[[1]]) <= 50))
  e <- tryCatch(restriktor:::stop(paste(rep("word", 30), collapse = " "), call. = FALSE),
                error = function(e) conditionMessage(e))
  expect_true(all(nchar(strsplit(e, "\n")[[1]]) <= 50))
  w <- tryCatch(restriktor:::warning(paste(rep("word", 30), collapse = " ")),
                warning = function(w) w)
  expect_true(all(nchar(strsplit(conditionMessage(w), "\n")[[1]]) <= 50))
  out <- capture.output(restriktor:::cat(paste(rep("word", 30), collapse = " "), "\n"))
  expect_true(all(nchar(out) <= 50))
  # conditions are passed on as they are
  cond <- simpleError("x")
  expect_error(restriktor:::stop(cond), "^x$")
})

test_that("the printed output of goric() is wrapped", {
  old <- options(restriktor.linewidth = 60)
  on.exit(options(old))
  set.seed(1)
  d <- data.frame(y = rnorm(60, rep(c(0, .3, .6), each = 20)), g = gl(3, 20))
  fit <- goric(lm(y ~ g - 1, data = d), hypotheses = list(H1 = "g1 < g2 < g3", H2 = "g1 < g2 = g3"))
  out <- capture.output(print(fit), type = "output")
  conclusion <- out[(grep("^Conclusion", out) + 1):(grep("^Conclusion", out) + 4)]
  expect_true(all(nchar(conclusion) <= 60))
})

test_that("goric() collects the messages and warnings; print() reports them", {
  est <- c(x1 = .5, x2 = .3, x3 = .1)
  m <- matrix(c(.5, .3, .1, .2), 2, dimnames = list(c("a", "b"), c("c", "d")))
  fit <- suppressMessages(goric(m, VCOV = diag(4) * .01, hypotheses = list(H1 = "a.c > b.c")))
  expect_true(any(grepl("converted into a vector", fit$run_messages)))
  expect_length(fit$run_warnings, 0)
  out <- capture.output(print(fit))
  expect_true(any(grepl("message was given when running goric", out)))
  expect_true(any(grepl("run_messages", out)))
  # no messages or warnings: no note
  fit2 <- goric(est, VCOV = diag(3) * .01, hypotheses = list(H1 = "x1 > x2"))
  expect_length(fit2$run_messages, 0)
  expect_false(any(grepl("given when running goric", capture.output(print(fit2)))))
  # warnings are repeated by print()
  fit2$run_warnings <- "a test warning"
  out2 <- capture.output(print(fit2))
  expect_true(any(grepl("1 warning was given", out2)))
  expect_true(any(grepl("a test warning", out2)))
})

test_that("goric() still finds objects in the caller's environment (e.g., arima data)", {
  f <- function() {
    y_local <- as.numeric(arima.sim(list(ar = .5), n = 100))
    fit <- arima(y_local, order = c(1, 0, 0))
    suppressMessages(goric(fit, hypotheses = list(H1 = "ar1 > 0"), type = "goric"))
  }
  set.seed(5)
  expect_identical(f()$type, "goric")
})
