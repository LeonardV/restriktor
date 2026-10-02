# =============================================================================
# Tests: regressietests voor de bevindingen van de audit (kern goric)
# A1, A2, A3, A5, A13, B1, B2, B3, B4, B5, B8, B9, B10, B14
# =============================================================================

est_a <- c(x1 = 5, x2 = 3, x3 = 1)
VCOV_a <- diag(3) * 0.1
est_b <- c(x1 = 0.5, x2 = 0.8, x3 = 0.3)
VCOV_b <- diag(3) * 0.05

cap_all <- function(expr) {
  msgs <- character(); warns <- character(); err <- NULL
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) { err <<- conditionMessage(e); NULL }),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(value = val, messages = msgs, warnings = warns, error = err)
}

# -----------------------------------------------------------------------------
# A1: summary() met Heq = TRUE rapporteert H1 vs. complement (op naam)
# -----------------------------------------------------------------------------

test_that("summary.con_goric met Heq: zin gaat over H1 vs. complement, niet Heq vs. H1", {
  g <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2 > x3"),
             Heq = TRUE, type = "gorica")
  expect_identical(g$result$model, c("Heq", "H1", "complement"))
  # onafhankelijk: ratio = exp((IC_c - IC_H1) / 2)
  ic <- setNames(g$result$gorica, g$result$model)
  ratio_hand <- exp((ic["complement"] - ic["H1"]) / 2)
  expect_equal(unname(g$ratio.gw["H1", "vs. complement"]), unname(ratio_hand))
  out <- capture.output(s <- summary(g))
  line <- grep("times more supported than its complement", out, value = TRUE)
  expect_length(line, 1L)
  expect_match(line, "'H1' is 5.07e\\+04 times more supported than its complement")
  expect_false(any(grepl("'Heq'.*times more", out)))
  # summary() geeft het object onzichtbaar terug; print(summary()) eindigt niet op "NULL"
  expect_s3_class(s, "con_goric")
  out2 <- capture.output(print(summary(g)))
  expect_false(any(out2 == "NULL"))
})

test_that("summary.con_goric met Heq: als Heq het beste is wordt dat gemeld (zoals print)", {
  g <- goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 > x2 > x3"), Heq = TRUE)
  expect_equal(g$result$model[which.max(g$result$gorica.weights)], "Heq")
  out <- capture.output(summary(g))
  expect_true(any(grepl("Heq\\) is the best in the set", out)))
  expect_false(any(grepl("'Heq' (is|has) .*times more", out)))
})

# -----------------------------------------------------------------------------
# A2 / B10: conclusietekst bij gelijke steun (ties), op naam geindexeerd
# -----------------------------------------------------------------------------

test_that("print.con_goric: ties worden als gelijke steun gemeld en ratio's staan bij de juiste hypothese", {
  g <- goric(c(x1 = 1, x2 = 1, x3 = 1), VCOV = diag(3),
             hypotheses = list(H1 = "x1 > 0", H2 = "x2 > 0", H3 = "x3 < 0"))
  # onafhankelijk: H1 en H2 identiek (ratio 1); H1 vs H3 = exp(0.5) (loglik
  # verschil 0.5 bij gelijke penalty)
  expect_equal(unname(g$ratio.gw["H1", "vs. H2"]), 1)
  expect_equal(unname(g$ratio.gw["H1", "vs. H3"]), exp(0.5))
  out <- capture.output(print(g))
  expect_true(any(grepl("'H1' and 'H2' have equal support", out)))
  expect_true(any(grepl("'H1' is 1.649 times more supported than 'H3'", out)))
  expect_false(any(grepl("NULL|NA times", out)))

  # beste hypothese niet op de eerste positie
  g2 <- goric(c(x1 = 1, x2 = 1, x3 = 1), VCOV = diag(3),
              hypotheses = list(H1 = "x3 < 0", H2 = "x1 > 0", H3 = "x2 > 0"))
  out2 <- capture.output(print(g2))
  expect_true(any(grepl("'H2' is 1.649 times more supported than 'H1'", out2)))
  expect_true(any(grepl("'H2' and 'H3' have equal support", out2)))

  # twee hypothesen met gelijke steun
  g3 <- goric(c(x1 = 1, x2 = 1), VCOV = diag(2),
              hypotheses = list(H1 = "x1 > 0", H2 = "x2 > 0"))
  out3 <- capture.output(print(g3))
  expect_true(any(grepl("'H1' and 'H2' have equal support", out3)))
  expect_false(any(grepl("NULL", out3)))
})

test_that("print.con_goric: comparison = 'none' met 2 hypothesen formuleert vanuit de beste hypothese", {
  g <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 < x2", H2 = "x1 > x2"),
             comparison = "none", type = "gorica")
  # onafhankelijk: loglik-verschil = (5-3)^2 / (2 * 2 * 0.1) = 10 -> ratio exp(10)
  expect_equal(unname(g$ratio.gw["H2", "vs. H1"]), exp(10), tolerance = 1e-8)
  out <- capture.output(print(g))
  expect_false(any(grepl("0.00 times", out)))
  expect_true(any(grepl("'H2' is 2.20e\\+04 times more supported than 'H1'", out)))
})

test_that("print.con_goric: prior 0 geeft 'infinitely more supported', geen NaN/NULL", {
  g <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2", H2 = "x1 > x3"),
             priorICweights = c(1, 0, 0), type = "gorica")
  out <- capture.output(print(g))
  expect_true(any(grepl("'H1' is infinitely more supported than 'H2'", out)))
  expect_false(any(grepl("NaN|NULL", out)))
  out2 <- capture.output(print(summary(g)))
  expect_false(any(out2 == "NULL"))
})

# -----------------------------------------------------------------------------
# A3: summary() met 1 hypothese en comparison = "none"
# -----------------------------------------------------------------------------

test_that("summary.con_goric werkt met 1 hypothese en comparison = 'none'", {
  g <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2 > x3"),
             comparison = "none", type = "gorica")
  expect_equal(dim(g$ratio.gw), c(1L, 1L))
  out <- NULL
  expect_no_error(out <- capture.output(summary(g)))
  expect_true(any(grepl("vs. H1", out)))
  expect_true(any(grepl("H1 \\(best\\)", out)))
  expect_no_error(capture.output(summary(g, brief = FALSE)))
})

# -----------------------------------------------------------------------------
# A5: format_numeric met NA/NaN/Inf
# -----------------------------------------------------------------------------

test_that("format_numeric: NA, NaN en Inf geven geen fout", {
  expect_identical(restriktor:::format_numeric(NA), "NA")
  expect_identical(restriktor:::format_numeric(NA_real_), "NA")
  expect_identical(restriktor:::format_numeric(NaN), "NaN")
  expect_identical(restriktor:::format_numeric(Inf), "Inf")
  expect_identical(restriktor:::format_numeric(-Inf), "-Inf")
  expect_identical(restriktor:::format_numeric(0.5), "0.500")
  expect_identical(restriktor:::format_numeric(0), "0.000")
})

# -----------------------------------------------------------------------------
# A13: ':=' gedefinieerde parameters in coef()/ormle$b.restr (numeric/lavaan)
# -----------------------------------------------------------------------------

test_that("goric numeric: gedefinieerde parameter (:=) in b.restr, zonder rbind-warning", {
  est <- c(a = 0.9, b = 0.5, c = 0.2)
  V <- diag(3) * 0.04
  r <- cap_all(goric(est, VCOV = V, hypotheses = list(H1 = "d := a - b; d > 0; c > d")))
  expect_length(r$warnings, 0L)
  b <- r$value$ormle$b.restr
  expect_identical(colnames(b), c("a", "b", "c", "d"))
  # onafhankelijk: projectie met actieve restrictie c - a + b >= 0 (gelijke
  # varianties): a = .9 - t, b = .5 + t, c = .2 + t met t = 1/15
  t <- 1 / 15
  expect_equal(unlist(b["H1", ]), c(a = 0.9 - t, b = 0.5 + t, c = 0.2 + t, d = 0.4 - 2 * t),
               tolerance = 1e-8)
  expect_equal(unlist(b["complement", ]), c(a = 0.9, b = 0.5, c = 0.2, d = 0.4),
               tolerance = 1e-8)
  expect_equal(coef(r$value), b)

  # unconstrained: hypothese zonder := krijgt NA voor d
  r2 <- cap_all(goric(est, VCOV = V, comparison = "unconstrained",
                      hypotheses = list(H1 = "d := a - b; d > 0; c > d", H2 = "a > b")))
  expect_length(r2$warnings, 0L)
  b2 <- r2$value$ormle$b.restr
  expect_identical(colnames(b2), c("a", "b", "c", "d"))
  expect_true(is.na(b2["H2", "d"]))
  expect_equal(unname(unlist(b2["unconstrained", ])), c(0.9, 0.5, 0.2, 0.4), tolerance = 1e-8)
})

test_that("goric lavaan: gedefinieerde parameter (:=) in de hypothese, zonder warning", {
  skip_if_not_installed("lavaan")
  set.seed(5)
  n <- 100
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 0.6 * d$x1 + 0.2 * d$x2 + rnorm(n)
  fit <- lavaan::sem("y ~ a*x1 + b*x2", data = d)
  r <- cap_all(goric(fit, hypotheses = list(H1 = "diff := a - b; diff > 0"), type = "gorica"))
  expect_null(r$error)
  expect_length(r$warnings, 0L)
  b <- r$value$ormle$b.restr
  expect_true("diff" %in% colnames(b))
  expect_equal(unname(b["complement", "diff"]), unname(b["complement", "a"] - b["complement", "b"]))
})

# -----------------------------------------------------------------------------
# B1: gereserveerde namen en namen die 'Heq'/'unconstrained' bevatten
# -----------------------------------------------------------------------------

test_that("goric: gereserveerde hypothesenamen geven een duidelijke fout", {
  for (nm in c("Heq", "complement", "unconstrained")) {
    hyp <- list("x1 > x2", "x2 > x3")
    names(hyp) <- c(nm, "H2")
    expect_error(goric(est_b, VCOV = VCOV_b, hypotheses = hyp), "reserved")
  }
})

test_that("goric: namen die 'Heq' of 'unconstrained' bevatten worden gewoon behandeld", {
  g <- goric(est_b, VCOV = VCOV_b, hypotheses = list(unconstrained1 = "x1 > x2", H2 = "x2 > x3"))
  expect_identical(g$result$model, c("unconstrained1", "H2", "unconstrained"))
  g_ref <- goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 > x2", H2 = "x2 > x3"))
  expect_identical(names(g$result), names(g_ref$result))
  expect_equal(g$result$gorica.weights, g_ref$result$gorica.weights)

  g <- goric(est_b, VCOV = VCOV_b, hypotheses = list(myHeq = "x1 > x2 > x3"), Heq = TRUE)
  expect_identical(g$result$model, c("Heq", "myHeq", "complement"))
  g <- goric(est_b, VCOV = VCOV_b, hypotheses = list(Heqx = "x1 > x2 > x3"), Heq = TRUE)
  expect_identical(g$result$model, c("Heq", "Heqx", "complement"))

  a <- goric(est_a, VCOV = VCOV_a, comparison = "none", type = "gorica",
             hypotheses = list(H1 = "x1 > x2", MyHeq = "x1 < x2", H3 = "x2 > x3"))
  expect_false(any(grepl("without", names(a$result))))
  a <- goric(est_a, VCOV = VCOV_a, comparison = "none", type = "gorica",
             hypotheses = list(H1 = "x1 > x2", unconstrained_alt = "x1 < x2", H3 = "x2 > x3"))
  expect_false(any(grepl("without", names(a$result))))
})

test_that("goric: een Heq-element in hypotheses (zoals uit benchmark) wordt op naam verwijderd", {
  g <- goric(est_b, VCOV = VCOV_b, Heq = TRUE,
             hypotheses = list(Heq = "x1 = x2 = x3", H1 = "x1 > x2 > x3"))
  g_ref <- goric(est_b, VCOV = VCOV_b, Heq = TRUE, hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_identical(g$result$model, c("Heq", "H1", "complement"))
  expect_equal(g$result, g_ref$result)
})

# -----------------------------------------------------------------------------
# B2: priorICweights als matrix, numeric(0), met namen
# -----------------------------------------------------------------------------

test_that("goric: priorICweights als 1 x k matrix wordt geaccepteerd", {
  g_m <- goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 > x2"),
               priorICweights = matrix(c(0.5, 0.5), 1))
  g_v <- goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 > x2"))
  expect_equal(g_m$result$gorica.weights, g_v$result$gorica.weights)
  expect_null(dim(g_m$priorICweights))
})

test_that("goric: priorICweights numeric(0) geeft de lengte-melding (correcte grammatica)", {
  expect_error(goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 > x2"),
                     priorICweights = numeric(0)),
               "should consist of 2 elements.*now consists of 0 elements")
  expect_error(goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 > x2"),
                     comparison = "none", priorICweights = c(0.5, 0.5)),
               "should consist of 1 element,.*now consists of 2 elements")
})

test_that("goric: benoemde priorICweights worden op naam toegepast en in die volgorde bewaard", {
  g <- goric(est_a, VCOV = VCOV_a, comparison = "none", type = "gorica",
             hypotheses = list(H1 = "x1 > x2 > x3", H2 = "x1 > x2 & x1 > x3"),
             priorICweights = c(H2 = 0.9, H1 = 0.1))
  expect_equal(g$priorICweights, c(H1 = 0.1, H2 = 0.9))
  # onafhankelijk: w_i = p_i exp(-IC_i/2) / som
  w <- exp(-g$result$gorica / 2)
  expect_equal(g$result$gorica.weights, c(0.1, 0.9) * w / sum(c(0.1, 0.9) * w))
  g_pos <- goric(est_a, VCOV = VCOV_a, comparison = "none", type = "gorica",
                 hypotheses = list(H1 = "x1 > x2 > x3", H2 = "x1 > x2 & x1 > x3"),
                 priorICweights = c(0.1, 0.9))
  expect_equal(g$result$gorica.weights, g_pos$result$gorica.weights)

  # Heq en complement op naam (in willekeurige volgorde)
  g2 <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2 > x3"), Heq = TRUE,
              type = "gorica", priorICweights = c(complement = 0.5, H1 = 0.2, Heq = 0.3))
  expect_equal(g2$priorICweights, c(Heq = 0.3, H1 = 0.2, complement = 0.5))
  # namen die niet overeenkomen
  expect_error(goric(est_a, VCOV = VCOV_a, comparison = "none", type = "gorica",
                     hypotheses = list(H1 = "x1 > x2 > x3", H2 = "x1 > x2 & x1 > x3"),
                     priorICweights = c(H3 = 0.9, H1 = 0.1)), "do not match")
})

# -----------------------------------------------------------------------------
# B3: Heq = TRUE bij een hypothese zonder (overblijvende) ongelijkheden
# -----------------------------------------------------------------------------

test_that("goric: Heq wordt genegeerd (met melding) bij een hypothese met alleen gelijkheden", {
  r <- cap_all(goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 = x2"), Heq = TRUE))
  expect_null(r$error)
  expect_true(any(grepl("'Heq' argument is ignored", r$messages)))
  g <- r$value
  expect_false(g$Heq)
  expect_identical(g$result$model, c("H1", "unconstrained"))
  # onafhankelijk: projectie x1 = x2 = 0.65, loglik-verschil -0.45, PT 2 vs 3
  # -> w_H1 = 1 / (1 + exp(-0.55))
  expect_equal(g$result$gorica.weights[1], 1 / (1 + exp(-0.55)), tolerance = 1e-8)
  expect_length(g$priorICweights, 2L)

  # ongelijkheid die door de gelijkheid wordt geimpliceerd
  r2 <- cap_all(goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 = 0; x1 > 0"), Heq = TRUE))
  expect_null(r2$error)
  expect_true(any(grepl("'Heq' argument is ignored", r2$messages)))
  expect_identical(r2$value$result$model, c("H1", "unconstrained"))

  # priors: lengte na de beslissing over Heq
  g3 <- suppressMessages(goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 = x2"),
                               Heq = TRUE, priorICweights = c(0.3, 0.7)))
  expect_equal(g3$priorICweights, c(0.3, 0.7))
})

# -----------------------------------------------------------------------------
# B4: Heq die niet gevormd kan worden (strijdig / range)
# -----------------------------------------------------------------------------

test_that("goric: Heq met strijdige gelijkheden geeft een duidelijke restriktor-fout", {
  expect_error(goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 > x2 > 0.1; x1 > 0.05"),
                     Heq = TRUE), "restriktor ERROR.*Heq.*cannot be formed")
  expect_error(goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "2*x1 > 1; x1 > 0.2; x2 > 0"),
                     Heq = TRUE), "restriktor ERROR.*Heq.*cannot be formed")
})

test_that("goric: Heq bij een range-restrictie geeft een duidelijke fout", {
  expect_error(goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "0.2 < x1 < 0.5"), Heq = TRUE),
               "restriktor ERROR.*Heq.*range restriction")
  expect_error(restriktor:::goric_heq_constraints(est_b, "x1 > 0.2; x1 < 0.5"),
               "range restriction")
  # dezelfde grens aan beide kanten is een gelijkheid: wel toegestaan
  expect_no_error(g <- goric(est_b, VCOV = VCOV_b, hypotheses = list(H1 = "x1 > 1; x1 < 1"), Heq = TRUE))
  expect_identical(g$result$model, c("Heq", "H1", "complement"))
  # werkende gevallen blijven werken
  expect_identical(restriktor:::goric_heq_constraints(est_b, "x1 > 0.5; x1 > 0.2; x2 > 0"),
                   "x1=0.5; x2=0")
  expect_null(restriktor:::goric_heq_constraints(est_b, "x1 = x2"))
  expect_null(restriktor:::goric_heq_constraints(est_b, "x1 = 0.5; x1 > 0.2; x1 > 0.1"))
})

# -----------------------------------------------------------------------------
# B5: mlm: melding over de coefficientnamen; goricc via summary.restriktor
# -----------------------------------------------------------------------------

set.seed(11)
n_mlm <- 40
d_mlm <- data.frame(x1 = rnorm(n_mlm), x2 = rnorm(n_mlm))
d_mlm$y1 <- 1 + 0.8 * d_mlm$x1 + rnorm(n_mlm)
d_mlm$y2 <- -0.5 + 0.3 * d_mlm$x1 + rnorm(n_mlm)
fit_mlm2 <- lm(cbind(y1, y2) ~ x1 + x2, d_mlm)

test_that("goric mlm: melding over de coefficientnamen (incl. intercept) bij elke comparison", {
  for (cmp in c("none", "unconstrained", "complement")) {
    r <- cap_all(goric(fit_mlm2, hypotheses = list(H1 = "y1..Intercept. > 0"), comparison = cmp))
    expect_null(r$error)
    msg <- grep("coefficient names", r$messages, value = TRUE)
    expect_length(msg, 1L)
    expect_match(msg, "y1..Intercept.", fixed = TRUE)
    expect_match(msg, "vcov")
  }
  expect_error(goric(fit_mlm2, hypotheses = list(H1 = "y1.(Intercept) > 0")))
})

test_that("summary.restriktor: goricc/goricac geven een fout voor mlm-objecten", {
  rs <- restriktor(fit_mlm2, constraints = "y1.x1 > 0")
  expect_error(summary(rs, goric = "goricc"), "not \\(yet\\) available for objects of class mlm")
  expect_error(summary(rs, goric = "GORICAC"), "not \\(yet\\) available for objects of class mlm")
  expect_no_error(s <- summary(rs, goric = "goric"))
  expect_true(is.finite(s$goric))
})

# -----------------------------------------------------------------------------
# B8 / B9: geen partial matching ($goric, se -> seed)
# -----------------------------------------------------------------------------

set.seed(1)
d_g <- data.frame(g = factor(rep(1:3, each = 10)))
d_g$y <- c(0, 0.5, 1)[d_g$g] + rnorm(30)
fit_g <- lm(y ~ -1 + g, d_g)

test_that("goric/print/summary geven geen partial-match warnings (alle types en comparisons)", {
  op <- options(warnPartialMatchDollar = TRUE, warnPartialMatchArgs = TRUE,
                warnPartialMatchAttr = TRUE)
  on.exit(options(op), add = TRUE)
  for (ty in c("goric", "goricc", "gorica", "goricac")) {
    for (cmp in c("unconstrained", "none")) {
      r <- cap_all({
        g <- goric(fit_g, hypotheses = list(H1 = "g1 < g2 < g3", H2 = "g1 > g2"),
                   type = ty, comparison = cmp)
        capture.output(print(g))
        capture.output(summary(g))
        capture.output(summary(g, brief = FALSE))
        g
      })
      expect_null(r$error, info = paste(ty, cmp))
      expect_length(r$warnings, 0L)
      expect_true(all(is.finite(r$value$result[[ty]])))
    }
    # complement (incl. con_gorica_est via compute_complement_likelihood: 'se')
    r <- cap_all({
      g <- goric(fit_g, hypotheses = list(H1 = "g1 < g2 < g3"), type = ty,
                 comparison = "complement", Heq = TRUE)
      capture.output(print(g))
      capture.output(summary(g))
      g
    })
    expect_null(r$error, info = ty)
    expect_length(r$warnings, 0L)
    expect_identical(r$value$result$model, c("Heq", "H1", "complement"))
  }
})

test_that("compute_complement_likelihood: 'se' wordt niet aan 'seed' gebonden (boot)", {
  op <- options(warnPartialMatchArgs = TRUE)
  on.exit(options(op), add = TRUE)
  r <- cap_all(goric(fit_g, hypotheses = list(H1 = "g1 < g2 < g3"), type = "gorica",
                     comparison = "complement", mix_weights = "boot",
                     control = list(mix_weights_bootstrap_limit = 200)))
  expect_null(r$error)
  expect_false(any(grepl("partial", r$warnings)))
})

# -----------------------------------------------------------------------------
# B14: ic_weights_log met -Inf; goric met oneindige penalty
# -----------------------------------------------------------------------------

test_that("ic_weights_log: IC = -Inf geeft gewicht 1 (gedeeld bij ties), geen NaN", {
  expect_equal(restriktor:::ic_weights_log(c(-Inf, 0, 1)), c(1, 0, 0))
  expect_equal(restriktor:::ic_weights_log(c(-Inf, -Inf, 1)), c(0.5, 0.5, 0))
  expect_equal(restriktor:::ic_weights_log(c(-Inf, 0), prior = c(0, 1)), c(0, 1))
  expect_equal(restriktor:::ic_weights_log(c(Inf, 0)), c(0, 1))
  expect_true(all(is.nan(restriktor:::ic_weights_log(c(NaN, 0)))))
  expect_true(all(is.nan(restriktor:::ic_weights_log(c(Inf, Inf)))))
  # namen op de prior worden niet doorgegeven aan de gewichten
  expect_null(names(restriktor:::ic_weights_log(c(1, 2), prior = c(a = 0.5, b = 0.5))))
})

test_that("goric: te kleine N voor goricac/goricc (N - p - 2 <= 0) geeft vooraf een duidelijke fout", {
  # N - p - 2 = 0: oneindige penalty; N - p - 2 < 0: negatieve penalty
  for (cmp in c("complement", "unconstrained", "none")) {
    for (N in c(5, 4)) {
      expect_error(goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2 > x3"),
                         type = "goricac", sample_nobs = N, comparison = cmp),
                   "restriktor ERROR.*sample size is too small.*goricac.*N = .*p = 3",
                   info = paste(cmp, N))
    }
  }
  # N - p - 2 = 1: wel toegestaan, eindige gewichten (ook bij comparison = 'none')
  r <- cap_all(goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2 > x3", H2 = "x1 > x2"),
                     type = "goricac", sample_nobs = 6, comparison = "none"))
  expect_null(r$error)
  expect_true(all(is.finite(r$value$result$goricac.weights)))
  # lm-route (goricc): N uit het model, p = aantal coefficienten
  expect_error(goric(lm(y ~ -1 + g, d_g[c(1, 2, 11, 12, 21), ]), hypotheses = list(H1 = "g1 < g2 < g3"),
                     type = "goricc"), "sample size is too small.*goricc.*N = 5 and p = 3")
})

# =============================================================================
# W1: residuele bevindingen van de her-verificatie
# =============================================================================

# -----------------------------------------------------------------------------
# W1-N-01: summary(brief = FALSE) met een gedefinieerde parameter (:=)
# -----------------------------------------------------------------------------

test_that("summary(brief = FALSE) werkt met ':=' (numeric, lm en lavaan)", {
  # numeric, alle comparisons
  for (cmp in c("none", "complement", "unconstrained")) {
    g <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "d := x1 - x2; d > 0"),
               comparison = cmp, type = "gorica")
    out <- NULL
    expect_no_error(out <- capture.output(summary(g, brief = FALSE)), 
                    message = paste("numeric", cmp))
    expect_true(any(grepl("^\\s+x1\\s+x2\\s+x3\\s+op\\s+rhs\\s+active", out)), 
                info = paste("numeric", cmp))
  }
  # lm, met een restrictie op de gedefinieerde parameter en op een gewone
  set.seed(3)
  n <- 50
  dd <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  dd$y <- dd$a + 0.5 * dd$b + rnorm(n)
  fit <- lm(y ~ a + b + c, dd)
  for (cmp in c("complement", "unconstrained", "none")) {
    g <- goric(fit, hypotheses = list(H1 = "d := a - b; d > 0; c > d"), comparison = cmp)
    expect_no_error(capture.output(summary(g, brief = FALSE)), message = paste("lm", cmp))
  }
  # lavaan
  skip_if_not_installed("lavaan")
  dl <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  dl$y <- dl$x1 + rnorm(60)
  fl <- lavaan::sem("y ~ a*x1 + b*x2\n d := a - b", dl)
  g <- suppressMessages(goric(fl, hypotheses = list(H1 = "d > 0"), type = "gorica"))
  expect_no_error(capture.output(summary(g, brief = FALSE)))
})

# -----------------------------------------------------------------------------
# W1-N-02: conclusiezinnen via support_sentence() (prior 0, 2 modellen)
# -----------------------------------------------------------------------------

test_that("print/summary: complement met IC-gewicht 0 geeft een leesbare zin", {
  g <- suppressMessages(goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2 > x3"),
                              priorICweights = c(0, 1), type = "gorica"))
  out <- capture.output(print(g))
  expect_true(any(grepl("^Its complement is infinitely more supported than 'H1' \\(which has an IC weight of 0\\)", out)))
  expect_false(any(grepl("hypothesis its complement", out)))
  out2 <- capture.output(summary(g))
  expect_true(any(grepl("^Its complement is infinitely more supported than 'H1'", out2)))
  # prior 0 op het complement: zin begint met de hypothese
  g2 <- suppressMessages(goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2 > x3"),
                               priorICweights = c(1, 0), type = "gorica"))
  out3 <- capture.output(print(g2))
  expect_true(any(grepl("^The order-restricted hypothesis 'H1' is infinitely more supported than its complement", out3)))
})

test_that("print: 2 modellen met comparison = 'unconstrained' gebruikt support_sentence", {
  g <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2"),
             comparison = "unconstrained", type = "gorica")
  out <- capture.output(print(g))
  r <- unname(g$ratio.gw["H1", "vs. unconstrained"])
  expect_equal(r, exp((g$result$gorica[2] - g$result$gorica[1]) / 2))
  expect_true(any(grepl(sprintf("^The order-restricted hypothesis 'H1' is %.3f times more supported than the unconstrained\\.$", r), out)))
  expect_false(any(grepl("> 1 times|< 1 times|= 1", out)))
  # hypothese slechter dan unconstrained: geformuleerd vanuit de unconstrained
  g2 <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 < x2"),
              comparison = "unconstrained", type = "gorica")
  out2 <- capture.output(print(g2))
  expect_true(any(grepl("^The order-restricted hypothesis 'H1' has less support than the unconstrained: the unconstrained is .* times more supported than 'H1'\\.$", out2)))
})

# -----------------------------------------------------------------------------
# W1-N-03: fit_hypothesis() geeft alleen inconsistentie-fouten een Heq-tekst
# -----------------------------------------------------------------------------

test_that("fit_hypothesis: niet-gerelateerde fouten bij Heq worden ongewijzigd doorgegeven", {
  set.seed(1)
  n <- 60
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- 1 + 0.6 * d$x1 + 0.3 * d$x2 + rnorm(n)
  fit <- MASS::rlm(y ~ x1 + x2 + x3, d)  # method = "M": niet ondersteund
  r_heq <- cap_all(goric(fit, hypotheses = list(H1 = "x1 > x2 > x3"), Heq = TRUE))
  r_no  <- cap_all(goric(fit, hypotheses = list(H1 = "x1 > x2 > x3"), Heq = FALSE))
  expect_identical(r_heq$error, r_no$error)
  expect_match(r_heq$error, "bisquare")
  expect_false(grepl("Heq", r_heq$error))
  # een inconsistente Heq krijgt wel de Heq-tekst
  expect_error(restriktor:::fit_hypothesis("Heq", function(...) stop("constraints are inconsistent, no solution!"),
                                           list(constraints = "x1 = 0")),
               "Heq.*cannot be formed.*constraints are inconsistent")
  expect_error(restriktor:::fit_hypothesis("Heq", function(...) stop("something else"),
                                           list(constraints = "x1 = 0")),
               "^something else$")
})

# -----------------------------------------------------------------------------
# W1-N-04: overlap-noot: exacte naam 'unconstrained', geen substring
# -----------------------------------------------------------------------------

test_that("print: overlap-noot en maximum-suffix ook voor een hypothese 'unconstrained1'", {
  est <- c(x1 = 1, x2 = 1, x3 = 1)
  for (nm in c("H1", "unconstrained1")) {
    hyp <- list("x1 > 0", "x2 > 0", "x3 < 0")
    names(hyp) <- c(nm, "H2", "H3")
    r <- cap_all({ g <- goric(est, VCOV = diag(3), hypotheses = hyp, type = "gorica"); capture.output(print(g)) })
    expect_null(r$error, info = nm)
    note <- grep("have equal log-likelihood", r$messages, value = TRUE)
    expect_length(note, 1L)
    expect_match(note, paste0("Hypotheses ('H2' and ", sQuote(nm), "|", sQuote(nm), " and 'H2') have equal log-likelihood"))
    expect_true(any(grepl(paste0(sQuote(nm), " and 'H2' have equal support \\(This relative support reached its maximum"), r$value)), info = nm)
  }
  # zonder overlap: geen fout (lege lijst van combinaties)
  expect_no_error(capture.output(print(goric(est_a, VCOV = VCOV_a, 
                                             hypotheses = list(H1 = "x1 > x2", H2 = "x2 > x3"), type = "gorica"))))
})

# -----------------------------------------------------------------------------
# W1-N-05: Heq = TRUE met een hypothese als constraint-matrix (lijst)
# -----------------------------------------------------------------------------

test_that("goric: Heq = TRUE met een constraint-matrix geeft dezelfde Heq als de tekstversie", {
  Amat <- rbind(c(1, -1, 0), c(0, 1, -1))
  g_mat <- goric(est_b, VCOV = VCOV_b, Heq = TRUE,
                 hypotheses = list(H1 = list(constraints = Amat, rhs = c(0, 0), neq = 0)))
  g_chr <- goric(est_b, VCOV = VCOV_b, Heq = TRUE, hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_identical(g_mat$result$model, c("Heq", "H1", "complement"))
  expect_equal(g_mat$result[, -1], g_chr$result[, -1], tolerance = 1e-6)
  expect_identical(g_mat$objectList$Heq$neq, 2L)
  expect_no_error(capture.output(summary(g_mat, brief = FALSE)))
  # lm-route
  g_lm <- goric(fit_g, Heq = TRUE,
                hypotheses = list(H1 = list(constraints = rbind(c(-1, 1, 0), c(0, -1, 1)), rhs = c(0, 0), neq = 0)))
  g_lm2 <- goric(fit_g, Heq = TRUE, hypotheses = list(H1 = "g1 < g2 < g3"))
  expect_equal(g_lm$result[, -1], g_lm2$result[, -1], tolerance = 1e-6)
  # redundantie: dubbele rij (grootste rhs blijft), rij geimpliceerd door gelijkheid
  h <- restriktor:::goric_heq_constraints(est_b, list(constraints = rbind(c(1, 0, 0), c(1, 0, 0), c(0, 1, 0)),
                                                       rhs = c(0.2, 0.5, 0), neq = 0))
  expect_equal(h, list(constraints = rbind(c(1, 0, 0), c(0, 1, 0)), rhs = c(0.5, 0), neq = 2L))
  expect_null(restriktor:::goric_heq_constraints(est_b, list(constraints = rbind(c(1, 0, 0), c(1, 0, 0)),
                                                              rhs = c(0.5, 0.2), neq = 1)))
  # alleen gelijkheden: Heq genegeerd met melding
  r <- cap_all(goric(est_b, VCOV = VCOV_b, Heq = TRUE,
                     hypotheses = list(H1 = list(constraints = rbind(c(1, -1, 0)), rhs = 0, neq = 1))))
  expect_null(r$error)
  expect_true(any(grepl("'Heq' argument is ignored", r$messages)))
  # range-restrictie en inconsistente Heq: duidelijke fouten
  expect_error(goric(est_b, VCOV = VCOV_b, Heq = TRUE,
                     hypotheses = list(H1 = list(constraints = rbind(c(1, 0, 0), c(-1, 0, 0)), rhs = c(0.2, -0.5), neq = 0))),
               "Heq.*range restriction")
  expect_error(goric(est_b, VCOV = VCOV_b, Heq = TRUE,
                     hypotheses = list(H1 = list(constraints = rbind(c(1, -1, 0), c(0, 1, 0), c(1, 0, 0)),
                                                 rhs = c(0, 0.1, 0.05), neq = 0))),
               "Heq.*cannot be formed.*constraint matrix")
})

# -----------------------------------------------------------------------------
# W1-N-06: glm en rlm: geen partial matching van $fitted
# -----------------------------------------------------------------------------

test_that("goric glm/rlm: geen partial-match warnings uit restriktor ($fitted)", {
  op <- options(warnPartialMatchDollar = TRUE, warnPartialMatchArgs = TRUE,
                warnPartialMatchAttr = TRUE)
  on.exit(options(op), add = TRUE)
  set.seed(5)
  n <- 50
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$yb <- rbinom(n, 1, plogis(d$x1))
  d$y <- d$x1 + rnorm(n)
  fit_glm <- glm(yb ~ x1 + x2, binomial, d)
  fit_rlm <- MASS::rlm(y ~ x1 + x2, d, method = "MM")
  for (cmp in c("complement", "unconstrained", "none")) {
    r <- cap_all({
      g <- goric(fit_glm, hypotheses = list(H1 = "x1 > x2"), comparison = cmp)
      capture.output(print(g)); capture.output(summary(g, brief = FALSE)); g
    })
    expect_null(r$error, info = paste("glm", cmp))
    expect_length(r$warnings, 0L)
    r <- cap_all({
      g <- goric(fit_rlm, hypotheses = list(H1 = "x1 > x2"), comparison = cmp)
      capture.output(print(g)); capture.output(summary(g, brief = FALSE)); g
    })
    expect_null(r$error, info = paste("rlm", cmp))
    # (MASS::vcov.rlm zelf geeft 'corr' -> 'correlation'; dat is geen restriktor-code)
    expect_false(any(grepl("fitted", r$warnings)))
  }
})

# -----------------------------------------------------------------------------
# W1-N-07: print/summary crashen niet op NaN-gewichten
# -----------------------------------------------------------------------------

test_that("print/summary: alle IC-gewichten NaN geeft een noot in plaats van een fout", {
  g <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2 > x3", H2 = "x1 > x2"),
             comparison = "none", type = "gorica")
  g$result$gorica.weights[] <- NaN
  r <- cap_all(capture.output(print(g)))
  expect_null(r$error)
  expect_true(any(grepl("IC weights are not available", r$messages)))
  g2 <- goric(est_a, VCOV = VCOV_a, hypotheses = list(H1 = "x1 > x2 > x3"), type = "gorica")
  g2$result$gorica.weights[] <- NaN
  r2 <- cap_all(capture.output(summary(g2)))
  expect_null(r2$error)
  expect_true(any(grepl("IC weights are not available", r2$messages)))
})

# -----------------------------------------------------------------------------
# W4: residuele bevindingen van de her-verificatie
# -----------------------------------------------------------------------------

test_that("W4-N-01: type in hoofdletters wordt niet stil naar gorica omgezet (numeric/lavaan-route)", {
  est <- c(a = 1, b = 2); V <- diag(2) * .1
  r_low <- suppressMessages(goric(est, VCOV = V, hypotheses = list(H1 = "a > b"),
                                  type = "goricac", sample_nobs = 10))
  for (ty in c("GORICAC", "Goricac")) {
    r_up <- suppressMessages(goric(est, VCOV = V, hypotheses = list(H1 = "a > b"),
                                   type = ty, sample_nobs = 10))
    expect_equal(r_up$type, "goricac")
    expect_equal(r_up$result$penalty, r_low$result$penalty)
    expect_true("goricac" %in% names(r_up$result))
  }
  # onafhankelijke referentie: PT_H1 = 0.5*N*2/(N-3) + 0.5*N*3/(N-4) - 1 (N = 10, p = 2)
  expect_equal(r_low$result$penalty[1], 0.5 * (10 * 2 / 7) + 0.5 * (10 * 3 / 6) - 1,
               tolerance = 1e-8)
  # de controle op een te kleine N geldt ook met hoofdletters
  expect_error(suppressMessages(goric(est, VCOV = V, hypotheses = list(H1 = "a > b"),
                                      type = "GORICAC", sample_nobs = 3)),
               "sample size is too small")
})

test_that("W4-N-02: summary.restriktor met goricc en te kleine N geeft een fout, geen Inf", {
  set.seed(3)
  n <- 5
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n)); d$y <- 1 + d$x1 + rnorm(n)
  fit <- lm(y ~ x1 + x2, data = d)
  r <- restriktor(fit, constraints = "x1 > x2")
  expect_error(summary(r, goric = "goricc"), "sample size is too small")
  # bij N - p - 2 = 1 werkt het wel en is de penalty eindig
  n <- 6
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n)); d$y <- 1 + d$x1 + rnorm(n)
  fit6 <- lm(y ~ x1 + x2, data = d)
  s6 <- summary(restriktor(fit6, constraints = "x1 > x2"), goric = "goricc")
  expect_true(is.finite(attr(s6$goric, "penalty")))
})

test_that("W4-N-04: summary(brief = FALSE) met Heq labelt het complement als 'not H1'", {
  set.seed(4)
  d <- data.frame(x1 = rnorm(40), x2 = rnorm(40)); d$y <- 1 + d$x1 + rnorm(40)
  fit <- lm(y ~ x1 + x2, data = d)
  g <- suppressMessages(goric(fit, hypotheses = list(H1 = "x1 > x2"), Heq = TRUE))
  out <- capture.output(summary(g, brief = FALSE))
  expect_true(any(grepl("not H1", out)))
  expect_false(any(grepl("not Heq", out)))
})
