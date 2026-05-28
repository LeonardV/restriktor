# =============================================================================
# Tests: print() en summary() methoden voor evSyn - alle input types
# =============================================================================

# --- Gedeelde testdata ---

# Estimates + VCOV voor 3 studies
est1 <- c(x1 = 5.0, x2 = 3.0)
est2 <- c(x1 = 4.5, x2 = 2.5)
est3 <- c(x1 = 6.0, x2 = 2.0)
VCOV1 <- diag(c(0.10, 0.10))
VCOV2 <- diag(c(0.15, 0.15))
VCOV3 <- diag(c(0.12, 0.12))

H1 <- list(H1 = "x1 > x2")

# Goric objecten
g1 <- goric(est1, VCOV = VCOV1, type = "gorica", hypotheses = H1)
g2 <- goric(est2, VCOV = VCOV2, type = "gorica", hypotheses = H1)
g3 <- goric(est3, VCOV = VCOV3, type = "gorica", hypotheses = H1)

# LL en PT uit de goric objecten
LL_list <- list(g1$result$loglik, g2$result$loglik, g3$result$loglik)
PT_list <- list(g1$result$penalty, g2$result$penalty, g3$result$penalty)

# IC values (gorica waarden)
IC_list <- list(g1$result$gorica, g2$result$gorica, g3$result$gorica)

# IC weights (gorica weights, summen tot 1)
W_list <- list(g1$result$gorica.weights, g2$result$gorica.weights, g3$result$gorica.weights)

# IC ratios (elke vector eindigt met 1 = genormaliseerd t.o.v. laatste hypothese)
R_list <- lapply(W_list, function(w) w / w[length(w)])


# =============================================================================
# 1. evSyn_gorica (input: list van goric objecten)
# =============================================================================

test_that("evSyn_gorica: print() werkt zonder fout", {
  es <- evSyn(object = list(g1, g2, g3))
  expect_output(print(es))
})

test_that("evSyn_gorica: summary() werkt zonder fout", {
  es <- evSyn(object = list(g1, g2, g3))
  expect_output(print(summary(es)))
})

test_that("evSyn_gorica: print() toont Evidence Synthesis", {
  es <- evSyn(object = list(g1, g2, g3))
  out <- capture.output(print(es))
  expect_true(any(grepl("Evidence Synthesis|evidence synthesis", out)))
})

test_that("evSyn_gorica: summary() toont study-specific en cumulative resultaten", {
  es <- evSyn(object = list(g1, g2, g3))
  out <- capture.output(print(summary(es)))
  expect_true(any(grepl("Study-specific|study-specific", out)))
  expect_true(any(grepl("Cumulative|cumulative", out)))
  expect_true(any(grepl("Final", out)))
})

test_that("evSyn_gorica: type_ev = 'equal' werkt met print en summary", {
  es <- evSyn(object = list(g1, g2, g3), type_ev = "equal")
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn_gorica: type_ev = 'average' werkt met print en summary", {
  es <- evSyn(object = list(g1, g2, g3), type_ev = "average")
  expect_output(print(es))
  expect_output(print(summary(es)))
})


# =============================================================================
# 2. evSyn_est (input: estimates + VCOV)
# =============================================================================

test_that("evSyn_est: print() werkt zonder fout", {
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H1)
  expect_output(print(es))
})

test_that("evSyn_est: summary() werkt zonder fout", {
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H1)
  expect_output(print(summary(es)))
})

test_that("evSyn_est: print() toont Evidence Synthesis", {
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H1)
  out <- capture.output(print(es))
  expect_true(any(grepl("Evidence Synthesis|evidence synthesis", out)))
})

test_that("evSyn_est: summary() toont study-specific en cumulative resultaten", {
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H1)
  out <- capture.output(print(summary(es)))
  expect_true(any(grepl("Study-specific|study-specific", out)))
  expect_true(any(grepl("Cumulative|cumulative", out)))
  expect_true(any(grepl("Final", out)))
})

test_that("evSyn_est: type_ev = 'equal' werkt met print en summary", {
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H1, type_ev = "equal")
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn_est: type_ev = 'average' werkt met print en summary", {
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H1, type_ev = "average")
  expect_output(print(es))
  expect_output(print(summary(es)))
})


# =============================================================================
# 3. evSyn_LL (input: log-likelihood + penalty)
# =============================================================================

test_that("evSyn_LL: print() werkt zonder fout", {
  es <- evSyn(object = LL_list, PT = PT_list)
  expect_output(print(es))
})

test_that("evSyn_LL: summary() werkt zonder fout", {
  es <- evSyn(object = LL_list, PT = PT_list)
  expect_output(print(summary(es)))
})

test_that("evSyn_LL: print() toont Evidence Synthesis", {
  es <- evSyn(object = LL_list, PT = PT_list)
  out <- capture.output(print(es))
  expect_true(any(grepl("Evidence Synthesis|evidence synthesis", out)))
})

test_that("evSyn_LL: summary() toont study-specific en cumulative resultaten", {
  es <- evSyn(object = LL_list, PT = PT_list)
  out <- capture.output(print(summary(es)))
  expect_true(any(grepl("Study-specific|study-specific", out)))
  expect_true(any(grepl("Cumulative|cumulative", out)))
  expect_true(any(grepl("Final", out)))
})

test_that("evSyn_LL: type_ev = 'equal' werkt met print en summary", {
  es <- evSyn(object = LL_list, PT = PT_list, type_ev = "equal")
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn_LL: type_ev = 'average' werkt met print en summary", {
  es <- evSyn(object = LL_list, PT = PT_list, type_ev = "average")
  expect_output(print(es))
  expect_output(print(summary(es)))
})


# =============================================================================
# 4. evSyn_ICvalues (input: IC waarden)
# =============================================================================

test_that("evSyn_ICvalues: print() werkt zonder fout", {
  es <- evSyn(object = IC_list)
  expect_output(print(es))
})

test_that("evSyn_ICvalues: summary() werkt zonder fout", {
  es <- evSyn(object = IC_list)
  expect_output(print(summary(es)))
})

test_that("evSyn_ICvalues: print() toont Evidence Synthesis", {
  es <- evSyn(object = IC_list)
  out <- capture.output(print(es))
  expect_true(any(grepl("Evidence Synthesis|evidence synthesis", out)))
})

test_that("evSyn_ICvalues: summary() toont cumulative resultaten", {
  es <- evSyn(object = IC_list)
  out <- capture.output(print(summary(es)))
  expect_true(any(grepl("Cumulative|cumulative|Final", out)))
})

test_that("evSyn_ICvalues: type_ev = 'average' werkt met print en summary", {
  es <- evSyn(object = IC_list, type_ev = "average")
  expect_output(print(es))
  expect_output(print(summary(es)))
})


# =============================================================================
# 5. evSyn_ICweights (input: IC gewichten, sommeert tot 1)
# =============================================================================

test_that("evSyn_ICweights: print() werkt zonder fout", {
  es <- evSyn(object = W_list)
  expect_output(print(es))
})

test_that("evSyn_ICweights: summary() werkt zonder fout", {
  es <- evSyn(object = W_list)
  expect_output(print(summary(es)))
})

test_that("evSyn_ICweights: print() toont Evidence Synthesis", {
  es <- evSyn(object = W_list)
  out <- capture.output(print(es))
  expect_true(any(grepl("Evidence Synthesis|evidence synthesis", out)))
})

test_that("evSyn_ICweights: summary() toont cumulative resultaten", {
  es <- evSyn(object = W_list)
  out <- capture.output(print(summary(es)))
  expect_true(any(grepl("Cumulative|cumulative|Final", out)))
})

test_that("evSyn_ICweights: type_ev = 'average' werkt met print en summary", {
  es <- evSyn(object = W_list, type_ev = "average")
  expect_output(print(es))
  expect_output(print(summary(es)))
})


# =============================================================================
# 6. evSyn_ICratios (input: ratio IC gewichten, elke vector eindigt op 1)
# =============================================================================

test_that("evSyn_ICratios: print() werkt zonder fout", {
  es <- evSyn(object = R_list)
  expect_output(print(es))
})

test_that("evSyn_ICratios: summary() werkt zonder fout", {
  es <- evSyn(object = R_list)
  expect_output(print(summary(es)))
})

test_that("evSyn_ICratios: print() toont Evidence Synthesis", {
  es <- evSyn(object = R_list)
  out <- capture.output(print(es))
  expect_true(any(grepl("Evidence Synthesis|evidence synthesis", out)))
})

test_that("evSyn_ICratios: summary() toont cumulative resultaten", {
  es <- evSyn(object = R_list)
  out <- capture.output(print(summary(es)))
  expect_true(any(grepl("Cumulative|cumulative|Final", out)))
})

test_that("evSyn_ICratios: type_ev = 'average' werkt met print en summary", {
  es <- evSyn(object = R_list, type_ev = "average")
  expect_output(print(es))
  expect_output(print(summary(es)))
})


# =============================================================================
# 7. Expliciete input_type parameter 
# =============================================================================

test_that("evSyn met input_type = 'gorica': print en summary werken", {
  es <- evSyn(object = list(g1, g2, g3), input_type = "gorica")
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn met input_type = 'est_vcov': print en summary werken", {
  es <- evSyn(object = list(est1, est2, est3), input_type = "est_vcov",
              VCOV = list(VCOV1, VCOV2, VCOV3), hypotheses = H1)
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn met input_type = 'll_pt': print en summary werken", {
  es <- evSyn(object = LL_list, input_type = "ll_pt", PT = PT_list)
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn met input_type = 'icvalues': print en summary werken", {
  es <- evSyn(object = IC_list, input_type = "icvalues")
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn met input_type = 'icweights': print en summary werken", {
  es <- evSyn(object = W_list, input_type = "icweights")
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn met input_type = 'icratios': print en summary werken", {
  es <- evSyn(object = R_list, input_type = "icratios")
  expect_output(print(es))
  expect_output(print(summary(es)))
})


# =============================================================================
# 8. Correcte class toewijzing
# =============================================================================

test_that("evSyn_gorica: class bevat 'evSyn' en 'evSyn_gorica'", {
  es <- evSyn(object = list(g1, g2, g3))
  expect_true("evSyn" %in% class(es))
  expect_true("evSyn_gorica" %in% class(es))
})

test_that("evSyn_est: class bevat 'evSyn' en 'evSyn_est'", {
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H1)
  expect_true("evSyn" %in% class(es))
  expect_true("evSyn_est" %in% class(es))
})

test_that("evSyn_LL: class bevat 'evSyn' en 'evSyn_LL'", {
  es <- evSyn(object = LL_list, PT = PT_list)
  expect_true("evSyn" %in% class(es))
  expect_true("evSyn_LL" %in% class(es))
})

test_that("evSyn_ICvalues: class bevat 'evSyn' en 'evSyn_ICvalues'", {
  es <- evSyn(object = IC_list)
  expect_true("evSyn" %in% class(es))
  expect_true("evSyn_ICvalues" %in% class(es))
})

test_that("evSyn_ICweights: class bevat 'evSyn' en 'evSyn_ICweights'", {
  es <- evSyn(object = W_list)
  expect_true("evSyn" %in% class(es))
  expect_true("evSyn_ICweights" %in% class(es))
})

test_that("evSyn_ICratios: class bevat 'evSyn' en 'evSyn_ICratios'", {
  es <- evSyn(object = R_list)
  expect_true("evSyn" %in% class(es))
  expect_true("evSyn_ICratios" %in% class(es))
})


# =============================================================================
# 9. summary retourneert correct object
# =============================================================================

test_that("summary.evSyn retourneert 'summary.evSyn' class", {
  es <- evSyn(object = list(g1, g2, g3))
  s <- summary(es)
  expect_true(inherits(s, "summary.evSyn"))
})

test_that("summary.evSyn bevat Final_Cumulative_results", {
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H1)
  s <- summary(es)
  expect_true(!is.null(s$Final_Cumulative_results))
})

test_that("summary.evSyn bevat Final_ratio_GORICA_weights", {
  es <- evSyn(object = list(g1, g2, g3))
  s <- summary(es)
  expect_true(!is.null(s$Final_ratio_GORICA_weights))
})


# =============================================================================
# 10. Meerdere hypothesen
# =============================================================================

test_that("evSyn_est: meerdere hypothesen - print en summary werken", {
  H2 <- list(H1 = "x1 > x2", H2 = "x1 < x2")
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H2)
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn_gorica: meerdere hypothesen - print en summary werken", {
  H2 <- list(H1 = "x1 > x2", H2 = "x1 < x2")
  g1b <- goric(est1, VCOV = VCOV1, type = "gorica", hypotheses = H2)
  g2b <- goric(est2, VCOV = VCOV2, type = "gorica", hypotheses = H2)
  g3b <- goric(est3, VCOV = VCOV3, type = "gorica", hypotheses = H2)
  es <- evSyn(object = list(g1b, g2b, g3b))
  expect_output(print(es))
  expect_output(print(summary(es)))
})


# =============================================================================
# 11. Edge case: 2 studies (minimum)
# =============================================================================

test_that("evSyn_gorica: 2 studies - print en summary werken", {
  es <- evSyn(object = list(g1, g2))
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn_LL: 2 studies - print en summary werken", {
  es <- evSyn(object = LL_list[1:2], PT = PT_list[1:2])
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn_ICvalues: 2 studies - print en summary werken", {
  es <- evSyn(object = IC_list[1:2])
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn_ICweights: 2 studies - print en summary werken", {
  es <- evSyn(object = W_list[1:2])
  expect_output(print(es))
  expect_output(print(summary(es)))
})

test_that("evSyn_ICratios: 2 studies - print en summary werken", {
  es <- evSyn(object = R_list[1:2])
  expect_output(print(es))
  expect_output(print(summary(es)))
})


# =============================================================================
# 12. study_names en hypo_names parameters
# =============================================================================

test_that("evSyn_gorica: custom study_names werken met print en summary", {
  es <- evSyn(object = list(g1, g2, g3), 
              study_names = c("Studie A", "Studie B", "Studie C"))
  out <- capture.output(print(summary(es)))
  expect_true(any(grepl("Studie A|Studie B|Studie C", out)))
})

test_that("evSyn_LL: custom hypo_names werken met print en summary", {
  es <- evSyn(object = LL_list, PT = PT_list, 
              hypo_names = c("Hypothese1", "Unconstrained"))
  out <- capture.output(print(summary(es)))
  expect_true(any(grepl("Hypothese1", out)))
})
