# =============================================================================
# Tests: penalty_factor in evSyn() en leave1studyout()
# =============================================================================

# --- Gedeelde testdata: LL en PT voor 4 studies ---
set.seed(1)
S <- 4
LL_list <- lapply(1:S, function(s) {
  c(H1 = rnorm(1, -50, 2), H2 = rnorm(1, -52, 2), Hu = rnorm(1, -49, 2))
})
PT_list <- lapply(1:S, function(s) c(H1 = 1.5, H2 = 2, Hu = 3))

LL_m <- do.call(rbind, LL_list)
PT_m <- do.call(rbind, PT_list)


test_that("evSyn_LL: default penalty_factor = 2 en wordt opgeslagen", {
  es <- evSyn(LL_list, PT = PT_list, type_ev = "equal")
  expect_identical(es$penalty_factor, 2)
  expect_equal(unname(es$GORICA_m), unname(-2 * LL_m + 2 * PT_m))
})

test_that("evSyn_LL: penalty_factor werkt door in GORICA_m en cumulatieve waarden", {
  pf <- 3
  es <- evSyn(LL_list, PT = PT_list, type_ev = "equal", penalty_factor = pf)
  expect_identical(es$penalty_factor, pf)
  expect_equal(unname(es$GORICA_m), unname(-2 * LL_m + pf * PT_m))
  # finale cumulatieve equal-evidence GORICA
  expect_equal(unname(es$Cumulative_GORICA[S, ]),
               unname(-2 * colSums(LL_m) + pf * colMeans(PT_m)))
})

test_that("evSyn_LL: ongeldige penalty_factor geeft fout", {
  expect_error(evSyn(LL_list, PT = PT_list, penalty_factor = -1),
               "penalty_factor")
  expect_error(evSyn(LL_list, PT = PT_list, penalty_factor = c(1, 2)),
               "penalty_factor")
})

test_that("leave1studyout: 'equal' gebruikt de penalty_factor uit het evSyn object", {
  pf <- 3
  es <- evSyn(LL_list, PT = PT_list, type_ev = "equal", penalty_factor = pf)
  l1o <- leave1studyout(es)
  for (s in 1:S) {
    keep <- seq_len(S) != s
    expected <- -2 * colSums(LL_m[keep, , drop = FALSE]) +
      pf * colMeans(PT_m[keep, , drop = FALSE])
    expect_equal(unname(l1o$OverallGorica[s, ]), unname(expected))
  }
})

test_that("leave1studyout: valt terug op penalty_factor = 2 voor oudere objecten", {
  es <- evSyn(LL_list, PT = PT_list, type_ev = "equal")
  es$penalty_factor <- NULL
  l1o <- leave1studyout(es)
  for (s in 1:S) {
    keep <- seq_len(S) != s
    expected <- -2 * colSums(LL_m[keep, , drop = FALSE]) +
      2 * colMeans(PT_m[keep, , drop = FALSE])
    expect_equal(unname(l1o$OverallGorica[s, ]), unname(expected))
  }
})

test_that("leave1studyout: 'added' en 'average' consistent met IC's (incl. penalty_factor)", {
  pf <- 3
  for (tev in c("added", "average")) {
    es <- evSyn(LL_list, PT = PT_list, type_ev = tev, penalty_factor = pf)
    l1o <- leave1studyout(es)
    agg <- if (tev == "added") colSums else colMeans
    for (s in 1:S) {
      keep <- seq_len(S) != s
      expected <- agg(es$GORICA_m[keep, , drop = FALSE])
      expect_equal(unname(l1o$OverallGorica[s, ]), unname(expected))
    }
  }
})

test_that("evSyn_est: penalty_factor werkt door in studie-specifieke en overall IC's", {
  est1 <- c(x1 = 5.0, x2 = 3.0)
  est2 <- c(x1 = 4.5, x2 = 2.5)
  est3 <- c(x1 = 6.0, x2 = 2.0)
  VCOV1 <- diag(c(0.10, 0.10))
  VCOV2 <- diag(c(0.15, 0.15))
  VCOV3 <- diag(c(0.12, 0.12))
  H1 <- list(H1 = "x1 > x2")

  pf <- 3
  es <- evSyn(object = list(est1, est2, est3), VCOV = list(VCOV1, VCOV2, VCOV3),
              hypotheses = H1, type_ev = "equal", penalty_factor = pf)
  expect_identical(es$penalty_factor, pf)
  # per studie: IC = -2*LL + pf*PT
  expect_equal(unname(es$GORICA_m), unname(-2 * es$LL_m + pf * es$PT_m))
  # leave1studyout gebruikt dezelfde penalty_factor
  l1o <- leave1studyout(es)
  for (s in 1:3) {
    keep <- seq_len(3) != s
    expected <- -2 * colSums(es$LL_m[keep, , drop = FALSE]) +
      pf * colMeans(es$PT_m[keep, , drop = FALSE])
    expect_equal(unname(l1o$OverallGorica[s, ]), unname(expected))
  }
})

test_that("evSyn_gorica: penalty_factor wordt uit de goric objecten overgenomen", {
  est1 <- c(x1 = 5.0, x2 = 3.0)
  est2 <- c(x1 = 4.5, x2 = 2.5)
  VCOV1 <- diag(c(0.10, 0.10))
  VCOV2 <- diag(c(0.15, 0.15))
  H1 <- list(H1 = "x1 > x2")

  pf <- 3
  g1 <- goric(est1, VCOV = VCOV1, type = "gorica", hypotheses = H1, penalty_factor = pf)
  g2 <- goric(est2, VCOV = VCOV2, type = "gorica", hypotheses = H1, penalty_factor = pf)

  es <- evSyn(object = list(g1, g2), type_ev = "equal")
  expect_identical(es$penalty_factor, pf)
  expect_equal(unname(es$GORICA_m[1, ]), unname(g1$result$gorica))
  expect_equal(unname(es$GORICA_m[2, ]), unname(g2$result$gorica))

  # verschillende penalty factors over studies is niet toegestaan
  g2b <- goric(est2, VCOV = VCOV2, type = "gorica", hypotheses = H1, penalty_factor = 2)
  expect_error(evSyn(object = list(g1, g2b)), "same penalty_factor")
})
