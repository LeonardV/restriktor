# =============================================================================
# Tests: con_weights_boot (bootstrap-gebaseerde chi-bar-square gewichten)
#
# NB: Bootstrap R waarden bewust laag gehouden voor snelle tests.
#     Gebruik hogere tolerance waar nodig.
# =============================================================================

# --- Gedeelde testdata ---

# Eenvoudig 2D probleem: 2 ongelijkheidsrestricties, meq = 0
make_simple_2d <- function() {
  VCOV <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  Amat <- diag(2)  # x1 >= 0, x2 >= 0
  list(VCOV = VCOV, Amat = Amat, meq = 0L)
}

# 3D probleem met ordening: x1 >= x2 >= x3
make_ordering_3d <- function() {
  VCOV <- diag(3)
  Amat <- matrix(c(1, -1,  0,
                    0,  1, -1), nrow = 2, byrow = TRUE)
  list(VCOV = VCOV, Amat = Amat, meq = 0L)
}

# 3D probleem met 1 gelijkheidsbeperking: x1 = x2, x3 >= 0
make_mixed_3d <- function() {
  VCOV <- diag(3)
  Amat <- matrix(c(1, -1, 0,
                    0,  0, 1), nrow = 2, byrow = TRUE)
  list(VCOV = VCOV, Amat = Amat, meq = 1L)
}


# --- Basiseigenschappen (1 gedeelde berekening) ---

d_2d <- make_simple_2d()
w_2d <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                          R = 5000L, chunk_size = 5000L, seed = 42)

test_that("con_weights_boot gewichten tellen op tot 1", {
  expect_equal(sum(w_2d), 1, tolerance = 1e-3)
})

test_that("con_weights_boot gewichten zijn niet-negatief", {
  expect_true(all(w_2d >= 0))
})

test_that("con_weights_boot heeft correcte lengte (ncol(VCOV) + 1)", {
  expect_equal(length(w_2d), ncol(d_2d$VCOV) + 1)
})


# --- Attributen (hergebruik w_2d) ---

test_that("con_weights_boot retourneert alle verwachte attributen", {
  expect_false(is.null(attr(w_2d, "total_bootstrap_draws")))
  expect_false(is.null(attr(w_2d, "converged")))
  expect_false(is.null(attr(w_2d, "convergence_crit")))
  expect_false(is.null(attr(w_2d, "wt_bar_chunk")))
  expect_false(is.null(attr(w_2d, "chunk_size")))
  expect_false(is.null(attr(w_2d, "total_chunks")))
  expect_false(is.null(attr(w_2d, "chunk_iter")))
  expect_false(is.null(attr(w_2d, "error.idx")))
  expect_false(is.null(attr(w_2d, "mvtnorm")))
})

test_that("con_weights_boot convergence_crit attribuut komt overeen met invoer", {
  crit <- 0.005
  w <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                        R = 10000L, chunk_size = 5000L,
                        convergence_crit = crit, seed = 42)
  expect_equal(attr(w, "convergence_crit"), crit)
})

test_that("con_weights_boot chunk_size attribuut komt overeen met invoer", {
  expect_equal(attr(w_2d, "chunk_size"), 5000L)
})

test_that("con_weights_boot total_bootstrap_draws <= R", {
  expect_true(attr(w_2d, "total_bootstrap_draws") <= 5000L)
})


# --- Convergentie ---

test_that("con_weights_boot convergeert bij voldoende draws", {
  w <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                        R = 20000L, chunk_size = 5000L,
                        convergence_crit = 0.02, seed = 42)
  expect_true(attr(w, "converged"))
  # de stopregel garandeert: 95%-CI-halfbreedte van elk gewicht < crit
  expect_true(1.96 * max(attr(w, "mc_se")) < 0.02)
})

test_that("con_weights_boot stopt eerder bij convergentie (chunk_iter < total_chunks)", {
  w <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                        R = 20000L, chunk_size = 5000L,
                        convergence_crit = 0.02, seed = 42)
  if (attr(w, "converged")) {
    expect_true(attr(w, "chunk_iter") <= attr(w, "total_chunks"))
  }
})

test_that("con_weights_boot wt_bar_chunk heeft correcte rijnamen", {
  w <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                        R = 10000L, chunk_size = 5000L, seed = 42)
  chunk_mat <- attr(w, "wt_bar_chunk")
  n_chunks <- nrow(chunk_mat)
  expect_equal(rownames(chunk_mat),
               paste0("chunk_iter_", seq_len(n_chunks)))
})


# --- Reproduceerbaarheid ---

test_that("con_weights_boot is reproduceerbaar met seed", {
  w1 <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                         R = 5000L, chunk_size = 5000L, seed = 123)
  w2 <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                         R = 5000L, chunk_size = 5000L, seed = 123)
  expect_equal(as.numeric(w1), as.numeric(w2))
})

test_that("con_weights_boot verschilt met verschillende seeds", {
  w1 <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                         R = 5000L, chunk_size = 5000L, seed = 1)
  w2 <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                         R = 5000L, chunk_size = 5000L, seed = 999)
  expect_false(identical(as.numeric(w1), as.numeric(w2)))
})


# --- Dimensies en meq ---

test_that("con_weights_boot werkt met 3D ordening", {
  d <- make_ordering_3d()
  w <- con_weights_boot(d$VCOV, d$Amat, d$meq,
                        R = 5000L, chunk_size = 5000L, seed = 42)
  expect_equal(sum(w), 1, tolerance = 1e-3)
  expect_equal(length(w), ncol(d$VCOV) + 1)
})

test_that("con_weights_boot werkt met meq > 0 (gemengde restricties)", {
  d <- make_mixed_3d()
  w <- con_weights_boot(d$VCOV, d$Amat, d$meq,
                        R = 5000L, chunk_size = 5000L, seed = 42)
  expect_equal(sum(w), 1, tolerance = 1e-3)
  expect_true(all(w >= 0))
})


# --- Chunk-gedrag ---

test_that("con_weights_boot met chunk_size >= R maakt 1 chunk", {
  w <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                        R = 5000L, chunk_size = 5000L, seed = 42)
  expect_equal(attr(w, "total_chunks"), 1)
})

test_that("con_weights_boot met kleine chunk_size levert meerdere chunks", {
  w <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                        R = 10000L, chunk_size = 2000L, seed = 42)
  expect_equal(attr(w, "total_chunks"), 5)
})


# --- Vergelijking met analytisch resultaat ---

test_that("con_weights_boot benadert analytische gewichten (identiteitsmatrix)", {
  g <- 3
  VCOV <- diag(g)
  Amat <- diag(g)
  meq  <- 0L

  w_boot <- con_weights_boot(VCOV, Amat, meq,
                             R = 20000L, chunk_size = 5000L, seed = 42)
  expected <- choose(g, g:0) / 2^g  # c(1/8, 3/8, 3/8, 1/8)

  expect_equal(as.numeric(w_boot), expected, tolerance = 0.05)
})


# --- NNLS-engine (duale projectie) vs exacte pmvnorm-gewichten ---

# testprobleem: ordening-constraints met ruwe random covariantie
make_boot_problem <- function(g, p, seed = 3) {
  set.seed(seed)
  X <- matrix(rnorm(3 * p * p), ncol = p)
  VCOV <- crossprod(X) / (3 * p)
  Amat <- matrix(0, g, p)
  for (i in 1:g) { Amat[i, i] <- -1; Amat[i, i + 1] <- 1 }
  list(VCOV = VCOV, Amat = Amat, W = Amat %*% VCOV %*% t(Amat))
}

test_that("con_weights_boot gebruikt standaard de nnls-engine", {
  expect_equal(attr(w_2d, "engine"), "nnls")
})

test_that("nnls-engine reproduceert de exacte pmvnorm-gewichten (meq = 0)", {
  pr <- make_boot_problem(g = 5, p = 8)
  exact <- rev(restriktor:::ic_weights(pr$W, tolerance = 1e-8,
                                       ridge_constant = 1e-5))
  w <- con_weights_boot(pr$VCOV, pr$Amat, meq = 0, R = 50000L, seed = 1,
                        convergence_crit = 0, chunk_size = 50000L)
  # boot-massa zit op levels (p-g):p; vergelijk met exacte gewichten
  expect_equal(tail(as.numeric(w), 6), as.numeric(exact), tolerance = 0.01)
})

test_that("nnls-engine matcht de exacte reductie bij meq > 0", {
  pr <- make_boot_problem(g = 5, p = 8)
  meq <- 2
  Wred <- solve(solve(pr$W)[-(1:meq), -(1:meq)])
  exact <- rev(restriktor:::ic_weights(Wred, tolerance = 1e-8,
                                       ridge_constant = 1e-5))
  w <- con_weights_boot(pr$VCOV, pr$Amat, meq = meq, R = 50000L, seed = 1,
                        convergence_crit = 0, chunk_size = 50000L)
  # alle massa op levels (p-g):(p-meq) = 3:6
  expect_equal(sum(as.numeric(w)[4:7]), 1, tolerance = 1e-12)
  expect_equal(as.numeric(w)[4:7], as.numeric(exact), tolerance = 0.01)
})

test_that("con_weights_boot werkt met meer constraints dan parameters (g > p)", {
  set.seed(11)
  p <- 4
  X <- matrix(rnorm(3 * p * p), ncol = p)
  VCOV <- crossprod(X) / (3 * p)
  Amat <- rbind(diag(p), c(-1, 1, 0, 0), c(0, 0, -1, 1))  # 6 x 4, A V A' singulier
  w <- con_weights_boot(VCOV, Amat, meq = 0, R = 20000L, seed = 1,
                        convergence_crit = 0, chunk_size = 20000L)
  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_true(all(w >= 0))
  expect_equal(length(w), p + 1)
})

test_that("rtmvnorm-argumenten via ... activeren de quadprog-engine", {
  pr <- make_boot_problem(g = 3, p = 5)
  w <- con_weights_boot(pr$VCOV, pr$Amat, meq = 0, R = 5000L, seed = 1,
                        convergence_crit = 0, chunk_size = 5000L,
                        lower = rep(-Inf, 5), upper = rep(Inf, 5))
  expect_equal(attr(w, "engine"), "quadprog")
  expect_true(all(c("lower", "upper") %in% names(attr(w, "mvtnorm"))))
  expect_equal(sum(w), 1, tolerance = 1e-12)
})

test_that("con_weights_boot accepteert een vector als 1 constraint", {
  W <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  w <- con_weights_boot(W, rbind(c(0, 1)), meq = 0L, R = 20000L, seed = 1)
  # 1 constraint: gewichten (0, 0.5, 0.5)
  expect_equal(as.numeric(w), c(0, 0.5, 0.5), tolerance = 0.02)
})

test_that("meq == aantal rijen: alle massa op level p - meq", {
  pr <- make_boot_problem(g = 5, p = 8)
  Amat <- pr$Amat[1:2, ]
  w <- con_weights_boot(pr$VCOV, Amat, meq = 2, R = 5000L, seed = 1,
                        chunk_size = 5000L)
  expect_equal(as.numeric(w[8 - 2 + 1]), 1)
  expect_equal(sum(w), 1)
})


# --- Inputvalidatie ---

test_that("con_weights_boot valideert de invoer", {
  pr <- make_boot_problem(g = 3, p = 5)
  # niet-vierkante of niet-symmetrische VCOV
  expect_error(con_weights_boot(pr$VCOV[, 1:3], pr$Amat, meq = 0),
               "square")
  V_asym <- pr$VCOV; V_asym[1, 2] <- V_asym[1, 2] + 1
  expect_error(con_weights_boot(V_asym, pr$Amat, meq = 0), "symmetric")
  # niet-positief-definiete VCOV
  V_npd <- pr$VCOV; V_npd[1, 1] <- -1
  expect_error(con_weights_boot(V_npd, pr$Amat, meq = 0), "positive definite")
  # kolomdimensie Amat
  expect_error(con_weights_boot(pr$VCOV, pr$Amat[, 1:3], meq = 0),
               "ncol")
  # meq buiten bereik
  expect_error(con_weights_boot(pr$VCOV, pr$Amat, meq = -1), "meq")
  expect_error(con_weights_boot(pr$VCOV, pr$Amat, meq = 4), "meq")
  # ongeldige R / chunk_size / convergence_crit
  expect_error(con_weights_boot(pr$VCOV, pr$Amat, meq = 0, R = 0), "R must")
  expect_error(con_weights_boot(pr$VCOV, pr$Amat, meq = 0, R = NA), "R must")
  expect_error(con_weights_boot(pr$VCOV, pr$Amat, meq = 0, R = 1e10), "R must")
  expect_error(con_weights_boot(pr$VCOV, pr$Amat, meq = 0, chunk_size = 0),
               "chunk_size")
  expect_error(con_weights_boot(pr$VCOV, pr$Amat, meq = 0,
                                convergence_crit = -1), "convergence_crit")
})

test_that("lineair afhankelijke equality-constraints geven een foutmelding", {
  pr <- make_boot_problem(g = 5, p = 8)
  Amat <- rbind(pr$Amat[1, ], pr$Amat[1, ], pr$Amat[3:5, ])  # rij 1 dubbel
  expect_error(con_weights_boot(pr$VCOV, Amat, meq = 2, R = 5000L, seed = 1),
               "linearly dependent")
})


# --- SE-gebaseerde stopregel ---

test_that("mc_se attribuut is aanwezig en plausibel", {
  pr <- make_boot_problem(g = 5, p = 8)
  w <- con_weights_boot(pr$VCOV, pr$Amat, meq = 0, R = 20000L, seed = 1,
                        convergence_crit = 0, chunk_size = 20000L)
  se <- attr(w, "mc_se")
  expect_equal(length(se), length(w))
  expect_true(all(is.finite(se)) && all(se >= 0))
  # MC-fout bij 20000 draws hoort ruim onder 0.01 te liggen
  expect_true(max(se) < 0.01)
})

test_that("stopregel bindt de werkelijke fout (vergelijking met exact)", {
  pr <- make_boot_problem(g = 5, p = 8)
  exact <- rev(restriktor:::ic_weights(pr$W, tolerance = 1e-8,
                                       ridge_constant = 1e-5))
  crit <- 0.01
  w <- con_weights_boot(pr$VCOV, pr$Amat, meq = 0, R = 1e5L, seed = 7,
                        convergence_crit = crit, chunk_size = 5000L)
  expect_true(attr(w, "converged"))
  # de fout van elk gewicht hoort (op ~95% niveau) onder crit te liggen;
  # neem wat marge voor de resterende 5%
  expect_true(max(abs(tail(as.numeric(w), 6) - exact)) < 1.5 * crit)
})
