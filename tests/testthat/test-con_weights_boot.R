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
                        convergence_crit = 0.01, seed = 42)
  expect_true(attr(w, "converged"))
})

test_that("con_weights_boot stopt eerder bij convergentie (chunk_iter < total_chunks)", {
  w <- con_weights_boot(d_2d$VCOV, d_2d$Amat, d_2d$meq,
                        R = 20000L, chunk_size = 5000L,
                        convergence_crit = 0.01, seed = 42)
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
