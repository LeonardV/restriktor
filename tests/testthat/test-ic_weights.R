# =============================================================================
# Tests: ic_weights (analytisch) en .orthant_prob
# =============================================================================


# --- Gedeelde testdata ---
tol <- 1e-8
ridge <- 1e-05

# Compound symmetry correlatiematrix
make_compound_sym <- function(g, rho = 0.5) {
  mat <- matrix(rho, g, g)
  diag(mat) <- 1
  mat
}

# AR(1) correlatiematrix
make_ar1 <- function(g, phi = 0.7) {
  mat <- matrix(0, g, g)
  for (i in 1:g) for (j in 1:g) mat[i, j] <- phi^abs(i - j)
  mat
}


# =============================================================================
# .orthant_prob
# =============================================================================

test_that(".orthant_prob geeft exact 0.5 voor dimensie 1", {
  sigma <- matrix(1, 1, 1)
  result <- restriktor:::.orthant_prob(sigma, 0, Inf)
  expect_equal(result, 0.5)
})

test_that(".orthant_prob geeft correcte waarde voor dimensie 2", {
  # Voor rho=0: P(X1>0, X2>0) = 1/4
  sigma <- diag(2)
  result <- restriktor:::.orthant_prob(sigma, rep(0, 2), rep(Inf, 2))
  expect_equal(result, 0.25, tolerance = 0.01)
  
  # Voor rho=1: P(X1>0, X2>0) = 1/2
  sigma_pos <- matrix(c(1, 0.999, 0.999, 1), 2, 2)
  result_pos <- restriktor:::.orthant_prob(sigma_pos, rep(0, 2), rep(Inf, 2))
  expect_lt(abs(result_pos - 0.5), 0.01)
})

test_that(".orthant_prob geeft correcte waarde voor dimensie 3 (Sheppard)", {
  # Onafhankelijke variabelen: P(X1>0, X2>0, X3>0) = 1/8
  sigma <- diag(3)
  result <- restriktor:::.orthant_prob(sigma, rep(0, 3), rep(Inf, 3))
  expect_equal(result, 0.125, tolerance = 1e-10)
})


# =============================================================================
# ic_weights: basiseigenschappen
# =============================================================================

test_that("ic_weights geeft gewichten die optellen tot 1", {
  for (g in c(2, 3, 4, 6)) {
    corr <- make_compound_sym(g)
    w <- restriktor:::ic_weights(corr, tolerance = tol, ridge_constant = ridge)
    expect_equal(sum(w), 1, tolerance = 1e-6,
                 info = paste("g =", g, "sum(weights) != 1"))
  }
})

test_that("ic_weights geeft niet-negatieve gewichten", {
  for (g in c(2, 3, 4, 6, 8)) {
    corr <- make_compound_sym(g)
    w <- restriktor:::ic_weights(corr, tolerance = tol, ridge_constant = ridge)
    expect_true(all(w >= -1e-10),
                info = paste("g =", g, "negatieve gewichten gevonden"))
  }
})

test_that("ic_weights heeft correcte namen (g:0)", {
  g <- 4
  corr <- make_compound_sym(g)
  w <- restriktor:::ic_weights(corr, tolerance = tol, ridge_constant = ridge)
  expect_equal(names(w), as.character(g:0))
})

test_that("ic_weights g=1 geeft (0.5, 0.5)", {
  corr <- matrix(1, 1, 1)
  w <- restriktor:::ic_weights(corr, tolerance = tol, ridge_constant = ridge)
  expect_equal(as.numeric(w), c(0.5, 0.5))
})

test_that("ic_weights identiteitsmatrix: symmetrische gewichten", {
  g <- 4
  corr <- diag(g)
  w <- restriktor:::ic_weights(corr, tolerance = tol, ridge_constant = ridge)
  # Voor identiteitsmatrix: w[k] = C(g,k) / 2^g
  expected <- choose(g, g:0) / 2^g
  expect_equal(as.numeric(w), expected, tolerance = 1e-6)
})


# =============================================================================
# ic_weights: foutafhandeling
# =============================================================================

test_that("ic_weights geeft fout bij niet-vierkante matrix", {
  mat <- matrix(1:6, 2, 3)
  expect_error(restriktor:::ic_weights(mat, tolerance = tol,
                                        ridge_constant = ridge),
               "square matrix")
})


# =============================================================================
# con_weights
# =============================================================================

test_that("con_weights werkt zonder gelijkheidsbeperkingen (meq=0)", {
  corr <- make_compound_sym(4)
  w <- restriktor:::con_weights(corr, meq = 0L, tolerance = tol,
                                 ridge_constant = ridge)
  expect_equal(sum(w), 1, tolerance = 1e-6)
  expect_equal(length(w), 5)
})

test_that("con_weights werkt met gelijkheidsbeperkingen (meq>0)", {
  corr <- make_compound_sym(4)
  w <- restriktor:::con_weights(corr, meq = 1L, tolerance = tol,
                                 ridge_constant = ridge)
  # Met 1 gelijkheidsbeperking: werkt op submatrix van dimensie g-meq=3
  expect_equal(sum(w), 1, tolerance = 1e-6)
  expect_equal(length(w), 4)  # g-meq+1 = 4
})


# =============================================================================
# ic_weights: diverse correlatiematrix-types
# =============================================================================

test_that("ic_weights werkt met AR(1) correlatiematrix", {
  corr <- make_ar1(6)
  w <- restriktor:::ic_weights(corr, tolerance = tol, ridge_constant = ridge)
  expect_equal(sum(w), 1, tolerance = 1e-6)
  expect_true(all(w >= -1e-10))
})

test_that("ic_weights werkt met hoge correlatie", {
  g <- 4
  corr <- matrix(0.9, g, g)
  diag(corr) <- 1
  w <- restriktor:::ic_weights(corr, tolerance = tol, ridge_constant = ridge)
  expect_equal(sum(w), 1, tolerance = 1e-6)
  # Bij hoge correlatie: meer gewicht op extremen (alle positief / alle negatief)
  expect_true(w[1] > 0.1)  # P(alle g positief) is substantieel
})
