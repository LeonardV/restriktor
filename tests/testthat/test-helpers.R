# =============================================================================
# Tests: Helper functies
# make_positive_semi_definite, nchoosek, GaussianElimination
# =============================================================================

# --- make_positive_semi_definite ---

test_that("make_positive_semi_definite behoudt positief-definiete matrix", {
  mat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  result <- restriktor:::make_positive_definite(mat)
  
  # Resultaat moet vierkant zijn met dezelfde dimensie
  expect_equal(nrow(result$corr), 2)
  expect_equal(ncol(result$corr), 2)
  
  # Eigenwaarden moeten positief zijn
  ev <- eigen(result$corr)$values
  expect_true(all(ev > 0))
  
  # Resultaat moet dicht bij origineel liggen
  expect_equal(result$corr, mat, tolerance = 1e-10)
})

test_that("make_positive_semi_definite corrigeert near-singular matrix", {
  # Maak een bijna-singuliere matrix
  mat <- matrix(c(1, 0.9999999, 0.9999999, 1), 2, 2)
  result <- restriktor:::make_positive_definite(mat, 
                                                tolerance = 1e-15, 
                                                ridge_constant = 1e-05)
  
  # Moet inverteerbaar zijn na correctie
  expect_no_error(solve(result$corr))
  
  # Eigenwaarden moeten strikt positief zijn
  ev <- eigen(result$corr)$values
  expect_true(all(ev > 0))
})

test_that("make_positive_semi_definite verwerkt identiteitsmatrix correct", {
  mat <- diag(5)
  result <- restriktor:::make_positive_definite(mat)
  expect_equal(result$corr, mat, tolerance = 1e-12)
})

test_that("make_positive_semi_definite past ridging toe bij singuliere matrix", {
  # Rank-deficiente matrix: kolom 3 = kolom 1 + kolom 2
  mat <- matrix(c(1, 0.8, 0.8,
                  0.8, 1, 0.8,
                  0.8, 0.8, 1), 3, 3)
  # Forceer singulariteit
  mat[3, ] <- mat[1, ] + mat[2, ] - mat[3, ]
  mat[, 3] <- mat[, 1] + mat[, 2] - mat[, 3]
  mat <- (mat + t(mat)) / 2
  
  result <- restriktor:::make_positive_definite(mat, ridge_constant = 1e-04)
  
  # Na correctie moet de matrix inverteerbaar zijn
  inv_result <- try(solve(result$corr), silent = TRUE)
  expect_false(inherits(inv_result, "try-error"))
})


# --- nchoosek ---

test_that("nchoosek geeft correcte combinaties voor n=4, k=2", {
  result <- restriktor:::nchoosek(4, 2)
  
  # C(4,2) = 6 combinaties, elk van lengte 2
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 6)
  
  # Elke kolom moet unieke elementen bevatten uit 1:4
  for (j in 1:ncol(result)) {
    expect_true(all(result[, j] %in% 1:4))
    expect_equal(length(unique(result[, j])), 2)
  }
  
  # Alle combinaties moeten uniek zijn
  combos <- apply(result, 2, paste, collapse = "-")
  expect_equal(length(unique(combos)), 6)
})

test_that("nchoosek geeft correcte combinaties voor n=5, k=3", {
  result <- restriktor:::nchoosek(5, 3)
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), choose(5, 3))
})

test_that("nchoosek geeft fout bij ongeldige invoer", {
  expect_error(restriktor:::nchoosek(3, 5), "0 <= k <= n")
  expect_error(restriktor:::nchoosek("a", 2), "non-NA numeric scalars")
  expect_error(restriktor:::nchoosek(NA, 2), "non-NA numeric scalars")
})

test_that("nchoosek werkt voor k=1 (triviale combinaties)", {
  result <- restriktor:::nchoosek(5, 1)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 5)
  expect_equal(as.vector(result), 1:5)
})


# --- GaussianElimination ---

test_that("GaussianElimination lost eenvoudig stelsel op", {
  # 2x + y = 5, x + 3y = 10 => x = 1, y = 3
  A <- matrix(c(2, 1, 1, 3), 2, 2)
  B <- matrix(c(5, 10), 2, 1)
  result <- restriktor:::GaussianElimination(A, B)
  
  # Check rang
  expect_equal(result$rank, 2)
  
  # Oplossing zit in de laatste kolom van E
  solution <- result$E[, ncol(result$E)]
  expect_equal(solution, c(1, 3), tolerance = 1e-10)
})

test_that("GaussianElimination detecteert rang correct", {
  # Rang-deficiente matrix
  A <- matrix(c(1, 2, 2, 4), 2, 2)
  result <- restriktor:::GaussianElimination(A)
  expect_equal(result$rank, 1)
})

test_that("GaussianElimination reduceert identiteit correct", {
  A <- diag(3)
  result <- restriktor:::GaussianElimination(A)
  expect_equal(result$rank, 3)
  expect_equal(result$E, diag(3))
})

test_that("GaussianElimination geeft pivot kolommen terug", {
  A <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), 3, 3)
  result <- restriktor:::GaussianElimination(A)
  expect_equal(result$pivot, c(1, 2, 3))
})

test_that("GaussianElimination verwerkt niet-vierkante matrix", {
  # 2 vergelijkingen, 3 onbekenden
  A <- matrix(c(1, 0, 0, 1, 1, 1), 2, 3)
  result <- restriktor:::GaussianElimination(A)
  expect_equal(result$rank, 2)
})
