# =============================================================================
# Tests: calculate_IC_weights (geexporteerde functie)
# =============================================================================

test_that("calculate_IC_weights geeft gewichten die optellen tot 1", {
  IC <- c(100, 102, 105)
  result <- restriktor:::calculate_IC_weights(IC)
  
  expect_equal(sum(result$IC_weights), 1, tolerance = 1e-10)
})

test_that("calculate_IC_weights laagste IC krijgt hoogste gewicht", {
  IC <- c(100, 102, 110)
  result <- restriktor:::calculate_IC_weights(IC)
  
  # H1 heeft laagste IC, dus hoogste gewicht
  expect_true(result$IC_weights[1] > result$IC_weights[2])
  expect_true(result$IC_weights[2] > result$IC_weights[3])
})

test_that("calculate_IC_weights gelijke IC geeft gelijke gewichten", {
  IC <- c(100, 100, 100)
  result <- restriktor:::calculate_IC_weights(IC)
  
  expect_equal(as.numeric(result$IC_weights), rep(1/3, 3), tolerance = 1e-10)
})

test_that("calculate_IC_weights correcte standaardnamen", {
  IC <- c(50, 55)
  result <- restriktor:::calculate_IC_weights(IC)
  
  expect_equal(names(result$IC), c("H1", "H2"))
  expect_equal(names(result$IC_weights), c("H1", "H2"))
})

test_that("calculate_IC_weights met aangepaste namen", {
  IC <- c(50, 55, 60)
  result <- restriktor:::calculate_IC_weights(IC, hypo_names = c("M1", "M2", "M3"))
  
  expect_equal(names(result$IC), c("M1", "M2", "M3"))
  expect_equal(names(result$IC_weights), c("M1", "M2", "M3"))
})

test_that("calculate_IC_weights ratio matrix is correct", {
  IC <- c(100, 102)
  result <- restriktor:::calculate_IC_weights(IC)
  
  # ratio[i,j] = weight[i] / weight[j]
  w <- result$IC_weights
  expect_equal(result$ratio_IC_weights[1, 1], 1, tolerance = 1e-10)
  expect_equal(result$ratio_IC_weights[2, 2], 1, tolerance = 1e-10)
  expect_equal(as.numeric(result$ratio_IC_weights[1, 2]),
               as.numeric(w[1] / w[2]), tolerance = 1e-10)
})

test_that("calculate_IC_weights klasse is goric_ICw", {
  IC <- c(100, 105)
  result <- restriktor:::calculate_IC_weights(IC)
  
  expect_s3_class(result, "goric_ICw")
  expect_true(is.list(result))
  expect_true(all(c("IC", "IC_weights", "ratio_IC_weights") %in% names(result)))
})

test_that("calculate_IC_weights foutmelding bij ongeldige invoer", {
  # Matrix met meerdere kolommen
  expect_error(restriktor:::calculate_IC_weights(matrix(1:4, 2, 2)),
               "vector or a matrix with one column")
  
  # Verkeerde aantal namen
  expect_error(restriktor:::calculate_IC_weights(c(1, 2), hypo_names = c("A")),
               "NrHypos")
})

test_that("calculate_IC_weights met grote IC-verschillen", {
  # Zeer groot verschil: practisch alle gewicht op H1
  IC <- c(100, 200)
  result <- restriktor:::calculate_IC_weights(IC)
  
  expect_true(result$IC_weights[1] > 0.999)
  expect_true(result$IC_weights[2] < 0.001)
})

test_that("calc_ICweights alias werkt identiek", {
  IC <- c(100, 105, 110)
  r1 <- restriktor:::calculate_IC_weights(IC)
  r2 <- restriktor:::calc_ICweights(IC)
  
  expect_equal(r1$IC_weights, r2$IC_weights)
  expect_equal(r1$ratio_IC_weights, r2$ratio_IC_weights)
})
