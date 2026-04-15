# =============================================================================
# Tests: con_solver functies (constraint solvers via quadprog)
# =============================================================================
# Zorg dat vereiste imports beschikbaar zijn in het huidige library-pad
# for (pkg in c("lavaan", "mvtnorm", "tmvtnorm", "quadprog", "norm")) {
#   if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# }

# # Laad het restriktor package vanuit broncode
# setwd("/Workspace/Users/l.vanbrabant@ggdwestbrabant.nl/restriktor")
# devtools::load_all(".")

# --- con_solver_gorica ---

test_that("con_solver_gorica vindt onbeperkt optimum als restricties inactief", {
  # Schattingen al in de feasible regio: restricties zijn niet-bindend
  est <- c(3, 2, 1)
  VCOV <- diag(3)
  # Restrictie: x1 >= 0, x2 >= 0, x3 >= 0 (reeds voldaan)
  Amat <- diag(3)
  bvec <- rep(0, 3)
  meq <- 0L
  
  result <- restriktor:::con_solver_gorica(est, VCOV, Amat, bvec, meq)
  
  # Oplossing moet (bijna) gelijk zijn aan originele schattingen
  expect_equal(result$solution, est, tolerance = 1e-6)
})

test_that("con_solver_gorica respecteert ongelijkheidsrestricties", {
  # Schattingen: (−1, 2) met VCOV = I
  # Restrictie: x1 >= 0
  est <- c(-1, 2)
  VCOV <- diag(2)
  Amat <- matrix(c(1, 0), nrow = 1)
  bvec <- 0
  meq <- 0L
  
  result <- restriktor:::con_solver_gorica(est, VCOV, Amat, bvec, meq)
  
  # x1 moet >= 0 zijn (restrictie actief)
  expect_true(result$solution[1] >= -1e-10)
  # x2 moet vrij zijn (onveranderd ≈ 2)
  expect_equal(result$solution[2], 2, tolerance = 1e-6)
})

test_that("con_solver_gorica respecteert gelijkheidsbeperkingen", {
  # Schattingen: (3, 1) met VCOV = I
  # Gelijkheidsbeperking: x1 = x2
  est <- c(3, 1)
  VCOV <- diag(2)
  Amat <- matrix(c(1, -1), nrow = 1)
  bvec <- 0
  meq <- 1L
  
  result <- restriktor:::con_solver_gorica(est, VCOV, Amat, bvec, meq)
  
  # x1 moet gelijk zijn aan x2
  expect_equal(result$solution[1], result$solution[2], tolerance = 1e-6)
  # Optimale waarde is het gemiddelde: (3+1)/2 = 2
  expect_equal(result$solution[1], 2, tolerance = 1e-6)
})

test_that("con_solver_gorica met ordening x1 >= x2 >= x3", {
  est <- c(1, 3, 2)  # volgorde geschonden
  VCOV <- diag(3)
  # x1 - x2 >= 0 en x2 - x3 >= 0
  Amat <- matrix(c(1, -1, 0,
                    0,  1, -1), nrow = 2, byrow = TRUE)
  bvec <- c(0, 0)
  meq <- 0L
  
  result <- restriktor:::con_solver_gorica(est, VCOV, Amat, bvec, meq)
  
  # Oplossing moet geordend zijn: x1 >= x2 >= x3
  expect_true(result$solution[1] >= result$solution[2] - 1e-10)
  expect_true(result$solution[2] >= result$solution[3] - 1e-10)
})


# --- con_solver_lm ---

test_that("con_solver_lm vindt beperkte oplossing bij eenvoudige regressie", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n), rnorm(n))
  beta_true <- c(1, 0.5, -0.3)
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  
  # Restrictie: beta2 >= 0 (al voldaan in beta_true)
  Amat <- matrix(c(0, 1, 0), nrow = 1)
  bvec <- 0
  meq <- 0L
  
  result <- restriktor:::con_solver_lm(X, y, w = NULL, Amat, bvec, meq)
  
  # Beperkte schatting beta2 >= 0
  expect_true(result$qp$solution[2] >= -1e-10)
  
  # s2 moet positief zijn
  expect_true(all(result$s2 > 0))
})

test_that("con_solver_lm convergeert bij gelijkheidsbeperkingen", {
  set.seed(123)
  n <- 100
  X <- cbind(1, rnorm(n), rnorm(n))
  y <- X %*% c(2, 1, 1) + rnorm(n)
  
  # Gelijkheid: beta2 = beta3
  Amat <- matrix(c(0, 1, -1), nrow = 1)
  bvec <- 0
  meq <- 1L
  
  result <- restriktor:::con_solver_lm(X, y, w = NULL, Amat, bvec, meq)
  
  # beta2 moet (bijna) gelijk zijn aan beta3
  expect_equal(result$qp$solution[2], result$qp$solution[3], tolerance = 1e-6)
})

test_that("con_solver_lm geeft zelfde resultaat als OLS zonder restricties", {
  set.seed(99)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- X %*% c(1, 2) + rnorm(n, sd = 0.5)
  
  # Triviale restrictie: altijd voldaan (geen echte beperking)
  Amat <- matrix(c(0, 1), nrow = 1)
  bvec <- -100  # beta1 >= -100 (altijd waar)
  meq <- 0L
  
  result <- restriktor:::con_solver_lm(X, y, w = NULL, Amat, bvec, meq)
  ols <- lm.fit(X, y)
  
  expect_equal(result$qp$solution, unname(coef(ols)), tolerance = 1e-4)
})
