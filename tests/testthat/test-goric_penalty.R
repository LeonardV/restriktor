# =============================================================================
# Tests: penalty_goric en penalty_complement_goric
# =============================================================================
# Zorg dat vereiste imports beschikbaar zijn in het huidige library-pad
# for (pkg in c("lavaan", "mvtnorm", "tmvtnorm", "quadprog", "norm")) {
#   if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# }

# # Laad het restriktor package vanuit broncode
# setwd("/Workspace/Users/l.vanbrabant@ggdwestbrabant.nl/restriktor")
# devtools::load_all(".")

# --- penalty_goric ---

test_that("penalty_goric met boot methode", {
  p <- 3
  Amat <- diag(p)
  LP <- c(0.1, 0.3, 0.4, 0.2)  # gewichten voor 0:p restricties
  attr(LP, "method") <- "boot"
  
  result <- restriktor:::penalty_goric(Amat, meq = 0L, LP, correction = FALSE)
  
  # PT = 1 + sum(0:p * LP) = 1 + (0*0.1 + 1*0.3 + 2*0.4 + 3*0.2) = 1 + 1.7 = 2.7
  expected <- 1 + sum(0:p * LP)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("penalty_goric met pmvnorm methode en restricties", {
  p <- 4
  q1 <- 2  # ongelijkheidsrestricties
  q2 <- 0  # gelijkheidsrestricties
  Amat <- matrix(0, nrow = q1, ncol = p)
  Amat[1, 1] <- 1
  Amat[2, 2] <- 1
  
  LP <- c(0.3, 0.5, 0.2)  # gewichten voor min_col:max_col
  attr(LP, "method") <- "pmvnorm"
  
  result <- restriktor:::penalty_goric(Amat, meq = q2, LP, correction = FALSE)
  
  min_col <- p - nrow(Amat)  # 4 - 2 = 2
  max_col <- p - q2           # 4 - 0 = 4
  expected <- 1 + sum(min_col:max_col * LP)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("penalty_goric met small sample correctie", {
  p <- 3
  Amat <- matrix(0, nrow = 1, ncol = p)
  LP <- rep(1 / (p + 1), p + 1)
  attr(LP, "method") <- "boot"
  N <- 100
  
  result <- restriktor:::penalty_goric(Amat, meq = 0L, LP, 
                                        correction = TRUE, 
                                        sample.nobs = N)
  
  # Met correctie: resultaat moet > standaard penalty zijn
  result_no_corr <- restriktor:::penalty_goric(Amat, meq = 0L, LP, 
                                                correction = FALSE)
  expect_true(result > result_no_corr)
})


# --- penalty_complement_goric ---

test_that("penalty_complement_goric berekent correcte complementpenalty (goric)", {
  p <- 3
  Amat <- diag(p)
  meq <- 0L
  
  wt.bar <- c(0.1, 0.3, 0.4, 0.2)
  attr(wt.bar, "method") <- "pmvnorm"
  
  result <- restriktor:::penalty_complement_goric(Amat, meq, type = "goric",
                                                   wt.bar = wt.bar)
  
  # PTc = 1 + p - wt.bar[last] * q1, met q1 = rank(Amat) - meq
  q1 <- qr(Amat)$rank - meq  # 3
  expected <- 1 + p - wt.bar[length(wt.bar)] * q1
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("penalty_complement_goric gorica trekt 1 af", {
  p <- 3
  Amat <- diag(p)
  meq <- 0L
  
  wt.bar <- c(0.1, 0.3, 0.4, 0.2)
  attr(wt.bar, "method") <- "pmvnorm"
  
  result_goric  <- restriktor:::penalty_complement_goric(Amat, meq, 
                                                          type = "goric", 
                                                          wt.bar = wt.bar)
  result_gorica <- restriktor:::penalty_complement_goric(Amat, meq, 
                                                          type = "gorica", 
                                                          wt.bar = wt.bar)
  
  # gorica = goric - 1
  expect_equal(result_gorica, result_goric - 1, tolerance = 1e-10)
})
