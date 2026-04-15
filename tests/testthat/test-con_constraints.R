# =============================================================================
# Tests: con_constraints (constraint parsing en matrix-opbouw)
# =============================================================================


# =============================================================================
# 1. Character constraints met numerieke vector
# =============================================================================

test_that("con_constraints parseert enkelvoudige ongelijkheidsrestrictie", {
  est <- c(x1 = 3, x2 = 1, x3 = 2)
  VCOV <- diag(3)
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = "x1 > x2"
  )
  
  # Output structuur controleren
  expect_type(result, "list")
  expect_named(result, c("CON", "parTable", "Amat", "bvec", "meq"))
  
  # Geen gelijkheidsbeperkingen
  expect_equal(result$meq, 0L)
  
  # Amat moet 1 rij hebben (1 restrictie)
  expect_equal(nrow(result$Amat), 1)
  expect_equal(ncol(result$Amat), 3)
  
  # bvec moet lengte 1 zijn met waarde 0
  expect_equal(length(result$bvec), 1)
  expect_equal(result$bvec, 0)
  
  # Amat rij: x1 - x2 >= 0 → coëfficiënten (1, -1, 0)
  expect_equal(as.numeric(result$Amat[1, ]), c(1, -1, 0))
  
  # CON object moet bestaan
  expect_false(is.null(result$CON))
})


test_that("con_constraints parseert gelijkheidsbeperkingen correct", {
  est <- c(x1 = 3, x2 = 1)
  VCOV <- diag(2)
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = "x1 == x2"
  )
  
  # Eén gelijkheidsbeperking
  expect_equal(result$meq, 1L)
  
  # Amat: x1 - x2 = 0 → coëfficiënten (1, -1)
  expect_equal(nrow(result$Amat), 1)
  expect_equal(as.numeric(result$Amat[1, ]), c(1, -1))
  
  # bvec = 0
  expect_equal(result$bvec, 0)
})


test_that("con_constraints parseert meerdere ongelijkheden", {
  est <- c(x1 = 3, x2 = 2, x3 = 1)
  VCOV <- diag(3)
  
  # Ordening: x1 > x2 > x3
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = "x1 > x2; x2 > x3"
  )
  
  expect_equal(result$meq, 0L)
  expect_equal(nrow(result$Amat), 2)
  
  # Rij 1: x1 - x2 >= 0 → (1, -1, 0)
  # Rij 2: x2 - x3 >= 0 → (0, 1, -1)
  expect_equal(as.numeric(result$Amat[1, ]), c(1, -1, 0))
  expect_equal(as.numeric(result$Amat[2, ]), c(0, 1, -1))
  
  expect_equal(result$bvec, c(0, 0))
})


test_that("con_constraints parseert gemengde gelijkheid en ongelijkheid", {
  est <- c(x1 = 3, x2 = 2, x3 = 1)
  VCOV <- diag(3)
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = "x1 == x2; x2 > x3"
  )
  
  # 1 gelijkheid + 1 ongelijkheid
  expect_equal(result$meq, 1L)
  expect_equal(nrow(result$Amat), 2)
})


test_that("con_constraints parseert constraint met constante rhs", {
  est <- c(x1 = 5, x2 = 2)
  VCOV <- diag(2)
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = "x1 > 2"
  )
  
  expect_equal(result$meq, 0L)
  expect_equal(nrow(result$Amat), 1)
  
  # bvec moet 2 zijn (x1 >= 2)
  expect_equal(result$bvec, 2)
  
  # Amat: alleen x1 coëfficiënt = 1
  expect_equal(as.numeric(result$Amat[1, ]), c(1, 0))
})


test_that("con_constraints bewaart oorspronkelijke constraints in CON", {
  est <- c(x1 = 3, x2 = 1)
  VCOV <- diag(2)
  constraints_text <- "x1 > x2"
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = constraints_text
  )
  
  expect_equal(result$CON$constraints, constraints_text)
})


test_that("con_constraints parseert compound constraints (x1 > x2 > x3)", {
  est <- c(x1 = 3, x2 = 2, x3 = 1)
  VCOV <- diag(3)
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = "x1 > x2 > x3"
  )
  
  expect_equal(result$meq, 0L)
  # Compound wordt gesplitst in 2 restricties
  expect_equal(nrow(result$Amat), 2)
  expect_equal(as.numeric(result$Amat[1, ]), c(1, -1, 0))
  expect_equal(as.numeric(result$Amat[2, ]), c(0, 1, -1))
})


# =============================================================================
# 2. Matrix constraints
# =============================================================================

test_that("con_constraints accepteert constraint matrix", {
  est <- c(x1 = 3, x2 = 1, x3 = 2)
  VCOV <- diag(3)
  
  # x1 >= 0, x2 >= 0, x3 >= 0
  Amat <- diag(3)
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = Amat
  )
  
  # CON moet NULL zijn bij matrix-invoer
  expect_null(result$CON)
  
  # Amat, bvec, meq moeten correct zijn
  expect_equal(nrow(result$Amat), 3)
  expect_equal(result$bvec, rep(0L, 3))
  expect_equal(result$meq, 0L)
})


test_that("con_constraints respecteert custom bvec en meq bij matrix", {
  est <- c(x1 = 3, x2 = 1)
  VCOV <- diag(2)
  
  Amat <- matrix(c(1, -1), nrow = 1)
  bvec_custom <- 2  # x1 - x2 >= 2
  meq_custom  <- 0L
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = Amat, bvec = bvec_custom, meq = meq_custom
  )
  
  expect_equal(result$bvec, 2)
  expect_equal(result$meq, 0L)
})


test_that("con_constraints converteert vector naar matrix", {
  est <- c(x1 = 3, x2 = 1)
  VCOV <- diag(2)
  
  # Vector in plaats van matrix
  constraints_vec <- c(1, 0)
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = constraints_vec
  )
  
  # Moet geconverteerd zijn naar matrix met 1 rij
  expect_true(is.matrix(result$Amat))
  expect_equal(nrow(result$Amat), 1)
})


# =============================================================================
# 3. lm model input
# =============================================================================

test_that("con_constraints werkt met lm model en character constraints", {
  set.seed(42)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 2 + 0.5 * x1 + 0.3 * x2 + rnorm(n, sd = 0.5)
  model <- lm(y ~ x1 + x2)
  
  VCOV <- vcov(model)
  est  <- coef(model)
  
  result <- restriktor:::con_constraints(
    model = model, VCOV = VCOV, est = est,
    constraints = "x1 > x2"
  )
  
  expect_type(result, "list")
  expect_equal(result$meq, 0L)
  expect_equal(nrow(result$Amat), 1)
  # Intercept niet betrokken: Amat rij moet 0 hebben op positie 1
  expect_equal(as.numeric(result$Amat[1, 1]), 0)
})


# =============================================================================
# 4. Foutafhandeling
# =============================================================================


test_that("con_constraints stopt bij lavaan-operatoren in constraints", {
  est <- c(x1 = 1, x2 = 2)
  VCOV <- diag(2)
  
  expect_error(
    restriktor:::con_constraints(
      model = est, VCOV = VCOV, est = est,
      constraints = "x1 =~ x2"
    ),
    "error in constraint syntax"
  )
  
  expect_error(
    restriktor:::con_constraints(
      model = est, VCOV = VCOV, est = est,
      constraints = "x1 ~~ x2"
    ),
    "error in constraint syntax"
  )
})


test_that("con_constraints stopt bij dubbele operatoren (>> of <<)", {
  est <- c(x1 = 1, x2 = 2)
  VCOV <- diag(2)
  
  expect_error(
    restriktor:::con_constraints(
      model = est, VCOV = VCOV, est = est,
      constraints = "x1 >> x2"
    ),
    "error in constraint syntax"
  )
})


test_that("con_constraints stopt wanneer meq > nrow(Amat)", {
  est <- c(x1 = 1, x2 = 2)
  VCOV <- diag(2)
  Amat <- matrix(c(1, -1), nrow = 1)
  
  expect_error(
    restriktor:::con_constraints(
      model = est, VCOV = VCOV, est = est,
      constraints = Amat, meq = 5L
    ),
    "meq.*cannot exceed.*nrow"
  )
})


# =============================================================================
# 5. Redundante constraints en randgevallen
# =============================================================================

test_that("con_constraints verwijdert redundante constraints", {
  est <- c(x1 = 5, x2 = 2)
  VCOV <- diag(2)
  
  # Twee identieke restricties: x1 >= 1 en x1 >= 2
  # Na remove_redundant_constraints blijft alleen x1 >= 1 over (strengste)
  Amat <- matrix(c(1, 0,
                   1, 0), nrow = 2, byrow = TRUE)
  bvec <- c(1, 2)
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = Amat, bvec = bvec
  )
  
  # Na deduplicatie: 1 rij over
  expect_equal(nrow(result$Amat), 1)
})


test_that("con_constraints output dimensies zijn consistent", {
  est <- c(x1 = 3, x2 = 2, x3 = 1)
  VCOV <- diag(3)
  
  result <- restriktor:::con_constraints(
    model = est, VCOV = VCOV, est = est,
    constraints = "x1 > x2; x2 > x3"
  )
  
  # nrow(Amat) moet gelijk zijn aan length(bvec)
  expect_equal(nrow(result$Amat), length(result$bvec))
  
  # ncol(Amat) moet gelijk zijn aan het aantal parameters
  expect_equal(ncol(result$Amat), length(est))
  
  # meq <= nrow(Amat)
  expect_true(result$meq <= nrow(result$Amat))
})
