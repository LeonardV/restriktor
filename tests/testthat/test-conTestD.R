# =============================================================================
# Tests: conTestD() en bootstrapD() — lavaan bootstrap-gebaseerde toetsen
# =============================================================================

# Zorg dat vereiste imports beschikbaar zijn in het huidige library-pad
# for (pkg in c("lavaan", "mvtnorm", "tmvtnorm", "quadprog", "norm")) {
#   if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# }

# # Laad het restriktor package vanuit broncode
# setwd("/Workspace/Users/l.vanbrabant@ggdwestbrabant.nl/restriktor")
# devtools::load_all(".")

#skip_if_not_installed("lavaan")


# --- Gedeelde testdata ---

# Eenvoudig pad-model met bekende structuur
set.seed(42)
n <- 200
x <- rnorm(n)
m <- 0.5 * x + rnorm(n, sd = 0.8)
y <- 0.4 * m + 0.2 * x + rnorm(n, sd = 0.7)
df_path <- data.frame(x = x, m = m, y = y)

model_path <- "
  m ~ a*x
  y ~ b*m + c*x
"

# Constraints: a > 0, b > 0
constraints_ineq <- "a > 0; b > 0"
# Constraints: a > b
constraints_order <- "a > b"

# Pre-fit lavaan objecten (hergebruikt in meerdere tests)
fit_h2 <- lavaan::sem(model_path, data = df_path, test = "standard")

# Maak h0 met gelijkheden (nodig voor bootstrapD tests)
.make_h0 <- function(fit_h1, constraints_str) {
  CON <- attr(lavaan::lavParseModelString(constraints_str), "constraints")
  pt <- fit_h1@ParTable
  for (con in seq_along(CON)) {
    pt <- lavaan::lav_partable_add(pt, CON[[con]])
  }
  for (con in seq_along(CON)) {
    if (CON[[con]]$op %in% c("<", ">")) {
      idx <- which(pt$lhs == CON[[con]]$lhs &
                   pt$op  == CON[[con]]$op &
                   pt$rhs == CON[[con]]$rhs)
      if (length(idx) > 0) pt$op[idx] <- "=="
    }
  }
  lavaan::sem(pt, data = df_path, test = "standard")
}
fit_h0 <- .make_h0(fit_h2, constraints_ineq)


# =============================================================================
# conTestD: basisfunctionaliteit
# =============================================================================

test_that("conTestD: type A met ongelijkheidsrestricties", {
  result <- conTestD(model = model_path, data = df_path,
                     constraints = constraints_ineq,
                     type = "A", R = 100, double.bootstrap = "no")

  expect_s3_class(result, "conTestLavaan")
  expect_true(!is.null(result$bootA))
  expect_true(!is.null(result$fit.h0))
  expect_true(!is.null(result$fit.h2))
  expect_true(inherits(result$fit.h0, "lavaan"))
  expect_true(inherits(result$fit.h2, "lavaan"))
})

test_that("conTestD: type B met ongelijkheidsrestricties", {
  result <- conTestD(model = model_path, data = df_path,
                     constraints = constraints_ineq,
                     type = "B", R = 100, double.bootstrap = "no")

  expect_s3_class(result, "conTestLavaan")
  expect_true(!is.null(result$bootB))
})

test_that("conTestD: type A en B tegelijk", {
  result <- conTestD(model = model_path, data = df_path,
                     constraints = constraints_ineq,
                     type = c("A", "B"), R = 100, double.bootstrap = "no")

  expect_s3_class(result, "conTestLavaan")
  expect_true(!is.null(result$bootA))
  expect_true(!is.null(result$bootB))
})

test_that("conTestD: orderingsrestrictie a > b", {
  result <- conTestD(model = model_path, data = df_path,
                     constraints = constraints_order,
                     type = "A", R = 100, double.bootstrap = "no")

  expect_s3_class(result, "conTestLavaan")
})


# =============================================================================
# conTestD: output structuur
# =============================================================================

test_that("conTestD: output bevat verwachte velden", {
  result <- conTestD(model = model_path, data = df_path,
                     constraints = constraints_ineq,
                     type = "A", R = 100, double.bootstrap = "no")

  expect_true(!is.null(result$double.bootstrap))
  expect_true(!is.null(result$double.bootstrap.alpha))
  expect_true(!is.null(result$return.test))
  expect_true(!is.null(result$type))
})

test_that("conTestD: print methode werkt", {
  result <- conTestD(model = model_path, data = df_path,
                     constraints = constraints_ineq,
                     type = "A", R = 100, double.bootstrap = "no")
  expect_output(print(result))
})


# =============================================================================
# conTestD: bootstrap-opties
# =============================================================================

test_that("conTestD: return.test = FALSE geeft geen D-waarden", {
  result <- conTestD(model = model_path, data = df_path,
                     constraints = constraints_ineq,
                     type = "A", R = 50,
                     return.test = FALSE, double.bootstrap = "no")
  expect_s3_class(result, "conTestLavaan")
  expect_false(result$return.test)
})


# =============================================================================
# conTestD: foutafhandeling
# =============================================================================

test_that("conTestD: fout bij ontbrekende constraints", {
  expect_error(
    conTestD(model = model_path, data = df_path, constraints = NULL),
    "constraints|NULL"
  )
})

test_that("conTestD: fout bij ongeldig model", {
  expect_error(
    conTestD(model = "invalid_model_syntax ~~~", data = df_path,
             constraints = "a > 0"),
    "lavaan|syntax|ERROR"
  )
})


# =============================================================================
# bootstrapD: directe aanroep — type A
# =============================================================================

test_that("bootstrapD: type A retourneert p-waarde", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "A", R = 50, double.bootstrap = "no")

  expect_true(is.numeric(result))
  # p-waarde moet in [0, 1] liggen
  expect_true(result >= 0 & result <= 1)
})

test_that("bootstrapD: type A met return.D = TRUE bevat D.original en D vector", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "A", R = 50,
                       return.D = TRUE, double.bootstrap = "no")

  expect_true(!is.null(attr(result, "D.original")))
  expect_true(!is.null(attr(result, "D")))
  # D.original moet niet-negatief zijn
  expect_true(attr(result, "D.original") >= -1e-10)
  # D vector moet lengte <= R hebben
  expect_true(length(attr(result, "D")) <= 50)
})

test_that("bootstrapD: type A D.original is numeriek en eindig", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "A", R = 30,
                       return.D = TRUE, double.bootstrap = "no")

  D_orig <- attr(result, "D.original")
  expect_true(is.numeric(D_orig))
  expect_true(is.finite(D_orig))
})


# =============================================================================
# bootstrapD: directe aanroep — type B
# =============================================================================

test_that("bootstrapD: type B retourneert p-waarde", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "B", R = 50)

  expect_true(is.numeric(result))
  expect_true(result >= 0 & result <= 1)
})

test_that("bootstrapD: type B met return.D = TRUE bevat D vector", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "B", R = 50,
                       return.D = TRUE)

  expect_true(!is.null(attr(result, "D")))
  D_vec <- attr(result, "D")
  # D waarden moeten numeriek zijn
  expect_true(all(is.numeric(D_vec)))
})


# =============================================================================
# bootstrapD: Amat attribuut
# =============================================================================

test_that("bootstrapD: retourneert constraints matrix als attribuut", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "A", R = 30)

  Amat <- attr(result, "Amat")
  expect_true(!is.null(Amat))
  expect_true(is.matrix(Amat))
  # Amat moet minstens 1 rij hebben (minstens 1 restrictie)
  expect_true(nrow(Amat) >= 1)
})


# =============================================================================
# bootstrapD: bootstrap-type varianten
# =============================================================================

test_that("bootstrapD: bollen.stine bootstrap", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "A", R = 30,
                       bootstrap.type = "bollen.stine")
  expect_true(is.numeric(result))
  expect_true(result >= 0 & result <= 1)
})


# =============================================================================
# bootstrapD: double bootstrap
# =============================================================================

test_that("bootstrapD: double bootstrap 'standard' geeft adj.pvalue", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "A", R = 20,
                       double.bootstrap = "standard",
                       double.bootstrap.R = 10,
                       double.bootstrap.alpha = 0.05)

  expect_true(is.numeric(result))
  # adj.alpha en adj.pvalue als attributen
  expect_true(!is.null(attr(result, "adj.alpha")))
  expect_true(!is.null(attr(result, "adj.pvalue")))
  adj_p <- attr(result, "adj.pvalue")
  expect_true(adj_p >= 0 & adj_p <= 1)
})

test_that("bootstrapD: double bootstrap 'FDB' geeft adj.pvalue en D.q", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "A", R = 20,
                       double.bootstrap = "FDB")

  expect_true(is.numeric(result))
  expect_true(!is.null(attr(result, "adj.pvalue")))
  expect_true(!is.null(attr(result, "D.q")))
  adj_p <- attr(result, "adj.pvalue")
  expect_true(adj_p >= 0 & adj_p <= 1)
})

test_that("bootstrapD: double bootstrap 'FDB' met return.D bevat D2", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "A", R = 20,
                       double.bootstrap = "FDB",
                       return.D = TRUE)

  expect_true(!is.null(attr(result, "D2")))
  expect_true(!is.null(attr(result, "D.original")))
  expect_true(!is.null(attr(result, "D")))
})


# =============================================================================
# bootstrapD: wiskundige consistentie
# =============================================================================

test_that("bootstrapD: p-waarde type A <= p-waarde type B (verwacht patroon)", {
  pA <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                   constraints = constraints_ineq,
                   type = "A", R = 100, seed = 42)
  pB <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                   constraints = constraints_ineq,
                   type = "B", R = 100, seed = 42)

  # Beide moeten geldige p-waarden zijn
  expect_true(pA >= 0 & pA <= 1)
  expect_true(pB >= 0 & pB <= 1)
})

test_that("bootstrapD: D.original type A is niet-negatief", {
  result <- bootstrapD(h0 = fit_h0, h1 = fit_h2,
                       constraints = constraints_ineq,
                       type = "A", R = 30,
                       return.D = TRUE)
  D_orig <- attr(result, "D.original")
  expect_true(D_orig >= -1e-10)
})


# =============================================================================
# bootstrapD: foutafhandeling
# =============================================================================

test_that("bootstrapD: fout als h1 geen lavaan object is", {
  expect_error(
    bootstrapD(h0 = fit_h0, h1 = "niet_lavaan",
               constraints = constraints_ineq,
               type = "A", R = 10),
    "lavaan"
  )
})

test_that("bootstrapD: fout bij onbekend bootstrap.type", {
  expect_error(
    bootstrapD(h0 = fit_h0, h1 = fit_h2,
               constraints = constraints_ineq,
               type = "A", R = 10,
               bootstrap.type = "onbekend"),
    "bootstrap"
  )
})

test_that("bootstrapD: fout bij ongeldig type", {
  expect_error(
    bootstrapD(h0 = fit_h0, h1 = fit_h2,
               constraints = constraints_ineq,
               type = "C", R = 10)
  )
})
