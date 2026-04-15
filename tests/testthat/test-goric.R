# =============================================================================
# Tests: goric() end-to-end
# =============================================================================

# Zorg dat vereiste imports beschikbaar zijn in het huidige library-pad
# for (pkg in c("lavaan", "mvtnorm", "tmvtnorm", "quadprog", "norm")) {
#   if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
# }
# 
# # Laad het restriktor package vanuit broncode
# setwd("/Workspace/Users/l.vanbrabant@ggdwestbrabant.nl/restriktor")
# devtools::load_all(".")

# --- Gedeelde testdata ---

# Groepsgemiddelden: mu1 > mu2 > mu3 (ware ordening)
set.seed(123)
n_per_group <- 30
group <- factor(rep(1:3, each = n_per_group))
y_ordered <- c(rnorm(n_per_group, mean = 5),
               rnorm(n_per_group, mean = 3),
               rnorm(n_per_group, mean = 1))
df_anova <- data.frame(y = y_ordered, group = group)
fit_lm <- lm(y ~ -1 + group, data = df_anova)

# Schattingen en VCOV voor gorica-pad
est_3 <- c(x1 = 5.0, x2 = 3.0, x3 = 1.0)
VCOV_3 <- diag(c(0.04, 0.04, 0.04))


# =============================================================================
# gorica via numeric (schattingen + VCOV)
# =============================================================================

test_that("goric.numeric: gorica met 1 hypothese en complement", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  # Output klasse
  expect_s3_class(result, "con_goric")

  # Result data.frame structuur
  expect_true(is.data.frame(result$result))
  expect_true("gorica" %in% names(result$result))
  expect_true("loglik" %in% names(result$result))
  expect_true("penalty" %in% names(result$result))

  # Bij 1 hypothese: default comparison = complement => 2 rijen
  expect_equal(nrow(result$result), 2)
  expect_true("complement" %in% result$result$model)

  # Gewichten tellen op tot 1
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(sum(result$result[[wt_col]]), 1, tolerance = 1e-6)
})

test_that("goric.numeric: gorica met meerdere hypothesen en unconstrained", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3",
                                    H2 = "x1 > x2 = x3"))

  # Bij >1 hypothese: default comparison = unconstrained => 3 rijen
  expect_equal(nrow(result$result), 3)
  expect_true("unconstrained" %in% result$result$model)

  # Gewichten optellen tot 1
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(sum(result$result[[wt_col]]), 1, tolerance = 1e-6)
})

test_that("goric.numeric: gorica comparison = 'none'", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  comparison = "none",
                  hypotheses = list(H1 = "x1 > x2 > x3",
                                    H2 = "x1 > x2 = x3"))

  # Geen unconstrained of complement rij
  expect_equal(nrow(result$result), 2)
  expect_false("unconstrained" %in% result$result$model)
  expect_false("complement" %in% result$result$model)
})

test_that("goric.numeric: ware hypothese krijgt hoogste gewicht", {
  # Data is gegenereerd met x1 > x2 > x3
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H_correct = "x1 > x2 > x3",
                                    H_fout    = "x3 > x2 > x1"))

  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  weights <- result$result[[wt_col]]
  names(weights) <- result$result$model

  # H_correct moet meer gewicht krijgen dan H_fout
  expect_true(weights["H_correct"] > weights["H_fout"])
})


# =============================================================================
# goric via lm object
# =============================================================================

test_that("goric.lm: goric met 1 hypothese en complement", {
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = "group1 > group2 > group3"))

  expect_s3_class(result, "con_goric")
  expect_true(is.data.frame(result$result))
  expect_true("goric" %in% names(result$result))

  # 1 hypothese => complement vergelijking
  expect_equal(nrow(result$result), 2)
  expect_true("complement" %in% result$result$model)
})

test_that("goric.lm: goric met meerdere hypothesen", {
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = "group1 > group2 > group3",
                                    H2 = "group1 > group2 = group3"))

  # 2 hypothesen => unconstrained vergelijking => 3 rijen
  expect_equal(nrow(result$result), 3)
  expect_true("unconstrained" %in% result$result$model)
})

test_that("goric.lm: ware ordening krijgt hoogste gewicht", {
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H_correct = "group1 > group2 > group3",
                                    H_reverse = "group3 > group2 > group1"))

  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  weights <- result$result[[wt_col]]
  names(weights) <- result$result$model

  expect_true(weights["H_correct"] > weights["H_reverse"])
})


# =============================================================================
# Output structuur
# =============================================================================

test_that("goric output bevat alle verwachte velden", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  # Kernvelden
  expect_true(!is.null(result$result))
  expect_true(!is.null(result$ratio.gw))
  expect_true(!is.null(result$ratio.pw))
  expect_true(!is.null(result$ratio.lw))
  expect_true(!is.null(result$objectList))
  expect_true(!is.null(result$objectNames))
  expect_true(!is.null(result$comparison))
  expect_true(!is.null(result$type))
  expect_true(!is.null(result$VCOV))
  expect_true(!is.null(result$b.unrestr))
  expect_true(!is.null(result$ormle))
  expect_true(!is.null(result$constraints))
  expect_true(!is.null(result$rhs))
  expect_true(!is.null(result$neq))
})

test_that("goric ratio.gw matrix is vierkant en benoemd", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3",
                                    H2 = "x1 > x2 = x3"))

  rw <- result$ratio.gw
  n_models <- nrow(result$result)
  expect_equal(nrow(rw), n_models)
  expect_equal(ncol(rw), n_models)
  # Diagonaal = 1
  expect_equal(diag(rw), rep(1, n_models), tolerance = 1e-10)
})

test_that("goric ormle bevat beperkte schattingen", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  coefs <- coef(result)
  expect_true(is.matrix(coefs) || is.data.frame(coefs))
  # Rijen = hypothesen + vergelijkingsmodel
  expect_equal(nrow(coefs), nrow(result$result))
})


# =============================================================================
# type varianten
# =============================================================================

test_that("goric type = 'gorica' geeft con_gorica klasse", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_true("con_gorica" %in% class(result))
  expect_equal(result$type, "gorica")
})

test_that("goric type = 'goricac' vereist sample_nobs", {
  expect_error(
    goric(est_3, VCOV = VCOV_3, type = "goricac",
          hypotheses = list(H1 = "x1 > x2 > x3")),
    "sample_nobs"
  )
})

test_that("goric type = 'goricac' met sample_nobs werkt", {
  result <- goric(est_3, VCOV = VCOV_3, type = "goricac",
                  sample_nobs = 90,
                  hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_s3_class(result, "con_goric")
  expect_equal(result$type, "goricac")
  expect_equal(result$sample_nobs, 90)
})


# =============================================================================
# Gelijkheidsbeperkingen
# =============================================================================

test_that("goric met alleen gelijkheidsbeperkingen", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 = x2 = x3"))

  expect_s3_class(result, "con_goric")
  # Alleen gelijkheden + complement => unconstrained ipv complement
  # (complement van pure gelijkheid = unconstrained)
  expect_true("unconstrained" %in% result$result$model)
})


# =============================================================================
# Heq parameter
# =============================================================================

test_that("goric met Heq = TRUE voegt gelijkheidsversie toe", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  comparison = "complement", Heq = TRUE,
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  # Heq rij moet aanwezig zijn
  expect_true("Heq" %in% result$result$model)
  expect_true(result$Heq)
})


# =============================================================================
# Foutafhandeling
# =============================================================================

test_that("goric fout bij ontbrekende hypotheses", {
  expect_error(
    goric(est_3, VCOV = VCOV_3, type = "gorica"),
    "hypotheses"
  )
})

test_that("goric fout bij ongeldige hypotheses (geen lijst)", {
  expect_error(
    goric(est_3, VCOV = VCOV_3, type = "gorica",
          hypotheses = "x1 > x2"),
    "hypotheses"
  )
})

test_that("goric fout bij niet-positief definiete VCOV", {
  bad_vcov <- matrix(c(1, 2, 2, 1), 2, 2)  # negatieve eigenwaarde
  est_2 <- c(x1 = 1, x2 = 2)
  expect_error(
    goric(est_2, VCOV = bad_vcov, type = "gorica",
          hypotheses = list(H1 = "x1 > x2")),
    "positive"
  )
})

test_that("goric fout bij ongeldig penalty_factor", {
  expect_error(
    goric(est_3, VCOV = VCOV_3, type = "gorica",
          penalty_factor = -1,
          hypotheses = list(H1 = "x1 > x2 > x3")),
    "penalty factor"
  )
})


# =============================================================================
# Wiskundige consistentie
# =============================================================================

test_that("goric: gorica = -2*loglik + 2*penalty", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  comparison = "none",
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  df_row <- result$result[1, ]
  expected_gorica <- -2 * df_row$loglik + 2 * df_row$penalty
  expect_equal(df_row$gorica, expected_gorica, tolerance = 1e-8)
})

test_that("goric: loglik onbeperkt model >= loglik beperkt model", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3",
                                    H2 = "x1 > x2 = x3"))

  df <- result$result
  ll_unc <- df$loglik[df$model == "unconstrained"]
  ll_models <- df$loglik[df$model != "unconstrained"]

  # Onbeperkte loglik moet >= alle beperkte logliks zijn
  for (ll in ll_models) {
    expect_true(ll_unc >= ll - 1e-6)
  }
})

test_that("goric: penalty onbeperkt > penalty beperkt (bij ongelijkheidsrestricties)", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  df <- result$result
  pt_unc <- df$penalty[df$model == "unconstrained" | df$model == "complement"]
  pt_h1 <- df$penalty[df$model == "H1"]

  # Beperkt model moet lagere penalty hebben dan onbeperkt
  expect_true(all(pt_h1 <= max(pt_unc) + 1e-6))
})

test_that("goric: gewichten zijn consistent met gorica waarden", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3",
                                    H2 = "x3 > x2 > x1"))

  df <- result$result
  wt_col <- grep("weights$", names(df), value = TRUE)[1]

  # Model met laagste gorica moet hoogste gewicht hebben
  min_gorica_idx <- which.min(df$gorica)
  max_weight_idx <- which.max(df[[wt_col]])
  expect_equal(min_gorica_idx, max_weight_idx)
})


# =============================================================================
# Matrixrestricties: goric.lm met list(constraints, rhs, neq)
# =============================================================================

# --- Basisfunctionaliteit ---

test_that("goric.lm matrix: enkele ongelijkheidsrestrictie (vector)", {
  # group1 > group2  =>  c(1, -1, 0) %*% beta >= 0
  h1 <- list(constraints = c(1, -1, 0))
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = h1))

  expect_s3_class(result, "con_goric")
  expect_equal(nrow(result$result), 2)  # H1 + complement
  expect_true("goric" %in% names(result$result))
})

test_that("goric.lm matrix: meerdere ongelijkheidsrestricties (matrix)", {
  # group1 > group2 > group3
  h1 <- list(constraints = rbind(c(1, -1,  0),
                                  c(0,  1, -1)),
             rhs = c(0, 0),
             neq = 0)
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = h1))

  expect_s3_class(result, "con_goric")
  expect_equal(nrow(result$result), 2)

  # Beperkte schattingen moeten geordend zijn
  b_restr <- coef(result)[1, ]
  expect_true(b_restr[1] >= b_restr[2] - 1e-6)
  expect_true(b_restr[2] >= b_restr[3] - 1e-6)
})

test_that("goric.lm matrix: gelijkheidsbeperkingen met neq", {
  # group1 = group2 (gelijkheid), group2 > group3 (ongelijkheid)
  h1 <- list(constraints = rbind(c(1, -1,  0),   # group1 = group2 (eerste neq rijen)
                                  c(0,  1, -1)),   # group2 > group3
             rhs = c(0, 0),
             neq = 1)
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = h1))

  expect_s3_class(result, "con_goric")
  b_restr <- coef(result)[1, ]
  # group1 moet (bijna) gelijk zijn aan group2
  expect_equal(unname(b_restr[1]), unname(b_restr[2]), tolerance = 1e-4)
  # group2 >= group3
  expect_true(b_restr[2] >= b_restr[3] - 1e-6)
})

test_that("goric.lm matrix: alleen gelijkheidsbeperkingen", {
  # group1 = group2 = group3
  h1 <- list(constraints = rbind(c(1, -1,  0),
                                  c(0,  1, -1)),
             rhs = c(0, 0),
             neq = 2)
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = h1))

  expect_s3_class(result, "con_goric")
  b_restr <- coef(result)[1, ]
  # Alle groepen gelijk
  expect_equal(unname(b_restr[1]), unname(b_restr[2]), tolerance = 1e-4)
  expect_equal(unname(b_restr[2]), unname(b_restr[3]), tolerance = 1e-4)
})


# --- Niet-nul rhs ---

test_that("goric.lm matrix: rhs ongelijk aan nul", {
  # group1 - group2 >= 0.5  (verschil minstens 0.5)
  h1 <- list(constraints = c(1, -1, 0),
             rhs = 0.5,
             neq = 0)
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H1 = h1))

  expect_s3_class(result, "con_goric")
  b_restr <- coef(result)[1, ]
  # group1 - group2 >= 0.5
  expect_true(b_restr[1] - b_restr[2] >= 0.5 - 1e-4)
})


# --- Meerdere hypothesen met matrixvorm ---

test_that("goric.lm matrix: meerdere hypothesen tegelijk", {
  h1 <- list(constraints = rbind(c(1, -1,  0),
                                  c(0,  1, -1)),
             rhs = c(0, 0),
             neq = 0)
  h2 <- list(constraints = rbind(c(-1, 1,  0),
                                  c( 0, -1, 1)),
             rhs = c(0, 0),
             neq = 0)
  result <- goric(fit_lm, type = "goric",
                  hypotheses = list(H_correct = h1, H_reverse = h2))

  expect_equal(nrow(result$result), 3)  # H_correct + H_reverse + unconstrained
  expect_true("unconstrained" %in% result$result$model)

  # Correcte hypothese (group1 > group2 > group3) moet meer gewicht krijgen
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  weights <- result$result[[wt_col]]
  names(weights) <- result$result$model
  expect_true(weights["H_correct"] > weights["H_reverse"])
})


# --- Equivalentie met string-hypothesen ---

test_that("goric.lm matrix: equivalent aan string-hypothese (complement)", {
  # Matrix: group1 > group2 > group3
  h_mat <- list(constraints = rbind(c(1, -1,  0),
                                     c(0,  1, -1)),
                rhs = c(0, 0),
                neq = 0)
  res_mat <- goric(fit_lm, type = "goric",
                   hypotheses = list(H1 = h_mat))

  # String: identiek
  res_str <- goric(fit_lm, type = "goric",
                   hypotheses = list(H1 = "group1 > group2 > group3"))

  # Loglik en penalty moeten gelijk zijn
  expect_equal(res_mat$result$loglik[1],  res_str$result$loglik[1],  tolerance = 1e-4)
  expect_equal(res_mat$result$penalty[1], res_str$result$penalty[1], tolerance = 1e-4)
  expect_equal(res_mat$result$goric[1],   res_str$result$goric[1],   tolerance = 1e-4)
})

test_that("goric.lm matrix: equivalent aan string-hypothese (unconstrained)", {
  # Matrix: twee hypothesen
  h1_mat <- list(constraints = rbind(c(1, -1,  0),
                                      c(0,  1, -1)),
                 rhs = c(0, 0), neq = 0)
  h2_mat <- list(constraints = c(1, -1, 0),
                 rhs = 0, neq = 1)  # group1 = group2
  res_mat <- goric(fit_lm, type = "goric",
                   hypotheses = list(H1 = h1_mat, H2 = h2_mat))

  # String: identiek
  res_str <- goric(fit_lm, type = "goric",
                   hypotheses = list(H1 = "group1 > group2 > group3",
                                     H2 = "group1 = group2"))

  # Gewichten moeten vergelijkbaar zijn
  wt_col_mat <- grep("weights", names(res_mat$result), value = TRUE)[1]
  wt_col_str <- grep("weights", names(res_str$result), value = TRUE)[1]
  expect_equal(res_mat$result[[wt_col_mat]], res_str$result[[wt_col_str]],
               tolerance = 1e-3)
})


# --- Matrixrestricties via gorica (numeric) ---

test_that("goric.numeric matrix: gorica met matrixrestricties", {
  h1 <- list(constraints = rbind(c(1, -1,  0),
                                  c(0,  1, -1)),
             rhs = c(0, 0),
             neq = 0)
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = h1))

  expect_s3_class(result, "con_goric")
  expect_true("con_gorica" %in% class(result))
  expect_equal(nrow(result$result), 2)  # H1 + complement
})

test_that("goric.numeric matrix: gorica equivalent aan string", {
  h_mat <- list(constraints = rbind(c(1, -1,  0),
                                     c(0,  1, -1)),
                rhs = c(0, 0), neq = 0)
  res_mat <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                   hypotheses = list(H1 = h_mat))
  res_str <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                   hypotheses = list(H1 = "x1 > x2 > x3"))

  expect_equal(res_mat$result$loglik[1],  res_str$result$loglik[1],  tolerance = 1e-6)
  expect_equal(res_mat$result$penalty[1], res_str$result$penalty[1], tolerance = 1e-6)
})

test_that("goric.numeric matrix: gorica met meerdere hypothesen en neq", {
  h1 <- list(constraints = rbind(c(1, -1,  0),
                                  c(0,  1, -1)),
             rhs = c(0, 0), neq = 0)
  h2 <- list(constraints = rbind(c(1, -1,  0),
                                  c(0,  1, -1)),
             rhs = c(0, 0), neq = 2)  # alle gelijkheden
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H_order = h1, H_equal = h2))

  expect_equal(nrow(result$result), 3)
  # Bij ware ordening x1>x2>x3 moet ordering-hypothese meer gewicht krijgen
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  weights <- result$result[[wt_col]]
  names(weights) <- result$result$model
  expect_true(weights["H_order"] > weights["H_equal"])
})


# --- Foutafhandeling matrixrestricties ---

test_that("goric.lm matrix: fout bij ongeldige lijstnamen", {
  h_bad <- list(wrong_name = c(1, -1, 0))
  expect_error(
    goric(fit_lm, type = "goric",
          hypotheses = list(H1 = h_bad)),
    "constraints.*rhs.*neq"
  )
})

test_that("goric.numeric matrix: fout bij ongeldige lijstnamen", {
  h_bad <- list(foo = rbind(c(1, -1, 0)))
  expect_error(
    goric(est_3, VCOV = VCOV_3, type = "gorica",
          hypotheses = list(H1 = h_bad)),
    "constraints.*rhs.*neq"
  )
})


# =============================================================================
# goric.CTmeta: meta-analyse (continuous-time) objecten
# =============================================================================

# CTmeta is niet op CRAN, daarom gebruiken we een mock-object dat dezelfde
# interface biedt: coef(), vcov(), en $N / $NrStudies attributen.

# Helper: maak een mock-CTmeta object
make_mock_CTmeta <- function(est, VCOV, N) {
  obj <- list(
    coefficients = est,
    vcov_matrix  = VCOV,
    N            = N,
    NrStudies    = length(N)
  )
  class(obj) <- "CTmeta"
  obj
}

# S3 methoden voor mock-CTmeta (nodig voor goric.CTmeta)
coef.CTmeta <- function(object, ...) object$coefficients
vcov.CTmeta <- function(object, ...) object$vcov_matrix

# Registreer S3 methoden zodat ze gevonden worden tijdens dispatch
# in de testthat child-omgeving (UseMethod zoekt niet in test-env)
registerS3method("coef", "CTmeta", coef.CTmeta)
registerS3method("vcov", "CTmeta", vcov.CTmeta)


# --- Gedeelde CTmeta testdata ---

# Simuleer een meta-analyse met 3 effectmaten
ct_est  <- c(phi1 = 0.5, phi2 = 0.3, phi3 = 0.1)
ct_VCOV <- matrix(c(0.02, 0.005, 0.001,
                     0.005, 0.03, 0.002,
                     0.001, 0.002, 0.025), 3, 3)
ct_N    <- c(50, 60, 45)  # steekproefgroottes per studie

ct_obj <- make_mock_CTmeta(ct_est, ct_VCOV, ct_N)


# --- Basisfunctionaliteit ---

test_that("goric.CTmeta: gorica met 1 hypothese en complement", {
  expect_message(
    result <- goric(ct_obj, type = "gorica",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3")),
    "covariance|vb"
  )

  expect_s3_class(result, "con_goric")
  expect_true("con_gorica" %in% class(result))
  expect_equal(nrow(result$result), 2)
  expect_true("complement" %in% result$result$model)
})

test_that("goric.CTmeta: gorica met meerdere hypothesen", {
  expect_message(
    result <- goric(ct_obj, type = "gorica",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3",
                                      H2 = "phi3 > phi2 > phi1")),
    "covariance|vb"
  )

  expect_equal(nrow(result$result), 3)
  expect_true("unconstrained" %in% result$result$model)
})


# --- Type-conversie ---

test_that("goric.CTmeta: type='goric' wordt omgezet naar gorica", {
  expect_message(
    result <- goric(ct_obj, type = "goric",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3")),
    "GORICA"
  )

  # check.type() converteert goric -> gorica voor CTmeta
  expect_equal(result$type, "gorica")
})

test_that("goric.CTmeta: type='goricc' wordt omgezet naar goricac", {
  expect_message(
    result <- goric(ct_obj, type = "goricc",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3")),
    "GORICAC"
  )

  expect_equal(result$type, "goricac")
})

test_that("goric.CTmeta: goricac gebruikt sum(N) als sample_nobs", {
  expect_message(
    result <- goric(ct_obj, type = "goricac",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3")),
    "vb|covariance"
  )

  expect_equal(result$type, "goricac")
  # N = sum(ct_N) = 155
  expect_equal(result$sample_nobs, sum(ct_N))
})


# --- Gewichten en wiskundige consistentie ---

test_that("goric.CTmeta: gewichten tellen op tot 1", {
  expect_message(
    result <- goric(ct_obj, type = "gorica",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3")),
    "vb|covariance"
  )

  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(sum(result$result[[wt_col]]), 1, tolerance = 1e-6)
})

test_that("goric.CTmeta: ware hypothese krijgt hoogste gewicht", {
  # ct_est: phi1 > phi2 > phi3
  expect_message(
    result <- goric(ct_obj, type = "gorica",
                    hypotheses = list(H_correct = "phi1 > phi2 > phi3",
                                      H_reverse = "phi3 > phi2 > phi1")),
    "vb|covariance"
  )

  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  weights <- result$result[[wt_col]]
  names(weights) <- result$result$model

  expect_true(weights["H_correct"] > weights["H_reverse"])
})

test_that("goric.CTmeta: gorica = -2*loglik + 2*penalty", {
  expect_message(
    result <- goric(ct_obj, type = "gorica", comparison = "none",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3")),
    "vb|covariance"
  )

  row <- result$result[1, ]
  expect_equal(row$gorica, -2 * row$loglik + 2 * row$penalty, tolerance = 1e-8)
})


# --- Comparison varianten ---

test_that("goric.CTmeta: comparison = 'none'", {
  expect_message(
    result <- goric(ct_obj, type = "gorica",
                    comparison = "none",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3",
                                      H2 = "phi1 > phi2 = phi3")),
    "vb|covariance"
  )

  expect_equal(nrow(result$result), 2)
  expect_false("unconstrained" %in% result$result$model)
})


# --- Output structuur ---

test_that("goric.CTmeta: output bevat VCOV en b.unrestr", {
  expect_message(
    result <- goric(ct_obj, type = "gorica",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3")),
    "vb|covariance"
  )

  expect_true(!is.null(result$VCOV))
  expect_true(!is.null(result$b.unrestr))
  expect_equal(length(result$b.unrestr), 3)
  expect_equal(names(result$b.unrestr), names(ct_est))
})


# --- Equivalentie met goric.numeric ---

test_that("goric.CTmeta: resultaten identiek aan goric.numeric met dezelfde input", {
  expect_message(
    res_ct <- goric(ct_obj, type = "gorica", comparison = "none",
                    hypotheses = list(H1 = "phi1 > phi2 > phi3")),
    "vb|covariance"
  )

  res_num <- goric(ct_est, VCOV = ct_VCOV, type = "gorica",
                   comparison = "none",
                   hypotheses = list(H1 = "phi1 > phi2 > phi3"))

  # Loglik, penalty, gorica moeten identiek zijn
  expect_equal(res_ct$result$loglik,  res_num$result$loglik,  tolerance = 1e-6)
  expect_equal(res_ct$result$penalty, res_num$result$penalty, tolerance = 1e-6)
  expect_equal(res_ct$result$gorica,  res_num$result$gorica,  tolerance = 1e-6)
})


# --- Gelijkheidsbeperkingen ---

test_that("goric.CTmeta: gelijkheidsbeperkingen werken", {
  expect_message(
    result <- goric(ct_obj, type = "gorica",
                    hypotheses = list(H1 = "phi1 = phi2 = phi3")),
    "vb|covariance"
  )

  expect_s3_class(result, "con_goric")
  expect_true("unconstrained" %in% result$result$model)
})


# --- Matrixrestricties ---

test_that("goric.CTmeta: matrixrestricties via list(constraints, rhs, neq)", {
  h_mat <- list(constraints = rbind(c(1, -1,  0),
                                     c(0,  1, -1)),
                rhs = c(0, 0),
                neq = 0)

  expect_message(
    result <- goric(ct_obj, type = "gorica",
                    hypotheses = list(H1 = h_mat)),
    "vb|covariance"
  )

  expect_s3_class(result, "con_goric")
  expect_equal(nrow(result$result), 2)
})


# =============================================================================
# Edge cases: 1D probleem (enkele parameter)
# =============================================================================

test_that("edge: gorica met enkele parameter (scalaire VCOV)", {
  est_1 <- c(mu = 3.0)
  vcov_1 <- matrix(0.05, 1, 1)
  result <- goric(est_1, VCOV = vcov_1, type = "gorica",
                  hypotheses = list(H1 = "mu > 0"))

  expect_s3_class(result, "con_goric")
  expect_equal(nrow(result$result), 2)
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(sum(result$result[[wt_col]]), 1, tolerance = 1e-6)
})

test_that("edge: gorica met scalaire VCOV (niet-matrix)", {
  est_1 <- c(mu = 3.0)
  # VCOV als scalar ipv matrix - goric.default converteert dit
  result <- goric(est_1, VCOV = 0.05, type = "gorica",
                  hypotheses = list(H1 = "mu > 0"))

  expect_s3_class(result, "con_goric")
  expect_equal(nrow(result$result), 2)
})


# =============================================================================
# Edge cases: VCOV validatie
# =============================================================================

test_that("edge: fout bij niet-vierkante VCOV", {
  est_2 <- c(x1 = 1, x2 = 2)
  bad_vcov <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(
    goric(est_2, VCOV = bad_vcov, type = "gorica",
          hypotheses = list(H1 = "x1 > x2")),
    "square"
  )
})

test_that("edge: fout bij VCOV met NA waarden", {
  est_2 <- c(x1 = 1, x2 = 2)
  na_vcov <- matrix(c(0.1, 0.01, 0.01, NA), 2, 2)
  expect_error(
    goric(est_2, VCOV = na_vcov, type = "gorica",
          hypotheses = list(H1 = "x1 > x2")),
    "NA|NaN"
  )
})

test_that("edge: fout bij (bijna-)singuliere VCOV", {
  est_2 <- c(x1 = 1, x2 = 2)
  # Exacte singulariteit: rij 2 = rij 1
  sing_vcov <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_error(
    goric(est_2, VCOV = sing_vcov, type = "gorica",
          hypotheses = list(H1 = "x1 > x2")),
    "singular"
  )
})


# =============================================================================
# Edge cases: Heq parameter
# =============================================================================

test_that("edge: Heq = TRUE met >1 hypothese geeft fout", {
  expect_error(
    goric(est_3, VCOV = VCOV_3, type = "gorica",
          Heq = TRUE,
          hypotheses = list(H1 = "x1 > x2", H2 = "x2 > x3")),
    "Heq.*one|one.*hypothesis"
  )
})

test_that("edge: Heq = TRUE met comparison='unconstrained' geeft waarschuwing", {
  expect_warning(
    result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                    comparison = "unconstrained", Heq = TRUE,
                    hypotheses = list(H1 = "x1 > x2 > x3")),
    "Heq.*ignored|ignored"
  )
  # Heq moet FALSE zijn in output
  expect_false(result$Heq)
})

test_that("edge: Heq = TRUE met comparison='complement' voegt Heq-rij toe", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  comparison = "complement", Heq = TRUE,
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  expect_true("Heq" %in% result$result$model)
  # Heq-rij bevat de gelijkheidsversie (x1 = x2 = x3)
  expect_equal(nrow(result$result), 3)  # Heq + H1 + complement
})


# =============================================================================
# Edge cases: onbenoemde en lege hypothesenamen
# =============================================================================

test_that("edge: hypothesen zonder namen krijgen automatische namen", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  comparison = "none",
                  hypotheses = list("x1 > x2", "x2 > x3"))

  # Automatisch genaamd H1, H2
  expect_true("H1" %in% result$result$model)
  expect_true("H2" %in% result$result$model)
})

test_that("edge: hypothesen met lege namen krijgen automatische namen", {
  hyps <- list("x1 > x2", "x2 > x3")
  names(hyps) <- c("", "")
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  comparison = "none",
                  hypotheses = hyps)

  expect_true("H1" %in% result$result$model)
  expect_true("H2" %in% result$result$model)
})


# =============================================================================
# Edge cases: onbekende argumenten
# =============================================================================

test_that("edge: onbekend argument in ... geeft fout", {
  expect_error(
    goric(est_3, VCOV = VCOV_3, type = "gorica",
          onbekend_arg = 42,
          hypotheses = list(H1 = "x1 > x2 > x3")),
    "Unknown|unknown"
  )
})


# =============================================================================
# Edge cases: penalty_factor varianten
# =============================================================================

test_that("edge: penalty_factor = NA geeft fout", {
  expect_error(
    goric(est_3, VCOV = VCOV_3, type = "gorica",
          penalty_factor = NA,
          hypotheses = list(H1 = "x1 > x2 > x3")),
    "penalty factor"
  )
})

test_that("edge: penalty_factor = Inf geeft grote IC-waarden", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  penalty_factor = Inf, comparison = "none",
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  expect_true(is.infinite(result$result$gorica[1]))
})


# =============================================================================
# Edge cases: sample_nobs als vector (groepsgroottes)
# =============================================================================

test_that("edge: sample_nobs als vector wordt gesommeerd (goricac)", {
  expect_message(
    result <- goric(est_3, VCOV = VCOV_3, type = "goricac",
                    sample_nobs = c(30, 30, 30),
                    hypotheses = list(H1 = "x1 > x2 > x3")),
    "group sizes|total"
  )

  expect_equal(result$sample_nobs, 90)
})


# =============================================================================
# Edge cases: complement bij alleen gelijkheidsbeperkingen
# =============================================================================

test_that("edge: complement + alleen gelijkheden schakelt over naar unconstrained", {
  expect_message(
    result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                    comparison = "complement",
                    hypotheses = list(H1 = "x1 = x2 = x3")),
    "unconstrained"
  )

  # Complement van pure gelijkheid = unconstrained
  expect_true("unconstrained" %in% result$result$model)
  expect_false("complement" %in% result$result$model)
})


# =============================================================================
# Edge cases: groot aantal hypothesen
# =============================================================================

test_that("edge: 5 hypothesen tegelijk", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(
                    H1 = "x1 > x2 > x3",
                    H2 = "x1 > x2 = x3",
                    H3 = "x1 = x2 > x3",
                    H4 = "x3 > x2 > x1",
                    H5 = "x1 = x2 = x3"
                  ))

  # 5 hypothesen + unconstrained = 6 rijen
  expect_equal(nrow(result$result), 6)
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(sum(result$result[[wt_col]]), 1, tolerance = 1e-6)

  # H1 (ware ordening) moet het hoogste gewicht hebben
  weights <- result$result[[wt_col]]
  names(weights) <- result$result$model
  expect_true(which.max(weights) == which(result$result$model == "H1"))
})


# =============================================================================
# Edge cases: identieke hypothesen (gedegeneerd geval)
# =============================================================================

test_that("edge: twee identieke hypothesen geven gelijke IC-waarden", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  comparison = "none",
                  hypotheses = list(H1 = "x1 > x2 > x3",
                                    H2 = "x1 > x2 > x3"))

  expect_equal(result$result$gorica[1], result$result$gorica[2], tolerance = 1e-8)
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(result$result[[wt_col]][1], result$result[[wt_col]][2], tolerance = 1e-6)
})


# =============================================================================
# Edge cases: 2D probleem met sterke correlatie
# =============================================================================

test_that("edge: sterk gecorreleerde VCOV werkt", {
  est_2 <- c(a = 3, b = 1)
  # Hoge correlatie (r = 0.9)
  vcov_hc <- matrix(c(0.1, 0.09, 0.09, 0.1), 2, 2)
  result <- goric(est_2, VCOV = vcov_hc, type = "gorica",
                  hypotheses = list(H1 = "a > b"))

  expect_s3_class(result, "con_goric")
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(sum(result$result[[wt_col]]), 1, tolerance = 1e-6)
})


# =============================================================================
# Edge cases: data op de grens (boundary solution)
# =============================================================================

test_that("edge: schattingen exact op de grens (x1 = x2)", {
  # x1 == x2, hypothese x1 > x2 => boundary
  est_boundary <- c(x1 = 3.0, x2 = 3.0, x3 = 1.0)
  result <- goric(est_boundary, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  expect_s3_class(result, "con_goric")
  # Model moet nog steeds geldige IC-waarden geven
  expect_true(all(is.finite(result$result$gorica)))
})

test_that("edge: schattingen die hypothese schenden (x3 > x1)", {
  # Schattingen schenden H1: x1 > x2 > x3
  est_violated <- c(x1 = 1.0, x2 = 2.0, x3 = 3.0)
  result <- goric(est_violated, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2 > x3"))

  expect_s3_class(result, "con_goric")
  # complement moet meer gewicht krijgen
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  weights <- result$result[[wt_col]]
  names(weights) <- result$result$model
  expect_true(weights["complement"] > weights["H1"])
})


# =============================================================================
# Edge cases: goric.lm specifiek
# =============================================================================

test_that("edge: goric.lm type='goricc' small sample correctie", {
  result_c  <- goric(fit_lm, type = "goricc",
                     hypotheses = list(H1 = "group1 > group2 > group3"))
  result_nc <- goric(fit_lm, type = "goric",
                     hypotheses = list(H1 = "group1 > group2 > group3"))

  # goricc penalty >= goric penalty (small sample correctie)
  expect_true(result_c$result$penalty[1] >= result_nc$result$penalty[1] - 1e-6)
  expect_equal(result_c$type, "goricc")
})

test_that("edge: goric.lm met model met intercept", {
  fit_int <- lm(y ~ group, data = df_anova)
  result <- goric(fit_int, type = "goric",
                  hypotheses = list(H1 = "group2 > 0; group3 > group2"))

  expect_s3_class(result, "con_goric")
})

test_that("edge: goric.lm comparison='complement' met >1 hypo -> unconstrained", {
  expect_warning(
    result <- goric(fit_lm, type = "goric",
                    comparison = "complement",
                    hypotheses = list(H1 = "group1 > group2 > group3",
                                      H2 = "group1 > group2 = group3")),
    "unconstrained"
  )

  expect_true("unconstrained" %in% result$result$model)
  expect_false("complement" %in% result$result$model)
})


# =============================================================================
# Edge cases: grote VCOV (veel parameters)
# =============================================================================

test_that("edge: 6 parameters, subset in hypothese", {
  est_6 <- c(a = 5, b = 4, c = 3, d = 2, e = 1, f = 0)
  vcov_6 <- diag(6) * 0.05

  # Hypothese over subset van parameters
  result <- goric(est_6, VCOV = vcov_6, type = "gorica",
                  hypotheses = list(H1 = "a > b > c"))

  expect_s3_class(result, "con_goric")
  wt_col <- grep("weights", names(result$result), value = TRUE)[1]
  expect_equal(sum(result$result[[wt_col]]), 1, tolerance = 1e-6)
})


# =============================================================================
# Edge cases: semicolon-gescheiden restricties
# =============================================================================

test_that("edge: hypothese met semicolons als scheidingsteken", {
  result <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                  hypotheses = list(H1 = "x1 > x2; x2 > x3"))

  expect_s3_class(result, "con_goric")
  # Moet equivalent zijn aan x1 > x2 > x3
  res_chain <- goric(est_3, VCOV = VCOV_3, type = "gorica",
                     hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_equal(result$result$gorica[1], res_chain$result$gorica[1], tolerance = 1e-6)
})
