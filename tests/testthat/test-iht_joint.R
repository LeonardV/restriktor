# =============================================================================
# Tests: iht_joint(), pchibar_joint(), pfbar_joint()
#
# Wolak (1989, J. Econometrics 41, Theorem 1): gezamenlijke verdeling van de
# type-A- (EI) en type-B- (IU) teststatistiek onder de gedeelde (least
# favorable) nulhypothese Rb = r:
#
#   Pr[IU >= c1, EI >= c2] = sum_k w(L,k,Omega) Pr[chi2_{q-k} >= c1]
#                                              Pr[chi2_k >= c2]
#
# Opzet volgt test-conTest.R: restriktor() met mix_weights = "pmvnorm" zodat
# wt.bar 1x berekend en gecachet wordt; alle assertions op gecachte resultaten.
# =============================================================================

# --- Gedeelde testdata ---
set.seed(42)
n <- 60
group <- factor(rep(1:3, each = n / 3))
x1 <- rnorm(n)
y <- 2 + 1.5 * (group == "2") + 3 * (group == "3") + rnorm(n)
df_test <- data.frame(y = y, group = group, x1 = x1)

fit_lm  <- lm(y ~ -1 + group, data = df_test)
fit_glm <- glm(y ~ -1 + group, family = gaussian, data = df_test)
fit_rlm <- MASS::rlm(y ~ -1 + group, data = df_test, method = "MM")

con_lm  <- restriktor(fit_lm,  constraints = "group1 < group2 < group3")
con_glm <- restriktor(fit_glm, constraints = "group1 < group2 < group3")
con_rlm <- restriktor(fit_rlm, constraints = "group1 < group2 < group3")
con_lm_ceq <- restriktor(fit_lm, constraints = "group1 = group2; group2 = group3")

out_joint <- iht_joint(con_lm)
wt_lm  <- out_joint$wt.bar
q_lm   <- nrow(out_joint$Amat)
df2_lm <- out_joint$df.residual


# =============================================================================
# (a) Numerieke consistentie: randgevallen en monotonie
# =============================================================================

test_that("pchibar_joint/pfbar_joint: randgevallen", {
  expect_equal(pchibar_joint(0, 0, q_lm, wt_lm), 1)
  expect_equal(pchibar_joint(Inf, Inf, q_lm, wt_lm), 0)
  expect_equal(pfbar_joint(0, 0, q_lm, df2_lm, wt_lm), 1)
  expect_equal(pfbar_joint(Inf, Inf, q_lm, df2_lm, wt_lm), 0)
  # df2 = Inf reduceert de F-bar-versie tot de chi-bar-versie
  expect_equal(pfbar_joint(2.5, 1.5, q_lm, Inf, wt_lm),
               pchibar_joint(2.5, 1.5, q_lm, wt_lm))
})

test_that("pfbar_joint: numeriek stabiel voor grote df2, convergeert naar chi-bar", {
  wt <- c(1/3, 1/2, 1/6)  # q = 2, corr(Omega) = -0.5
  for (df2 in c(1000, 1e5, 1e7)) {
    expect_equal(pfbar_joint(0, 4.9, 2, df2, wt),
                 1 - restriktor:::pfbar(4.9, df1 = 0:2, df2 = df2, wt.bar = wt),
                 tolerance = 1e-08)
  }
  expect_equal(pfbar_joint(5.4, 4.9, 2, 1e7, wt),
               pchibar_joint(5.4, 4.9, 2, wt), tolerance = 1e-05)
})

test_that("pchibar_joint/pfbar_joint: monotoon dalend in beide argumenten", {
  grid <- c(0, 0.5, 1, 2, 5, 10)
  for (fixed in c(0, 1.3)) {
    pc1 <- sapply(grid, function(cc) pchibar_joint(cc, fixed, q_lm, wt_lm))
    pc2 <- sapply(grid, function(cc) pchibar_joint(fixed, cc, q_lm, wt_lm))
    pf1 <- sapply(grid, function(cc) pfbar_joint(cc, fixed, q_lm, df2_lm, wt_lm))
    pf2 <- sapply(grid, function(cc) pfbar_joint(fixed, cc, q_lm, df2_lm, wt_lm))
    expect_true(all(diff(pc1) <= 1e-12))
    expect_true(all(diff(pc2) <= 1e-12))
    expect_true(all(diff(pf1) <= 1e-10))
    expect_true(all(diff(pf2) <= 1e-10))
  }
})


# =============================================================================
# (b) Marginalisatie: joint met een grens op 0 = bestaande marginale functies
# =============================================================================

test_that("marginalisatie: pchibar_joint consistent met pchibar", {
  meq <- con_lm$neq  # 0 hier
  L <- q_lm - meq
  for (cc in c(0.3, 1.2, 3.7, 8, 15)) {
    # cIU = 0: joint reduceert tot marginale type-A-staartkans
    expect_equal(pchibar_joint(0, cc, q_lm, wt_lm),
                 1 - restriktor:::pchibar(cc, df = 0:L, wt.bar = wt_lm),
                 tolerance = 1e-12)
    # cEI = 0: joint reduceert tot marginale type-B-staartkans
    expect_equal(pchibar_joint(cc, 0, q_lm, wt_lm),
                 1 - restriktor:::pchibar(cc, df = meq:q_lm, wt.bar = rev(wt_lm)),
                 tolerance = 1e-12)
  }
})

test_that("marginalisatie: pfbar_joint consistent met pfbar", {
  meq <- con_lm$neq
  L <- q_lm - meq
  # tolerance is relatief; de kwadratuur in pfbar_joint (rel.tol 1e-09) geeft
  # in verre staarten (kansen ~1e-4) relatieve afwijkingen tot ~1e-07
  for (cc in c(0.3, 1.2, 3.7, 8, 15)) {
    expect_equal(pfbar_joint(0, cc, q_lm, df2_lm, wt_lm),
                 1 - restriktor:::pfbar(cc, df1 = 0:L, df2 = df2_lm,
                                        wt.bar = wt_lm),
                 tolerance = 1e-06)
    expect_equal(pfbar_joint(cc, 0, q_lm, df2_lm, wt_lm),
                 1 - restriktor:::pfbar(cc, df1 = meq:q_lm, df2 = df2_lm,
                                        wt.bar = rev(wt_lm)),
                 tolerance = 1e-06)
  }
})

test_that("iht_joint: statistieken en marginale p-waarden identiek aan iht", {
  outA <- iht(con_lm, type = "A", test = "F")
  outB <- iht(con_lm, type = "B", test = "F")
  expect_equal(out_joint$Ts.A, outA$Ts, tolerance = 1e-10)
  expect_equal(out_joint$Ts.B, outB$Ts, tolerance = 1e-10)
  expect_equal(as.numeric(out_joint$pvalue.A), as.numeric(outA$pvalue),
               tolerance = 1e-10)
  expect_equal(as.numeric(out_joint$pvalue.B), as.numeric(outB$pvalue),
               tolerance = 1e-10)
  # joint p-waarde is de staartkans op de waargenomen (Ts.B, Ts.A)
  expect_equal(out_joint$pvalue.joint,
               pfbar_joint(out_joint$Ts.B, out_joint$Ts.A, q_lm, df2_lm, wt_lm),
               tolerance = 1e-12)
  expect_true(out_joint$pvalue.joint >= 0 && out_joint$pvalue.joint <= 1)
  # joint staartkans kan niet groter zijn dan elk van de marginale staartkansen
  expect_lte(out_joint$pvalue.joint, as.numeric(out_joint$pvalue.A) + 1e-08)
  expect_lte(out_joint$pvalue.joint, as.numeric(out_joint$pvalue.B) + 1e-08)
})

test_that("iht_joint: werkt ook direct op een ongefit lm-object", {
  out2 <- iht_joint(fit_lm, constraints = "group1 < group2 < group3")
  expect_equal(out2$Ts.A, out_joint$Ts.A, tolerance = 1e-10)
  expect_equal(out2$Ts.B, out_joint$Ts.B, tolerance = 1e-10)
  expect_equal(out2$pvalue.joint, out_joint$pvalue.joint, tolerance = 1e-10)
})


# =============================================================================
# (d) Onafhankelijke referentie-implementatie: P = 3, niet-diagonale Omega
#     Gewichten via de subset-decompositie met gesloten-vorm orthantkansen
#     (arcsin-formule, Kudo 1963), onafhankelijk van ic_weights().
# =============================================================================

# positieve-orthantkans P(Z > 0), Z ~ N(0, S), gesloten vorm voor dim 1-3
orthant_ref <- function(S) {
  k <- NROW(S)
  if (k == 0L) return(1)
  if (k == 1L) return(0.5)
  r <- cov2cor(S)
  if (k == 2L) return(0.25 + asin(r[1, 2]) / (2 * pi))
  0.125 + (asin(r[1, 2]) + asin(r[1, 3]) + asin(r[2, 3])) / (4 * pi)
}

# level-kansen w_j = Pr[EI-component ~ chi2_j], j = 0..g (Kudo 1963;
# Silvapulle & Sen 2005, eq. 3.23), W = Amat %*% Sigma %*% t(Amat):
#   w_j = sum_{|S| = g - j} P(Z_S > 0 | (W_SS)^-1) *
#         P(Z_Sbar > 0 | W_SbarSbar - W_SbarS (W_SS)^-1 W_SSbar)
wt_bar_ref <- function(W) {
  g <- nrow(W)
  sapply(0:g, function(j) {
    k <- g - j  # subsetgrootte
    if (k == 0L) return(orthant_ref(W))
    subsets <- utils::combn(g, k)
    total <- 0
    for (s in seq_len(ncol(subsets))) {
      S  <- subsets[, s]
      Sc <- setdiff(seq_len(g), S)
      WSS_inv <- solve(W[S, S, drop = FALSE])
      schur <- W[Sc, Sc, drop = FALSE] -
        W[Sc, S, drop = FALSE] %*% WSS_inv %*% W[S, Sc, drop = FALSE]
      total <- total + orthant_ref(WSS_inv) * orthant_ref(schur)
    }
    total
  })
}

test_that("P = 3, niet-diagonale Omega: wt.bar en joint vs referentie", {
  set.seed(123)
  n3 <- 80
  # gecorreleerde predictoren -> niet-diagonale Omega = Amat Sigma Amat'
  z1 <- rnorm(n3)
  z2 <- 0.6 * z1 + sqrt(1 - 0.36) * rnorm(n3)
  z3 <- -0.4 * z1 + 0.3 * z2 + rnorm(n3)
  y3 <- 0.4 * z1 + 0.3 * z2 + 0.2 * z3 + rnorm(n3)
  d3 <- data.frame(y3, z1, z2, z3)
  fit3 <- lm(y3 ~ z1 + z2 + z3, data = d3)
  con3 <- restriktor(fit3, constraints = "z1 > 0; z2 > 0; z3 > 0")

  Amat3 <- con3$constraints
  Sigma3 <- vcov(fit3)
  W <- Amat3 %*% Sigma3 %*% t(Amat3)
  expect_true(max(abs(cov2cor(W)[upper.tri(W)])) > 0.05)  # echt niet-diagonaal

  # 1. gewichten: gecachete wt.bar vs onafhankelijke referentie
  wt_ref <- wt_bar_ref(W)
  expect_equal(sum(wt_ref), 1, tolerance = 1e-10)
  expect_equal(as.numeric(con3$wt.bar), wt_ref, tolerance = 1e-06)

  # 2. joint staartkans: pchibar_joint vs directe som met referentiegewichten
  joint_ref <- function(cIU, cEI) {
    g <- nrow(W)
    sum(sapply(0:g, function(k) {
      tIU <- if (g - k == 0) as.numeric(cIU <= 0) else
        pchisq(cIU, g - k, lower.tail = FALSE)
      tEI <- if (k == 0) as.numeric(cEI <= 0) else
        pchisq(cEI, k, lower.tail = FALSE)
      wt_ref[k + 1] * tIU * tEI
    }))
  }
  for (cc in list(c(0, 0), c(1.5, 0.8), c(4.2, 2.9), c(9, 6), c(0.3, 11))) {
    expect_equal(pchibar_joint(cc[1], cc[2], 3, con3$wt.bar),
                 joint_ref(cc[1], cc[2]), tolerance = 1e-06)
  }

  # 3. iht_joint op dezelfde fit gebruikt dezelfde gewichten
  outj3 <- iht_joint(con3)
  expect_equal(outj3$wt.bar, wt_ref, tolerance = 1e-06)
  expect_true(outj3$pvalue.joint >= 0 && outj3$pvalue.joint <= 1)
})


# =============================================================================
# glm / rlm: rooktests
# =============================================================================

test_that("iht_joint werkt voor conGLM en conRLM", {
  outg <- iht_joint(con_glm)
  outr <- iht_joint(con_rlm, test = "Wald")
  for (out in list(outg, outr)) {
    expect_s3_class(out, "conTestJoint")
    expect_true(out$pvalue.joint >= 0 && out$pvalue.joint <= 1)
    expect_true(out$Ts.A >= -1e-10 && out$Ts.B >= -1e-10)
    expect_true(out$decision %in% c("a", "b", "c"))
  }
})

test_that("print.conTestJoint geeft leesbare output", {
  txt <- capture.output(print(out_joint))
  expect_true(any(grepl("Wolak", txt)))
  expect_true(any(grepl("Type A test", txt)))
  expect_true(any(grepl("Type B test", txt)))
  expect_true(any(grepl("Joint tail probability", txt)))
  expect_true(any(grepl("Decision rule", txt)))
})


# =============================================================================
# Integratie in conTest_summary / iht(type = "summary")
# =============================================================================

test_that("iht(type = 'summary') bevat de joint p-waarde (Wolak)", {
  s <- iht(con_lm)  # type = "summary" is de default
  expect_false(is.null(s$joint))
  expect_equal(s$joint$pvalue, out_joint$pvalue.joint, tolerance = 1e-10)
  expect_equal(s$joint$Ts.A, out_joint$Ts.A, tolerance = 1e-10)
  expect_equal(s$joint$Ts.B, out_joint$Ts.B, tolerance = 1e-10)
  # print toont de joint-regel
  txt <- capture.output(print(s))
  expect_true(any(grepl("Joint Type A and Type B test \\(Wolak", txt)))
})

test_that("summary zonder chi-bar-gewichten of met bootstrap: geen joint", {
  con_nw <- restriktor(fit_lm, constraints = "group1 < group2 < group3",
                       mix_weights = "none")
  s_nw <- iht(con_nw, boot = "parametric", R = 9)
  expect_true(is.null(s_nw$joint))
})


# =============================================================================
# Beslisregel en foutafhandeling
# =============================================================================

test_that("beslisregel: regio's a/b/c", {
  # duidelijk oplopende gemiddelden: type A verwerpt, type B niet -> regio (b)
  expect_identical(out_joint$decision, "b")
  # geen effect: type A verwerpt niet -> regio (a)
  set.seed(7)
  y0 <- rnorm(n)
  fit0 <- lm(y0 ~ -1 + group, data = data.frame(y0, group))
  out0 <- iht_joint(fit0, constraints = "group1 < group2 < group3")
  expect_identical(out0$decision, "a")
  # sterk dalende gemiddelden: type B verwerpt -> regio (c)
  set.seed(8)
  yc <- 6 - 3 * as.numeric(group) + rnorm(n)
  fitc <- lm(yc ~ -1 + group, data = data.frame(yc, group))
  outc <- iht_joint(fitc, constraints = "group1 < group2 < group3")
  expect_identical(outc$decision, "c")
})

test_that("iht_joint foutafhandeling", {
  # alleen gelijkheidsrestricties
  expect_error(iht_joint(con_lm_ceq), "equality restrictions only")
  # neq.alt > 0 niet ondersteund
  expect_error(iht_joint(con_lm, neq.alt = 1), "neq.alt")
  # bootstrap p-waarden niet ondersteund
  expect_error(iht_joint(con_lm, boot = "parametric"), "boot")
  # geen constraints
  expect_error(iht_joint(fit_lm), "no constraints")
  # geen gewichten beschikbaar
  con_nw <- restriktor(fit_lm, constraints = "group1 < group2 < group3",
                       mix_weights = "none")
  expect_error(iht_joint(con_nw), "mix_weights")
  # ongeldig alpha
  expect_error(iht_joint(con_lm, alpha = 1.5), "alpha")
})


# =============================================================================
# (c) Monte-Carlo type-I-fout onder de gedeelde nul Rb = r
#
# >= 50.000 simulaties; duurt ~10-15 minuten en wordt daarom alleen gedraaid
# als RESTRIKTOR_RUN_MC=yes is gezet (en nooit op CRAN):
#   RESTRIKTOR_RUN_MC=yes Rscript -e 'testthat::test_local(filter="iht_joint")'
# =============================================================================

test_that("Monte-Carlo: marginale en gezamenlijke type-I-fout volgen Wolak", {
  skip_on_cran()
  skip_if(Sys.getenv("RESTRIKTOR_RUN_MC") != "yes",
          "zet RESTRIKTOR_RUN_MC=yes voor de Monte-Carlo type-I-fout test")

  nsim <- 50000L
  set.seed(2026)
  n_mc <- 60L
  group_mc <- factor(rep(1:3, each = n_mc / 3))
  Amat_mc <- rbind(c(-1, 1, 0),
                   c( 0,-1, 1))
  q_mc <- nrow(Amat_mc)

  TsA <- TsB <- pA <- pB <- numeric(nsim)
  wt_mc <- NULL
  df2_mc <- NULL
  for (i in seq_len(nsim)) {
    # DGP onder Rb = r: alle groepsgemiddelden gelijk (least favorable nul)
    y_mc <- rnorm(n_mc)
    fit_mc <- lm(y_mc ~ -1 + group_mc)
    out <- iht_joint(fit_mc, constraints = Amat_mc, rhs = c(0, 0), neq = 0)
    TsA[i] <- out$Ts.A
    TsB[i] <- out$Ts.B
    pA[i]  <- out$pvalue.A
    pB[i]  <- out$pvalue.B
    if (i == 1L) {
      wt_mc  <- out$wt.bar   # vaste X -> Omega en dus wt.bar identiek per sim
      df2_mc <- out$df.residual
    }
  }

  mc_se <- function(p) sqrt(p * (1 - p) / nsim)

  # -- marginale verwerpingsfrequenties = nominaal alpha (binnen 4 MC-se) --
  for (a in c(0.10, 0.05, 0.01)) {
    expect_lt(abs(mean(pA < a) - a), 4 * mc_se(a))
    expect_lt(abs(mean(pB < a) - a), 4 * mc_se(a))
  }

  # -- gezamenlijke verwerpingsfrequentie = pfbar_joint, != alpha_A * alpha_B --
  for (a in c(0.10, 0.05)) {
    # marginale kritieke waarden bij nominale alpha
    cA <- uniroot(function(cc) {
      (1 - restriktor:::pfbar(cc, df1 = 0:q_mc, df2 = df2_mc, wt.bar = wt_mc)) - a
    }, interval = c(1e-08, 500))$root
    cB <- uniroot(function(cc) {
      (1 - restriktor:::pfbar(cc, df1 = 0:q_mc, df2 = df2_mc,
                              wt.bar = rev(wt_mc))) - a
    }, interval = c(1e-08, 500))$root

    emp_joint <- mean(TsB >= cB & TsA >= cA)
    theo_joint <- pfbar_joint(cB, cA, q_mc, df2_mc, wt_mc)
    naive_prod <- a * a

    # empirisch consistent met Wolak Theorem 1 (F-bar-vorm)
    expect_lt(abs(emp_joint - theo_joint), 4 * mc_se(max(theo_joint, 1e-04)))
    # het naieve product van de marginale niveaus wijkt systematisch af ...
    expect_gt(abs(emp_joint - naive_prod), 8 * mc_se(max(theo_joint, 1e-04)))
    # ... en de joint-functie zit er aantoonbaar dichter op
    expect_lt(abs(emp_joint - theo_joint), abs(emp_joint - naive_prod))
  }
})
