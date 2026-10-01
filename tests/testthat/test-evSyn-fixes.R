# Tests voor bugfixes in evSyn() (comparison-default, priorWeights-deprecatie,
# priorICweights per studie, study_weights & volgorde, LL weights,
# leave1studyout, order_studies voor ICweights/ICratios, Href, foutmeldingen).

est_fx <- list(c(x1 = 5.0, x2 = 3.0, x3 = 1), c(x1 = 4.5, x2 = 2.5, x3 = 3),
               c(x1 = 6.0, x2 = 2.0, x3 = 2.2), c(x1 = 3, x2 = 3.2, x3 = 1))
V_fx   <- list(diag(3) * .10, diag(3) * .15, diag(3) * .12, diag(3) * .3)
H2_fx  <- list(H1 = "x1 > x2 > x3", H2 = "x1 > x2 & x1 > x3")
H1_fx  <- list(H1 = "x1 > x2 > x3")

g_fx  <- lapply(1:4, function(i) goric(est_fx[[i]], VCOV = V_fx[[i]],
                                       hypotheses = H2_fx, type = "gorica"))
IC_fx <- lapply(g_fx, function(x) x$result$gorica)
W_fx  <- lapply(g_fx, function(x) x$result$gorica.weights)
R_fx  <- lapply(W_fx, function(w) w / w[length(w)])
sw_fx <- c(.4, .3, .2, .1)
pr_fx <- c(.6, .3, .1)

final_w <- function(x) unname(x$Cumulative_GORICA_weights["Final", ])

# 1. comparison default --------------------------------------------------------
test_that("evSyn_est: één hypothese -> standaard vergelijking met complement", {
  res <- evSyn(list(c(x1 = 3, x2 = 2, x3 = 1), c(x1 = 2.5, x2 = 2, x3 = .5)),
               VCOV = list(diag(3) * .1, diag(3) * .1),
               hypotheses = list(H1 = "x1 > x2 > x3"))
  expect_equal(colnames(res$Cumulative_GORICA_weights), c("H1", "Complement"))
  # waarden zoals in master
  expect_equal(final_w(res), c(0.991769685385915, 0.00823031461408453),
               tolerance = 1e-8)
  expect_equal(unname(res$Cumulative_LL_weights["Final", ]),
               c(0.957912272084381, 0.0420877279156189), tolerance = 1e-8)
})

test_that("evSyn_est: lijst-in-lijst met één hypothese per set -> complement", {
  res_ll <- evSyn(est_fx, VCOV = V_fx, hypotheses = rep(list(H1_fx), 4))
  res    <- evSyn(est_fx, VCOV = V_fx, hypotheses = H1_fx)
  expect_equal(colnames(res_ll$Cumulative_GORICA_weights), c("H1", "Complement"))
  expect_equal(res_ll$Cumulative_GORICA_weights, res$Cumulative_GORICA_weights)
})

test_that("evSyn_est: expliciete comparison wordt gerespecteerd", {
  res <- evSyn(est_fx, VCOV = V_fx, hypotheses = H1_fx, comparison = "unconstrained")
  expect_equal(colnames(res$Cumulative_GORICA_weights), c("H1", "Unconstrained"))
  res2 <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx)
  expect_equal(colnames(res2$Cumulative_GORICA_weights), c("H1", "H2", "Unconstrained"))
})

# 2. priorWeights deprecatie --------------------------------------------------
test_that("priorWeights geeft deprecatie-waarschuwing en werkt als priorICweights", {
  expect_warning(
    res_old <- evSyn(W_fx, priorWeights = pr_fx),
    "argument 'priorWeights' is deprecated; use 'priorICweights'", fixed = TRUE
  )
  res_new <- evSyn(W_fx, priorICweights = pr_fx)
  expect_equal(res_old, res_new)

  expect_warning(res_old_r <- evSyn(R_fx, priorWeights = pr_fx), "deprecated")
  expect_equal(res_old_r, evSyn(R_fx, priorICweights = pr_fx))

  expect_warning(res_old_ic <- evSyn(IC_fx, priorWeights = pr_fx), "deprecated")
  expect_equal(res_old_ic, evSyn(IC_fx, priorICweights = pr_fx))
})

test_that("priorWeights en priorICweights samen geeft een fout", {
  expect_error(evSyn(W_fx, priorWeights = pr_fx, priorICweights = pr_fx),
               "priorWeights")
})

# 3. priorICweights per studie (S != H) -----------------------------------------
test_that("studiespecifieke gewichten met priorICweights zijn correct als S != H", {
  expected <- t(vapply(W_fx, function(w) pr_fx * w / sum(pr_fx * w), numeric(3)))
  resW  <- evSyn(W_fx, priorICweights = pr_fx)
  resIC <- evSyn(IC_fx, priorICweights = pr_fx)
  expect_equal(unname(resW$GORICA_weight_m), unname(expected))
  expect_equal(unname(resIC$GORICA_weight_m), unname(expected))
  expect_equal(unname(resW$GORICA_weight_m[2, ]), c(0.61645448, 0.33501344, 0.04853208),
               tolerance = 1e-6)
  expect_equal(unname(rowSums(resW$GORICA_weight_m)), rep(1, 4))
})

# 4. study_weights & volgorde -------------------------------------------------
test_that("resultaat met study_weights hangt niet af van de volgorde van studies", {
  for (obj in list(IC_fx, W_fx, R_fx)) {
    a <- evSyn(obj, study_weights = sw_fx)
    b <- evSyn(obj, study_weights = sw_fx, order_studies = 4:1)
    cc <- evSyn(obj[4:1], study_weights = sw_fx[4:1])
    expect_equal(final_w(a), final_w(b))
    expect_equal(final_w(a), final_w(cc))
    expect_equal(b$study_weights, sw_fx[4:1])
  }
  a <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx, study_weights = sw_fx)
  b <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx, study_weights = sw_fx,
             order_studies = "descending")
  expect_equal(final_w(a), final_w(b))
  expect_equal(unname(a$Cumulative_LL_weights["Final", ]),
               unname(b$Cumulative_LL_weights["Final", ]))
  a <- evSyn(g_fx, study_weights = sw_fx)
  b <- evSyn(g_fx, study_weights = sw_fx, order_studies = c(2, 4, 1, 3))
  expect_equal(final_w(a), final_w(b))
})

test_that("ratio_GORICA_weight_mu wordt mee herordend", {
  a <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx)
  b <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx, order_studies = 4:1)
  expect_equal(unname(b$ratio_GORICA_weight_mu), unname(a$ratio_GORICA_weight_mu[4:1, ]))
})

# 5. LL weights en leave1studyout met study_weights -----------------------------
test_that("LL weights reageren op study_weights", {
  e0 <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx)
  e1 <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx, study_weights = sw_fx)
  expect_false(isTRUE(all.equal(e0$Cumulative_LL_weights, e1$Cumulative_LL_weights)))
  expect_false(isTRUE(all.equal(e0$Final_ratio_LL_weights, e1$Final_ratio_LL_weights)))
  # gewogen som van LL (gewichten sommeren tot S)
  ll <- colSums(e1$LL_m * 4 * sw_fx)
  expect_equal(unname(e1$Cumulative_LL_weights["Final", ]),
               unname(exp(ll - max(ll)) / sum(exp(ll - max(ll)))))
  # LL-route en gorica-route geven hetzelfde
  LL <- lapply(g_fx, function(x) x$result$loglik)
  PT <- lapply(g_fx, function(x) x$result$penalty)
  eL <- evSyn(LL, PT = PT, study_weights = sw_fx)
  expect_equal(unname(eL$Cumulative_LL_weights), unname(e1$Cumulative_LL_weights))
})

test_that("leave1studyout is consistent met evSyn zonder de weggelaten studie", {
  for (te in c("added", "equal", "average")) {
    e1 <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx, study_weights = sw_fx,
                priorICweights = pr_fx, type_ev = te)
    l1 <- leave1studyout(e1)
    for (s in 1:4) {
      ref <- suppressMessages(
        evSyn(est_fx[-s], VCOV = V_fx[-s], hypotheses = H2_fx,
              study_weights = sw_fx[-s], priorICweights = pr_fx, type_ev = te)
      )
      expect_equal(unname(l1$OverallGoricaWeights[s, ]), final_w(ref))
      expect_equal(unname(l1$OverallGorica[s, ]), unname(ref$Cumulative_GORICA["Final", ]))
    }
  }
})

test_that("leave1studyout zonder study_weights is ongewijzigd", {
  e0 <- evSyn(IC_fx)
  l0 <- leave1studyout(e0)
  expect_equal(unname(l0$OverallGorica[1, ]), unname(colSums(do.call(rbind, IC_fx)[-1, ])))
})

# 6. order_studies ascending/descending voor ICweights / ICratios ----------------
test_that("order_studies = 'ascending'/'descending' werkt voor ICweights en ICratios", {
  for (o in c("ascending", "descending")) {
    rW <- evSyn(W_fx, order_studies = o)
    rR <- evSyn(R_fx, order_studies = o)
    rI <- evSyn(IC_fx, order_studies = o)
    expect_equal(rW$order_studies, rI$order_studies)
    expect_equal(rR$order_studies, rI$order_studies)
    expect_equal(final_w(rW), final_w(evSyn(W_fx)))
  }
  w1 <- vapply(W_fx, function(w) w[1], numeric(1))
  expect_equal(evSyn(W_fx, order_studies = "ascending")$order_studies, order(w1))
})

# 7. Href zonder gemeenschappelijke referentiehypothese --------------------------
test_that("ICratios zonder gemeenschappelijke referentie: Href gevonden en input herschaald", {
  expect_message(
    res <- evSyn(list(c(1, 1, .5), c(.5, .2, 1), c(2, 1, 3))),
    "Not all studies used the same reference hypothesis"
  )
  expect_equal(unname(res$Href), 1)
  expect_equal(names(res$Href), "H1")
  # ratio's zijn nu t.o.v. Href = H1 in alle studies
  expect_equal(unname(res$GwRatio_m[, 1]), rep(1, 3))
  expect_equal(unname(res$GwRatio_m[2, ]), c(1, .4, 2))
  expect_equal(unname(res$ICdiff_m[, 1]), rep(0, 3))
  expect_false(anyNA(res$Cumulative_GORICA_weights))
  # finale gewichten onafhankelijk van referentie (zelfde als via IC weights)
  W <- lapply(list(c(1, 1, .5), c(.5, .2, 1), c(2, 1, 3)), function(x) x / sum(x))
  expect_equal(final_w(res), final_w(evSyn(W)))
})

test_that("ICratios: gemeenschappelijke referentie wordt gevonden (ook niet de eerste 1)", {
  res <- evSyn(list(c(1, 2, 1), c(.5, 3, 1), c(2, 1, 1)))
  expect_equal(unname(res$Href), 3)
})

# 8. Labels in de summary van ICratios -------------------------------------------
test_that("ICratios summary: juiste referentiehypothese en studienamen", {
  res <- suppressMessages(evSyn(list(c(1, 1, .5), c(.5, .2, 1), c(2, 1, 3)),
                                study_names = c("A", "B", "C")))
  expect_equal(rownames(res$GORICA_weight_m), c("A", "B", "C"))
  out <- capture.output(print(summary(res)))
  expect_true(any(grepl("(versus reference hypothesis H1)", out, fixed = TRUE)))
  # benoemde input: label volgt kolomnamen (H1, H2, ...)
  resn <- evSyn(lapply(R_fx, function(x) setNames(x, c("a", "b", "c"))))
  outn <- capture.output(print(summary(resn)))
  expect_true(any(grepl("(versus reference hypothesis H3)", outn, fixed = TRUE)))
})

# 9. Foutmeldingen ------------------------------------------------------------
test_that("ICweights validatie geeft de juiste foutmelding", {
  expect_error(restriktor:::evSyn_ICweights(list(c(.5, .6, -.1), c(.2, .3, .5))),
               "Not all values are >= 0")
  err <- tryCatch(restriktor:::evSyn_ICweights(list(c(.5, .6, -.1), c(.2, .3, .5))),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("<= 1", err))
  expect_false(grepl("do not sum to 1", err))
  expect_error(restriktor:::evSyn_ICweights(list(c(.5, .6, .1), c(.2, .3, .5))),
               "do not sum to 1")
  expect_error(restriktor:::evSyn_ICweights(list(c(1.5, -.6, .1), c(.2, .3, .5))),
               "Not all values are <= 1")
})

test_that("ICvalues hypo_names foutmelding noemt het juiste aantal", {
  expect_error(evSyn(IC_fx, hypo_names = c("a", "b")),
               "should consist of 3 names")
})

# =============================================================================
# Fixes n.a.v. externe review: identiteit van hypothesen, validatie van
# gewichten (incl. studiegewicht 0), underflow bij producten van IC weights,
# S = 1, voorkeurshypothese o.b.v. finale gewichten, penalty_factor.
# =============================================================================

# 10. Identiteit van hypothesen bij goric-objecten ------------------------------
test_that("evSyn_gorica: hypothesen worden op naam gekoppeld, niet op positie", {
  g1 <- goric(c(x = 2), VCOV = matrix(1), hypotheses = list(Hp = "x > 0", Hn = "x < 0"),
              type = "gorica")
  g2 <- goric(c(x = 2), VCOV = matrix(1), hypotheses = list(Hn = "x < 0", Hp = "x > 0"),
              type = "gorica")
  expect_message(res <- evSyn(list(g1, g2)), "matched by name")
  ref <- evSyn(list(g1, g1))
  expect_equal(colnames(res$Cumulative_GORICA_weights), c("Hp", "Hn", "unconstrained"))
  expect_equal(res$Cumulative_GORICA_weights, ref$Cumulative_GORICA_weights)
  expect_equal(unname(res$LL_m[2, ]), unname(g2$result$loglik[c(2, 1, 3)]))
  # zonder koppeling op naam zou Hp van studie 1 met Hn van studie 2 zijn gecombineerd
  expect_true(final_w(res)[1] > 0.7)
})

test_that("evSyn_gorica: ongenoemde hypothesen met verschillende tekst geven een fout", {
  g3 <- goric(c(x = 2), VCOV = matrix(1), hypotheses = list("x > 0", "x < 0"), type = "gorica")
  g4 <- goric(c(x = 2), VCOV = matrix(1), hypotheses = list("x < 0", "x > 0"), type = "gorica")
  expect_error(evSyn(list(g3, g4)), "ambiguous")
  # identieke ongenoemde hypothesen zijn prima
  expect_silent(evSyn(list(g3, g3)))
})

test_that("evSyn_gorica: andere namen, andere comparison of andere penalty_factor geven een fout", {
  g1 <- goric(c(x = 2), VCOV = matrix(1), hypotheses = list(Hp = "x > 0", Hn = "x < 0"),
              type = "gorica")
  g5 <- goric(c(x = 2), VCOV = matrix(1), hypotheses = list(Ha = "x > 0", Hn = "x < 0"),
              type = "gorica")
  g6 <- goric(c(x = 2), VCOV = matrix(1), hypotheses = list(Hp = "x > 0"),
              type = "gorica", comparison = "complement")
  g7 <- goric(c(x = 2), VCOV = matrix(1), hypotheses = list(Hp = "x > 0", Hn = "x < 0"),
              type = "gorica", penalty_factor = 5)
  expect_error(evSyn(list(g1, g5)), "identical across the goric objects")
  expect_error(evSyn(list(g1, g6)), "same comparison")
  expect_error(evSyn(list(g1, g7)), "same 'penalty_factor'")
  # genoemde hypothesen met andere tekst (bv. andere parameternamen): waarschuwing
  g8 <- goric(c(y = 2), VCOV = matrix(1), hypotheses = list(Hp = "y > 0", Hn = "y < 0"),
              type = "gorica")
  expect_warning(evSyn(list(g1, g8)), "assumed to represent the same theory")
})

test_that("lijst-input: namen van de invoervectoren worden gecontroleerd", {
  expect_message(res <- evSyn(list(c(a = 0, b = 2), c(b = 0, a = 2)), input_type = "icvalues"),
                 "matched by name")
  expect_equal(unname(res$GORICA_m), rbind(c(0, 2), c(2, 0)))
  expect_error(evSyn(list(c(a = 0, b = 2), c(c = 0, a = 2)), input_type = "icvalues"),
               "identical across studies")
  expect_warning(evSyn(list(c(a = 0, b = 2), c(0, 2)), input_type = "icvalues"),
                 "matched across studies by position")
  expect_error(evSyn(list(c(a = .5, b = .5), c(c = .5, a = .5)), input_type = "icweights"),
               "identical across studies")
  expect_error(evSyn(list(c(a = 1, b = .5), c(c = 1, a = .5)), input_type = "icratios"),
               "identical across studies")
  expect_error(evSyn(list(c(a = -1, b = -2), c(c = -1, a = -2)), PT = list(c(1, 2), c(1, 2))),
               "identical across studies")
  # consistente namen: identiek aan ongenoemde invoer
  expect_equal(final_w(evSyn(lapply(W_fx, function(w) setNames(w, c("a", "b", "c"))))),
               final_w(evSyn(W_fx)))
})

# 11. Validatie van priorICweights en study_weights ------------------------------
test_that("priorICweights worden in elke route gevalideerd", {
  expect_error(evSyn(list(c(0, 2), c(0, 2)), input_type = "icvalues", priorICweights = c(-1, 2)),
               "non-negative")
  expect_error(evSyn(W_fx, priorICweights = c(-1, 1, 1)), "non-negative")
  expect_error(evSyn(R_fx, priorICweights = c(0, 0, 0)), "at least one positive")
  expect_error(evSyn(IC_fx, priorICweights = c(1, NA, 1)), "finite")
  expect_error(evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx, priorICweights = c(1, -1, 1)),
               "non-negative")
  expect_error(evSyn(g_fx, priorICweights = c(1, 2)), "should consist of 3 elements")
  LL <- lapply(g_fx, function(x) x$result$loglik)
  PT <- lapply(g_fx, function(x) x$result$penalty)
  expect_error(evSyn(LL, PT = PT, priorICweights = c(1, 2, Inf)), "finite")
  # prior gewicht 0 -> IC gewicht 0
  res <- evSyn(IC_fx, priorICweights = c(1, 1, 0))
  expect_equal(unname(final_w(res)[3]), 0)
  expect_equal(unname(res$GORICA_weight_m[, 3]), rep(0, 4))
})

test_that("study_weights worden in elke route gevalideerd", {
  expect_error(evSyn(IC_fx, study_weights = c(-1, 1, 1, 1)), "non-negative")
  expect_error(evSyn(W_fx, study_weights = c(1, 1, 1)), "should consist of 4 elements")
  expect_error(evSyn(R_fx, study_weights = rep(0, 4)), "at least one positive")
  expect_error(evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx, study_weights = c(1, NA, 1, 1)),
               "finite")
  expect_error(evSyn(g_fx, study_weights = c(1, 1, 1, Inf)), "finite")
})

test_that("studiegewicht 0: studie draagt niets bij, geen NaN", {
  # reproducer: gaf NaN
  res <- evSyn(list(c(0, 2), c(0, 2)), input_type = "icvalues", study_weights = c(0, 1))
  expect_false(anyNA(res$Cumulative_GORICA_weights))
  # studie 1 heeft gewicht 0: nog geen evidentie, dus prior (uniform) gewichten
  expect_equal(unname(res$Cumulative_GORICA_weights[1, ]), c(.5, .5))
  # finale resultaat = resultaat van studie 2 alleen (studie 1 telt niet mee)
  expect_equal(unname(res$Cumulative_GORICA["Final", ]), c(0, 2))
  expect_equal(unname(res$Cumulative_GORICA_weights["Final", ]),
               final_w(evSyn(list(c(0, 2)), input_type = "icvalues")))
  res_p <- evSyn(list(c(0, 2), c(0, 2)), input_type = "icvalues", study_weights = c(0, 1),
                 priorICweights = c(.2, .8))
  expect_equal(unname(res_p$Cumulative_GORICA_weights[1, ]), c(.2, .8))
  # laatste studie gewicht 0: finale resultaat = cumulatief t/m de vorige studie
  for (obj in list(IC_fx, W_fx, R_fx)) {
    a <- evSyn(obj, study_weights = c(1, 1, 1, 0))
    b <- evSyn(obj[1:3])
    expect_equal(unname(a$Cumulative_GORICA_weights["Final", ]),
                 unname(b$Cumulative_GORICA_weights["Final", ]))
    expect_equal(unname(a$Cumulative_GORICA_weights[3, ]),
                 unname(a$Cumulative_GORICA_weights[4, ]))
    expect_false(anyNA(a$Cumulative_GORICA_weights))
  }
  a <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx, study_weights = c(0, 1, 1, 1), type_ev = "equal")
  b <- evSyn(est_fx[2:4], VCOV = V_fx[2:4], hypotheses = H2_fx, type_ev = "equal")
  expect_equal(unname(a$Cumulative_GORICA_weights["Final", ]), final_w(b))
  expect_equal(unname(a$Cumulative_LL_weights["Final", ]), unname(b$Cumulative_LL_weights["Final", ]))
  expect_equal(unname(a$Cumulative_GORICA_weights[1, ]), rep(1/3, 3))
  expect_equal(unname(a$Cumulative_LL_weights[1, ]), rep(1/3, 3))
  expect_false(anyNA(a$Cumulative_GORICA))
  expect_output(print(a)); expect_output(print(summary(a)))
  # ook ICweights/ICratios met gewicht 0 voor de eerste studie
  w0 <- evSyn(W_fx, study_weights = c(0, 1, 1, 1), priorICweights = pr_fx)
  expect_equal(unname(w0$Cumulative_GORICA_weights[1, ]), pr_fx)
  expect_equal(unname(w0$Cumulative_GORICA_weights["Final", ]),
               final_w(evSyn(W_fx[2:4], priorICweights = pr_fx)))
  # leave1studyout met een gewicht 0
  l0 <- leave1studyout(evSyn(IC_fx, study_weights = c(0, 1, 1, 1)))
  expect_false(anyNA(l0$OverallGoricaWeights))
})

# 12. Underflow bij producten van IC weights ------------------------------------
test_that("veel studies: geen underflow bij IC weights en ratios", {
  res <- evSyn(rep(list(c(.5, .5)), 1100), input_type = "icweights")
  expect_equal(unname(final_w(res)), c(.5, .5))
  res_a <- evSyn(rep(list(c(.5, .5)), 1100), input_type = "icweights", type_ev = "average")
  expect_equal(unname(final_w(res_a)), c(.5, .5))
  res_r <- evSyn(rep(list(c(1, 1)), 1100), input_type = "icratios")
  expect_equal(unname(final_w(res_r)), c(.5, .5))
  res_o <- evSyn(rep(list(c(.4, .6)), 1100), input_type = "icweights", order_studies = "descending")
  expect_false(anyNA(res_o$Cumulative_GORICA_weights))
  expect_equal(unname(final_w(res_o)), c(0, 1))
  # normale invoer: identiek aan de directe berekening met producten
  res_n <- evSyn(W_fx, priorICweights = pr_fx, study_weights = sw_fx)
  W <- do.call(rbind, W_fx)
  CumW <- pr_fx * apply(W^(4 * sw_fx), 2, prod)
  expect_equal(unname(final_w(res_n)), unname(CumW / sum(CumW)), tolerance = 1e-10)
  res_n2 <- evSyn(W_fx, type_ev = "average", priorICweights = pr_fx)
  CumW2 <- pr_fx * apply(W^(1/4), 2, prod)
  expect_equal(unname(final_w(res_n2)), unname(CumW2 / sum(CumW2)), tolerance = 1e-10)
  # IC weight van 0 in een studie
  resz <- evSyn(list(c(0, 1), c(.5, .5)), input_type = "icweights")
  expect_equal(unname(final_w(resz)), c(0, 1))
})

# 13. Eén studie --------------------------------------------------------------
test_that("S = 1 werkt in elke route (cumulatief == de enkele studie)", {
  res <- evSyn(list(c(.7, .3)), input_type = "icweights")
  expect_equal(unname(res$Cumulative_GORICA_weights["Final", ]), c(.7, .3))
  expect_equal(unname(res$Cumulative_GORICA_weights[1, ]), c(.7, .3))
  res <- evSyn(list(c(.7, .3)), input_type = "icweights", type_ev = "average")
  expect_equal(unname(final_w(res)), c(.7, .3))
  res <- evSyn(list(c(1, .3)), input_type = "icratios")
  expect_equal(unname(final_w(res)), c(1, .3) / 1.3)
  res <- evSyn(list(c(0, 3)), input_type = "icvalues", order_studies = "ascending")
  expect_equal(unname(final_w(res)), unname(exp(-.5 * c(0, 3)) / sum(exp(-.5 * c(0, 3)))))
  expect_equal(unname(res$GORICA_weight_m[1, ]), unname(final_w(res)))
  res <- evSyn(list(c(-1, -2)), PT = list(c(1, 2)), order_studies = 1)
  expect_equal(unname(final_w(res)), unname(exp(-.5 * c(4, 8)) / sum(exp(-.5 * c(4, 8)))))
  res <- evSyn(est_fx[1], VCOV = V_fx[1], hypotheses = H2_fx, order_studies = "descending")
  expect_equal(unname(final_w(res)), unname(g_fx[[1]]$result$gorica.weights))
  expect_equal(unname(res$Cumulative_LL_weights["Final", ]), unname(g_fx[[1]]$result$loglik.weights))
  res <- evSyn(g_fx[1])
  expect_equal(unname(final_w(res)), unname(g_fx[[1]]$result$gorica.weights))
  # print, summary en plot
  for (r in list(evSyn(list(c(.7, .3)), input_type = "icweights"),
                 evSyn(list(c(1, .3)), input_type = "icratios"),
                 evSyn(IC_fx[1]), evSyn(g_fx[1]),
                 evSyn(est_fx[1], VCOV = V_fx[1], hypotheses = H2_fx))) {
    expect_output(print(r))
    expect_output(print(summary(r)))
    pdf(NULL); expect_no_error(print(plot(r))); dev.off()
  }
})

# 14. Voorkeurshypothese o.b.v. de finale (prior-gewogen) gewichten --------------
test_that("leave1studyout: OverallPrefHypo volgt de finale gewichten (incl. priorICweights)", {
  res <- evSyn(list(c(0, 2), c(0, 2), c(0, 2)), input_type = "icvalues",
               priorICweights = c(.01, .99))
  l1 <- leave1studyout(res)
  expect_equal(unname(l1$OverallPrefHypo[, 1]), rep("H2", 3))
  expect_true(all(l1$OverallGoricaWeights[, 2] > .9))
  expect_equal(unname(l1$OverallGorica[, 1]), rep(0, 3))
  # zonder prior: H1
  l0 <- leave1studyout(evSyn(list(c(0, 2), c(0, 2), c(0, 2)), input_type = "icvalues"))
  expect_equal(unname(l0$OverallPrefHypo[, 1]), rep("H1", 3))
})

test_that("order_studies ascending/descending: voorkeurshypothese o.b.v. finale gewichten", {
  # H1 heeft de laagste IC, maar met de prior wint H2
  IC <- list(c(0, 1, 5), c(0, 1, 5), c(0, 2, 5), c(0, .5, 5))
  pr <- c(.01, .99, 0)
  r <- evSyn(IC, priorICweights = pr, order_studies = "descending")
  w2 <- vapply(IC, function(x) exp(-.5 * x)[2] / sum(exp(-.5 * x)), numeric(1))
  expect_equal(r$order_studies, order(w2, decreasing = TRUE))
  r0 <- evSyn(IC, order_studies = "descending")
  w1 <- vapply(IC, function(x) exp(-.5 * x)[1] / sum(exp(-.5 * x)), numeric(1))
  expect_equal(r0$order_studies, order(w1, decreasing = TRUE))
  W <- lapply(IC, function(x) exp(-.5 * x) / sum(exp(-.5 * x)))
  expect_equal(evSyn(W, priorICweights = pr, order_studies = "descending")$order_studies,
               r$order_studies)
  R <- lapply(W, function(w) w / w[1])
  expect_equal(evSyn(R, priorICweights = pr, order_studies = "descending")$order_studies,
               r$order_studies)
})

test_that("leave1studyout werkt voor icweights- en icratios-objecten", {
  W <- list(c(.5, .5), c(.1, .9), c(.2, .8))
  l1 <- leave1studyout(evSyn(W, input_type = "icweights", priorICweights = c(.3, .7)))
  # weglaten van studie 1: product van de overige gewichten (x prior)
  expect_equal(unname(l1$OverallGoricaWeights[1, ]), c(.3 * .1 * .2, .7 * .9 * .8) / (.3 * .1 * .2 + .7 * .9 * .8))
  expect_true(isTRUE(l1$IC_is_diff))
  l2 <- leave1studyout(evSyn(list(c(1, 1), c(1, 9), c(1, 4)), input_type = "icratios"))
  expect_equal(unname(l2$OverallGoricaWeights[1, ]), c(1, 36) / 37)
  expect_equal(unname(l2$OverallPrefHypo[, 1]), rep("H2", 3))
  # consistent met evSyn zonder de weggelaten studie
  for (s in 1:4) {
    lw <- leave1studyout(evSyn(W_fx, priorICweights = pr_fx, study_weights = sw_fx))
    ref <- evSyn(W_fx[-s], priorICweights = pr_fx, study_weights = sw_fx[-s])
    expect_equal(unname(lw$OverallGoricaWeights[s, ]), final_w(ref))
    lr <- leave1studyout(evSyn(R_fx, priorICweights = pr_fx, study_weights = sw_fx))
    expect_equal(unname(lr$OverallGoricaWeights[s, ]), final_w(ref))
  }
})

# 15. penalty_factor ------------------------------------------------------------
test_that("evSyn_gorica neemt penalty_factor over van de goric-objecten", {
  gp <- lapply(1:4, function(i) goric(est_fx[[i]], VCOV = V_fx[[i]], hypotheses = H2_fx,
                                      type = "gorica", penalty_factor = 5))
  res <- evSyn(gp)
  expect_equal(res$penalty_factor, 5)
  expect_equal(unname(res$GORICA_m), unname(t(vapply(gp, function(x) x$result$gorica, numeric(3)))))
  expect_equal(unname(res$GORICA_weight_m[1, ]), unname(gp[[1]]$result$gorica.weights))
  expect_equal(unname(final_w(evSyn(gp[1]))), unname(gp[[1]]$result$gorica.weights))
  # LL-route met penalty_factor geeft hetzelfde
  LL <- lapply(gp, function(x) x$result$loglik)
  PT <- lapply(gp, function(x) x$result$penalty)
  for (te in c("added", "equal", "average")) {
    expect_equal(unname(evSyn(LL, PT = PT, penalty_factor = 5, type_ev = te)$Cumulative_GORICA_weights),
                 unname(evSyn(gp, type_ev = te)$Cumulative_GORICA_weights))
  }
  expect_false(isTRUE(all.equal(final_w(evSyn(LL, PT = PT)), final_w(res))))
  expect_error(evSyn(LL, PT = PT, penalty_factor = -1), "penalty_factor")
  # leave1studyout (type_ev = 'equal') gebruikt de penalty_factor
  l1 <- leave1studyout(evSyn(gp, type_ev = "equal"))
  ref <- evSyn(gp[-1], type_ev = "equal")
  expect_equal(unname(l1$OverallGoricaWeights[1, ]), final_w(ref))
  # est-route met penalty_factor via ...
  re <- evSyn(est_fx, VCOV = V_fx, hypotheses = H2_fx, penalty_factor = 5)
  expect_equal(re$penalty_factor, 5)
  expect_equal(unname(re$Cumulative_GORICA_weights), unname(res$Cumulative_GORICA_weights))
})
