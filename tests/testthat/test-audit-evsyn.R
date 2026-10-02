# Regressietests n.a.v. de audit (review4) van evSyn(): uitlijning van
# hypothesesets op naam (A6), ongelijke aantallen hypothesen (A7),
# auto-detectie van het invoertype (A8/B20), hypo_names als permutatie (B13),
# NA/NaN in de invoer (B14), priorWeights-shim (B15), "average" (B16),
# plot(ll_weights) voor IC-routes (B18), benoemde priorICweights/study_weights
# (B21), order_studies als studienamen (B22), type_ev = "equal" bij IC-routes (B23).

final_w <- function(x) unname(x$Cumulative_GORICA_weights["Final", ])
icw <- function(ic) { e <- exp(-(ic - min(ic)) / 2); e / sum(e) }

# A8 / B20: auto-detectie van het invoertype ------------------------------------
test_that("auto-detectie: IC-waarden met een 1 (V4) worden gemeld; expliciet icvalues geeft handberekening", {
  IC <- list(c(1, 3, 5), c(1, 2.5, 4))
  # handberekening (added): som IC = (2, 5.5, 9)
  hand <- icw(c(2, 5.5, 9))
  expect_equal(round(hand, 4), c(0.8306, 0.1443, 0.0251))
  res <- suppressMessages(evSyn(IC, input_type = "icvalues"))
  expect_s3_class(res, "evSyn_ICvalues")
  expect_equal(final_w(res), hand, tolerance = 1e-6)
  # auto-detectie: alle waarden > 0 en een gemeenschappelijke 1 -> icratios, met melding
  expect_message(r <- evSyn(IC), "treated as ratios of IC weights")
  expect_message(evSyn(IC), "input_type = 'icvalues'")
  expect_s3_class(r, "evSyn_ICratios")
})

test_that("auto-detectie: delta-IC met een 0 en een 1 is ambigu -> fout; met input_type correct", {
  dIC <- list(c(0, 1, 2.3), c(0, 1, .5))
  expect_error(evSyn(dIC), "ambiguous")
  expect_error(evSyn(dIC), "input_type")
  res <- suppressMessages(evSyn(dIC, input_type = "icvalues"))
  expect_equal(final_w(res), icw(c(0, 2, 2.8)))
  # beste hypothese (H1) krijgt het hoogste gewicht, niet 0
  expect_equal(which.max(final_w(res)), 1L)
})

test_that("auto-detectie: 1 op verschillende posities, of niet-positieve waarden, is ambigu", {
  expect_error(evSyn(list(c(1, 3, 5), c(2, 1, 4))), "not at the same position")
  expect_error(evSyn(list(c(1, 1, .5), c(.5, .2, 1), c(2, 1, 3))), "ambiguous")
  expect_error(evSyn(list(c(1, 0, 2), c(1, 2, 3))), "not all values are positive")
  # met input_type werkt het wel (ratio's zonder gemeenschappelijke referentie)
  expect_message(res <- evSyn(list(c(1, 1, .5), c(.5, .2, 1), c(2, 1, 3)), input_type = "icratios"),
                 "same reference hypothesis")
  expect_equal(unname(res$Href), 1)
})

test_that("auto-detectie: gewichten, ratio's en IC-waarden worden herkend en gemeld", {
  expect_message(r <- evSyn(list(c(.6, .3, .1), c(.2, .5, .3))), "treated as IC weights")
  expect_s3_class(r, "evSyn_ICweights")
  expect_message(r <- evSyn(list(c(1, .5, .2), c(1, .7, .1))), "treated as ratios of IC weights")
  expect_s3_class(r, "evSyn_ICratios")
  expect_equal(unname(r$Href), 1)
  expect_message(r <- evSyn(list(c(3, 5, 1), c(2.5, 4, 1))), "hypothesis 3 has a value of exactly 1")
  expect_message(r <- evSyn(list(c(0, 2, 4), c(1, 0, 3))), "treated as IC values")
  expect_s3_class(r, "evSyn_ICvalues")
  expect_message(evSyn(list(c(0, 2, 4), c(1, 0, 3))), "override")
  # afgeronde gewichten (B20): ambigu, geen stille icvalues-route
  expect_error(evSyn(list(c(.3333, .3333, .3333), c(.5, .25, .2499))), "approximately")
  expect_error(evSyn(list(c(.3333, .3333, .3333), c(.5, .25, .2499))), "input_type")
  # est- en LL-route: geen melding over het invoertype
  expect_silent(evSyn(list(c(-1, -2), c(-1, -2)), PT = list(c(1, 2), c(1, 2))))
})

# B14: NA/NaN in de invoer -------------------------------------------------------
test_that("NA/NaN in de invoer geeft een duidelijke restriktor-fout, Inf is toegestaan", {
  expect_error(evSyn(list(c(Inf, NaN, Inf), c(1, 1, 2))), "restriktor ERROR")
  expect_error(evSyn(list(c(Inf, NaN, Inf), c(1, 1, 2))), "NA or NaN")
  expect_error(evSyn(list(c(0, NA), c(1, 2)), input_type = "icvalues"), "NA or NaN")
  expect_error(evSyn(list(c(-1, -2), c(-1, -2)), PT = list(c(1, NA), c(1, 2))), "'PT' contains NA or NaN")
  expect_error(restriktor:::evSyn_ICweights(list(c(NaN, 1), c(.5, .5))), "NA or NaN")
  expect_error(restriktor:::evSyn_ICratios(list(c(1, NaN), c(1, .5))), "NA or NaN")
  # oneindige waarden blijven toegestaan (F-05: ratio Inf -> gewichten (1, 0))
  res <- suppressMessages(evSyn(list(c(Inf, 1), c(2, 1)), input_type = "icratios"))
  expect_equal(final_w(res), c(1, 0))
})

# A7: ongelijk aantal hypothesen per studie --------------------------------------
test_that("verschillend aantal hypothesen per studie geeft een fout (geen recycling)", {
  expect_error(evSyn(list(c(0, 2, 4), c(0, 2)), input_type = "icvalues"),
               "number of hypotheses must be identical across studies")
  expect_error(evSyn(list(c(.2, .3, .5), c(.5, .5)), input_type = "icweights"),
               "number of hypotheses must be identical across studies")
  expect_error(evSyn(list(c(1, 2, 4), c(1, 2)), input_type = "icratios"),
               "number of hypotheses must be identical across studies")
  expect_error(evSyn(list(c(-1, -2, -3), c(-1, -2)), PT = list(c(1, 2, 3), c(1, 2))),
               "number of hypotheses must be identical across studies")
  expect_error(evSyn(list(c(-1, -2, -3), c(-1, -2, -3)), PT = list(c(1, 2, 3), c(1, 2))),
               "'PT' have lengths")
  expect_error(evSyn(list(c(-1, -2, -3), c(-1, -2, -3)), PT = list(c(1, 2, 3))),
               "number of elements in 'PT'")
  expect_error(evSyn(list(c(-1, -2, -3), c(-1, -2, -3)), PT = list(c(1, 2, 3), c(1, 2, 3), c(1, 2, 3))),
               "number of elements in 'PT'")
  expect_error(restriktor:::evSyn_ICvalues(list(c(0, 2, 4), c(0, 2))),
               "number of hypotheses must be identical across studies")
})

# A6: est-route, hypothesesets per studie op naam uitgelijnd ---------------------
test_that("evSyn_est: hypothesesets met dezelfde namen in andere volgorde worden op naam gekoppeld", {
  V <- diag(2) * .3
  est <- list(c(x = 2, y = 1), c(x = 2, y = 1))
  sets <- list(list(Hp = "x > y", Hn = "x < y"), list(Hn = "x < y", Hp = "x > y"))
  ref <- evSyn(est, VCOV = list(V, V), hypotheses = list(Hp = "x > y", Hn = "x < y"))
  expect_message(r1 <- evSyn(est, VCOV = list(V, V), hypotheses = sets, hypo_names = c("Hp", "Hn")),
                 "matched by name")
  expect_message(r2 <- evSyn(est, VCOV = list(V, V), hypotheses = sets), "matched by name")
  for (r in list(r1, r2)) {
    expect_equal(colnames(r$LL_m), c("Hp", "Hn", "Unconstrained"))
    expect_equal(unname(r$LL_m[1, ]), unname(r$LL_m[2, ]))
    expect_equal(r$Cumulative_GORICA_weights, ref$Cumulative_GORICA_weights)
    expect_equal(round(final_w(r), 3), c(0.642, 0.121, 0.236))
  }
  # drie studies, alleen set 3 gepermuteerd
  sets3 <- list(list(Hp = "x > y", Hn = "x < y"), list(Hp = "x > y", Hn = "x < y"), list(Hn = "x < y", Hp = "x > y"))
  r3 <- suppressMessages(evSyn(est[c(1, 2, 1)], VCOV = list(V, V, V), hypotheses = sets3))
  expect_equal(unname(r3$LL_m[3, ]), unname(r3$LL_m[1, ]))
  # verschillende namen -> fout
  expect_error(evSyn(est, VCOV = list(V, V), hypotheses = list(list(Ha = "x > y"), list(Hb = "x > y"))),
               "identical across the hypothesis sets")
  expect_error(evSyn(est, VCOV = list(V, V), hypotheses = list(list(Hp = "x > y", Hn = "x < y"), list(Hp = "x > y", Hz = "x < y"))),
               "identical across the hypothesis sets")
  # deels benoemd -> waarschuwing en positie
  expect_warning(r4 <- evSyn(est, VCOV = list(V, V), hypotheses = list(list(Hp = "x > y", Hn = "x < y"), list("x > y", "x < y"))),
                 "matched across studies by position")
  expect_equal(colnames(r4$LL_m), c("H1", "H2", "Unconstrained"))
  # onbenoemde sets: op positie, zonder melding (labels H1, H2)
  expect_silent(r5 <- evSyn(est, VCOV = list(V, V), hypotheses = list(list("x > y", "x < y"), list("x > y", "x < y"))))
  expect_equal(colnames(r5$LL_m), c("H1", "H2", "Unconstrained"))
  expect_equal(unname(r5$Cumulative_GORICA_weights), unname(ref$Cumulative_GORICA_weights))
})

# B13: hypo_names als permutatie van de invoernamen ------------------------------
test_that("hypo_names die een permutatie van de invoernamen is, geeft een waarschuwing (labels positioneel)", {
  LL <- list(c(a = -1, b = -5, c = -3), c(a = -2, b = -4, c = -3))
  PT <- list(c(a = 1, b = 2, c = 3), c(a = 1, b = 2, c = 3))
  expect_warning(r <- evSyn(LL, PT = PT, hypo_names = c("c", "b", "a")), "in another order")
  expect_warning(evSyn(LL, PT = PT, hypo_names = c("c", "b", "a")), "NOT re-ordered")
  # labels positioneel (gedocumenteerd): kolom 'c' bevat de waarde van 'a'
  expect_equal(unname(r$LL_m[1, "c"]), -1)
  # niet-permutatie waarschuwt nog steeds; gelijke volgorde niet
  expect_warning(evSyn(LL, PT = PT, hypo_names = c("x", "y", "z")), "differ from the hypothesis names")
  expect_silent(evSyn(LL, PT = PT, hypo_names = c("a", "b", "c")))
  # icweights en icvalues
  expect_warning(evSyn(list(c(a = .7, b = .2, c = .1), c(a = .5, b = .3, c = .2)), input_type = "icweights",
                       hypo_names = c("c", "b", "a")), "in another order")
  expect_warning(evSyn(list(c(a = 0, b = 2, c = 4), c(a = 0, b = 2, c = 4)), input_type = "icvalues",
                       hypo_names = c("c", "b", "a")), "in another order")
  # goric-objecten
  g1 <- goric(c(x = 2), VCOV = matrix(1), hypotheses = list(Hp = "x > 0", Hn = "x < 0"), type = "gorica")
  expect_warning(evSyn(list(g1, g1), hypo_names = c("Hn", "Hp", "unconstrained")), "in another order")
  # herlabelen met andere namen is prima (geen waarschuwing)
  expect_silent(evSyn(list(g1, g1), hypo_names = c("H1.1", "H1.2", "Hu")))
  # est-route
  V <- diag(2) * .3
  expect_warning(evSyn(list(c(x = 2, y = 1), c(x = 2, y = 1)), VCOV = list(V, V),
                       hypotheses = list(Hp = "x > y", Hn = "x < y"), hypo_names = c("Hn", "Hp")),
                 "in another order")
})

# B15: priorWeights-shim in de directe aanroepen ---------------------------------
test_that("priorWeights werkt (met waarschuwing) in evSyn_LL, evSyn_ICvalues, evSyn_gorica en evSyn_est", {
  IC <- list(c(0, 2, 4), c(0, 2, 4)); p <- c(.1, .1, .8)
  expect_warning(r1 <- restriktor:::evSyn_ICvalues(IC, priorWeights = p), "deprecated")
  expect_equal(r1$priorICweights, p)
  expect_equal(r1, restriktor:::evSyn_ICvalues(IC, priorICweights = p))
  LL <- list(c(-1, -2, -3), c(-1, -2, -3)); PT <- list(c(1, 2, 3), c(1, 2, 3))
  expect_warning(r2 <- restriktor:::evSyn_LL(LL, PT = PT, priorWeights = p), "deprecated")
  expect_equal(r2, restriktor:::evSyn_LL(LL, PT = PT, priorICweights = p))
  g <- lapply(1:2, function(i) goric(c(x = 2, y = 1), VCOV = diag(2),
                                     hypotheses = list(Hp = "x > y", Hn = "x < y"), type = "gorica"))
  expect_warning(r3 <- restriktor:::evSyn_gorica(g, priorWeights = p), "deprecated")
  expect_equal(r3, restriktor:::evSyn_gorica(g, priorICweights = p))
  est <- list(c(x = 2, y = 1), c(x = 2, y = 1)); V <- list(diag(2), diag(2))
  H <- list(Hp = "x > y", Hn = "x < y")
  expect_warning(r4 <- restriktor:::evSyn_est(est, VCOV = V, hypotheses = H, priorWeights = p), "deprecated")
  expect_equal(r4, restriktor:::evSyn_est(est, VCOV = V, hypotheses = H, priorICweights = p))
  expect_equal(r4$priorICweights, p)
  expect_error(restriktor:::evSyn_est(est, VCOV = V, hypotheses = H, priorWeights = p, priorICweights = p),
               "priorWeights")
  # overige ...-argumenten gaan nog steeds naar goric()
  expect_equal(suppressWarnings(evSyn(est, VCOV = V, hypotheses = H, priorWeights = p, penalty_factor = 5))$penalty_factor, 5)
})

# B16: type_ev = "average" -------------------------------------------------------
test_that("average: cumulatieve LL is een gemiddelde, Final-tabel reconstrueerbaar, LL-gewichten consistent", {
  LL <- list(c(-1, -5, -3), c(-2, -4, -3), c(-1.5, -4.5, -3))
  PT <- list(c(1, 2, 3), c(1, 2, 3), c(1, 2, 3))
  LLm <- do.call(rbind, LL)
  for (te in c("added", "equal", "average")) {
    r <- evSyn(LL, PT = PT, type_ev = te)
    ft <- summary(r)$Final_Cumulative_results
    expect_equal(unname(ft["GORICA values", ]),
                 unname(-2 * ft["Log-likelihood values", ] + 2 * ft["Penalty term values", ]))
  }
  r <- evSyn(LL, PT = PT, type_ev = "average")
  expect_equal(unname(r$Cumulative_LL["Study nr.s 1-3   ", ]), colMeans(LLm))
  expect_equal(unname(r$Cumulative_LL[2, ]), colMeans(LLm[1:2, ]))
  m <- colMeans(LLm)
  expect_equal(unname(r$Cumulative_LL_weights["Final", ]), unname(exp(m - max(m)) / sum(exp(m - max(m)))))
  expect_equal(unname(summary(r)$Final_Cumulative_results["Log-likelihood values", ]), unname(colMeans(LLm)))
  # added: som (ongewijzigd)
  ra <- evSyn(LL, PT = PT, type_ev = "added")
  expect_equal(unname(ra$Cumulative_LL[3, ]), colSums(LLm))
  # studiegewicht 0 telt niet mee in het gemiddelde
  r0 <- evSyn(LL, PT = PT, type_ev = "average", study_weights = c(1, 0, 1))
  expect_equal(unname(r0$Cumulative_LL[3, ]), colMeans(LLm[c(1, 3), ]))
  expect_equal(unname(r0$Cumulative_LL[2, ]), unname(r0$Cumulative_LL[1, ]))
  # est-route en goric-route identiek aan LL-route
  est <- list(c(x1 = 5, x2 = 3, x3 = 1), c(x1 = 4.5, x2 = 2.5, x3 = 3))
  V <- list(diag(3) * .1, diag(3) * .15)
  H <- list(H1 = "x1 > x2 > x3", H2 = "x1 > x2")
  g <- lapply(1:2, function(i) goric(est[[i]], VCOV = V[[i]], hypotheses = H, type = "gorica"))
  re <- evSyn(est, VCOV = V, hypotheses = H, type_ev = "average")
  rg <- evSyn(g, type_ev = "average")
  rl <- evSyn(lapply(g, function(x) x$result$loglik), PT = lapply(g, function(x) x$result$penalty), type_ev = "average")
  expect_equal(unname(re$Cumulative_LL), unname(rl$Cumulative_LL))
  expect_equal(unname(rg$Cumulative_LL_weights), unname(rl$Cumulative_LL_weights))
  expect_equal(unname(re$Cumulative_LL[2, ]), colMeans(rbind(g[[1]]$result$loglik, g[[2]]$result$loglik)))
})

# B18: plot(output_type = "ll_weights") voor IC-routes --------------------------
test_that("plot(ll_weights) geeft een duidelijke fout voor icvalues/icweights/icratios", {
  pdf(NULL); on.exit(dev.off())
  for (it in c("icvalues", "icweights", "icratios")) {
    obj <- switch(it, icvalues = list(c(0, 2, 4), c(1, 0, 3)),
                  icweights = list(c(.5, .3, .2), c(.2, .3, .5)),
                  icratios = list(c(1, .6, .4), c(1, 2, 3)))
    r <- suppressMessages(evSyn(obj, input_type = it))
    expect_error(plot(r, output_type = "ll_weights"), "not available")
    expect_error(plot(r, output_type = "ll_weights"), "restriktor ERROR")
    expect_no_error(print(plot(r)))
  }
  r <- evSyn(list(c(-1, -2), c(-1, -2)), PT = list(c(1, 2), c(1, 2)))
  expect_no_error(print(plot(r, output_type = "ll_weights")))
})

# B21: benoemde priorICweights en study_weights ----------------------------------
test_that("benoemde priorICweights worden op naam gekoppeld; niet-passende namen geven een fout", {
  LL <- list(c(H1 = -1, H2 = -5, H3 = -3), c(H1 = -2, H2 = -4, H3 = -3))
  PT <- list(c(H1 = 1, H2 = 2, H3 = 3), c(H1 = 1, H2 = 2, H3 = 3))
  r0 <- evSyn(LL, PT = PT, priorICweights = c(.1, .1, .8))
  r1 <- evSyn(LL, PT = PT, priorICweights = c(H1 = .1, H2 = .1, H3 = .8))
  r2 <- evSyn(LL, PT = PT, priorICweights = c(H3 = .8, H1 = .1, H2 = .1))
  expect_equal(r1, r0)
  expect_equal(r2, r0)
  expect_equal(r2$priorICweights, c(.1, .1, .8))
  expect_error(evSyn(LL, PT = PT, priorICweights = c(a = .8, b = .1, c = .1)), "names of 'priorICweights'")
  expect_error(evSyn(LL, PT = PT, priorICweights = c(H1 = .8, H2 = .1, H2 = .1)), "names of 'priorICweights'")
  # IC-routes: kolomnamen H1, H2, H3 (of hypo_names)
  IC <- list(c(0, 2, 4), c(1, 0, 3))
  expect_equal(suppressMessages(evSyn(IC, priorICweights = c(H3 = .8, H1 = .1, H2 = .1))),
               suppressMessages(evSyn(IC, priorICweights = c(.1, .1, .8))))
  expect_equal(suppressMessages(evSyn(IC, priorICweights = c(c = .8, a = .1, b = .1), hypo_names = c("a", "b", "c"))),
               suppressMessages(evSyn(IC, priorICweights = c(.1, .1, .8), hypo_names = c("a", "b", "c"))))
  W <- list(c(.6, .3, .1), c(.2, .5, .3))
  expect_equal(evSyn(W, input_type = "icweights", priorICweights = c(H2 = .3, H3 = .6, H1 = .1)),
               evSyn(W, input_type = "icweights", priorICweights = c(.1, .3, .6)))
  R <- lapply(W, function(w) w / w[3])
  expect_equal(evSyn(R, input_type = "icratios", priorICweights = c(H2 = .3, H3 = .6, H1 = .1)),
               evSyn(R, input_type = "icratios", priorICweights = c(.1, .3, .6)))
  # est-route: incl. de naam van de failsafe-hypothese
  est <- list(c(x = 2, y = 1), c(x = 2, y = 1)); V <- list(diag(2) * .3, diag(2) * .3)
  H <- list(Hp = "x > y", Hn = "x < y")
  expect_equal(evSyn(est, VCOV = V, hypotheses = H, priorICweights = c(Unconstrained = .5, Hn = .3, Hp = .2)),
               evSyn(est, VCOV = V, hypotheses = H, priorICweights = c(.2, .3, .5)))
  # het failsafe-label mag ook in kleine letters (zoals in goric()) worden gegeven
  expect_equal(evSyn(est, VCOV = V, hypotheses = H, priorICweights = c(unconstrained = .5, Hn = .3, Hp = .2)),
               evSyn(est, VCOV = V, hypotheses = H, priorICweights = c(.2, .3, .5)))
  expect_error(evSyn(est, VCOV = V, hypotheses = H, priorICweights = c(unconstr = .5, Hn = .3, Hp = .2)),
               "names of 'priorICweights'")
  # goric-route
  g <- lapply(1:2, function(i) goric(est[[i]], VCOV = V[[i]], hypotheses = H, type = "gorica"))
  expect_equal(evSyn(g, priorICweights = c(unconstrained = .5, Hn = .3, Hp = .2)),
               evSyn(g, priorICweights = c(.2, .3, .5)))
})

test_that("benoemde study_weights worden op naam gekoppeld aan study_names", {
  LL <- list(c(-1, -5, -3), c(-2, -4, -3)); PT <- list(c(1, 2, 3), c(1, 2, 3))
  r0 <- evSyn(LL, PT = PT, study_names = c("A", "B"), study_weights = c(.1, .9))
  r1 <- evSyn(LL, PT = PT, study_names = c("A", "B"), study_weights = c(B = .9, A = .1))
  expect_equal(r1, r0)
  expect_equal(r1$study_weights, c(.1, .9))
  # zonder study_names: namen "1", "2", ...
  r2 <- evSyn(LL, PT = PT, study_weights = c("2" = .9, "1" = .1))
  expect_equal(r2$study_weights, c(.1, .9))
  expect_error(evSyn(LL, PT = PT, study_names = c("A", "B"), study_weights = c(X = .9, Y = .1)),
               "names of 'study_weights'")
  expect_error(evSyn(LL, PT = PT, study_weights = c(A = .9, B = .1)), "names of 'study_weights'")
  # IC- en est-route
  IC <- list(c(0, 2, 4), c(1, 0, 3))
  expect_equal(suppressMessages(evSyn(IC, study_names = c("A", "B"), study_weights = c(B = .9, A = .1))),
               suppressMessages(evSyn(IC, study_names = c("A", "B"), study_weights = c(.1, .9))))
  est <- list(c(x = 2, y = 1), c(x = 3, y = 1)); V <- list(diag(2) * .3, diag(2) * .3)
  H <- list(Hp = "x > y", Hn = "x < y")
  expect_equal(evSyn(est, VCOV = V, hypotheses = H, study_names = c("A", "B"), study_weights = c(B = .9, A = .1)),
               evSyn(est, VCOV = V, hypotheses = H, study_names = c("A", "B"), study_weights = c(.1, .9)))
  # in combinatie met herordenen: gewichten volgen de studies
  r3 <- evSyn(LL, PT = PT, study_names = c("A", "B"), study_weights = c(B = .9, A = .1), order_studies = c(2, 1))
  expect_equal(r3$study_weights, c(.9, .1))
  expect_equal(r3$study_names, c("B", "A"))
})

# B22: order_studies als studienamen / ongeldige invoer --------------------------
test_that("order_studies: studienamen worden vertaald naar posities; lijst/logical geeft een duidelijke fout", {
  IC <- list(c(0, 2, 4), c(1, 0, 3), c(2, 2, 0))
  r <- suppressMessages(evSyn(IC, input_type = "icvalues", study_names = c("A", "B", "C"), order_studies = c("C", "A", "B")))
  expect_equal(r$order_studies, c(3L, 1L, 2L))
  expect_equal(r$study_names, c("C", "A", "B"))
  expect_equal(r, suppressMessages(evSyn(IC, input_type = "icvalues", study_names = c("A", "B", "C"), order_studies = c(3, 1, 2))))
  for (os in list(list(1, 2, 3), TRUE, c("C", "A", "B"), c("A", "B", "X"), c("A", "B"), 1i)) {
    expect_error(evSyn(IC, input_type = "icvalues", order_studies = os), "restriktor ERROR")
    expect_error(evSyn(IC, input_type = "icvalues", order_studies = os), "order_studies")
  }
  expect_error(evSyn(IC, input_type = "icvalues", study_names = c("A", "B", "C"), order_studies = c("A", "B", "X")),
               "order_studies")
  expect_error(evSyn(IC, input_type = "icvalues", order_studies = "nonsense"), "order_studies")
  # numerieke validatie ongewijzigd
  expect_error(evSyn(IC, input_type = "icvalues", order_studies = c(1, 1, 2)), "permutation of 1:3")
  expect_error(evSyn(IC, input_type = "icvalues", order_studies = c(1, 2)), "must equal the number of studies")
  # alle routes
  LL <- list(c(-1, -5, -3), c(-2, -4, -3), c(-1.5, -4.5, -3)); PT <- rep(list(c(1, 2, 3)), 3)
  expect_equal(evSyn(LL, PT = PT, study_names = c("A", "B", "C"), order_studies = c("B", "C", "A"))$order_studies, c(2L, 3L, 1L))
  W <- lapply(IC, function(x) exp(-x / 2) / sum(exp(-x / 2)))
  expect_equal(evSyn(W, input_type = "icweights", study_names = c("A", "B", "C"), order_studies = c("B", "C", "A"))$order_studies, c(2L, 3L, 1L))
  R <- lapply(W, function(w) w / w[1])
  expect_equal(evSyn(R, input_type = "icratios", study_names = c("A", "B", "C"), order_studies = c("B", "C", "A"))$order_studies, c(2L, 3L, 1L))
  est <- list(c(x = 2, y = 1), c(x = 3, y = 1), c(x = 1, y = 2)); V <- rep(list(diag(2) * .3), 3)
  H <- list(Hp = "x > y", Hn = "x < y")
  re <- evSyn(est, VCOV = V, hypotheses = H, study_names = c("A", "B", "C"), order_studies = c("B", "C", "A"))
  expect_equal(re$order_studies, c(2L, 3L, 1L))
  expect_equal(re$LL_m, evSyn(est, VCOV = V, hypotheses = H, study_names = c("A", "B", "C"), order_studies = c(2, 3, 1))$LL_m)
  g <- lapply(1:3, function(i) goric(est[[i]], VCOV = V[[i]], hypotheses = H, type = "gorica"))
  expect_equal(evSyn(g, study_names = c("A", "B", "C"), order_studies = c("B", "C", "A"))$order_studies, c(2L, 3L, 1L))
  # gedeeltelijke string blijft werken
  expect_equal(suppressMessages(evSyn(IC, input_type = "icvalues", order_studies = "desc"))$order_studies,
               suppressMessages(evSyn(IC, input_type = "icvalues", order_studies = "descending"))$order_studies)
})

# B23: expliciet input_type + type_ev = "equal" op IC-routes --------------------
test_that("type_ev = 'equal' valt bij IC-routes terug op 'added' (met melding), ook met expliciet input_type", {
  W <- list(c(.6, .3, .1), c(.2, .5, .3))
  expect_message(a <- evSyn(W, input_type = "icweights", type_ev = "equal"), "added-evidence approach is used instead")
  expect_equal(a$type_ev, "added")
  expect_equal(a, evSyn(W, input_type = "icweights", type_ev = "added"))
  expect_message(b <- evSyn(list(c(0, 2, 4), c(1, 0, 3)), input_type = "icvalues", type_ev = "equal"),
                 "added-evidence approach is used instead")
  expect_equal(b$type_ev, "added")
  expect_message(d <- evSyn(list(c(1, .5, .2), c(1, 2, 3)), input_type = "icratios", type_ev = "equal"),
                 "added-evidence approach is used instead")
  expect_equal(d$type_ev, "added")
  expect_message(restriktor:::evSyn_ICvalues(list(c(0, 2, 4), c(1, 0, 3)), type_ev = "equal"), "added-evidence")
  # ook via auto-detectie
  expect_message(e <- evSyn(W, type_ev = "equal"), "added-evidence approach is used instead")
  expect_equal(e$type_ev, "added")
  # andere waarden blijven ongeldig
  expect_error(evSyn(W, input_type = "icweights", type_ev = "nonsense"))
})

# W2-N-01: Heq = TRUE in de est-route ---------------------------------------------
test_that("evSyn_est met Heq = TRUE: kolommen Heq, hypothese, Complement; gelijk aan goric() per studie", {
  V <- diag(2) * .3; e1 <- c(x = 2, y = 1); e2 <- c(x = 1, y = 2)
  g1 <- goric(e1, VCOV = V, hypotheses = list(Hp = "x > y"), comparison = "complement", Heq = TRUE, type = "gorica")
  g2 <- goric(e2, VCOV = V, hypotheses = list(Hp = "x > y"), comparison = "complement", Heq = TRUE, type = "gorica")
  r <- evSyn(list(e1, e2), VCOV = list(V, V), hypotheses = list(Hp = "x > y"), Heq = TRUE)
  expect_equal(colnames(r$LL_m), c("Heq", "Hp", "Complement"))
  expect_equal(unname(r$LL_m[1, ]), g1$result$loglik)
  expect_equal(unname(r$LL_m[2, ]), g2$result$loglik)
  expect_equal(unname(r$PT_m[1, ]), g1$result$penalty)
  expect_equal(unname(r$ratio_GORICA_weight_mu[, 1]), c(g1$ratio.gw[2, 3], g2$ratio.gw[2, 3]))
  expect_equal(colnames(r$ratio_GORICA_weight_mu), "Hp vs. Complement")
  # handberekening (added): som van de IC-waarden
  ic <- colSums(-2 * r$LL_m + 2 * r$PT_m)
  expect_equal(final_w(r), unname(icw(ic)))
  # gelijk aan de goric-objectroute
  rg <- evSyn(list(g1, g2))
  expect_equal(unname(r$Cumulative_GORICA_weights), unname(rg$Cumulative_GORICA_weights))
  # een hypothese 'Heq' in de set wordt (zoals in goric) verwijderd en opnieuw gegenereerd
  r2 <- evSyn(list(e1, e2), VCOV = list(V, V), hypotheses = list(Heq = "x = y", Hp = "x > y"),
              comparison = "complement", Heq = TRUE)
  expect_equal(r2$LL_m, r$LL_m)
  # hypo_names (1 label) en benoemde priorICweights (incl. Heq en Complement)
  r3 <- evSyn(list(e1, e2), VCOV = list(V, V), hypotheses = list(Hp = "x > y"), Heq = TRUE,
              hypo_names = "Pos", priorICweights = c(Complement = .5, Hp = .3, Heq = .2))
  expect_equal(colnames(r3$LL_m), c("Heq", "Pos", "Complement"))
  expect_equal(r3$priorICweights, c(.2, .3, .5))
  # leave1studyout werkt met de extra kolom
  expect_s3_class(leave1studyout(r), "leave1studyout.evSyn")
  # meer dan een hypothese: Heq wordt (met een waarschuwing) genegeerd
  expect_warning(r4 <- evSyn(list(e1, e2), VCOV = list(V, V), hypotheses = list(Hp = "x > y", Hn = "x < y"), Heq = TRUE),
                 "'Heq' argument is ignored")
  expect_equal(colnames(r4$LL_m), c("Hp", "Hn", "Unconstrained"))
  # hypothese zonder ongelijkheden: goric negeert Heq -> duidelijke fout (geen subscript out of bounds)
  expect_error(suppressMessages(evSyn(list(e1, e2), VCOV = list(V, V), hypotheses = list(Hp = "x = y"), Heq = TRUE)),
               "'Heq' option cannot be used")
})

# W2-N-02: numerieke study_names + benoemde study_weights -------------------------
test_that("benoemde study_weights met numerieke study_names (jaartallen) worden op naam gekoppeld", {
  LL <- list(c(-1, -5, -3), c(-2, -4, -3), c(-1.5, -4.5, -3)); PT <- rep(list(c(1, 2, 3)), 3)
  r0 <- evSyn(LL, PT = PT, study_names = c(2001, 2005, 2010), study_weights = c(.2, .3, .5))
  r1 <- evSyn(LL, PT = PT, study_names = c(2001, 2005, 2010), study_weights = c("2010" = .5, "2001" = .2, "2005" = .3))
  expect_equal(r1$study_weights, c(.2, .3, .5))
  expect_equal(r1, r0)
  expect_equal(rownames(r1$LL_m), c("2001", "2005", "2010"))
  expect_error(evSyn(LL, PT = PT, study_names = c(2001, 2005, 2010), study_weights = c("2011" = .5, "2001" = .2, "2005" = .3)),
               "names of 'study_weights'")
  # ook order_studies als (character) studienamen bij numerieke study_names
  expect_equal(evSyn(LL, PT = PT, study_names = c(2001, 2005, 2010), order_studies = c("2010", "2001", "2005"))$order_studies,
               c(3L, 1L, 2L))
})

# W2-N-03: penalty_factor bij een lijst van goric-objecten --------------------------
test_that("evSyn(<goric-lijst>, penalty_factor = ...) geeft een duidelijke fout bij een afwijkende waarde", {
  V <- diag(2) * .3; H <- list(Hp = "x > y", Hn = "x < y")
  g <- lapply(list(c(x = 2, y = 1), c(x = 1, y = 2)), function(e) goric(e, VCOV = V, hypotheses = H, type = "gorica"))
  expect_error(evSyn(g, penalty_factor = 3), "taken from the goric objects")
  expect_error(restriktor:::evSyn_gorica(g, penalty_factor = 3), "taken from the goric objects")
  # dezelfde waarde als in de objecten is toegestaan
  expect_equal(evSyn(g, penalty_factor = 2), evSyn(g))
  expect_equal(evSyn(g)$penalty_factor, 2)
})

# W2-N-04: benoemde invoervectoren + priorICweights met de invoernamen ------------------
test_that("benoemde priorICweights mogen de namen van de invoervectoren gebruiken (ook met hypo_names)", {
  LL <- list(c(a = -1, b = -5, c = -3), c(a = -2, b = -4, c = -3)); PT <- rep(list(c(a = 1, b = 2, c = 3)), 2)
  r0 <- evSyn(LL, PT = PT, priorICweights = c(.1, .1, .8))
  r1 <- evSyn(LL, PT = PT, priorICweights = c(c = .8, a = .1, b = .1))
  expect_equal(colnames(r1$LL_m), c("H1", "H2", "H3"))
  expect_equal(r1, r0)
  # met hypo_names: labels of invoernamen, niet gemengd
  # (de waarschuwing dat hypo_names afwijkt van de invoernamen is hier verwacht)
  r2 <- suppressWarnings(evSyn(LL, PT = PT, hypo_names = c("A", "B", "C"), priorICweights = c(c = .8, a = .1, b = .1)))
  r3 <- suppressWarnings(evSyn(LL, PT = PT, hypo_names = c("A", "B", "C"), priorICweights = c(C = .8, A = .1, B = .1)))
  expect_equal(colnames(r2$LL_m), c("A", "B", "C"))
  expect_equal(r2$priorICweights, c(.1, .1, .8))
  expect_equal(r3, r2)
  expect_error(suppressWarnings(evSyn(LL, PT = PT, hypo_names = c("A", "B", "C"), priorICweights = c(C = .8, a = .1, b = .1))),
               "names of the input")
  # IC-routes
  IC <- list(c(a = 1, b = 5, c = 3), c(a = 2, b = 4, c = 3))
  expect_equal(evSyn(IC, input_type = "icvalues", priorICweights = c(c = .8, a = .1, b = .1))$priorICweights, c(.1, .1, .8))
  W <- list(c(a = .2, b = .5, c = .3), c(a = .1, b = .4, c = .5))
  expect_equal(suppressWarnings(evSyn(W, input_type = "icweights", hypo_names = c("A", "B", "C"),
                                      priorICweights = c(c = .8, a = .1, b = .1)))$priorICweights, c(.1, .1, .8))
  R <- list(c(a = 1, b = 2, c = 3), c(a = 1, b = 4, c = 5))
  expect_equal(evSyn(R, input_type = "icratios", priorICweights = c(c = .8, a = .1, b = .1))$priorICweights, c(.1, .1, .8))
  # est-route (namen van de hypothesen + failsafe) en goric-route (modelnamen) met hypo_names
  V <- diag(2) * .3; est <- list(c(x = 2, y = 1), c(x = 1, y = 2)); H <- list(Hp = "x > y", Hn = "x < y")
  re <- evSyn(est, VCOV = list(V, V), hypotheses = H, hypo_names = c("P", "N"),
              priorICweights = c(Unconstrained = .5, Hn = .3, Hp = .2))
  expect_equal(re$priorICweights, c(.2, .3, .5))
  expect_equal(colnames(re$LL_m), c("P", "N", "Unconstrained"))
  g <- lapply(est, function(e) goric(e, VCOV = V, hypotheses = H, type = "gorica"))
  rg <- evSyn(g, hypo_names = c("P", "N", "U"), priorICweights = c(unconstrained = .5, Hn = .3, Hp = .2))
  expect_equal(rg$priorICweights, c(.2, .3, .5))
  expect_equal(rg, evSyn(g, hypo_names = c("P", "N", "U"), priorICweights = c(U = .5, N = .3, P = .2)))
  expect_error(evSyn(g, priorICweights = c(X = .5, Hn = .3, Hp = .2)), "names of 'priorICweights'")
})

# Klein: exacte studienamen in order_studies; type_ev = "eq" bij IC-routes ---------
test_that("order_studies: studienaam die een afkorting van een keuze is, wordt exact als studienaam genomen; 'eq' = 'equal'", {
  expect_equal(restriktor:::.evSyn_order_studies("a", 1, "a"), 1L)
  expect_equal(restriktor:::.evSyn_order_studies("asc", 1, "a"), "ascending")
  expect_equal(restriktor:::.evSyn_order_studies(c("in", "desc", "asc"), 3, c("asc", "desc", "in")), c(3L, 2L, 1L))
  W <- list(c(.6, .3, .1), c(.2, .5, .3))
  expect_message(a <- evSyn(W, input_type = "icweights", type_ev = "eq"), "added-evidence approach is used instead")
  expect_equal(a$type_ev, "added")
  expect_equal(suppressMessages(evSyn(list(c(0, 2, 4), c(1, 0, 3)), input_type = "icvalues", type_ev = "av"))$type_ev, "average")
})

test_that("W4-N-03: prior-naam 'complement' (kleine letters) wordt geaccepteerd in de est-route met Heq", {
  e1 <- c(x = 2, y = 1); e2 <- c(x = 2.5, y = 1); V <- diag(2) * .3
  base <- suppressMessages(evSyn(list(e1, e2), VCOV = list(V, V),
                                 hypotheses = list(Hp = "x > y"), Heq = TRUE,
                                 priorICweights = c(.2, .5, .3)))
  low <- suppressMessages(evSyn(list(e1, e2), VCOV = list(V, V),
                                hypotheses = list(Hp = "x > y"), Heq = TRUE,
                                priorICweights = c(complement = .3, Heq = .2, Hp = .5)))
  up  <- suppressMessages(evSyn(list(e1, e2), VCOV = list(V, V),
                                hypotheses = list(Hp = "x > y"), Heq = TRUE,
                                priorICweights = c(Complement = .3, Heq = .2, Hp = .5)))
  expect_equal(low$Cumulative_GORICA_weights, base$Cumulative_GORICA_weights)
  expect_equal(up$Cumulative_GORICA_weights, base$Cumulative_GORICA_weights)
  # zonder Heq: 'complement' en 'Complement' beide goed, verkeerde naam fout
  b2 <- suppressMessages(evSyn(list(e1, e2), VCOV = list(V, V), hypotheses = list(Hp = "x > y"),
                               priorICweights = c(complement = .7, Hp = .3)))
  expect_equal(b2$priorICweights, c(.3, .7))
  expect_error(suppressMessages(evSyn(list(e1, e2), VCOV = list(V, V), hypotheses = list(Hp = "x > y"),
                                      priorICweights = c(compl = .7, Hp = .3))), "names of 'priorICweights'")
})

# FX6-A: leave1studyout() herberekent op de logschaal (geen underflow) ------------
test_that("leave1studyout: icweights met extreme gewichten en priors geeft dezelfde uitkomst als evSyn op S-1 studies", {
  W <- list(c(1e-250, 1), c(1, 1e-200), c(1, 1e-200), c(.5, .5))
  p <- c(1e-100, 1)
  e <- suppressMessages(evSyn(W, input_type = "icweights", priorICweights = p))
  # de oorspronkelijke log-gewichten (zonder priors) zitten in het object
  expect_equal(unname(e$logW_m), log(do.call(rbind, W)))
  l <- leave1studyout(e)
  # handberekening (logschaal, in eenheden log(10)): weglaten van studie 4
  #   H1: log(1e-250) + log(1e-100) = -350 log(10); H2: 2 log(1e-200) = -400 log(10)
  #   -> H1 wint met een factor 1e50; IC-verschillen (-2 log W): (500, 800) log(10)
  expect_equal(unname(l$OverallGorica[4, ]), c(500, 800) * log(10))
  expect_equal(unname(l$OverallPrefHypo[, 1]), c("H1", "H2", "H2", "H1"))
  expect_equal(unname(l$OverallGoricaWeights[4, 1]), 1)
  expect_equal(log10(l$OverallGoricaWeights[4, 2] / l$OverallGoricaWeights[4, 1]), -50)
  # weglaten van studie 1: H1: -100 log(10) + log(.5); H2: -400 log(10) + log(.5)
  expect_equal(log10(l$OverallGoricaWeights[1, 2] / l$OverallGoricaWeights[1, 1]), -300)
  # weglaten van studie 2 of 3: H1: -350 log(10) + log(.5); H2: -200 log(10) + log(.5)
  expect_equal(log10(l$OverallGoricaWeights[2, 1] / l$OverallGoricaWeights[2, 2]), -150)
  expect_equal(log10(l$OverallGoricaWeights[3, 1] / l$OverallGoricaWeights[3, 2]), -150)
  # gelijk aan evSyn() op de overige studies (vergeleken op de logschaal)
  for (s in 1:4) {
    d <- suppressMessages(evSyn(W[-s], input_type = "icweights", priorICweights = p))
    expect_equal(log(unname(l$OverallGoricaWeights[s, ])), log(final_w(d)))
  }
  # een (oud) object zonder de log-gewichten wordt niet uit de genormaliseerde
  # (prior-gewogen) gewichten gereconstrueerd, maar geeft een duidelijke fout
  e_old <- e
  e_old$logW_m <- NULL
  expect_error(leave1studyout(e_old), "re-create")
})

test_that("leave1studyout: alle IC-routes en type_ev komen overeen met evSyn op S-1 studies (ook met study_weights)", {
  W  <- list(c(1e-250, 1, 1e-10), c(1, 1e-200, 1e-5), c(.2, .3, .5), c(1e-20, 1, 1e-20))
  W  <- lapply(W, function(w) w / sum(w))
  IC <- lapply(W, function(w) -2 * log(w))
  R  <- lapply(W, function(w) w / w[2])
  sw <- c(.5, 1.5, 1, 1)
  pr <- c(1e-30, .5, .5)
  inputs <- list(icweights = W, icvalues = IC, icratios = R)
  for (nm in names(inputs)) {
    for (te in c("added", "average")) {
      e <- suppressMessages(evSyn(inputs[[nm]], input_type = nm, type_ev = te,
                                  priorICweights = pr, study_weights = sw))
      l <- leave1studyout(e)
      for (s in 1:4) {
        d <- suppressMessages(evSyn(inputs[[nm]][-s], input_type = nm, type_ev = te,
                                    priorICweights = pr, study_weights = sw[-s]))
        expect_equal(log(unname(l$OverallGoricaWeights[s, ])), log(final_w(d)),
                     info = paste(nm, te, s))
        expect_equal(unname(l$OverallPrefHypo[s, 1]),
                     colnames(d$Cumulative_GORICA_weights)[which.max(final_w(d))],
                     info = paste(nm, te, s))
      }
    }
  }
  # icweights, icvalues en icratios geven onderling dezelfde leave-one-out gewichten
  lw <- leave1studyout(suppressMessages(evSyn(W, input_type = "icweights", priorICweights = pr, study_weights = sw)))
  lv <- leave1studyout(suppressMessages(evSyn(IC, input_type = "icvalues", priorICweights = pr, study_weights = sw)))
  lr <- leave1studyout(suppressMessages(evSyn(R, input_type = "icratios", priorICweights = pr, study_weights = sw)))
  expect_equal(log(lw$OverallGoricaWeights), log(lv$OverallGoricaWeights))
  expect_equal(log(lw$OverallGoricaWeights), log(lr$OverallGoricaWeights))
  # est-route, alle type_ev
  est <- list(c(x = 2, y = 1), c(x = 1.5, y = 1.2), c(x = .5, y = .7))
  V <- list(diag(2) * .1, diag(2) * .2, diag(2) * .05)
  H <- list(H1 = "x > y")
  for (te in c("added", "equal", "average")) {
    e <- suppressMessages(evSyn(est, VCOV = V, hypotheses = H, type_ev = te,
                                priorICweights = c(.2, .8), study_weights = c(1, 2, 1)))
    l <- leave1studyout(e)
    for (s in 1:3) {
      d <- suppressMessages(evSyn(est[-s], VCOV = V[-s], hypotheses = H, type_ev = te,
                                  priorICweights = c(.2, .8), study_weights = c(1, 2, 1)[-s]))
      expect_equal(unname(l$OverallGoricaWeights[s, ]), final_w(d), info = paste(te, s))
      expect_equal(unname(l$OverallGorica[s, ]), unname(d$Cumulative_GORICA["Final", ]), info = paste(te, s))
    }
  }
})
