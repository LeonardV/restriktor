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
