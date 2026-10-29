# group_means(): group (cell) means in terms of the lm parameters

tg <- ToothGrowth
tg$dose <- factor(tg$dose)
tg <- tg[-c(1:3, 35), ] # unbalanced
cell_means <- aggregate(len ~ supp + dose, tg, mean)
cell_means <- cell_means[order(cell_means$supp, cell_means$dose), "len"]
implied <- function(fit) {
  gm <- suppressWarnings(group_means(fit))
  unname(drop(gm$L %*% coef(fit)[!is.na(coef(fit))]))
}

test_that("group means equal the cell means for any coding", {
  expect_equal(implied(lm(len ~ supp * dose, tg)), cell_means)
  expect_equal(implied(lm(len ~ 0 + supp:dose, tg)), cell_means)
  expect_equal(implied(lm(len ~ supp * dose, tg,
                          contrasts = list(supp = "contr.sum", dose = "contr.sum"))),
               cell_means)
  expect_equal(implied(lm(len ~ supp * dose, tg,
                          contrasts = list(supp = "contr.helmert", dose = "contr.helmert"))),
               cell_means)
})

test_that("notation and group sizes", {
  gm <- group_means(lm(len ~ supp * dose, tg))
  expect_equal(gm$table$name, c("mu_11", "mu_12", "mu_13", "mu_21", "mu_22", "mu_23"))
  expect_equal(sum(gm$table$n), nrow(tg))
  expect_equal(group_means(lm(len ~ supp * dose, tg), labels = "level")$table$name[1],
               "mu_OJ_0.5")
})

test_that("goric with group means: same results for dummy and effects coding", {
  H <- list(H1 = "mu_11 < mu_12 < mu_13; mu_21 < mu_22 < mu_23")
  g_d <- goric(lm(len ~ supp * dose, tg), hypotheses = H)
  g_e <- goric(lm(len ~ supp * dose, tg,
                  contrasts = list(supp = "contr.sum", dose = "contr.sum")),
               hypotheses = H)
  expect_equal(g_d$result$goric, g_e$result$goric, tolerance = 1e-4)
})

test_that("warnings for factors coded by hand", {
  d <- tg
  d$d_e1 <- ifelse(d$dose == "0.5", 1, ifelse(d$dose == "2", -1, 0))
  d$d_e2 <- ifelse(d$dose == "1", 1, ifelse(d$dose == "2", -1, 0))
  d$d1 <- as.numeric(d$dose == "1")
  d$d2 <- as.numeric(d$dose == "2")
  d$sex <- rep(0:1, length.out = nrow(d))
  expect_warning(group_means(lm(len ~ supp + d_e1 + d_e2, d)), "effects-coded")
  expect_warning(group_means(lm(len ~ d1 + d2, d)), "dummy-coded")
  # no warning for factors (any coding) and a single 0/1 dummy
  expect_no_warning(group_means(lm(len ~ supp * dose + sex, d,
                                   contrasts = list(dose = "contr.sum"))))
})
