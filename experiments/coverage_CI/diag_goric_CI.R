# Diagnostics for the conditions with a (strict) coverage below .90 in the goric_CI() simulation:
# strict coverage, coverage with a small tolerance, and the share of degenerate (point) CIs.
suppressPackageStartupMessages(library(restriktor))
set.seed(4)
R <- 200; B <- 1000; tol <- 0.1   # tolerance on the scale of the difference in GORICA values
thetas <- list(inside = c(b1 = .4, b2 = .2, b3 = 0), boundary = c(b1 = .2, b2 = 0, b3 = 0),
               outside = c(b1 = 0, b2 = .2, b3 = .4), null = c(b1 = 0, b2 = 0, b3 = 0))
sets <- list(A = list(hypotheses = list(H1 = "b1 > b2 > b3", H2 = "b1 > b2 = b3"), comparison = "unconstrained"),
             C = list(hypotheses = list(H1 = "b1 > b2"), comparison = "complement"),
             D = list(hypotheses = list(H1 = "b1 = b2 = b3"), comparison = "unconstrained"))
cases <- data.frame(
  set      = c("A", "A", "A", "A", "A", "C", "C", "D"),
  theta    = c("boundary", "null", "outside", "inside", "inside", "boundary", "inside", "null"),
  n        = c(100, 100, 100, 400, 25, 100, 100, 25),
  quantity = c("H1 vs. H2", "H1 vs. H2", "H1 vs. H2", "H1 vs. failsafe", "H1 vs. failsafe",
               "value H1", "value H1", "H1 vs. failsafe"),
  method   = c("boot", "boot", "boot", "boot", "boot", "boot", "boot", "invert"),
  stringsAsFactors = FALSE)
out <- list()
for (i in seq_len(nrow(cases))) {
  cs <- cases[i, ]; hs <- sets[[cs$set]]; th <- thetas[[cs$theta]]; V <- diag(1 / cs$n, 3)
  g0 <- goric(th, VCOV = V, hypotheses = hs$hypotheses, comparison = hs$comparison)
  mods <- as.character(g0$result$model)
  EIC <- colMeans(restriktor:::goric_CI_boot(g0, mods, B = 2e5))
  k <- if (cs$quantity == "H1 vs. H2") 2 else length(mods)
  target <- if (cs$quantity == "value H1") EIC[1] else EIC[1] - EIC[k]
  b <- t(replicate(R, {
    g <- goric(th + rnorm(3, 0, sqrt(1 / cs$n)), VCOV = V, hypotheses = hs$hypotheses, comparison = hs$comparison)
    ci <- suppressMessages(goric_CI(g, method = cs$method, B = B))
    if (cs$quantity == "value H1") c(ci$IC$lower[1], ci$IC$upper[1])
    else { j <- which(ci$difference$model == mods[1] & ci$difference$vs == mods[k])
           c(ci$difference$lower[j], ci$difference$upper[j]) }
  }))
  out[[i]] <- cbind(cs, target = target,
                    strict = mean(b[, 1] <= target & target <= b[, 2]),
                    tolerant = mean(b[, 1] - tol <= target & target <= b[, 2] + tol),
                    degenerate = mean(b[, 2] - b[, 1] < 1e-6),
                    below = mean(target < b[, 1] - tol), above = mean(target > b[, 2] + tol))
  print(out[[i]], row.names = FALSE)
}
saveRDS(do.call(rbind, out), "diag_goric_CI.rds")
