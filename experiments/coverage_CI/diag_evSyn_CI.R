# Diagnostic for the evSyn_CI() condition with a low (strict) coverage for H1 versus H2
# (parameters in the opposite order, S = 3, n = 100): strict versus tolerant coverage and degenerate CIs.
suppressPackageStartupMessages(library(restriktor))
set.seed(5)
R <- 200; B <- 500; tol <- 0.1; S <- 3; n <- 100
th <- c(b1 = 0, b2 = .2, b3 = .4); V <- diag(1 / n, 3)
H <- list(H1 = "b1 > b2 > b3", H2 = "b1 > b2 = b3")
g0 <- goric(th, VCOV = V, hypotheses = H)
EIC <- colMeans(restriktor:::goric_CI_boot(g0, as.character(g0$result$model), B = 2e5))
target <- S * (EIC[1] - EIC[2])
b <- t(replicate(R, {
  est <- lapply(1:S, function(s) setNames(th + rnorm(3, 0, sqrt(1 / n)), names(th)))
  ci <- evSyn_CI(suppressMessages(evSyn(est, VCOV = rep(list(V), S), hypotheses = H)), B = B)$difference
  j <- which(ci$model == "H1" & ci$vs == "H2"); c(ci$lower[j], ci$upper[j])
}))
d <- data.frame(S = S, n = n, comparison = "H1 vs. H2", theta = "outside", target = target,
                strict = mean(b[, 1] <= target & target <= b[, 2]),
                tolerant = mean(b[, 1] - tol <= target & target <= b[, 2] + tol),
                degenerate = mean(b[, 2] - b[, 1] < 1e-6))
print(d, row.names = FALSE)
saveRDS(d, "diag_evSyn_CI.rds")
