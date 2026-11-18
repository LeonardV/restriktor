suppressPackageStartupMessages(library(restriktor))
set.seed(11)
R <- 250; B <- 500
thetas <- list(null = c(b1 = 0, b2 = 0, b3 = 0), boundary = c(b1 = .2, b2 = 0, b3 = 0),
               outside = c(b1 = 0, b2 = .2, b3 = .4), inside = c(b1 = .4, b2 = .2, b3 = 0))
hs <- list(A = list(H1 = "b1 > b2 > b3", H2 = "b1 > b2 = b3"), B = list(H1 = "b1 > b2 > b3"))
out <- list()
for (hn in names(hs)) for (tn in names(thetas)) for (S in c(3, 10)) {
  n <- 25; V <- diag(1 / n, 3); th <- thetas[[tn]]
  g0 <- goric(th, VCOV = V, hypotheses = hs[[hn]]); mods <- as.character(g0$result$model); M <- length(mods)
  EIC <- colMeans(restriktor:::goric_CI_boot(g0, mods, B = 1e5))
  target <- S * (EIC[1] - EIC[M])                         # pooled difference H1 vs failsafe (added)
  cov <- matrix(NA, R, 3, dimnames = list(NULL, c("percentile", "centered_perc", "centered_basic")))
  for (r in 1:R) {
    Tobs <- 0; Tb <- 0
    for (s in 1:S) {
      g <- goric(th + rnorm(3, 0, sqrt(1 / n)), VCOV = V, hypotheses = hs[[hn]])
      ICb <- restriktor:::goric_CI_boot(g, mods, B = B)
      Tobs <- Tobs + (g$result$gorica[1] - g$result$gorica[M])
      Tb <- Tb + (ICb[, 1] - ICb[, M])
    }
    q <- quantile(Tb, c(.025, .975)); Z <- Tb - mean(Tb)
    qz <- quantile(Z, c(.025, .975))
    cov[r, ] <- c(q[1] <= target & target <= q[2],
                  Tobs + qz[1] <= target & target <= Tobs + qz[2],
                  Tobs - qz[2] <= target & target <= Tobs - qz[1])
  }
  out[[length(out) + 1]] <- data.frame(set = hn, theta = tn, S = S, t(colMeans(cov)))
  print(tail(out, 1)[[1]], row.names = FALSE)
}
saveRDS(do.call(rbind, out), "proto_centered.rds")
