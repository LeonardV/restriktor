## Onderzoek: kunnen de level probabilities (chi-bar-kwadraat gewichten) in
## ic_weights() sneller en schaalbaar berekend worden?
##
## Dit script bevat het prototype con_weights_nnls() en reproduceert de
## benchmarks uit README.md. Draaien vanuit de root van de repo:
##   Rscript benchmarks/ic_weights_level_probabilities/bench_ic_weights.R

suppressPackageStartupMessages({
  library(mvtnorm)
  library(quadprog)
  library(tmvtnorm)
  library(nnls)
  library(future)
  library(future.apply)
})
plan(sequential)

source("R/ext_ic_infer.R")
source("R/compute_chiBarSquare_weights.R")

## ---------------------------------------------------------------------------
## Prototype: level probabilities via NNLS-simulatie
##
## De level probabilities zijn
##   w_k(W) = P( projectie van s ~ N(0, W) op {b >= 0} in metriek W^{-1}
##               heeft precies k positieve componenten ),  W = A V A' (g x g).
## Met Cholesky W = U'U, M = U^{-T} en z = M s ~ N(0, I) geldt
##   min_{b>=0} (s-b)' W^{-1} (s-b)  ==  min_{b>=0} || z - M b ||^2,
## een standaard NNLS-probleem met vaste designmatrix M. Elke draw kost een
## Lawson-Hanson solve (gecompileerd Fortran) in g dimensies, onafhankelijk
## van het aantal parameters p. Antithetische paren (z, -z) reduceren de
## variantie. Voor meq > 0: pas toe op solve(solve(W)[-(1:meq), -(1:meq)]),
## zoals nu in con_weights().
## ---------------------------------------------------------------------------
con_weights_nnls <- function(W, R = 1e5L, seed = NULL, antithetic = TRUE,
                             tol = 1e-10) {
  if (!is.null(seed)) set.seed(seed)
  g <- ncol(W)
  U <- chol(W)                       # W = U'U
  M <- t(backsolve(U, diag(g)))      # M = U^{-T}

  if (antithetic) {
    half <- ceiling(R / 2)
    Z <- matrix(rnorm(g * half), nrow = g)
    Z <- cbind(Z, -Z)[, seq_len(R), drop = FALSE]
  } else {
    Z <- matrix(rnorm(g * R), nrow = g)
  }

  counts <- integer(g + 1L)
  for (i in seq_len(R)) {
    k <- sum(nnls::nnls(M, Z[, i])$x > tol)
    counts[k + 1L] <- counts[k + 1L] + 1L
  }
  w <- counts / R
  names(w) <- 0:g                    # oplopend, zoals rev(ic_weights(W))
  w
}

## testprobleem: ordering-constraints b1 < b2 < ... met ruwe random covariantie
make_problem <- function(g, p = g + 3, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(3 * p * p), ncol = p)
  VCOV <- crossprod(X) / (3 * p)
  Amat <- matrix(0, g, p)
  for (i in 1:g) { Amat[i, i] <- -1; Amat[i, i + 1] <- 1 }
  W <- Amat %*% VCOV %*% t(Amat)
  list(VCOV = VCOV, Amat = Amat, W = W)
}

## ---------------------------------------------------------------------------
## 1) Correctheid: NNLS-schatter vs exacte ic_weights
## ---------------------------------------------------------------------------
cat("== correctheid (R = 2e5, antithetic) ==\n")
for (g in c(4, 6, 8, 10)) {
  pr <- make_problem(g)
  exact <- rev(ic_weights(pr$W, tolerance = 1e-15, ridge_constant = 1e-5))
  est <- con_weights_nnls(pr$W, R = 2e5L, seed = 42)
  cat(sprintf("g=%2d  max|diff| = %.5f   |E[k] diff| = %.5f\n", g,
              max(abs(exact - est)),
              abs(sum(0:g * exact) - sum(0:g * est))))
}

## ---------------------------------------------------------------------------
## 2) Vergelijking met bestaande con_weights_boot (quadprog) bij g = 8
## ---------------------------------------------------------------------------
cat("\n== NNLS vs bestaande con_weights_boot (g = 8) ==\n")
pr <- make_problem(8)
t_boot <- system.time(
  wb <- con_weights_boot(pr$VCOV, pr$Amat, meq = 0, R = 1e5L, seed = 42,
                         convergence_crit = 1e-3, chunk_size = 5000L)
)["elapsed"]
t_nnls <- system.time(wn <- con_weights_nnls(pr$W, R = 1e5L, seed = 42))["elapsed"]
exact <- rev(ic_weights(pr$W, tolerance = 1e-15, ridge_constant = 1e-5))
cat(sprintf("boot:  %5.2fs (draws=%d)  max|diff vs exact| = %.5f\n", t_boot,
            attr(wb, "total_bootstrap_draws"),
            max(abs(tail(as.numeric(wb), 9) - exact))))
cat(sprintf("nnls:  %5.2fs (draws=%d)  max|diff vs exact| = %.5f\n",
            t_nnls, 1e5L, max(abs(wn - exact))))

## ---------------------------------------------------------------------------
## 3) Schaalbaarheid: tijd vs aantal constraints g
## ---------------------------------------------------------------------------
cat("\n== timings exact vs NNLS ==\n")
for (g in c(6, 8, 10, 12, 14)) {
  pr <- make_problem(g)
  t_exact <- system.time(
    we <- ic_weights(pr$W, tolerance = 1e-15, ridge_constant = 1e-5)
  )["elapsed"]
  t_nnls <- system.time(wn <- con_weights_nnls(pr$W, R = 1e5L, seed = 1))["elapsed"]
  cat(sprintf("g=%2d  exact=%8.2fs   nnls(R=1e5)=%5.2fs   max|diff|=%.5f\n",
              g, t_exact, t_nnls, max(abs(rev(we) - wn))))
}

cat("\n== NNLS bij grote g (exact niet meer haalbaar) ==\n")
for (g in c(20, 30, 50, 75, 100)) {
  pr <- make_problem(g)
  t_nnls <- system.time(wn <- con_weights_nnls(pr$W, R = 1e5L, seed = 1))["elapsed"]
  cat(sprintf("g=%3d  nnls(R=1e5)=%6.2fs   E[k]=%.3f\n", g, t_nnls, sum(0:g * wn)))
}

## ---------------------------------------------------------------------------
## 4) p-afhankelijkheid: boot werkt in theta-ruimte (p), NNLS in
##    constraint-ruimte (g)
## ---------------------------------------------------------------------------
cat("\n== p-afhankelijkheid (g = 8, R = 2e4, boot zonder vroege stop) ==\n")
for (p in c(11, 25, 50, 100)) {
  pr <- make_problem(8, p)
  t_boot <- system.time(
    con_weights_boot(pr$VCOV, pr$Amat, meq = 0, R = 2e4L, seed = 42,
                     convergence_crit = 0, chunk_size = 5000L)
  )["elapsed"]
  t_nnls <- system.time(con_weights_nnls(pr$W, R = 2e4L, seed = 42))["elapsed"]
  cat(sprintf("p=%3d  boot=%6.2fs (%5.1f us/draw)   nnls=%5.2fs (%5.1f us/draw)\n",
              p, t_boot, 1e6 * t_boot / 2e4, t_nnls, 1e6 * t_nnls / 2e4))
}

## ---------------------------------------------------------------------------
## 5) Precisie van de goric-penalty PT = 1 + sum(k * w_k) bij g = 10
## ---------------------------------------------------------------------------
cat("\n== precisie penalty (g = 10, 20 herhalingen per R) ==\n")
pr <- make_problem(10, 13)
exact <- rev(ic_weights(pr$W, tolerance = 1e-15, ridge_constant = 1e-5))
PT_exact <- 1 + sum(0:10 * exact)
for (R in c(1e4, 5e4, 1e5)) {
  pts <- replicate(20, 1 + sum(0:10 * con_weights_nnls(pr$W, R = as.integer(R))))
  cat(sprintf("R=%6d  PT_exact=%.4f  mean=%.4f  sd=%.4f  max|err|=%.4f\n",
              R, PT_exact, mean(pts), sd(pts), max(abs(pts - PT_exact))))
}
