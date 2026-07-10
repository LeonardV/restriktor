## Randgevallen voor het vervangen van con_weights_boot (quadprog) door NNLS:
##   1. duale vorm == primale vorm == exact (meq = 0, W regulier)
##   2. g > p: W = A V A' singulier, chol(W) faalt; duale vorm matcht boot
##   3. meq > 0 via reductie matcht exact en boot
## Draaien vanuit de root van de repo:
##   Rscript benchmarks/ic_weights_level_probabilities/check_edge_cases.R

suppressPackageStartupMessages({
  library(mvtnorm); library(quadprog); library(tmvtnorm); library(nnls)
  library(future); library(future.apply)
})
plan(sequential)
source("R/ext_ic_infer.R")
source("R/compute_chiBarSquare_weights.R")

## primale NNLS (zoals prototype): vereist W positief-definiet
w_nnls_primal <- function(W, R = 2e5L, seed = 1, tol = 1e-10) {
  set.seed(seed); g <- ncol(W)
  U <- chol(W); M <- t(backsolve(U, diag(g)))
  half <- R %/% 2L
  Z <- matrix(rnorm(g * half), nrow = g); Z <- cbind(Z, -Z)
  counts <- integer(g + 1L)
  for (i in seq_len(2L * half)) {
    k <- sum(nnls::nnls(M, Z[, i])$x > tol)
    counts[k + 1L] <- counts[k + 1L] + 1L
  }
  counts / (2 * half)
}

## duale NNLS in theta-ruimte: werkt voor elke Amat (ook g > p, redundante
## rijen) en telt actieve constraints precies zoals solve.QP in
## con_weights_boot: dimL = p - #actieve ongelijkheden.
## projectie: min (z-th)'V^-1(z-th) s.t. A th >= 0 heeft duale vorm
##   min_{lam>=0} || ztilde - D lam ||^2, D = L' A' (p x g), V = L L'
w_nnls_dual <- function(VCOV, Amat, R = 2e5L, seed = 1, tol = 1e-10) {
  set.seed(seed)
  p <- ncol(VCOV); g <- nrow(Amat)
  L <- t(chol(VCOV))          # V = L L'
  D <- t(L) %*% t(Amat)       # p x g
  half <- R %/% 2L
  Z <- matrix(rnorm(p * half), nrow = p); Z <- cbind(Z, -Z)
  counts <- integer(p + 1L)
  for (i in seq_len(2L * half)) {
    iact <- sum(nnls::nnls(D, Z[, i])$x > tol)
    dimL <- p - iact
    counts[dimL + 1L] <- counts[dimL + 1L] + 1L
  }
  counts / (2 * half)
}

cat("== check 1: duale vorm == primale vorm == exact (g=6, p=9, meq=0) ==\n")
set.seed(7)
p <- 9; g <- 6
X <- matrix(rnorm(3*p*p), ncol = p); V <- crossprod(X)/(3*p)
A <- matrix(0, g, p); for (i in 1:g) { A[i,i] <- -1; A[i,i+1] <- 1 }
W <- A %*% V %*% t(A)
exact <- rev(ic_weights(W, tolerance = 1e-15, ridge_constant = 1e-5))
wp <- w_nnls_primal(W)
wd <- w_nnls_dual(V, A)
cat("max|primal - exact| =", max(abs(wp - exact)), "\n")
cat("max|dual(tail) - exact| =", max(abs(tail(wd, g+1) - exact)), "\n")
cat("dual: massa onder level p-g (moet 0 zijn):", sum(head(wd, p-g)), "\n")

cat("\n== check 2: g > p (W singulier): chol(W) faalt, dual + boot werken ==\n")
set.seed(11)
p <- 4; g <- 6
X <- matrix(rnorm(3*p*p), ncol = p); V <- crossprod(X)/(3*p)
A <- rbind(diag(p), c(-1,1,0,0), c(0,0,-1,1))   # 6 constraints, 4 params
W <- A %*% V %*% t(A)
cat("rank(W) =", qr(W)$rank, "van", g, "; chol(W):",
    class(try(chol(W), silent = TRUE))[1], "\n")
wb <- con_weights_boot(V, A, meq = 0, R = 2e5L, seed = 42,
                       convergence_crit = 0, chunk_size = 2e5L)
wd <- w_nnls_dual(V, A)
cat("boot: ", round(as.numeric(wb), 4), "\n")
cat("dual: ", round(wd, 4), "\n")
cat("max|dual - boot| =", max(abs(wd - as.numeric(wb))), "(MC-fout, 2 onafh. runs)\n")

cat("\n== check 3: meq > 0 via reductie, vs boot en exact ==\n")
set.seed(13)
p <- 8; g <- 5; meq <- 2
X <- matrix(rnorm(3*p*p), ncol = p); V <- crossprod(X)/(3*p)
A <- matrix(0, g, p); for (i in 1:g) { A[i,i] <- -1; A[i,i+1] <- 1 }
W <- A %*% V %*% t(A)
## exacte route (con_weights): ic_weights op gereduceerde matrix
Wred <- solve(solve(W)[-(1:meq), -(1:meq)])
exact <- rev(ic_weights(Wred, tolerance = 1e-15, ridge_constant = 1e-5)) # levels 0:(g-meq)
wp <- w_nnls_primal(Wred)
cat("max|nnls(reductie) - exact| =", max(abs(wp - exact)), "\n")
## boot in theta-ruimte, levels 0:p; bij meq=2 zit alle massa op levels
## (p-g):(p-meq) = 3:6
wb <- con_weights_boot(V, A, meq = meq, R = 2e5L, seed = 42,
                       convergence_crit = 0, chunk_size = 2e5L)
cat("boot levels 0:p     :", round(as.numeric(wb), 4), "\n")
cat("exact op levels 3:6 :", round(exact, 4), "\n")
cat("max|boot[4:7] - exact| =", max(abs(as.numeric(wb)[4:7] - exact)), "\n")
