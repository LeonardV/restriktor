## Per-draw equivalentie van de duale NNLS-projectie en solve.QP, met
## GEDEELDE trekkingen. Let op de tekenconventie: NNLS-respons z_tilde hoort
## bij QP-draw z = -U' z_tilde (de lineaire term van het duale probleem
## wisselt van teken); distributioneel is dat irrelevant omdat z_tilde
## symmetrisch verdeeld is en antithetisch wordt getrokken, maar voor een
## per-draw vergelijking moet de koppeling exact zijn.
## Resultaat: 100% overeenstemming in alle scenario's (20.000 draws elk).
## Draaien vanuit de root van de repo:
##   Rscript benchmarks/ic_weights_level_probabilities/check_shared_draws.R

suppressPackageStartupMessages({ library(quadprog); library(nnls) })

count_qp <- function(z, V, A, meq) {
  invW <- solve(V)
  QP <- try(solve.QP(Dmat = 2*invW, dvec = 2*(invW %*% z), Amat = t(A),
                     bvec = rep(0, nrow(A)), meq = meq), silent = TRUE)
  if (inherits(QP, "try-error")) return(NA_integer_)
  if (QP$iact[1] == 0L) 0L else length(QP$iact)
}
count_nnls <- function(zt, V, A, meq, tol = 1e-10) {
  U <- chol(V); D <- U %*% t(A)
  if (meq > 0) {
    qr_eq <- qr(D[, 1:meq, drop = FALSE])
    D2 <- qr.resid(qr_eq, D[, -(1:meq), drop = FALSE])
    zt <- qr.resid(qr_eq, zt)
  } else D2 <- D
  as.integer(meq + sum(nnls(D2, zt)$x > tol))
}

## correcte koppeling: NNLS-respons zt hoort bij QP-draw z = -U' zt
## (dual: min_{lam>=0} lam' W lam + 2 lam' A z  ==  || zt - D lam ||^2
##  met D = U A' en A z = -D' zt)
run_shared <- function(V, A, meq, n = 20000, label) {
  p <- ncol(V); U <- chol(V)
  set.seed(99)
  Zt <- matrix(rnorm(p * n), nrow = p)
  Z  <- -crossprod(U, Zt)                 # z = -U' zt
  cq <- integer(n); cn <- integer(n)
  for (i in 1:n) {
    cq[i] <- count_qp(Z[, i], V, A, meq)
    cn[i] <- count_nnls(Zt[, i], V, A, meq)
  }
  qp_na <- sum(is.na(cq))
  agree <- mean(cq[!is.na(cq)] == cn[!is.na(cq)])
  cat(sprintf("%-38s  overeenstemming = %.4f%%  (QP-fouten: %d/%d)\n",
              label, 100*agree, qp_na, n))
  if (agree < 1) {
    idx <- which(!is.na(cq) & cq != cn)
    cat("   aantal afwijkende draws:", length(idx),
        "| QP vs NNLS:", paste(head(cq[idx], 6), head(cn[idx], 6), sep="/", collapse=" "), "\n")
  }
}

set.seed(3)
p <- 8; g <- 5
X <- matrix(rnorm(3*p*p), ncol = p); V <- crossprod(X)/(3*p)
A <- matrix(0, g, p); for (i in 1:g) { A[i,i] <- -1; A[i,i+1] <- 1 }

run_shared(V, A, 0, label = "generiek, meq = 0")
run_shared(V, A, 2, label = "gelijkheden, meq = 2")
run_shared(V, rbind(A, A[1, ] + 1e-6 * rnorm(p)), 0, label = "bijna-afhankelijke rij (1e-6)")
run_shared(V, rbind(A, A[1, ]), 0, label = "gedupliceerde rij (exact redundant)")
run_shared(V, diag(p), 0, label = "p orthant-restricties (veel actief)")
