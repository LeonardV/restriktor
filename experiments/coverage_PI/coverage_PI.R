# Coverage of the prediction interval (evSyn_PI) for the ratio of GORICA weights in a new study.
# Resumable per condition (results saved per condition).
suppressPackageStartupMessages(library(restriktor))
args <- commandArgs(TRUE); R <- if (length(args)) as.integer(args[1]) else 200
B <- 200
outdir <- "/home/claude/work/sim/pi/res"; dir.create(outdir, showWarnings = FALSE)
sets <- list(
  comp = list(H = list(H1 = "b1 > b2 > b3"), comparison = "complement"),
  unc  = list(H = list(H1 = "b1 > b2 > b3", H2 = "b1 > b2 = b3"), comparison = "unconstrained"))
theta0 <- c(b1 = .25, b2 = .12, b3 = 0)
conds <- expand.grid(set = names(sets), S = c(5, 10, 20), tau = c(0, .10), stringsAsFactors = FALSE)
draw_study <- function(n, tau) {
  th <- theta0 + rnorm(3, 0, tau)                    # study-specific true values (heterogeneity)
  V <- diag(1 / n, 3)
  list(est = setNames(th + rnorm(3, 0, sqrt(1 / n)), names(theta0)), V = V, n = n)
}
for (ci in seq_len(nrow(conds))) {
  cd <- conds[ci, ]; f <- file.path(outdir, sprintf("%s_S%d_tau%.2f.rds", cd$set, cd$S, cd$tau))
  if (file.exists(f)) next
  hs <- sets[[cd$set]]; set.seed(1000 + ci)
  res <- matrix(NA, R, 6, dimnames = list(NULL, c("covered", "lower", "upper", "obs", "Pfav", "fav_new")))
  for (r in seq_len(R)) {
    ns <- sample(c(50, 100, 200, 400), cd$S, replace = TRUE)
    st <- lapply(ns, draw_study, tau = cd$tau)
    ev <- tryCatch(suppressMessages(evSyn(lapply(st, `[[`, "est"), VCOV = lapply(st, `[[`, "V"),
                        hypotheses = hs$H, comparison = hs$comparison, study_sample_nobs = ns)),
                   error = function(e) NULL)
    if (is.null(ev)) next
    pi <- tryCatch(suppressMessages(evSyn_PI(ev, B = B)), error = function(e) NULL)
    if (is.null(pi)) next
    best <- pi$best; row <- pi$PI[1, ]                # best vs. first other hypothesis
    n_new <- row$n_new
    nw <- draw_study(n_new, cd$tau)
    g <- suppressMessages(goric(nw$est, VCOV = nw$V, hypotheses = hs$H, comparison = hs$comparison))
    w <- g$result$gorica.weights; mods <- as.character(g$result$model)
    lab <- tolower(c(row$model, row$vs)); im <- match(lab, tolower(mods))
    obs <- w[im[1]] / w[im[2]]
    eps <- 1e-8
    res[r, ] <- c(row$lower * (1 - eps) <= obs & obs <= row$upper * (1 + eps), row$lower, row$upper, obs,
                  row$P_favours_best, obs > 1)
  }
  saveRDS(list(cond = cd, res = res), f)
  cat(sprintf("%s S=%d tau=%.2f: coverage %.3f (n=%d); mean P(fav) %.3f vs observed %.3f\n", cd$set, cd$S, cd$tau,
              mean(res[, 1], na.rm = TRUE), sum(!is.na(res[, 1])), mean(res[, 5], na.rm = TRUE), mean(res[, 6], na.rm = TRUE)))
}
