con_wt_lm <- function(cov, meq) {
  if (meq == 0L) {
    wt.bar <- con_weights_lm(cov)
  }
  else if (meq > 0) {
    wt.bar <- con_weights_lm(solve(solve(cov)[-(1:meq), -(1:meq)]))
  }
  wt.bar
}

