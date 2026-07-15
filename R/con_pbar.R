pfbar <- function(x, df1, df2, wt.bar) {
  if (x <= 0) {
    return(0)
  }
  zed <- df1 == 0
  cdf <- ifelse(any(zed), wt.bar[zed], 0)
  cdf <- cdf + sum(pf(x/df1[!zed], df1[!zed], df2) * wt.bar[!zed])
  return(cdf)
}
  
pchibar <- function(x, df, wt.bar) {
  if (x <= 0) {
    return(0)
  }
  zed <- df == 0
  cdf <- ifelse(any(zed), wt.bar[zed], 0)
  cdf <- cdf + sum(pchisq(x, df[!zed]) * wt.bar[!zed])
  return(cdf)
}


# joint upper-tail probability of the type B (IU) and type A (EI) statistics
# under the least favorable null hypothesis Rb = r:
#
#   Pr[IU >= cIU, EI >= cEI] = sum_{k=0}^{L} wt.bar[k+1] *
#                              Pr[chi2_{df_total - k} >= cIU] * Pr[chi2_k >= cEI]
#
# with L = length(wt.bar) - 1 = (number of inequality constraints) and
# df_total = total number of constraints (equalities + inequalities), so that
# conditional on k "positive" components EI ~ chi2_k and IU ~ chi2_{df_total - k},
# and EI and IU are conditionally independent.
# REF: Wolak, F. (1989). Local and global testing of linear and nonlinear
# inequality constraints in nonlinear econometric models. Journal of
# Econometrics, 41, 205-235 (Theorem 1).
#
# wt.bar follows the type A convention as stored on a restriktor object
# (wt.bar[k+1] is the weight of chi2_k in the type A mixture), so that
# pchibar_joint(cIU, 0, ...) equals the marginal type B tail probability and
# pchibar_joint(0, cEI, ...) the marginal type A tail probability.
pchibar_joint <- function(cIU, cEI, df_total, wt.bar) {
  L <- length(wt.bar) - 1L
  k <- 0:L
  df.IU <- df_total - k
  df.EI <- k
  # chi2_0 is a point mass at zero: Pr[chi2_0 >= c] = 1 if c <= 0, else 0
  tail.IU <- ifelse(df.IU == 0, as.numeric(cIU <= 0),
                    pchisq(cIU, df.IU, lower.tail = FALSE))
  tail.EI <- ifelse(df.EI == 0, as.numeric(cEI <= 0),
                    pchisq(cEI, df.EI, lower.tail = FALSE))
  sum(wt.bar * tail.IU * tail.EI)
}


# F-bar analogue of pchibar_joint for test statistics whose quadratic forms are
# weighted by an estimated error variance (Wolak, 1987, Theorem 4.4). Both
# statistics share the same estimate s2 = sigma^2 * V / df2 with V ~ chi2_{df2}
# independent of the numerator quadratic forms, so the exact joint tail is the
# expectation over V of the conditional chi-square product:
#
#   Pr[IU >= cIU, EI >= cEI] = sum_k wt.bar[k+1] *
#      E_V[ Pr(chi2_{df_total-k} >= cIU * V/df2) * Pr(chi2_k >= cEI * V/df2) ]
#
# Note: because the denominator is shared, this is NOT the product of the two
# marginal F-bar tails; the product form is only the limit df2 -> Inf, in which
# case pfbar_joint reduces to pchibar_joint. The marginals are unaffected:
# pfbar_joint(cIU, 0, ...) equals the marginal type B F-bar tail probability
# (consistent with pfbar/con_pvalue_Fbar) and pfbar_joint(0, cEI, ...) the
# marginal type A F-bar tail probability.
# cIU and cEI are on the scale of the Ts statistics in conTest (i.e., the
# quadratic forms themselves; the division by the degrees of freedom that
# pfbar applies is done internally per mixture component).
pfbar_joint <- function(cIU, cEI, df_total, df2, wt.bar) {
  if (cIU <= 0 && cEI <= 0) {
    return(1)
  }
  if (is.infinite(df2)) {
    return(pchibar_joint(cIU, cEI, df_total, wt.bar))
  }
  L <- length(wt.bar) - 1L
  k <- 0:L
  df.IU <- df_total - k
  df.EI <- k

  # E_V[g(V)] with V ~ chi2_{df2} is computed via the probability integral
  # transform, integrate g(qchisq(t, df2)) over t in (0, 1): for large df2 the
  # chi-square density is a sharp peak that adaptive quadrature on (0, Inf)
  # can miss entirely, whereas the transformed integrand is smooth and O(1)
  integrand <- function(t) {
    s <- qchisq(t, df2) / df2
    vapply(s, function(si) {
      tail.IU <- ifelse(df.IU == 0, as.numeric(cIU <= 0),
                        pchisq(cIU * si, df.IU, lower.tail = FALSE))
      tail.EI <- ifelse(df.EI == 0, as.numeric(cEI <= 0),
                        pchisq(cEI * si, df.EI, lower.tail = FALSE))
      sum(wt.bar * tail.IU * tail.EI)
    }, numeric(1))
  }

  out <- integrate(integrand, lower = 0, upper = 1, rel.tol = 1e-09,
                   subdivisions = 500L)
  # tail probabilities live in [0, 1]; guard against tiny numerical overshoot
  min(max(out$value, 0), 1)
}


  
  