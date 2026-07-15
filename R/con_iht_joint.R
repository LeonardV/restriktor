# joint informative hypothesis test for the type A and type B statistics.
#
# Under the shared (least favorable) null hypothesis Rb = r, Wolak (1989,
# Theorem 1) gives the joint distribution of the type B (IU) and type A (EI)
# statistics:
#
#   Pr[IU >= c1, EI >= c2] = sum_{k=0}^{L} w(L, k, Omega) *
#                            Pr[chi2_{q-k} >= c1] * Pr[chi2_k >= c2],
#
# where q = total number of constraints, L = q - meq the number of inequality
# constraints, and w(L, k, Omega) the chi-bar-square weights that restriktor
# already computes (con_weights()/ic_weights(), cached as object$wt.bar).
# Conditional on k, IU and EI are independent. For test statistics that are
# weighted by an estimated error variance the F-bar analogue applies
# (Wolak, 1987, Theorem 4.4); see pfbar_joint() in con_pbar.R, which accounts
# for the shared variance estimate in both statistics.
#
# The type A and type B statistics themselves are computed by the existing
# conTest() machinery (conTest_conLM.R / conTest_conGLM.R / conTest_conRLM.R),
# so iht_joint() works for lm, glm and rlm fits and for all test statistics
# ("F", "LRT", "score", "Wald") that the marginal tests support. For lm the
# F-bar distribution is exact; for glm and rlm it holds asymptotically, as for
# the marginal tests.
#
# REFS:
# Wolak, F. (1989). Local and global testing of linear and nonlinear
#   inequality constraints in nonlinear econometric models. Journal of
#   Econometrics, 41, 205-235.
# Wolak, F. (1987). An exact test for multiple inequality and equality
#   constraints in the linear regression model. Journal of the American
#   Statistical Association, 82, 782-793.
iht_joint <- function(object, constraints = NULL, test = "F",
                      rhs = NULL, neq = 0L, alpha = 0.05, ...) {

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("restriktor ERROR: alpha must be a single number in (0, 1).")
  }

  # unfitted lm/glm/rlm object: fit restriktor first (same dispatch as conTest)
  if (inherits(object, c("lm", "rlm", "glm")) && !inherits(object, "restriktor")) {
    if (is.null(constraints)) {
      stop("restriktor ERROR: no constraints found.")
    }
    ldots <- list(...)
    if (is.null(ldots$se)) {
      ldots$se <- "none"
    }
    m.restr <- match(names(ldots), c("se", "B", "mix_weights",
                                     "parallel", "ncpus", "cl", "seed", "control",
                                     "verbose", "debug"), 0L)
    CALL.restr <- c(list(object = object, constraints = constraints, rhs = rhs,
                         neq = neq), ldots[m.restr > 0L])
    object <- do.call("restriktor", CALL.restr)
  }

  if (!inherits(object, "restriktor")) {
    stop("restriktor ERROR: object must be of class restriktor, lm, glm or rlm.")
  }
  if (!inherits(object, c("conLM", "conGLM", "conRLM"))) {
    stop("restriktor ERROR: iht_joint() is only available for conLM, conGLM and conRLM objects.")
  }

  Amat <- object$constraints
  bvec <- object$rhs
  meq  <- object$neq

  if (all(Amat == 0)) {
    stop("restriktor ERROR: no constraints specified!")
  }
  if (meq == nrow(Amat)) {
    stop("restriktor ERROR: test not applicable for object with equality restrictions only.")
  }

  ldots <- list(...)
  if (!is.null(ldots$neq.alt) && ldots$neq.alt > 0L) {
    stop("restriktor ERROR: iht_joint() requires neq.alt = 0; the joint distribution",
         " (Wolak, 1989, Theorem 1) is derived under the full null hypothesis Rb = r.")
  }
  if (!is.null(ldots$boot) && ldots$boot != "no") {
    stop("restriktor ERROR: bootstrapped p-values are not available for iht_joint();",
         " the joint p-value requires the chi-bar-square weights (mix_weights).")
  }

  wt.bar <- object$wt.bar
  if (is.null(wt.bar) || is.null(attr(wt.bar, "method")) ||
      attr(wt.bar, "method") == "none") {
    stop("restriktor ERROR: no chi-bar-square weights available. Refit the model",
         " with mix_weights = \"pmvnorm\" or mix_weights = \"boot\".")
  }

  # marginal type A and type B tests via the existing machinery; this reuses
  # the cached wt.bar for the marginal p-values
  m.test <- match(names(ldots), c("control", "verbose"), 0L)
  CALL.test <- c(list(object = object, test = test), ldots[m.test > 0L])
  CALL.test$type <- "A"
  outA <- do.call("conTest", CALL.test)
  CALL.test$type <- "B"
  outB <- do.call("conTest", CALL.test)

  df.residual <- outA$df.residual
  q <- nrow(Amat)

  joint <- con_pvalue_joint(object, Ts.A = outA$Ts, Ts.B = outB$Ts)
  wt.bar.joint <- joint$wt.bar
  pvalue.joint <- joint$pvalue

  # combined decision rule at level alpha, based on the marginal tests:
  #  (a) H0: Rb = r is retained          (type A not rejected)
  #  (b) Rb >= r and Rb != r             (type A rejected, type B not rejected)
  #  (c) H0: Rb >= r is rejected         (type B rejected, regardless of type A)
  decision <- if (!is.na(outB$pvalue) && outB$pvalue < alpha) {
    "c"
  } else if (!is.na(outA$pvalue) && outA$pvalue < alpha) {
    "b"
  } else {
    "a"
  }

  OUT <- list(CON          = object$CON,
              Amat         = Amat,
              bvec         = bvec,
              meq          = meq,
              meq.alt      = 0L,
              iact         = object$iact,
              test         = outA$test,
              Ts.A         = outA$Ts,
              Ts.B         = outB$Ts,
              pvalue.A     = outA$pvalue,
              pvalue.B     = outB$pvalue,
              pvalue.joint = pvalue.joint,
              wt.bar       = wt.bar.joint,
              df.residual  = df.residual,
              b.eqrestr    = outA$b.eqrestr,
              b.restr      = object$b.restr,
              b.unrestr    = object$b.unrestr,
              alpha        = alpha,
              decision     = decision,
              R2.org       = object$R2.org,
              R2.reduced   = object$R2.reduced,
              model.org    = object$model.org)

  class(OUT) <- "conTestJoint"

  OUT
}


# internal workhorse shared by iht_joint() and conTest_summary(): the joint
# tail probability Pr[IU >= Ts.B, EI >= Ts.A] under the least favorable null
# Rb = r (Wolak, 1989, Theorem 1; F-bar form of Wolak, 1987, Theorem 4.4),
# computed from the cached chi-bar-square weights on the restriktor object.
con_pvalue_joint <- function(object, Ts.A, Ts.B) {
  Amat <- object$constraints
  meq  <- object$neq
  q <- nrow(Amat)

  # align wt.bar with the type A convention wt.bar[k+1] <-> chi2_k, k = 0:L.
  # mix_weights = "pmvnorm" stores exactly L + 1 weights; mix_weights = "boot"
  # stores level probabilities for dimensions 0:p (p = number of parameters),
  # of which the window p-q ... p-meq corresponds to k = 0:L.
  wt.bar <- object$wt.bar
  L <- q - meq
  p <- length(object$b.unrestr)
  if (length(wt.bar) == L + 1L) {
    wt.bar.joint <- as.numeric(wt.bar)
  } else if (length(wt.bar) == p + 1L) {
    wt.bar.joint <- as.numeric(wt.bar[(p - q + 1L):(p - meq + 1L)])
  } else {
    stop("restriktor ERROR: length of wt.bar (", length(wt.bar),
         ") does not match the number of (inequality) constraints.")
  }
  if (abs(sum(wt.bar.joint) - 1) > 1e-04) {
    warning("restriktor WARNING: chi-bar-square weights do not sum to 1 (sum = ",
            round(sum(wt.bar.joint), 6), "); the joint p-value may be inaccurate.")
  }

  pvalue <- pfbar_joint(cIU = Ts.B, cEI = Ts.A, df_total = q,
                        df2 = object$df.residual, wt.bar = wt.bar.joint)

  list(pvalue = pvalue, wt.bar = wt.bar.joint, Ts.A = Ts.A, Ts.B = Ts.B,
       df.residual = object$df.residual)
}


print.conTestJoint <- function(x, digits = max(3, getOption("digits") - 2), ...) {

  if (!inherits(x, "conTestJoint")) {
    stop("x must be of class \"conTestJoint\"")
  }

  Amat <- x$Amat
  meq  <- x$meq
  bvec <- x$bvec
  iact <- x$iact

  fmt.p <- function(p) {
    if (is.na(p)) {
      as.character(NA)
    } else if (p < 1e-04) {
      "<0.0001"
    } else {
      sprintf("%.4f", p)
    }
  }

  cat("\nRestriktor: joint restricted hypothesis test (Wolak, 1989, Theorem 1):\n")

  if (!inherits(x$model.org, "glm")) {
    if (all((x$R2.org - x$R2.reduced) < 1e-08)) {
      cat("\nMultiple R-squared remains", sprintf("%5.3f", x$R2.org), "\n")
    } else {
      cat("\nMultiple R-squared reduced from", sprintf("%5.3f", x$R2.org), "to",
          sprintf("%5.3f", x$R2.reduced), "\n")
    }
  }

  colnames(Amat) <- names(x$b.unrestr)
  out.rest <- cbind(round(Amat, 4), c(rep("   ==", meq),
                                      rep("   >=", nrow(Amat) - meq)), bvec, " ")
  rownames(out.rest) <- paste(seq_len(nrow(out.rest)), ":", sep = "")
  colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs", "active")
  idx <- ncol(out.rest)
  out.rest[, idx] <- "no"
  out.rest[iact, idx] <- "yes"
  out.rest <- as.data.frame(out.rest)

  cat("\nConstraint matrix:\n")
  print(out.rest, quote = FALSE, scientific = FALSE)

  cat("\n\nType A test: H0: all restrictions are equalities (==)", "\n",
      "        vs. HA: at least one inequality restriction is strictly true (>)\n")
  cat("       ", x$test, "-test statistic: ", sprintf("%.4f", x$Ts.A),
      ",   marginal p-value: ", fmt.p(x$pvalue.A), "\n\n", sep = "")

  cat("Type B test: H0: all restrictions hold in the population", "\n",
      "        vs. HA: at least one restriction is violated\n")
  cat("       ", x$test, "-test statistic: ", sprintf("%.4f", x$Ts.B),
      ",   marginal p-value: ", fmt.p(x$pvalue.B), "\n\n", sep = "")

  cat("Joint tail probability under H0: Rb = r (least favorable null):\n")
  cat("        Pr[Type B >= ", sprintf("%.4f", x$Ts.B),
      ", Type A >= ", sprintf("%.4f", x$Ts.A), "] = ",
      fmt.p(x$pvalue.joint), "\n", sep = "")
  cat("        (mixture of products of F-distributions with a shared\n")
  cat("         denominator on", x$df.residual, "residual degrees of freedom)\n\n")

  decision.txt <- switch(x$decision,
    "a" = "H0: Rb = r is retained (neither test rejects).",
    "b" = "Rb >= r and Rb != r: the order-constrained hypothesis is supported\n        (type A rejected, type B not rejected).",
    "c" = "H0: Rb >= r is rejected (type B rejected).")
  cat("Decision rule at alpha = ", format(x$alpha), ":\n        ",
      decision.txt, "\n\n", sep = "")

  invisible(x)
}
