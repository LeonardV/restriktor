\name{goric_CI}
\alias{goric_CI}
\alias{print.goric_CI}
\title{Confidence Intervals for GORIC(A) Values, Differences, Weights, and Ratios of Weights}
\description{
  \code{goric_CI} computes confidence intervals (CIs) for the GORIC(A) values,
  for the differences in GORIC(A) values between all pairs of hypotheses, for
  the GORIC(A) weights, and for the ratios of GORIC(A) weights, based on a
  fitted \code{goric} object. Analytic (normal-based) CIs and
  parametric-bootstrap CIs are available.
}
\usage{
goric_CI(object, level = 0.95, method = c("invert", "boot", "delta", "mc"),
         approach = c("auto", "casewise", "estimates"),
         nsim = 10000L, B = 1000L, seed = NULL, tol = 1e-8)

\method{print}{goric_CI}(x, digits = max(3, getOption("digits") - 4),
      all_pairs = FALSE, \dots)
}
\arguments{
  \item{object}{an object of class \code{con_goric} (i.e., output from the
  \code{\link{goric}} function).}
  \item{level}{the confidence level (default 0.95).}
  \item{method}{the type of CIs: \code{"delta"} gives analytic
  (normal-based) CIs, where the CIs for the weights are obtained via the delta
  method on the logit scale; \code{"mc"} gives the same analytic CIs, but the
  CIs for the weights are obtained via Monte Carlo draws from the (multivariate
  normal) distribution of the GORIC(A) values; \code{"boot"} gives
  parametric-bootstrap percentile CIs (see Details), which are recommended when
  hypotheses are close or overlap; \code{"invert"} (default; recommended
  in general, see 'Choosing a method' in Details) gives, for
  a hypothesis consisting of only equality constraints or of a single inequality
  constraint versus the unconstrained hypothesis or versus its complement,
  CIs obtained by inverting the sampling distribution of the difference in
  GORIC(A) values (see Details), and bootstrap CIs otherwise.}
  \item{approach}{the approach used for the analytic CIs: \code{"casewise"}
  uses the casewise (individual) log-likelihood contributions (only for
  \code{type = "goric"} or \code{"goricc"}); \code{"estimates"} uses the
  estimates and their covariance matrix (as in the GORICA). By default
  (\code{"auto"}), \code{"casewise"} is used for the GORIC(C) and
  \code{"estimates"} for the GORICA(C).}
  \item{nsim}{the number of Monte Carlo draws when \code{method = "mc"}.}
  \item{B}{the number of draws when \code{method = "boot"} or \code{"invert"} (default 1000; also used in \code{\link{goric}}).}
  \item{seed}{an optional seed (for \code{method = "mc"} or \code{"boot"}).
  Note that this sets the seed of the global random number generator.}
  \item{tol}{tolerance used to determine whether two hypotheses have the same
  order-restricted maximum likelihood estimates (or-mle's) in the sample, that is,
  whether the standard error of the difference in their
  GORIC(A) values is (numerically) zero.}
  \item{x}{an object of class \code{goric_CI}.}
  \item{digits}{the number of digits to be printed.}
  \item{all_pairs}{if \code{FALSE} (default), the differences in GORIC(A)
  values and the ratios of GORIC(A) weights are printed for the best hypothesis
  versus the others; if \code{TRUE}, for all pairs of hypotheses.}
  \item{\dots}{not used.}
}
\details{
  Let \eqn{IC_m = -2 \ell_m + c\, PT_m} be the GORIC(A) value of hypothesis
  \eqn{m}, with \eqn{\ell_m} the (order-restricted) log-likelihood, \eqn{PT_m}
  the penalty, and \eqn{c} the penalty factor (\code{penalty_factor},
  default 2). With prior model probabilities \eqn{p_m} (\code{priorICweights}),
  let \eqn{\eta_m = \log p_m - IC_m/2}; then the GORIC(A) weights are
  \eqn{w_m = \exp(\eta_m) / \sum_k \exp(\eta_k)} and the ratio of weights is
  \eqn{w_m / w_k = \exp(\eta_m - \eta_k)}. The penalties do not depend on the
  outcome values (and the priors are fixed), so all sampling variability
  comes from the log-likelihood values; the penalty factor and priors as used
  in \code{object} are taken into account.

  \emph{Weights and ratios of weights.} Since the weights and the ratios of
  weights are functions of \eqn{\eta} (i.e., of the GORIC(A) values), their CIs
  follow from the distribution of the GORIC(A) values, for each method:
  \itemize{
    \item The ratio of weights \eqn{w_m / w_k = \exp(\eta_m - \eta_k)} is a monotone
    function of the difference \eqn{IC_m - IC_k}, so its CI is the
    transformation \eqn{\exp(\log(p_m/p_k) - \cdot/2)} of the CI for that
    difference (for all methods).
    \item A weight depends on all differences jointly. Its CI is obtained via
    the delta method on the logit scale (\code{"delta"}), via Monte Carlo
    draws from the (multivariate normal) distribution of the GORIC(A) values
    (\code{"mc"}), or as the percentiles of the weights over the bootstrap
    draws (\code{"boot"} and \code{"invert"}). With only two hypotheses, a
    weight is a monotone function of their difference,
    \eqn{w_1 = 1 / (1 + (p_2/p_1) \exp((IC_1 - IC_2)/2))}, so then its CI is
    the transformation of the CI for the difference; for \code{"invert"}, this
    gives the inversion CI for the weights (when the pair is eligible, see
    below). The column \code{CI_method} shows which CI is given.
  }
  The population counterparts of the weights and ratios of weights are those
  based on the expected GORIC(A) values.

  \emph{Choosing a method.} Based on the simulation results described below:
  \tabular{ll}{
    \strong{Quantity} \tab \strong{Recommended method} \cr
    Single GORIC(C) value \tab analytic casewise (\code{"delta"}); single values \cr
    \tab are rarely of interest, only their differences \cr
    Single GORICA(C) value \tab \code{"boot"} (the analytic CI is degenerate for \cr
    \tab hypotheses in agreement with the estimates) \cr
    Difference / ratio of weights: \tab \code{"invert"} (coverage .98 to 1.00) \cr
    \eqn{H_m} with only equalities or one inequality \tab \cr
    versus the unconstrained hypothesis or its complement \tab \cr
    Difference / ratio of weights: other pairs \tab \code{"boot"} (coverage about .91 to .99) \cr
    Weights, two hypotheses \tab as for their difference: \code{"invert"} if \cr
    \tab eligible (coverage .96 to .98), otherwise \code{"boot"} \cr
    Weights, more than two hypotheses \tab \code{"boot"} (coverage about .93 to .97) \cr
  }
  Hence, in general, \code{method = "invert"} is recommended: it gives the
  inversion CIs where these apply and the bootstrap CIs otherwise. The
  analytic CIs (\code{"delta"}, \code{"mc"}) need no simulation and are fine for
  clearly non-overlapping hypotheses in large samples, but they are too
  narrow when hypotheses are close or overlap, which is common for
  order-restricted hypotheses (coverage for differences about .83 to .92,
  for weights as low as .29 to .89).

  \emph{Analytic CIs.} Let \eqn{\Sigma} be the covariance matrix of
  \eqn{(\ell_1, \ldots, \ell_M)}.
  \itemize{
    \item For the GORIC(C) (\code{approach = "casewise"}),
    \eqn{\ell_m = \sum_i \ell_{im}} is a sum of casewise contributions, and
    \eqn{\Sigma_{mk}} is estimated by \eqn{\sum_i (\ell_{im} - \bar\ell_m)(\ell_{ik} - \bar\ell_k)},
    that is, \eqn{n} times the sample covariance of the casewise contributions
    (cf. Vuong, 1989). This gives
    \eqn{SE(IC_m) = 2 \sqrt{n}\, s_m} (cf. Vehtari, Gelman, & Gabry, 2017, Eq. 23) and
    \eqn{SE(IC_m - IC_k) = 2 \sqrt{n}\, \hat\omega_{mk}} (Merkle, You, & Preacher, 2016,
    Eq. 16; Vehtari et al., 2017, Eq. 24), with \eqn{s_m} the
    standard deviation of the casewise contributions and \eqn{\hat\omega_{mk}}
    that of their differences. Merkle et al. (2016) first apply Vuong's (1989)
    formal test of distinguishability; that test assumes estimates that are
    not on the boundary of the parameter space, which does not hold for
    order-restricted maximum likelihood estimates (or-mle's), so here pairs are only flagged when their
    standard error is (numerically) zero (see below). Casewise
    contributions are available for \code{lm} (also with prior weights),
    \code{rlm}, \code{mlm}, and \code{glm} objects (gaussian, poisson, binomial,
    Gamma, and inverse.gaussian family), and for \code{lm} objects with
    \code{missing = "fiml"} (observed-data contributions).
    \item For the GORICA(C) (\code{approach = "estimates"}),
    \eqn{-2\ell_m = \log|2\pi V| + D_m}, with
    \eqn{D_m = (\hat\theta - \tilde\theta_m)' V^{-1} (\hat\theta - \tilde\theta_m)}
    the squared Mahalanobis distance from the estimates to hypothesis \eqn{m}.
    The delta method gives \eqn{\Sigma_{mk} = r_m' V^{-1} r_k} with
    \eqn{r_m = \hat\theta - \tilde\theta_m}, so that
    \eqn{SE(IC_m - IC_k) = 2 \sqrt{(\tilde\theta_m - \tilde\theta_k)' V^{-1} (\tilde\theta_m - \tilde\theta_k)}}.
  }
  Then, the CI for \eqn{IC_m - IC_k} is \eqn{IC_m - IC_k \pm z\, SE(IC_m - IC_k)},
  the CI for the ratio of weights \eqn{w_m / w_k} is its exponential
  transformation \eqn{\exp(\eta_m - \eta_k \pm z\, SE(IC_m - IC_k)/2)}, and the
  CI for a weight is obtained on the logit scale via the delta method
  (\code{method = "delta"}) or via Monte Carlo draws (\code{method = "mc"}).

  \emph{Inversion} (\code{method = "invert"}; cf. the noncentral chi-square based
  CIs for the expected BIC and differences therein in Preacher & Merkle, 2012).
  For a hypothesis \eqn{H_m} versus the unconstrained hypothesis or versus its
  complement, the difference \eqn{T = IC_m - IC_k} depends on the true parameters
  via their signed distance \eqn{s} to the boundary of \eqn{H_m} (in the
  \eqn{V^{-1}} metric), when \eqn{H_m} consists of only equality constraints
  (then via the noncentrality \eqn{s^2}; for the GORICA, \eqn{T} then is a
  shifted noncentral chi-square variable) or of a single inequality constraint.
  On the line through the estimates perpendicular to the boundary of \eqn{H_m},
  the set of \eqn{s} for which the observed \eqn{T} lies in between the
  \eqn{(1 \pm level)/2} quantiles of its (simulated) distribution is determined,
  and the CI for the expected difference \eqn{E[T]} (and, via
  \eqn{\exp(\log(p_m/p_k) - E[T]/2)}, for the ratio of weights) follows. For
  only equality constraints, this is exactly the noncentrality-parameter
  inversion. For hypotheses with two or more inequality constraints (e.g.,
  \eqn{\mu_1 < \mu_2 < \mu_3}), the distribution also depends on the location
  along the boundary and this reduction does not hold (in a simulation, the
  coverage was then very poor), so the bootstrap CI is used for those pairs;
  the column \code{CI_method} shows which CI is given. Note that this CI
  concerns the \emph{expected} difference, so the sample difference can lie
  outside it. In a small simulation (three groups, \eqn{H: \mu_1 < \mu_2}
  versus the unconstrained hypothesis, \eqn{n} = 60 or 150), the coverage of
  95\% inversion CIs was .98 to 1.00 in all settings, whereas the bootstrap
  CIs ranged from .58 to 1.00 (the bootstrap undercovers when the estimates
  are far inside such a hypothesis, because its distribution is then
  degenerate).

  \emph{Hypotheses with the same or-mle's.} The analytic CIs require that the
  hypotheses do not have the same or-mle's in the sample (i.e., are distinguishable, in
  Vuong's terms). This can only occur for overlapping hypotheses, when the estimates
  lie in their common part; overlapping hypotheses can also have
  different or-mle's. When two hypotheses lead to the same order-restricted maximum likelihood estimates
  (or-mle's; e.g., a hypothesis and the unconstrained hypothesis
  when the estimates are in agreement with that hypothesis), the standard
  error of their difference is zero, while the difference in GORIC(A) values
  is in fact not fixed (it follows a non-normal, chi-bar-square-type
  distribution). Such pairs are flagged (\code{same_ormle = TRUE}) and get
  no analytic CI. Also when the true parameters are close to the boundaries of
  the hypotheses (relative to the sample size), the analytic CIs tend to be
  too narrow. In these cases, use \code{method = "boot"}.

  \emph{Parametric bootstrap} (\code{method = "boot"}). Estimates
  \eqn{\theta^*_b \sim N(\hat\theta, V)}, \eqn{b = 1, \ldots, B}, are drawn,
  with \eqn{V} the covariance matrix of the estimates in \code{object} (these
  are the same draws as those for the observed population in
  \code{\link{benchmark_asymp}}). For each draw, the squared Mahalanobis
  distance to each hypothesis is computed (for the complement: to the closest
  boundary of the hypothesis), giving GORICA-type values
  \eqn{IC^*_{mb} = \log|2\pi V| + D^*_{mb} + c\, PT_m}. The percentile CIs are
  the corresponding quantiles of the differences, ratios of weights, and
  weights over the draws. For a GORIC(C) object, the bootstrap uses this
  GORICA-type log-likelihood (asymptotically equivalent for differences), and
  the CIs for the single GORIC(C) values are the analytic (casewise) ones.

  In a small simulation study (three groups, \eqn{H_1: \mu_1 < \mu_2 < \mu_3}
  versus \eqn{H_2: \mu_1 < \mu_2}, true means (0, .5, .2), \eqn{n} = 60 or 150),
  the coverage of 95\% CIs for the difference in GORIC(A) values was about .83
  to .86 for the analytic CIs (approaching .95 only for large samples or
  clearly separated hypotheses) and about .91 to .93 for the bootstrap CIs.
  For the weights of \eqn{H_1}, \eqn{H_2}, and the unconstrained hypothesis
  (true means (0, .5, .2), (0, .3, .6), and (0, 1, .4), \eqn{n} = 60 or 150),
  the coverage was .29 to .89 for the analytic CIs and .93 to .97 for the
  bootstrap CIs, except for \eqn{H_2} when \eqn{H_1} (a subset of \eqn{H_2}) is
  true (.73). With two hypotheses (\eqn{H: \mu_1 < \mu_2} versus the
  unconstrained hypothesis or versus its complement, \eqn{n} = 60), the
  coverage for the weight of \eqn{H} was .05 to .82 for \code{"delta"}, .61 to
  .95 for \code{"boot"}, and .96 to .98 for \code{"invert"}.

  The CIs concern the sampling variability of the GORIC(A) values (and their
  differences, weights, and ratios) as estimates of their population
  counterparts; they are not posterior intervals. The penalties (and the
  chi-bar-square weights / level probabilities on which they are based) are
  treated as known.
}
\value{
  An object of class \code{goric_CI}, a list with, among others:
  \item{IC}{a data frame with, per hypothesis, the GORIC(A) value, its
  standard error (\code{SE}), and the lower and upper bound of the CI.}
  \item{difference}{a data frame with, for all ordered pairs of hypotheses
  (\code{model} versus \code{vs}), the difference in GORIC(A) values, its
  standard error, the CI, whether the hypotheses have the same or-mle's in the sample (\code{same_ormle}), and
  the method used for the CI (\code{CI_method}).}
  \item{weights}{a data frame with the GORIC(A) weights, their CIs, and the
  method used for the CI (\code{CI_method}).}
  \item{ratio_weights}{a data frame with, for all ordered pairs of hypotheses,
  the ratio of GORIC(A) weights and its CI.}
  \item{Sigma}{the estimated covariance matrix of the log-likelihood values.}
  \item{best}{the name of the hypothesis with the highest GORIC(A) weight.}
  \item{any_same_ormle}{whether any pair of hypotheses has the same or-mle's in this
  sample.}
  Also returned are \code{type}, \code{level}, \code{method}, \code{approach},
  \code{n} (the number of cases used in the casewise approach), \code{B},
  \code{penalty_factor}, and \code{priorICweights}.
}
\references{
  Merkle, E. C., You, D., & Preacher, K. J. (2016). Testing nonnested
  structural equation models. \emph{Psychological Methods, 21}(2), 151--163.
  \doi{10.1037/met0000038}

  Preacher, K. J., & Merkle, E. C. (2012). The problem of model selection
  uncertainty in structural equation modeling. \emph{Psychological Methods,
  17}(1), 1--14. \doi{10.1037/a0026804}

  Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
  evaluation using leave-one-out cross-validation and WAIC. \emph{Statistics
  and Computing, 27}(5), 1413--1432. \doi{10.1007/s11222-016-9696-4}

  Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
  non-nested hypotheses. \emph{Econometrica, 57}(2), 307--333.
  \doi{10.2307/1912557}
}
\author{Rebecca Kuiper}
\seealso{
  \code{\link{goric}}, \code{\link{benchmark}}
}
\examples{
## Example: three groups (lm object)
set.seed(123)
n <- 90
g <- factor(rep(c("a", "b", "c"), each = n / 3))
y <- rnorm(n, mean = c(0, 0.5, 0.2)[g])
fit <- lm(y ~ -1 + g)

out <- goric(fit, hypotheses = list(H1 = "ga < gb < gc", H2 = "ga < gb"))

## Analytic CIs (casewise approach, since type = "goric")
goric_CI(out)

## Parametric-bootstrap CIs
\donttest{
goric_CI(out, method = "boot", B = 1000, seed = 42)

## Inversion CIs (here for H2 versus the unconstrained hypothesis;
## bootstrap CIs for the other pairs)
goric_CI(out, method = "invert", B = 1000, seed = 42)
}

## GORICA based on estimates and their covariance matrix
est  <- coef(fit)
VCOV <- vcov(fit)
out_a <- goric(est, VCOV = VCOV,
               hypotheses = list(H1 = "ga < gb < gc", H2 = "ga < gb"))
ci_a <- goric_CI(out_a)
print(ci_a, all_pairs = TRUE)
ci_a$weights
}
