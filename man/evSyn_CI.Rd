\name{evSyn_CI}
\alias{evSyn_CI}
\alias{print.evSyn_CI}
\title{Confidence intervals for the pooled results of a GORICA evidence synthesis}
\description{
Computes confidence intervals (CIs) for the pooled (i.e., final) GORICA values,
differences in GORICA values, GORICA weights, and ratios of GORICA weights of an
evidence synthesis (\code{\link{evSyn}}) based on (standardized) estimates and their
covariance matrices (input type \code{"est_vcov"} or \code{"escalc"}) or on goric objects
(input type \code{"gorica"}, i.e., a list of \code{\link{goric}} objects).
It extends \code{\link{goric_CI}} from one study to the pooled results over multiple studies.
}
\usage{
evSyn_CI(object, level = 0.95, method = c("boot", "delta", "mc"),
         nsim = 10000L, B = 1000L, seed = NULL, tol = 1e-8)

\method{print}{evSyn_CI}(x, digits = max(3, getOption("digits") - 4),
      all_pairs = FALSE, \dots)
}
\arguments{
  \item{object}{an object of class \code{evSyn} from \code{\link{evSyn}}, with
  (standardized) estimates and their covariance matrices or goric objects as input and
  \code{type_ev = "added"} or \code{type_ev = "equal"}.}
  \item{level}{the confidence level; by default, 0.95.}
  \item{method}{\code{"boot"} (default): parametric-bootstrap percentile CIs,
  based on \code{B} draws from N(estimates, VCOV) per study, independently over the
  studies; \code{"delta"}: normal-based CIs (delta method) based on the estimates and
  their covariance matrices, with the CIs for the weights on the logit scale;
  \code{"mc"}: as \code{"delta"}, but the CIs for the weights are obtained by Monte
  Carlo draws. See Details. The inversion method of \code{\link{goric_CI}}
  (\code{"invert"}) is not available for pooled results; then, the bootstrap is used.}
  \item{nsim}{the number of Monte Carlo draws for \code{method = "mc"}.}
  \item{B}{the number of bootstrap draws (per study) for \code{method = "boot"}.}
  \item{seed}{optional seed for the random number generator.}
  \item{tol}{tolerance used to determine whether the standard error of a difference is zero.}
  \item{x}{an object of class \code{evSyn_CI}.}
  \item{digits}{the number of decimals to print.}
  \item{all_pairs}{logical; if \code{FALSE} (default), the differences and ratios are
  printed for the best hypothesis (i.e., the one with the highest pooled weight) versus
  the others only; if \code{TRUE}, for all pairs of hypotheses.}
  \item{\dots}{not used.}
}
\details{
Since the studies are independent, the sampling (co)variance of the pooled
log-likelihood values is the sum of the study-specific ones (obtained as in
\code{\link{goric_CI}}); the penalties and the prior weights (\code{priorICweights})
are fixed. As in \code{\link{evSyn}}:
\itemize{
  \item the pooled GORICA values (and their differences) do not use the study weights;
    in the added-evidence approach, they are -2 x the sum of the log-likelihood values
    + penalty_factor x the sum of the penalty values; in the equal-evidence approach,
    -2 x the sum of the log-likelihood values + penalty_factor x the average of the penalty values;
  \item the pooled GORICA weights (and their ratios) do use the study weights
    (\code{study_weights} in \code{\link{evSyn}}): the log-likelihood values (and penalties)
    are then weighted by the study weights (rescaled to sum to the number of studies).
    Hence, the study weights also affect the CIs of the weights and ratios of weights.
}
The parametric bootstrap is recommended when hypotheses have the same
order-restricted maximum likelihood estimates (or-mle's) in (some of) the studies
(e.g., a hypothesis and the unconstrained hypothesis when the estimates are in
agreement with that hypothesis; reported as \code{n_same_ormle}, the number of studies
in which this holds): the differences in GORICA values are then not normally
distributed and the normal-based CIs are too narrow. Notably, when a hypothesis is
contained in the other one (e.g., in the unconstrained hypothesis), the ratio of
weights has a maximum (determined by the difference in penalties), which then also
bounds the CI.

For goric objects of type GORIC(C) (e.g., based on a fitted \code{lm} object), the
covariance matrices are based on the casewise log-likelihood contributions and the bootstrap
uses the GORICA-type log-likelihood (as in \code{\link{goric_CI}}), with its draws centered
at the observed GORIC(C) log-likelihood values; the CIs of the pooled GORIC(C) values are then
the analytic (casewise) ones.

The CIs are not available for \code{type_ev = "average"}: the analyses are then based
on the same data and are thus dependent, so that their sampling variability cannot be
combined from the analysis-specific results (this would require re-sampling the data
and re-running all analyses). They are also not available for the other input types
of \code{\link{evSyn}} (log-likelihood and penalty values, IC values, IC weights, or
ratios of IC weights), since these contain no information about the sampling variability
(for goric objects, the estimates and their covariance matrix are stored in the objects).
}
\value{
An object of class \code{evSyn_CI}, a list with, among others:
\item{IC}{a data frame with the pooled GORICA values, their standard errors, and CIs.}
\item{difference}{a data frame with the differences in pooled GORICA values for all
  ordered pairs of hypotheses, their standard errors, CIs, and \code{n_same_ormle}.}
\item{weights}{a data frame with the pooled GORICA weights and their CIs.}
\item{ratio_weights}{a data frame with the ratios of pooled GORICA weights for all
  ordered pairs of hypotheses, their CIs, and \code{n_same_ormle}.}
\item{Sigma_values, Sigma_eta}{the covariance matrices of the pooled log-likelihood values,
  without and with study weights, respectively.}
\item{best}{the hypothesis with the highest pooled GORICA weight.}
}
\seealso{\code{\link{evSyn}} (with the argument \code{CI = TRUE}, the CIs are computed and
printed by \code{evSyn} itself), \code{\link{evSyn_forest}} (forest-type plot), \code{\link{goric_CI}}}
\author{Rebecca M. Kuiper}
\examples{
# Three studies with estimates of three (standardized) parameters
est <- list(c(b1 = 0.32, b2 = 0.12, b3 = 0.01),
            c(b1 = 0.25, b2 = 0.15, b3 = 0.04),
            c(b1 = 0.38, b2 = 0.09, b3 = -0.02))
vcov <- rep(list(diag(c(0.010, 0.012, 0.011))), 3)
H <- list(H1 = "b1 > b2 > b3", H2 = "b1 > b2 = b3")
res <- evSyn(est, VCOV = vcov, hypotheses = H)
evSyn_CI(res, B = 200, seed = 1)       # parametric bootstrap
evSyn_CI(res, method = "delta")       # delta method
}
