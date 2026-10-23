\name{evSyn_CI}
\alias{evSyn_CI}
\alias{print.evSyn_CI}
\title{Confidence intervals for the aggregated results of a GORICA evidence synthesis}
\description{
Computes confidence intervals (CIs) for the aggregated (i.e., final) GORICA values,
differences in GORICA values, GORICA weights, and ratios of GORICA weights of an
evidence synthesis (\code{\link{evSyn}}) based on (standardized) estimates and their
covariance matrices (input type \code{"est_vcov"} or \code{"escalc"}) or on goric objects
(input type \code{"gorica"}, i.e., a list of \code{\link{goric}} objects).
It extends \code{\link{goric_CI}} from one study to the aggregated results over multiple studies.

Notably, this regards the confidence interval of the aggregated evidence, not of pooled parameter 
estimates (as done in meta-analysis). Since aggregated evidence changes with the combined sample size, the 
aggregated evidence and its interval may be at another location than the study-/analysis-specific 
evidence and their intervals (it is not a weighted average, as it is for pooled parameter estimates; 
so, unlike in a meta-analysis, the aggregated evidence does not need to lie in between the 
study-specific evidence).
}
\usage{
evSyn_CI(object, level = 0.95, method = c("boot", "delta", "mc"),
         boot_type = c("auto", "percentile", "centered"),
         nsim = 10000L, B = 1000L, seed = NULL, tol = 1e-8, rho = 1,
         MI = isTRUE(object$MI))

\method{print}{evSyn_CI}(x, digits = max(3, getOption("digits") - 4),
      all_pairs = FALSE, \dots)
}
\arguments{
  \item{object}{an object of class \code{evSyn} from \code{\link{evSyn}}, with
  (standardized) estimates and their covariance matrices or goric objects as input (for 
  \code{dependence = "rve"} in \code{\link{evSyn}}, also log-likelihood and penalty values).}
  \item{level}{the confidence level; by default, 0.95.}
  \item{method}{\code{"boot"} (default): parametric-bootstrap CIs (see \code{boot_type}),
  based on \code{B} draws from N(estimates, VCOV) per study, independently over the
  studies; \code{"delta"}: normal-based CIs (delta method) based on the estimates and
  their covariance matrices, with the CIs for the weights on the logit scale;
  \code{"mc"}: as \code{"delta"}, but the CIs for the weights are obtained by Monte
  Carlo draws. See Details. The inversion method of \code{\link{goric_CI}}
  (\code{"invert"}) is not available for aggregated results; then, the bootstrap is used.}
  \item{boot_type}{the type of bootstrap CI (for \code{method = "boot"}): \code{"percentile"}, 
  the percentiles of the aggregated bootstrap draws; \code{"centered"}, the bootstrap draws are first 
  centered at the observed value (i.e., their mean shift is removed), after which the basic bootstrap 
  CI is used (on the scale of the GORICA values and their differences, of the log of the ratios of 
  weights, and of the logit of the weights); \code{"auto"} (default): \code{"centered"} from 5 studies 
  onwards and \code{"percentile"} otherwise. See Details.}
  \item{nsim}{the number of Monte Carlo draws for \code{method = "mc"}.}
  \item{B}{the number of bootstrap draws (per study) for \code{method = "boot"}.}
  \item{seed}{optional seed for the random number generator.}
  \item{tol}{tolerance used to determine whether the standard error of a difference is zero.}
  \item{rho}{for \code{type_ev = "average"} and for clusters of dependent studies 
  (\code{cluster} in \code{\link{evSyn}}, with \code{dependence = "cluster"}): the assumed correlation 
  (in between 0 and 1) between the sampling errors of the log-likelihood values of two analyses (or of 
  two studies within a cluster). The default, \code{rho = 1}, treats them as fully dependent, which is 
  the safe (conservative) choice. See Details for guidance; a sensitivity analysis over several values 
  (e.g., 0, .5, .8, and 1) is recommended. Not used with \code{MI = TRUE}.}
  \item{MI}{logical; only for \code{type_ev = "average"}: \code{TRUE} when the analyses are the analyses of 
  \strong{multiple imputed data sets} (one analysis per imputed data set); then, the CIs are based on Rubin's 
  rules (see Details). This is \strong{not} appropriate for different analyses of the same data (e.g., a 
  multiverse analysis, or other operationalizations or models), for which \code{MI = FALSE} (default) with 
  \code{rho} should be used. By default, the value used in \code{\link{evSyn}} (\code{MI = TRUE} there). 
  Alternatively, the MI estimates (i.e., the estimates pooled over the imputed data sets) and their 
  covariance matrix (the total covariance matrix from Rubin's rules) can be used as input to 
  \code{\link{goric}}.}
  \item{x}{an object of class \code{evSyn_CI}.}
  \item{digits}{the number of decimals to print.}
  \item{all_pairs}{logical; if \code{FALSE} (default), the differences and ratios are
  printed for the best hypothesis (i.e., the one with the highest aggregated weight) versus
  the others only; if \code{TRUE}, for all pairs of hypotheses.}
  \item{\dots}{not used.}
}
\details{
Since the studies are independent, the sampling (co)variance of the aggregated
log-likelihood values is the sum of the study-specific ones (obtained as in
\code{\link{goric_CI}}); the penalties and the prior IC weights (\code{priorICweights})
are fixed. As in \code{\link{evSyn}}:
\itemize{
  \item the aggregated GORICA values (and their differences) do not use the study weights;
    in the added-evidence approach, they are -2 x the sum of the log-likelihood values
    + penalty_factor x the sum of the penalty values; in the equal-evidence approach,
    -2 x the sum of the log-likelihood values + penalty_factor x the average of the penalty values;
  \item the aggregated GORICA weights (and their ratios) do use the study weights
    (\code{study_weights} in \code{\link{evSyn}}): the log-likelihood values (and penalties)
    are then weighted by the study weights (rescaled to sum to the number of studies).
    Hence, the study weights also affect the CIs of the weights and ratios of weights.
}
The parametric bootstrap is recommended when hypotheses have the same
maximum order-restricted log likelihood (i.e., LL in \code{\link{evSyn}}) in (some of) the studies
(e.g., a hypothesis and the unconstrained hypothesis when the estimates are in
agreement with that hypothesis; reported as \code{n_same_loglik}, the number of studies
in which this holds): the differences in GORICA values are then not normally
distributed and the normal-based CIs are too narrow. Notably, when a hypothesis is
contained in the other one (e.g., in the unconstrained hypothesis), the ratio of
weights has a maximum (determined by the difference in penalties), which then also
bounds the CI.

The bootstrap distribution of a difference in GORICA values is shifted relative to its sampling 
distribution (by roughly the number of active constraints per study). In one study, this shift is 
small relative to the spread, but in the aggregated results the shifts add up over the studies, whereas 
the spread grows only with the square root of the number of studies. In a simulation study (available 
in the folder \code{experiments/coverage_CI} of the source code), the percentile CIs had a coverage of 
about .95 or higher with 3 studies, but undercovered with 10 studies (down to about .75, mainly for 
the best hypothesis versus the failsafe hypothesis), which the centered CIs largely corrected (coverage 
.96 to 1.00 with 10 studies). With few studies, the centered CIs can be somewhat too narrow in some 
conditions (about .89 to .92 with 3 studies). Therefore, by default, the centered CIs are used from 5 
studies onwards. The delta method undercovers when hypotheses have the same maximum order-restricted 
log likelihood in (some of) the studies.

For goric objects of type GORIC(C) (e.g., based on a fitted \code{lm} object), the
covariance matrices are based on the casewise log-likelihood contributions and the bootstrap
uses the GORICA-type log-likelihood (as in \code{\link{goric_CI}}), with its draws centered
at the observed GORIC(C) log-likelihood values; the CIs of the aggregated GORIC(C) values are then
the analytic (casewise) ones.

\subsection{Average approach (\code{type_ev = "average"})}{
The analyses are then based on the same data and are thus dependent, so their sampling 
variability cannot simply be added (and they cannot be bootstrapped independently). Their 
covariances are unknown, and the estimates of the analyses need not be comparable; therefore, the 
correlation between the sampling errors of the log-likelihood values of two analyses is set to 
\code{rho}: the covariance matrix of the (weighted) mean of the log-likelihood values is then 
\deqn{(1 - \rho) \sum_s w_s^2 \Sigma_s + \rho (\sum_s w_s \Sigma_s^{1/2})^2,}
with \eqn{\Sigma_s} the covariance matrix in analysis \eqn{s} (as in \code{\link{goric_CI}}) and 
\eqn{w_s} the (rescaled) study weights. With \code{rho = 1} (default), the analyses are treated as 
fully dependent; this gives the maximal variance and thus conservative CIs. The CIs are based on the 
delta method (\code{method = "mc"} or \code{"delta"}). Besides the CIs, the between-analysis variability 
is reported (\code{SD_between}: the standard deviation over the analyses of the GORICA values, of their 
differences, and, \code{SD_between_log}, of the log ratios of weights). This variability reflects the 
analysis choices, not sampling variability, and is not part of the CIs. (The analyses are not regarded 
as multiple imputations, so Rubin's rules are not used.)

Guidance for \code{rho} (the correlation between the sampling errors, not the similarity of the 
results): analyses that differ only slightly (e.g., an extra covariate or another coding of the same 
variables) have strongly correlated sampling errors (about .8 to .95); analyses with different measures 
of the same constructs, typically about .5 to .8; analyses based on different (sub)sets of variables or 
partly different subsamples, lower (about .3 to .5). Hence, \code{rho = 1} is the worst case rather than 
the most likely value.

\emph{Multiple imputation} (\code{MI = TRUE}). When the analyses are the analyses of \eqn{m} multiple imputed 
data sets (and only then; not for a multiverse analysis), the aggregated GORICA values and weights are those of 
the average approach, but their (co)variance is obtained with Rubin's rules: the total covariance matrix is 
\eqn{T = \bar W + (1 + 1/m) B}, with \eqn{\bar W} the mean of the within-imputation covariance matrices of 
the log-likelihood values (as in \code{\link{goric_CI}}) and \eqn{B} the between-imputation covariance matrix 
(of -1/2 times the GORICA values). The CIs use t quantiles with Rubin's degrees of freedom, 
\eqn{\nu = (m - 1)(1 + 1/r)^2} with \eqn{r = (1 + 1/m) b / \bar w} for the quantity of interest (reported as 
\code{df}). Here, the between-imputation variability reflects the uncertainty due to the missing data and is 
therefore part of the CIs (whereas, for different analyses of the same data, the between-analysis variability 
reflects the analysis choices and is not part of the CIs). As elsewhere, the normal approximation is 
inaccurate when hypotheses have the same maximum order-restricted log likelihood in (some of) the imputed 
data sets. Alternatively, the MI estimates (i.e., the estimates pooled over the imputed data sets with 
Rubin's rules) and their covariance matrix (the total covariance matrix from Rubin's rules) can be used as 
input to \code{\link{goric}} (via \code{VCOV}); this evaluates the hypotheses once, on the pooled estimates.

A more exact alternative, which would not need \code{rho}, is a joint case bootstrap (re-running all 
analyses on the same bootstrap samples of the cases), which requires the data (e.g., fitted models on the 
same data frame); this is not (yet) available.
}

\subsection{Dependent studies (\code{cluster} in \code{evSyn})}{
For clusters of dependent studies (e.g., studies with overlapping samples or from the same lab or 
cohort), see \code{cluster} and \code{dependence} in \code{\link{evSyn}}:
\itemize{
  \item \code{dependence = "cluster"} (two-level evidence synthesis): the evidence is averaged within each 
    cluster and combined over the (independent) clusters; the covariance matrix of a cluster is obtained as 
    in the average approach above, with \code{rho} the correlation between the sampling errors of the 
    studies within a cluster. For overlapping samples, \code{rho} is roughly the proportion of shared 
    participants, \eqn{n_{shared} / \sqrt{n_1 n_2}}; for studies from the same lab or cohort but with 
    different participants, the sampling errors are (nearly) independent (\code{rho} about 0 to .2), and the 
    dependence mainly regards shared heterogeneity.
  \item \code{dependence = "rve"} (robust variance estimation; cf. Hedges, Tipton, & Johnson, 2010): the 
    aggregated evidence is that of the studies as such, and its covariance matrix is the cluster-robust one 
    (based on the residuals of the study-specific log-likelihood values around their (weighted) mean, 
    summed per cluster, with the small-sample factor C/(C - 1) and t(C - 1) quantiles, with C the number of 
    clusters). These CIs are robust to unknown dependence within the clusters, but also include the 
    heterogeneity of the studies, and require a sufficient number of clusters (about 10 or more).
}
In both cases, the bootstrap cannot be used (the studies are not independent), and \code{method = "mc"} 
is used instead.
}

The CIs are not available for the other input types
of \code{\link{evSyn}} (log-likelihood and penalty values, IC values, IC weights, or
ratios of IC weights), since these contain no information about the sampling variability
(for goric objects, the estimates and their covariance matrix are stored in the objects).
}
\value{
An object of class \code{evSyn_CI}, a list with, among others:
\item{IC}{a data frame with the aggregated GORICA values, their standard errors, and CIs.}
\item{difference}{a data frame with the differences in aggregated GORICA values for all
  ordered pairs of hypotheses, their standard errors, CIs, and \code{n_same_loglik}.}
\item{weights}{a data frame with the aggregated GORICA weights and their CIs.}
\item{ratio_weights}{a data frame with the ratios of aggregated GORICA weights for all
  ordered pairs of hypotheses, their CIs, and \code{n_same_loglik}.}
\item{Sigma_values, Sigma_eta}{the covariance matrices of the aggregated log-likelihood values,
  without and with study weights, respectively.}
\item{MI, df_MI}{whether Rubin's rules (multiple imputation) are used, and the corresponding degrees of 
  freedom (also given as the column \code{df} in the tables).}
\item{rho, dependence, n_clusters, df}{for the average approach and for dependent studies: the 
  correlation used, the type of dependence (\code{"cluster"} or \code{"rve"}), the number of clusters, 
  and the degrees of freedom of the t quantiles (\code{Inf}: normal quantiles).}
\item{best}{the hypothesis with the highest aggregated GORICA weight.}
}
\references{
Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance estimation in meta-regression with 
dependent effect size estimates. \emph{Research Synthesis Methods, 1}, 39--65.
}
\seealso{\code{\link{evSyn}} (with the argument \code{CI = TRUE}, the CIs are computed and
printed by \code{evSyn} itself), \code{\link{evSyn_forest}} (forest-type plot), 
\code{\link{evSyn_PI}} (prediction interval for a new study), \code{\link{goric_CI}}}
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

# Average approach (analyses of the same data): sensitivity to rho
res_av <- evSyn(est, VCOV = vcov, hypotheses = H, type_ev = "average")
evSyn_CI(res_av, rho = 1, seed = 1)
evSyn_CI(res_av, rho = 0.5, seed = 1)
}
