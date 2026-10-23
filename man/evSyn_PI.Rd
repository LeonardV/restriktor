\name{evSyn_PI}
\alias{evSyn_PI}
\alias{print.evSyn_PI}
\title{Prediction interval for the ratio of GORICA weights in a new study}
\description{
Computes a prediction interval (PI) for the ratio of GORICA weights of the overall best
hypothesis versus each other hypothesis in a new (comparable) study, based on an evidence
synthesis (\code{\link{evSyn}}). Where the confidence interval of \code{\link{evSyn_CI}}
regards the aggregated evidence, the prediction interval regards the evidence to expect in
a next study; it is the frequentist counterpart of the Bayes factor for the next person
(Klaassen et al., 2018), and summarizes the heterogeneity of the studies in an interpretable way.
For the average approach (\code{type_ev = "average"}), the interval regards another analysis
of the same data.
}
\usage{
evSyn_PI(object, vs = NULL, best = NULL, level = 0.95, n_new = NULL,
         type = c("observed", "true"), B = 1000L, seed = NULL, tol = 1e-8)

\method{print}{evSyn_PI}(x, digits = max(3, getOption("digits") - 4), \dots)
}
\arguments{
  \item{object}{an object of class \code{evSyn} from \code{\link{evSyn}}; for the added and equal
  approach with (standardized) estimates and their covariance matrices or goric objects as input
  (input type \code{"est_vcov"}, \code{"escalc"}, or \code{"gorica"}), and at least three studies
  (or clusters).}
  \item{vs}{the hypotheses to compare the best hypothesis with; by default, all other hypotheses.}
  \item{best}{the hypothesis of interest; by default, the one with the highest final GORICA weight.}
  \item{level}{the level of the prediction interval; by default, 0.95.}
  \item{n_new}{the sample size of the new study; by default, the median sample size of the studies.
  Only used when the sample sizes of the studies are known.}
  \item{type}{\code{"observed"} (default): the interval for the ratio observed in a new study (which
  includes its sampling variability); \code{"true"}: the interval for the true (population) ratio in a
  new study.}
  \item{B}{the number of bootstrap draws (per study) used for the study-specific standard errors.}
  \item{seed}{optional seed for the random number generator.}
  \item{tol}{tolerance used to determine whether two hypotheses have the same maximum
  order-restricted log likelihood in a study.}
  \item{x}{an object of class \code{evSyn_PI}.}
  \item{digits}{the number of decimals to print.}
  \item{\dots}{not used.}
}
\details{
The log of the ratio of GORICA weights of the best hypothesis (b) versus another one (k) in
study \eqn{s} equals \eqn{d_s + p_s}, with \eqn{d_s = LL_{b,s} - LL_{k,s}} the difference in the
maximum order-restricted log likelihoods (i.e., loglik in \code{\link{goric}} and LL in
\code{\link{evSyn}}) and \eqn{p_s = -(c/2)(PT_{b,s} - PT_{k,s})} the penalty part (with \eqn{c} the
penalty factor). The log-likelihood difference grows (roughly linearly) with the sample size, the
penalty part does not. Therefore, the per-observation difference \eqn{z_s = d_s / n_s} is modelled
with a random-effects model, \eqn{z_s \sim N(\mu, v_s + \tau^2)}, with \eqn{v_s} based on the
parametric-bootstrap standard error of \eqn{d_s} (as in \code{\link{goric_CI}}) and \eqn{\tau^2}
estimated with REML. For a new study of size \eqn{n_{new}}, the prediction interval of the log ratio
is \eqn{n_{new} (\hat\mu \pm t_{k-2} \sqrt{SE(\hat\mu)^2 + \hat\tau^2 + v_{new}}) + p_{new}}
(cf. Higgins, Thompson, & Spiegelhalter, 2009), with \eqn{v_{new}} the sampling variance in the new
study (omitted for \code{type = "true"}), \eqn{p_{new}} the median penalty part, and \eqn{k} the number
of studies used (at least 1 degree of freedom). When the sample sizes are unknown (specify
\code{study_sample_nobs} in \code{\link{evSyn}} or \code{sample_nobs} in \code{\link{goric}}), the
model is fitted to \eqn{d_s} itself, which assumes studies of comparable size.

Ceiling: when the best hypothesis and the other one have the same maximum order-restricted log
likelihood in a study (\eqn{d_s = 0}; the hypotheses overlap there and, when nested, the ratio is at its
bound), \eqn{d_s} is at its bound and not normally distributed. Then, a two-part model is used: the
probability that a new study is at the ceiling is estimated by the proportion of studies at the
ceiling, and the random-effects model is fitted to the other studies and truncated at the bound. The
prediction interval consists of the quantiles of this mixture (and is thus capped at the bound), and the
probability that a new study is at the ceiling is reported (\code{P(ceiling)}).

Heterogeneity measures (for the added and equal approach; based on the same random-effects model, and thus 
on the studies not at the ceiling): \eqn{\tau}, the between-study standard deviation of the per-observation 
evidence; \eqn{n_{med} \tau}, the same on the scale of the log ratio of weights for a study of median sample 
size \eqn{n_{med}} (and the corresponding factor \eqn{\exp(1.96\, n_{med} \tau)} by which the true ratios of 
such studies differ from the typical ratio); \eqn{I^2 = \tau^2 / (\tau^2 + v_{med})}, the proportion of the 
variance in the study-specific evidence that is due to true differences between the studies rather than 
sampling error, with \eqn{v_{med}} the typical within-study variance for a study of median sample size; and 
Cochran's \eqn{Q} test of homogeneity (with \eqn{k - 1} degrees of freedom; with few studies, it has little 
power). Since the evidence grows with the sample size, these are computed on the per-observation scale, on 
which the sample size cancels out; the median sample size is only used to express \eqn{\tau} on the ratio 
scale and to compute \eqn{I^2}. These measures regard heterogeneity in the evidence (which combines direction 
and size of the effects); for heterogeneity in the parameters themselves, a (multivariate) meta-analysis of 
the estimates (e.g., with the \pkg{metafor} package) is the appropriate tool.

Since the measures above are based on the studies not at the ceiling only, they can understate the 
heterogeneity when many studies are at the ceiling. Therefore, heterogeneity measures that use all studies are 
reported as well, based on the \emph{signed distance} of the estimates to the best hypothesis: the 
(Mahalanobis) distance to its boundary, positive inside the hypothesis (the depth, i.e., the distance to its 
complement) and negative outside (minus the distance to the hypothesis), that is, 
\eqn{s = \mathrm{sign}(D) \sqrt{2 |D|}} with \eqn{D} the difference in log likelihood of the best hypothesis 
and its complement (for a hypothesis with only equality constraints, the unconstrained hypothesis). At the 
ceiling, the log-likelihood difference is 0 whether the estimates lie just inside or far inside the 
hypothesis, whereas the signed distance distinguishes these. It is continuous, without a bound, and roughly 
normally distributed (its standard error is obtained with the parametric bootstrap); \eqn{s / \sqrt{n}} is a 
standardized margin per observation, on which the sample size cancels out. On this scale, \eqn{\tau}, 
\eqn{\sqrt{n_{med}} \tau} (the between-study SD of the signed distance for a study of median sample size, in 
standard-error units), \eqn{I^2} (for a study of median sample size), and Cochran's \eqn{Q} are computed. 
These regard how clearly the estimates lie inside (or outside) the best hypothesis, not the ratio of weights 
versus a specific other hypothesis. They are not available when the best hypothesis is the unconstrained 
hypothesis.

Also reported is the probability that a new study favours the best hypothesis (i.e., a ratio larger than
1), which relates to the evidence rate (ER) of the measures of variability (see \code{evSyn_variability}
in \code{\link{evSyn}}).

As the study-specific weights in \code{\link{evSyn}}, the study-specific ratios (and thus the prediction
interval) do not include the prior IC weights (\code{priorICweights}).

For clusters of dependent studies (\code{cluster} in \code{\link{evSyn}}), the clusters are the units:
the cluster-specific evidence is the (weighted) mean over its studies, which are treated as fully
dependent (the standard error of the mean is the weighted mean of the standard errors).

For the average approach, the analyses are based on the same data; the interval then regards another
analysis of the same data: \eqn{\bar y \pm t_{S-1} \sqrt{(1 + 1/S) B}}, with \eqn{y_s} the
analysis-specific log ratios and \eqn{B} their (between-analysis) variance. The sampling variance is
shared by all analyses and is not part of this interval.

With few studies, the heterogeneity (\eqn{\tau}) is estimated imprecisely, and so is the prediction
interval. In a simulation study (available in the folder \code{experiments/coverage_PI} of the source
code; 5, 10, or 20 studies of varying size, without and with heterogeneity, and with and without studies at
the ceiling), the coverage of the 95\% prediction interval for the observed ratio in a new study ranged from
.90 to 1.00 (somewhat conservative without heterogeneity). The probability that a new study favours the best
hypothesis is an approximation: it deviated up to about .15 from the observed proportions (too low for a
hypothesis versus its complement and too high with studies at the ceiling), partly because the best
hypothesis is selected based on the same studies.
}
\value{
An object of class \code{evSyn_PI}, a list with, among others, \code{PI}: a data frame with, per
comparison, the typical ratio (\code{ratio_typical}; the median of the predictive distribution), the
bounds of the prediction interval (\code{lower}, \code{upper}), the probability that a new study favours
the best hypothesis (\code{P_favours_best}), the number of studies at the ceiling and the probability
that a new study is at the ceiling (\code{n_ceiling}, \code{P_ceiling}), the estimated heterogeneity
(\code{tau}, on the per-observation scale), \code{n_new}, the number of studies used (\code{k}), the
degrees of freedom (\code{df}), and the method; and \code{heterogeneity}: a data frame with, per 
comparison, \code{k}, \code{n_ceiling}, \code{tau}, \code{tau_logratio} (\eqn{n_{med} \tau}), \code{I2}, 
\code{Q}, \code{Q_df}, \code{Q_p}, and \code{n_med} (for the added and equal approach); and 
\code{heterogeneity_sd}: the corresponding measures based on the signed distance to the best hypothesis (all 
studies), with also the numbers of studies with estimates inside and outside it, the mean margin 
(\code{mean_margin}, per observation), and \code{tau_distance} (\eqn{\sqrt{n_{med}} \tau}).
}
\references{
Higgins, J. P. T., Thompson, S. G., & Spiegelhalter, D. J. (2009). A re-evaluation of random-effects
meta-analysis. \emph{Journal of the Royal Statistical Society: Series A, 172}, 137--159.

Klaassen, F., Zedelius, C. M., Veling, H., Aarts, H., & Hoijtink, H. (2018). All for one or some for all?
Evaluating informative hypotheses using multiple N = 1 studies. \emph{Behavior Research Methods, 50},
2276--2291.
}
\seealso{\code{\link{evSyn}}, \code{\link{evSyn_CI}}, \code{\link{evSyn_forest}}, \code{\link{goric_CI}}}
\author{Rebecca M. Kuiper}
\examples{
set.seed(1)
est <- lapply(1:6, function(i) c(b1 = 0.30, b2 = 0.15, b3 = 0) + rnorm(3, 0, 0.1))
vcov <- rep(list(diag(0.01, 3)), 6)
res <- evSyn(est, VCOV = vcov, hypotheses = list(H1 = "b1 > b2 > b3"),
             study_sample_nobs = rep(300, 6))
evSyn_PI(res, B = 200, seed = 1)
}
