\name{evSyn_forest}
\alias{evSyn_forest}
\title{Forest-type plot of the study-specific and final ratios of GORICA weights}
\description{
Creates a forest-type plot for a GORICA evidence synthesis (\code{\link{evSyn}}) based on
(standardized) estimates and their covariance matrices (input type \code{"est_vcov"} or
\code{"escalc"}) or on goric objects (input type \code{"gorica"}): per study, the ratio of GORICA weights of the best hypothesis versus another
hypothesis with its confidence interval (CI; from \code{\link{goric_CI}}), and, at the bottom,
the final (pooled) ratio with its CI (from \code{\link{evSyn_CI}}).
}
\usage{
evSyn_forest(object, vs = NULL, best = NULL, level = 0.95,
             method_study = c("invert", "boot", "delta"),
             method_pooled = c("boot", "delta", "mc"),
             B = 1000L, seed = NULL, \dots)
}
\arguments{
  \item{object}{an object of class \code{evSyn} from \code{\link{evSyn}}, with
  (standardized) estimates and their covariance matrices or goric objects as input.}
  \item{vs}{character vector with the hypotheses that the best hypothesis is compared
  with (one panel per hypothesis); by default, all other hypotheses in the set.}
  \item{best}{the hypothesis in the numerator of the ratios; by default, the hypothesis
  with the highest final GORICA weight.}
  \item{level}{the confidence level; by default, 0.95.}
  \item{method_study}{the method used for the study-specific CIs, see \code{\link{goric_CI}};
  by default, its recommended method \code{"invert"}.}
  \item{method_pooled}{the method used for the CI of the final (pooled) ratio, see
  \code{\link{evSyn_CI}}; by default, \code{"boot"}. When the evSyn object already contains
  CIs (\code{evSyn(..., CI = TRUE)}) with the same level and method, these are used.}
  \item{B}{the number of bootstrap draws (per study).}
  \item{seed}{optional seed for the random number generator.}
  \item{\dots}{not used.}
}
\details{
The studies are shown in the same order and with the same labels as in the evSyn object
(see \code{order_studies} and \code{study_names} in \code{\link{evSyn}}). The x-axis has a
log scale, and a dashed line denotes a ratio of 1 (equal support). With unequal study weights,
the size of the study points reflects the study weights. The study-specific ratios include the
ratio of the prior weights (\code{priorICweights}) used in \code{evSyn}, as do the study-specific
weights in \code{evSyn} (for goric objects as input, the prior weights used in these objects are
replaced by those).

Open points denote studies in which the two hypotheses have the same fit (i.e., the same
order-restricted MLEs). If one of the two hypotheses is contained in the other (e.g., a
hypothesis in the unconstrained hypothesis), the ratio is then at its bound (determined by the
difference in penalties): its maximum if the best hypothesis is contained in the other one, and
its minimum otherwise; this bound also bounds the CI.

For \code{type_ev = "average"}, the final ratio is shown without a CI, since the analyses
(of the same data) are dependent (see \code{\link{evSyn_CI}}).
}
\value{
A \code{ggplot} object, which can be printed or adjusted further. The plotted data are
available as \code{attr(plot, "forest_data")}.
}
\seealso{\code{\link{evSyn}}, \code{\link{evSyn_CI}}, \code{\link{goric_CI}}}
\author{Rebecca M. Kuiper}
\examples{
est <- list(c(b1 = 0.32, b2 = 0.12, b3 = 0.01),
            c(b1 = 0.25, b2 = 0.15, b3 = 0.04),
            c(b1 = 0.38, b2 = 0.09, b3 = -0.02))
vcov <- rep(list(diag(c(0.010, 0.012, 0.011))), 3)
H <- list(H1 = "b1 > b2 > b3", H2 = "b1 > b2 = b3")
res <- evSyn(est, VCOV = vcov, hypotheses = H)
evSyn_forest(res, B = 200, seed = 1)
}
