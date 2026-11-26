\name{benchmark_tsme}
\alias{benchmark_tsme}
\alias{print.benchmark_tsme}
\title{Parametric-Bootstrap Benchmark for the GORIC(A) of Time-Series Models}
\description{
  \code{benchmark_tsme} benchmarks the GORIC(A) weights of a \code{goric()}
  analysis of an arima or arme object (e.g., about measurement error, ME). Data
  sets are simulated from the fitted model under a population restriction
  (by default: no ME), the model is refitted and the analysis is redone. The
  observed weights are then located in this bootstrap distribution. Unlike
  \code{\link{benchmark}}, this does not assume that the estimates are
  normally distributed, which matters when the ME parameters are weakly
  identified (e.g., a small autoregressive parameter).
}
\usage{
benchmark_tsme(object, pop = NULL, B = 500L, seed = NULL,
               probs = c(0.05, 0.35, 0.50, 0.65, 0.95))

\method{print}{benchmark_tsme}(x, digits = 3, \dots)
}
\arguments{
  \item{object}{output of \code{goric()} for an object of class \code{Arima}
  (\code{\link{goric.Arima}}) or \code{arme} (\code{\link{arme}}).}
  \item{pop}{the population restriction(s), in the syntax of the hypotheses,
  e.g., \code{"pme = 0"} or \code{"pme = 0.1"}. The model is fitted under these
  restrictions (order-restricted maximum likelihood) and the data are
  simulated from this fit. Default: no ME (\code{"pme = 0"}, or all
  \code{pmei = 0} for a VAR(1) model).}
  \item{B}{the number of bootstrap samples.}
  \item{seed}{a seed (optional).}
  \item{probs}{the quantiles of the bootstrap distribution that are reported.}
  \item{x}{an object of class \code{benchmark_tsme}.}
  \item{digits}{the number of digits to print.}
  \item{...}{not used.}
}
\details{
  Per bootstrap sample, the same hypotheses, comparison, type (GORIC or
  GORICA) and \code{delta} are used as in \code{object}; the covariance matrix
  is obtained via the delta method. The \code{percentile} is the percentage of
  bootstrap samples in which the weight is smaller than or equal to the
  observed weight. With the default population (no ME), a percentile of, say,
  97 for the hypothesis with ME means that such a high weight occurs in only 3
  percent of the samples without ME. This can be used to calibrate a decision
  rule (e.g., conclude 'relevant ME' only above the 95th percentile), which
  keeps the rate of false conclusions at about 5 percent also when the ME
  parameters are weakly identified, at the cost of power.

  The data are needed: for an arima object, they are taken from the call of
  \code{arima()}, or from the arguments \code{y} (and \code{xreg}) of
  \code{goric()}.
}
\value{
  An object of class \code{benchmark_tsme}, a list with the population
  restriction and estimates, the number of (usable) bootstrap samples, the
  observed weights, their percentiles, the quantiles of the bootstrap
  distribution, and the bootstrap weights (\code{draws}).
}
\seealso{\code{\link{goric.Arima}}, \code{\link{arme}}, \code{\link{benchmark}}}
\author{Rebecca M. Kuiper and Leonard Vanbrabant}
\examples{
\donttest{
set.seed(123)
n <- 200; phi <- 0.3; lambda <- 0.7
r <- (1 / (1 - phi^2)) * (1 - lambda) / lambda
y <- as.numeric(arima.sim(model = list(ar = phi), n = n)) + rnorm(n, sd = sqrt(r))
fit <- arima(y, order = c(1, 0, 1), method = "ML")

out <- goric(fit, default_hypotheses = "minimum")   # H_ME: pme > 0.1 vs complement
out$result
benchmark_tsme(out, B = 200, seed = 1)              # population: no ME (pme = 0)
benchmark_tsme(out, pop = "pme = 0.1", B = 200, seed = 1)
}
}
\keyword{models}
\keyword{ts}
