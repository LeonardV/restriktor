\name{arme}
\alias{arme}
\alias{goric.arme}
\alias{print.arme}
\alias{coef.arme}
\alias{vcov.arme}
\alias{logLik.arme}
\alias{nobs.arme}
\title{AR(p) and VAR(1) Models with Measurement Error, and their GORIC(A)}
\description{
  \code{arme} fits an AR(p) model (one variable) or a VAR(1) model (two or more
  variables) with measurement error (ME) by exact Gaussian maximum likelihood.
  The ME variances are not restricted to be positive, such that hypotheses
  like 'the ME variance is positive' can be evaluated against their complement
  with \code{goric()}. For one variable and p = 1, the model is equivalent to
  an ARMA(1,1) model; see \code{\link{goric.Arima}}.
}
\usage{
arme(y, p = 1L, start = NULL, control = list())

\method{goric}{arme}(object, ..., hypotheses = NULL, comparison = NULL,
      type = "gorica",
      default_hypotheses = c("relevance", "minimum", "zero", "sign"),
      delta = 0.1, sensitivity = FALSE,
      delta_values = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50),
      VCOV_method = c("delta", "bootstrap"),
      B = 1000L, seed = NULL, boot_trim = 0.025,
      sample_nobs = NULL, debug = FALSE)

\method{print}{arme}(x, digits = max(3, getOption("digits") - 4), \dots)
}
\arguments{
  \item{y}{a numeric vector or univariate time series (AR(p) model), or a
  numeric matrix or data frame with one column per variable (VAR(1) model). No
  missing values are allowed; the time points are assumed to be equally
  spaced.}
  \item{p}{the order of the autoregressive part (only \code{p = 1} for more
  than one variable).}
  \item{start}{optional starting values (in the order of the parameters, see
  'Details').}
  \item{control}{a list of control parameters passed to
  \code{\link[stats]{optim}} (method \code{"BFGS"}).}
  \item{object}{an object of class \code{arme}.}
  \item{x}{an object of class \code{arme}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{...}{further arguments passed to \code{\link{goric}}.}
  \item{hypotheses, comparison, type, default_hypotheses, delta, sensitivity,
  delta_values, VCOV_method, B, seed, boot_trim, sample_nobs, debug}{as in
  \code{\link{goric.Arima}}.}
}
\details{
  \bold{Model.} For one variable:
  \deqn{x_t = \phi_1 x_{t-1} + \dots + \phi_p x_{t-p} + \epsilon_t, \quad
        \epsilon_t \sim N(0, s2), \qquad y_t = \mu + x_t + u_t, \quad
        u_t \sim N(0, r).}
  For k variables:
  \deqn{x_t = \Phi x_{t-1} + \epsilon_t, \quad \epsilon_t \sim N(0, \Sigma),
        \qquad y_t = \mu + x_t + u_t, \quad u_t \sim N(0, diag(r_1, \dots, r_k)).}
  The log-likelihood is computed with the Kalman filter, started from the
  stationary distribution of the signal. The ME variances may be negative, as
  long as the implied covariance matrix of the data is positive definite and
  \code{pme >= -1} (i.e., the variance of the signal is at most twice the
  variance of y). Without the latter bound, a large negative ME variance
  combined with a large signal variance can mimic other dynamics and give
  degenerate solutions, especially for a VAR(1) model. The
  covariance matrix of the estimates is the inverse of the (numerical) Hessian.

  \bold{Parameters} that can be used in hypotheses (\code{goric()}):
  \itemize{
    \item one variable: \code{mu}, \code{phi} (p = 1) or \code{phi1}, ...,
          \code{phip}, \code{s2}, \code{r}; and the derived \code{sx2} (variance
          of the signal), \code{gamma0} (variance of y), \code{lambda}
          (reliability, \code{sx2 / gamma0}) and \code{pme} (proportion of ME
          variance, \code{r / gamma0}).
    \item k variables: \code{mu1}, ..., \code{Phi11}, \code{Phi12}, ...
          (\code{Phiij}: the effect of variable j at t - 1 on variable i at t),
          \code{Sigma11}, \code{Sigma21}, \code{Sigma22}, ..., \code{r1}, ...;
          and the derived \code{sx2_1}, ..., \code{gamma0_1}, ...,
          \code{lambda1}, ... and \code{pme1}, ....
  }

  \bold{Default hypotheses} (see also \code{\link{goric.Arima}}). For one
  variable, the same sets as for an ARMA(1,1) arima object. For k variables:
  \itemize{
    \item \code{"relevance"}: per pattern of variables with relevant ME
          (\code{pmei > delta}) and negligible ME (\code{-delta < pmei < delta}),
          e.g., for two variables \code{H_negl}, \code{H_ME1}, \code{H_ME2} and
          \code{H_ME12}, with the unconstrained hypothesis as failsafe;
    \item \code{"minimum"}: \code{H_ME}: all \code{pmei > delta}, versus its
          complement;
    \item \code{"zero"}: \code{H0}: all \code{ri = 0} and \code{Hpos}: all
          \code{ri > 0}, with the unconstrained hypothesis as failsafe;
    \item \code{"sign"}: \code{Hpos}: all \code{ri > 0}, versus its complement.
  }

  The GORIC (\code{type = "goric"}) uses the order-restricted maximum
  log-likelihood of the model (see \code{\link{goric.Arima}}); the default is
  the GORICA. For a VAR(1) model, the GORIC takes some time (a constrained
  optimization per hypothesis).
}
\value{
  \code{arme} returns an object of class \code{arme}, a list with, a.o.,
  \code{coefficients}, \code{vcov}, \code{loglik}, \code{nobs}, \code{derived}
  (the derived parameters) and the data. \code{goric.arme} returns an object
  of class \code{goric_arme} (see \code{\link{goric.Arima}} for its
  elements).
}
\seealso{\code{\link{goric.Arima}}, \code{\link{benchmark_tsme}}, \code{\link{goric}}}
\author{Rebecca M. Kuiper and Leonard Vanbrabant}
\examples{
## One variable: AR(1) with ME (reliability 0.7)
set.seed(123)
n <- 300; phi <- 0.6; lambda <- 0.7
r <- (1 / (1 - phi^2)) * (1 - lambda) / lambda
y <- 5 + as.numeric(arima.sim(model = list(ar = phi), n = n)) + rnorm(n, sd = sqrt(r))

fit1 <- arme(y)
fit1
goric(fit1)                      # default 'relevance' set, GORICA
goric(fit1, type = "goric")      # likelihood-based GORIC

## AR(2) with ME, own hypotheses
fit2 <- arme(y, p = 2)
goric(fit2, hypotheses = list(H1 = "phi1 > 0; phi2 < 0; pme > delta"))

\donttest{
## Two variables: VAR(1) with ME in the first variable only
Phi   <- matrix(c(0.5, 0.2, -0.1, 0.4), 2, byrow = TRUE)
Sigma <- matrix(c(1, 0.3, 0.3, 1), 2)
n <- 300
x <- matrix(0, n + 100, 2)
for (t in 2:(n + 100)) x[t, ] <- Phi \%*\% x[t - 1, ] + t(chol(Sigma)) \%*\% rnorm(2)
x <- x[-(1:100), ]
Y <- x + cbind(rnorm(n, sd = 0.8), 0)
colnames(Y) <- c("stress", "mood")

fit_v <- arme(Y)
fit_v
goric(fit_v)                     # which variables have relevant ME?
goric(fit_v, hypotheses = list(H1 = "Phi12 > 0; Phi21 < 0"))
}
}
\keyword{models}
\keyword{ts}
