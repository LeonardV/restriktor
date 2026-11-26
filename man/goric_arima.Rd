\name{goric.Arima}
\alias{goric.Arima}
\alias{goric_arima}
\alias{print.goric_tsme}
\title{GORIC(A) for (AR)(I)MA Models, incl. an AR(1) Model with Measurement Error}
\description{
  \code{goric.Arima} evaluates (order-restricted) hypotheses about the
  parameters of a model fitted with \code{\link[stats]{arima}}, using the
  GORICA or the (likelihood-based) GORIC. For an ARMA(1,1) model, which is equivalent to an AR(1) model with
  measurement error (ME), the parameters of the latter model are derived, such
  that hypotheses about the ME can be evaluated, for example, whether the ME is
  negligible or relevant. Default hypotheses are available for this, together
  with a sensitivity analysis for the bound that defines 'negligible'.
}
\usage{
\method{goric}{Arima}(object, ..., hypotheses = NULL, comparison = NULL,
      type = "gorica",
      default_hypotheses = c("relevance", "minimum", "zero", "sign"),
      delta = 0.1, sensitivity = FALSE,
      delta_values = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50),
      VCOV_method = c("delta", "bootstrap"),
      B = 1000L, seed = NULL, boot_trim = 0.025,
      y = NULL, xreg = NULL,
      sample_nobs = NULL, debug = FALSE)

\method{print}{goric_tsme}(x, digits = max(3, getOption("digits") - 4), \dots)
}
\arguments{
  \item{object}{an object of class \code{Arima}, i.e., output from
  \code{\link[stats]{arima}}.}
  \item{...}{further arguments passed to \code{\link{goric}} (e.g.,
  \code{posthoc}, \code{priorICweights}, \code{CI}).}
  \item{hypotheses}{a named list of hypotheses, specified as in
  \code{\link{goric}}, using the parameter names listed under 'Details'. The
  word \code{delta} in a hypothesis is replaced by the value of the argument
  \code{delta} (e.g., \code{"pme < delta"}). If \code{NULL} (default), the
  default hypotheses chosen via \code{default_hypotheses} are used; these are
  only available for an ARMA(1,1) model.}
  \item{comparison}{as in \code{\link{goric}}. If \code{NULL} and the default
  hypotheses are used, the comparison that belongs to that set is used (see
  'Details'); otherwise, the default of \code{\link{goric}}.}
  \item{type}{\code{"gorica"} (default), \code{"goricac"}, \code{"goric"} or
  \code{NULL}. The GORIC uses the order-restricted maximum log-likelihood of
  the arima model (see 'Details'); it needs the data (found via the call of
  \code{arima()}, or given via \code{y}) and is available for non-seasonal
  models without differencing and without fixed coefficients; when it is not
  available, the GORICA is used (with a message). With \code{NULL}, the GORIC is
  used when available and the GORICA otherwise.}
  \item{default_hypotheses}{which default set of hypotheses is used when
  \code{hypotheses = NULL}: \code{"relevance"} (default), \code{"minimum"},
  \code{"zero"} or \code{"sign"}; see 'Details'.}
  \item{delta}{a number > 0 that is inserted for the word \code{delta} in the
  hypotheses. In the default \code{"relevance"} set, it is the bound for the
  proportion of ME variance (\code{pme}) below which the ME is called
  negligible. Default 0.1: then ignoring the ME attenuates the
  autoregressive parameter by at most 10 percent (see 'Details').}
  \item{sensitivity}{logical; if \code{TRUE}, the hypotheses (that contain
  \code{delta}) are also evaluated for all values in \code{delta_values}.
  Default \code{FALSE}.}
  \item{delta_values}{the values of \code{delta} used in the sensitivity
  analysis.}
  \item{VCOV_method}{how the covariance matrix of the estimates is obtained:
  \code{"delta"} (default; \code{var.coef} of the arima object, with the delta
  method for the derived parameters) or \code{"bootstrap"} (a parametric
  bootstrap; see 'Details').}
  \item{B}{the number of bootstrap samples (for \code{VCOV_method =
  "bootstrap"}).}
  \item{seed}{a seed for the bootstrap (optional).}
  \item{boot_trim}{each parameter's bootstrap estimates are winsorized at the
  \code{boot_trim} and \code{1 - boot_trim} quantiles before the covariance
  matrix is computed (default 0.025; 0 = no winsorizing).}
  \item{y, xreg}{the data (series and, if used in \code{arima}, the
  regressors). Only needed for \code{type = "goric"}, \code{VCOV_method =
  "bootstrap"} and \code{\link{benchmark_tsme}}, and only when they cannot be
  retrieved from the call of \code{arima()} (e.g., when \code{goric()} is
  called in another environment). It is checked that the log-likelihood of the
  arima object is reproduced.}
  \item{sample_nobs}{not used: the number of observations is taken from the
  arima object (needed for \code{type = "goricac"}).}
  \item{debug}{as in \code{\link{goric}}.}
  \item{x}{an object of class \code{goric_tsme} (\code{goric_arima} or
  \code{goric_arme}).}
  \item{digits}{the number of significant digits to use when printing.}
}
\details{
  \bold{Parameters.} The coefficients of the arima object can be used by their
  names as given by \code{coef(object)}: \code{ar1}, \code{ar2}, ...,
  \code{ma1}, ..., \code{sar1}, ..., \code{sma1}, ..., \code{intercept} (the
  mean, when \code{include.mean = TRUE} and there is no differencing) and the
  names of the \code{xreg} columns. Additionally, \code{sigma2} (the innovation
  variance) can be used. Fixed coefficients (argument \code{fixed} of
  \code{arima}) cannot be used.

  \bold{AR(1) model with measurement error.} An AR(1) process observed with
  white-noise ME,
  \deqn{x_t = \phi x_{t-1} + \epsilon_t, \quad \epsilon_t \sim N(0, s2), \qquad
        y_t = \mu + x_t + u_t, \quad u_t \sim N(0, r),}
  is an ARMA(1,1) process, \eqn{(1 - \phi L)(y_t - \mu) = (1 + \theta L) a_t}
  with \eqn{Var(a_t) = \sigma^2_a}, where \eqn{\phi} = \code{ar1},
  \eqn{\theta} = \code{ma1} and \eqn{\sigma^2_a} = \code{sigma2}. For an
  ARMA(1,1) model (\code{arima(y, order = c(1, 0, 1))}), the following
  parameters are derived and can be used in the hypotheses:
  \itemize{
    \item \code{r} \eqn{= -\theta \sigma^2_a / \phi}: the ME variance;
    \item \code{s2} \eqn{= \sigma^2_a (\phi + \theta)(1 + \theta\phi) / \phi}:
          the innovation variance of the signal;
    \item \code{sx2} \eqn{=} \code{s2} \eqn{/ (1 - \phi^2)}: the variance of
          the signal;
    \item \code{gamma0} \eqn{=} \code{sx2 + r}: the variance of y;
    \item \code{lambda} \eqn{=} \code{sx2 / gamma0}: the reliability;
    \item \code{pme} \eqn{=} \code{r / gamma0} \eqn{= 1 -} \code{lambda}: the
          proportion of the variance of y that is ME variance.
  }
  These are unconstrained estimates: \code{r} and \code{pme} can be negative
  (when \eqn{\theta\phi > 0}), and \code{s2} can be negative (when
  \eqn{|\theta| > |\phi|}). This is needed to evaluate, for example, the
  complement of \code{r > 0}. Support for a negative ME variance indicates that
  the AR(1) model with ME does not describe the data well (e.g., omitted
  dynamics).

  When the model is ignoring the ME, the autoregressive parameter is
  attenuated: the lag-1 autocorrelation of y is \code{lambda} \eqn{\phi}
  \eqn{= (1 -} \code{pme}\eqn{) \phi}. Hence, \code{pme} < \code{delta} means
  that ignoring the ME biases \eqn{\phi} by less than \code{100 * delta}
  percent. The value of \code{delta} should be chosen before seeing the data,
  based on what is considered negligible in the field (e.g., 0.1).

  \bold{Default hypotheses} (for an ARMA(1,1) model):
  \itemize{
    \item \code{"relevance"} (default): \code{H_negl: -delta < pme < delta}
          (negligible ME), \code{H_ME: pme > delta} (relevant ME) and
          \code{H_neg: pme < -delta} (negative ME variance; misfit), with
          \code{comparison = "none"} (together, these cover all possible
          values). Since \code{H_negl} is a range and not a point, it can
          receive (very) strong support, unlike \code{r = 0}.
    \item \code{"minimum"}: \code{H_ME: pme > delta} (at least relevant ME)
          versus its complement (\code{pme < delta}: negligible or negative
          ME variance). This is the same as \code{r > delta * gamma0} versus
          its complement, that is, a minimum ME variance relative to the
          variance of y.
    \item \code{"zero"}: \code{H0: r = 0}, \code{Hpos: r > 0} and
          \code{Hneg: r < 0}, with \code{comparison = "none"}. Note that the
          support for \code{H0} is bounded, since it lies on the boundary of
          \code{Hpos} (and \code{Hneg}).
    \item \code{"sign"}: \code{Hpos: r > 0} versus its complement. Without
          ME, \code{Hpos} is preferred in about half of the samples, so this
          set cannot show that there is no (relevant) ME; use it with a cut-off
          (e.g., via \code{\link{benchmark_tsme}}).
  }

  \bold{GORIC versus GORICA.} The GORICA approximates the likelihood by a
  normal distribution for the estimates of the parameters in the hypotheses.
  The GORIC (\code{type = "goric"}) uses the log-likelihood of the arima model
  itself (computed with \code{\link[stats]{KalmanLike}}), maximized under the
  restrictions of each hypothesis (order-restricted maximum likelihood; for
  restrictions on derived parameters such as \code{pme}, a penalty method is
  used, started from several points). The penalties are the same as those of
  the GORICA. When the ME parameters are weakly identified, the likelihood is
  far from normal. In a simulation (AR(1) with and without ME, scored on the
  conclusion 'ME' versus 'no ME'), the GORIC was less often convincingly wrong
  when the autoregressive parameter was small, but also gave less convincing
  support to the truth; overall, the GORICA performed as well or slightly
  better (and is much faster), which is why it is the default.

  \bold{Covariance matrix.} With \code{VCOV_method = "delta"}, the covariance
  matrix of the coefficients is \code{var.coef} of the arima object (based on
  the observed information), with \eqn{Var(\sigma^2_a) = 2\sigma^4_a / n}
  (asymptotically independent of the coefficients); for the derived
  parameters, the delta method is used. With \code{VCOV_method =
  "bootstrap"}, \code{B} series are simulated from the fitted model, the model
  is refitted and the parameters are derived; the covariance matrix of these
  (winsorized, see \code{boot_trim}) estimates is used. The bootstrap needs
  the data and is available for non-seasonal models without differencing and
  without fixed coefficients. For a benchmark of the weights that does not
  assume normally distributed estimates, see \code{\link{benchmark_tsme}}.

  \bold{Identification.} The ME parameters are only identified via the
  autocorrelations of y, i.e., via the dynamics: when \eqn{\phi} is close to 0
  (or when \eqn{\phi + \theta} is close to 0, i.e., y is close to white noise),
  the signal and the ME are hardly separable. Their estimates are then far
  from normally distributed and the GORICA (which approximates the likelihood
  by a normal distribution) can be overconfident. A message is given when
  \eqn{|\hat\phi| / SE < 2} or \eqn{|\hat\phi + \hat\theta| / SE < 2}. Then
  consider \code{type = "goric"} and/or \code{\link{benchmark_tsme}}, more
  time points, or a shorter time interval between the measurements (which
  increases \eqn{\phi} for a continuous-time process). For AR(p) and VAR(1)
  models with ME, see \code{\link{arme}}.
}
\value{
  An object of class \code{goric_arima} (and \code{goric_tsme}), which inherits
  from the class of the \code{\link{goric}} output (\code{con_goric}), so
  that, e.g., \code{summary()} and \code{\link{benchmark}} can be used. For
  \code{type = "goric"}, the log-likelihoods in \code{result} are those of the
  arima model and \code{ormle} contains the order-restricted maximum likelihood
  estimates of the parameters in the hypotheses. It contains an extra element
  \code{tsme}, a list with:
  \item{estimates}{a data frame with all available parameters: estimate,
  standard error, and whether it is involved in the hypotheses.}
  \item{VCOV}{the covariance matrix of all available parameters (singular when
  derived parameters are included).}
  \item{VCOV_method}{the method used for the covariance matrix.}
  \item{bootstrap}{for \code{VCOV_method = "bootstrap"}: the number of
  (usable) bootstrap samples and the bootstrap estimates.}
  \item{me}{whether the model is an ARMA(1,1) model, i.e., whether the ME
  parameters are available.}
  \item{label}{a description of the model.}
  \item{hypotheses, hypotheses_used}{the hypotheses as specified and with the
  value of \code{delta} inserted.}
  \item{delta}{the value of \code{delta} (if used in the hypotheses).}
  \item{default_hypotheses}{the default set used (if any).}
  \item{identification}{a note on weak identification (if applicable).}
  \item{sensitivity}{for \code{sensitivity = TRUE}: a data frame with, per
  value of \code{delta}, the GORIC(A) weights and the preferred hypothesis.}
  \item{loglik_unrestricted}{for the GORIC: the maximum log-likelihood.}
  \item{fit, mi}{the arima object and the internal model description (used by
  \code{\link{benchmark_tsme}}).}
}
\seealso{\code{\link{goric}}, \code{\link{arme}}, \code{\link{benchmark_tsme}},
\code{\link{benchmark}},
\code{\link[stats]{arima}}}
\author{Rebecca M. Kuiper and Leonard Vanbrabant}
\examples{
## Simulate an AR(1) process with measurement error (reliability 0.7)
set.seed(123)
n      <- 300
phi    <- 0.6
lambda <- 0.7                                  # reliability
sx2    <- 1 / (1 - phi^2)                      # variance of the signal
r      <- sx2 * (1 - lambda) / lambda          # ME variance
x <- arima.sim(model = list(ar = phi), n = n)  # signal
y <- 5 + as.numeric(x) + rnorm(n, sd = sqrt(r))

## Fit the AR(1) model with ME as an ARMA(1,1) model
fit <- arima(y, order = c(1, 0, 1), method = "ML")
coef(fit)

## Default: is the ME negligible (pme < 0.1), relevant or negative? (GORICA)
out <- goric(fit)
out
out$tsme$estimates       # all parameters, incl. the derived ME parameters

## The likelihood-based GORIC (needs the data, found via the call of arima())
goric(fit, type = "goric")

## Another bound for 'negligible', and a sensitivity analysis for delta
out_s <- goric(fit, delta = 0.2, sensitivity = TRUE)
out_s$tsme$sensitivity

## Other default sets
goric(fit, default_hypotheses = "minimum")   # pme > delta versus its complement
goric(fit, default_hypotheses = "zero")
goric(fit, default_hypotheses = "sign")

## Own hypotheses (may contain 'delta'), e.g., about the reliability and
## the autoregressive parameter
H <- list(H1 = "lambda > 1 - delta; ar1 > 0",
          H2 = "lambda < 1 - delta; ar1 > 0")
goric(fit, hypotheses = H, delta = 0.2)

## The same 'relevance' hypotheses, but linear in r and gamma0
H_lin <- list(H_negl = "r - delta*gamma0 < 0; r + delta*gamma0 > 0",
              H_ME   = "r - delta*gamma0 > 0",
              H_neg  = "r + delta*gamma0 < 0")
goric(fit, hypotheses = H_lin, comparison = "none")

\donttest{
## Parametric-bootstrap covariance matrix (e.g., when ar1 is small)
goric(fit, VCOV_method = "bootstrap", B = 500, seed = 1)

## Benchmark: what would the weights look like without ME (pme = 0)?
benchmark_tsme(out, B = 100, seed = 1)             # parametric bootstrap
benchmark(out, model_type = "asymp", pop_est = c(pme = 0), seed = 1)  # normal
}

## Other arima models: hypotheses about the coefficients
fit_ar2 <- arima(lh, order = c(2, 0, 0))
goric(fit_ar2, hypotheses = list(H1 = "ar1 > 0; ar2 < 0"))

year  <- time(LakeHuron) - 1920
fit_x <- arima(LakeHuron, order = c(2, 0, 0), xreg = year)
names(coef(fit_x))
goric(fit_x, hypotheses = list(H1 = "year < 0; ar1 > 0; ar2 < 0"))
}
\keyword{models}
\keyword{ts}
