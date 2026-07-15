\name{iht_joint}
\alias{iht_joint}
\alias{pchibar_joint}
\alias{pfbar_joint}
\alias{print.conTestJoint}

\title{joint informative hypothesis test for the type A and type B statistics}
\description{\code{iht_joint} computes the type A and the type B test
statistic for the same set of order restrictions and evaluates them
\emph{jointly} under the shared (least favorable) null hypothesis
\eqn{R\theta = rhs}, using the joint distribution of the two statistics
derived by Wolak (1989, Theorem 1). In addition to the two marginal
p-values it reports the joint tail probability
\eqn{Pr[T_B \ge t_B, T_A \ge t_A]} and a combined decision rule.
}

\usage{
iht_joint(object, constraints = NULL, test = "F",
          rhs = NULL, neq = 0L, alpha = 0.05, \dots)

pchibar_joint(cIU, cEI, df_total, wt.bar)

pfbar_joint(cIU, cEI, df_total, df2, wt.bar)
}

\arguments{
  \item{object}{an object of class \code{lm}, \code{glm} or \code{rlm}.
  In this case, the constraint syntax needs to be specified

  OR

  an object of class \code{restriktor} (conLM, conGLM or conRLM). The
  constraints are inherited from the fitted restriktor object and do not
  need to be specified again.}

  \item{constraints}{constraint syntax (text-based or a matrix
  \eqn{R}); see \code{\link{restriktor}} for all options.}

  \item{test}{test statistic; \code{"F"} (default), \code{"LRT"} or
  \code{"score"} for lm/glm objects, \code{"F"}, \code{"Wald"} or
  \code{"score"} for rlm objects; identical to the \code{test} argument
  of \code{\link{iht}}.}

  \item{rhs}{vector on the right-hand side of the constraints. Only used
  if the constraints input is a matrix or vector; see \code{\link{iht}}.}

  \item{neq}{integer (default = 0): the first \code{neq} rows of the
  constraints matrix are treated as equality constraints. Only used if
  the constraints input is a matrix or vector; see \code{\link{iht}}.}

  \item{alpha}{significance level (default = 0.05) used by the combined
  decision rule; see details.}

  \item{cIU}{critical value / observed value of the type B
  (inequality-unrestricted, IU) statistic.}

  \item{cEI}{critical value / observed value of the type A
  (equality-inequality, EI) statistic.}

  \item{df_total}{total number of constraints \eqn{q} (equalities plus
  inequalities), i.e. \code{nrow(Amat)}.}

  \item{df2}{residual degrees of freedom of the unrestricted model
  (\code{df2 = Inf} gives the chi-bar-square version).}

  \item{wt.bar}{chi-bar-square weights in the type A convention as
  stored on a restriktor object fitted with
  \code{mix_weights = "pmvnorm"}: \code{wt.bar[k+1]} is the weight of
  \eqn{\chi^2_k} in the type A mixture, \eqn{k = 0, \dots, L} with
  \eqn{L} the number of inequality constraints.}

  \item{\dots}{further options for the \code{\link{restriktor}} function
  (e.g. \code{mix_weights}, \code{control}) when \code{object} is an
  unrestricted model. \code{neq.alt} and \code{boot} are not available
  for the joint test; see details.}
}

\details{
The type A test evaluates H0: \eqn{R\theta = rhs} against
HA: \eqn{R\theta \ge rhs} with test statistic
\eqn{T_A = (\tilde{\theta} - \bar{\theta})' \Sigma^{-1} (\tilde{\theta} - \bar{\theta})}
(EI, distance between the equality-restricted and the order-restricted
estimate), and the type B test evaluates H0: \eqn{R\theta \ge rhs}
against an unrestricted alternative with
\eqn{T_B = (\hat{\theta} - \tilde{\theta})' \Sigma^{-1} (\hat{\theta} - \tilde{\theta})}
(IU, distance between the unrestricted and the order-restricted
estimate). Both statistics are computed by the same machinery as
\code{\link{iht}}.

Under the least favorable null configuration \eqn{R\theta = rhs}, which
the two hypotheses share, Wolak (1989, Theorem 1) shows that

\deqn{Pr[T_B \ge c_1, T_A \ge c_2] = \sum_{k=0}^{L} w(L, k, \Omega)
      Pr[\chi^2_{q-k} \ge c_1] Pr[\chi^2_k \ge c_2],}

where \eqn{q} is the total number of constraints, \eqn{L = q - neq} the
number of inequality constraints, and \eqn{w(L, k, \Omega)} the usual
chi-bar-square (level probability) weights that restriktor computes for
the marginal p-values. Conditional on \eqn{k}, \eqn{T_A \sim \chi^2_k}
and \eqn{T_B \sim \chi^2_{q-k}} are independent, but marginally the two
statistics are dependent: the naive product of the two marginal
p-values/levels is \emph{not} the joint rejection probability.
\code{pchibar_joint} implements this formula.

For linear models the error variance is estimated and both statistics
share the same estimate \eqn{s^2}; the exact finite-sample analogue
replaces the chi-square distributions by a mixture involving
F-distributions (Wolak, 1987, Theorem 4.4). \code{pfbar_joint}
implements the exact joint version by integrating the conditional
chi-square product over the distribution of the shared variance
estimate; because the denominator is shared, the joint tail probability
is not the product of the two marginal F-bar tails (the product form is
recovered in the limit \code{df2 = Inf}). \code{iht_joint} uses
\code{pfbar_joint} with \code{df2} equal to the residual degrees of
freedom, consistent with the marginal p-values of \code{\link{iht}}.
For lm objects this is exact under normality; for glm and rlm objects
it holds asymptotically, as for the marginal tests.

The combined decision rule at level \code{alpha} distinguishes three
regions based on the marginal tests:
\itemize{
  \item (a) type A not rejected: H0: \eqn{R\theta = rhs} is retained;
  \item (b) type A rejected and type B not rejected:
  \eqn{R\theta \ge rhs} and \eqn{R\theta \ne rhs}, i.e. the
  order-constrained hypothesis is supported;
  \item (c) type B rejected (regardless of type A): H0:
  \eqn{R\theta \ge rhs} is rejected.
}
The joint tail probability \eqn{Pr[T_B \ge t_B, T_A \ge t_A]} quantifies
how exceptional the observed \emph{pair} of statistics is under the
shared null, which is the appropriate error measure when conclusions
depend on both tests simultaneously.

Restrictions: the joint distribution is derived under the full null
hypothesis, so \code{neq.alt} (partial alternatives for type B) is not
available; bootstrapped p-values (\code{boot}) are not available; the
restriktor object must carry chi-bar-square weights
(\code{mix_weights = "pmvnorm"} (default) or \code{"boot"}), and at
least one constraint must be an inequality.
}

\value{
  \code{iht_joint} returns an object of class \code{conTestJoint}, for
  which a print function is available, with among others the elements
  \code{Ts.A}, \code{Ts.B} (the two test statistics), \code{pvalue.A},
  \code{pvalue.B} (marginal p-values, identical to \code{\link{iht}}),
  \code{pvalue.joint} (the joint tail probability
  \eqn{Pr[T_B \ge Ts.B, T_A \ge Ts.A]}), \code{wt.bar},
  \code{b.eqrestr}, \code{b.restr}, \code{b.unrestr}, \code{alpha} and
  \code{decision} (\code{"a"}, \code{"b"} or \code{"c"}, see details).

  \code{pchibar_joint} and \code{pfbar_joint} return the joint upper
  tail probability (a single number in [0, 1]).

  The joint tail probability is also reported by
  \code{iht(..., type = "summary")} (see \code{\link{conTest_summary}}),
  next to the marginal type A and type B results.
}

\references{
  Kudo, A. (1963). A multivariate analogue of the one-sided test.
  \emph{Biometrika}, \bold{50}, 403--418.

  Silvapulle, M.J. and Sen, P.K. (2005). \emph{Constrained Statistical
  Inference}. Wiley, New York.

  Wolak, F. (1987). An exact test for multiple inequality and
  equality constraints in the linear regression model.
  \emph{Journal of the American Statistical Association},
  \bold{82}, 782--793.

  Wolak, F. (1989). Local and global testing of linear and nonlinear
  inequality constraints in nonlinear econometric models.
  \emph{Journal of Econometrics}, \bold{41}, 205--235.
}

\author{Leonard Vanbrabant}

\seealso{
  \code{\link{iht}},
  \code{\link{restriktor}}
}

\examples{
# the data consist of ages (in months) at which an
# infant starts to walk alone.
DATA <- subset(ZelazoKolb1972, Group != "Control")

# fit unrestricted linear model
fit.lm <- lm(Age ~ -1 + Group, data = DATA)

# joint test of the ordered hypothesis
myConstraints <- ' GroupActive < GroupPassive < GroupNo '
iht_joint(fit.lm, myConstraints)

# or based on a fitted restriktor object
fit.con <- restriktor(fit.lm, constraints = myConstraints)
iht_joint(fit.con)

# the joint tail probability directly:
out <- iht_joint(fit.con)
pfbar_joint(cIU = out$Ts.B, cEI = out$Ts.A,
            df_total = nrow(out$Amat), df2 = out$df.residual,
            wt.bar = out$wt.bar)
}
