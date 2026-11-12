\name{goric_posthoc}
\alias{goric_posthoc}
\alias{print.goric_posthoc}
\title{Post-hoc Analysis: the Best Hypothesis versus its Complement}
\description{
  \code{goric_posthoc} evaluates the best hypothesis of a GORIC(A) analysis
  (i.e., the hypothesis with the highest GORIC(A) weight in the set) against
  its complement. This gives the ratio of GORIC(A) weights of the best
  hypothesis versus its complement and the conditional error probability,
  that is, the probability (given the data) that choosing the best hypothesis
  is incorrect. The same analysis is performed by \code{goric()} (argument
  \code{posthoc}, default \code{TRUE}) when there is more than one specified
  hypothesis and one of these is the best.
}
\usage{
goric_posthoc(object, priorICweights = NULL, ...)

\method{print}{goric_posthoc}(x, digits = max(3, getOption("digits") - 4), \dots)
}
\arguments{
  \item{object}{an object of class \code{con_goric} (i.e., output from the
  \code{\link{goric}} function).}
  \item{priorICweights}{\code{NULL} (default) or a vector of two prior
  GORIC(A) weights, for the best hypothesis and its complement, respectively.
  By default, equal prior weights are used if the prior weights in the
  original set were equal; otherwise, the best hypothesis keeps its prior
  weight from the original set and its complement gets 1 minus that value.}
  \item{x}{an object of class \code{goric_posthoc}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{...}{currently not used.}
}
\details{
  The post-hoc analysis is useful when a set of hypotheses is evaluated (e.g.,
  several order-restricted hypotheses together with the unconstrained
  hypothesis as a safeguard) and one of the specified hypotheses is the best.
  The GORIC(A) weights then tell which hypothesis is the best in the set, but
  not how much support it has compared to all other possible orderings of the
  parameters. That is what the evaluation of the best hypothesis against its
  complement shows.

  The analysis is not performed (\code{$applicable = FALSE}, with the reason
  in \code{$reason}) when the best hypothesis is the unconstrained hypothesis
  (then none of the hypotheses is supported), a complement, or the
  equality-restricted hypothesis Heq; nor when the output already contains the
  (single) hypothesis versus its complement (\code{comparison = "complement"}).
  If the best hypothesis contains only equality restrictions, its complement
  is the unconstrained hypothesis.

  Since the best hypothesis is selected based on the data, the result is
  exploratory (post-hoc) and can be used to inform future research, but it is
  not a confirmatory test.
}
\value{
  An object of class \code{goric_posthoc}, a list with:
  \item{applicable}{whether the post-hoc analysis was performed.}
  \item{reason}{if not applicable, the reason why.}
  \item{best}{the name of the best hypothesis in the original set.}
  \item{hypothesis}{the best hypothesis (as specified), if available.}
  \item{comparison}{\code{"complement"} (or \code{"unconstrained"} for a
  hypothesis with only equality restrictions).}
  \item{result}{the result table of the best hypothesis versus its complement.}
  \item{ratio.gw}{the ratio of GORIC(A) weights of the best hypothesis
  versus its complement.}
  \item{error_prob}{the conditional error probability, that is, the GORIC(A)
  weight of the complement.}
  \item{priorICweights}{the prior weights used.}
  \item{goric}{the full \code{con_goric} object of the best hypothesis versus
  its complement.}
}
\seealso{\code{\link{goric}}}
\author{Rebecca M. Kuiper and Leonard Vanbrabant}
\examples{
est  <- c(x1 = 0.5, x2 = 0.3, x3 = 0.1)
VCOV <- diag(c(0.01, 0.01, 0.01))
fit <- goric(est, VCOV = VCOV,
             hypotheses = list(H1 = "x1 > x2 > x3", H2 = "x1 < x2 < x3"))
goric_posthoc(fit)

# this is also done by goric() by default (posthoc = TRUE):
fit$posthoc

# switch it off (e.g., in simulations):
goric(est, VCOV = VCOV,
      hypotheses = list(H1 = "x1 > x2 > x3", H2 = "x1 < x2 < x3"),
      posthoc = FALSE)
}
