\name{calculate_IC_weights}
\alias{calculate_IC_weights}
\alias{calc_ICweights}
\alias{IC_weights}
\alias{print.goric_ICw}
\title{Calculating IC weights based on IC values (AIC, ORIC, GORIC(A), BIC, SIC, ...)}

\usage{
calculate_IC_weights(IC, hypo_names = NULL)
calc_ICweights(IC, hypo_names = NULL)

\method{print}{goric_ICw}(x, digits = max(3, getOption("digits") - 4), use_scientific, \dots)

}
\arguments{
  \item{IC}{A vector or one-column matrix with information criteria (AIC, ORIC, GORIC(A), BIC, SIC,
            \ldots) values of length 'NrHypos', where 'NrHypos' stands for the number of hypotheses/
            models.}
  \item{x}{an object of class \code{con_goric}.}
  \item{hypo_names}{Optional. Vector containing 'NrHypos' characters which will be used for labeling the
                    hypothesis. Default: H1, H2, \ldots}
  \item{use_scientific}{If TRUE (default), the IC weights and ratio of IC weights 
  will be formatted using scientific notation. If FALSE, standard numeric formatting is used.}                  
  \item{digits}{the number of significant digits to use when printing.}
  \item{\ldots}{no additional arguments for now.}
}
\value{
IC weights, which quantify the relative strength of hypotheses/models.
}
\description{
This function transforms IC values into IC weights: IC values denote the ordering 
of hypotheses/models, while IC weights quantify the relative strength of hypotheses/models.
}

\examples{
IC <- c(1,2,3)
calculate_IC_weights(IC)


## PT weights 
# This examples shows how to calculate PT weights.
# Notably, one is interested in PT weights when the log-likelihood for two or more
# hypotheses are (approximately) equal.
# Then, the comparison between those hypotheses is solely based on the PT values.
# The IC weights will then equal the PT weights.
# In that case, there is support for the overlap (boundary) of these hypotheses.
# Thus, when the IC weights equal the PT weights for a (sub)set of hypotheses,
# then there is support for the overlap (boundary) of these hypotheses.

y <- rnorm(30)
group <- factor(rep(c("A","B","C"), each = 10))
fit.lm <- lm(y ~ -1 + group)
est <- coef(fit.lm)
VCOV_est <- vcov(fit.lm)
H1 <- "groupA < groupB < groupC"
results <- goric(est, VCOV = VCOV_est, hypotheses = list(H1), 
                 comparison = "complement", type = "gorica")
calculate_IC_weights(results$result[,3])
}