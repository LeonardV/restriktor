\docType{data}
\name{Kuiper2012estimates}
\alias{Kuiper2012estimates}
\title{Estimates and standard errors from four studies on past experience and buyer trust}
\format{
  An object of class \code{data.frame} with 4 rows and 4 columns.
  \describe{
    \item{study}{Character. Study identifier (e.g., authors).}
    \item{year}{Integer. Publication year of the study.}
    \item{estimate}{Numeric. Estimate of \eqn{\beta_{past}} (effect of previous experience on trust).}
    \item{se}{Numeric. Standard error of \code{estimate}.}
  }
}
\usage{
data(Kuiper2012estimates)
}
\description{
This dataset contains effect estimates for the association between previous experience and buyer trust
from four independent studies. The studies used different statistical models (linear, probit, and
three-level logistic regression), making the raw parameter estimates not directly comparable across studies.
Therefore, GORICA evidence aggregation can be performed using only the study-specific estimates of
\eqn{\beta_{past}} and their uncertainty (standard errors, or equivalently variances).
}
\details{
Each row corresponds to one study. For each study, the hypothesis of interest is
\eqn{H_{+}: \beta_{past} > 0}. Study-specific GORICA weights quantify support for this hypothesis,
and GORICA evidence aggregation combines this support across the four studies to synthesize evidence
for the central theory that previous experience positively relates to trust.
}
\examples{
data(Kuiper2012estimates)

# Inspect the data
Kuiper2012estimates

# Optional: compute variances from standard errors
# Kuiper2012estimates$var <- Kuiper2012estimates$se^2
}
\keyword{datasets}
