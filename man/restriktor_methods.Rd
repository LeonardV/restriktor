\name{restriktor-methods}
\alias{print.conLM}
\alias{summary.conLM}
\alias{print.summary.conLM}
\alias{coef.conLM}
\alias{model.matrix.conLM}

\title{Methods for restriktor}
\description{Restrikted estimation and confidence intervals for 
(robust) linear (in)equality restrikted hypotheses.
}
\usage{\method{print}{conLM and conRLM}(x, digits = max(3, getOption("digits") - 3), ...)

\method{summary}{conLM and conRLM}(object, bootCIs = TRUE, 
        bty = "basic", level = 0.95, GORIC = TRUE, ...)

\method{print.summary.conLM}{summary.conLM and summary.conRLM}(x, digits = max(3, getOption("digits") - 2),
                    signif.stars = getOption("show.signif.stars"), ...)     

\method{coef}{conLM and conRLM}(object, ...)

\method{model.matrix}{conLM and conRLM}(object, ...)
  
}
\arguments{
  \item{object}{an object of class \code{conLM} or \code{conRLM}.}
  \item{x}{for \code{print} an object of class \code{conLM} or 
  \code{conRLM}. For \code{print.summary.conLM} an object of class
  \code{summary.conLM} or \code{summary.conRLM}.}
  \item{bootCIs}{if TRUE (default), nonparametric bootstrap 
  confidence intervals are generated. Only available if \code{object}
  contains \code{bootout} object.}
  \item{bty}{a character string representing the type of interval
  required. The value should be any of the values \code{"norm",
  "basic","perc","bca"}. The value \code{"stud"} is not supported.
  For more details see \code{\link{boot.ci}}.}
  \item{level}{the confidence level of the interval (default = 0.95).}
  \item{GORIC}{if TRUE (default), the Generalized Order-Restrikted 
  Information Criterion value is computed, the penalty term and the
  log-likelihood. The GORIC is a modification of the AIC and the 
  ORIC such that it can be applied to a (in)equality restrikted
  hypotheses in multivariate normal linear models. For more 
  information see \link[goric:goric]{goric}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{If TRUE, ‘significance stars’ are printed 
  for each coefficient.}
  \item{\ldots}{no additional arguments for now.}
}
\details{
  The function \code{print} returns the restrikted coefficients.
  The output from the \code{print.summary.conLM} function provides 
  information that is comparable with the output from 
  \code{print.summary.lm}. Additional information is provided 
  about the unrestrikted and restrikted R-square and by default 
  the output of the GORIC. If bootstrapped standard errors are 
  requested (e.g., option \code{se = "boot.model.based"} in the 
  \code{restriktor} function and \code{bootCI = TRUE} in the 
  summary function) standard errors and confidence intervals 
  are provided.
}

\value{
  The function \code{summary} computes and returns a list of 
  summary statistics of the fitted unrestrikted and restrikted 
  (robust) linear model given in \code{object}, plus
  
  \item{se.type}{type of standard error computed, equal to input 
  \code{se} in the \code{restriktor} function.}
  
  \item{residuals}{the weighted residuals.}
  
  \item{coefficients}{a p x 4 matrix with columns for the 
  estimated coefficient, its standard error, t-statistic and 
  corresponding p-value. If \code{bootCIs = TRUE} and the 
  \code{bootout} object is available in the object, bootstrapped
  standard errors and confidence intervals are produced.}
  
  \item{rdf}{residual degrees of freedom.}
  
  \item{R2.org}{unrestrikted R-squared.}
  
  \item{R2.reduced}{restrikted R-squared.}
  
  \item{goric}{goric value and attributed its penalty term and
  log-likelihood.}
}

\examples{
# unrestricted linear model for ages (in months) at which an 
# infant starts to walk alone.
DATA <- ZelazoKolb1972
idx <- which(DATA$Group == 3)
DATA <- DATA[-idx, ]
DATA$Group <- factor(DATA$Group)

# fit unrestrikted linear model
fit.lm <- lm(Age ~ Group, data = DATA)

# restrikted linear model with restriktions that the walking 
# exercises would not have a negative effect of increasing the 
# mean age at which a child starts to walk. 
fit.con <- restriktor(fit.lm, constraints = "Group2 > 0; Group2 < Group4")

summary(fit.con)  
}
