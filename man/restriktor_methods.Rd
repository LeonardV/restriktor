\name{restriktor-methods}
\alias{restriktor-methods}
\alias{print.restriktor}
\alias{summary.restriktor}
\alias{print.summary.restriktor}
\alias{coef.restriktor}
\alias{model.matrix.restriktor}
\alias{logLik.restriktor}

\title{Methods for restriktor}
\description{restricted estimation and confidence intervals for  
(robust) linear (in)equality restricted hypotheses.
}
\usage{

\method{print}{restriktor}(x, digits = max(3, getOption("digits") - 2), \dots)

\method{summary}{restriktor}(object, bootCIs = TRUE, 
        bty = "perc", level = 0.95, goric = "goric", \dots)

\method{print}{summary.restriktor}(x, digits = max(3, getOption("digits") - 2),
              signif.stars = getOption("show.signif.stars"), \dots)     

\method{coef}{restriktor}(object, \dots)

\method{model.matrix}{restriktor}(object, \dots)

\method{logLik}{restriktor}(object, \dots)
  
}
\arguments{
  \item{object}{an object of class \code{restriktor}.}
  
  \item{x}{an object of class \code{restriktor}.}
  
  \item{bootCIs}{if TRUE (default), nonparametric bootstrap 
  confidence intervals are generated. Only available if \code{object}
  contains \code{bootout} object.}
  
  \item{bty}{a character string representing the type of interval
  required. The value should be any of the values \code{"norm",
  "basic","perc","bca"}. The value \code{"stud"} is not supported.
  For more details see \code{\link{boot.ci}}.}
  
  \item{level}{the confidence level of the interval (default = 0.95).}
  
  \item{goric}{if \code{"goric"} (default), the generalized order-restricted 
  information criterion value is computed. If \code{"gorica"} the 
  log-likihood is computed using the multivariate normal distribution 
  function. If \code{"goricc" or "goricca"}, a small sample version
  of the \code{"goric" or "gorica"} is computed.}
  
  \item{digits}{the number of significant digits to use when printing.}
  
  \item{signif.stars}{If TRUE, "significance stars are printed 
  for each coefficient.}
  
  \item{\ldots}{no additional arguments for now.}
}
\details{
  The function \code{print} returns the restricted coefficients.
  The output from the \code{print.summary.conLM} function provides 
  information that is comparable with the output from 
  \code{print.summary.lm}. Additional information is provided 
  about the unrestricted and restricted R-square and by default 
  the output of the GORIC. If bootstrapped standard errors are 
  requested (e.g., option \code{se = "boot.model.based"} in the 
  \code{restriktor} function and \code{bootCI = TRUE} in the 
  summary function) standard errors and confidence intervals 
  are provided.
}

\value{
  The function \code{summary} computes and returns a list of 
  summary statistics of the fitted unrestricted and restricted 
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
  
  \item{R2.org}{unrestricted R-squared.}
  
  \item{R2.reduced}{restricted R-squared.}
  
  \item{goric}{goric value and attributed its penalty term and
  log-likelihood.}
}

\examples{
  # unrestricted linear model for ages (in months) at which an 
  # infant starts to walk alone.
  
  # prepare data
  DATA <- subset(ZelazoKolb1972, Group != "Control")
  
  # fit unrestricted linear model
  fit.lm <- lm(Age ~ -1 + Group, data = DATA)
  
  # restricted linear model with restrictions that the walking 
  # exercises would not have a negative effect of increasing the 
  # mean age at which a child starts to walk. 
  fit.con <- restriktor(fit.lm, constraints = ' GroupActive  < GroupPassive; 
                                                GroupPassive < GroupNo ')
  
  summary(fit.con)
}
