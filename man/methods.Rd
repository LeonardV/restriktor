\name{restriktor-methods}
\alias{summary.conLM}
\alias{coef.conLM}
\alias{model.matrix.conLM}

\title{Methods for restriktor}
\description{Restrikted estimation and confidence intervals for 
linear (in)equality restrikted hypotheses.
}
\usage{
  \method{summary}{conLM and conRLM}(object, bootCIs = TRUE, 
        bty = "basic", level = 0.95,
        GORIC = TRUE, digits = max(3, getOption("digits") - 2),
        signif.stars = getOption("show.signif.stars"), ...)
        
  \method{coef}{conLM and conRLM}(object, ...)
  \method{model.matrix}{conLM and conRLM}(object, ...)
  
}
\arguments{
  \item{object}{an object of class \code{conLM} or \code{conRLM}.}
  \item{bootCIs}{if TRUE (default), nonparametric bootstrap 
  confidence intervals are generated. Only available if \code{object}
  contains \code{bootout} object.}
  \item{bty}{a character string representing the type of interval
  required. The value should be any of the values \code{"norm",
  "basic","perc","bca"}. The value \code{"stud"} is not supported.
  For more details see \code{\link{boot.ci}}.}
  \item{level}{the confidence level of the interval (default = 0.95).}
  \item{GORIC}{if TRUE (default), a Generalized Order-Restrikted 
  Information Criterion value is computed, a penalty term and the
  log-likelihood. The GORIC is a modification of the AIC and the 
  ORIC such that it can be applied to a (in)equality restrikted
  hypotheses in multivariate normal linear models. For more 
  information see \code{\link{goric}}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{If TRUE, ‘significance stars’ are printed for each coefficient.}
  \item{\ldots}{no additional arguments for now.}
}
\details{
  The methods for general linear hypotheses as described by objects returned
  by \code{\link{glht}} can be used to actually test the global
  null hypothesis, each of the partial hypotheses and for
  simultaneous confidence intervals for the linear function \eqn{K \theta}.
  
  The \code{\link{coef}} and \code{\link{vcov}} methods compute the linear
  function \eqn{K \hat{\theta}} and its covariance, respectively.
  
  The \code{test} argument to \code{summary} takes a function specifying
  the type of test to be applied. Classical Chisq (Wald test) or F statistics
  for testing the global hypothesis \eqn{H_0} are implemented in functions
  \code{Chisqtest} and \code{Ftest}. Several approaches to multiplicity adjusted p 
  values for each of the linear hypotheses are implemented 
  in function \code{adjusted}. The \code{type}
  argument to \code{adjusted} specifies the method to be applied:
    \code{"single-step"} implements adjusted p values based on the joint
  normal or t distribution of the linear function, and
  \code{"Shaffer"} and \code{"Westfall"} implement logically constraint 
  multiplicity adjustments (Shaffer, 1986; Westfall, 1997). 
  \code{"free"} implements multiple testing procedures under free 
  combinations (Westfall et al, 1999).
  In addition, all adjustment methods
  implemented in \code{\link{p.adjust}} are available as well.
  
  Simultaneous confidence intervals for linear functions can be computed
  using method \code{\link{confint}}. Univariate confidence intervals
  can be computed by specifying \code{calpha = univariate_calpha()}
  to \code{confint}. The critical value can directly be specified as a scalar 
  to \code{calpha} as well. Note that \code{plot(a)} for some object \code{a} of class
  \code{glht} is equivalent to \code{plot(confint(a))}.
  
  All simultaneous inference procedures implemented here control
  the family-wise error rate (FWER). Multivariate
  normal and t distributions, the latter one only for models of 
  class \code{\link{lm}}, are evaluated using the procedures
  implemented in package \code{mvtnorm}.
  
  A more detailed description of the underlying methodology is
  available from Hothorn et al. (2008) and Bretz et al. (2010).
  
  
}
\value{
  
  \code{summary} computes (adjusted) p values for general linear hypotheses,
  \code{\link{confint}} computes (adjusted) confidence intervals. 
  \code{\link{coef}} returns estimates of the linear function \eqn{K \theta}
  and \code{\link{vcov}} its covariance. 
  
}
\references{
  
  
}
\examples{
  
}
