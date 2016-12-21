\name{conTestC}
\alias{conTestC.conLM}
\alias{conTestC.conRLM}

\title{Hypothesis test Type C
}
\description{\code{conTestC} tests linear inequality restricted hypotheses for (robust) 
linear models by a one-sided t-test. This method is based on the union-intersection
principle. For more information see details. Robust linear models are also supported
via the \code{conTestC.conLM} function.
}

\usage{

\method{conTestC}{conLM}(object, \dots)

}

\arguments{
  \item{object}{an object of class \code{conLM} or \code{conRLM}.} 
  
  \item{\dots}{no additional arguments for now.}
}


\details{Hypothesis test Type C:
\itemize{
  \item Test H0: at least one restriction false ("<") 
  against HA: all constraints strikty true (">"). This test is 
  based on the intersection-union principle. Note that, this test only makes 
  sense in case of no equality constraints.
}

The null-distribution of hypothesis test Type C is based on a 
t-distribution (one-sided). Its power can be poor in case of many 
inequalty constraints. Its main role is to prevent wrong 
conclusions from significant results from hypothesis test Type A.

}

\value{
  An object of class conTest, for which a print is available. 
  More specifically, it is a list with the following items:
  \item{CON}{a list with useful information about the constraints.}
  \item{Amat}{constraints matrix.}
  \item{bvec}{vector of right-hand side elements.}
  \item{meq}{number of equality constraints.}
  \item{test}{same as input.}
  \item{Ts}{test-statistic value.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{pvalue}{tail probability for \code{Ts}.}
  \item{b.unrestr}{unrestricted regression coefficients.}
  \item{b.restr}{restricted regression coefficients.}
  \item{Sigma}{variance-covariance matrix of unrestricted model.}
  \item{R2.org}{unrestricted R-squared.}
  \item{R2.reduced}{restricted R-squared.}
  \item{boot}{"no", not used (yet).}
  \item{model.org}{original model.}
  
}

\references{
Silvapulle, M.J. and Sen, P.K. (2005, chapter 5.). \emph{Constrained 
Statistical Inference}. Wiley, New York

}

\author{Leonard Vanbrabant and Yves Rosseel}


\seealso{ 
  \link[quadprog:solve.QP]{quadprog}, 
  \code{\link{conTest}}
}

\examples{
## example 1:
# unrestricted linear model for ages (in months) at which an 
# infant starts to walk alone.

# prepare data
DATA1 <- subset(ZelazoKolb1972, Group != "Control")

# fit unrestricted linear model
fit1_lm <- lm(Age ~ -1 + Group, data = DATA1)

# the variable names can be used to impose constraints on
# the corresponding regression parameters.
coef(fit1_lm)

# restricted linear model with constraints that the walking 
# exercises would not have a negative effect of increasing the 
# mean age at which a child starts to walk. 
fit1_con <- restriktor(fit1_lm, constraints = ' GroupActive  < GroupPassive; 
                                                GroupPassive < GroupNo ')

# compute F-test
conTest(fit1_con, type = "C")

\dontrun{ 
  # Or in matrix notation.
  Amat1 <- rbind(c(-1, 1, 0),
                 c( 0,-1, 1))
  myRhs1 <- rep(0L, nrow(Amat1)) 
  myNeq1 <- 0
  
  fit1_con <- restriktor(fit1_lm, constraints = Amat1,
                         rhs = myRhs1, neq = myNeq1)
  conTest(fit1_con, type = "C")
}            



}