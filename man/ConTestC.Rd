\name{conTestC}
\alias{conTestC}
\alias{conTestC.restriktor}

\title{one-sided t-test for iht}
\description{\code{conTestC} tests linear inequality restricted 
hypotheses for (robust) linear models by a one-sided t-test. This 
method is based on the union-intersection principle. It is called 
by the \code{conTest} function if all restrictions are equalities.
For more information see details. 
}

\usage{

\method{conTestC}{restriktor}(object, \dots)

}

\arguments{
  \item{object}{an object of class \code{restriktor}.} 
  
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
  \code{\link{iht}}
}

\examples{
## example 1:
# the data consist of ages (in months) at which an 
# infant starts to walk alone.

# prepare data
DATA1 <- subset(ZelazoKolb1972, Group != "Control")

# fit unrestricted linear model
fit1.lm <- lm(Age ~ -1 + Group, data = DATA1)

# the variable names can be used to impose constraints on
# the corresponding regression parameters.
coef(fit1.lm)

# constraint syntax: assuming that the walking 
# exercises would not have a negative effect of increasing the 
# mean age at which a child starts to walk. 
myConstraints1 <- ' GroupActive  < GroupPassive < GroupNo '

iht(fit1.lm, myConstraints1, type = "C")


# another way is to first fit the restricted model
fit.restr1 <- restriktor(fit1.lm, constraints = myConstraints1)

iht(fit.restr1, type = "C")

\dontrun{ 
  # Or in matrix notation.
  Amat1 <- rbind(c(-1, 0,  1),
                 c( 0, 1, -1))
  myRhs1 <- rep(0L, nrow(Amat1)) 
  myNeq1 <- 0
  
  fit1.con <- restriktor(fit1.lm, constraints = Amat1,
                         rhs = myRhs1, neq = myNeq1)
  iht(fit1.con, type = "C")
}            

}