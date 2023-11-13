\name{conTest_summary}
\alias{conTest_summary}
\alias{iht_summary}
\alias{conTest_summary.restriktor}

\title{function for computing all available hypothesis tests

}
\description{conTest_summary computes all available hypothesis tests and returns
and object of class \code{conTest} for which a print function is available. The 
\code{conTest_summary} can be used directly and is called by the \code{conTest} 
function if \code{type = "summary"}. 
}

\usage{

\method{conTest_summary}{restriktor}(object, test = "F", ...)

}

\arguments{
  \item{object}{an object of class \code{restriktor}.} 
  
  \item{test}{test statistic; for information about the 
   null-distribution see details.
    \itemize{
      \item for object of class lm; if "F" (default), the classical F-statistic
      is computed. If "Wald", the classical Wald-statistic is computed. If "score", 
      the classical score test statistic is computed. 
      
      \item for object of class rlm; if "F" (default), a robust 
      likelihood ratio type test statistic (Silvapulle, 1992a) is 
      computed. If "Wald", a robust Wald test statistic (Silvapulle, 1992b) 
      is computed. If "score", a score test statistic (Silvapulle, 1996) is 
      computed. 
    }
  }
  
  \item{\dots}{the same arguments as passed to the \code{\link{iht}} function,
  except for \code{type}, of course.}
}


\value{
  An object of class conTest, for which a print is available. 
  More specifically, it is a list with the following items:
  \item{CON}{a list with useful information about the constraints.}
  \item{Amat}{constraints matrix.}
  \item{bvec}{vector of right-hand side elements.}
  \item{meq}{number of equality constraints.}
  \item{meq.alt}{same as input neq.alt.}
  \item{iact}{number of active constraints.}
  \item{type}{same as input.}
  \item{test}{same as input.}
  \item{Ts}{test-statistic value.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{pvalue}{tail probability for \code{Ts}.}
  \item{b.eqrestr}{equality restricted regression coefficients. 
  Only available for \code{type = "A"} and \code{type = "global"}, 
  else \code{b.eqrestr = NULL}.}
  \item{b.unrestr}{unrestricted regression coefficients.}
  \item{b.restr}{restricted regression coefficients.}
  \item{b.restr.alt}{restricted regression coefficients under HA 
  if some equality constraints are maintained.}
  \item{Sigma}{variance-covariance matrix of unrestricted model.}
  \item{R2.org}{unrestricted R-squared.}
  \item{R2.reduced}{restricted R-squared.}
  \item{boot}{same as input.}
  \item{model.org}{original model.}
}


\references{
Shapiro, A. (1988). Towards a unified theory of inequality-constrained 
testing in multivariate analysis. \emph{International Statistical 
Review} \bold{56}, 49--62.

Silvapulle, M. (1992a). Robust tests of inequality constraints and 
one-sided hypotheses in the linear model. \emph{Biometrika}, 
\bold{79}, 621--630.

Silvapulle, M. (1992b). Robust Wald-Type Tests of One-Sided Hypotheses 
in the Linear Model. \emph{Journal of the American Statistical Association}, 
\bold{87}, 156--161.

Silvapulle, M. and Silvapulle, P. (1995). A score test against 
one-sided alternatives. \emph{American statistical association}, 
\bold{90}, 342--349.

Silvapulle, M. (1996) On an F-type statistic for testing one-sided 
hypotheses and computation of chi-bar-squared weights. 
\emph{Statistics and probability letters}, \bold{28}, 137--141.

Silvapulle, M. (1996) Robust bounded influence tests against 
one-sided hypotheses in general parametric models. 
\emph{Statistics and probability letters}, \bold{31}, 45--50.

Silvapulle, M.J. and Sen, P.K. (2005). \emph{Constrained 
Statistical Inference}. Wiley, New York

Wolak, F. (1987). An exact test for multiple inequality and 
equality constraints in the linear regression model. 
\emph{Journal of the American statistical association}, 
\bold{82}, 782--793.
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

iht(fit1.lm, myConstraints1)


# another way is to first fit the restricted model
fit.restr1 <- restriktor(fit1.lm, constraints = myConstraints1)

iht(fit.restr1)

# Or in matrix notation.
Amat1 <- rbind(c(-1, 0,  1),
               c( 0, 1, -1))
myRhs1 <- rep(0L, nrow(Amat1)) 
myNeq1 <- 0

fit1.con <- restriktor(fit1.lm, constraints = Amat1,
                       rhs = myRhs1, neq = myNeq1)
iht(fit1.con)

}