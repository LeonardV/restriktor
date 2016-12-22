\name{conTest_ceq}
\alias{conTest_ceq}
\alias{conTest_ceq.conLM}
\alias{conTest_ceq.conRLM}

\title{Tests for iht with equality constraints only}
\description{\code{conTest_ceq} tests linear equality restricted hypotheses for
(robust) linear models by F-, Wald-, and score-tests. It can be used directly 
and is called by the \code{conTest} function if all restrictions are equalities.
}

\usage{

\method{conTest_ceq}{conLM}(object, test = "F", boot = "no", 
            R = 9999, p.distr = "N", df = 7, parallel = "no", 
            ncpus = 1L, cl = NULL, seed = 1234, verbose = FALSE, \dots)

\method{conTest_ceq}{conRLM}(object, test = "F", boot = "no", 
            R = 9999, p.distr = "N", df = 7, parallel = "no", 
            ncpus = 1L, cl = NULL, seed = 1234, verbose = FALSE, \dots)
}

\arguments{
  \item{object}{an object of class \code{conLM} or \code{conRLM}.} 
  
  \item{test}{test statistic; for information about the 
   null-distribution see details.
    \itemize{
      \item for object of class lm; if "F" (default), the classical F-statistic
      is computed. If "Wald", the classical Wald-statistic is computed. If "score", 
      the classical score test statistic is computed. 
      
      \item for object of class rlm; if "F" (default), a robust 
      likelihood ratio type test statistic (Silvapulle, 1992a) is 
      computed. If "Wald", a robust Wald test statistic 
      (Silvapulle, 1996) is computed. If "Wald2", a robust Wald test 
      statistic (Silvapulle, 1992b) is computed. If "score", a 
      score test statistic (Silvapulle, 1996) is computed. 
    }
  }
  
  \item{boot}{if \code{"parametric"}, the p-value is computed based on the 
  parametric bootstrap. If \code{"model.based"}, a model-based bootstrap method 
  is used. 
  }
  
  \item{R}{integer; number of bootstrap draws for \code{boot}. 
  The default value is set to 9999.}
  
  \item{p.distr}{If \code{"N"}, samples for the parametric 
  bootstrap are drawn from the normal distribution (default) with 
  mean zero and variance one. If \code{"t"}, samples are drawn 
  from a t-distribution with df = 7 (degrees of freedom) by 
  default. If \code{"Chi"}, samples are drawn from a chi-square 
  distribution with df = 7 by default. The df can be adjusted via 
  the \code{df} argument.}
  
  \item{df}{degrees of freedom for \code{p.distr}.}
  
  \item{parallel}{the type of parallel operation to be used 
  (if any). If missing, the default is set "no".}
  
  \item{ncpus}{integer: number of processes to be used in parallel 
  operation: typically one would chose this to the number of 
  available CPUs.}
  
  \item{cl}{an optional parallel or snow cluster for use if 
  parallel = "snow". If not supplied, a cluster on the local 
  machine is created for the duration of the conTest call.}
  
  \item{seed}{seed value. The default value is set to 1234.}
  
  \item{verbose}{logical; if TRUE, information is shown at each 
  bootstrap draw.}
  
  \item{\dots}{no additional arguments for now.}
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
  \item{R2.org}{unrestricted R-squared.}
  \item{R2.reduced}{restricted R-squared.}
}

\references{
Silvapulle, M. (1992a). Robust tests of inequality constraints and 
one-sided hypotheses in the linear model. \emph{Biometrika}, 
\bold{79}, 621--630.

Silvapulle, M. (1996) Robust bounded influence tests against 
one-sided hypotheses in general parametric models. 
\emph{Statistics & probability letters}, \bold{31}, 45--50.

Silvapulle, M. (1992b). Robust Wald-Type Tests of One-Sided Hypotheses 
in the Linear Model. \emph{Journal of the American Statistical Association}, 
\bold{87}, 156--161.

Silvapulle, M. (1996) Robust bounded influence tests against one-sided hypotheses
in general parametric models. \emph{Statistics & probability letters}, 
\bold{31}, 45--50.

}

\author{Leonard Vanbrabant and Yves Rosseel}


\seealso{ 
  \link[quadprog:solve.QP]{quadprog}, 
  \code{\link{conTest}}
}

\examples{
## example 1:
# the data consist of ages (in months) at which an 
# infant starts to walk alone.

# prepare data
DATA1 <- subset(ZelazoKolb1972, Group != "Control")

# fit unrestricted linear model
fit1_lm <- lm(Age ~ -1 + Group, data = DATA1)

# the variable names can be used to impose constraints on
# the corresponding regression parameters.
coef(fit1_lm)

# constraint syntax: assuming that the walking 
# exercises would not have a negative effect of increasing the 
# mean age at which a child starts to walk. 
myConstraints1 <- ' GroupActive  == GroupPassive; 
                    GroupPassive == GroupNo '

conTest(fit1_lm, myConstraints1)


# another way is to first fit the restricted model
fit_restr1 <- restriktor(fit1_lm, constraints = myConstraints1)

conTest(fit_restr1)

\dontrun{ 
  # Or in matrix notation.
  Amat1 <- rbind(c(-1, 0,  1),
                 c( 0, 1, -1))
  myRhs1 <- rep(0L, nrow(Amat1)) 
  myNeq1 <- 2
  
  conTest(fit1_lm, constraints = Amat1,
          rhs = myRhs1, neq = myNeq1)
}            

}