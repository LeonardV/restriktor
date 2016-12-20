\name{conTest_ceq}
\alias{conTest_ceq.conLM}
\alias{conTest_ceq.conRLM}

\title{F-, Wald-, and score-tests for (robust) iht with equality constraints only}
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
fit1_con <- restriktor(fit1_lm, constraints = ' GroupActive  == GroupPassive; 
                                                GroupPassive == GroupNo ')

# compute F-test
conTest(fit1_con, test = "F")

\dontrun{ 
  # Or in matrix notation.
  Amat1 <- rbind(c(-1, 1, 0),
                 c( 0,-1, 1))
  myRhs1 <- rep(0L, nrow(Amat1)) 
  myNeq1 <- 2
  
  fit1_con <- restriktor(fit1_lm, constraints = Amat1,
                         rhs = myRhs1, neq = myNeq1)
  conTest(fit1_con, test = "F")
}            

#########################
## Artificial examples ##
#########################
# generate data
n <- 10
means <- c(1,2,1,3)
nm <- length(means)
group <- as.factor(rep(1:nm, each = n))
y <- rnorm(n * nm, rep(means, each = n))
DATA2 <- data.frame(y, group)

# fit unrestricted linear model
fit2_lm <- lm(y ~ -1 + group, data = DATA2)
coef(fit2_lm)

## example 2: increasing means
myConstraints2 <- ' group1 == group2
                    group2 == group3
                    group3 == group4 '

# fit restricted robust linear model and compute 
# Huber-White (robust) standard errors.
fit2_con <- restriktor(fit2_lm, constraints = myConstraints2, se = "HC0")

# compute F test for hypothesis test Type A and compute the tail 
# probability based on the parametric bootstrap. We only generate 9 
# bootstrap samples in this example; in practice you may wish to 
# use a much higher number.
conTest(fit2_con, type = "A", test = "F", boot = "parametric", R = 9)

\dontrun{ 
  # increasing means in matrix notation.
  Amat2 <- rbind(c(-1, 1, 0, 0),
                 c( 0,-1, 1, 0),
                 c( 0, 0,-1, 1))
  myRhs2 <- rep(0L, nrow(Amat2)) 
  myNeq2 <- 3
  
  fit2_con <- restriktor(fit2_lm, constraints = Amat2,
                         rhs = myRhs2, neq = myNeq2)

  conTest(fit2_con, type = "A", test = "F", boot = "parametric", R = 9)
}            


}