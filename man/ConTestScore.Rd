\name{conTestScore}
\alias{conTestScore}
\alias{conTestScore.conLM}
\alias{conTestScore.conRLM}
\alias{conTestScore.conGLM}

\title{Score-bar test for iht}
\description{\code{conTestScore} tests linear equality and/or 
inequality restricted hypotheses for (robust) linear models by score-tests. It can be 
used directly and is called by the \code{conTest} function if \code{test = "score"}.

}

\usage{

\method{conTestScore}{conLM}(object, type = "A", neq.alt = 0, 
           boot = "no", R = 9999, p.distr = rnorm, 
           parallel = "no", ncpus = 1L, cl = NULL, seed = 1234, 
           verbose = FALSE, control = NULL, \dots)

\method{conTestScore}{conRLM}(object, type = "A", neq.alt = 0, 
           boot = "no", R = 9999, p.distr = rnorm,  
           parallel = "no", ncpus = 1L, cl = NULL, seed = 1234, 
           verbose = FALSE, control = NULL, \dots)

\method{conTestScore}{conGLM}(object, type = "A", neq.alt = 0, 
           boot = "no", R = 9999, p.distr = rnorm,  
           parallel = "no", ncpus = 1L, cl = NULL, seed = 1234, 
           verbose = FALSE, control = NULL, \dots)

}

\arguments{
  \item{object}{an object of class \code{conLM}, \code{conRLM} or \code{conGLM}.} 
  
  \item{type}{hypothesis test type "A", "B", "C", "global", or 
   "summary" (default). See details for more information.}
  
  \item{neq.alt}{integer: number of equality constraints that 
  are maintained under the alternative hypothesis (for hypothesis 
  test type "B"), see example 3.}
  
  \item{boot}{the null-distribution of these test-statistics 
  (except under type "C", see details) takes the form of a mixture of 
  F-distributions. The tail probabilities can be computed directly 
  via bootstrapping; if \code{"parametric"}, the p-value is 
  computed based on the parametric bootstrap. By default, samples 
  are drawn from a normal distribution with mean zero and varance 
  one. See \code{p.distr} for other distributional options. If 
  \code{"model.based"}, a model-based bootstrap method is used. 
  Instead of computing the p-value via simulation, the p-value
  can also be computed using the chi-bar-square weights. If 
  \code{"no"}, the p-value is computed based on the weights 
  obtained via simulation (\code{mix.weights = "boot"}) or using the 
  multivariate normal distribution function (\code{mix.weights = "pmvnorm"}).
  Note that, these weights are already available in the restriktor objected and
  do not need to be estimated again. However, there are two exception for objects 
  of class \code{conRLM}, namely for computing the p-value for the robust 
  \code{test} = \code{"Wald"} and the robust \code{"score"}. In these cases 
  the weights need to be recalculated.
  }
  
  \item{R}{integer; number of bootstrap draws for \code{boot}. 
  The default value is set to 9999.}
  
  \item{p.distr}{random generation distribution for the parametric bootstrap. For 
  all available distributions see \code{?distributions}. For example, 
  if \code{rnorm}, samples are drawn from the normal distribution (default) 
  with mean zero and variance one. If \code{rt}, samples are drawn from 
  a t-distribution. If \code{rchisq}, samples are drawn from a chi-square 
  distribution. The random generation distributional parameters will 
  be passed in via \ldots.}
  
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
  
  \item{control}{a list of control arguments: 
    \itemize{
      \item \code{absval} tolerance criterion for convergence 
      (default = sqrt(.Machine$double.eps)). Only used for model 
      of class lm.
      \item \code{maxit} the maximum number of iterations for the 
      optimizer (default = 10000). Only used for model of class 
      mlm (not yet supported). 
      \item \code{tol} numerical tolerance value. Estimates 
      smaller than \code{tol} are set to 0.
    }
  }
  \item{\dots}{additional arguments to be passed to the p.distr function.}
}

\details{The following hypothesis tests are available:
\itemize{
  \item Type A: Test H0: all constraints with equalities ("=") 
  active against HA: at least one inequality restriction (">") 
  strictly true.
  \item Type B: Test H0: all constraints with inequalities (">") 
  (including some equalities ("=")) active against HA: at least 
  one restriction false (some equality constraints may be 
  maintained).
  \item Type C: Test H0: at least one restriction false ("<") 
  against HA: all constraints strikty true (">"). This test is 
  based on the intersection-union principle (Silvapulle and Sen, 
  2005, chp 5.3). Note that, this test only makes sense in case 
  of no equality constraints.
  \item Type global: equal to Type A but H0 contains additional 
  equality constraints. This test is analogue to the global 
  F-test in lm, where all coefficients but the intercept equal 0.
}

The null-distribution of hypothesis test Type C is based on a 
t-distribution (one-sided). Its power can be poor in case of many 
inequalty constraints. Its main role is to prevent wrong 
conclusions from significant results from hypothesis test Type A.

The exact finite sample distributions of the non-robust F-, 
score- and LR-test statistics based on restricted OLS estimates 
and normally distributed errors, are a mixture of F-distributions 
under the null hypothesis (Wolak, 1987). In agreement with 
Silvapulle (1992), we found that the results based on these 
mixtures of F-distributions approximate the tail probabilities of 
the robust tests better than their asymptotic distributions. 
Therefore, all p-values for hypothesis test Type \code{"A"}, 
\code{"B"} and \code{"global"} are computed based on mixtures of 
F-distributions. 

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
  if some equality constraints are maintained. Only available for 
  \code{type = "B"} else \code{b.restr.alt = NULL}.}
  \item{Sigma}{variance-covariance matrix of unrestricted model.}
  \item{R2.org}{unrestricted R-squared, not available for objects of class \code{conGLM}.}
  \item{R2.reduced}{restricted R-squared, not available for objects of class \code{conGLM}.}
  \item{boot}{same as input.}
  \item{model.org}{original model.}
}

\references{
Silvapulle, M. and Silvapulle, P. (1995). A score test against 
one-sided alternatives. \emph{American statistical association}, 
\bold{90}, 342--349.

Silvapulle, M. (1996) Robust bounded influence tests against 
one-sided hypotheses in general parametric models. 
\emph{Statistics & probability letters}, \bold{31}, 45--50.

Silvapulle, M.J. and Sen, P.K. (2005). \emph{Constrained 
Statistical Inference}. Wiley, New York

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
fit1.lm <- lm(Age ~ -1 + Group, data = DATA1)

# the variable names can be used to impose constraints on
# the corresponding regression parameters.
coef(fit1.lm)

# constraint syntax: assuming that the walking 
# exercises would not have a negative effect of increasing the 
# mean age at which a child starts to walk. 
myConstraints1 <- ' GroupActive  < GroupPassive; 
                    GroupPassive < GroupNo '

conTest(fit1.lm, myConstraints1, test = "score")


# another way is to first fit the restricted model
fit.restr1 <- restriktor(fit1.lm, constraints = myConstraints1)

conTest(fit.restr1, test = "score")

\dontrun{ 
  # Or in matrix notation.
  Amat1 <- rbind(c(-1, 0,  1),
                 c( 0, 1, -1))
  myRhs1 <- rep(0L, nrow(Amat1)) 
  myNeq1 <- 0
  
  conTest(fit1.lm, constraints = Amat1, test = "score",
          rhs = myRhs1, neq = myNeq1)
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
fit2.lm <- lm(y ~ -1 + group, data = DATA2)
coef(fit2.lm)

## example 2: increasing means
myConstraints2 <- ' group1 < group2
                    group2 < group3
                    group3 < group4 '

# compute F-test for hypothesis test Type A and compute the tail 
# probability based on the parametric bootstrap. We only generate 9 
# bootstrap samples in this example; in practice you may wish to 
# use a much higher number.
conTest(fit2.lm, constraints = myConstraints2, type = "A", test = "score",
        boot = "parametric", R = 9)


# or fit restricted linear model
fit2.con <- restriktor(fit2.lm, constraints = myConstraints2)

conTest(fit2.con, test = "score")

\dontrun{ 
  # increasing means in matrix notation.
  Amat2 <- rbind(c(-1, 1, 0, 0),
                 c( 0,-1, 1, 0),
                 c( 0, 0,-1, 1))
  myRhs2 <- rep(0L, nrow(Amat2)) 
  myNeq2 <- 0

  conTest(fit2.con, constraints = Amat2, rhs = myRhs2, neq = myNeq2, 
          type = "A", test = "score", boot = "parametric", R = 9)
}            

## example 3:
# combination of equality and inequality constraints.
myConstraints3 <- ' group1 == group2
                    group3  < group4 '

conTest(fit2.lm, constraints = myConstraints3, type = "B", test = "score", neq.alt = 1)

# fit resticted model and compute model-based bootstrapped 
# standard errors. We only generate 9 bootstrap samples in this 
# example; in practice you may wish to use a much higher number.
# Note that, a warning message may be thrown because the number of 
# bootstrap samples is too low.
fit3.con <- restriktor(fit2.lm, constraints = myConstraints3, 
                       se = "boot.model.based", B = 9)
conTest(fit3.con, type = "B", test = "score", neq.alt = 1)


## example 4:
# restriktor can also be used to define effects using the := operator 
# and impose constraints on them. For example, is the 
# average effect (AVE) larger than zero?
# generate data
n <- 30
b0 <- 10; b1 = 0.5; b2 = 1; b3 = 1.5
X <- c(rep(c(0), n/2), rep(c(1), n/2))
set.seed(90) 
Z <- rnorm(n, 16, 5)
y <- b0 + b1*X + b2*Z + b3*X*Z + rnorm(n, 0, sd = 10) 
DATA3 = data.frame(cbind(y, X, Z))

# fit linear model with interaction
fit4.lm <- lm(y ~ X*Z, data = DATA3)

# constraint syntax
myConstraints4 <- ' AVE := X + 16.86137*X.Z; 
                    AVE > 0 '

conTest(fit4.lm, constraints = myConstraints4, test = "score")

# or
fit4.con <- restriktor(fit4.lm, constraints = ' AVE := X + 16.86137*X.Z; 
                                                AVE > 0 ')
conTest(fit4.con, test = "score")
}