\name{conTest}
\alias{conTest}
\alias{print.conTest}
\alias{summary.conTest.conLM}
\alias{conTestF.conLM}
\alias{conTestLRT.conLM}
\alias{conTestScore.conLM}
\alias{conTestF.conRLM}
\alias{conTestWald.conRLM}
\alias{conTestScore.conRLM}
\alias{conTestC.conLM}
\alias{conTestEq.conLM}
\alias{conTestEq.conRLM}


\title{Function for testing equality and/or inequality restrikted 
 hypotheses for (robust) linear models}
\description{\code{conTest} tests linear equality and/or 
inequality restrikted hypotheses for (robust) linear models by F-, 
Wald-, likelihood ratio-, and score-tests.}

\usage{
conTest(object, type = "summary", ...)

\method{print}{conTest}(x, digits = max(3, getOption("digits") - 4), ...)
\method{summary.conTest.conLM}{conLM and conRLM}(object, test = "F", ...)

\method{conTestF}{conLM and conRLM}(object, type = "A", neq.alt = 0, 
         boot = "no", B = 9999, p.distr = "N", df = 7, 
         parallel = "no", ncpus = 1L, cl = NULL, seed = 1234, 
         verbose = FALSE, control = NULL, ...)

\method{conTestScore}{conLM and conRLM}(object, type = "A", neq.alt = 0, 
             boot = "no", B = 9999, p.distr = "N", df = 7, 
             parallel = "no", ncpus = 1L, cl = NULL, seed = 1234, 
             verbose = FALSE, control = NULL, ...)

\method{conTestLRT}{conLM}(object, type = "A", neq.alt = 0, boot = "no", 
           B = 9999, p.distr = "N", df = 7, parallel = "no", 
           ncpus = 1L, cl = NULL, seed = 1234, verbose = FALSE, 
           control = NULL, ...)

\method{conTestWald}{conRLM}(object, type = "A", neq.alt = 0, boot = "no", 
            B = 9999, p.distr = "N", df = 7, parallel = "no", 
            ncpus = 1L, cl = NULL, seed = 1234, verbose = FALSE, 
            control = NULL, ...)

\method{conTestC}{conLM and conRLM}(object, type = "C", ...)

\method{conTestEq}{conLM and conRLM}(object, test = "F", boot = "no", 
          B = 9999, p.distr = "N", df = 7, parallel = "no", 
          ncpus = 1L, cl = NULL, seed = 1234, verbose = FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{conLM} or \code{conRLM}.} 
  \item{type}{hypothesis test type "A", "B", "C", "global", or 
   "summary" (default). See details for more information.}
  \item{test}{test statistic; for information about the 
   null-distribution see details.
    \itemize{
      \item for object of class lm; if "F" (default), the F-bar 
      statistic (Silvapulle, 1996) is computed. If "LRT", a 
      likelihood ratio test statistic (Silvapulle and Sen, 2005, 
      chp 3.) is computed. If "score", a one-sided score test 
      statistic (Silvapulle, M.J. and Silvapulle, P., 1995) is 
      computed. Note that, in case of equality constraints only, 
      the usual unconstrained F-, Wald-, LR- and score-test 
      statistic is computed. 
      
      \item for object of class rlm; if "F" (default), a robust 
      likelihood ratio type test statistic (Silvapulle, 1992) is 
      computed. If "Wald", a robust Wald test statistic 
      (Silvapulle, 1996) is computed. If "score", a one-sided 
      score test statistic (Silvapulle, 1996) is computed. Note 
      that, in case of equality constraints only, unconstrained 
      robust F-, Wald-, score-test statistics are computed.
    }
  }
  \item{neq.alt}{integer: number of equality restriktions that 
  are maintained under the alternative hypothesis (for hypothesis 
  test type "B"), see example 4.}
  
  \item{boot}{the null-distribution of these test-statistics 
  (except under type "C") takes the form of a mixture of 
  F-distributions. The tail probabilities can be computed directly 
  via bootstrapping; if \code{"parametric"}, the p-value is 
  computed based on the parametric bootstrap. By default, samples 
  are drawn from a normal distribution with mean zero and varance 
  one. See \code{p.distr} for other distributional options. If 
  \code{"model.based"}, a model-based bootstrap method is used. 
  Instead of computing the p-value via simulation, the p-value
  can also be computed based on the chi-bar-square weights. If 
  \code{"no"}, the p-value is computed based on the weights 
  obtained via simulation (\code{bootWt = TRUE}) or using the 
  multivariate normal distribution function (\code{bootWt = FALSE}).}
  
  \item{B}{integer; number of bootstrap draws for \code{boot}. 
  The default value is set to 9999.}
  
  \item{p.distr}{If \code{"N"}, samples for the parametric 
  bootstrap are drawn from the normal distribution (default) with 
  mean zero and variance one. If \code{"t"}, samples are drawn 
  from a t-distribution with df = 7 (degrees of freedom) by 
  default. If \code{"Chi"}, samples are drawn from a chi-square 
  distribution with df = 7 by default. The df can be adjusted via 
  the \code{df} argument.}
  
  \item{df}{degrees of freedom for \code{p.distr}.}
  
  \item{digits}{the number of significant digits to use when printing.}
  
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
  \item{\dots}{no additional arguments for now.}
}

\details{The following hypothesis tests are available:
\itemize{
  \item Type A: Test H0: all restriktions with equalities ("=") 
  active against HA: at least one inequality restriktion (">") 
  strictly true.
  \item Type B: Test H0: all restriktions with inequalities (">") 
  (including some equalities ("=")) active against HA: at least 
  one restriktion false (some equality restriktions may be 
  maintained).
  \item Type C: Test H0: at least one restriktion false ("<") 
  against HA: all restriktions strikty true (">"). This test is 
  based on the intersection-union principle (Silvapulle and Sen, 
  2005, chp 5.3). The null-distribution is based on a 
  t-distribution (one-sided). Note that, this test only makes 
  sense in case of no equality constraints.
  \item Type global: equal to Type A but H0 contains additional 
  equality constraints. This test is analogue to the global 
  F-test in lm, where all coefficients but the intercept equal 0.
}

The exact finite sample distributions of the non-robust F-, 
score- and LR-test statistics based on restrikted OLS estimates 
and normally distributed errors, are a mixture of F-distributions 
under the null hypothesis (Wolak, 1987). In agreement with 
Silvapulle (1992), we found that the results based on these 
mixtures of F-distributions approximate the tail probabilities of 
the robust tests better than their asymptotic distributions. 
Therefore, all p-values for hypothesis test Type \code{"A"}, 
\code{"B"} and \code{"global"} are computed based on mixtures of 
F-distributions. 

Note that, in case of equality constraints only, the 
null-distribution of the (robust) F-test statistics is based on 
an F-distribution. The (robust) Wald- and (robust) score-test 
statistics are based on chi-square distributions.

The \code{summary.conTest.conLM} function is a wrapper function, 
which works for both objects of class \code{conLM} and \code{conRLM}.
The function computes all available hypohtesis tests and prints
a nice summary. 
}

\value{
  An object of class conTest, for which a print is available. 
  More specifically, it is a list with the following items:
  \item{CON}{a list with useful information about the restriktions.}
  \item{type}{same as input.}
  \item{b.eqrestr}{equality restrikted regression coefficients. 
  Only available for \code{type = "A"} and \code{type = "global"}, 
  else \code{b.eqrestr = NULL}.}
  \item{b.unrestr}{unrestrikted regression coefficients.}
  \item{b.restr}{restrikted regression coefficients.}
  \item{b.restr.alt}{restrikted regression coefficients under HA 
  if some equality constraints are maintained.}
  \item{constraints}{matrix with restriktions.}
  \item{rhs}{vector of right-hand side elements.}
  \item{neq}{number of equality restriktions.}
  \item{neq.alt}{same as input.}
  \item{iact}{active restriktions.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{Sigma}{variance-covariance matrix of unrestrikted model.}
  \item{Ts}{test-statistic value.}
  \item{pvalue}{tail probability for \code{Ts}.}
  \item{model.org}{original model.}
}

\references{
Shapiro, A. (1988). Towards a unified theory of inequality-constrained 
testing in multivariate analysis. \emph{International Statistical 
Review} \bold{56}, 49--62.

Silvapulle, M. (1992). Robust tests of inequality constraints and 
one-sided hypotheses in the linear model. \emph{Biometrika}, 
\bold{79}, 621--630.

Silvapulle, M. and Silvapulle, P. (1995). A score test against 
one-sided alternatives. \emph{American statistical association}, 
\bold{90}, 342--349.

Silvapulle, M. (1996) On an F-type statistic for testing one-sided 
hypotheses and computation of chi-bar-squared weights. 
\emph{Statistics & probability letters}, \bold{28}, 137--141.

Silvapulle, M. (1996) Robust bounded influence tests against 
one-sided hypotheses in general parametric models. 
\emph{Statistics & probability letters}, \bold{31}, 45--50.

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
  \code{\link{conTest}},
  \link[sandwich:vcovHC]{sandwich}
}

\examples{
## example 1:
# unrestricted linear model for ages (in months) at which an 
# infant starts to walk alone.
DATA <- ZelazoKolb1972
idx <- which(DATA$Group == 3)
DATA <- DATA[-idx, ]
DATA$Group <- factor(DATA$Group)

# fit unrestrikted linear model
fit1.lm <- lm(Age ~ Group, data = DATA)

# the variable names can be used to impose restriktions on
# the corresponding regression parameters.
coef(fit1.lm)

# restrikted linear model with restriktions that the walking 
# exercises would not have a negative effect of increasing the 
# mean age at which a child starts to walk. 
fit1.con <- restriktor(fit1.lm, constraints = "Group2 > 0; Group2 < Group4")

# compute F-test
conTest(fit1.con)

\dontrun{ 
  # Or in matrix notation.
  myConstraints_mat1 <- rbind(c(0, 1, 0),
                              c(0,-1, 1))
  myRhs1 <- rep(0L, nrow(myConstraints_mat1)) 
  myNeq1 <- 0
  
  fit1.con <- restriktor(fit1.lm, constraints = myConstraints_mat1,
                         rhs = myRhs1, neq = myNeq1)
  conTest(fit1.con)
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

# fit unrestrikted linear model
fit2.lm <- lm(y ~ group, data = DATA2)
coef(fit2.lm)

## example 2: increasing means
myConstraints2 <- ' group2 > 0
                    group2 < group3
                    group3 < group4'

# fit restrikted robust linear model and compute 
# Huber-White (robust) standard errors.
fit2.con <- restriktor(fit2.lm, constraints = myConstraints2, se = "HC0")

# compute score test type A and compute the tail probability based 
# on the parametric bootstrap. We only generate 9 bootstrap samples 
# in this example; in practice you may wish to use a much higher 
# number.
conTest(fit2.con, type = "A", test = "score", boot = "parametric", B = 9)

\dontrun{ 
  # increasing means in matrix notation.
  myConstraints_mat2 <- rbind(c(0, 1, 0, 0),
                              c(0,-1, 1, 0),
                              c(0, 0,-1, 1))
  myRhs2 <- rep(0L, nrow(myConstraints_mat2)) 
  myNeq2 <- 0
  
  fit2.con <- restriktor(fit2.lm, constraints = myConstraints_mat2,
                         rhs = myRhs2, neq = myNeq2)

  conTest(fit2.con, type = "A", test = "score", boot = "parametric", B = 9)
}            

## example 3: equality restriktions only.
myConstraints3 <- 'group2 == 0
                   group2 == group3
                   group3 == group4'
fit3.con <- restriktor(fit2.lm, constr = myConstraints3)
conTest(fit3.con)


## example 4:
# combination of equality and inequality restriktions.
myConstraints4 <- ' group2 == 0
                    group3 < group4
'

# fit restikted model and compute model-based bootstrapped 
# standard errors. We only generate 9 bootstrap samples in this 
# example; in practice you may wish to use a much higher number.
# Note that, a warning message may be thrown because the number of 
# bootstrap samples is too low.
fit4.con <- restriktor(fit2.lm, constr = myConstraints4, 
                       se = "boot.model.based", B = 9)
conTest(fit4.con, type = "B", neq.alt = 1)


## example 5:
# restriktor can also be used to define effects using the := operator 
# and impose restriktions on them. For example, is the 
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
fit5.lm <- lm(y ~ X*Z, data = DATA3)

fit5.con <- restriktor(fit5.lm, constraints = "AVE := X + 16.86137*X.Z; 
                                               AVE > 0")
conTest(fit5.con, test = "LRT")
}