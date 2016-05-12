\name{conTest}
\alias{conTest}
\alias{conTestF.lm}
\alias{conTestLRT.lm}
\alias{conTestScore.lm}
\alias{conTestF.rlm}
\alias{conTestWald.rlm}
\alias{conTestScore.rlm}

\title{Function for testing (in)equality restrikted hypotheses for
linear models}
\description{\code{conTest} tests linear (in)equality constrained hypotheses for
linear models by F-, Wald-, likelihood ratio-, and score-tests.}

\usage{
conTest(object, type = "A", ...)
\method{conTest}{lm and rlm}(object, type = "A", ...)
}
\arguments{
  \item{object}{an object of class \code{conLM} or \code{conRLM}.} 
  \item{type}{hypothesis test type "A", "B", or "C". See details for more information.}
  \item{test}{test statistic; for information about the null-distributions see details.
    \itemize{
      \item for object of class lm; if "F" (default), the F-bar statistic (Kudo, 1963) is
      computed. If "LRT", a likelihood ratio test statistic (Silvapulle and Sen, 
      2005, chp 3.) is computed. If "score", a one-sided score test statistic 
      (Silvapulle, M.J. and Silvapulle, P., 1995) is computed. 
      
      \item for object of class rlm; if "F" (default), a robust likelihood ratio 
      type test statistic (Silvapulle, 1992) is computed. If "Wald", a robust 
      Wald test statistic (Silvapulle, 1996) is computed. If "score", a one-sided 
      score test statistic (Silvapulle, 1995) is computed.
    }
  }
  \item{neq.alt}{integer: number of equality restriktions that are maintained 
  under the alternative hypothesis (for hypothesis test type "B"), see Xrd example.}
  \item{boot}{the null-distribution of these test-statistics (except under type 
  "C") takes the form of a mixture of F-distributions. The tail probabilities can 
  be computed in several ways; if \code{"parametric"}, the p-value is computed 
  based on the parametric bootstrap. By default, samples are drawn from a normal 
  distribution with mean zero and varance one. See \code{p.distr} for other 
  distributional options. If \code{"model.based"}, a model-based bootstrap method 
  is used. Instead of computing the p-value directly, two approaches are available 
  to compute the mixing weights. If \code{"mix.weights"}, mixing weights are 
  computed based on a simulation approach. If \code{"no"} (default), no bootstrap 
  method is used and the mixing weights are computed using the multivariate normal 
  distribution function.}
  \item{B}{integer; number of bootstrap draws. The default value is set to 9999.}
  \item{p.distr}{If \code{"N"}, samples for the parametric bootstrap are drawn 
  from the normal distribution (default) with mean zero and variance one. If 
  \code{"t"}, samples are drawn from a t-distribution with df = 7 (degrees of 
  freedom) by default. If \code{"chi"}, samples are drawn from a chi-square 
  distribution with df = 7 by default. The df can be adjusted via the \code{df} 
  argument.}
  \item{df}{degrees of freedom for \code{p.distr}.}
  \item{parallel}{the type of parallel operation to be used (if any). If missing, 
  the default is set "no".}
  \item{ncpus}{integer: number of processes to be used in parallel operation: 
  typically one would chose this to the number of available CPUs.}
  \item{cl}{an optional parallel or snow cluster for use if parallel = "snow". 
  If not supplied, a cluster on the local machine is created for the duration of 
  the conTest call.}
  \item{seed}{seed value.}
  \item{verbose}{logical; if TRUE, information is shown at each bootstrap draw.}
  \item{control}{a list of control arguments: 
    \itemize{
      \item \code{absval} tolerance criterion for convergence 
      (default = sqrt(.Machine$double.eps)). Only used for model of class lm.
      \item \code{maxit} the maximum number of iterations for the optimizer 
      (default = 10000). Only used for model of class mlm (not yet available). 
      \item \code{tol} numerical tolerance value. Estimates smaller than \code{tol} 
      are set to 0.
    }
  }
  \item{\dots}{no additional arguments for now.}
}

\details{
The following hypothesis tests are implemented:
\itemize{
  \item Type A: H0: restrictions valid with equality vs. HA: at least one inequality true
  \item Type B: H0: all restrictions true vs. HA: at least one restriction false
  \item Type C: H0: restrictions false vs. HA: restrictions true (with inequality)
}


TP=11: H0: restriction valid with equality and further linear equalities vs. H1: at least one equality from H0 violated, restriction valid

TP=21: H0: restrictions valid (including some equality restrictions) vs. H1: at least one restriction from H0 violated, some equality restrictions are maintained

Note that TPs 1 and 11 can reject H0 even if H1 is violated by the data. Rejection of H0 does not provide evidence for H1 (but only against H0) in these TPs because H1 is not the opposite of H0. The tests concentrate their power in H1, but are only guaranteed to observe their level for the stated H0.

Also note that TP 3 does not make sense if obj involves equality restrictions (obj$meq>0).

Under TPs 1, 2, 11, and 21, the distributions of test statistics are mixtures of chi-square distributions (df.error=Inf) or beta-distributions (df.error finite) with different degrees of freedom (chi-square) or parameter combinations (beta). Shapiro (1988) gives detailed information on the mixing weights for the different scenarios. Basically, there are two different situations:

If meq=0, the weights are probabilities that a random variable with covariance matrix ui%*%cov%*%t(ui) is realized in the positive orthant or its lower-dimensional faces, respectively (if ui has too few columns, blow up by columns of 0s in appropriate positions) (Shapiro, formulae (5.5) or (5.10), respectively).

If meq > 0 (but not all restrictions are equality restrictions), the weights are probabilities that a random variable with covariance matrix the inverse of the lower right corner of solve(ui%*%cov%*%t(ui)) is realized in the positive orthant or its lower-dimensional faces, respectively (Shapiro, formula (5.9)).

These weights must then be combined with the appropriate degrees of freedom - these can be worked out by realizing that either the null hypothesis or the alternative hypothesis has fixed dimension and the respective mixing degrees of freedom are obtained by taking the difference to the dimension of the respective other hypothesis, which is correct because - given a certain dimension of the inequality-restricted estimate, the inequality-restricted estimate is a projection onto a linear space of that dimension.

The test for TP 3 (cf. e.g. Sasabuchi 1980) is based on the intersection-union principle and simply obtains its p-value as the maximum p-value from testing the individual restrictions.


The exact finite sample distributions of the non-robust F-, score- and 
LR-test statistics based on OLS estimates and normally distributed errors, are 
a mixture of F-distributions under the null hypothesis (Wolak, 1987). In 
agreement with Silvapulle (1992b), we found that the results based on these 
mixtures of F-distributions approximate the tail probabilities of the robust 
tests better than their asymptotic distributions. Therefore, all p-values are 
computed based on mixtures of F-distributions.


}

\value{
  An object of class conTest, for which a print and a 
  summary method is available. More specifically, it is a list 
  with the following items:
    
  \item{CON}{a list with useful information about the restriktions.}
  \item{type}{}
  \item{b.eqrestr}{}
  \item{b.unrestr}{}
  \item{b.restr}{}
  \item{b.restr.alt}{}
  \item{constraints}{}
  \item{rhs}{}
  \item{neq}{}
  \item{neq.alt}{}
  \item{iact}{}
  \item{df.residual}{}
  \item{Sigma}{}
  \item{Ts}{}
  \item{pvalue}{}
  \item{model.org}{}
}

\references{
Kudo, A. (1963). A multivariate analogue of the one-sided test. \emph{Biometrika}, 
\bold{50}, 403--418.

Silvapulle, M. (1992a). Robust tests of inequality constraints and one-sided 
hypotheses in the linear model. \emph{Biometrika}, \bold{79}, 621--630.

Silvapulle, M. (1992b). Robust wald-type tests of one-sided hypotheses in 
the linear model. \emph{Journal of the American Statistical Association}, 
\bold{87}, 156--161.

Silvapulle, M. and Silvapulle, P. (1995). A score test against one-sided 
alternatives. \emph{American statistical association}, \bold{90}, 342--349.

Silvapulle, M. (1996). Robust bounded influence tests against one-sided 
hypotheses in general parametric models. \emph{Statistics & probability 
letters}, \bold{31}, 45--50.

Silvapulle, M.J. and Sen, P.K. (2005). \emph{Constrained Statistical Inference}. 
Wiley, New York

Wolak, F. (1987). An exact test for multiple inequality and equality constraints 
in the linear regression model. \emph{Journal of the American statistical 
association}, \bold{82}, 782--793.
}

\author{Leonard Vanbrabant and Yves Rosseel}


\seealso{ 
  \code{\link{solve.QP}}, 
  \code{\link{conTest}},
  \code{\link{vcovHC}}
}

\examples{
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
  summary(fit1.con)
  
  \dontrun{ 
    # Or in matrix notation.
    myConstraints_mat1 <- rbind(c(0, 1, 0),
                                c(0,-1, 1))
    myRhs1 <- rep(0L, nrow(myConstraints_mat1)) 
    myNeq1 <- 0
    
    fit1.con <- restriktor(fit1.lm, constraints = myConstraints_mat1,
                           rhs = myRhs1, neq = myNeq1)
    summary(fit1.con)
  }            
  
  #########################
  ## Artificial examples ##
  #########################
  # generate data
  n <- 10
  means <- c(1,2,1,3)
  nm <- length(means)
  group <- as.factor(rep(1:nm, each = n))
  y <- cbind(c(MASS:::mvrnorm(n, mu = means, Sigma = diag(nm)), nrow(n)))
  DATA2 <- data.frame(y, group)
  
  # fit unrestrikted robust linear model
  fit2.rlm <- MASS:::rlm(y ~ group, data = DATA2, method = "MM")
  coef(fit2.rlm)
  
  # increasing means
  myConstraints2 <- ' group2 > 0
  group2 < group3
  group3 < group4
  '
  
  # fit restrikted robust linear model and compute 
  # Huber-White (robust) standard errors.
  fit2.con <- restriktor(fit2.rlm, constraints = myConstraints2, se = "HC0")
  summary(fit2.con)
  
  \dontrun{ 
    # increasing means in matrix notation.
    myConstraints_mat2 <- rbind(c(0, 1, 0, 0),
                                c(0,-1, 1, 0),
                                c(0, 0,-1, 1))
    myRhs2 <- rep(0L, nrow(myConstraints_mat2)) 
    myNeq2 <- 0
    
    fit2.con <- restriktor(fit2.rlm, constraints = myConstraints_mat2,
                           rhs = myRhs2, neq = myNeq2)
    summary(fit2.con)
  }            
  
  # equality restriktions only.
  myConstraints3 <- ' group2 == 0
  group2 == group3
  group3 == group4
  '
  fit3.con <- restriktor(fit2.lm, constr = myConstraints3)
  summary(fit3.con)
  
  
  # combination of equality and inequality restriktions.
  myConstraints4 <- ' group2 == 0
  group3 < group4
  '
  
  # fit restikted model and compute model-based bootstrapped 
  # standard errors. We only generate 9 bootstrap samples in this 
  # example; in practice you may wish to use a much higher number.
  fit4.con <- restriktor(fit2.lm, constr = myConstraints4, 
                         se = "boot.model.based", B = 9)
  summary(fit4.con)
  
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
  summary(fit5.con)
}