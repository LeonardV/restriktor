\name{conTest}
\alias{conTest}
\alias{iht}

\title{function for informative hypothesis testing (iht)}
\description{\code{conTest} tests linear equality and/or inequality 
restricted hypotheses for linear models. 
}

\usage{
conTest(object, constraints = NULL, type = "summary", test = "F", 
        rhs = NULL, neq = 0, \dots)
}

\arguments{
  \item{object}{an object of class \code{lm} or \code{rlm}. In this 
  case, the constraint syntax needs to be specified
  
  OR
  
  an object of class \code{restriktor}. The constraints are inherited 
  from the fitted restriktor object and do not to be specified again.
  } 
  
  \item{constraints}{there are two ways to constrain parameters. 
  First, the constraint syntax consists of one or more text-based
  descriptions, where the syntax can be specified as a literal 
  string enclosed by single quotes. Only the names of \code{coef(model)}
  can be used as names. See details \code{\link{restriktor}} for more information. 
  
  Second, the constraint syntax consists of a matrix \eqn{R} (or a vector in 
  case of one constraint) and defines the left-hand side of the 
  constraint \eqn{R\theta \ge rhs}, where each row represents one 
  constraint. The number of columns needs to correspond to the 
  number of parameters estimated (\eqn{\theta}) by model. The rows 
  should be linear independent, otherwise the function gives an 
  error. For more information about constructing the matrix \eqn{R} and 
  \eqn{rhs} see the details in the \code{\link{restriktor}} function.}
  
  \item{type}{hypothesis test type "A", "B", "C", "global", or 
   "summary" (default). See details for more information.}
  
  \item{test}{test statistic; for information about the 
   null-distribution see details.
    \itemize{
      \item for object of class lm; if "F" (default), the F-bar 
      statistic (Silvapulle, 1996) is computed. If "LRT", a 
      likelihood ratio test statistic (Silvapulle and Sen, 2005, 
      chp 3.) is computed. If "score", a global score test 
      statistic (Silvapulle and Silvapulle, 1995) is 
      computed. Note that, in case of equality constraints only, 
      the usual unconstrained F-, Wald-, LR- and score-test 
      statistic is computed. 
      
      \item for object of class rlm; if "F" (default), a robust 
      likelihood ratio type test statistic (Silvapulle, 1992a) is 
      computed. If "Wald", a robust Wald test statistic (Silvapulle, 1992b) 
      is computed. If "score", a global score test statistic (Silvapulle, 
      and Silvapulle, 1995) is computed. Note that, in case of equality 
      constraints only, unconstrained robust F-, Wald-, score-test 
      statistics are computed.
      
      \item for object of class glm; if "F" (default), the F-bar 
      statistic (Silvapulle, 1996) is computed. If "LRT", a 
      likelihood ratio test statistic (Silvapulle and Sen, 2005, 
      chp 4.) is computed. If "score", a global score test 
      statistic (Silvapulle and Silvapulle, 1995) is 
      computed. Note that, in case of equality constraints only, 
      the usual unconstrained F-, Wald-, LR- and score-test 
      statistic is computed. 
    }
  }
 
  \item{rhs}{vector on the right-hand side of the constraints; 
  \eqn{R\theta \ge rhs}. The length of this vector equals the 
  number of rows of the constraints matrix \eqn{R} and consists of 
  zeros by default. Note: only used if constraints input is a 
  matrix or vector.}
  
  \item{neq}{integer (default = 0) treating the number of 
  constraints rows as equality constraints instead of inequality 
  constraints. For example, if \code{neq = 2}, this means that the 
  first two rows of the constraints matrix \eqn{R} are treated as 
  equality constraints. Note: only used if constraints input is a 
  matrix or vector.}

  
  \item{\dots}{futher options for the \code{conTest} and/or 
  \code{restriktor} function. See details for more information.}

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
under the null hypothesis (Wolak, 1987). For the robust tests, we 
found that the results based on these mixtures of F-distributions 
approximate the tail probabilities better than their asymptotic 
distributions. 

Note that, in case of equality constraints only, the 
null-distribution of the (non-)robust F-test statistics are 
based on an F-distribution. The (non-)robust Wald- and (non-)robust 
score-test statistics are based on chi-square distributions.

If object is of class \code{lm} or \code{rlm}, the \code{conTest} function
internally calls the \code{restriktor} function. Arguments for the 
\code{\link{restriktor}} function can be passed on via the \code{\ldots}. Additional
arguments for the \code{conTest} function can also passed on via the \code{\ldots}.
See for example \code{\link{conTestF}} for all available arguments.

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
Robertson, T., Wright, F.T. and Dykstra, R.L. (1988). \emph{Order Restricted 
Statistical Inference} New York: Wiley.

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

conTest(fit1.lm, myConstraints1)


# another way is to first fit the restricted model
fit.restr1 <- restriktor(fit1.lm, constraints = myConstraints1)

conTest(fit.restr1)

\dontrun{ 
  # Or in matrix notation.
  Amat1 <- rbind(c(-1, 0,  1),
                 c( 0, 1, -1))
  myRhs1 <- rep(0L, nrow(Amat1)) 
  myNeq1 <- 0
  
  conTest(fit1.lm, constraints = Amat1,
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
conTest(fit2.lm, constraints = myConstraints2, type = "A", 
        boot = "parametric", R = 9)


# or fit restricted linear model
fit2.con <- restriktor(fit2.lm, constraints = myConstraints2)

conTest(fit2.con)

\dontrun{ 
  # increasing means in matrix notation.
  Amat2 <- rbind(c(-1, 1, 0, 0),
                 c( 0,-1, 1, 0),
                 c( 0, 0,-1, 1))
  myRhs2 <- rep(0L, nrow(Amat2)) 
  myNeq2 <- 0

  conTest(fit2.con, constraints = Amat2, rhs = myRhs2, neq = myNeq2, 
          type = "A", boot = "parametric", R = 9)
}            

## example 3: equality constraints only.
myConstraints3 <- ' group1 == group2
                    group2 == group3
                    group3 == group4 '

conTest(fit2.lm, constraints = myConstraints3)

# or
fit3.con <- restriktor(fit2.lm, constraints = myConstraints3)
conTest(fit3.con)


## example 4:
# combination of equality and inequality constraints.
myConstraints4 <- ' group1 == group2
                    group3  < group4 '

conTest(fit2.lm, constraints = myConstraints4, type = "B", neq.alt = 1)

# fit resticted model and compute model-based bootstrapped 
# standard errors. We only generate 9 bootstrap samples in this 
# example; in practice you may wish to use a much higher number.
# Note that, a warning message may be thrown because the number of 
# bootstrap samples is too low.
fit4.con <- restriktor(fit2.lm, constraints = myConstraints4, 
                       se = "boot.model.based", B = 9)
conTest(fit4.con, type = "B", neq.alt = 1)


## example 5:
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
fit5.lm <- lm(y ~ X*Z, data = DATA3)

# constraint syntax
myConstraints5 <- ' AVE := X + 16.86137*X.Z; 
                    AVE > 0 '

conTest(fit5.lm, constraints = myConstraints5)

# or
fit5.con <- restriktor(fit5.lm, constraints = ' AVE := X + 16.86137*X.Z; 
                                                AVE > 0 ')
conTest(fit5.con)
}