\name{iht}
\alias{iht}
\alias{conTest}
\alias{contest}
\alias{conTestD}
\alias{contestD}

\title{function for informative hypothesis testing (iht)}
\description{\code{iht} tests linear equality and/or inequality 
restricted hypotheses for linear models. 
}

\usage{
iht(\ldots)

conTest(object, constraints = NULL, type = "summary", test = "F", 
        rhs = NULL, neq = 0, \dots)

conTestD(model = NULL, data = NULL, constraints = NULL, 
         type = c("A","B"), R = 1000L, bootstrap.type = "bollen.stine", 
         return.test = TRUE, neq.alt = 0, 
         double.bootstrap = "standard", double.bootstrap.R = 249, 
         double.bootstrap.alpha = 0.05, 
         parallel = c("no", "multicore", "snow"), 
         ncpus = 1L, cl = NULL, verbose = FALSE, \ldots)        
}

\arguments{
  \item{object}{an object of class \code{lm} or \code{rlm}. In this 
  case, the constraint syntax needs to be specified
  
  OR
  
  an object of class \code{restriktor}. The constraints are inherited 
  from the fitted restriktor object and do not to be specified again.
  } 
  
  \item{model}{lavaan model syntax specifying the model. See \code{\link[lavaan]{model.syntax}} 
    for more information. }

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
  
  \item{data}{the data frame containing the observed variables being used to 
    fit the lavaan model.}

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
  \item{neq.alt}{integer: number of equality constraints that are maintained under 
    the alternative hypothesis (for hypothesis test type "B").}

  \item{R}{Integer; number of bootstrap draws. The default value is set to 1000.}
  \item{bootstrap.type}{If \code{"parametric"}, the parametric bootstrap is used. 
    If \code{"bollen.stine"}, the semi-nonparametric Bollen-Stine bootstrap 
    is used. The default is set to \code{"bollen.stine"}.}
  \item{return.test}{Logical; if \code{TRUE}, the function returns bootstrapped         
    test-values.}
  \item{double.bootstrap}{If \code{"standard"} (default) the genuine double bootstrap is 
    used to compute an additional set of plug-in p-values for each bootstrap       
    sample. If \code{"no"}, no double bootstrap is used. If \code{"FDB"}, 
    the fast double bootstrap is used to compute second level LRT-values for 
    each bootstrap sample. Note that the \code{"FDB"} is experimental and should 
    not be used by inexperienced users.}    
  \item{double.bootstrap.R}{Integer; number of double bootstrap draws. The default 
    value is set to 249.}
  \item{double.bootstrap.alpha}{The significance level to compute the adjusted 
    alpha based on the plugin p-values. Only used if \code{double.bootstrap = "standard"}. 
    The default value is set to 0.05.}
  \item{parallel}{The type of parallel operation to be used (if any). If missing, 
    the default is set "no".}
  \item{ncpus}{Integer: number of processes to be used in parallel operation: 
      typically one would chose this to the number of available CPUs.} 
  \item{cl}{An optional parallel or snow cluster for use if 
    \code{parallel = "snow"}. If not supplied, a cluster on the local machine 
    is created for the duration of the \code{InformativeTesting} call.}
  \item{verbose}{Logical; if \code{TRUE}, information is shown at each bootstrap         
    draw.}
  \item{\dots}{futher options for the \code{iht} and/or 
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
  An object of class conTest or conTestLavaan for which a print is available. 
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
  \emph{Statistics and probability letters}, \bold{28}, 137--141.
  
  Silvapulle, M. (1996) Robust bounded influence tests against 
  one-sided hypotheses in general parametric models. 
  \emph{Statistics and probability letters}, \bold{31}, 45--50.
  
  Silvapulle, M.J. and Sen, P.K. (2005). \emph{Constrained 
  Statistical Inference}. Wiley, New York
  
  Vanbrabant, L., Van de Schoot, R., Van Loey, N.E.E. and Rosseel, Y. (2017). 
  A General Procedure for Testing Inequality Constrained Hypotheses in SEM.
  \emph{Methodology European Journal of Research Methods for the Behavioral and Social Sciences},
  \bold{13}, 61-70.
  
  Van de Schoot, R., Hoijtink, H., and Dekovic, M. (2010). 
  Testing inequality constrained hypotheses in SEM models. 
  \emph{Structural Equation Modeling}, \bold{17}, 443-463.
  
  Van de Schoot, R., Strohmeier, D. (2011). 
  Testing informative hypotheses in SEM increases power: An 
  illustration contrasting classical. \emph{International Journal 
  of Behavioral Development}, \bold{35}, 180-190.
  
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
myConstraints1 <- ' GroupActive < GroupPassive < GroupNo '

iht(fit1.lm, myConstraints1)


# another way is to first fit the restricted model
fit.restr1 <- restriktor(fit1.lm, constraints = myConstraints1)

iht(fit.restr1)

\dontrun{ 
  # Or in matrix notation.
  Amat1 <- rbind(c(-1, 0,  1),
                 c( 0, 1, -1))
  myRhs1 <- rep(0L, nrow(Amat1)) 
  myNeq1 <- 0
  
  iht(fit1.lm, constraints = Amat1,
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
myConstraints2 <- ' group1 < group2 < group3 < group4 '

# compute F-test for hypothesis test Type A and compute the tail 
# probability based on the parametric bootstrap. We only generate 9 
# bootstrap samples in this example; in practice you may wish to 
# use a much higher number.
iht(fit2.lm, constraints = myConstraints2, type = "A", 
    boot = "parametric", R = 9)


# or fit restricted linear model
fit2.con <- restriktor(fit2.lm, constraints = myConstraints2)

iht(fit2.con)

\dontrun{ 
  # increasing means in matrix notation.
  Amat2 <- rbind(c(-1, 1, 0, 0),
                 c( 0,-1, 1, 0),
                 c( 0, 0,-1, 1))
  myRhs2 <- rep(0L, nrow(Amat2)) 
  myNeq2 <- 0

  iht(fit2.con, constraints = Amat2, rhs = myRhs2, neq = myNeq2, 
      type = "A", boot = "parametric", R = 9)
}            

## example 3: equality constraints only.
myConstraints3 <- ' group1 = group2 = group3 = group4 '

iht(fit2.lm, constraints = myConstraints3)

# or
fit3.con <- restriktor(fit2.lm, constraints = myConstraints3)
iht(fit3.con)


## example 4:
# combination of equality and inequality constraints.
myConstraints4 <- ' group1 = group2
                    group3 < group4 '

iht(fit2.lm, constraints = myConstraints4, type = "B", neq.alt = 1)

# fit resticted model and compute model-based bootstrapped 
# standard errors. We only generate 9 bootstrap samples in this 
# example; in practice you may wish to use a much higher number.
# Note that, a warning message may be thrown because the number of 
# bootstrap samples is too low.
fit4.con <- restriktor(fit2.lm, constraints = myConstraints4, 
                       se = "boot.model.based", B = 9)
iht(fit4.con, type = "B", neq.alt = 1)


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

iht(fit5.lm, constraints = myConstraints5)

# or
fit5.con <- restriktor(fit5.lm, constraints = ' AVE := X + 16.86137*X.Z; 
                                                AVE > 0 ')
iht(fit5.con)


\dontrun{
  # testing equality and/or inequality restrictions in SEM:
  
  #########################
  ### real data example ###
  #########################
  # Multiple group path model for facial burns example.
  
  # model syntax with starting values.
  burns.model <- 'Selfesteem ~ Age + c(m1, f1)*TBSA + HADS +
                             start(-.10, -.20)*TBSA  
               HADS ~ Age + c(m2, f2)*TBSA + RUM +
                      start(.10, .20)*TBSA '
  
  
  # constraints syntax
  burns.constraints <- 'f2 > 0  ; m1 < 0
                        m2 > 0  ; f1 < 0
                        f2 > m2 ; f1 < m1'
  
  # we only generate 2 bootstrap samples in this example; in practice
  # you may wish to use a much higher number. 
  # the double bootstrap was switched off; in practice you probably 
  # want to set it to "standard".
  example6 <- iht(model = burns.model, data = FacialBurns,
                  R = 2, constraints = burns.constraints,
                  double.bootstrap = "no", group = "Sex")
  
  example6
  
  ##########################
  ### artificial example ###
  ##########################
  # Simple ANOVA model with 3 groups (N = 20 per group)
  set.seed(1234)
  Y <- cbind(c(rnorm(20,0,1), rnorm(20,0.5,1), rnorm(20,1,1)))
  grp <- c(rep("1", 20), rep("2", 20), rep("3", 20))
  Data <- data.frame(Y, grp)
  
  #create model matrix
  fit.lm <- lm(Y ~ grp, data = Data)
  mfit <- fit.lm$model
  mm <- model.matrix(mfit)
  
  Y <- model.response(mfit)
  X <- data.frame(mm[,2:3])
  names(X) <- c("d1", "d2")
  Data.new <- data.frame(Y, X)
  
  # model
  model <- 'Y ~ 1 + a1*d1 + a2*d2'
  
  # fit without constraints
  fit <- lavaan::sem(model, data = Data.new)
  
  # constraints syntax: mu1 < mu2 < mu3
  constraints <- ' a1 > 0
                   a1 < a2 '
  
  # we only generate 10 bootstrap samples in this example; in practice
  # you may wish to use a much higher number, say > 1000. The double 
  # bootstrap is not necessary in case of an univariate ANOVA model.
  example7 <- iht(model = model, data = Data.new, 
                  start = lavaan::parTable(fit),
                  R = 10L, double.bootstrap = "no",
                  constraints = constraints)
  example7
  }
}