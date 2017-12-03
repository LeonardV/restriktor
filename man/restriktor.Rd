\name{restriktor}
\alias{restriktor}
\alias{conLM.lm}
\alias{conRLM.rlm}
\alias{conGLM.glm}
\alias{conMLM.mlm}

\title{Estimating linear regression models with (in)equality restrictions}
\description{Function \code{restriktor} estimates the parameters 
  of an univariate and a multivariate linear model (\code{lm}), a 
  robust estimation of the linear model (\code{rlm}) and a generalized 
  linear model (\code{glm}) subject to linear equality and linear 
  inequality restrictions. It is a convenience function. The real work 
  horses are the \code{conLM}, \code{conMLM}, \code{conRLM} and 
  the \code{conGLM} functions.}

\usage{
restriktor(object, constraints = NULL, \dots)

\method{conLM}{lm}(object, constraints = NULL, se = "standard", 
      B = 999, rhs = NULL, neq = 0L, mix.weights = "pmvnorm", 
      mix.bootstrap = 99999L, parallel = "no", ncpus = 1L, 
      cl = NULL, seed = NULL, control = list(), 
      verbose = FALSE, debug = FALSE, \dots)
      
\method{conRLM}{rlm}(object, constraints = NULL, se = "standard", 
       B = 999, rhs = NULL, neq = 0L, mix.weights = "pmvnorm", 
       mix.bootstrap = 99999L, parallel = "no", ncpus = 1L, 
       cl = NULL, seed = NULL, control = list(), 
       verbose = FALSE, debug = FALSE, \dots)
       
\method{conGLM}{glm}(object, constraints = NULL, se = "standard", 
       B = 999, rhs = NULL, neq = 0L, mix.weights = "pmvnorm", 
       mix.bootstrap = 99999L, parallel = "no", ncpus = 1L, 
       cl = NULL, seed = NULL, control = list(), 
       verbose = FALSE, debug = FALSE, \dots)

\method{conMLM}{mlm}(object, constraints = NULL, se = "none", 
       B = 999, rhs = NULL, neq = 0L, mix.weights = "pmvnorm", 
       mix.bootstrap = 99999L, parallel = "no", ncpus = 1L, 
       cl = NULL, seed = NULL, control = list(), 
       verbose = FALSE, debug = FALSE, \dots)
}

\arguments{
  \item{object}{a fitted linear model object of class "lm", "mlm",
  "rlm" or "glm". For class "rlm" only the loss function \code{bisquare} 
  is supported for now, otherwise the function gives an error.} 
  
  \item{constraints}{there are two ways to constrain parameters. 
  First, the constraint syntax consists of one or more text-based
  descriptions, where the syntax can be specified as a literal 
  string enclosed by single quotes. Only the names of \code{coef(model)}
  can be used as names. See details for more information. Note that
  objects of class "mlm" do not (yet) support this method. 
  
  Second, the constraint syntax consists of a matrix \eqn{R} (or a vector in 
  case of one constraint) and defines the left-hand side of the 
  constraint \eqn{R\theta \ge rhs}, where each row represents one 
  constraint. The number of columns needs to correspond to the 
  number of parameters estimated (\eqn{\theta}) by model. The rows 
  should be linear independent, otherwise the function gives an 
  error. For more information about constructing the matrix \eqn{R} and 
  \eqn{rhs} see details.}
  
  \item{se}{if "\code{standard}" (default), conventional standard errors 
  are computed based on inverting the observed augmented information 
  matrix. If "const", homoskedastic standard errors are computed. 
  If "\code{HC0}" or just "\code{HC}", heteroskedastic robust standard 
  errors are computed (a.k.a Huber White). The options "\code{HC1}", 
  "\code{HC2}", "\code{HC3}", "\code{HC4}", "\code{HC4m}", and 
  "\code{HC5}" are refinements of "\code{HC0}". For more details about 
  heteroskedastic robust standard errors see the \pkg{sandwich} 
  package. If "\code{boot.standard}", bootstrapped standard 
  errors are computed using standard bootstrapping. If "\code{boot.model.based}" 
  or "\code{boot.residual}", bootstrapped standard errors are computed 
  using model-based bootstrapping. If "\code{none}", no standard errors 
  are computed. Note that for objects of class "mlm" no standard errors 
  are available (yet).}
  \item{B}{integer; number of bootstrap draws for \code{se}. The 
  default value is set to 999. Parallel support is available.}
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
  \item{mix.weights}{if \code{"pmvnorm"} (default), the chi-bar-square 
  weights are computed based on the multivariate normal distribution 
  function with additional Monte Carlo steps. If \code{"boot"}, the 
  chi-bar-square weights are computed using parametric bootstrapping. 
  If \code{"none"}, no chi-bar-square weights are computed. The 
  weights are necessary in the \code{restriktor.summary} function 
  for computing the GORIC. Moreover, the weights are re-used in the 
  \code{\link{iht}} function for computing the p-value for the 
  test-statistic, unless the p-value is computed directly via bootstrapping.}
  \item{mix.bootstrap}{integer; number of bootstrap draws for 
  \code{mix.weights = "boot"}. The default value is set to 99999. 
  Parallel support is available.}
  \item{parallel}{the type of parallel operation to be used (if any). 
  If missing, the default is set "no".}
  \item{ncpus}{integer: number of processes to be used in parallel 
  operation: typically one would chose this to the number of 
  available CPUs.}
  \item{cl}{an optional parallel or snow cluster for use if 
  parallel = "snow". If not supplied, a cluster on the local machine 
  is created for the duration of the \code{restriktor} call.}
  \item{seed}{seed value.}
  \item{control}{a list of control arguments: 
    \itemize{
      \item \code{absval} tolerance criterion for convergence 
      (default = sqrt(.Machine$double.eps)). 
      \item \code{maxit} the maximum number of iterations for the 
      optimizer (default = 10000). 
      \item \code{tol} numerical tolerance value. Estimates smaller 
      than \code{tol} are set to 0.
    }
  }
  \item{verbose}{logical; if TRUE, information is shown at each 
  bootstrap draw.}
  \item{debug}{if TRUE, debugging information about the constraints
  is printed out.
  }
  \item{\ldots}{no additional arguments for now.}
}

\details{
  The constraint syntax can be specified in two ways. First as a 
  literal string enclosed by single quotes as shown below: 
  
  \preformatted{myConstraints <- '
    # 1. inequality constraints
      x1 > 0
      x1 < x2
    
    ! 2. equality constraints  
      x3 == x4; x4 == x5 '  
  }
  The variable names x1 to x5 refer to the corresponding regression
  coefficient. Thus, constraints are impose on regression coefficients
  and not on the data.
  
  Second, the above constraints syntax can also be written in 
  matrix/vector notation as:
  
  (The first column refers to the intercept, the remaining five
  columns refer to the regression coefficients x1 to x5.)
  \preformatted{myConstraints <-
    rbind(c(0, 0, 0,-1, 1, 0), #equality constraint x3 == x4
          c(0, 0, 0, 0,-1, 1), #equality constraint x4 == x5
          c(0, 1, 0, 0, 0, 0), #inequality constraint x1 > rhs
          c(0,-1, 1, 0, 0, 0)) #inequality constraint x1 < x2
  
# the length of rhs is equal to the number of myConstraints rows.      
myRhs <- c(0,0,0,0) 
    
# the first two rows should be considered as equality constraints
myNeq <- 2  
}
  
Blank lines and comments can be used in between the constraints, 
and constraints can be split over multiple lines. Both the 
hashtag (#) and the exclamation (!) characters can be used to 
start a comment. Multiple constraints can be placed on a single 
line if they are separated by a semicolon (;).

There can be three types of text-based descriptions in the constraints 
syntax:
  \enumerate{
      \item Equality constraints: The "\code{==}" operator can be 
      used to define equality constraints (e.g., \code{x1 == 1} or 
      \code{x1 == x2}).
      
      \item Inequality constraints: The "\code{<}" or "\code{>}" 
      operator can be used to define inequality constraints 
      (e.g., \code{x1 > 1} or \code{x1 < x2}).
      
      \item Newly defined parameters: The "\code{:=}" operator can 
      be used to define new parameters, which take on values that 
      are an arbitrary function of the original model parameters. 
      The function must be specified in terms of the parameter names 
      in \code{coef(model)} (e.g., \code{new := x1 + 2*x2}). By 
      default, the standard errors for these defined parameters are 
      computed by using the so-called Delta method.
  }

  Variable names of interaction effects in objects of class lm, 
  rlm and glm contain a semi-colon (:) between the variables. To impose 
  constraints on parameters of interaction effects, the semi-colon 
  must be replaced by a dot (.) (e.g., \code{x3:x4} becomes 
  \code{x3.x4}). In addition, the intercept variable names is shown 
  as "\code{(Intercept)}". To impose restrictions on the intercept 
  both parentheses must be replaced by a dot "\code{.Intercept.}" 
  (e.g.,\code{.Intercept. > 10}). Note: in most practical situations 
  we do not impose restrictions on the intercept because we do not 
  have prior knowledge about the intercept. Moreover, the sign of 
  the intercept can be changed arbitrarily by shifting the response 
  variable \eqn{y}.
  
  Each element can be modified using arithmetic operators. For example, 
  if \code{x2} is expected to be twice as large as \code{x1}, 
  then "\code{2*x2 == x1}". 
  
  If \code{constraints = NULL}, the unrestricted model is fitted.
}

\value{
An object of class restriktor, for which a print and a 
summary method are available. More specifically, it is a list 
with the following items:

\item{CON}{a list with useful information about the restrictions.}
\item{call}{the matched call.}
\item{timing}{how much time several tasks take.}
\item{parTable}{a parameter table with information about the 
 observed variables in the model and the imposed restrictions.}
\item{b.unrestr}{unrestricted regression coefficients.}
\item{b.restr}{restricted regression coefficients.}
\item{residuals}{restricted residuals.}
\item{wresid}{a working residual, weighted for "inv.var" weights 
only (rlm only)}
\item{fitted}{restricted fitted mean values.}
\item{weights}{(only for weighted fits) the specified weights.}
\item{wgt}{the weights used in the IWLS process (rlm only).}
\item{scale}{the robust scale estimate used (rlm only).}
\item{stddev}{a scale estimate used for the standard errors.}
\item{R2.org}{unrestricted R-squared.}
\item{R2.reduced}{restricted R-squared.}
\item{df.residual}{the residual degrees of freedom}
\item{s2.unrestr}{mean squared error of unrestricted model.}
\item{s2.restr}{mean squared error of restricted model.}
\item{loglik}{restricted log-likelihood.}
\item{Sigma}{variance-covariance matrix of unrestricted model.}
\item{constraints}{matrix with restrictions.}
\item{rhs}{vector of right-hand side elements.}
\item{neq}{number of equality restrictions.}
\item{wt.bar}{chi-bar-square mixing weights or a.k.a. level probabilities.}
\item{iact}{active restrictions.}
\item{converged}{did the IWLS converge (rlm only)?}
\item{iter}{number of iteration needed for convergence (rlm only).}
\item{bootout}{object of class boot. Only available if bootstrapped
standard errors are requested, else bootout = NULL.} 
\item{control}{list with control options.}
\item{model.org}{original model.}
\item{se}{as input. This information is needed in the summary 
function.}
\item{information}{observed information matrix with the inverted 
information matrix and the augmented information matrix as attributes.}
}

\references{
Schoenberg, R. (1997). Constrained Maximum Likelihood. \emph{Computational 
Economics}, \bold{10}, 251--266.

Shapiro, A. (1988). Towards a unified theory of inequality-constrained 
testing in multivariate analysis. \emph{International Statistical Review} 
\bold{56}, 49--62.

Silvapulle, M.J. and Sen, P.K. (2005). \emph{Constrained Statistical Inference}. 
Wiley, New York
}

\author{Leonard Vanbrabant and Yves Rosseel}


\seealso{ 
  \code{\link{iht}},
  \code{\link{goric}}
}

\examples{
## lm
## unrestricted linear model for ages (in months) at which an 
## infant starts to walk alone.

# prepare data
DATA1 <- subset(ZelazoKolb1972, Group != "Control")

# fit unrestricted linear model
fit1.lm <- lm(Age ~ -1 + Group, data = DATA1)

# the variable names can be used to impose restrictions on
# the corresponding regression parameters.
coef(fit1.lm)

# restricted linear model with restrictions that the walking 
# exercises would not have a negative effect of increasing the 
# mean age at which a child starts to walk. 
fit1.con <- restriktor(fit1.lm, constraints = ' GroupActive  < GroupPassive; 
                                                GroupPassive < GroupNo ')
summary(fit1.con)

\dontrun{ 
# Or in matrix notation.
myConstraints1 <- rbind(c(-1, 1, 0),
                        c( 0,-1, 1))
myRhs1 <- rep(0L, nrow(R1)) 
myNeq1 <- 0

fit1.con <- restriktor(fit1.lm, constraints = myConstraints1,
                       rhs = myRhs1, neq = myNeq1)
summary(fit1.con)
}            

#########################
## Artificial examples ##
#########################
library(MASS)

## mlm
# generate data
n <- 30
mu <- rep(0, 4)
Sigma <- matrix(5,4,4)
  diag(Sigma) <- c(10,10,10,10)
# 4 Y's.
Y <- mvrnorm(n, mu, Sigma)

# fit unrestricted multivariate linear model
fit.mlm <- lm(Y ~ 1)

# constraints
myConstraints2 <- diag(0,4)
  diag(myConstraints2) <- 1

# fit restricted multivariate linear model
fit2.con <- restriktor(fit.mlm, constraints = myConstraints2)

summary(fit2.con)


## rlm
# generate data
n <- 10
means <- c(1,2,1,3)
nm <- length(means)
group <- as.factor(rep(1:nm, each = n))
y <- rnorm(n * nm, rep(means, each = n))
DATA2 <- data.frame(y, group)

# fit unrestricted robust linear model
fit3.rlm <- rlm(y ~ -1 + group, data = DATA2, method = "MM")
coef(fit3.rlm)

## increasing means
myConstraints3 <- ' group1 < group2
                    group2 < group3
                    group3 < group4 '

# fit restricted robust linear model and compute 
# Huber-White (robust) standard errors.
fit3.con <- restriktor(fit3.rlm, constraints = myConstraints2, 
                       se = "HC0")
summary(fit3.con)

\dontrun{ 
## increasing means in matrix notation.
myConstraints3 <- rbind(c(-1, 1, 0, 0),
                        c( 0,-1, 1, 0),
                        c( 0, 0,-1, 1))
myRhs3 <- rep(0L, nrow(myConstraints3)) 
myNeq2 <- 0

fit3.con <- restriktor(fit3.rlm, constraints = myConstraints3,
                       rhs = myRhs2, neq = myNeq2, se = "HC0")
summary(fit3.con)
}            

## equality restrictions only.
myConstraints4 <- ' group1 == group2
                    group2 == group3
                    group3 == group4 '
                    
fit4.con <- restriktor(fit3.rlm, constraints = myConstraints4)
summary(fit4.con)


## combination of equality and inequality restrictions.
myConstraints5 <- ' group1 == group2
                    group3  < group4 '

# fit restricted model and compute model-based bootstrapped 
# standard errors. We only generate 9 bootstrap samples in this 
# example; in practice you may wish to use a much higher number.
fit5.con <- restriktor(fit3.rlm, constraints = myConstraints4, 
                       se = "boot.model.based", B = 9)
# an error is probably thrown, due to a too low number of bootstrap draws.
summary(fit5.con)

# restriktor can also be used to define effects using the := operator 
# and impose restrictions on them. For example, compute the average 
# effect (AVE) and impose the restriction AVE > 0.
# generate data
n <- 30
b0 <- 10; b1 = 0.5; b2 = 1; b3 = 1.5
X <- c(rep(c(0), n/2), rep(c(1), n/2))
set.seed(90) 
Z <- rnorm(n, 16, 5)
y <- b0 + b1*X + b2*Z + b3*X*Z + rnorm(n, 0, sd = 10) 
DATA3 = data.frame(cbind(y, X, Z))

# fit linear model with interaction
fit6.lm <- lm(y ~ X*Z, data = DATA3)

fit6.con <- restriktor(fit6.lm, constraints = ' AVE := X + 16.86137*X.Z; 
                                                AVE > 0 ')
summary(fit6.con)
}