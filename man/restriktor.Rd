\name{restriktor}
\alias{restriktor}
\alias{conLM.lm}
\alias{conRLM.rlm}

\title{Estimating linear regression models with equality and/or 
inequality restriktions}
\description{Function \code{restriktor} estimates the parameters of a linear model (\code{lm}) 
  or robust linear model (\code{rlm}) subject to equality and/or inequality 
  restriktions. The function \code{restriktor} is a convenience function that can 
  take a symbolic description as constraint input and translates it to right 
  constraint matrices and vectors. The real work horses are the \code{conLM} and 
  the \code{conRLM} functions for linear and robust linear models respectively.}

\usage{
restriktor(model, constraints, se = "standard", 
           rhs = NULL, neq = NULL, control = NULL, debug = FALSE, ...)
\method{conLM}{lm}(model, constraints, se = "default", 
      rhs = NULL, neq = NULL, control = NULL, ...)
\method{conRLM}{rlm}(model, constraints, se = "default", 
       rhs = NULL, neq = NULL, control = NULL, ...)
}
\arguments{
  \item{model}{a fitted linear model object of class lm, or a fitted robust 
  linear model object of class rlm. For class rlm, both \code{M} and
  \code{MM} estimation are supported. For the loss function only the 
  \code{bisquare} is supported for now, otherwise the function gives an error.} 
  
  \item{constraints}{there are two ways to constrain parameters. 
  First, the constraint syntax consists of one or more symbolic
  descriptions, where the syntax can be specified as a literal 
  string enclosed by single quotes. Only the names of coef(model) 
  can be used as constraints names. See details for more information. 
  
  Second, the constraints consists of a matrix \eqn{R} (or a vector in 
  case of one constraint) and defines the left-hand side of the 
  constraint \eqn{R\theta \ge rhs}, where each row represents one 
  constraint. The number of columns needs to correspond to the 
  number of parameters estimated (\eqn{\theta}) by model. The rows 
  should be linear independent, otherwise the function gives an 
  error. For more information about constructing the matrix \eqn{R} 
  see details. 
  }
  \item{se}{if "\code{standard}" (default), conventional standard errors 
  are computed based on inverting the observed augmented information 
  matrix. If "const", homoskedastic standard errors are computed. 
  If "\code{HC0}" or just "\code{HC}", heteroskedastic robust standard 
  errors are computed (a.k.a Huber White). The options "\code{HC1}", 
  "\code{HC2}", "\code{HC3}", "\code{HC4}", "\code{HC4m}", and 
  "\code{HC5}" are refinements of "\code{HC0}". For more details about 
  heteroskedastic robust standard errors see \code{\link{vcovHC}} in 
  the \pkg{sandwich} package. If "\code{boot.standard}", bootstrapped standard 
  errors are computed using standard bootstrapping. If "\code{boot.model.based}" 
  or "\code{boot.residual}", bootstrapped standard errors are computed 
  using model-based bootstrapping. If "\code{none}", no standard errors 
  are computed. 
  }
  \item{rhs}{vector on the right-hand side of the constraints; 
  \eqn{R\theta \ge rhs}. The length of this vector equals the 
  number of rows of the constraints matrix \eqn{R} and consists of zeros 
  by default. Note: only used if constraints input is a matrix or 
  vector.
  }
  \item{neq}{integer (default = 0) treating the number of 
  constraints rows as equality constraints instead of inequality 
  constraints. For example, if \code{neq = 2}, this means that the first 
  two rows of the constraints matrix \eqn{R} are treated as equality 
  constraints. Note: only used if constraints input is a matrix or 
  vector.
  }
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
  \item{debug}{if TRUE, debugging information about the constraints
  is printed out.
  }
  \item{\ldots}{not used for now.}
}

\details{
  \code{restriktor} is a convenience function that takes a symbolic
  description as constraint input and translates it to the constraint
  matrix/vector \eqn{R}, the \eqn{rhs} and \eqn{neq}. The real work 
  horses are the \code{conLM} and the \code{conRLM} functions. The 
  constraint syntax can be specified in two ways. First as a 
  literal string enclosed by single quotes as shown below: 
  
  \preformatted{myConstraints <- '
    # 1. inequality constraints
      x1 > 0
      x1 < x2
    
    ! 2. equality constraints  
      x3 == x4; x4 == x5    
  '  
  }
  
  Second, the above constraints syntax can also be written in 
  matrix/vector notation as:
  \preformatted{myConstraints <-
    rbind(c(0, 0, 0,-1, 1, 0), #equality constraint x3 == x4
          c(0, 0, 0, 0,-1, 1), #equality constraint x4 == x5
          c(0, 1, 0, 0, 0, 0), #inequality constraint x1 > rhs
          c(0,-1, 1, 0, 0, 0)) #inequality constraint x1 < x2
  
  The first column refers to the intercept, the remaining five
  columns refer to the regression coefficients x1 to x5.
  
  #rhs is equal to the number of myConstraints rows.      
  myRhs <- c(0,0,0,0) 
    
  #the first two rows should be considered as equality constraints
  myNeq <- 2  
  }
  
  Blank lines and comments can be used in between the constraints, 
  and constraints can be split over multiple lines. Both the sharp 
  (#) and the exclamation (!) characters can be used to start a 
  comment. Multiple constraints can be placed on a single line if 
  they are separated by a semicolon (;).
  
  There can be three types of symoblic descriptions in the 
  constraints syntax:
  \enumerate{
  
      \item Equality constraints: The "\code{==}" operator can be used to 
      define equality constraints (e.g., \code{x1 == 1} or \code{x1 == x2}).
      
      \item Inequality constraints: The "\code{<}" or "\code{>}" operator can be 
      used to define inequality constraints (e.g., \code{x1 > 1} or \code{x1 < x2}).
      
      \item Newly defined parameters: The "\code{:=}" operator can be used 
      to define new parameters, which take on values that are an 
      arbitrary function of the original model parameters. The 
      function must be specified in terms of the parameter names 
      in \code{coef(model)} (e.g., \code{new := x1 + 2*x2}). By default, the 
      standard errors for these defined parameters are computed by 
      using the so-called Delta method.
  }

  Variable names of interaction effects in an lm object and an 
  rlm object are shown as semi-colons (:) between the variables. 
  In order to impose restriktions on parameters of interaction 
  effects, the semi-colon must be replaced by a dot (.). Thus, 
  an interaction effect "\code{x3:x4}" becomes "\code{x3.x4}" 
  (e.g., \code{x3.x4 > 0}). In addition, the intercept is shown as 
  "\code{(Intercept)}". In order to impose restriktions on the intercept 
  both parentheses must be replaced by a dot "\code{.Intercept.}" 
  (e.g.,\code{.Intercept. > 10}). Note: in most practical situations 
  we do not impose restriktions on the intercept because we do not have 
  prior knowledge about the intercept. Moreover, the sign of the 
  intercept can be changed arbitrarily by shifting the response 
  variable \eqn{y}.
  
  Each element can be modified using arithmetic operators. For example, 
  if \code{x2} is expected to be twice as large as \code{x1}, 
  then "\code{2*x2 == x1}". 
}
\value{
An object of class conLM or conRLM, for which a print and a 
summary method are available. More specifically, it is a list 
with the following items:

\item{CON}{a list with useful information about the restriktions.}

\item{parTable}{a parameter table with information about the 
observed variables in the model and the imposed restriktions.}

\item{b.restr}{restrikted regression coefficients.}

\item{b.unrestr}{unrestrikted regression coefficients.}

\item{residuals}{restrikted residuals.}

\item{fitted}{restrikted fitted mean values.}

\item{weights}{(only for weighted fits) the specified weights.}

\item{df.residual}{the residual degrees of freedom. Note: the 
df is corrected in case of equality constraints.}

\item{R2.org}{unrestrikted R-squared.}

\item{R2.reduced}{restrikted R-squared.}

\item{s2.unc}{mean squared error of unrestrikted model.}

\item{s2.restr}{mean squared error of restrikted model.}

\item{loglik}{restrikted log-likelihood.}

\item{Sigma}{variance-covariance matrix of unrestrikted model.}

\item{Amat}{matrix with restriktions.}

\item{bvec}{vector of right-hand side elements.}

\item{meq}{number of equality restriktions.}

\item{iact}{active restriktions.}

\item{bootout}{object of class boot. Only available if bootstrapped
standard errors are requested, else bootout = NULL.} 

\item{call}{the matched call.}

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