\name{restriktor}
\alias{restriktor}
\alias{conLM.lm}
\alias{conRLM.rlm}

\title{Fitting linear regression models with equality and/or 
inequality restriktions}
\description{}

\usage{
restriktor(model, constraints, se = "standard", 
   rhs = NULL, neq = NULL, control = NULL, debug = FALSE, ...)
           
\method{restriktor}{restriktor}(model, constraints, se = "default", 
   rhs = NULL, neq = NULL, control = NULL, debug = FALSE, ...)
}
\arguments{
  \item{model}{a fitted linear model object of class lm, or a 
  fitted robust linear model object of class rlm.}
  \item{constraints}{there are two ways to constrain parameters. 
  First, the constraint syntax consists of one or more formula-like 
  expressions, where the syntax can be specified as a literal 
  string enclosed by single quotes. Only the names of coef(model) 
  can be used as constraints names. See details for more information. 
  
  Second, the constraints consists of a matrix R (or a vector in 
  case of one constraint) and defines the left-hand side of the 
  constraint \eqn{R\theta \ge rhs}, where each row represents one 
  constraint. The number of columns needs to correspond to the 
  number of parameters estimated (\eqn{\theta}) by model. The rows 
  should be linear independent, otherwise the function gives an 
  error. For more information about constructing the matrix R see
  details. 
  }
  \item{se}{if "standard" (default), conventional standard errors 
  are computed based on inverting the observed augmented information 
  matrix. If "const", homoskedastic standard errors are computed. 
  If "HC0" or just "HC", heteroskedastic robust standard errors 
  are computed (a.k.a Huber White). The options "HC1", "HC2", 
  "HC3", "HC4", "HC4m", and "HC5" are refinements of "HC0". 
  For more details see \code{\link{vcovHC}} in the \pkg{sandwich} 
  package. If "boot.standard" bootstrapped, standard errors are 
  computed using standard bootstrapping. If "boot.model.based", 
  bootstrapped standard errors are computed using model-based 
  bootstrapping. If "none", no standard errors are computed. 
  }
  \item{rhs}{vector on the right-hand side of the constraints; 
  \eqn{R\theta \ge rhs}. The length of this vector equals the 
  number of rows of the constraints matrix and consists of zeros 
  by default. Note: only used if constraints input is a matrix or 
  vector.
  }
  \item{neq}{integer (default = 0) treating the number of 
  constraints rows as equality constraints instead of inequality 
  constraints. For example, if neq = 2, this means that the first 
  two rows of the constraints matrix are treated as equality 
  constraints. Note: only used if constraints input is a matrix or 
  vector.
  }
  \item{control}{a list of control arguments; see details.
  }
  \item{debug}{if TRUE, debugging information is printed out.
  }
}

\details{
The constraint syntax consists of one or more formula-like 
expressions. The constraint syntax can be specified as a literal 
string enclosed by single quotes as shown below: 

\preformatted{myConstraints <- '
  # 1. inequality constraints
    x1 > 0
    x1 < x2
  
  ! 2. equality constraints  
    x3 == x4; x4 == x5    
'  
}


The above constraints syntax can also be written as:
\preformatted{myConstraints <-
  rbind(c(0, 0, 0,-1, 1, 0), #equality constraint x3 == x4
        c(0, 0, 0, 0,-1, 1), #equality constraint x4 == x5
        c(0, 1, 0, 0, 0, 0), #inequality constraint x1 > rhs
        c(0,-1, 1, 0, 0, 0)) #inequality constraint x1 < x2

The first column refers to the intercept, the remaining five
columns refer to the regression coefficients x1 to x5.

#rhs equal to the number of myConstraints rows.      
myRhs <- c(0,0,0,0) 
  
#the first two rows should be considered as equality constraints
myNeq <- 2  
}

Blank lines and comments can be used in between the constraints, 
and constraints can be split over multiple lines. Both the sharp 
(#) and the exclamation (!) characters can be used to start a 
comment. Multiple constraints can be placed on a single line if 
they are separated by a semicolon (;).

There can be three types of formula-like expressions in the 
constraints syntax:
\enumerate{

    \item Equality constraints: The "==" operator can be used to 
    define equality constraints. 
    
    \item Inequality constraints: The "<" or ">" operator can be 
    used to define inequality constraints.
    
    \item Newly defined parameters: The ":=" operator can be used 
    to define new parameters, which take on values that are an 
    arbitrary function of the original model parameters. The 
    function must be specified in terms of the parameter names 
    in coef(model). By default, the standard errors for these 
    defined parameters are computed by using the so-called Delta 
    method.
}

  Variable names of interaction effects in an lm object and an 
  rlm object are shown as semi-colons (:) between the variables. 
  In order to impose restriktions on parameters of interaction 
  effects, the semi-colon must be replaced by a dot (.). Thus, 
  an interaction effect "x3:x4" becomes "x3.x4". In addition, the 
  intercept is shown as "(Intercept)". In order to impose 
  restriktions on the intercept both parentheses must be replaced 
  by a dot (.Intercept.).
  
  Each element can be modified using arithmetic operators ("*", 
  "+"", "-" and "/"). For example, if x2 is expected to be twice 
  as large as x1, then "2*x2 == x1". 
  
  
}
\value{
An object of class conLM or conRLM, for which a print and a 
summary method are available. More specifically, it is a list 
with the following items:

\item{CON}{a list with useful information about the restriktions.}

\item{parTable}{a parameter table with information about the 
observed variables in the model and the imposed restriktions.}

\item{constraints}{a specification of the linear hypotheses. 
Depending on the input, the value is a formula-like expression 
or a matrix.}

\item{b.restr}{restrikted regression coefficients.}

\item{b.unrestr}{unrestrikter regression coefficients.}

\item{residuals}{restrikted residuals.}

\item{fitted}{restrikted fitted mean values.}

\item{weights}{(only for weighted fits) the specified weights.}

\item{df.residual}{the residual degrees of freedom (df). The df 
is corrected in case of equality constraints.}

\item{R2.org}{unrestrikted R-squared.}

\item{R2.reduced}{restrikted R-squared.}

\item{s2.unc}{mean squared error of unrestrikted model.}

\item{s2.restr}{mean squared error of restrikted model.}

\item{loglik}{restrikted log-likelihood.}

\item{Sigma}{variance-covariance matrix of unrestrikted model.}

\item{Amat}{matrix with constraints.}

\item{bvec}{vector of right-hand side elements.}

\item{meq}{number of equality constraints.}

\item{iact}{active restriktions.}

\item{bootout}{object of class boot. Only available if bootstrapped
standard errors are requested, else bootout = NULL.} 

\item{call}{the matched call.}
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

fit.lm <- lm(Age ~ Group, data = DATA)
summary(fit.lm)

# restrikted linear model with restriktions that the walking 
# exercises would not have a negative effect of increasing the 
# mean age at which a child starts to walk. 

fit.con <- restriktor(fit.lm, constraints = "Group2 > 0; Group2 < Group4")
summary(fit.con)

#########################
## Artificial examples ##
#########################
n <- 10
m <- c(1,2,1,3)

DATA <- data.frame(grp = as.factor(rep(1:nm, each = n)),
                   y = rnorm(n * nm, rep(m, each = n), 1))

# unrestrikted linear model
fit.lm <- lm(y ~ grp, data = DATA)

# increasing means
myConstraints <- ' grp2 > 0
  grp2 < grp3
  grp3 < grp4
'

fit1.con <- restriktor(fit.lm, constraints = myConstraints)
summary(fit1.con)


#myConstraints_mat <- rbind(c(0, 1, 0, 0),
#                           c(0,-1, 1, 0),
#                           c(0, 0,-1, 1))
#rhs = rep(0L, nrow(myConstraints_mat)) 
#neq = 0
            

# equality restriktions (all parameters equal)
fit2.con <- restriktor(fit.lm, constr = myConstraints,
                       rhs = rep(0, nrow(A1)), neq = nrow(A1))

#restrikting the 1st and the 4th parameter
#to their true values, and the 2nd and 3rd between them

}

