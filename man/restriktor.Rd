\name{restriktor}
\alias{restriktor}

\title{Fitting linear regression models with equality and/or inequality constraints}
\description{}

\usage{
restriktor(model, constraints, se = "standard", rhs = NULL, neq = NULL, 
           control = NULL, debug = FALSE, ...)
           
\method{restriktor}{restriktor}(model, constraints, se = "default", rhs = NULL, 
           neq = NULL, control = NULL, debug = FALSE, ...)
}
\arguments{
  \item{model}{A fitted linear model object of class lm, or a fitted robust linear 
  model object of class rlm.}
  \item{constraints}{There are two ways to constrain parameters. First, the 
  constraint syntax consists of one or more formula-like expressions, where the 
  syntax can be specified as a literal string enclosed by single quotes. Only
  the names of coef(model) can be used as constraints names. See details for 
  more information. 
  
  Second, the constraints consists of a matrix R (or a vector in case of one 
  constraint) and defines the left-hand side of the constraint \eqn{R\theta \ge rhs}, 
  where each row represents one constraint. The number of columns needs to 
  correspond to the number of parameters estimated (\eqn{\theta}) by model. The 
  rows should be linear independent, otherwise the function gives an error. 
  }
  \item{se}{If "standard" (default), conventional standard errors are computed 
  based on inverting the observed information matrix. if "const", homoskedastic
  standard errors are computed. If "HC0", heteroskedastic robust standard errors 
  are computed (a.k.a Huber White). The options "HC1", "HC2", "HC3", "HC4", 
  "HC4m", and "HC5" are refinements of "HC0". For more details see \code{\link{vcovHC}} 
  in the \pkg{sandwich} package. If "boot.standard" bootstrapped, standard errors 
  are computed using standard bootstrapping. If "boot.model.based", bootstrapped 
  standard errors are computed using model-based bootstrapping. If "none", no 
  standard errors are computed. 
  }
  \item{rhs}{Vector on the right-hand side of the constraints; \eqn{R\theta \ge rhs}. 
  The length of this vector equals the number of rows of the constraints matrix
  and consists of zeros by default. Note: only used if constraints input is a 
  matrix or vector.
  }
  \item{neq}{Integer (default = 0) treating the number of constraints rows 
  as equality constraints instead of inequality constraints. For example, if 
  neq = 2, this means that the first two rows of the constraints matrix are 
  treated as equality constraints. Note: only used if constraints input is a 
  matrix or vector.
  }
  \item{control}{A list of control arguments; see details.
  }
  \item{debug}{If TRUE, debugging information is printed out.
  }
}

\details{
The constraint syntax consists of one or more formula-like expressions. The 
constraint syntax can be specified as a literal string enclosed by single quotes 
as shown below: 

\preformatted{myConstraints <- '
  # 1. inequality constraints
    x1 > 0
    x1 < x2
  
  ! 2. equality constraints  
    x3 == x4
'  
}

Blank lines and comments can be used in between the formulas, and formulas can 
be split over multiple lines. Both the sharp (#) and the exclamation (!) 
characters can be used to start a comment. Multiple formulas can be placed on a 
single line if they are separated by a semicolon (;).

There can be three types of formula-like expressions in the constraints syntax:
\enumerate{

    \item Equality constraints: The "==" operator can be used to define equality
    constraints. 
    
    \item Inequality constraints: The "<" or ">" operator can be used to define
    inequality constraints.
    
    \item Newly defined parameters: The ":=" operator can be used to define new
    parameters, which take on values that are an arbitrary function of the 
    original model parameters. The function must be specified in terms of the 
    parameter names in coef(model). By default, the standard errors for these 
    defined parameters are computed by using the so-called Delta method.
}

  Variable names of interaction effects in an lm object and an rlm object are 
  shown as semi-colons (:) between the variables. In order to impose restriktions 
  on interaction variables, the semi-colon must be replaced by a dot (.). Thus,
  an interaction effect "x3:x4" becomes "x3.x4". In addition, the intercepts is 
  shown as "(Intercept)". In order to impose restriktions on intercepts both 
  parentheses must be replaced by a dot (.Intercept.).
  
  Each element can be modified using arithmetic operators ("*", "+"", "-" and "/"). 
  For example, if x2 is expected to be twice as large as x1, "2*x2 == x1". 
}
\value{An object of class conLM or conRLM, for which a print and a summary method
are available.}

\references{
Schoenberg, R. (1997). Constrained Maximum Likelihood. \emph{Computational 
Economics}, \bold{10}, 251--266.
}

\author{Leonard Vanbrabant and Yves Rosseel}


\seealso{ 
  \code{\link{solve.QP}}, 
  \code{\link{conTest}},
  \code{\link{vcovHC}}
}

\examples{
# unrestricted linear model for ages (in months) at which an infant starts to 
# walk alone.
DATA <- ZelazoKolb1972
idx <- which(DATA$Group == 3)
DATA <- DATA[-idx, ]
DATA$Group <- factor(DATA$Group)

fit.lm <- lm(Age ~ Group, data = DATA)
summary(fit.lm)

# restrikted linear model with restriktions that the walking exercises would 
# not have a negative effect of increasing the mean age at which a child starts 
# to walk. 

fit.con <- restriktor(fit.lm, constraints = "Group2 > 0; Group2 < Group4")
summary(fit.con)

#########################
## Artificial examples ##
#########################
n <- 10
m <- c(1,2,1,5)
nm <- length(m)
DATA <- data.frame(grp = as.factor(rep(1:nm, each = n)),
                   y = rnorm(n * nm, rep(m, each = n), 1))

# unrestricted linear model
fit.lm <- lm(y ~ grp, data = DATA)

# increasing means
A1 <- rbind(c(0, -1, 1, 0, 0),
            c(0, 0, -1, 1, 0),
            c(0, 0, 0, -1, 1))
             
fit1.con <- restriktor(fit.lm, constraints = A1,
                      rhs = rep(0L, nrow(A1)), neq = 0)
summary(fit.con)

# equality constraints (all parameters equal)
fit2.con <- restriktor(fit.lm, constr = A1,
                       rhs = rep(0, nrow(A1)), neq = nrow(A1))

# constraining the 1st and the 4th parameter
# to their true values, and the 2nd and 3rd between them
A3 <- rbind(c(0, 1, 0, 0, 0),
            c(0,-1, 1, 0, 0),
            c(0, 0, -1, 0, 1),
            c(0,-1, 0, 1, 0),
            c(0, 0, 0, -1, 1),
            c(0, 0, 0, 0, 1))

fit3.con <- restriktor(fit.lm, constraints = A3,
                       rhs = c(1,rep(0,4),5), neq = c(TRUE,rep(FALSE,4),TRUE))

}

\keyword{models}
