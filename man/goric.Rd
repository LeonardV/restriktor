\name{goric}
\alias{goric}
\alias{goric.default}
\alias{print.con_goric}
\alias{summary.con_goric}
\alias{coef.con_goric}

\title{Generalized Order-Restricted Information Criterion (Approximation) Weights}
\description{The \code{goric} function computes GORIC(A) weights, which are 
comparable to the Akaike weights. 
}
\usage{goric(object, ...)

\method{goric}{default}(object, ..., comparison = c("unconstrained", "complement", "none"), 
      VCOV = NULL, sample.nobs = NULL, type = "goric", bound = NULL, debug = FALSE)

\method{print}{con_goric}(x, digits = max(3, getOption("digits") - 4), \dots)

\method{summary}{con_goric}(object, brief = TRUE, digits = max(3, getOption("digits") - 4), \dots)

\method{coef}{con_goric}(object, \dots)
}

\arguments{
  \item{object}{an object containing the outcome of a statistical analysis.
  Currently, the following objects can be processed:
  \itemize{
    \item a fitted object of class restriktor. 
    \item a fitted unconstrained object of class \code{lm}, \code{rlm} or \code{glm}. 
    \item a numeric vector containing the unconstrained estimates resulting from any 
    statistical analysis.
    \item a fitted object of class \code{CTmeta}. 
  }
  }
  
  \item{x}{an object of class \code{con_goric}.}

  \item{\ldots}{this depends on the class of the object. If object is of class 
  restriktor, further objects of class restriktor can be passed. If object is of
  class \code{lm}, \code{rlm} or \code{glm}, the constraints can be passed. If
  object is of class \code{lavaan}, the standardized or unstandardized vcov can
  be used, using setting \code{standardized = TRUE}. If object is of 
  class numeric, the constraints can be passed. See details for more
  information.}
  
  \item{comparison}{if "\code{unconstrained}" (default) the unconstrained model is 
  included in the set of models. If "\code{complement}" then the restricted object 
  is compared against its complement. Note that the complement can only be computed 
  for one model/hypothesis at a time (for now). If "\code{none}" the model is only compared
  against the models provided by the user.}
  
  \item{VCOV}{variance-coviance matrix. Only needed if object is of class numeric and
  \code{type = "gorica"}.}
  
  \item{sample.nobs}{not used for now.}
  
  \item{type}{if \code{"goric"} (default), the generalized order-restricted 
  information criterion value is computed. If \code{"gorica"} the 
  log-likihood is computed using the multivariate normal distribution 
  function. }
  
  \item{bound}{not used yet.}
  
  \item{digits}{the number of significant digits to use when printing.}
  
  \item{debug}{if TRUE, debugging information is printed out.}
  
  \item{brief}{if FALSE, an extended overview is printed.}
}

\details{The GORIC(A) values themselves are not interpretable and 
the interest lie in their differences. The GORIC(A) weights reflect 
the support of each hypothesis in the set. To compare two hypotheses 
(and not one to the whole set), one can examine the ratio of the 
two corresponding GORIC(A) weights. To avoid selecting a weakly supported 
hypothesis as the best one, the unconstrained hypothesis is usually 
included as safeguard.

In case of one order-constrained hypothesis, say H1, the complement 
Hc can be computed as competing hypothesis. The complement is defined 
as Hc = not H1.

If the object(s) is of class restriktor the constraints are automatically extracted.
Otherwise, the constraint syntax can be parsed via the \ldots. If the object is 
an unconstrained model of class \code{lm}, \code{rlm} or \code{glm}, then the 
constraints can be specified in two ways, see \code{\link{restriktor}}. Note 
that if the constraints are written in matrix notation, then the constraints
for each model/hypothesis is put in a named list. For example, 
h1 <- list(constraints = "x1 > 0", rhs = 0, neq = 0). The \code{rhs} and \code{neq}
are not required if they are equal to 0. If \code{type = "gorica"},
then the object might be a (named) numeric vector. The constraints can again be 
specified in two ways, see \code{\link{restriktor}}. For examples, see below.

To determine the penalty term values, the chi-bar-square weights (a.k.a. level
probabilities) must be computed. If \code{"mix.weights = "pmvnorm" "} (default), 
the chi-bar-square weights are computed based on the multivariate normal distribution 
function with additional Monte Carlo steps. If \code{"mix.weights = "boot" "}, the 
chi-bar-square weights are computed using parametric bootstrapping (see \code{\link{restriktor}}).
}

\value{The function returns a dataframe with the log-likelihood,
penalty term, GORIC(A) values and the GORIC(A) weights. Furthermore, a dataframe
is returned with the relative GORIC(A) weights. 
}

\references{
Kuiper, R.M., Hoijtink, H., and Silvapulle, M.J. (2011). An Akaike-type
information criterion for model selection under inequality constraints.
\emph{Biometrika}, \bold{98}, 2, 495--501.

Vanbrabant, L. and Kuiper, R. (2020). Evaluating a theory-based hypothesis against 
its complement using an AIC-type information criterion with an application to 
facial burn injury. Psychological Methods. 
}

\author{Leonard Vanbrabant and Rebecca Kuiper}


\examples{
library(MASS)
## lm
## unrestricted linear model for ages (in months) at which an 
## infant starts to walk alone.

# prepare data
DATA <- subset(ZelazoKolb1972, Group != "Control")
  
# fit unrestrikted linear model
fit1.lm <- lm(Age ~ Group, data = DATA)

# some artificial restrictions
fit1.con <- restriktor(fit1.lm, constraints = "GroupPassive > 0;  GroupPassive < GroupNo")
fit2.con <- restriktor(fit1.lm, constraints = "GroupPassive > 0;  GroupPassive > GroupNo")
fit3.con <- restriktor(fit1.lm, constraints = "GroupPassive == 0; GroupPassive < GroupNo")
fit4.con <- restriktor(fit1.lm) # unrestricted model

goric(fit1.con, fit2.con, fit3.con, fit4.con)

# fit1.con versus the complement
goric(fit1.con, comparison = "complement")


## GORICA
# generate data
n <- 10
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + x1 + x2 + rnorm(n)
# fit unconstrained linear model
fit.lm <- lm(y ~ x1 + x2)

# extract unconstrained estimates
est <- coef(fit.lm)
# unconstrained variance-covariance matrix
VCOV <- vcov(fit.lm)

## constraint syntax (character)
h1 <- "x1 > 0"
h2 <- "x1 > 0; x2 > 0"
# use fitted unconstrained linear model
out <- goric(fit.lm, h1, h2, type = "gorica")
# use unconstrained estimates
out <- goric(est, VCOV = VCOV, h1, h2, type = "gorica")

## constraint syntax (matrix notation)
h1 <- list(constraints = c(0,1,0))
h2 <- list(constraints = rbind(c(0,1,0), c(0,0,1)))
out <- goric(fit.lm, h1, h2, type = "gorica")
out <- goric(est, VCOV = VCOV, h1, h2, type = "gorica")


## mlm
# generate data
n <- 30
mu <- c(1,2,3,4)
Sigma <- matrix(5,4,4)
  diag(Sigma) <- c(10,10,10,10)
# 4 Y's.
Y <- mvrnorm(n, mu, Sigma)

# fit unrestricted multivariate linear model
fit2.mlm <- lm(Y ~ 1)

# constraints
myConstraints2 <- rbind(c(-1,1,0,0), c(0,-1,1,0), c(0,0,-1,1))

# fit restricted multivariate linear model
fit5.con <- restriktor(fit2.mlm, constraints = myConstraints2)
}
