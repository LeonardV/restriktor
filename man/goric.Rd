\name{goric}
\alias{goric}
\alias{goric.default}
\alias{goric.lm}
\alias{goric.numeric}
\alias{goric.lavaan}
\alias{goric.CTmeta}
\alias{goric.rma}
\alias{print.con_goric}
\alias{summary.con_goric}
\alias{coef.con_goric}

\title{Generalized Order-Restricted Information Criterion (Approximation) Weights}
\description{The \code{goric} function computes GORIC(A) weights, which are 
comparable to the Akaike weights. 
}

\usage{goric(object, ...)

\method{goric}{default}(object, ..., hypotheses = NULL,
      comparison = NULL, type = "goric",
      VCOV = NULL, sample_nobs = NULL,
      penalty_factor = 2, Heq = FALSE, control = list(), debug = FALSE)

\method{goric}{lm}(object, ..., hypotheses = NULL,
      comparison = NULL, type = "goric",
      missing = "none", auxiliary = c(), emControl = list(),
      debug = FALSE)

\method{goric}{numeric}(object, ..., hypotheses = NULL,
      comparison = NULL, type = "gorica", 
      VCOV = NULL, sample_nobs = NULL,
      debug = FALSE)

\method{goric}{lavaan}(object, ..., hypotheses = NULL,
      comparison = NULL, type = "gorica",
      standardized = FALSE, debug = FALSE)

\method{goric}{CTmeta}(object, ..., hypotheses = NULL,
      comparison = NULL, type = "gorica", 
      sample_nobs = NULL, debug = FALSE)

\method{goric}{rma}(object, ..., hypotheses = NULL,
      comparison = NULL, type = "gorica",
      VCOV = NULL, sample_nobs = NULL, 
      debug = FALSE)


\method{print}{con_goric}(x, digits = max(3, getOption("digits") - 4), \dots)

\method{summary}{con_goric}(object, brief = TRUE, digits = max(3, getOption("digits") - 4), \dots)
      
\method{coef}{con_goric}(object, \dots)
}

\arguments{
  \item{object}{an object containing the outcome of an unconstrained statistical analysis.
  Currently, the following objects can be processed:
  \itemize{
    \item a fitted unconstrained object of class \code{aov}, \code{lm}, \code{glm}, \code{rlm} or \code{mlm}. 
    \item a numeric (named) vector containing the unconstrained estimates resulting from any 
    statistical analysis.
    \item a fitted object of class \code{lavaan}. See examples on how to specify the hypotheses.
    \item a fitted object of class \code{CTmeta}.
    \item a fitted object of class \code{rma}.
    \item a fitted object of class \code{nlmerMod}.
    \item a fitted object of class \code{glmerMod}.
    \item a fitted object of class \code{lmerMod}.
  }
  }
  
  \item{x}{an object of class \code{con_goric}.}

  \item{\ldots}{this depends on the class of the object. If object is of class 
  \code{lavaan}, the standardized or unstandardized vcov can be used, using 
  setting \code{standardized = TRUE}. See details for more information.
  
  \emph{Options for calculating the chi-bar-square weights:}

  Parameters passed to the truncated multivariate normal distribution. By default, 
  restriktor (i.e. \code{con_weights_boot} function) uses no truncation points 
  for calculating the chi-bar-square weights, which renders to the multivariate 
  normal distribution. See the manual page of the \code{rtmvnorm} function from 
  the \pkg{rtmvnorm} to see how to specify a truncated mvnorm distribution and 
  the possible arguments.}
  
  \item{hypotheses}{a named list; Please note that the hypotheses argument in the given 
  context serves the same purpose as the constraints argument utilized in the 
  restriktor function. The distinction between them is solely semantic.
  
  There are two ways to constrain parameters. 
  First, the hypothesis syntax consists of one or more text-based
  descriptions, where the syntax can be specified as a literal 
  string enclosed by single quotes. Only the names of \code{coef(model)} or
  \code{names(vector)} can be used as names. See details for more information. 
  Note that objects of class "mlm" do not (yet) support this method.
  % TO DO ik gebruik nu rijnamen van vcov en dan moet gebruiker evt : vervangen door .. 
  %       Werkt dat? Zo ja, dan kan deze opmerking gewijzigd worden.
  
  Second, the hypothesis syntax consists of a matrix \eqn{R} (or a vector in 
  case of one constraint) and defines the left-hand side of the constraint 
  \eqn{R\theta \ge rhs}, where each row represents one constraint. The number of 
  columns needs to correspond to the number of parameters estimated (\eqn{\theta}) 
  by model. The rows should be linear independent, otherwise the function gives an 
  error. For more information about constructing the matrix \eqn{R} and \eqn{rhs} 
  see details.}
  
  \item{comparison}{The default behavior depends on the number of user-specified 
  order-restricted hypotheses. If a single hypothesis is specified, it is compared 
  against its complement by default. When multiple order-restricted hypotheses are 
  specified, the unconstrained model is added as a safeguard to the set of hypotheses.

  These default settings can be overridden. Use "\code{unconstrained}" to include 
  the unconstrained model in the set of models. Use "\code{complement}" to compare 
  the order-restricted object against its complement; note that the complement can 
  only be computed for one model/hypothesis at a time for now. If "\code{none}" 
  is chosen, the model is compared only against those provided by the user.}
  
  \item{type}{if \code{"goric"}, the generalized order-restricted 
  information criterion value is computed. 
  If \code{"gorica"}, a proxi for the log-likihood is calculated by assuming that 
  the parameter estimates (of interest) are multivariate normally distributed. 
  Then, either the unbiased covariance matrix of the estimates (based on 'N'; e.g., for an lm object) is used
  or the biased one obtained via the vcov function (e.g., for a lavaan object). 
  In the latter case, a corresponding message will be printed. 
  If \code{"goricc"} or \code{"goricac"}, 
  then the small-sample-size correction is applied.}
  
  \item{VCOV}{variance-covariance matrix. Only needed if object is of class 'numeric' (and, therefore,
  \code{type = "gorica"} or \code{type = "goricac"}).
  For many fit objects, one can use the function vcov() to obtain the VCOV. 
  This renders the (biased) restricted sample covariance matrix and not the unbiased sample covariance matrix (based on 'N').
  To obtain the unbiased one, one should make an adjustment based on samples size and the rank of the model (i.e., number of predictors including intercept); 
  as shwon in an example below.
  }
  
  \item{sample_nobs}{the number of observations, only needed if \code{type = "goricac"}. Note that, 
  for many fit objects and/or if \code{type = "goricc"},
  the number of observations are inherited from the fitted object.
  In case a vector is given, it is assumed that those are group sizes; 
  then, their sum is used as the number of observations.}
  
  \item{penalty_factor}{The penalty factor adjusts the penalty in the GORIC(A) formula
  (GORIC(A) = -2 x log-likelihood + penalty_factor x penalty). 
  By default, GORIC(A) uses a penalty factor of 2, but penalty factor allows this 
  to be customized. Higher values of penalty factor place a stronger emphasis on 
  model simplicity, helping to prevent overfitting by penalizing complex models 
  more heavily.}
  
  \item{Heq}{If \code{TRUE}, the null hypothesis is added to the set of models. 
  This means that all inequality constraints are replaced with equality constraints, 
  effectively testing the hypothesis that the parameters satisfy exact equality 
  rather than inequality. Defaults to \code{FALSE}.}
  
  \item{missing}{the default setting for objects of class "lm" is listwise:
  all cases with missing values are removed from the data before the analysis. This
  is only valid if the data are missing completely at random (MCAR). Another 
  option is to use "two.stage". In this approach, missing data are imputed using 
  an EM algorithm. However, we cannot use the complete data as input for futher 
  analyses, because the resulting complete data variance-covariance matrix will 
  not be correct. Therefore, we compute the correct aymptotic covariance (Savalei and Bentler, 2009)
  and use it as input for the \code{goric.numeric} function to compute a GORICA(C)
  value. Note that, the parameter estimates are also recomputed using the complete data.} 
  
  \item{auxiliary}{Vector. The inclusion of auxiliary variables can improve the
  imputation model. These auxiliary variables are not part of the target model.}
  
  \item{emControl}{a list of control arguments for the \code{emnorm} function
  from the \pkg{norm} package.}
  
  \item{standardized}{if TRUE, standardized parameter estimates are used.
  % TO DO Only used for class lavaan?
  }
  
  \item{digits}{the number of significant digits to use when printing.}
  
  \item{debug}{if TRUE, debugging information is printed out.}
  
  Control options for calculating the chi-bar-square weights:
  \item{control}{
    \itemize{
      \item \code{chunk_size} integer; the chi-bar-square weights are computed for samples
      of size \code{chunk_size = 5000L}. This process is repeated iteratively until the 
      weights converges (see \code{convergenge_crit}) or the maximum is reached, i.e., 
      \code{mix_weights_bootstrap_limit}.
      
      \item \code{mix_weights_bootstrap_limit} integer; maximum number of bootstrap draws.
      The default value is set to 1e5.
      
      \item \code{convergence_crit} the convergence criterion for the iterative bootstrap
      process. Default is 1e-03.
    }
  }
  \item{brief}{if TRUE, a short overview is printed.}
}

\details{The GORIC(A) can be used to evaluate theory-based, informative hypotheses,
  that is, quantify the relative support of theory-based, informative hypotheses. 
  The GORIC(A) values themselves are not interpretable and 
  the interest lie in their differences. The GORIC(A) weights reflect 
  the support of each hypothesis in the set. To compare two hypotheses 
  (and not one to the whole set), one can examine the ratio of the 
  two corresponding GORIC(A) weights. To avoid selecting a weakly supported 
  hypothesis as the best one, the unconstrained hypothesis is usually 
  included as safeguard.
  
  In case of one order-constrained hypothesis, say H1, the complement 
  Hc can be used as competing hypothesis. The complement is defined 
  as 'Hc: not H1'.
  
  The hypothesis syntax can be parsed via the hypotheses argument. 
  If the object is an unconstrained model of class 
  \code{aov}, \code{lm}, \code{glm}, \code{rlm}, or \code{mlm}, 
  then the hypotheses can be specified in two ways, see \code{\link{restriktor}}. Note 
  that if the hypotheses are written in matrix notation, then the hypotheses
  for each model/hypothesis is put in a named list with specific names constraints, rhs, 
  and neq. For example with three parameters x1, x2, x3, and x1 > 0: 
  list(model1 = list(constraints = rbind(c(1, 0, 0)), rhs = 0, neq = 0))). The \code{rhs} 
  and \code{neq} are not required if they are equal to 0. If \code{type = "gorica"},
  then the object might be a (named) numeric vector. The hypotheses can again be 
  specified in two ways, see \code{\link{restriktor}}. For examples, see below.
  
  To determine the penalty term values, the chi-bar-square weights (a.k.a. level
  probabilities) must be computed. If \code{"mix_weights = "pmvnorm" "} (default), 
  the chi-bar-square weights are computed based on the multivariate normal distribution 
  function with additional Monte Carlo steps. If \code{"mix_weights = "boot" "}, the 
  chi-bar-square weights are computed using parametric bootstrapping (see \code{\link{restriktor}}).
  
  
  The "two.stage" approach for missing data uses the EM algorithm from the 
  \code{norm} package. The response variables are assumed to be jointly 
  normal. In practice, missing-data procedures designed for variables that are 
  normal are sometimes applied to variables that are not. Binary and ordinal 
  variables are sometimes imputed under a normal model, and the imputed values 
  may be classified or rounded. This is also how restriktor handles (ordered)
  factors for now. 
  
  A better strategy (not implemented yet) would be to convert (ordered) factors 
  into a pair of dummy variables. If the (ordered) factors have missing values, 
  the dummy variables could be included as columns of Y and imputed, but then you 
  have to decide how to convert the continuously distributed imputed values for 
  these dummy codes back into categories. 
  
  ### Note on not full row-rank ###
  
  If the restriction matrix is not of full row-rank, this means one of the following:
  \itemize{
   \item There is at least one redundant restriction specified in the hypothesis. Then, either 
    \itemize{
      \item[a] Leave the redundant one out 
      \item[b] Use another (more time-consuming) way of obtaining the level probabilities 
      for the penalty term (goric function does this by default): Bootstrapping, as discussed above.
      }
   \item There is at least one range restriction (e.g., -2 < group1 < 2). 
         Such a restriction can be evaluated but there is a sensitivity (of a scaling 
         factor in the covariance matrix, like with a prior in a Bayes factor) which 
         currently cannot be checked for.
   \item There is at least one conflicting restriction (e.g., 2 < group1 < -2).
  }
    Such a restriction can evidently never hold and is thus impossible to evaluate. 
    To prevent this type of error delete the one that is incorrect.
}

\value{The function returns a dataframe with the log-likelihood,
penalty term, GORIC(A) values and the GORIC(A) weights. Furthermore, a dataframe
is returned with the relative GORIC(A) weights. 
}


\references{
Kuiper, R.M., Hoijtink, H., and Silvapulle, M.J. (2011). An Akaike-type
information criterion for model selection under inequality constraints.
\emph{Biometrika}, 98, 2, 495--501.
\url{https://doi.org/10.1093/biomet/asr002}

Vanbrabant, L. and Kuiper, R. (2020). Evaluating a theory-based hypothesis against 
its complement using an AIC-type information criterion with an application to 
facial burn injury. Psychological Methods.
Vanbrabant, L., Van Loey, N., and Kuiper, R. M. (2020). 
Evaluating a Theory-Based Hypothesis Against Its Complement Using an AIC-Type Information Criterion With an Application to Facial Burn Injury. 
\emph{Psychological Methods}, 25(2), 129-142. 
\url{https://doi.org/10.1037/met0000238}

Victoria Savalei and Peter M. Bentler (2009) A Two-Stage Approach to
Missing Data: Theory and Application to Auxiliary Variables, Structural Equation 
Modeling: A Multidisciplinary Journal, 16:3, 477-497, DOI: 10.1080/10705510903008238
}

\author{Leonard Vanbrabant and Rebecca M. Kuiper}


\examples{
## By following these examples, you can appropriately specify hypotheses based on 
## your research questions and analytical framework.

# The hypotheses (i.e., constraints) have to be in a list. It is recommended to name
# each hypothesis in the list. Otherwise the hypotheses are named accordingly 'H1', 'H2', \ldots.
# Another option is to use the \code{llist()} function from the \pkg{Hmisc} package, where.

# text-based syntax (the labels x1, x2, and x2 are the names of coef(model) or names(vector))
h1 <- '(x1, x2, x3) > 0'
h2 <- '(x1, x3) > 0; x2 = 0'
h3 <- 'x1 > 0; x2 < 0; x3 = 0'
hypotheses = list(hypo1 = h1, hypo2 = h2, hypo3 = h3)

# define constraints matrix directly (note that the constraints have to be specified pairwise).
# the element names (i.e., constraints, rhs, neq) must be used. 
h1 <- list(constraints = c(0,1,0)) 
h2 <- list(constraints = rbind(c(0,1,0), c(0,0,1)), rhs = c(0.5, 1), neq = 0) 
hypotheses = list(H1 = h1, H2 = h2)

# mixed syntax:  
hypotheses = list(Ha = h1, Hb = 'x1 = x2 > x3')

# lavaan object syntax:
# the recommended option for objects of class lavaan is to use labels (here a, b and c) 
# to define our hypothesis.
model.lav <- "y ~ 1 + a*x1 + b*x2 + c*x3 + x4"
# fit lavaan model, for example
# library(lavaan)
# fit.lav <- sem(model, data = DATA)
# define hypothesis syntax
hypotheses = list(h1 = 'a > b > c')
# if needed absolute values can be used.
hypotheses = list(h1 = 'abs(a) > abs(b) > abs(c)')

library(MASS)
## lm
## unrestricted linear model for ages (in months) at which an 
## infant starts to walk alone.

# prepare data
DATA <- subset(ZelazoKolb1972, Group != "Control")
  
# fit unrestrikted linear model
fit1.lm <- lm(Age ~ Group, data = DATA)

# some artificial restrictions
H1 <- "GroupPassive > 0; GroupPassive < GroupNo"
H2 <- "GroupPassive > 0; GroupPassive > GroupNo"
H3 <- "GroupPassive = 0; GroupPassive < GroupNo"


# object is of class lm
goric(fit1.lm, hypotheses = list(H1 = H1, H2 = H2, H3 = H3))

# same result, but using the parameter estimates and covariance matrix as input
# Note, that in case of a numeric input only the gorica(c) can be computed. 
goric(coef(fit1.lm), VCOV = vcov(fit1.lm), hypotheses = list(H1 = H1, H2 = H2, H3 = H3))


# hypothesis H1 versus the complement (i.e., not H1)
goric(fit1.lm, hypotheses = list(H1 = H1), comparison = "complement")


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
#biased VCOV
VCOV <- vcov(fit.lm)
#unbiased VCOV
N_min_k <- fit.lm$df.residual
N <- N_min_k + fit.lm$rank
VCOV_unbiased <- vcov(fit.lm) * N_min_k / N

## constraint syntax (character)
h1 <- "x1 > 0"
h2 <- "x1 > 0; x2 > 0"
# use fitted unconstrained linear model
goric(fit.lm, hypotheses = list(h1 = h1, h2 = h2), type = "gorica")
# use unconstrained estimates
# unbiased VCOV - then same result as with fit object above
goric(est, VCOV = VCOV_unbiased, hypotheses = list(h1 = h1, h2 = h2), type = "gorica")
# biased VCOV
goric(est, VCOV = VCOV, hypotheses = list(h1 = h1, h2 = h2), type = "gorica")

## constraint syntax (matrix notation)
h1 <- list(constraints = c(0,1,0))
h2 <- list(constraints = rbind(c(0,1,0), c(0,0,1)), rhs = c(0.5, 1), neq = 0)
goric(fit.lm, hypotheses = list(h1 = h1, h2 = h2), type = "gorica")
goric(est, VCOV = VCOV, hypotheses = list(h1 = h1, h2 = h2), type = "gorica")

}
