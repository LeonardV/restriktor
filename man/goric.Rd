\name{goric}
\alias{goric}
\alias{print.goric}

\title{Generalized Order-Restrikted Information Criterion Weights}
\description{The \code{goric} function computes GORIC 
weights, which are comparable to the Akaike weights. 
}
\usage{goric(object, \dots, complement = FALSE, bound = NULL, 
             digits = max(3, getOption("digits") - 2), debug = FALSE)

\method{print}{goric}(x, digits = max(3, getOption("digits") - 2), \dots)

}

\arguments{
  \item{object}{a fitted restriktor object of class conLM, conMLM, conRLM and conGLM.}
  
  \item{x}{an object of class \code{goric}.}

  \item{\ldots}{further fitted restriktor objects of class conLM, conMLM, conRLM and conGLM.}
  
  \item{complement}{logical; if TRUE the complement is computed. The complement can only
  be computed for one restriktor object at a time.}
  
  \item{bound}{not used yet.}
  
  \item{digits}{the number of significant digits to use when printing.}
  
  \item{debug}{if TRUE, debugging information is printed out.}
}

\details{The GORIC values themselves are not interpretable and 
the interest lie in their differences. The GORIC weights reflect 
the support of each hypothesis in the set. To compare two hypotheses 
(and not one to the whole set), one can examine the ratio of the 
two corresponding GORIC weights. To avoid selecting a weakly supported 
hypothesis as the best one, the unconstrained hypothesis is usually 
included as safeguard.

In case of one order-constrained hypothesis, say H1, the complement 
Hc can be computed as competing hypothesis. The complement is defined 
as Hc = not H1.
}

\value{The function returns a dataframe with the log-likelihood,
penalty term, goric, and the GORIC weights values.
}

\references{
Kuiper, R.M., Hoijtink, H., and Silvapulle, M.J. (2011). An Akaike-type
information criterion for model selection under inequality constraints.
\emph{Biometrika}, \bold{98}, 2, 495--501.
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

# fit1.con verus the complement
goric(fit1.con, complement = TRUE)


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

# fit5.con verus the complement
goric(fit5.con, complement = TRUE)
}
