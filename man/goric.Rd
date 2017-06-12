\name{goric}
\alias{goric}

\title{Generalized Order-Restrikted Information Criterion Weights}
\description{The \code{goric} function computes GORIC 
weights, which are comparable to the Akaike weights. 
}
\usage{goric(object, \dots, complement = FALSE, digits = max(3, getOption("digits") - 2))
}

\arguments{
  \item{object}{a fitted restriktor object of class conLM or conRLM.}
  
  \item{\ldots}{further fitted restriktor objects of class conLM or conRLM.}
  
  \item{complement}{logical; if TRUE the complement is computed. The complement can only
  be computed for one restriktor object at a time.}
  
  \item{digits}{the number of significant digits to use when printing.}
}

\details{The GORIC values themselves are not interpretable and 
the interest lie in their differences. The GORIC weights reflect 
the support of each hypothesis in the set. To compare two hypotheses 
(and not one to the whole set), one can examine the ratio of the 
two corresponding GORIC weights. To avoid drawing wrong conclusions
it is recommended to include the unconstrained hypothesis, see examples.

In case of one hypothesis, say H1, the complement Hc can be computed 
as competing hypothesis. The complement is defined as Hc = not H1.
}

\value{The function returns a dataframe with the log-likelihood,
penalty term, goric, and the GORIC weights values.
}

\references{
Kuiper, R.M., Hoijtink, H., and Silvapulle, M.J. (2011). An Akaike-type
information criterion for model selection under inequality constraints.
\emph{Biometrika}, \bold{98}, 2, 495--501.
}

\author{Leonard Vanbrabant and Yves Rosseel}


\examples{
## unrestricted linear model for ages (in months) at which an 
## infant starts to walk alone.

# prepare data
DATA <- subset(ZelazoKolb1972, Group != "Control")
  
# fit unrestrikted linear model
fit1.lm <- lm(Age ~ Group, data = DATA)

# some artificial restrictions
restr1 <- restriktor(fit1.lm, constraints = "GroupPassive > 0;  GroupPassive < GroupNo")
restr2 <- restriktor(fit1.lm, constraints = "GroupPassive > 0;  GroupPassive > GroupNo")
restr3 <- restriktor(fit1.lm, constraints = "GroupPassive == 0; GroupPassive < GroupNo")
restr4 <- restriktor(fit1.lm) # unrestricted model

goric(restr1, restr2, restr3, restr4)

## use the complement as competing hypothesis for the first hypothesis

goric(restr1, complement = TRUE)

}
