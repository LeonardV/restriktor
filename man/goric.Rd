\name{goric}
\alias{goric}

\title{Generalized Order-Restrikted Information Criterion Weights}
\description{The \code{goricWeights} function computes GORIC 
weights, which are comparable to the Akaike weights. 
}
\usage{goric(x, \dots)
}

\arguments{
  \item{x}{a list with objects of class \code{summary.conLM} or
  \code{summary.conRLM}.}
  
  \item{\ldots}{no additional arguments for now.}
}

\details{The GORIC values themselves are not interpretable and 
the interest lie in their differences. The GORIC weights reflect 
the support of each hypothesis in the set. To compare two hypotheses 
(and not one to the whole set), one should examine the ratio of the 
two corresponding GORIC weights. To avoid draweing wrong conclusions
it is recommended to include the unrestrikted model, see examples.
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
# unrestricted linear model for ages (in months) at which an 
# infant starts to walk alone.
DATA <- ZelazoKolb1972
idx <- which(DATA$Group == "Control")
DATA <- DATA[-idx, ]

# fit unrestrikted linear model
fit1.lm <- lm(Age ~ Group, data = DATA)

# some artificial restriktions
fit1.con <- restriktor(fit1.lm, constraints = "GroupPassive > 0;  GroupPassive < GroupNo")
fit2.con <- restriktor(fit1.lm, constraints = "GroupPassive > 0;  GroupPassive > GroupNo")
fit3.con <- restriktor(fit1.lm, constraints = "GroupPassive == 0; GroupPassive < GroupNo")
fitu.con <- restriktor(fit1.lm) # unrestrikted model

sum1.con <- summary(fit1.con)
sum2.con <- summary(fit2.con)
sum3.con <- summary(fit3.con)
sum4.con <- summary(fitu.con)

goric(list(sum1.con, sum2.con, sum3.con, sum4.con))
}
