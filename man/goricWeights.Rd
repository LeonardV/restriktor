\name{goricWeights}
\alias{goricWeights}

\title{Generalized Order-Restrikted Information Criterion Weights}
\description{The \code{goricWeights} function computes GORIC 
weights, which are comparable to the Akaike weights. 
}
\usage{\method{goricWeights}{conLM and conRLM}(x, ...)
}

\arguments{
  \item{x}{a list with objects of class \code{summary.conLM} or
  \code{summary.conRLM}.}
  
  \item{\ldots}{no additional arguments for now.}
}

\details{The GORIC values themselves are not interpretable and 
the interest lie in their differences. The GORIC can be seen as a 
likelihood for hypothesis \emph{m}, the GORIC weights represent 
the weights of evidence of hypothesis \emph{m} give the data and 
the set of M hypotheses. For example, if hypothesis H1 has a 
GORIC weight of w_1 = 0.54, and hypothesis H2 has a GORIC weight 
of w_2 = 0.32. The relative evidence for hypothesis H1 is 
0.54 / 0.32 = 1.67 more likely than H2. Thus, although H1 is 
preferred, there is no compelling evidence, since H2 receives 
quite some supports as well. 
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

\seealso{ 
  \code{\link{goric}} 
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

# some artificial restriktions
fit1.con <- restriktor(fit1.lm, constraints = "Group2 > 0; Group2 < Group4")
fit2.con <- restriktor(fit1.lm, constraints = "Group2 > 0; Group2 > Group4")
fit3.con <- restriktor(fit1.lm, constraints = "Group2 == 0; Group2 < Group4")
fit4.con <- restriktor(fit1.lm, constraints = "Group2 == 0; Group2 == Group4")

sum1.con <- summary(fit1.con)
sum2.con <- summary(fit2.con)
sum3.con <- summary(fit3.con)
sum4.con <- summary(fit4.con)

goricWeights(list(sum1.con, sum2.con, sum3.con, sum4.con))
}
