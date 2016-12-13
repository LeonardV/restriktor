\name{AngerManagement}
\alias{AngerManagement}
\docType{data}
\title{
  Reduction of aggression levels Dataset (4 treatment groups)
}
\description{
  The anger management dataset consists of reduction of aggression levels
  between week 1 (intake) and week 8 (end of training) from four 
  different treatment groups (No-exercises, Physical-exercises, 
  Behavioral-exercises, combination of physical and behavioral
  exercises). 
}

\usage{data(AngerManagement)}
\format{
  A data frame of 40 observations of 4 treatment variables and covariate age.
  \describe{
    \item{\code{Anger}}{reduction in aggression levels}
    \item{\code{Group}}{No, Physical, Behavioral, Both}
    \item{\code{Age}}{persons' age}
  }
}

\references{Hoijtink, H. Informative Hypotheses: Theory and Practice 
for Behavioral and Social Scientists Boca Raton, FL: Taylor & Francis, 2012.
}

\examples{
  head(AngerManagement)
}
