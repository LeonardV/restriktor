\name{Hurricanes}
\alias{Hurricanes}
\docType{data}
\title{
The Hurricanes Dataset
}
\description{
The data provide information on the effect of El Nino (Cold, Neutral, 
Warm) on the number of hurricanes from 1950 to 1995. 
}

\usage{data(Hurricanes)}
\format{
  A data frame of 46 observations of 3 variables.
  \describe{
    \item{\code{Year}}{}
    \item{\code{Hurricanes}}{Number of Hurricanes}
    \item{\code{ElNino}}{1=Cold, 2=Neutral, 3=Warm}
  }
}
\references{The original source of these data is the National Hurricane Center
in Australia. The dataset was extracted from the table on page 5 in Silvapulle 
and Sen (2005).
}
\examples{
head(Hurricanes)
}
