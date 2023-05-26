\docType{data}
\name{myLLs}
\alias{myLLs}
\title{An example of log likelihood (LL) values}
\format{
  An object of class \code{matrix} (inherits from \code{array}) with 4 rows and 3 columns.
}
\usage{
  data(myLLs)
}
\description{
  An example of LL: A matrix with log likelihood values of size 4 x 3 (in general: S x 'NrHypos+1', where 'NrHypos+1' stands for the number of theory-based hypotheses plus a safeguard hypothesis (the complement or unconstrained)).
}
\examples{
  data(myLLs)
  
}
\keyword{datasets}