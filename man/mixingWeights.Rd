\name{con_weights_boot}
\alias{con_weights_boot}

\title{function for computing the chi-bar-square weights based on 
Monte Carlo simulation.}
\description{The null-distribution of the test statistics under 
inequality constraints takes the form of mixtures of F-distributions. 
This function computes these mixing weights (a.k.a chi-bar-square weights
and level probabilities). It can be used directly and is called by 
the \code{conTest} function.
}
\usage{con_weights_boot(VCOV, Amat, meq, 
                 R = 99999L, parallel = c("no", "multicore", "snow"),
                 ncpus = 1L, cl = NULL, seed = NULL, verbose = FALSE, ...)
}

\arguments{
  \item{VCOV}{variance-covariance matrix of the data for which
  the weights are to be calculated.}
  
  \item{Amat}{constraints matrix \eqn{R} (or a vector in 
  case of one constraint) and defines the left-hand side of the 
  constraint \eqn{R\theta \ge rhs}, where each row represents one 
  constraint. The number of columns needs to correspond to the 
  number of parameters estimated (\eqn{\theta}). The rows should 
  be linear independent, otherwise the function gives an 
  error. For more information about constructing the matrix \eqn{R} 
  and \eqn{rhs} see \code{\link{restriktor}}.}
  
  \item{meq}{integer (default = 0) treating the number of 
  constraints rows as equality constraints instead of inequality 
  constraints. For example, if \code{meq = 2}, this means that the 
  first two rows of the constraints matrix \eqn{R} are treated as 
  equality constraints. }
  
  \item{R}{integer; number of bootstrap draws for \code{mix.bootstrap}. 
  The default value is set to 99999.}
  
  \item{parallel}{the type of parallel operation to be used 
  (if any). If missing, the default is set "no".}
  
  \item{ncpus}{integer: number of processes to be used in parallel 
  operation: typically one would chose this to the number of 
  available CPUs.}
  
  \item{cl}{an optional parallel or snow cluster for use if 
  parallel = "snow". If not supplied, a cluster on the local 
  machine is created for the duration of the conTest call.}
  
  \item{seed}{seed value.}
  
  \item{verbose}{logical; if TRUE, information is shown at each 
  bootstrap draw.}
  
  \item{\ldots}{no additional arguments for now.}
}

\value{The function returns a vector with the mixing weights}

\references{
  Silvapulle, M.J. and Sen, P.K. (2005, p.79). \emph{Constrained 
  Statistical Inference}. Wiley, New York.
}

\author{Leonard Vanbrabant and Yves Rosseel}


\examples{
W <- matrix(c(1,0.5,0.5,1),2,2)
Amat <- rbind(c(0,1))
meq <- 0L
# we only generate 99 bootstrap samples in this 
# example; in practice you may wish to use a much higher number.
wt.bar <- con_weights_boot(W, Amat, meq, R = 99)
wt.bar
}
