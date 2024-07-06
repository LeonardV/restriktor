\name{con_weights_boot}
\alias{con_weights_boot}

\title{function for computing the chi-bar-square weights based on 
Monte Carlo simulation.}
\description{The null-distribution of the test statistics under 
inequality constraints takes the form of mixtures of F-distributions. 
This function computes these mixing weights (a.k.a chi-bar-square weights
and level probabilities). 
}
\usage{con_weights_boot(VCOV, Amat, meq, R = 1e5L, 
                 chunk_size = 5000L, convergence_crit = 1e-03, 
                 seed = NULL, verbose = FALSE, ...)
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
  
  \item{R}{integer; the maximum number of bootstrap draws for 
  \code{mix_weights_bootstrap_limit}. 
  The default value is set to 1e5. See details for more information.}
  
  \item{chunk_size}{integer; the chi-bar-square weights are computed for samples of
  size \code{chunk_size = 5000L}. This process is repeated iteratively until the 
  weights converges (see \code{convergenge_crit}) or the maximum is reached, i.e.,
  mix_weights_bootstrap_limit.}
  
  \item{convergence_crit}{the convergence criterion for the iterative process. 
  The default is 1e-03. See details for more information.}

  \item{seed}{seed value.}

  \item{verbose}{logical; if TRUE, information is shown at each 
  bootstrap draw.}
  
  \item{\ldots}{additional parameters for the \code{rtmvnorm} function. }
}

\details{
## Iterative Weight Updating and Convergence Checking ##
The function adds in each run chunks of 5000 samples (default) to compute 
the chi-bar-square weights. After each iteration, the function checks if the 
weights have converged. This is determined by the convergence_crit parameter.

Convergence is assessed by comparing the absolute difference between the 
current and previous iteration's weights against the convergence_crit. If the 
change in weights is smaller than the convergence criterion, it indicates that 
the weights have stabilized, suggesting convergence.

If the weights have not converged and the \code{mix_weights_bootstrap_limit} has 
not been reached, the function proceeds with adding another set of 5000 samples 
and updates the weights accordingly.If the maximum number of iterations is 
reached without convergence, the function returns the (non-converged) weights. 
In this situation, it is advisible to increase the number of 
\code{mix_weights_bootstrap_limit}.
}


\value{If convergence is reached, the function returns a vector with the mixing 
weights with the following attributes:

\item{total_bootstrap_draws}{total number of bootstrap draws}
\item{converged}{have the chi-bar-square weights converged}
\item{convergence_crit}{convergence criterium}
\item{wt_bar_chunk}{matrix with the chi-bar-square weights for each iteration}
\item{chunk_size}{how many samples are added in each iteration}
\item{total_chunks}{what is the maximum number of chunks based on 
\code{mix_weights_bootstrap_limit} and \code{chunk_size}}
\item{chunk_iter}{number of iterations run}
\item{error.idx}{which bootstrap samples were not succesful}
\item{mix_weights_bootstrap_limit}{the maximum number of bootstrap draws}
}

\references{
  Silvapulle, M.J. and Sen, P.K. (2005, p.79). \emph{Constrained 
  Statistical Inference}. Wiley, New York.
}

\author{Leonard Vanbrabant and Yves Rosseel}


\examples{
\donttest{
W <- matrix(c(1,0.5,0.5,1),2,2)
Amat <- rbind(c(0,1))
meq <- 0L
wt.bar <- con_weights_boot(W, Amat, meq, R = 99999)
wt.bar

# in practice you want to use are more conservative convergence criterion
wt.bar2 <- con_weights_boot(W, Amat, meq, R = 99999, convergence_crit = 1e-02)
wt.bar2
}
}
