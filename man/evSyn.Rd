\name{evSyn}
\alias{evSyn}
\alias{evSyn.default}
\alias{evSyn.est}
\alias{evSyn.ICweights}
\alias{evSyn.ICvalues}
\alias{evSyn.LL}
\alias{plot.evSyn}

\title{GORIC(A) Evidence synthesis}
\description{GORIC(A) evidence synthesis aggregates the evidence for theory-based 
hypotheses from multiple studies that may use diverse designs to investigate the 
same central theory.}

\usage{evSyn(object, \ldots)

\method{evSyn}{est}(object, VCOV = NULL, constraints = NULL,
      type = c("equal", "added"), 
      comparison = c("unconstrained", "complement", "none"))

\method{evSyn}{LL}(object, PT, type = c("added", "equal"),
      hypo_names = NULL)

\method{evSyn}{ICvalues}(object, hypo_names = NULL)

\method{evSyn}{ICweights}(object, PriorWeights = NULL, hypo_names = NULL)

\method{plot}{evSyn}(x, \dots)
}

\arguments{
  \item{object}{Currently, the following objects can be processed:
    \itemize{
      \item a list with (standardized) parameter estimates;   
      \item a list with log-likelihood values;
      \item a list with GORIC(A) weights;
      \item a list wiht GORIC(A) values;
    }
  }
  \item{VCOV}{a list with covariance matrices of the (standardized) parameter estimates 
  of interest.} 
  
  \item{PT}{a list with the penalty values.}
  
  \item{x}{an object of class \code{evSyn}}

  \item{type}{type of evidence-synthesis approach: Equal-evidence approach 
  (\code{type = "equal"}) or Added-evidence approach (\code{type = "added"}).
  In case of an equal-evidence approach, aggregating evidence from, say, 5 studies 
  with n = 100 observations is the same as obtaining evidence from 1 study (as if it was
  possible) with n = 500 observations (like meta-analysis does). In the added-evidence 
  approach, the aggregated evidence from, says, 5 studies is stronger than as if 
  the data were combined (as if that was possible).}

  \item{constraints}{see \code{\link{restriktor}}.}

  \item{comparison}{see \code{\link{goric}}.} 
  
  \item{PriorWeights}{vector that represents the prior belief for this model. By default,
  equal prior weights are used (i.e., 1/(#hypotheses)). Notably, in case the prior
  weights do not sum to 1, it will be rescaled such that it does; which implies that
  relative importance can be used and not per se weights.}
  
  \item{hypo_names}{character vector for labelling the hypotheses. By default the
  names are set to H1, H2, ...}

  \item{\ldots}{This depends on the class of the object.}
}


\value{The output comprises, among other things, the cumulative and final evidence 
for the theory-based hypotheses.}


\author{Leonard Vanbrabant and Rebecca Kuiper}

\examples{
# Number of primary studies for which the evidence for a central theory should be
# determined and synthesized.
S <- 4 

# beta values from the analyses
est_1 <- c(beta1 = 0.09)
est_2 <- c(beta1 = 0.14)
est_3 <- c(beta1 = 1.09)
est_4 <- c(beta1 = 1.781)
Param_studies <- list(est_1, est_2, est_3, est_4)

# standard error of the beta's (from the primary studies)
vcov_est_1 <- matrix(c(0.029^2), nrow = 1)
vcov_est_2 <- matrix(c(0.054^2), nrow = 1)
vcov_est_3 <- matrix(c(0.093^2), nrow = 1)
vcov_est_4 <- matrix(c(0.179^2), nrow = 1)
CovMx_studies <- list(vcov_est_1, vcov_est_2, vcov_est_3, vcov_est_4)

# Set of hypotheses for each study
H0   <- "beta1 == 0"
Hpos <- "beta1 > 0"
Hneg <- "beta1 < 0"
Hypo_studies <- list(H0, Hpos, Hneg)

# Since this covers the whole space / covers all theories, we do not need a safeguard
# hypothesis:
comparison <- "none"

out1 <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
              constraints = Hypo_studies,
              type = "added", 
              comparison = "none") 
\dontrun{
  plot(out1)
}

}


