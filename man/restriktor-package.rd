\name{restriktor-package}
\alias{restriktor-package}
\title{Package for equality and inequality restricted estimation and hypothesis testing}
\description{
Package \code{restriktor} implements estimation, testing and evaluating of linear equality and 
inequality restriktions about parameters and effects for univariate and multivariate 
normal models and generalized linear models.}

\details{
  \tabular{ll}{
  Package: \tab restriktor\cr
  Type: \tab Package\cr
  Version: \tab 0.2-15\cr
  Date: \tab 2019-05-24\cr
  License: \tab GPL (>=2)\cr
  LazyLoad: \tab yes\cr
  }
  
  Function \code{restriktor} estimates the parameters of an univariate
  and multivariate linear model (\code{lm}), robust estimation of the 
  linear model (\code{rlm}) or a generalized linear model (\code{glm}) 
  subject to linear equality and/or inequality restriktions. The real 
  work horses are the \code{conLM}, \code{conMLM}, the \code{conRLM}, 
  and the \code{conGLM} functions. A major advantage of \pkg{restriktor} 
  is that the constraints can be specified by a text-based description. 
  This means that users do not have to specify the complex constraint matrix 
  (comparable with a contrast matrix) themselves. 
  
  The function \code{restriktor} offers the possibility to compute 
  (model robust) standard errors under the restriktions. The 
  parameter estimates can also be bootstrapped, where bootstrapped 
  standard errors and confidence intervals are available via the 
  summary function. Moreover, the function computes the Generalized 
  Order-restricted Information Criterion (GORIC), which is a 
  modification of the AIC and a generalization of the ORIC.
  
  The function \code{iht} (alias \code{conTest}) conducts restricted 
  hypothesis tests. F, Wald/LRT and score test-statistics are available. 
  The null-distribution of these test-statistics takes the form of a 
  mixture of F-distributions. The mixing weights (a.k.a. chi-bar-square 
  weights or level probabilities) can be computed using the multivariate 
  normal distribution function with additional Monte Carlo steps or via 
  a simulation approach. Bootstrap methods are available to calculate the 
  mixing weights and to compute the p-value directly. Parameters estimates 
  under the null- and alternative-hypothesis are available from the 
  summary function. 
  
  The function \code{goric} (generalized order-restricted information
  criterion) computes GORIC values, weights and relative-weights or GORICA
  (generalized order-restricted information crittion approximation) values,
  weights and relative weights. The GORIC(A) values are comparable to the AIC 
  values. The function offers the possibility to evaluate an order-restricted 
  hypothesis against its complement, the unconstrained hypothesis or against
  a set of hypotheses. For now, only one order-restricted hypothesis can be 
  evaluated against its complement but work is in progress to evaluate a set 
  of order-restricted hypothesis against its complement. 
  
  The package makes use of various other R packages: \pkg{quadprog} 
  is used for restricted estimation, \pkg{boot} for bootstrapping, 
  \pkg{ic.infer} for computing the mixing weights based on the 
  multivariate normal distribution, \pkg{lavaan} for parsing the 
  constraint syntax. 
}

\value{
  The output of function \code{restriktor} belongs to S3 class 
  \code{conLM}, \code{conMLM}, \code{conRLM} or \code{conGLM}.  
  
  The output of function \code{conTest} belongs to S3 class \code{conTest}. 
  
  These classes offer print and summary methods. 
}

\section{Acknowledgements}{
  This package uses as an internal function the function 
  \code{nchoosek} from \pkg{ic.infer}, which is originally from 
  \pkg{vsn}, authored by Wolfgang Huber, available under LGPL. 
  
  The output style of the \code{iht} print function is strongly 
  inspired on the summary of the \code{ic.test} function from the 
  \pkg{ic.infer} package.
}


\examples{
## Data preparation
## Ages (in months) at which an infant starts to walk alone.
DATA <- ZelazoKolb1972
idx <- which(DATA$Group == "Control")
DATA <- DATA[-idx, ]

## unrestricted linear model 
fit.lm <- lm(Age ~ -1 + Group, data = DATA)
summary(fit.lm)

## restricted linear model with restriktions that the walking 
## exercises would not have a negative effect of increasing the 
## mean age at which a child starts to walk. 

myConstraints <- ' GroupActive  < GroupPassive; 
                   GroupPassive < GroupNo '
                   
fit.con <- restriktor(fit.lm, constraints = myConstraints)
summary(fit.con)

}

\references{
    Groemping, U. (2010). Inference With Linear Equality And Inequality
    Constraints Using R: The Package ic.infer. \emph{Journal of Statistical 
    Software}, Forthcoming. 

    Kuiper R.M., Hoijtink H., Silvapulle M.J. (2011). An Akaike-type Information
    Criterion for Model Selection Under Inequality Constraints. \emph{Biometrika}, 
    \bold{98}, 495--501.

    Kuiper R.M., Hoijtink H., Silvapulle M.J. (2012). Generalization of the 
    Order-Restricted Information Criterion for Multivariate Normal Linear Models. 
    \emph{Journal of Statistical Planning and Inference}, \bold{142}, 2454--2463. 
    doi:10.1016/j.jspi.2012.03.007.
    
    Robertson T, Wright F, Dykstra R (1988). \emph{Order-Restricted Inference}. 
    Wiley, New York.

    Schoenberg, R. (1997). Constrained Maximum Likelihood. \emph{Computational 
    Economics}, \bold{10}, 251--266.

    Shapiro, A. (1988). Towards a unified theory of inequality-constrained 
    testing in multivariate analysis. \emph{International Statistical Review} 
    \bold{56}, 49--62.
    
    Silvapulle, M. (1992a). Robust tests of inequality constraints and one-sided 
    hypotheses in the linear model. \emph{Biometrika}, \bold{79}, 621--630.

    Silvapulle, M. (1992b). Robust wald-type tests of one-sided hypotheses in 
    the linear model. \emph{Journal of the American Statistical Association}, 
    \bold{87}, 156--161.

    Silvapulle, M. (1996). Robust bounded influence tests against one-sided 
    hypotheses in general parametric models. \emph{Statistics & probability 
    letters}, \bold{31}, 45--50.
    
    Silvapulle, M.J. and Sen, P.K. (2005). \emph{Constrained Statistical Inference}. 
    Wiley, New York
    
    Vanbrabant, L. and Kuiper, R. (n.d.). Giving the complement a compliment: Evaluating 
    a theory-based hypothesis against its complement using the GORIC. 
}

\author{ Leonard Vanbrabant and Yves Rosseel - Ghent University}
\seealso{ 
See also \code{\link{restriktor}}, \code{\link{iht}}, 
          packages \pkg{boot}, \pkg{goric}, \pkg{ic.infer}, 
          \pkg{mvtnorm}, and \pkg{quadprog}.
}
