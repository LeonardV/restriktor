\name{bootstrapD}
\alias{bootstrapD}
\alias{bootstrapD}
\alias{print.conTestLavaan}
\title{Bootstrapping a Lavaan Model}
\description{Bootstrap the D statistic.}
\usage{
  bootstrapD(h0 = NULL, h1 = NULL, constraints, type = "A", 
             bootstrap.type = "bollen.stine", R = 1000L,  
             return.D = FALSE, double.bootstrap = "no", 
             double.bootstrap.R = 500L, double.bootstrap.alpha = 0.05, 
             verbose = FALSE, warn = -1L, 
             parallel = c("no", "multicore", "snow"), ncpus = 1L, cl = NULL, 
             seed = NULL)
             
  \method{print}{conTestLavaan}(x, digits = max(3, getOption("digits") - 2), \dots)
}
\arguments{ 
  \item{h0}{An object of class \code{\linkS4class{lavaan}}. The restricted
    model.}
  \item{h1}{An object of class \code{\linkS4class{lavaan}}. The unrestricted
    model.}
  \item{x}{an object of class \code{conTestLavaan}.}
  \item{constraints}{The imposed (in)equality constraints on the model.}  
  \item{type}{hypothesis test type "A", "B".}
  \item{bootstrap.type}{If \code{"parametric"}, the parametric bootstrap is used. 
    If \code{"bollen.stine"}, the semi-nonparametric Bollen-Stine bootstrap 
    is used. The default is set to \code{"bollen.stine"}.}
  \item{R}{Integer. The number of bootstrap draws.}
  \item{return.D}{Logical; if \code{TRUE}, the function returns bootstrapped         
    D-values.}
  \item{double.bootstrap}{If \code{"standard"} (default) the genuine double bootstrap is 
    used to compute an additional set of plug-in p-values for each bootstrap       
    sample. If \code{"no"}, no double bootstrap is used. If \code{"FDB"}, 
    the fast double bootstrap is used to compute second level LRT-values for 
    each bootstrap sample. Note that the \code{"FDB"} is experimental and should 
    not be used by inexperienced users.}    
  \item{double.bootstrap.R}{Integer; number of double bootstrap draws. The default 
    value is set to 249.}
  \item{double.bootstrap.alpha}{The significance level to compute the adjusted 
    alpha based on the plugin p-values. Only used if \code{double.bootstrap = "standard"}. 
    The default value is set to 0.05.}
  \item{verbose}{If \code{TRUE}, show information for each bootstrap draw.}
  \item{warn}{Sets the handling of warning messages. See \code{\link{options}}.}
  \item{parallel}{The type of parallel operation to be used (if any).  If
    missing, the default is \code{"no"}.  }
  \item{ncpus}{Integer: number of processes to be used in parallel operation:
      typically one would chose this to the number of available CPUs.  }
  \item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
    \code{parallel = "snow"}.  If not supplied, a cluster on the local machine is
    created for the duration of the \code{bootstrapLavaan} or \code{bootstrapLRT}
    call.}
  \item{digits}{the number of significant digits to use when printing.}  
  \item{\ldots}{no additional arguments for now.}
  \item{seed}{An integer to set the seed. Or NULL if no reproducible seeds are
    needed.}
}
\author{Leonard Vanbrabant}
\references{
  Bollen, K. and Stine, R. (1992) Bootstrapping Goodness of Fit Measures in
  Structural Equation Models. Sociological Methods and Research, 21,
  205--229.
  
  Silvapulle, M.J. and Sen, P.K. (2005). \emph{Constrained 
  Statistical Inference}. Wiley, New York

  Yuan, K.-H., Hayashi, K., and Yanagihara, H. (2007). A class of population
  covariance matrices in the bootstrap approach to covariance structure analysis.
  Multivariate Behavioral Research, 42, 261--281.
}

\examples{
 \dontrun{
    #########################
    ### real data example ###
    #########################
    # Multiple group path model for facial burns example.
    
    # model syntax with starting values.
    burns.model <- 'Selfesteem ~ Age + c(m1, f1)*TBSA + HADS +
                               start(-.10, -.20)*TBSA  
                 HADS ~ Age + c(m2, f2)*TBSA + RUM +
                        start(.10, .20)*TBSA '
    
    
    # constraints syntax
    burns.constraints <- 'f2 > 0  ; m1 < 0
                          m2 > 0  ; f1 < 0
                          f2 > m2 ; f1 < m1'
    
    # we only generate 2 bootstrap samples in this example; in practice
    # you may wish to use a much higher number. 
    # the double bootstrap was switched off; in practice you probably 
    # want to set it to "standard".
    example1 <- conTestD(model = burns.model, data = FacialBurns,
                         R = 2, constraints = burns.constraints,
                         double.bootstrap = "no", group = "Sex")
    
    example1
    
    ##########################
    ### artificial example ###
    ##########################
    # Simple ANOVA model with 3 groups (N = 20 per group)
    set.seed(1234)
    Y <- cbind(c(rnorm(20,0,1), rnorm(20,0.5,1), rnorm(20,1,1)))
    grp <- c(rep("1", 20), rep("2", 20), rep("3", 20))
    Data <- data.frame(Y, grp)
    
    #create model matrix
    fit.lm <- lm(Y ~ grp, data = Data)
    mfit <- fit.lm$model
    mm <- model.matrix(mfit)
    
    Y <- model.response(mfit)
    X <- data.frame(mm[,2:3])
    names(X) <- c("d1", "d2")
    Data.new <- data.frame(Y, X)
    
    # model
    model <- 'Y ~ 1 + a1*d1 + a2*d2'
    
    # fit without constraints
    fit <- lavaan::sem(model, data = Data.new)
    
    # constraints syntax: mu1 < mu2 < mu3
    constraints <- ' a1 > 0
                     a1 < a2 '
    
    # we only generate 10 bootstrap samples in this example; in practice
    # you may wish to use a much higher number, say > 1000. The double 
    # bootstrap is not necessary in case of an univariate ANOVA model.
    example2 <- conTestD(model = model, data = Data.new, 
                         start = lavaan::parTable(fit),
                         R = 10L, double.bootstrap = "no",
                         constraints = constraints)
    example2
  }
}
