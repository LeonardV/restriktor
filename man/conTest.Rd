\name{conTest}
\alias{conTest}
\alias{conTestF.lm}
\alias{conTestLRT.lm}
\alias{conTestScore.lm}

\title{Function for testing equality and/or inequality restrikted hypotheses for
linear models}
\description{\code{conTest} tests linear (in)equality constrained hypotheses for
linear models by Wald-, likelihood ratio-, and score-tests.}

\usage{
conTest(object, type = "A", ...)
\method{conTest}{lm and rlm}(object, type = "A", ...)
}
\arguments{
  \item{object}{an object of class \code{conLM} or \code{conRLM}.} 
  \item{type}{}
  \item{test}{}
  \item{boot}{}
  \item{neq.alt}{}
  \item{control}{}
  \item{...}{}
  
  object, type = "A", neq.alt = 0, boot = "no", B = 9999, 
                            p.distr = "N", df = 7, parallel = "no", ncpus = 1L,
                            cl = NULL, seed = 1234, verbose = FALSE,
                            control = NULL, ...
}

\details{}
\value{
  An object of class conTest or conRLM, for which a print and a 
  summary method is available. More specifically, it is a list 
  with the following items:
    
  \item{}{}
}

\references{
  
}

\author{Leonard Vanbrabant and Yves Rosseel}


\seealso{ 
  \code{\link{solve.QP}}, 
  \code{\link{conTest}},
  \code{\link{vcovHC}}
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
  
  # the variable names can be used to impose restriktions on
  # the corresponding regression parameters.
  coef(fit1.lm)
  
  # restrikted linear model with restriktions that the walking 
  # exercises would not have a negative effect of increasing the 
  # mean age at which a child starts to walk. 
  fit1.con <- restriktor(fit1.lm, constraints = "Group2 > 0; Group2 < Group4")
  summary(fit1.con)
  
  \dontrun{ 
    # Or in matrix notation.
    myConstraints_mat1 <- rbind(c(0, 1, 0),
                                c(0,-1, 1))
    myRhs1 <- rep(0L, nrow(myConstraints_mat1)) 
    myNeq1 <- 0
    
    fit1.con <- restriktor(fit1.lm, constraints = myConstraints_mat1,
                           rhs = myRhs1, neq = myNeq1)
    summary(fit1.con)
  }            
  
  #########################
  ## Artificial examples ##
  #########################
  # generate data
  n <- 10
  means <- c(1,2,1,3)
  nm <- length(means)
  group <- as.factor(rep(1:nm, each = n))
  y <- cbind(c(MASS:::mvrnorm(n, mu = means, Sigma = diag(nm)), nrow(n)))
  DATA2 <- data.frame(y, group)
  
  # fit unrestrikted robust linear model
  fit2.rlm <- MASS:::rlm(y ~ group, data = DATA2, method = "MM")
  coef(fit2.rlm)
  
  # increasing means
  myConstraints2 <- ' group2 > 0
  group2 < group3
  group3 < group4
  '
  
  # fit restrikted robust linear model and compute 
  # Huber-White (robust) standard errors.
  fit2.con <- restriktor(fit2.rlm, constraints = myConstraints2, se = "HC0")
  summary(fit2.con)
  
  \dontrun{ 
    # increasing means in matrix notation.
    myConstraints_mat2 <- rbind(c(0, 1, 0, 0),
                                c(0,-1, 1, 0),
                                c(0, 0,-1, 1))
    myRhs2 <- rep(0L, nrow(myConstraints_mat2)) 
    myNeq2 <- 0
    
    fit2.con <- restriktor(fit2.rlm, constraints = myConstraints_mat2,
                           rhs = myRhs2, neq = myNeq2)
    summary(fit2.con)
  }            
  
  # equality restriktions only.
  myConstraints3 <- ' group2 == 0
  group2 == group3
  group3 == group4
  '
  fit3.con <- restriktor(fit2.lm, constr = myConstraints3)
  summary(fit3.con)
  
  
  # combination of equality and inequality restriktions.
  myConstraints4 <- ' group2 == 0
  group3 < group4
  '
  
  # fit restikted model and compute model-based bootstrapped 
  # standard errors. We only generate 9 bootstrap samples in this 
  # example; in practice you may wish to use a much higher number.
  fit4.con <- restriktor(fit2.lm, constr = myConstraints4, 
                         se = "boot.model.based", B = 9)
  summary(fit4.con)
  
  # restriktor can also be used to define effects using the := operator 
  # and impose restriktions on them. For example, is the 
  # average effect (AVE) larger than zero?
  # generate data
  n <- 30
  b0 <- 10; b1 = 0.5; b2 = 1; b3 = 1.5
  X <- c(rep(c(0), n/2), rep(c(1), n/2))
  set.seed(90) 
  Z <- rnorm(n, 16, 5)
  y <- b0 + b1*X + b2*Z + b3*X*Z + rnorm(n, 0, sd = 10) 
  DATA3 = data.frame(cbind(y, X, Z))
  
  # fit linear model with interaction
  fit5.lm <- lm(y ~ X*Z, data = DATA3)
  
  fit5.con <- restriktor(fit5.lm, constraints = "AVE := X + 16.86137*X.Z; 
                                               AVE > 0")
  summary(fit5.con)
}