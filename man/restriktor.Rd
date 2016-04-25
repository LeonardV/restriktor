\name{restriktor}
\alias{restriktor}

\title{Fitting linear regression models with equality and/or inequality constraints}
\description{}

\usage{
restriktor(model, constraints, se = "default", rhs = NULL, neq = NULL, 
           control = NULL, debug = FALSE, ...)
\method{restriktor}{restriktor}(model, constraints, se = "default", rhs = NULL, 
                    neq = NULL, control = NULL, debug = FALSE, ...)
}
\arguments{
  \item{model}{}
  \item{constraints}{}
  \item{se}{}
  \item{rhs}{vector of right hand side elements; \eqn{Constr \; \theta \geq
      rhs}; number should equal the number of rows of the constr matrix}
  \item{neq}{number of equality constraints; a numeric value treating
    the first nec constr rows as equality constraints, or a logical
    vector with \code{TRUE} for equality- and \code{FALSE} for
    inequality constraints.}
  \item{control}{a list of control arguments; see details.}
  \item{debug}{}
  }
\details{
The contraints in the hypothesis of interest are defined by 
\eqn{Constr}, \eqn{rhs}, and \eqn{nec}.
The first \eqn{nec} constraints are the equality contraints:
\eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]};
and the remaing ones are the inequality contraints:
\eqn{Constr[nec+1:c_m, 1:tk] \theta \geq rhs[nec+1:c_m]}.

}
\value{an object of class conLM}
\references{}
\author{Leonard Vanbrabant and Yves Rosseel}


\seealso{\code{\link{solve.QP}}, \code{\link{conTest}}}

\examples{
# unrestricted linear model for ages (in months) at which an infant starts to 
# walk alone.
DATA <- ZelazoKolb1972
idx <- which(DATA$Group == 3)
DATA <- DATA[-idx, ]
DATA$Group <- factor(DATA$Group)

fit.lm <- lm(Age ~ Group, data = DATA)
summary(fit.lm)

# restrikted linear model with restriktions that the walking exercises would 
# not have a negative effect of increasing the mean age at which a child starts 
# to walk. 

fit.con <- restriktor(fit.lm, constraints = "Group2 > 0; Group2 < Group4")
summary(fit.con)

########################
## Artificial example ##
########################
n <- 10
m <- c(1,2,1,5)
nm <- length(m)
dat <- data.frame(grp=as.factor(rep(1:nm, each=n)),
                  y=rnorm(n*nm, rep(m, each=n), 1))

# unrestricted linear model
fit.lm <- lm(y ~ grp-1, data = dat)

# increasing means
A1 <- rbind(c(-1,1,0,0),
            c(0,-1,1,0),
            c(0,0,-1,1))
             
fm1 <- restriktor(fit.lm, constraints = A1,
                  rhs = rep(0L, nrow(A1)), neq = 0)

# larger than average of the neighboring first 2 parameters
cm4 <- rbind(c(-0.5,-0.5,1,0),
             c(0,-0.5,-0.5,1))
fm4 <- restriktor(y ~ grp-1, data=dat, constr=cm4,
            rhs=rep(0,nrow(cm4)), nec=0)

# equality constraints (all parameters equal)
fm5 <- restriktor(fit.lm, constr = A1,
                 rhs = rep(0, nrow(A1)), neq = nrow(A1))

# alternatively
fm5 <- restriktor(y ~ grp-1, data=dat, constr=cm2,
            rhs=rep(0,nrow(cm2)), nec=c(TRUE,TRUE,TRUE))

# constraining the 1st and the 4th parameter
# to their true values, and the 2nd and 3rd between them
cm6 <- rbind(c( 1,0,0,0),
             c(-1,1,0,0),
             c(0,-1,0,1),
             c(-1,0,1,0),
             c(0,0,-1,1),
             c(0,0, 0,1))
fm6 <- restriktor(y ~ grp-1, data=dat, constr=cm6,
            rhs=c(1,rep(0,4),5), nec=c(TRUE,rep(FALSE,4),TRUE))


# constraint definition
cmat <- cbind(diag(3), 0) + cbind(0, -diag(3))
constr <- kronecker(diag(3), cmat)

# no effect model
fm0 <- retriktor(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
            constr=constr, rhs=rep(0, nrow(constr)), nec=nrow(constr))
fm0

# order constrained model (increasing serum levels with increasing doses)
fm1 <- retriktor(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
            constr=constr, rhs=rep(0, nrow(constr)), nec=0)
summary(fm1)

# unconstrained model
fmunc <- retriktor(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
              constr=matrix(0, nrow=1, ncol=12), rhs=0, nec=0)
fmunc
}

\keyword{models}
