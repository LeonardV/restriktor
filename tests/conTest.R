library(MASS)
library(restriktor)
n = 100
type = "A"
meq.alt = 0
betas = c(1,0,0,0,0)
cont = 0#.10
eta = 12
set.seed(3013073)

X <- cbind(1, mvtnorm:::rmvnorm(n, mean=rep(0,4), sigma=diag(4)))
colnames(X) <- c("intercept","x1","x2","x3","x4")
y <- data.frame(y=X%*%betas + rnorm(n))
  
idx <- sample(1:nrow(y), n*cont, replace=FALSE)                             
X[idx,2] <- rnorm(length(idx), 5, 0.1)  
y[idx,1] <- rnorm(length(idx), eta, 0.1)  
DATA <- cbind(y,X)
  
fit_rlm   <- rlm(y~x1+x2+x3+x4, data = DATA, method="MM", maxit = 5000)
fit_restr <- restriktor(fit_rlm, constraints = "x2 > 0; x3 > 0; x4 > 0", se = "none")

# robust tests with inequalities
conTest(fit_restr, test = "F", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "Wald", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "Wald2", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "score", type = type, meq.alt = meq.alt)[[1]]$pvalue

# robust tests with equalities
fit_restr <- restriktor(fit_rlm, constraints = "x2 == 0; x3 == 0; x4 == 0", se = "none")
conTest(fit_restr, test = "F", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "Wald", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "Wald2", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "score", type = type, meq.alt = meq.alt)[[1]]$pvalue


# non-robust tests with inequalities
fit_lm <- lm(y ~ x1+x2+x3+x4, data = DATA)
fit_restr <- restriktor(fit_lm, constraints = "x2 > 0; x3 > 0; x4 > 0", se = "none")
conTest(fit_restr, test = "F", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "LRT", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "score", type = type, meq.alt = meq.alt)[[1]]$pvalue

# non-robust tests with equalities
fit_restr <- restriktor(fit_lm, constraints = "x2 == 0; x4 == 0", se = "none")
conTest(fit_restr, test = "F", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "Wald", type = type, meq.alt = meq.alt)[[1]]$pvalue
conTest(fit_restr, test = "score", type = type, meq.alt = meq.alt)[[1]]$pvalue


