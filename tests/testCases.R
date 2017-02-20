## test constrained lm, rlm, glm
# mixing weights
# bootstrap standard errors
# models with and without intercept
# case weights
# constraints on factors and interactions

library(restriktor)
n <- 100
p <- 4
betas = c(.1,.2,.3,.4,.5)
set.seed(3013073)
X <- cbind(mvtnorm:::rmvnorm(n, mean=rep(0,p), sigma=diag(p)), rbinom(n,1,0.5))
colnames(X) <- c("x1","x2","x3","x4","f1")
z <- X %*% betas        
y <- z + rnorm(n)
DATA <- data.frame(y, X)

# intercept model
model1 <- y ~  1 + x1 + x2 + x3 + x4
# no intercept model
model2 <- y ~ -1 + x1 + x2 + x3 + x4
# intercept model with interaction
model3 <- y ~ 1 + x1*f1 + x2*f1 + x3*f1 + x4*f1
# intercept model with interaction
model4 <- y ~ -1 + x1*f1 + x2*f1 + x3*f1 + x4*f1

############################ lm #################################
linmod1 <- lm(model1, data = DATA)
linmod2 <- lm(model2, data = DATA)
linmod3 <- lm(model3, data = DATA)
linmod1wt <- lm(model1, data = DATA, weights = abs(rnorm(n)))
linmod2wt <- lm(model2, data = DATA, weights = abs(rnorm(n)))
linmod1fac <- lm(model3, data = DATA)
linmod2fac <- lm(model4, data = DATA)

# check normal functionality mixing weights and compare with bootstrapped mixing weights
restr1a <- restriktor(linmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none")
restr1b <- restriktor(linmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none", 
                      mix.weights = "boot", mix.bootstrap = 19999, seed = 123)
if (!all(abs(restr1a$wt.bar - restr1b$wt.bar) < .01)) {
  stop("mixing weights are not approx. equal")
}

restr1a <- restriktor(linmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none")
restr1b <- restriktor(linmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none", 
                      mix.weights = "boot", mix.bootstrap = 19999, seed = 123)
if (!all(abs(restr1a$wt.bar - restr1b$wt.bar) < .01)) {
  stop("mixing weights are not approx. equal")
}

restr1a <- restriktor(linmod3, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none")
restr1b <- restriktor(linmod3, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none", 
                      mix.weights = "boot", mix.bootstrap = 19999, seed = 123)
if (!all(abs(restr1a$wt.bar - restr1b$wt.bar) < .01)) {
  stop("mixing weights are not approx. equal")
}

restr1a <- restriktor(linmod1wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none")
restr1b <- restriktor(linmod1wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none", 
                      mix.weights = "boot", mix.bootstrap = 19999, seed = 123)
if (!all(abs(restr1a$wt.bar - restr1b$wt.bar) < .01)) {
  stop("mixing weights are not approx. equal")
}

restr1a <- restriktor(linmod2wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none")
restr1b <- restriktor(linmod2wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "none", 
                      mix.weights = "boot", mix.bootstrap = 19999, seed = 123)
if (!all(abs(restr1a$wt.bar - restr1b$wt.bar) < .01)) {
  stop("mixing weights are not approx. equal")
}


# check normal functionality robust standard errors
restr1 <- restriktor(linmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "HC0", B = 100, mix.weights = "pmvnorm")
summary(restr1)
out1 <- iht(restr1, test = "F")
round(out1$global$pvalue, 5) == .0437
round(out1$A$pvalue, 5) == .02607
round(out1$B$pvalue, 5) == .00414  
out1 <- iht(restr1, test = "LRT")
round(out1$global$pvalue, 5) == .07524
round(out1$A$pvalue, 5) == .0427
round(out1$B$pvalue, 5) == .00464
out1 <- iht(restr1, test = "score")
round(out1$global$pvalue, 5) == .08896
round(out1$A$pvalue, 5) == .05128
round(out1$B$pvalue, 5) == .01154

restr2 <- restriktor(linmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "HC0", B = 100, mix.weights = "pmvnorm")
summary(restr2)
out2 <- iht(restr2, test = "F")
round(out2$global$pvalue, 5) == .05882
round(out2$A$pvalue, 5) == .03988
round(out2$B$pvalue, 5) == .00402
out2 <- iht(restr2, test = "LRT")
round(out2$global$pvalue, 5) == .09875
round(out2$A$pvalue, 5) == .06254
round(out2$B$pvalue, 5) == .00495
out2 <- iht(restr2, test = "score")
round(out2$global$pvalue, 5) == .10913
round(out2$A$pvalue, 5) == .06982
round(out2$B$pvalue, 5) == .01107

# check normal functionality standard bootstrap
restr1 <- restriktor(linmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr1)
out1 <- iht(restr1, test = "F")
round(out1$global$pvalue, 5) == 4e-5
round(out1$A$pvalue, 5) == .00011
round(out1$B$pvalue, 5) == .49232
out1 <- iht(restr1, test = "LRT")
round(out1$global$pvalue, 5) == 8e-05
round(out1$A$pvalue, 5) == .00016
round(out1$B$pvalue, 5) == .47838
out1 <- iht(restr1, test = "score")
round(out1$global$pvalue, 5) == .2e-04
round(out1$A$pvalue, 5) == .00031
round(out1$B$pvalue, 5) == .49369

restr2 <- restriktor(linmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
out2 <- iht(restr2, test = "F")
round(out2$global$pvalue, 5) == .00039
round(out2$A$pvalue, 5) == .00133
round(out2$B$pvalue, 5) == .48783
out2 <- iht(restr2, test = "LRT")
round(out2$global$pvalue, 5) == .00059
round(out2$A$pvalue, 5) == .00157
round(out2$B$pvalue, 5) == .47723
out2 <- iht(restr2, test = "score")
round(out2$global$pvalue, 5) == .00101
round(out2$A$pvalue, 5) == .00219
round(out2$B$pvalue, 5) == .4892

restr1 <- restriktor(linmod1wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr1)
out1 <- iht(restr1, test = "F")
round(out1$global$pvalue, 5) == 0
round(out1$A$pvalue, 5) == 0
round(out1$B$pvalue, 5) == 0.03582
out1 <- iht(restr1, test = "LRT")
round(out1$global$pvalue, 5) == 0
round(out1$A$pvalue, 5) == 0
round(out1$B$pvalue, 5) == .03421
out1 <- iht(restr1, test = "score")
round(out1$global$pvalue, 5) == 2e-05
round(out1$A$pvalue, 5) == 1e-05
round(out1$B$pvalue, 5) == .09424

restr2 <- restriktor(linmod2wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
out2 <- iht(restr2, test = "F")
round(out2$global$pvalue, 5) == 0
round(out2$A$pvalue, 5) == 3e-05
round(out2$B$pvalue, 5) == .25947
out2 <- iht(restr2, test = "LRT")
round(out2$global$pvalue, 5) == 1e-05
round(out2$A$pvalue, 5) == 6e-05
round(out2$B$pvalue, 5) == .25084
out2 <- iht(restr2, test = "score")
round(out2$global$pvalue, 5) == .00221
round(out2$A$pvalue, 5) == .00248
round(out2$B$pvalue, 5) == .46014

restr1 <- restriktor(linmod1fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr1)
out1 <- iht(restr1, test = "F")
round(out1$global$pvalue, 5) == .00054
round(out1$A$pvalue, 5) == .00012
round(out1$B$pvalue, 5) == .25609
out1 <- iht(restr1, test = "LRT")
round(out1$global$pvalue, 5) == .00099
round(out1$A$pvalue, 5) == .00015
round(out1$B$pvalue, 5) == .22567
out1 <- iht(restr1, test = "score")
round(out1$global$pvalue, 5) == .00291
round(out1$A$pvalue, 5) == .00045
round(out1$B$pvalue, 5) == .26551

restr2 <- restriktor(linmod2fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
out2 <- iht(restr2, test = "F")
round(out2$global$pvalue, 5) == .00025
round(out2$A$pvalue, 5) == .00015
round(out2$B$pvalue, 5) == .23206
out2 <- iht(restr2, test = "LRT")
round(out2$global$pvalue, 5) == .00058
round(out2$A$pvalue, 5) == 2e-04
round(out2$B$pvalue, 5) == .20667
out2 <- iht(restr2, test = "score")
round(out2$global$pvalue, 5) == .00178
round(out2$A$pvalue, 5) == .00054
round(out2$B$pvalue, 5) == .24195


# check normal functionality model-based bootstrap
restr1 <- restriktor(linmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(linmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")
restr1 <- restriktor(linmod1wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(linmod2wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")
restr1 <- restriktor(linmod1fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(linmod2fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")


# check functionality computation p-value methods
out <- iht(restr1, test = "F") 
out.bootpar   <- iht(restr1, test = "F", boot = "parametric", R = 999, parallel = "multicore", ncpus = 2)
out.bootmodel <- iht(restr1, test = "F", boot = "model.based", R = 999, parallel = "multicore", ncpus = 2)

if (!(out$global$pvalue[1] - out.bootpar$global$pvalue[1] - out.bootmodel$global$pvalue[1] < 1e-03)) {
  stop("check calculation pvalue global test")
}
if (!(out$A$pvalue[1] - out.bootpar$A$pvalue[1] - out.bootmodel$A$pvalue[1] < 1e-03)) {
  stop("check calculation pvalue global test")
}
if (!(out$B$pvalue[1] - out.bootpar$B$pvalue[1] - out.bootmodel$B$pvalue[1] < 1e-03)) {
  stop("check calculation pvalue global test")
}


########################### rlm #################################
library(MASS)
rlinmod1 <- rlm(model1, data = DATA, method = "MM", maxit = 5000)
rlinmod2 <- rlm(model2, data = DATA, method = "MM", maxit = 5000)
rlinmod3 <- rlm(model3, data = DATA, method = "MM", maxit = 5000)
rlinmod1wt <- rlm(model1, data = DATA, weights = abs(rnorm(n)), method = "MM", maxit = 5000)
rlinmod2wt <- rlm(model2, data = DATA, weights = abs(rnorm(n)), method = "MM", maxit = 5000)
rlinmod1fac <- rlm(model3, data = DATA, method = "MM", maxit = 5000)
rlinmod2fac <- rlm(model4, data = DATA, method = "MM", maxit = 5000)


# check normal functionality robust standard errors
restr1 <- restriktor(rlinmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "HC0", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "Wald") 
iht(restr1, test = "Wald2")
iht(restr1, test = "score")
restr2 <- restriktor(rlinmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "HC0", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "Wald") 
iht(restr2, test = "Wald2")
iht(restr2, test = "score")
# check normal functionality standard bootstrap
restr1 <- restriktor(rlinmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "Wald") 
iht(restr1, test = "Wald2")
iht(restr1, test = "score")
restr2 <- restriktor(rlinmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "Wald") 
iht(restr2, test = "Wald2")
iht(restr2, test = "score")
restr1 <- restriktor(rlinmod1wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "Wald") 
iht(restr1, test = "Wald2")
iht(restr1, test = "score")
restr2 <- restriktor(rlinmod2wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "Wald") 
iht(restr2, test = "Wald2")
iht(restr2, test = "score")
restr1 <- restriktor(rlinmod1fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "Wald") 
iht(restr1, test = "Wald2")
iht(restr1, test = "score")
restr2 <- restriktor(rlinmod2fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "Wald") 
iht(restr2, test = "Wald2")
iht(restr2, test = "score")

# check normal functionality model-based bootstrap
restr1 <- restriktor(rlinmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "Wald") 
iht(restr1, test = "Wald2")
iht(restr1, test = "score")
restr2 <- restriktor(rlinmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "Wald") 
iht(restr2, test = "Wald2")
iht(restr2, test = "score")
restr1 <- restriktor(rlinmod1wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "Wald") 
iht(restr1, test = "Wald2")
iht(restr1, test = "score")
restr2 <- restriktor(rlinmod2wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "Wald") 
iht(restr2, test = "Wald2")
iht(restr2, test = "score")
restr1 <- restriktor(rlinmod1fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "Wald") 
iht(restr1, test = "Wald2")
iht(restr1, test = "score")
restr2 <- restriktor(rlinmod2fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "Wald") 
iht(restr2, test = "Wald2")
iht(restr2, test = "score")

########################### glm #################################
glinmod1 <- glm(model1, data = DATA)
glinmod2 <- glm(model2, data = DATA)
glinmod3 <- glm(model3, data = DATA)
glinmod1wt <- glm(model1, data = DATA, weights = abs(rnorm(n)))
glinmod2wt <- glm(model2, data = DATA, weights = abs(rnorm(n)))
glinmod1fac <- glm(model3, data = DATA)
glinmod2fac <- glm(model4, data = DATA)

# check normal functionality robust standard errors
restr1 <- restriktor(glinmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "HC0", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(glinmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "HC0", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")

# check normal functionality standard bootstrap
restr1 <- restriktor(glinmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(glinmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")
restr1 <- restriktor(glinmod1wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(glinmod2wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")
restr1 <- restriktor(glinmod1fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(glinmod2fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")

# check normal functionality model-based bootstrap
restr1 <- restriktor(glinmod1, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(glinmod2, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")
restr1 <- restriktor(glinmod1wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(glinmod2wt, constraints = "x2 > 0; x3 > 0; x4 == 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")
restr1 <- restriktor(glinmod1fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.model.based", B = 100, mix.weights = "pmvnorm")
summary(restr1)
iht(restr1, test = "F") 
iht(restr1, test = "LRT") 
iht(restr1, test = "score")
restr2 <- restriktor(glinmod2fac, constraints = "x2 > 0; x3 > 0; x4 == 0; f1.x3 < 0", se = "boot.standard", B = 100, mix.weights = "pmvnorm")
summary(restr2)
iht(restr2, test = "F") 
iht(restr2, test = "LRT") 
iht(restr2, test = "score")

