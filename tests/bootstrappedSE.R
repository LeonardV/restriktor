# test bootstrap standard errors
# with and without intercept
# include weights
# lm and rlm

library(restriktor)
n = 50
X <- cbind(1, mvtnorm:::rmvnorm(n, mean=rep(0,3), sigma=diag(3)))
colnames(X) <- c("intercept","x1","x2","x3")
betas <- c(1,0,0,0) 
y <- data.frame(y=X%*%betas + rnorm(n))
DATA <- cbind(y,X)

# lm - boot.standard
fit.ic.lm  <- lm(y ~ 1 + x1 + x2 +x3, data = DATA)
fit.nic.lm <- lm(y ~ -1 + x1 + x2 +x3, data = DATA)
fit.ic.w.lm  <- lm(y ~ 1 + x1 + x2 +x3, data = DATA, weights = abs(rnorm(n)))
fit.nic.w.lm <- lm(y ~ -1 + x1 + x2 +x3, data = DATA, weights = abs(rnorm(n)))

fit.ic.con  <- restriktor(fit.ic.lm, constraints = "x2 > 0; x3 > 0", se = "boot.standard", B = 99)
fit.nic.con <- restriktor(fit.nic.lm, constraints = "x2 > 0; x3 > 0", se = "boot.standard", B = 99)
fit.ic.w.con  <- restriktor(fit.ic.w.lm, constraints = "x2 > 0; x3 > 0", se = "boot.standard", B = 99)
fit.nic.w.con <- restriktor(fit.nic.w.lm, constraints = "x2 > 0; x3 > 0", se = "boot.standard", B = 99)

summary(fit.ic.con, bty = "basic")
summary(fit.nic.con, bty = "basic")
summary(fit.ic.w.con, bty = "basic")
summary(fit.nic.w.con, bty = "basic")

# lm - boot.model.based
fit.ic.con  <- restriktor(fit.ic.lm, constraints = "x2 > 0; x3 > 0", se = "boot.model.based", B = 99)
fit.nic.con <- restriktor(fit.nic.lm, constraints = "x2 > 0; x3 > 0", se = "boot.model.based", B = 99)
fit.ic.w.con  <- restriktor(fit.ic.w.lm, constraints = "x2 > 0; x3 > 0", se = "boot.model.based", B = 99)
fit.nic.w.con <- restriktor(fit.nic.w.lm, constraints = "x2 > 0; x3 > 0", se = "boot.model.based", B = 99)

summary(fit.ic.con, bty = "basic")
summary(fit.nic.con, bty = "basic")
summary(fit.ic.w.con, bty = "basic")
summary(fit.nic.w.con, bty = "basic")

################### rlm ################

# rlm - boot.standardlibrary(MASS)
fit.ic.rlm  <- rlm(y ~ 1 + x1 + x2 + x3, data = DATA, method = "MM")
fit.nic.rlm <- rlm(y ~ -1 + x1 + x2 + x3, data = DATA, method = "MM")
fit.ic.w.rlm  <- rlm(y ~ 1 + x1 + x2 + x3, data = DATA, weights = abs(rnorm(n)), method = "MM")
fit.nic.w.rlm <- rlm(y ~ -1 + x1 + x2 + x3, data = DATA, weights = abs(rnorm(n)), method = "MM")

fit.ic.con  <- restriktor(fit.ic.rlm, constraints = "x2 > 0; x3 > 0", se = "boot.standard", B = 99)
fit.nic.con <- restriktor(fit.nic.rlm, constraints = "x2 > 0; x3 > 0", se = "boot.standard", B = 99)
fit.ic.w.con  <- restriktor(fit.ic.rlm, constraints = "x2 > 0; x3 > 0", se = "boot.standard", B = 99)
fit.nic.w.con <- restriktor(fit.nic.rlm, constraints = "x2 > 0; x3 > 0", se = "boot.standard", B = 99)

summary(fit.ic.con, bty = "basic")
summary(fit.nic.con, bty = "basic")
summary(fit.ic.w.con, bty = "basic")
summary(fit.nic.w.con, bty = "basic")

# rlm - boot.model.based
fit.ic.con  <- restriktor(fit.ic.rlm, constraints = "x2 > 0; x3 > 0", se = "boot.model.based", B = 99)
fit.nic.con <- restriktor(fit.nic.rlm, constraints = "x2 > 0; x3 > 0", se = "boot.model.based", B = 99)
fit.ic.w.con  <- restriktor(fit.ic.rlm, constraints = "x2 > 0; x3 > 0", se = "boot.model.based", B = 99)
fit.nic.w.con <- restriktor(fit.nic.rlm, constraints = "x2 > 0; x3 > 0", se = "boot.model.based", B = 99)

summary(fit.ic.con, bty = "basic")
summary(fit.nic.con, bty = "basic")
summary(fit.ic.w.con, bty = "basic")
summary(fit.nic.w.con, bty = "basic")

