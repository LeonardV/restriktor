library(restriktor)
n <- 100
means <- c(1,2,1,3)
nm <- length(means)
group <- as.factor(rep(1:nm, each = n))
y <- rnorm(n * nm, rep(means, each = n))
DATA <- data.frame(Age = y, Group = group)

myConstraints <- ' Group1 < Group2
                   Group2 < Group3
                   Group3 < Group4 '
## lm
fit_lm <- lm(Age ~ -1 + Group, data = DATA)

fit1a_con <- restriktor(fit_lm, constraints = myConstraints,
                       Wt = "mvnorm", bootWt.R = 999,
                       se = "boot.model.based", B = 99)
summary(fit1a_con)
conTest(fit1a_con)


fit1b_con <- restriktor(fit_lm, constraints = myConstraints,
                       Wt = "boot", bootWt.R = 999,
                       se = "boot.model.based", B = 99)
summary(fit1b_con)
conTest(fit1b_con)


fit1c_con <- restriktor(fit_lm, constraints = myConstraints,
                       Wt = "none", bootWt.R = 999,
                       se = "boot.model.based", B = 99)
summary(fit1c_con)
conTest(fit1c_con)


fit2a_con <- restriktor(fit_lm, constraints = myConstraints,
                       Wt = "mvnorm", bootWt.R = 999,
                       se = "boot.standard", B = 99)

summary(fit2a_con)
conTest(fit2a_con)

fit2b_con <- restriktor(fit_lm, constraints = myConstraints,
                       Wt = "boot", bootWt.R = 999,
                       se = "boot.standard", B = 99)

summary(fit2b_con)
conTest(fit2b_con)

fit2c_con <- restriktor(fit_lm, constraints = myConstraints,
                       Wt = "none", bootWt.R = 999,
                       se = "boot.standard", B = 99)

summary(fit2c_con)
conTest(fit2c_con)


fit3a_con <- restriktor(fit_lm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)
summary(fit3a_con)
conTest(fit3a_con)


fit3b_con <- restriktor(fit_lm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)
summary(fit3b_con)
conTest(fit3b_con)


fit3c_con <- restriktor(fit_lm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)
summary(fit3c_con)
conTest(fit3c_con)

## rlm
library(MASS)
fit_rlm <- rlm(Age ~ -1 + Group, data = DATA, method = "MM", maxit = 5000)

fit1a_con <- restriktor(fit_rlm, constraints = myConstraints,
                        Wt = "mvnorm", bootWt.R = 999,
                        se = "boot.model.based", B = 99)
summary(fit1a_con)
conTest(fit1a_con)


fit1b_con <- restriktor(fit_rlm, constraints = myConstraints,
                        Wt = "boot", bootWt.R = 999,
                        se = "boot.model.based", B = 99)
summary(fit1b_con)
conTest(fit1b_con)


fit1c_con <- restriktor(fit_rlm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.model.based", B = 99)
summary(fit1c_con)
conTest(fit1c_con)


fit2a_con <- restriktor(fit_rlm, constraints = myConstraints,
                        Wt = "mvnorm", bootWt.R = 999,
                        se = "boot.standard", B = 99)

summary(fit2a_con)
conTest(fit2a_con)

fit2b_con <- restriktor(fit_rlm, constraints = myConstraints,
                        Wt = "boot", bootWt.R = 999,
                        se = "boot.standard", B = 99)

summary(fit2b_con)
conTest(fit2b_con)

fit2c_con <- restriktor(fit_rlm, constraints = myConstraints,
                       Wt = "none", bootWt.R = 999,
                       se = "boot.standard", B = 99)

summary(fit2c_con)
conTest(fit2c_con)


fit3a_con <- restriktor(fit_rlm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)
summary(fit3a_con)
conTest(fit3a_con)


fit3b_con <- restriktor(fit_rlm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)
summary(fit3b_con)
conTest(fit3b_con)


fit3c_con <- restriktor(fit_rlm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)
summary(fit3c_con)
conTest(fit3c_con)


## glm
library(MASS)
fit_glm <- glm(Age ~ -1 + Group, data = DATA, family = gaussian)

fit1a_con <- restriktor(fit_glm, constraints = myConstraints,
                        Wt = "mvnorm", bootWt.R = 999,
                        se = "boot.model.based", B = 99)
summary(fit1a_con)
conTest(fit1a_con)


fit1b_con <- restriktor(fit_glm, constraints = myConstraints,
                        Wt = "boot", bootWt.R = 999,
                        se = "boot.model.based", B = 99)
summary(fit1b_con)
conTest(fit1b_con)


fit1c_con <- restriktor(fit_glm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.model.based", B = 99)
summary(fit1c_con)
conTest(fit1c_con)


fit2a_con <- restriktor(fit_glm, constraints = myConstraints,
                        Wt = "mvnorm", bootWt.R = 999,
                        se = "boot.standard", B = 99)

summary(fit2a_con)
conTest(fit2a_con)

fit2b_con <- restriktor(fit_glm, constraints = myConstraints,
                        Wt = "boot", bootWt.R = 999,
                        se = "boot.standard", B = 99)

summary(fit2b_con)
conTest(fit2b_con)

fit2c_con <- restriktor(fit_glm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)

summary(fit2c_con)
conTest(fit2c_con)


fit3a_con <- restriktor(fit_glm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)
summary(fit3a_con)
conTest(fit3a_con)


fit3b_con <- restriktor(fit_glm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)
summary(fit3b_con)
conTest(fit3b_con)


fit3c_con <- restriktor(fit_glm, constraints = myConstraints,
                        Wt = "none", bootWt.R = 999,
                        se = "boot.standard", B = 99)
summary(fit3c_con)
conTest(fit3c_con)

