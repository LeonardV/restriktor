n <- 10
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 1 + x1 + x2 + rnorm(n)
# fit unconstrained linear model
fit.lm <- lm(y ~ x1 + x2)


## constraint syntax (character)
h1 <- "x1 > 0"
h2 <- "x1 > 0; x2 > 0"
# use fitted unconstrained linear model
out <- goric(fit.lm, h1, h2, type = "gorica")
## constraint syntax (character)
h1 <- "x1 > 0"
h2 <- "x1 > 0, x2 > 0"
# use fitted unconstrained linear model
out1 <- goric(fit.lm, h1, h2, type = "gorica")
## constraint syntax (character)
h1 <- "x1 > 0"
h2 <- "x1 > 0 & x2 > 0"
# use fitted unconstrained linear model
out2 <- goric(fit.lm, h1, h2, type = "gorica")

## constraint syntax (character)
h1 <- "x1 > 0"
h2 <- "x1 > 0 and x2 > 0"
# use fitted unconstrained linear model
out3 <- goric(fit.lm, h1, h2, type = "gorica")

#those must be identical
out$constraints
out1$constraints
out2$constraints
out3$constraints

#those must be identical
out$rhs
out1$rhs
out2$rhs
out3$rhs

#those must be identical
out$neq
out1$neq
out2$neq
out3$neq