library(testthat)
library(restriktor)

test_check("restriktor")

set.seed(3013073)
n = 100
betas <- c(1,2,3,4)
p <- length(betas)
mu <- rep(0L, p)
X <- MASS:::mvrnorm(n, mu, Sigma = diag(p), empirical = TRUE) 
y <- 1 + X%*%betas + rnorm(n)
fit.lm <- lm(y~X)
s2 <- summary(fit.lm)$sigma^2
I <- 1/s2 * crossprod(cbind(1,X))  
I <- round(I, 10)

test_that("augmented information matrix", {
  
  I.aug <- con_augmented_information(information = I, 
                                     X = cbind(1,X), 
                                     b.unconstr = coef(fit.lm), 
                                     b.constr = coef(fit.lm),  
                                     Amat = rbind(c(0,1,0,0,0),
                                                  c(0,0,-1,0,0)), 
                                     bvec = c(0), 
                                     meq = 0)
          
  expect_that( I.aug, is_a("matrix") )Ttt
  expect_that( all(I.aug[,3]), equals(0) )
  expect_that( all(diag(I.aug)), equals(0.009138493) )
})

