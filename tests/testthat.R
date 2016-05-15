library(testthat)
library(restriktor)

test_check("restriktor")

DATA <- restriktor::ZelazoKolb1972
idx <- which(DATA$Group == 3)
DATA <- DATA[-idx, ]
DATA$Group <- factor(DATA$Group)
fit.lm <- lm(Age ~ Group, data = DATA)
X <- model.matrix(fit.lm)[,,drop=FALSE]
s2 <- summary(fit.lm)$sigma^2
I <- 1/s2 * crossprod(X)
Amat <- rbind(c(0,1,0), c(0,-1,1))
bvec <- rep(0, nrow(Amat))
meq = 0

test_that("augmented information matrix", {
  I.aug <- restriktor:::con_augmented_information(information = I, 
                                     X = X, 
                                     b.unrestr = coef(fit.lm), 
                                     b.restr = coef(fit.lm),  
                                     Amat = Amat, 
                                     bvec = bvec, 
                                     meq = meq)
          
  expect_that( I.aug, is_a("matrix") )
  expect_that( mean(I.aug), equals(0.1360317, tolerance = 0.000001) )
})


