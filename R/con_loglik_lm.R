# compute the loglikelihood based on the (un)constrained coefficients.
con_loglik_lm <- function(X, y, b, detU = 1) {
  n <- dim(X)[1]
  Sigma  <- (t(y - X%*%matrix(b, ncol=ncol(y))) %*%
               (y - X%*%matrix(b, ncol=ncol(y)))) / n #ML
  loglik <- (-n/2)*log(2*pi) + (-1/2)*(n*log(det(Sigma)) +
                                             ncol(y)*log(detU)) - (1/2)*n

  OUT <- list(loglik = loglik, Sigma = Sigma)

  OUT
}
