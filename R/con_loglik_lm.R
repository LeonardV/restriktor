# compute the (weighted) loglikelihood based on the regression coefficients.
con_loglik_lm <- function(X, y, b, w) {
  n <- dim(X)[1]
  if (is.null(w)) {
    w <- rep(1, n)
  }
  W <- diag(w)
  S  <- ( t(y - X %*% matrix(b, ncol = ncol(y)) ) %*% W %*% 
            (y - X %*% matrix(b, ncol = ncol(y))) ) 
  ll <- 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(det(S))))

  OUT <- list(loglik = ll, Sigma = S)

  OUT
}