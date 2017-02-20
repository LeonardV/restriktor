# compute the (weighted) loglikelihood based on the regression coefficients.
con_loglik_lm <- function(object, ...) {
  res <- object$residuals
  N <- length(res)
  if (is.null(w <- object$weights)) {
    w <- rep.int(1, N)
  }
  else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }
  OUT <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
  
  OUT
}


con_loglik_glm <- function(x) {
  p <- x$rank
  if (x$family$family %in% c("gaussian", "Gamma", "inverse.gaussian")) {
    p <- p + 1
  }
  
  OUT <- p - x$aic / 2
  
  OUT
}


