# compute the (weighted) loglikelihood
con_loglik_lm <- function(object, ...) {
  res <- as.matrix(object$residuals)
  s2 <- crossprod(res) / nrow(res)
  n <- length(res)
  
  w <- object$weights
  if (is.null(w)) {
    w <- rep.int(1, n)
  }
  else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }
  
  # weights are not supported for mlm
  if (ncol(res) == 1L) {
    OUT <- 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
  } else if (ncol(res) > 1L) {
    # mlm
    OUT <- (-n/2)*log(2*pi) + (-1/2)*(nrow(res)*log(det(s2)) + ncol(res)*log(1)) - (1/2)*n
  }
    
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


