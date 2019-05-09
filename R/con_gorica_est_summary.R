summary.gorica_est <- function(object, ...) {
  z <- object
  
  if (!inherits(z, "gorica_est")) {
    stop("object of class ", sQuote(class(z)), " is not supported.")
  }
  
  Amat <- z$constraints
  meq <- z$neq
  b.restr <- z$b.restr
  
  ans <- z$model.org
  
  wt.bar <- z$wt.bar
  ## compute goric
  if (!(attr(wt.bar, "method") == "none")) {
    ## REF: Kuiper, R.M.; Hoijtink, H.J.A.; Silvapulle, M. J. (2012) 
    ## Journal of statistical planning and inference, volume 142, pp. 2454 - 2463
    
    # compute penalty term based on simulated level probabilities (wt.bar)
    if (all(c(Amat) == 0)) {
      # unconstrained case
      PT <- length(b.restr)
    } else if (attr(wt.bar, "method") == "boot") { 
      PT <- sum(0 : ncol(Amat) * wt.bar)  
    } else if (attr(wt.bar, "method") == "pmvnorm") {
      min.C <- ncol(Amat) - nrow(Amat)
      max.C <- ncol(Amat) - meq
      PT <- sum(min.C : max.C * wt.bar) 
    } else {
      stop("restriktor ERROR: unable to compute penalty for GORIC.")  
    }
    
    ans$goric <- -2*(z$loglik - PT)
    attr(ans$goric, "penalty") <- PT
    attr(ans$goric, "loglik")  <- z$loglik 
  }
  
  class(ans) <- c("summary.gorica_est")
  
  
  ans
}
