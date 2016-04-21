summary.conLM <- function(object, bootCIs = TRUE, bty = "basic", level = 0.95,
                          bootWt = FALSE, R = 9999, 
                          digits = max(3, getOption("digits") - 2),
                          signif.stars = getOption("show.signif.stars"), ...) {
  
  if (!inherits(object, "conLM")) {
    stop("object of class ", sQuote(class(object)), " is not supported.")
  }
  
  z <- object
  
  # bty = "stud" needs bootstrap variances
  if (bootCIs & !(bty %in% c("norm", "basic", "perc", "bca"))) {
    stop("bty is invalid.")
  }
  if (bootCIs & (level < 0.5 | level > 1)) {
    stop("invalid confidence level")
  }
  
  Amat <- z$Amat
  meq <- z$meq
  p <- z$model.org$rank
  rdf <- z$df.residual
  r <- c(z$residuals)
  #f <- z$fitted
  est <- z$b.constr
  #weights <- z$weights
  r <- c(weighted.residuals(z))
  
  ans <- z[c("call", if (!is.null(z$weights)) "weights")]
  ans$model.org <- z$model.org
  
  se.type <- z$se
  ans$se.type <- se.type
    attr(ans$se.type, "bootCIs") <- bootCIs    
    attr(ans$se.type, "level") <- level    
    attr(ans$se.type, "bty") <- bty

  ans$residuals <- r
  if (is.null(z$bootout) && se.type != "none") {
    vcovHC <- sandwich(z, bread.=bread(z), meat.=meatHC(z, type = se.type))
    se <- sqrt(diag(vcovHC))
    tval <- ifelse(se != 0, est/se, 0L)
    ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
                                                    rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(est),
                                       c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    # add new defined parameters 
    if (any(z$parTable$op == ":=")) {
      b.def <- z$CON$def.function(est)
      JAC <- lav_func_jacobian_complex(func = z$CON$def.function, x = est)
      def.cov <- JAC %*% vcovHC %*% t(JAC)
      diag.def.cov <- diag(def.cov)
      diag.def.cov[ diag.def.cov < 0 ] <- as.numeric(NA)
      se.def <- sqrt(diag.def.cov)
      tval.def <- ifelse(se.def != 0, b.def/se.def, 0L)
      ans$coefficients <- rbind(ans$coefficients, cbind(b.def, se.def, tval.def, 
                                                        2 * pt(abs(tval.def), 
                                                               rdf, lower.tail = FALSE)))
    }  
  } else if (bootCIs && (ans$se.type %in% c("boot.model.based", "boot.standard"))) {
    cis <- matrix(0, length(z$b.constr), 2)
    colnames(cis) <- c("lower", "upper")
    for (i in 1:length(z$b.constr)) {
      if (!bty %in% c("norm", "perc")) {
        cis[i, ] <- boot.ci(z$bootout, conf = level,
                            type = bty, index = i)[[bty]][4:5]
      } else if (bty == "perc") {
        cis[i, ] <- boot.ci(z$bootout, conf = level,
                            type = bty, index = i)[["percent"]][4:5]
      } else if (bty == "norm") {
        cis[i, ] <- boot.ci(z$bootout, conf = level,
                            type = bty, index = i)[["normal"]][2:3]
      }  
    }
    se <- apply(z$bootout$t, 2, sd)
    ans$coefficients <- cbind(est, se, cis)
    colnames(ans$coefficients) <- c("Estimate", "Std. Error", "Lower", "Upper")
  } else {
    ans$coefficients <- cbind(est)
    colnames(ans$coefficients) <- "Estimate"
  }
  
  ans$rdf <- rdf
  ans$R2.org <- z$R2.org
  if (attr(z$model.org$terms, "intercept") != p) {
    ans$R2.reduced <- z$R2.reduced
  } else {
    ans$R2.reduced <- 0
  }
  # compute goric
  # REF: Kuiper, R.M.; Hoijtink, H.J.A.; Silvapulle, M. J. (2012) 
  # Journal of statistical planning and inference, volume 142, pp. 2454 - 2463
  ## TO DO: add small samples correction 
  if (inherits(object, "conRLM")) {
    s2ml.unc <- c(summary_rlm(z$model.org, ml = TRUE)$stddev^2)
  } else { 
    s2ml.unc <- c(z$s2ml.unc)
  }
  X <- model.matrix(z$model.org)[,,drop=FALSE]
  invW <- kronecker(solve(s2ml.unc), t(X) %*% X)
  W <- solve(invW)
  # compute penalty term
  if (bootWt) {
    wt <- mix.boot(object, R = R)
    wt <- wt[-1]
    PT <- 1 + sum( (1:ncol(W)) * wt)
  } else if (!bootWt & (meq < nrow(Amat))) {
    wt <- rev(con_wt(Amat %*% W %*% t(Amat), meq = meq))
    start.idx <- 1 + (ncol(Amat)-nrow(Amat) - 1)
    end.idx <- ncol(Amat) - meq
    PT <- 1 + sum(start.idx:end.idx * wt)      
  } else if (!bootWt & (meq == nrow(Amat))) {
    stop("set bootWt = TRUE.")  
    }
  ans$goric <- -2*(z$loglik - PT)
  attr(ans$goric, "weights") <- wt
  attr(ans$goric, "penalty") <- PT
  attr(ans$goric, "loglik")  <- z$loglik 
  
  if (class(object)[1] == "conLM") {
    class(ans) <- "summary.conLM"
  } else if (class(object)[1] == "conRLM") {
    class(ans) <- c("summary.conLM","summary.conRLM")
  } 
    
  ans
}
