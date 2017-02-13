summary.restriktor <- function(object, bootCIs = TRUE, bty = "perc", 
                               level = 0.95, GORIC = TRUE, ...) {
  z <- object
  
  if (!inherits(z, "restriktor")) {
    stop("object of class ", sQuote(class(z)), " is not supported.")
  }
  # bty = "stud" needs bootstrap variances
  if (bootCIs & !(bty %in% c("norm", "basic", "perc", "bca"))) {
    if (bty == "stud") {
      stop("Restriktor ERROR: studentized intervals not implemented.")
    } else {
      stop("bty is invalid.")
    }
  }
  if (bootCIs & (level < 0.5 | level > 1)) {
    stop("invalid confidence level")
  }
  
  Amat <- z$constraints
  meq <- z$neq
  p <- z$model_org$rank
  rdf <- z$df.residual
  b_restr <- z$b_restr
  r <- c(weighted.residuals(z))
  
  ans <- z[c("call", if (!is.null(z$weights)) "weights")]
  ans$model_org <- z$model_org
  
  se_type <- z$se
  ans$se_type <- se_type
    attr(ans$se_type, "bootCIs") <- bootCIs    
    attr(ans$se_type, "level") <- level    
    attr(ans$se_type, "bty") <- bty

  ans$residuals <- r
  if (is.null(z$bootout) && se_type != "none") {
    if (se_type == "standard") {
      V <- attr(z$information, "inverted")
      se <- sqrt(diag(V))
    } else {
      V <- sandwich(z, 
                    bread. = bread(z), 
                    meat.  = meatHC(z, type = se_type))
      se <- sqrt(diag(V))
    }
    ans$V <- V
    tval <- ifelse(se != 0, b_restr/se, 0L)
    ans$coefficients <- cbind(b_restr, se, tval, 2 * pt(abs(tval), 
                                                    rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(b_restr),
                                       c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    # add new defined parameters 
    if (any(z$parTable$op == ":=")) {
      b.def <- z$CON$def.function(b_restr)
      JAC <- lav_func_jacobian_complex(func = z$CON$def.function, x = b_restr)
      def.cov <- JAC %*% V %*% t(JAC)
      diag.def.cov <- diag(def.cov)
      diag.def.cov[ diag.def.cov < 0 ] <- as.numeric(NA)
      se.def <- sqrt(diag.def.cov)
      tval.def <- ifelse(se.def != 0, b.def / se.def, 0)
      ans$coefficients <- rbind(ans$coefficients, cbind(b.def, se.def, tval.def, 
                                                        2 * pt(abs(tval.def), 
                                                               rdf, lower.tail = FALSE)))
    }
  } else if (bootCIs && (ans$se_type %in% c("boot.model.based", "boot.standard"))) {
      cis <- matrix(0, length(z$b_restr), 2)
      colnames(cis) <- c("lower", "upper")
      for (i in 1:length(z$b_restr)) {
        if (!bty %in% c("norm", "perc")) { # basic and adjusted percentile
          cis[i, ] <- boot.ci(z$bootout, conf = level,
                              type = bty, index = i)[[bty]][4:5]
        } else if (bty == "perc") { # percentile method 
          #bty name "perc" is different from name in list "percent"
          cis[i, ] <- boot.ci(z$bootout, conf = level,
                              type = bty, index = i)[["percent"]][4:5] 
        } else if (bty == "norm") { # normal approximation
          cis[i, ] <- boot.ci(z$bootout, conf = level,
                              type = bty, index = i)[["normal"]][2:3]
        }  
      }
      se <- apply(z$bootout$t, 2, sd)
      ans$coefficients <- cbind(b_restr, se, cis)
      colnames(ans$coefficients) <- c("Estimate", "Std. Error", "Lower", "Upper")
      # bootstrapped standard errors for newly defined parameters
      if (any(z$parTable$op == ":=")) {
        b.def <- z$CON$def.function(b_restr)
        JAC <- lav_func_jacobian_complex(func = z$CON$def.function, x = b_restr)
        bootout.def <- matrix(apply(z$bootout$t, 1, function(x) JAC %*% x), nrow(JAC))
        se.def <- apply(bootout.def, 1, function(x) sd(x))
        cis.def <- matrix(0, length(b.def), 2)
        colnames(cis) <- c("lower", "upper")
        for (i in 1:length(b.def)) {
          if (!bty %in% c("norm", "perc")) { 
            cis.def[i, ] <- boot.ci(z$bootout, conf = level, type = bty, 
                                    t0 = b.def[i], t = bootout.def[i,])[[bty]][4:5]
          } else if (bty == "perc") { 
            cis.def[i, ] <- boot.ci(z$bootout, conf = level, type = bty, 
                                    t0 = b.def[i], t = bootout.def[i,])[["percent"]][4:5] 
          } else if (bty == "norm") { 
            cis.def[i, ] <- boot.ci(z$bootout, conf = level, type = bty, 
                                    t0 = b.def[i], t = bootout.def[i,])[["normal"]][2:3]
          }  
        } 
        ans$coefficients <- rbind(ans$coefficients, cbind(b.def, se.def, cis.def))
      }
  } else if (is.null(z$bootout) && se_type == "none" && !any(z$parTable$op == ":=")) {
    ans$coefficients <- cbind(b_restr)
    colnames(ans$coefficients) <- "Estimate"
  } else if (is.null(z$bootout) && se_type == "none" && any(z$parTable$op == ":=")) {
      b.def <- z$CON$def.function(b_restr)
      ans$coefficients <- cbind(c(b_restr, b.def))
      colnames(ans$coefficients) <- "Estimate"
  } else {
      stop("restriktor ERROR")
    }
  
  ans$s2_unrestr <- z$s2_unrestr
  ans$s2_restr   <- z$s2_restr
  ans$rdf <- rdf
  
  if (inherits(z, c("conLM", "conRLM"))) {
    ans$R2_org <- z$R2_org
    if (attr(z$model_org$terms, "intercept") != p) {
      ans$R2_reduced <- z$R2_reduced
    } else {
      ans$R2_reduced <- 0
    }
  }
  
  if (inherits(z, "conGLM")) {
    ans$family <- z$family
    ans$dispersion <- z$dispersion  
    ans$deviance_null <- z$deviance_null
    ans$deviance <- z$deviance
    ans$df.residual_null <- z$df.residual_null
  }
    
  wt <- z$wt
  ## compute goric
  if (GORIC && !(attr(wt, "method") == "none")) {
    ## REF: Kuiper, R.M.; Hoijtink, H.J.A.; Silvapulle, M. J. (2012) 
    ## Journal of statistical planning and inference, volume 142, pp. 2454 - 2463
    
    # compute penalty term based on simulated level probabilities (wt)
    # The value 1 is the penalty for estimating the variance/dispersion parameter.
    if (attr(wt, "method") == "boot") {
      if (!(ncol(Amat) + 1 == length(wt))) {
        PT <- 1 + sum(0 : ncol(Amat) * wt)
      } else {
        warning("restriktor WARNING: unable to compute penalty for GORIC.")
        PT <- as.numeric(NA)
      }
      # unconstrained case
    } else if (attr(wt, "method") == "pmvnorm" && all(c(Amat) == 0)) {
      PT <- p + 1
    } else if (attr(wt, "method") == "pmvnorm") {
      min_C <- ncol(Amat) - nrow(Amat)
      max_C <- ncol(Amat) - meq
      PT <- 1 + sum(min_C : max_C * wt) 
    } else {
      stop("restriktor ERROR: unable to compute penalty for GORIC.")  
    }
    
    if (inherits(z, "conLM")) {
      ans$goric <- -2*(z$loglik - PT)
    } else if (inherits(z, "conGLM")) {
      if (!(z$model_org$family$family %in% c("gaussian", "Gamma", "inverse.gaussian"))) {
        PT <- PT - 1
      }
      ans$goric <- -2*z$loglik / 1 + 2*PT
    }
    attr(ans$goric, "penalty") <- PT
    attr(ans$goric, "loglik")  <- z$loglik 
  }
  
  if (inherits(z, "conRLM")) {
    class(ans) <- c("summary.restriktor", "summary.conRLM")
  } else if (inherits(z, "conGLM")) {
    class(ans) <- c("summary.restriktor", "summary.conGLM")
  } else if (inherits(z, "conLM")) {
    class(ans) <- c("summary.restriktor", "summary.conLM")
  } 
    
  ans
}
