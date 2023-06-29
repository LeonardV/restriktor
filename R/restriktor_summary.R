## to do
# gebruik con_goric() functie om LL, Penalty and Goric te berekenen.

summary.restriktor <- function(object, bootCIs = TRUE, bty = "perc", 
                               level = 0.95, 
                               goric = "goric", ...) {
  z <- object
  
  if (!inherits(z, "restriktor")) {
    stop("object of class ", sQuote(class(z)), " is not supported.")
  }
  
  goric <- tolower(goric)
  stopifnot(goric %in% c("goric", "goricc", "gorica", "goricac", "none"))
  
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
  meq  <- z$neq
  p    <- z$model.org$rank
  rdf  <- z$df.residual
  b.restr <- z$b.restr
  r <- weighted.residuals(z)
  
  ans <- z[c("call", if (!is.null(z$weights)) "weights")]
  ans$model.org <- z$model.org
  
  se.type <- z$se
  ans$se.type <- se.type
    attr(ans$se.type, "bootCIs") <- bootCIs    
    attr(ans$se.type, "level")   <- level    
    attr(ans$se.type, "bty")     <- bty

  ans$residuals <- r
  if (is.null(z$bootout) && se.type != "none") {
    if (se.type == "standard") {
      V <- attr(z$information, "inverted")
      se <- sqrt(diag(V))
    } else {
      V <- sandwich(z, 
                    bread. = bread(z), 
                    meat.  = meatHC(z, type = se.type))
      se <- sqrt(diag(V))
    }
    ans$V <- V
    tval <- ifelse(se != 0, b.restr/se, 0L)
    ans$coefficients <- cbind(b.restr, se, tval, 2 * pt(abs(tval), 
                                                    rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(b.restr),
                                       c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    # add new defined parameters 
    if (any(z$parTable$op == ":=")) {
      b.def <- z$CON$def.function(b.restr)
      JAC <- lav_func_jacobian_complex(func = z$CON$def.function, x = b.restr)
      def.cov <- JAC %*% V %*% t(JAC)
      diag.def.cov <- diag(def.cov)
      diag.def.cov[ diag.def.cov < 0 ] <- as.numeric(NA)
      se.def <- sqrt(diag.def.cov)
      tval.def <- ifelse(se.def != 0, b.def / se.def, 0)
      ans$coefficients <- rbind(ans$coefficients, cbind(b.def, se.def, tval.def, 
                                                        2 * pt(abs(tval.def), 
                                                               rdf, lower.tail = FALSE)))
    }
  } else if (bootCIs && (ans$se.type %in% c("boot.model.based", "boot.standard"))) {
      cis <- matrix(0, length(z$b.restr), 2)
      colnames(cis) <- c("lower", "upper")
      for (i in 1:length(z$b.restr)) {
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
      ans$coefficients <- cbind(b.restr, se, cis)
      colnames(ans$coefficients) <- c("Estimate", "Std. Error", "Lower", "Upper")
      # bootstrapped standard errors for newly defined parameters
      if (any(z$parTable$op == ":=")) {
        b.def <- z$CON$def.function(b.restr)
        JAC <- lav_func_jacobian_complex(func = z$CON$def.function, x = b.restr)
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
  } else if (is.null(z$bootout) && se.type == "none" && !any(z$parTable$op == ":=")) {
    ans$coefficients <- cbind(b.restr)
    if (!(inherits(z,"conMLM"))) {
      colnames(ans$coefficients) <- "Estimate"
    }
  } else if (is.null(z$bootout) && se.type == "none" && any(z$parTable$op == ":=")) {
      b.def <- z$CON$def.function(b.restr)
      ans$coefficients <- cbind(c(b.restr, b.def))
      colnames(ans$coefficients) <- "Estimate"
  } else {
      stop("restriktor ERROR: you may have found a bug, please contact me at info@restriktor.org")
    }
  
  ny <- ncol(coef(object$model.org))
  if (!is.null(ny) && ny > 1L) {
    ynames <- colnames(ans$coefficients)
    if (is.null(ynames)) {
      lhs <- object$model.org$terms[[2L]]
      if (mode(lhs) == "call" && lhs[[1L]] == "cbind") {
        ynames <- as.character(lhs)[-1L]
      } else {
        ynames <- paste0("Y", seq_len(ny))
      }
    }
    ind <- ynames == ""
    if (any(ind)) 
      ynames[ind] <- paste0("Y", seq_len(ny))[ind]
      colnames(ans$coefficients) <- ynames
  }
  
  ans$rdf <- rdf
  if (inherits(z, "conRLM")) {
    ans$wgt <- z$wgt
    ans$scale <- z$scale
    ans$stddev <- z$stddev
    ans$iter <- z$iter
  } else {
    ans$s2 <- z$s2
  }
  
  if (inherits(z, c("conLM", "conRLM"))) {
    ans$R2.org <- z$R2.org
    if (attr(z$model.org$terms, "intercept") != p) {
      ans$R2.reduced <- z$R2.reduced
    } else {
      ans$R2.reduced <- as.numeric(NA)
    }
  }
  
  if (inherits(z, "conGLM")) {
    ans$family <- z$family
    ans$dispersion <- z$dispersion  
    ans$deviance.null <- z$deviance.null
    ans$deviance <- z$deviance
    ans$df.residual.null <- z$df.residual.null
  }
    
  wt.bar <- z$wt.bar
  ## compute goric
  if (goric != "none" && !(attr(wt.bar, "method") == "none")) {
    ## REF: Kuiper, R.M.; Hoijtink, H.J.A.; Silvapulle, M. J. (2012) 
    ## Journal of statistical planning and inference, volume 142, pp. 2454 - 2463
    
    # compute penalty term based on simulated level probabilities (wt.bar)
    # The value 1 is the penalty for estimating the variance/dispersion parameter.
    if (goric %in% c("goric", "gorica")) {
      PT <- penalty_goric(Amat        = Amat, 
                          meq         = meq, 
                          LP          = wt.bar, 
                          correction  = FALSE, 
                          sample.nobs = NULL)
      if (goric == "gorica") {
        PT <- PT - 1 
      }
    } else if (goric %in% c("goricc", "goricac")) {
      PT <- penalty_goric(Amat        = Amat, 
                          meq         = meq, 
                          LP          = wt.bar, 
                          correction  = TRUE, 
                          sample.nobs = length(r))
      if (goric == "goricac") {
        PT <- PT - 1 
      }
    } else {
      stop("Restriktor ERROR: ", sQuote(goric), ": unknown goric method.")  
    }
    
    # compute log-likelihood value
    if (goric %in% c("goric", "goricc")) {
      ll   <- z$loglik
    } else if (goric %in% c("gorica", "goricac")) {
      # unconstrained vcov
      VCOV <- z$Sigma
      ll <- dmvnorm(c(z$b.unrestr - z$b.restr), sigma = VCOV, log = TRUE)  
    }
    
    if (inherits(z, c("conLM", "conMLM"))) {
      ans$goric <- -2*(ll - PT) 
    } else if (inherits(z, "conGLM")) {
      if (!(z$model.org$family$family %in% c("gaussian", "Gamma", "inverse.gaussian"))) {
        PT <- PT - 1
      }
      ans$goric <- -2*ll / 1 + 2*PT
    }
    attr(ans$goric, "type")    <- goric
    attr(ans$goric, "penalty") <- PT
    attr(ans$goric, "loglik")  <- ll
  }
  
  ans$messages <- z$messages
  
  if (inherits(z, "conRLM")) {
    class(ans) <- c("summary.restriktor", "summary.conRLM")
  } else if (inherits(z, "conGLM")) {
    class(ans) <- c("summary.restriktor", "summary.conGLM")
  } else if (inherits(z, "conLM")) {
    class(ans) <- c("summary.restriktor", "summary.conLM")
  } else if (inherits(z, "conMLM")) {
    class(ans) <- c("summary.restriktor", "summary.conMLM")
  }
    
  ans
}
