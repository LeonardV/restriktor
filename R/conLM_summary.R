summary.conLM <- function(object, bootCIs = TRUE, bty = "basic", level = 0.95,
                          GORIC = TRUE, digits = max(3, getOption("digits") - 2),
                          signif.stars = getOption("show.signif.stars"), ...) {
  
  z <- object
  
  if (!inherits(z, "conLM")) {
    stop("object of class ", sQuote(class(z)), " is not supported.")
  }
  # bty = "stud" needs bootstrap variances
  if (bootCIs & !(bty %in% c("norm", "basic", "perc", "bca"))) {
    if (bty == "stud") {
      stop("bootstrap variances needed for studentized intervals.")
    } else {
      stop("bty is invalid.")
    }
  }
  if (bootCIs & (level < 0.5 | level > 1)) {
    stop("invalid confidence level")
  }
  
  Amat <- z$constraints
  meq <- z$neq
  p <- z$model.org$rank
  rdf <- z$df.residual
  r <- c(z$residuals)
  b.restr <- z$b.restr
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
    if (se.type == "standard") {
      V <- attr(z$information, "inverted.information")
      se <- sqrt(diag(V))
    } else {
      V <- sandwich(z, bread. = bread(z), meat. = meatHC(z, type = se.type))
      se <- sqrt(diag(V))
    }
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
    } else {
    ans$coefficients <- cbind(b.restr)
    colnames(ans$coefficients) <- "Estimate"
    }
  
  ans$rdf <- rdf
  ans$R2.org <- z$R2.org
  if (attr(z$model.org$terms, "intercept") != p) {
    ans$R2.reduced <- z$R2.reduced
  } else {
    ans$R2.reduced <- 0
  }
  
  if (GORIC) {
    wt <- z$wt
    bootWt <- attr(wt, "bootWt")
    # compute goric
    # REF: Kuiper, R.M.; Hoijtink, H.J.A.; Silvapulle, M. J. (2012) 
    # Journal of statistical planning and inference, volume 142, pp. 2454 - 2463
    ## TO DO: add small sample corrections
    # compute penalty term
    if (bootWt) { # compute mixing weights based on simulation
      PT <- 1 + sum( (0 : ncol(Amat)) * wt)
    } else if (!bootWt & (meq < nrow(Amat))) { # compute mixing weights based on mvnorm
    #  wt <- rev(con_wt(Amat %*% W %*% t(Amat), meq = meq))
      start.idx <- 1 + (ncol(Amat) - nrow(Amat) - 1)
      end.idx <- ncol(Amat) - meq
      PT <- 1 + sum(start.idx:end.idx * wt)      
    } else if (!bootWt & (meq == nrow(Amat))) { # only equality constraints
      PT <- 1 + sum( (0 : ncol(Amat)) * wt)
    } else {
      stop("restriktor ERROR: unable to compute penalty for GORIC.")  
    }
    ans$goric <- -2*(z$loglik - PT)
    attr(ans$goric, "penalty") <- PT
    attr(ans$goric, "loglik")  <- z$loglik 
  }
  
  if (class(object)[1] == "conLM") {
    class(ans) <- "summary.conLM"
  } else if (class(object)[1] == "conRLM") {
    class(ans) <- c("summary.conLM","summary.conRLM")
  } 
    
  ans
}
