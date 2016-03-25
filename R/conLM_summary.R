summary.conLM <- function(x, digits = max(3, getOption("digits") - 2),
                          bootCIs = TRUE, bty = "basic", level = 0.95, 
                          signif.stars = getOption("show.signif.stars")) {

  # bty = "stud" needs bootstrap variances
  if (bootCIs & !is.null(x$bootout) & !bty %in% c("norm", "basic", "perc", "bca")) {
    stop("bty is invalid.")
  }
  if (bootCIs & !is.null(x$bootout) & (level < 0.5 | level > 1)) {
    stop("invalid confidence level")
  }

  if (class(x)[1] == "conLM") {
    cat("\nRestriktor: constrained linear model:\n\n")
  } else if (class(x)[1] == "conRLM") {
    cat("\nRestriktor: constrained robust linear model:\n\n")
  }
  
  cat("Residuals:\n")
  resid <- x$residuals
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- quantile(resid)
   names(rq) <- nam
  print(rq, digits)
  cat("\n")

  se <- x$se
  if (length(x$b.constr) && is.null(x$bootout) && !(se == "no")) {
    cat("Coefficients:\n")
    vcovHC <- sandwich(x, bread.=bread(x), meat.=meatHC(x, type = se))
    SE <-
#    if (se == "const" | se == "default") {
#      sqrt(diag(x$information.inverted))
#    } else {
      sqrt(diag(vcovHC))  
#    }  
    
    ##########
    tval <- ifelse(SE != 0, x$b.constr/SE, 0L)
    coefficients <- cbind(x$b.constr, SE, tval, 2 * pt(abs(tval),
                                                    x$df.residual, lower.tail = FALSE))
    dimnames(coefficients) <- list(names(x$model.org$coefficients),
                                   c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    
    #cat("\nDefined new paramters:\n")
    
    ######################### defined parameters ###############################
    coefficients.def.ciq <- coefficients.def.ceq <- NULL
    if (!(is.null(x$partable)) && any(x$partable$op == ":=")) {
      #equality constraints defined parameters
      if (x$meq > 0L) {
        eq.idx <- which(x$partable$op == "==")
        def.idx <- which(x$partable$op == ":=")
        lhs.labels <- all.vars(parse(file = "", text = x$partable$lhs[eq.idx]))
        rhs.labels <- all.vars(parse(file = "", text = x$partable$rhs[eq.idx]))
        eq.labels <- unique(c(lhs.labels, rhs.labels))
        def.names <- as.character(x$partable$lhs[def.idx])
        d.idx <- which(eq.labels%in%def.names)
        
        H.ceq <- x$Amat[1:x$meq,,drop=FALSE]  
        bvec.ceq <- x$bvec[1:x$meq]  
        #"new > 0; x2 > 0; new2 := 2*x3; new2 == 1; new := 2*x4 + x1"
        b.def.ceq <- as.numeric(H.ceq[d.idx,,drop=FALSE] %*% coef(x)) #- bvec.ceq[d.idx]                          #is this correct?
          b.def.ceq[abs(b.def.ceq) < 1e-14] <- 0L
          n.idx <- which(eq.labels%in%def.names)
          names(b.def.ceq) <- eq.labels[n.idx]
        cov.diag <- diag(H.ceq[d.idx,,drop=FALSE] %*% vcovHC %*% t(H.ceq[d.idx,,drop=FALSE]))
          cov.diag[cov.diag < 0] <- 0L
        SE.def.ceq <- sqrt(cov.diag)            
        tval.def.ceq <- ifelse(SE.def.ceq != 0, b.def.ceq/SE.def.ceq, 0L)
        coefficients.def.ceq <- rbind(cbind(b.def.ceq, SE.def.ceq, tval.def.ceq, 2 * pt(abs(tval.def.ceq),
                                            x$df.residual, lower.tail = FALSE)))
      }  
      
      if (any(x$partable$op == "<" | x$partable$op == ">")) {
        #inequality constraints defined parameters
        ineq.idx <- which(x$partable$op == ">" | x$partable$op == "<")
        def.idx <- which(x$partable$op == ":=")
        lhs.labels <- all.vars(parse(file = "", text = x$partable$lhs[ineq.idx]))
        rhs.labels <- all.vars(parse(file = "", text = x$partable$rhs[ineq.idx]))
        ineq.labels <- unique(c(lhs.labels, rhs.labels))
        def.names <- as.character(x$partable$lhs[def.idx])
        d.idx <- which(ineq.labels%in%def.names)
        
        H.ciq <- x$Amat[-c(1:x$meq),,drop=FALSE]
        bvec.ciq <- x$bvec[-c(1:x$meq)]  
        b.def.ciq <- as.numeric(H.ciq[d.idx,,drop=FALSE] %*% coef(x)) #- bvec.ciq[d.idx] 
          b.def.ciq[abs(b.def.ciq) < 1e-14] <- 0L
          n.idx <- which(ineq.labels%in%def.names)
          names(b.def.ciq) <- ineq.labels[n.idx]
        cov.diag <- diag(H.ciq[d.idx,,drop=FALSE] %*% vcovHC %*% t(H.ciq[d.idx,,drop=FALSE]))
          cov.diag[cov.diag < 0] <- 0L
        SE.def.ciq <- sqrt(cov.diag)            
        tval.def.ciq <- ifelse(SE.def.ciq != 0, b.def.ciq/SE.def.ciq, 0L)
        
        coefficients.def.ciq <- rbind(cbind(b.def.ciq, SE.def.ciq, tval.def.ciq, 2 * pt(abs(tval.def.ciq),
                                            x$df.residual, lower.tail = FALSE)))  
      }
    }  
    ############################################################################
    coefficients <- rbind(coefficients, coefficients.def.ceq, coefficients.def.ciq)
      coefficients[,4][coefficients[,4] < 2e-16] <- 2e-16
    printCoefmat(coefficients, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA")
    
    cat("\n")
    if (se == "const") {
      cat("Homoskedastic standard errors\n")
    } else {
      cat("Heteroskedastic robust standard errors:", se ,"\n")
    }
    
    if (round(x$R2.org,3) != round(x$R2.reduced, 3)) {
      cat("Constrained model: R2 reduced from", round(x$R2.org,3), "to", round(x$R2.reduced, 3),"\n")  
    } else {
      cat("Constrained model: R2 remains", round(x$R2.org,3),"\n")
    }
    
#    if (class(x)[1] == "conRLM") {
#      cat("Convergence in ")
#    }
    
  } else if (length(x$b.constr) && !is.null(x$bootout)) {
    if (bootCIs) {
      cis <- matrix(0, length(x$b.constr), 2)
      colnames(cis) <- c("lower", "upper")
      for (i in 1:length(x$b.constr)) {
        if (!bty %in% c("norm", "perc")) {
          cis[i, ] <- boot.ci(x$bootout, conf = level,
                              type = bty, index = i)[[bty]][4:5]
        } else if (bty == "perc") {
          cis[i, ] <- boot.ci(x$bootout, conf = level,
                              type = bty, index = i)[["percent"]][4:5]
        } else if (bty == "norm") {
          cis[i, ] <- boot.ci(x$bootout, conf = level,
                              type = bty, index = i)[["normal"]][2:3]
        }  
      }
      cat("\nCoefficients from constrained model\nwith",
          100 * level, "pct bootstrap confidence intervals (",bty,"):", "\n")
    }

    se <- apply(x$bootout$t, 2, sd)
    est <- round(coef(x), 9)
    icc <- cbind(est, se, round(cis, 9))
    colnames(icc) <- c("Estimate", "Std. Error", "Lower", "Upper")
    print(icc, quote = FALSE, digits = digits)

    cat("\nConstrained model: R2 reduced from", round(x$R2.org,3), "to", round(x$R2.reduced, 3),"\n")
  } else if (se == "no") {
    cat("Coefficients:\n")
    print(coef(x), digits = digits, scientific = FALSE, print.gap = 2L,
          quote = FALSE)
    cat("\nConstrained model: R2 reduced from", round(x$R2.org,3), "to", round(x$R2.reduced, 3),"\n")
  } else {
    cat("No coefficients\n")
  }  
  cat("\n")
  invisible(x)
}
