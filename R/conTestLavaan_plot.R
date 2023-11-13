
plot.conTestLavaan <- function(x, ...,
                               type       = c("test", "ppv"),
                               main       = "main",
                               xlab       = "xlabel",
                               ylab       = "Frequency",
                               freq       = TRUE,
                               breaks     = 15,
                               cex.main   = 1,
                               cex.lab    = 1,
                               cex.axis   = 1,
                               col        = "grey",
                               border     = par("fg"),
                               vline      = TRUE,
                               vline.col  = c("red", "blue"),
                               lty        = c(1,2),
                               lwd        = 1,
                               legend     = TRUE,
                               bty        = "o",
                               cex.legend = 1,
                               loc.legend = "topright") {
  object <- x
  return.test <- object$return.test
  double.bootstrap <- object$double.bootstrap
  double.bootstrap.alpha <- object$double.bootstrap.alpha
  pvalue <- c(object$bootA[1], object$bootB[1])
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  
  par(mfrow = c(1, 2))
  if (length(type) == 2) {
    par(mfrow = c(2, 2))
  }
  
  if (return.test && (type == "test" || length(type) == 2)) {
    D.obs  <- c(attr(object$bootA, "D.original"),
                attr(object$bootB, "D.original"))
    D.A <- attr(object$bootA, "D")
    D.B <- attr(object$bootB, "D")
    if (length(D.A) - length(D.B) < 0L) {
      D <- as.data.frame(cbind(c(D.A, rep(as.numeric(NA), length(D.B) -
                                                length(D.A))), D.B))
    } else {
      D <- as.data.frame(cbind(D.A, c(D.B, rep(as.numeric(NA),
                                                     length(D.A) -
                                                       length(D.B)))))
    }
    names(D) <- c("D.A", " D.B")
    
    if (xlab == "xlabel") {
      xlab.D <- c("Bootstrapped D values")
    }
    if (main == "main") {
      main.D <- c("Distr. of D values - Type A",
                    "Distr. of D values - Type B")
    }
    
    for (i in 1:2) {
      plot <- hist(D[,i], plot = FALSE, breaks = breaks)
      plot(plot, ...,
           freq     = freq,
           main     = main.D[i],
           xlab     = xlab.D,
           ylab     = ylab,
           cex.axis = cex.axis,
           cex.main = cex.main,
           cex.lab  = cex.lab,
           col      = col,
           border   = border,
           axes     = FALSE,
           xaxt     = "n")
      
      axis(side = 1)
      axis(side = 2)
      box(lty = 1, col = "black")
      
      if (vline) {
        abline(v   = D.obs[i],
               col = vline.col[1],
               lty = lty[1],
               lwd = lwd)
      }
      if (legend) {
        ppvalue  <- sprintf("%.2f", pvalue[i])
        obs.D  <- sprintf("%.2f", D.obs[i])
        ppval    <- paste0("plug-in p value = ", ppvalue)
        obs.D  <- paste0("observed D = ", obs.D)
        legend.obj <- c(obs.D, ppval)
        if (!vline) {
          legend(loc.legend, legend.obj,
                 lty = c(0, 0),
                 lwd = lwd,
                 cex = cex.legend,
                 bty = bty)
        } else {
          legend(loc.legend, legend.obj,
                 lty = c(lty[1], 0),
                 col = vline.col[1],
                 lwd = lwd,
                 cex = cex.legend,
                 bty = bty)
        }
      }
    }
  }
  
  if (double.bootstrap == "standard" && (type == "ppv" || length(type) == 2)) {
    ppvalue.A <- attr(object$bootA, "plugin.pvalues")
    ppvalue.B <- attr(object$bootB, "plugin.pvalues")
    adj.a <- c(quantile(ppvalue.A, double.bootstrap.alpha),
               quantile(ppvalue.B, double.bootstrap.alpha))
    adj.ppv <- c(attr(object$bootA, "adj.pvalue"),
                 attr(object$bootB, "adj.pvalue"))
    if (length(ppvalue.A) - length(ppvalue.B) < 0L) {
      ppv <- as.data.frame(cbind(c(ppvalue.A, rep(NA, length(ppvalue.B) -
                                                    length(ppvalue.A))), ppvalue.B))
    } else {
      ppv <- as.data.frame(cbind(ppvalue.A, c(ppvalue.B, rep(NA, length(ppvalue.A) -
                                                               length(ppvalue.B)))))
    }
    names(ppv) <- c("ppA", "ppB")
    
    if (xlab == "xlabel") {
      xlab.ppv  <- c("Bootstrapped plug-in p-values")
    }
    if (main == "main") {
      main.ppv  <- c("Distr. of plug-in p-values - Type A",
                     "Distr. of plug-in p-values - Type B")
    }
    
    for (i in 1:2) {
      plot <- hist(ppv[,i], plot = FALSE, breaks=breaks)
      plot(plot, ...,
           freq     = freq,
           main     = main.ppv[i],
           xlab     = xlab.ppv,
           ylab     = ylab,
           cex.axis = cex.axis,
           cex.main = cex.main,
           cex.lab  = cex.lab,
           col      = col,
           border   = border,
           axes     = FALSE,
           xaxt     = "n")
      
      axis(side = 1, at = seq(0,1,0.1))
      axis(side = 2)
      box(lty = 1, col = "black")
      if (vline) {
        abline(v   = adj.a[i],
               col = vline.col[1],
               lty = lty[1],
               lwd = lwd)
        abline(v   = adj.ppv[i],
               col = vline.col[2],
               lty = lty[2],
               lwd = lwd)
      }
      if (legend) {
        adj.alpha  <- sprintf("%.2f", adj.a[i])
        adj.pval   <- sprintf("%.2f", adj.ppv[i])
        adja <- paste0("Adjusted alpha = ", adj.alpha)
        adjp <- paste0("Adjusted p-value = ", adj.pval)
        legend.obj <- c(adja, adjp)
        if (!vline) {
          legend(loc.legend, legend.obj,
                 lty = 0,
                 col = vline.col,
                 lwd = lwd,
                 cex = cex.legend,
                 bty = bty)
        } else {
          legend(loc.legend, legend.obj,
                 lty = lty,
                 col = vline.col,
                 lwd = lwd,
                 cex = cex.legend,
                 bty = bty)
        }
      }
    }
  }
}
