print.conTest <- function(x, digits = max(3, getOption("digits") - 4), brief = FALSE, ...) {

  if (!("conTest" %in% class(x))) {
    stop("x must be of class \"conTest\"")
  }
  if (nrow(x$Amat) == x$meq) {
    cat("\nRestriktor: restrikted hypothesis test\n")
  } else {
    cat("\nRestriktor: restrikted hypothesis test type", x$type, "\n")
    if (x$type != "C") {
      if (x$boot != "no") {
        cat("( Number of successful bootstrap draws:", attr(x$pvalue, "B"),")\n")
      }
    }
  }
  vnames <- names(x$b.unconstr)
  Amat <- x$Amat
  colnames(Amat) <- vnames
  out.rest <- cbind(Amat, c(rep("   ==", x$meq), rep("   >=", nrow(x$Amat) -
                                                     x$meq)), x$bvec)
  
  rownames(out.rest) <- paste(1:nrow(out.rest), ":", sep = "")
  
  colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs")
  out.rest <- cbind(rep(" ", nrow(out.rest)), out.rest)
  out.rest[x$iact, 1] <- "A"
  if (nrow(x$Amat) == x$meq) {
    out.rest[1:nrow(Amat), 1] <- "A"
  }  
  out.rest <- as.data.frame(out.rest)
  names(out.rest)[1] <- ""

  out.test <- c(sprintf("%.4f", x$Ts), 
                                   if (x$pvalue < 1e-04) { "<0.0001" } 
                                   else { sprintf("%.4f", x$pvalue) })
  names(out.test) <- c(" Test statistic", "p-value")

  if (nrow(x$Amat) > x$meq) {
    if (x$type == "global") {
      cat("\nGlobal model test:\n\n")
      print(out.test, quote = FALSE, scientific = FALSE)
      cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
      cat("\nrestrikted estimate under H0:\n")
      print.default(format(x$b.eqconstr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\nrestrikted estimate under union of H0 and HA:\n")
      print.default(format(x$b.constr, digits = digits),
                    print.gap = 2, quote = FALSE)
    
    } else if (x$type == "A") {
      cat("\n H0: all restriktions active (=)",
          "\n HA: at least one restriktion strictly true (>)","\n\n")
      print(out.test, quote = FALSE, scientific = FALSE)        
      cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
      cat("\nrestrikted estimate under H0:\n")
      print.default(format(x$b.eqconstr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\nrestrikted estimate under union of H0 and HA:\n")
      print.default(format(x$b.constr, digits = digits),
                    print.gap = 2, quote = FALSE)
    } else if (x$type == "B" && x$meq.alt == 0L) {
      cat("\n H0: all restriktions true (>=)",
          "\n HA: at least one restriktion violated (<)", "\n\n")
      print(out.test, quote = FALSE)
      cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
      cat("\nrestrikted estimate under H0:\n")
      print.default(format(x$b.constr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\nUnrestrikted estimate:\n")
      print.default(format(x$b.unconstr, digits = digits),
                    print.gap = 2, quote = FALSE)
    } else if (x$type == "B" && x$meq.alt > 0L) {
        cat("\n H0: all restriktions true (>= or =)",
            "\n HA: at least one restriktion violated (<), some =-restriktions maintained",
            "\n\n")
        print(out.test, quote = FALSE)
        cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestrikted estimate under H0:\n")
        print.default(format(x$b.constr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nrestrikted estimate under HA:\n")
        print.default(format(x$b.constr.alt, digits = digits),
                      print.gap = 2, quote = FALSE)
      } else if (x$type == "C") {
        cat("\n H0: at least one restriktion not strictly true (<=)",
            "\n HA: all restriktions strictly true (>)",
            "\n\n")
        print(out.test, quote = FALSE)
        cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nUnrestrikted estimate:\n")
        print.default(format(x$b.unconstr, digits = digits),
                      print.gap = 2, quote = FALSE)
      }
  } else { #equality constraints only
    cat("\n HA: at least one restriktion violated (=)",
        "\n\n")
    print(out.test, quote = FALSE)
    if (!brief) {
      cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
      cat("\nrestrikted estimate under union of H0 and HA:\n")
      print.default(format(x$b.constr, digits = digits),
                    print.gap = 2, quote = FALSE)
    }
  }
}

