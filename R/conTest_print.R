print.conTest <- function(x, digits = max(3, getOption("digits") - 2), ...) {

  if (!("conTest" %in% class(x))) {
    stop("x must be of class \"conTest\"")
  }
  
  cat("\nRestriktor: restrikted hypothesis test\n")
  
  if (x$type == "Ax" && !is.null(x$typ3)) { x$type <- "global" }
  
  if (x$type != "C" && nrow(x$Amat) > x$meq) {
    if (x$boot != "no") {
      cat("( Number of successful bootstrap draws:", attr(x$pvalue, "R"),")\n")
    }
  }

  vnames <- names(x$b.unrestr)
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
      cat("\n\nGlobal test: H0: all parameters are restrikted to be equal", "\n", 
          "        vs. HA: at least one restriktion strictly true", "\n\n")
      print(out.test, quote = FALSE, scientific = FALSE)
      cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
      cat("\nrestrikted estimate under H0:\n")
      print.default(format(x$b.eqrestr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\nrestrikted estimate under union of H0 and HA:\n")
      print.default(format(x$b.restr, digits = digits),
                    print.gap = 2, quote = FALSE)
    } else if (x$type == "A") {
      cat("\n\nType A test: H0: all restriktions active (=)", "\n", 
          "        vs. HA: at least one inequality restriktion strictly true", "\n\n")
      print(out.test, quote = FALSE, scientific = FALSE)        
      cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
      cat("\nrestrikted estimate under H0:\n")
      print.default(format(x$b.eqrestr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\nrestrikted estimate under union of H0 and HA:\n")
      print.default(format(x$b.restr, digits = digits),
                    print.gap = 2, quote = FALSE)
    } else if (x$type == "B" && x$meq.alt == 0L) {
      cat("\n\nType B test: H0: all restriktions true", "\n", 
          "        vs. HA: at least one restriktion violated ", "\n\n")
      print(out.test, quote = FALSE)
      cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
      cat("\nrestrikted estimate under H0:\n")
      print.default(format(x$b.restr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\nUnrestrikted estimate:\n")
      print.default(format(x$b.unrestr, digits = digits),
                    print.gap = 2, quote = FALSE)
    } else if (x$type == "B" && x$meq.alt > 0L) {
        cat("\n H0: all restriktions true (>= or =)",
            "\n HA: at least one restriktion violated (<), some =-restriktions maintained",
            "\n\n")
        print(out.test, quote = FALSE)
        cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nrestrikted estimate under H0:\n")
        print.default(format(x$b.restr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nrestrikted estimate under HA:\n")
        print.default(format(x$b.restr.alt, digits = digits),
                      print.gap = 2, quote = FALSE)
      } else if (x$type == "C") {
        cat("\n\nType C test: H0: at least one restriktion false or active (=)", 
            "\n", "        vs. HA: all restriktions strictly true (>)", "\n\n")
        print(out.test, quote = FALSE)
        cat("\n\n(rows indicated with an \"A\" are active restriktions)\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nunrestrikted estimate:\n")
        print.default(format(x$b.unrestr, digits = digits),
                      print.gap = 2, quote = FALSE)
      }
  } else { #equality constraints only
    cat("\n\nClassical test: H0: all restriktions active (=)", 
        "\n", "           vs. HA: at least one equality restriktion violated", "\n\n")
    print(out.test, quote = FALSE)
    cat("\n\n(all rows are active restriktions under H0, H1 is unrestrikted!)\n")
    print(out.rest, quote = FALSE, scientific = FALSE)
    cat("\nrestrikted estimate under H0:\n")
    print.default(format(x$b.restr, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\nunrestrikted estimate:\n")
    print.default(format(x$b.unrestr, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
}

