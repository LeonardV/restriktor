print.conTest <- function(x, digits = max(3, getOption("digits") - 4), brief = FALSE, ...) {

  if (!("conTest" %in% class(x))) {
    stop("x must be of class \"conTest\"")
  }
  if (nrow(x$Amat) == x$meq) {
    cat("\nConstrained hypothesis test\n")
  } else {
    cat("\nConstrained hypothesis test type", x$type, "\n")
    if (x$boot != "no") {
      cat("( Number of successful bootstrap draws:", attr(x$pvalue, "B"),")\n")
    }
  }
  vnames <- names(x$b.unconstr)
  if (is.null(vnames)) { 
    vnames <- paste("m", 1:length(x$b.unconstr), sep = "") 
  }
  Amat <- x$Amat
  colnames(Amat) <- vnames
  out.rest <- cbind(Amat, c(rep("   ==", x$meq), rep("   >=", nrow(x$Amat) -
                                                     x$meq)), x$bvec)
  
  rownames(out.rest) <- paste(1:nrow(out.rest), ":", sep = "")
  
  colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs")
  out.rest <- cbind(rep(" ", nrow(out.rest)), out.rest)
  out.rest[x$iact, 1] <- "A"
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
      if (!brief) {
        #cat("\n\nConstraints on", vnames[colSums(!x$Amat == 0) > 0], fill = TRUE)
        cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nConstrained estimate under H0:\n")
        print.default(format(x$b.eqconstr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nConstrained estimate under union of H0 and HA:\n")
        print.default(format(x$b.constr, digits = digits),
                      print.gap = 2, quote = FALSE)
      }
    } else if (x$type == "A") {
      cat("\n H0: all constraints active (=)",
          "\n HA: at least one restriction strictly true (>)",
          "\n\n")
      print(out.test, quote = FALSE, scientific = FALSE)
      if (!brief) {
        #cat("\n\nConstraints on", vnames[colSums(!x$Amat == 0) > 0], fill = TRUE)
        cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
        cat("\nConstrained estimate under H0:\n")
        print.default(format(x$b.eqconstr, digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\nConstrained estimate under union of H0 and HA:\n")
        print.default(format(x$b.constr, digits = digits),
                      print.gap = 2, quote = FALSE)
      }
    }
    if (x$type == "B" && x$meq.alt == 0L) {
        cat("\n H0: all constraints true (>=)",
            "\n HA: at least one constraint violated (<)",
            "\n\n")
        print(out.test, quote = FALSE)
        if (!brief) {
          #cat("\nConstraints on", vnames[colSums(!x$Amat == 0) > 0], "\n")
          cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
          print(out.rest, quote = FALSE, scientific = FALSE)
          cat("\nConstrained estimate under H0:\n")
          print.default(format(x$b.constr, digits = digits),
                        print.gap = 2, quote = FALSE)
          cat("\nUnconstrained estimate:\n")
          print.default(format(x$b.unconstr, digits = digits),
                        print.gap = 2, quote = FALSE)
        }
      }
      else if (x$type == "B" && x$meq.alt > 0L) {
        cat("\n H0: all constraints true (>= or =)",
            "\n HA: at least one constraint violated (<), some =-constraints maintained",
            "\n\n")
        print(out.test, quote = FALSE)
        if (!brief) {
          #cat("\n\nConstraints on", vnames[colSums(!x$Amat == 0) > 0], "\n")
          cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
          print(out.rest, quote = FALSE, scientific = FALSE)
          cat("\nConstrained estimate under H0:\n")
          print.default(format(x$b.constr, digits = digits),
                        print.gap = 2, quote = FALSE)
          cat("\nConstrained estimate under HA:\n")
          print.default(format(x$b.constr.alt, digits = digits),
                        print.gap = 2, quote = FALSE)
        }
      }
      if (x$type == "C") {
        cat("\n H0: at least one constraint not strictly true (<=)",
            "\n HA: all constraints strictly true (>)",
            "\n\n")
        print(out.test, quote = FALSE)
        if (!brief) {
          #cat("\n\nConstraints on", vnames[colSums(!x$Amat == 0) > 0], "\n")
          cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
          print(out.rest, quote = FALSE, scientific = FALSE)
          cat("\nUnconstrained estimate:\n")
          print.default(format(x$b.unconstr, digits = digits),
                        print.gap = 2, quote = FALSE)
        }
      }
  } else {
    cat("\n H0: all constraints active (=)",
        "\n HA: at least one constraint violated (=)",
        "\n\n")
    print(out.test, quote = FALSE)
    if (!brief) {
      #cat("\n\nConstraints on", vnames[colSums(!x$Amat == 0) > 0], fill = TRUE)
      cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
      cat("\nConstrained estimate under H0:\n")
      print.default(format(x$b.constr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\nConstrained estimate under union of H0 and HA:\n")
      print.default(format(x$b.unconstr, digits = digits),
                    print.gap = 2, quote = FALSE)
    }
  }
}

