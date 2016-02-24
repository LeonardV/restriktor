print.conTest <- function(x, digits = max(3, getOption("digits") - 4), brief = FALSE, ...) {

  if (!("conTest" %in% class(x))) {
    stop("x must be of class \"conTest\"")
  }
  object <- x
  boot <- object$boot
  cat("Constrained hypothesis test type", object$type, "\n")
  vnames <- names(object$b.constr)
  if (is.null(vnames)) { 
    vnames <- paste("m", 1:length(object$b.constr), sep = "") 
  }
  Amat <- object$Amat
  colnames(Amat) <- vnames
  out.rest <- cbind(Amat, c(rep("   ==", object$meq), rep("   >=", nrow(object$Amat) -
                                                     object$meq)), object$bvec)
  
  rownames(out.rest) <- paste(1:nrow(out.rest), ":", sep = "")
  
  colnames(out.rest)[(ncol(Amat) + 1):ncol(out.rest)] <- c("op", "rhs")
  out.rest <- cbind(rep(" ", nrow(out.rest)), out.rest)
  out.rest[object$iact, 1] <- "A"
  out.rest <- as.data.frame(out.rest)
  names(out.rest)[1] <- ""

  if (object$type != "C" && boot == "none") {
    out.test <- c(format(object$Ts, digits = digits), format(round(object$pvalue$pf.value, 4), digits = digits))
  } else if (object$type != "C" && boot != "none") {
    out.test <- c(format(object$Ts, digits = digits), format(if (object$pvalue < 1e-04) "<0.0001" else round(object$pvalue, 4), digits = 4))
  } else {
    out.test <- c(format(object$Ts, digits = digits), format(if (object$pvalue < 1e-04) "<0.0001" else round(object$pvalue, 4), digits = 4))
  }
  names(out.test) <- c(" Test statistic", "p-value")
  
  if (object$type == "A") {
    cat("\n H0: all constraints active (=)",
        "\n HA: at least one restriction strictly true (>)",
        "\n\n")
    print(out.test, quote = FALSE)
    if (!brief) {
      #cat("\n\nConstraints on", vnames[colSums(!object$Amat == 0) > 0], fill = TRUE)
      cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
      print(out.rest, quote = FALSE, scientific = FALSE)
    }
    cat("\n Constrained estimate under H0:\n")
    print.default(format(object$b.eqconstr, digits = digits),
                  print.gap = 2, quote = FALSE)
    cat("\n Constrained estimate under union of H0 and H1:\n")
    print.default(format(object$b.constr, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
  if (object$type == "B" && object$meq.alt == 0L) {
      cat("\n H0: all constraints true (>=)",
          "\n HA: at least one constraint violated (<)",
          "\n\n")
      print(out.test, quote = FALSE)
      if (!brief) {
        #cat("\nConstraints on", vnames[colSums(!object$Amat == 0) > 0], "\n")
        cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
      }
      cat("\n Constrained estimate under H0:\n")
      print.default(format(object$b.constr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\n Unconstrained estimate:\n")
      print.default(format(object$b.unconstr, digits = digits),
                    print.gap = 2, quote = FALSE)
    }
    else if (object$type == "B" && object$meq.alt > 0L) {
      cat("\n H0: all constraints true (>= or =)",
          "\n HA: at least one constraint violated (<), some =-constraints maintained",
          "\n\n")
      print(out.test, quote = FALSE)
      if (!brief) {
        #cat("\n\nConstraints on", vnames[colSums(!object$Amat == 0) > 0], "\n")
        cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
      }
      cat("\n Constrained estimate under H0:\n")
      print.default(format(object$b.constr, digits = digits),
                    print.gap = 2, quote = FALSE)
      cat("\n Constrained estimate under H1:\n")
      print.default(format(object$b.constr.alt, digits = digits),
                    print.gap = 2, quote = FALSE)
    }
    if (object$type == "C") {
      cat("\n H0: at least one constraint not strictly true (<=)",
          "\n HA: all constraints strictly true (>)",
          "\n\n")
      print(out.test, quote = FALSE)
      if (!brief) {
        #cat("\n\nConstraints on", vnames[colSums(!object$Amat == 0) > 0], "\n")
        cat("\n\n(rows indicated with an \"A\" are active constraints)\n")
        print(out.rest, quote = FALSE, scientific = FALSE)
      }
      cat("\n Unconstrained estimate:\n")
      print.default(format(object$b.unconstr, digits = digits),
                    print.gap = 2, quote = FALSE)
    }
  }


# add overall test
#summary.iclm <- function() {
#
#
#}
