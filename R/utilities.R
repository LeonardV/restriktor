## utility functions
coef.restriktor <- function(object, ...)  {
  
  b.def <- c()
  b.restr <- object$b.restr
  
  if (any(object$parTable$op == ":=")) {
    b.def <- object$CON$def.function(object$b.restr)
  }
  
  if (inherits(object, "conMLM")) {
    OUT <- rbind(b.restr, b.def)
  } else {
    OUT <- c(b.restr, b.def)
  }
  
  OUT
}


logLik.restriktor <- function(object, ...) {
  object$loglik
}


model.matrix.restriktor <- function(object, ...) {
  model.matrix(object$model.org)
}


tukeyChi <- function(x, c = 4.685061, deriv = 0, ...) {
  u <- x / c
  out <- abs(x) > c
  if (deriv == 0) { # rho function
    r <- 1 - (1 - u^2)^3
    r[out] <- 1
  } else if (deriv == 1) { # rho' = psi function
    r <- 6 * x * (1 - u^2)^2 / c^2
    r[out] <- 0
  } else if (deriv == 2) { # rho'' 
    r <- 6 * (1 - u^2) * (1 - 5 * u^2) / c^2
    r[out] <- 0
  } else {
    stop("deriv must be in {0,1,2}")
  }
  r
}


# code taken from robustbase package.
# slightly addapted by LV (3-12-2017).
robWeights <- function(w, eps = 0.1/length(w), eps1 = 0.001, ...) {
  stopifnot(is.numeric(w))
  cat("Robustness weights:", "\n")
  cat0 <- function(...) cat("", ...)
  n <- length(w)
  if (n <= 10) 
    print(w, digits = 5, ...)
  else {
    n1 <- sum(w1 <- abs(w - 1) < eps1)
    n0 <- sum(w0 <- abs(w) < eps)
    if (any(w0 & w1)) 
      warning("weights should not be both close to 0 and close to 1!\n", 
              "You should use different 'eps' and/or 'eps1'")
    if (n0 > 0 || n1 > 0) {
      if (n0 > 0) {
        formE <- function(e) formatC(e, digits = max(2, 
                                                     5 - 3), width = 1)
        i0 <- which(w0)
        maxw <- max(w[w0])
        c3 <- paste0("with |weight| ", if (maxw == 0) 
          "= 0"
          else paste("<=", formE(maxw)), " ( < ", formE(eps), 
          ");")
        cat0(if (n0 > 1) {
          cc <- sprintf("%d observations c(%s)", n0, 
                        strwrap(paste(i0, collapse = ",")))
          c2 <- " are outliers"
          paste0(cc, if (nchar(cc) + nchar(c2) + nchar(c3) > 
                         getOption("width")) 
            "\n\t", c2)
        }
        else sprintf("observation %d is an outlier", 
                     i0), c3, "\n")
      }
      if (n1 > 0) 
        cat0(ngettext(n1, "one weight is", sprintf("%s%d weights are", 
                                                   if (n1 == n) 
                                                     "All "
                                                   else "", n1)), "~= 1.")
      n.rem <- n - n0 - n1
      if (n.rem <= 0) {
        if (n1 > 0) 
          cat("\n")
        return(invisible())
      }
    }
  }
}
