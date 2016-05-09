# acknowledgement: lqs function taken from the MASS package.
Sestimator <- function(x, y, lqs.control= NULL) {
  out <- do.call("lqs",
                  c(list(x, y, intercept = FALSE, method = "S",
                  k0 = 1.54764), lqs.control))
  out
}
