#acknowledgement: code taken from robustbase package.
tukeyChi <- function (x, cc = 4.685061, deriv = 0) {
  x <- x/cc
  x2 <- x * x
  out <- x2 > 1
  switch(deriv + 1, {
    r <- x2 * (3 + x2 * (-3 + x2))
    r[out] <- 1
  }, {
    r <- 6/cc * x * (1 - x2)^2
    r[out] <- 0
  }, {
    r <- 6/(cc^2) * (1 - x2) * (1 - 5 * x2)
    r[out] <- 0
  }, stop("deriv must be in {0,1,2}"))
  r
}


tukeyPsi <- function (x, cc = 4.685061, deriv = 0) {
  x2 <- (x/cc)^2
  if (deriv < 0) 
    out <- x2 > 1
  else in. <- x2 < 1
  switch(deriv + 2, {
    c. <- cc^2/6
    r <- c. * (1 - (1 - x2)^3)
    r[out] <- c.
    r
  }, {
    in. * x * (1 - x2)^2
  }, {
    in. * (1 - x2) * (1 - 5 * x2)
  }, {
    in. * 4 * x/cc^2 * (5 * x2 - 3)
  }, stop("deriv must be in {-1,0,1,2}"))
}