tukeyChi <- function (x, cc, deriv = 0) 
{
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
