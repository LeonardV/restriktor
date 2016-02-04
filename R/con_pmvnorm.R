con_pmvnorm <- function(Ts = NULL, df1 = NULL, df2, weights = NULL, ...) {
  
  #acknowlegment: taken from ic.infer package (Ulrike Groemping)
  pfbar <- function(x, df1, df2, wt) {
    if (x <= 0) {
      return(0)
    }
    zed <- df1 == 0
    cdf <- ifelse(any(zed), wt[zed], 0)
    cdf <- cdf + sum(pf(x/df1[!zed], df1[!zed], df2) * wt[!zed])
    return(cdf)
  }
  # p value based on F-distribution
  pf.value <- 1 - pfbar(x = Ts, df1 = df1, df2 = df2, wt = weights)
  
  # acknowlegment: taken from ic.infer package (Ulrike Groemping)
  pchibar <- function (x, df1, wt) 
  {
    if (x <= 0) {
      return(0)
    }
    zed <- df1 == 0
    cdf <- ifelse(any(zed), wt[zed], 0)
    cdf <- cdf + sum(pchisq(x, df1[!zed]) * wt[!zed])
    return(cdf)
  }
  # p value based on chi^2-distribution
  pchi.value <- 1 - pchibar(x = Ts, df1 = df1, wt = weights)

  OUT <- list(pf.value = pf.value, pchi.value = pchi.value)

  return(OUT)
}



