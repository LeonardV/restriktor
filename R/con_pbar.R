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
  



