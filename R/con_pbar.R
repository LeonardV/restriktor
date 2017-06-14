## acknowlegment: taken from package ic.infer
## pfbar added by LV.

pfbar <- function(x, df1, df2, wt.bar) {
  if (x <= 0) {
    return(0)
  }
  zed <- df1 == 0
  cdf <- ifelse(any(zed), wt.bar[zed], 0)
  cdf <- cdf + sum(pf(x/df1[!zed], df1[!zed], df2) * wt.bar[!zed])
  return(cdf)
}
  

# pchibar <- function (x, df1, wt.bar) {
#   if (x <= 0) {
#     return(0)
#   }
#   zed <- df1 == 0
#   cdf <- ifelse(any(zed), wt.bar[zed], 0)
#   cdf <- cdf + sum(pchisq(x, df1[!zed]) * wt.bar[!zed])
#   return(cdf)
# }

  
  