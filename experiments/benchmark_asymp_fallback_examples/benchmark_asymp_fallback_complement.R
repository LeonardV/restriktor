# Demonstrates benchmark_asymp()'s pop_est auto-derivation for a
# comparison = "complement" object, where the complement safeguard row
# wins overall.
#
# For this case, the fallback silently uses the one specified hypothesis's
# own constraint boundary (again, no message -- unlike the
# comparison = "unconstrained" case, which does message, see
# benchmark_asymp_fallback_unconstrained.R).

library(restriktor)

est <- c(x1 = 5.0, x2 = -5.0, x3 = 0.0)
VCOV <- diag(c(0.5, 0.5, 0.5))

H1 <- "x1 > 10"   # badly violated -> complement should win easily

fit <- goric(est, VCOV = VCOV, hypotheses = list(H1 = H1), comparison = "complement", type = "gorica")
print(fit$result)
preferred <- fit$result$model[which.max(fit$result[, 7])]
cat("\npreferred overall:", preferred, "\n\n")
stopifnot(identical(preferred, "complement"))

msgs <- character(0)
bm <- withCallingHandlers({
  benchmark(fit, iter = 10)
}, message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })

cat("messages captured during benchmark():\n")
print(msgs)
# the complement fallback itself is silent -- any message(s) above are unrelated
stopifnot(!any(grepl("No inequalities.* population is based on", msgs)))

cat("\npop_est used:\n")
print(bm$pop_est)
stopifnot(identical(rownames(bm$pop_est), c("No inequalities", "Observed")))
# the fallback landed on H1's own boundary (x1 = 10) -- the only specified hypothesis
stopifnot(abs(bm$pop_est["No inequalities", "x1"] - 10) < 1e-6)

cat("\nComplement-safeguard case correctly fell back to H1's own constraint boundary,\n")
cat("silently (no fallback-specific message). OK.\n")
