# Demonstrates benchmark_asymp()'s pop_est auto-derivation when the add_Hc
# safeguard row ("Complement of <ref>") wins overall.
#
# Background: when pop_est is not supplied, benchmark_asymp() derives a
# "No inequalities" population by projecting the sample estimates onto the
# *overall-preferred* hypothesis's constraint boundary (replacing its
# inequality constraints by equalities). If the overall-preferred row is
# itself a safeguard row with no constraint matrix of its own, it falls back
# to a well-defined specified hypothesis instead. For a
# goric(..., add_Hc = ...) object, that fallback is the add_Hc *reference*
# hypothesis itself -- and it happens silently (no message).
#
# This example constructs data where H1 (the add_Hc reference) is badly
# violated, so "Complement of H1" wins overall, and confirms the fallback
# lands on H1's own constraint boundary (x1 = 10), not H2's or H3's.

library(restriktor)

est <- c(x1 = -5.0, x2 = 5.0, x3 = 0.0)
VCOV <- diag(c(0.5, 0.5, 0.5))

H1 <- "x1 > 10"   # badly violated (est x1 = -5)
H2 <- "x2 > 10"   # also badly violated
H3 <- "x3 > 10"   # also badly violated

fit <- goric(est, VCOV = VCOV, hypotheses = list(H1 = H1, H2 = H2, H3 = H3),
             comparison = "unconstrained", add_Hc = "H1", type = "gorica")
print(fit$result)
preferred <- fit$result$model[which.max(fit$result[, 7])]
cat("\npreferred overall:", preferred, "\n\n")
stopifnot(identical(preferred, "Complement of H1"))

msgs <- character(0)
bm <- withCallingHandlers({
  benchmark(fit, iter = 10)
}, message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })

cat("messages captured during benchmark():\n")
print(msgs)
# the add_Hc fallback itself is silent -- any message(s) above are unrelated
# (e.g. about benchmark draw/percentile stability), never about the fallback
stopifnot(!any(grepl("add_Hc|No inequalities.* population is based on", msgs)))

cat("\npop_est used:\n")
print(bm$pop_est)
stopifnot(identical(rownames(bm$pop_est), c("No inequalities", "Observed")))
# the fallback landed on H1's own boundary (x1 = 10), not H2's or H3's
stopifnot(abs(bm$pop_est["No inequalities", "x1"] - 10) < 1e-6)

cat("\nadd_Hc fallback correctly used H1's own constraint boundary, silently\n")
cat("(no fallback-specific message), when the add_Hc safeguard row won overall. OK.\n")
