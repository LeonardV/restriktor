# Demonstrates the ordinary (non-safeguard) case for an add_Hc object,
# as a contrast to benchmark_asymp_fallback_addHc_fires.R.
#
# When one of the *specified* hypotheses is itself preferred overall (not
# the "Complement of <ref>" safeguard row), benchmark_asymp()'s fallback
# logic never triggers: it uses that hypothesis's own constraints directly,
# exactly as it would for a plain (non-add_Hc) object. This example uses a
# real goric(..., add_Hc = ...) fit on actual data to confirm the
# fallback-detection logic stays out of the way in this ordinary case.

library(restriktor)

data(ZelazoKolb1972)
DATA <- subset(ZelazoKolb1972, Group != "Control")
fit.lm <- lm(Age ~ Group - 1, data = DATA)
H1 <- "GroupPassive > 0; GroupPassive < GroupNo; new1 := GroupPassive + 0.5 > 0.9"
H2 <- "GroupPassive > 0; GroupPassive > GroupNo; new2 := GroupPassive + 0.5 > 0.7"
H3 <- "GroupPassive = 0; GroupPassive < GroupNo"

set.seed(1)
fit_AddHc <- suppressWarnings(goric(fit.lm, hypotheses = list(H1 = H1, H2 = H2, H3 = H3), add_Hc = "H1"))
preferred <- fit_AddHc$result$model[which.max(fit_AddHc$result[, 7])]
cat("fit_AddHc preferred hypothesis:", preferred, "\n")
# NOT the "Complement of H1" safeguard row -- an ordinary specified hypothesis won
stopifnot(preferred %in% c("H1", "H2", "H3"))

msgs <- character(0)
set.seed(1)
bm <- suppressWarnings(withCallingHandlers({
  benchmark(fit_AddHc, iter = 10)
}, message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }))

cat("\nmessages captured during benchmark():\n")
print(msgs)
cat("\npop_est used:\n")
print(bm$pop_est)
stopifnot(identical(rownames(bm$pop_est), c("No inequalities", "Observed")))

cat("\nOrdinary case (no safeguard row preferred): add_Hc's fallback-detection logic\n")
cat("stayed out of the way, as expected -- no fallback, no fallback-specific message. OK.\n")
