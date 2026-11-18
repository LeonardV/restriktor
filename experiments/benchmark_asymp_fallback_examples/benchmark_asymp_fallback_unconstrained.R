# Demonstrates benchmark_asymp()'s pop_est auto-derivation for a
# comparison = "unconstrained" object, where the unconstrained model wins
# overall.
#
# Unlike the add_Hc and complement safeguard cases, this fallback DOES
# message -- since it may not be obvious to the user that the "No
# inequalities" population ends up derived from one of the *specified*
# hypotheses (the best-supported of them) rather than from the preferred
# model itself. The message is deliberately generic: it reports that the
# preferred hypothesis is the unconstrained model and that the "No
# inequalities" population is based on the best of the specified
# hypotheses, but it does NOT name which one that is.

library(restriktor)

est <- c(x1 = 5.0, x2 = -5.0, x3 = 0.0)
VCOV <- diag(c(0.5, 0.5, 0.5))

H1 <- "x1 > 10"   # badly violated
H2 <- "x2 > -4.0" # best of the specified (moderate violation)
H3 <- "x3 > 8"    # badly violated

fit <- goric(est, VCOV = VCOV, hypotheses = list(H1 = H1, H2 = H2, H3 = H3),
             comparison = "unconstrained", type = "gorica")
cat("Confirmed setup: preferred overall =", fit$result$model[which.max(fit$result[, 7])],
    "; best among specified =", fit$objectNames[which.max(fit$result[seq_along(fit$constraints), 7])], "\n\n")

cat("=== Calling benchmark() (pop_est = NULL -> triggers fallback) ===\n")
msgs <- character(0)
bm <- withCallingHandlers({
  benchmark(fit, iter = 10)
}, message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })

cat("\n--- messages seen ---\n")
cat(paste(msgs, collapse = "\n---\n"), "\n")

fallback_msg <- msgs[grepl("preferred hypothesis is the unconstrained", msgs)]
stopifnot(length(fallback_msg) == 1)
stopifnot(grepl("'No inequalities' population is based on the best hypothesis of the specified ones", fallback_msg))
# the message is deliberately generic and does not name which specified
# hypothesis (H1/H2/H3) was actually used as the fallback
stopifnot(!grepl("'H1'", fallback_msg), !grepl("'H2'", fallback_msg), !grepl("'H3'", fallback_msg))
cat("\nFallback message correctly reports that the preferred hypothesis is unconstrained\n")
cat("and that the 'No inequalities' population is based on the best hypothesis of the\n")
cat("specified ones (without naming it). OK.\n")

cat("\npop_est used (No inequalities row should come from projecting onto H2's constraint,\n")
cat("since H2 is the best-supported of H1/H2/H3):\n")
print(bm$pop_est)
stopifnot(identical(rownames(bm$pop_est), c("No inequalities", "Observed")))
