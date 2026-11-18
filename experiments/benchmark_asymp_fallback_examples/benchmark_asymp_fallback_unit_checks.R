# Regression test for the pop_est auto-derivation fallback in
# benchmark_asymp() (goric_benchmark.R).
#
# Background: when the user doesn't supply pop_est, benchmark_asymp() derives
# a 'No inequalities' population by projecting the observed estimates onto the
# *preferred* hypothesis's constraint boundary (theta_restricted(), treating
# all its restrictions as equalities). This needs object$constraints[[NE_hypo]]
# / object$rhs[[NE_hypo]] -- but those lists only have one entry per
# *specified* hypothesis, never one for an implicit safeguard row. So this
# breaks (subscript out of bounds) whenever the preferred row is a safeguard:
# the complement (an ordinary comparison = "complement" object, always
# exactly 1 specified hypothesis), the unconstrained model (an ordinary
# comparison = "unconstrained" object), or the "Complement of <add_Hc>" row
# of a goric(..., add_Hc = ...) object. The fix detects NE_hypo pointing
# past object$constraints and falls back to a well-defined reference
# hypothesis:
#  - add_Hc: the add_Hc reference hypothesis itself (no message).
#  - comparison = "complement": the (only) specified hypothesis (no message).
#  - comparison = "unconstrained": the BEST-SUPPORTED of the specified
#    hypotheses (highest GORIC(A) weight among just those, excluding the
#    unconstrained row) -- NOT simply "the first" one. This case alone gets
#    a message(), since it may not be obvious/interesting to the user that
#    the "No inequalities" population is derived from one of the specified
#    hypotheses rather than the (preferred) unconstrained model; the message
#    is deliberately generic and does NOT name which specified hypothesis was
#    used (confirmed design: add_Hc/complement fall back silently, only the
#    unconstrained case messages, and that message never names a hypothesis).
#
# See the sibling scripts in this folder for the same three cases exercised
# through the real goric() + benchmark() pipeline on actual/simulated data:
#   benchmark_asymp_fallback_addHc_fires.R
#   benchmark_asymp_fallback_addHc_not_triggered.R  (contrast: no fallback needed)
#   benchmark_asymp_fallback_complement.R
#   benchmark_asymp_fallback_unconstrained.R
#
# This script instead extracts the pop_est-derivation block verbatim from the
# real, currently-shipped goric_benchmark.R by text-slicing between two
# unique anchor lines, then eval()'s it as a small standalone function -- so
# it tests the actual shipped code directly against hand-crafted edge cases,
# without needing data that happens to produce each scenario.
#
# Run this script with its own folder as the working directory (e.g. via
# RStudio's "Source" button, or `Rscript benchmark_asymp_fallback_unit_checks.R`
# from inside this folder) -- it locates the package's R/goric_benchmark.R
# two directories up (experiments/<this folder>/../../R).

pkg_root <- normalizePath(file.path("..", ".."))
bm_lines <- readLines(file.path(pkg_root, "R", "goric_benchmark.R"))
# benchmark_asymp() is defined a second time (after its S3-generic stub at
# the top of the file); "NE_hypo <- which.max(...)" also appears once in
# benchmark_means() and a second time later in benchmark_asymp() itself (for
# the simulation setup, unrelated to pop_est) -- so anchor the search to
# start after benchmark_asymp()'s real definition and take the first match
# after that point, which is the one inside the pop_est-derivation block.
asymp_def_idx <- grep("^benchmark_asymp <- function\\(object, pop_est", bm_lines)[1]
stopifnot(!is.na(asymp_def_idx))
candidate_starts <- grep("NE_hypo <- which.max\\(object\\$result\\[, 7\\]\\)", bm_lines)
start_idx <- min(candidate_starts[candidate_starts > asymp_def_idx])
candidate_ends <- grep("row.names\\(pop_est\\) <- c\\(\"No inequalities\", \"Observed\"\\)", bm_lines)
end_idx <- min(candidate_ends[candidate_ends > start_idx])
stopifnot(!is.na(start_idx), !is.na(end_idx), end_idx > start_idx)
block <- bm_lines[start_idx:end_idx]

# Sanity: confirm the extracted block actually contains the new fallback
# logic (not an unrelated match, and not a stale pre-update copy), before
# trusting the eval'd test below.
stopifnot(any(grepl("NE_hypo > length\\(object\\$constraints\\)", block)))
stopifnot(any(grepl("add_Hc", block)))
stopifnot(any(grepl("best_among_specified", block)))
stopifnot(any(grepl("identical\\(object\\$comparison, \"complement\"\\)", block)))

fn_src <- c(
  "derive_pop_est <- function(object, est_sample, VCOV, theta_restricted) {",
  block,
  "  pop_est",
  "}"
)
eval(parse(text = fn_src))

cat("=== Part 0: extracted block contains the new guard (sanity check) ===\n")
cat("OK -- testing the actual shipped code block.\n\n")

# a trivial theta_restricted() stand-in: records which (R, rhs) it was
# called with (by name, via a side channel) and returns a recognizable value
last_used_constraints <- NULL
mock_theta_restricted <- function(theta, V, R, rhs) {
  last_used_constraints <<- list(R = R, rhs = rhs)
  theta - 1  # arbitrary, just needs to be distinguishable
}

cat("=== Part 1: ordinary case -- NE_hypo within range, no fallback needed ===\n")
object1 <- list(
  result = data.frame(goric.weights = c(0.6, 0.3, 0.1)),  # col 7 stand-in: see note below
  constraints = list(H1 = "R1", H2 = "R2"),
  rhs = list(H1 = "rhs1", H2 = "rhs2"),
  objectNames = c("H1", "H2"),
  comparison = "unconstrained",
  add_Hc = NULL
)
# object$result[, 7] is accessed positionally; build a 7-col df so column 7
# holds the weights used to pick NE_hypo (mirrors the real $result shape).
object1$result <- as.data.frame(matrix(NA, nrow = 2, ncol = 7))
object1$result[, 7] <- c(0.3, 0.7)  # H2 preferred (NE_hypo = 2), in range

out1 <- withCallingHandlers(
  derive_pop_est(object1, est_sample = c(1, 2), VCOV = diag(2), theta_restricted = mock_theta_restricted),
  message = function(m) stop("unexpected message (no fallback should fire here): ", conditionMessage(m))
)
stopifnot(identical(last_used_constraints$R, "R2"), identical(last_used_constraints$rhs, "rhs2"))
cat("H2 preferred (in range): used H2's own constraints directly, no fallback message. OK.\n\n")

cat("=== Part 2: comparison = 'unconstrained', 2+ hypotheses -- falls back to the BEST\n")
cat("    among the specified ones (here H2), not simply hypothesis 1, with a generic\n")
cat("    message that does NOT name which hypothesis was used ===\n")
object2 <- list(
  result = as.data.frame(matrix(NA, nrow = 4, ncol = 7)),
  constraints = list(H1 = "R1", H2 = "R2", H3 = "R3"),  # 3 entries -- no entry for "unconstrained"
  rhs = list(H1 = "rhs1", H2 = "rhs2", H3 = "rhs3"),
  objectNames = c("H1", "H2", "H3"),
  comparison = "unconstrained",
  add_Hc = NULL
)
# row 4 ("unconstrained") preferred overall -> NE_hypo = 4, out of range;
# among the 3 specified rows, H2 (row 2) has the highest weight (0.35) --
# NOT H1 (row 1, 0.10) -- so the fallback must pick H2, not "the first".
object2$result[, 7] <- c(0.10, 0.35, 0.05, 0.50)

msgs2 <- character(0)
out2 <- withCallingHandlers(
  derive_pop_est(object2, est_sample = c(1, 2), VCOV = diag(2), theta_restricted = mock_theta_restricted),
  message = function(m) { msgs2 <<- c(msgs2, conditionMessage(m)); invokeRestart("muffleMessage") }
)
stopifnot(identical(last_used_constraints$R, "R2"), identical(last_used_constraints$rhs, "rhs2"))
stopifnot(length(msgs2) == 1)
stopifnot(grepl("preferred hypothesis is the unconstrained", msgs2[1]))
stopifnot(grepl("'No inequalities' population is based on the best hypothesis of the specified ones", msgs2[1]))
# the message is deliberately generic: it never names which specified
# hypothesis (H1/H2/H3) was actually used as the fallback
stopifnot(!grepl("'H1'", msgs2[1]), !grepl("'H2'", msgs2[1]), !grepl("'H3'", msgs2[1]))
cat("Unconstrained preferred overall (NE_hypo=4, out of range): fell back to H2's\n")
cat("constraints (the best-supported of H1/H2/H3, verified directly via the mock), with\n")
cat("the generic message that does not name a specific hypothesis. OK.\n\n")

cat("=== Part 3: comparison = 'complement' (single hypothesis) -- falls back to that\n")
cat("    one hypothesis silently (no message) ===\n")
object2b <- list(
  result = as.data.frame(matrix(NA, nrow = 2, ncol = 7)),
  constraints = list(H1 = "R1"),  # exactly one entry -- no entry for "complement"
  rhs = list(H1 = "rhs1"),
  objectNames = c("H1"),
  comparison = "complement",
  add_Hc = NULL
)
object2b$result[, 7] <- c(0.05, 0.95)  # row 2 ("complement") preferred -> NE_hypo = 2, out of range

msgs2b <- character(0)
out2b <- withCallingHandlers(
  derive_pop_est(object2b, est_sample = c(1, 2), VCOV = diag(2), theta_restricted = mock_theta_restricted),
  message = function(m) { msgs2b <<- c(msgs2b, conditionMessage(m)); invokeRestart("muffleMessage") }
)
stopifnot(identical(last_used_constraints$R, "R1"), identical(last_used_constraints$rhs, "rhs1"))
stopifnot(length(msgs2b) == 0)
cat("Complement preferred (NE_hypo=2, out of range): fell back to H1 (the only\n")
cat("specified hypothesis) silently, with no message. OK.\n\n")

cat("=== Part 4: add_Hc complement wins -- falls back to the add_Hc reference hypothesis,\n")
cat("    silently (no message) ===\n")
object3 <- list(
  result = as.data.frame(matrix(NA, nrow = 4, ncol = 7)),
  constraints = list(H1 = "R_H1", H2 = "R_H2", H3 = "R_H3"),  # 3 entries, no entry for the complement row
  rhs = list(H1 = "rhs_H1", H2 = "rhs_H2", H3 = "rhs_H3"),
  objectNames = c("H1", "H2", "H3"),
  comparison = "none",
  add_Hc = "H2"  # reference hypothesis is H2, not H1 -- exercises the non-default case too
)
object3$result[, 7] <- c(0.1, 0.2, 0.15, 0.55)  # row 4 ("Complement of H2") preferred -> NE_hypo = 4

msgs3 <- character(0)
out3 <- withCallingHandlers(
  derive_pop_est(object3, est_sample = c(1, 2), VCOV = diag(2), theta_restricted = mock_theta_restricted),
  message = function(m) { msgs3 <<- c(msgs3, conditionMessage(m)); invokeRestart("muffleMessage") }
)
# must fall back to H2 (the add_Hc reference), NOT H1, and NOT "best among
# specified" (that reasoning is for the plain 'unconstrained' case only)
stopifnot(identical(last_used_constraints$R, "R_H2"), identical(last_used_constraints$rhs, "rhs_H2"))
stopifnot(length(msgs3) == 0)
cat("add_Hc = 'H2' complement preferred (NE_hypo=4, out of range): fell back to H2's (the\n")
cat("reference hypothesis's) own constraints -- not H1, not 'best among specified' --\n")
cat("silently, with no message. OK.\n\n")

cat("=== Part 5: pop_est row names and shape are still correct in all fallback cases ===\n")
for (out in list(out2, out2b, out3)) {
  stopifnot(identical(rownames(out), c("No inequalities", "Observed")))
  stopifnot(nrow(out) == 2)
}
cat("pop_est output shape/row names unaffected by the fallback path. OK.\n\n")

cat("BENCHMARK_ASYMP POP_EST FALLBACK TEST OK.\n")
