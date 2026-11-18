# Forest-type plots (evSyn_forest()) for evidence syntheses with mixed evidence:
# some studies support the overall best hypothesis, others do not.
library(restriktor)
library(ggplot2)

n <- 200                       # observations per group in each study
V <- diag(1 / n, 3)            # covariance matrix of the three estimated group means
draw_est <- function(theta) setNames(theta + rnorm(3, 0, sqrt(1 / n)), c("b1", "b2", "b3"))

# -----------------------------------------------------------------------------
# Example 1: H1 versus its complement
# Four studies with true means in line with H1 (b1 > b2 > b3) and two studies
# with true means in the opposite order. These two studies give strong evidence
# for the complement, so that the complement is the overall best hypothesis,
# while the four other studies (weakly) favour H1.
# -----------------------------------------------------------------------------
set.seed(123)
theta_1 <- list(c(.45, .25, 0), c(.40, .20, 0), c(.35, .20, .05), c(.40, .15, 0),
                c(0, .30, .60), c(0, .25, .55))
est_1 <- lapply(theta_1, draw_est)
names_1 <- paste("Study", 1:6)
res_1 <- evSyn(est_1, VCOV = rep(list(V), 6), hypotheses = list(H1 = "b1 > b2 > b3"),
               study_names = names_1, CI = TRUE, variability = TRUE)
print(res_1)
p1 <- evSyn_forest(res_1, seed = 1)
ggsave("forest_mixed_H1_vs_complement.png", p1, width = 8, height = 5, dpi = 120)

# -----------------------------------------------------------------------------
# Example 2: H1 and H2, with the unconstrained hypothesis as failsafe
# Three studies with strictly ordered true means (in line with H1, not H2),
# two studies with b2 = b3 (in line with H2, and on the boundary of H1), and
# one study with true means outside both hypotheses.
# -----------------------------------------------------------------------------
set.seed(2024)
theta_2 <- list(c(.45, .25, 0), c(.50, .25, 0), c(.40, .20, 0),
                c(.40, 0, 0), c(.35, 0, 0),
                c(0, .30, .10))
est_2 <- lapply(theta_2, draw_est)
names_2 <- paste("Study", 1:6)
res_2 <- evSyn(est_2, VCOV = rep(list(V), 6),
               hypotheses = list(H1 = "b1 > b2 > b3", H2 = "b1 > b2 = b3"),
               study_names = names_2, CI = TRUE, variability = TRUE)
print(res_2)
p2 <- evSyn_forest(res_2, seed = 1)
ggsave("forest_mixed_H1_H2_unconstrained.png", p2, width = 10, height = 5, dpi = 120)
