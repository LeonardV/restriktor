Coverage simulation of evSyn_PI() (prediction interval for the ratio of GORICA weights in a new study).
- coverage_PI.R: the simulation script (run: Rscript coverage_PI.R 200; results per condition in res/).
- coverage_PI_summary.csv / coverage_PI.log: results (200 replications per condition).
Design: S = 5, 10, 20 studies with sample sizes drawn from {50, 100, 200, 400}; true values
theta0 = (.25, .12, 0) plus study-specific deviations N(0, tau^2), tau = 0 or .10; hypothesis sets
'comp' (H1: b1 > b2 > b3 vs. its complement) and 'unc' (H1 and H2: b1 > b2 = b3, with the unconstrained
hypothesis; with studies at the ceiling). Per replication, a new study (n = median n) is drawn from the
same population and its observed ratio (best vs. the first other hypothesis) is compared to the PI.
