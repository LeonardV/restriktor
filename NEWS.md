# restriktor 0.7 (development, branch claude/lucid-bell-7cfntx)

## goric()

* Multivariate linear models (`mlm`, e.g. `lm(cbind(y1, y2) ~ x)`) are
  supported: the unrestricted log-likelihood is the multivariate-normal
  log-likelihood, the sample size is `N` (not `N` times the number of
  responses), and the coefficients are handled as a named vector in the
  order of `vcov()`. Hypotheses use the names of `vcov()` with `:` replaced
  by `.` (e.g. `y1.x`, `y1..Intercept.`). `type = "goricc"`/`"goricac"` is
  not available for `mlm` objects until the multivariate small-sample
  correction has been derived.
* New argument `priorICweights`: prior weights for the models in the result
  table, validated (finite, non-negative, positive sum, correct length),
  matched by name when named, and applied via a numerically stable
  log-sum-exp (a zero prior gives weight 0, not `NaN`).
* IC weights, log-likelihood weights and penalty weights are computed on the
  log scale (no underflow for large IC differences).
* `Heq = TRUE`: the equality-restricted hypothesis is formed after removing
  redundant inequalities; `Heq` is ignored (with a message) when more than
  one hypothesis is given or when no inequality remains; infeasible or
  range-type `Heq` constructions give a clear error.
* Equality restrictions are always kept first in the constraint matrix
  (previously an inequality could be imposed as an equality after
  de-duplication); duplicated and conflicting restrictions are detected,
  conflicting ones give an error.
* Defined parameters (`:=`) are included in the coefficient table for all
  input routes.
* Standardized lavaan estimates are merged on parameter identity (defined
  parameters are no longer dropped).
* Hypothesis names `"Heq"`, `"complement"` and `"unconstrained"` are
  reserved and refused.
* `print()`/`summary()`: the "times more supported" sentences index the
  ratio matrix by name (correct with ties and with `Heq`), `summary()` works
  for a single hypothesis with `comparison = "none"`, and the best model is
  marked only in the printed output (row names of the ratio matrices are
  unchanged).
* Removed the unimplemented arguments `add_Hc` and `posthoc`.

## evSyn()

* New arguments `priorICweights` and `study_weights` (study weights are
  reordered with `order_studies`, zero weights are supported, and both are
  validated and matched by name when named). `priorWeights` is deprecated.
* One hypothesis is compared with its complement again (also for a list of
  one-hypothesis sets).
* Hypotheses are aligned across studies by name (goric objects, LL/PT and
  the est route); differing sets, differing numbers of hypotheses,
  differing `comparison` or `penalty_factor` give an error.
* The input type is detected with a stricter rule (IC ratios only when all
  studies have a 1 at the same position; ambiguous input requires
  `input_type`), and the chosen route is always reported.
* Products of IC weights are formed on the log scale (no underflow), a
  single study works in all routes, `type_ev = "average"` is consistent
  across routes, and `leave1studyout()` uses the same weighting and the
  prior-weighted preferred hypothesis (also for IC weights/ratios objects).
* `order_studies` accepts study names.

## benchmark()

* The simulation uses the same criterion (GORICA/GORICAC), sample size,
  prior weights and penalty factor as the benchmarked object.
* Population means are generated from the centred observed pattern (or
  `ratio_pop_means`) with a positive scale factor, so the ordering of the
  observed means is preserved; the 'Observed' population uses the observed
  estimates exactly.
* Cohen's f uses the residual variance `RSS/N` (consistent with the
  simulated covariance matrix); covariates and additional factors are kept
  at their observed values.
* Ratios of weights are formed on the log scale (`rgw_log`, `rlw_log`);
  non-finite draws are handled explicitly (overlap is `NA` with a note).
* Warnings inside a draw are muffled and the draw is kept; failed draws are
  counted and reported (`iter` is the number of successful draws,
  `iter_requested` the requested number).
* Adaptive `iter` stops after two consecutive stable rounds; the `iter_*`
  arguments and thresholds are validated and documented.
* Duplicate population names, vector sample sizes with
  `alt_sample_size`, ANCOVA models with a second factor, and intercept
  models are handled or refused with a clear message.
* `benchmark()` no longer changes the global `progressr` handlers.

## Other

* `restriktor()`: standard errors via HAC covariance matrices (`se =
  "HAC"`, `"kernHAC"`, `"NeweyWest"`).
* Many documentation corrections; the benchmark vignette builds again.
