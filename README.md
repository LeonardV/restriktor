# restriktor

**Restricted statistical estimation and inference for linear models**

[![CRAN status](https://www.r-pkg.org/badges/version/restriktor)](https://CRAN.R-project.org/package=restriktor)
[![License: GPL v2+](https://img.shields.io/badge/License-GPL%20v2%2B-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

`restriktor` is a free, open-source R package for estimating, testing, and
evaluating linear equality and inequality constraints on parameters in linear
models. Typical applications are informative (theory-based) hypotheses such as
`mu1 < mu2 < mu3`, where the substantive theory predicts an ordering of means,
effects, or regression coefficients.

Website: [restriktor.org](https://restriktor.org) · CRAN:
[restriktor](https://CRAN.R-project.org/package=restriktor)

## Features

- **Restricted estimation** — maximum likelihood estimates of linear model
  parameters subject to equality and/or inequality constraints, for `lm`,
  `glm`, `rlm` (robust), and `mlm` (multivariate) models, with a user-friendly
  text-based constraint syntax.
- **Informative hypothesis testing (IHT)** — F-bar, LR-bar, score-bar, and
  Wald-bar tests (`iht()` / `conTest()`), including standard-error types for
  heteroskedasticity (HC) and mixtures of F-distributions for the
  null-distribution (chi-bar-square weights).
- **GORIC and GORICA** — the Generalized Order-Restricted Information
  Criterion (`goric()`) and its approximation for model selection among
  informative hypotheses, including comparison against the complement or the
  unconstrained model, and support for `lavaan` structural equation models.
- **Evidence synthesis** — aggregate GORIC(A) evidence for a central theory
  over multiple, possibly heterogeneous studies (`evSyn()`), with
  leave-one-study-out sensitivity analysis (`leave1studyout()`).
- **Benchmarks** — put GORIC(A) weights into perspective using
  simulation-based benchmarks (`benchmark()`).
- **Utilities** — IC weights from any information criterion
  (`calc_ICweights()`), chi-bar-square weights via simulation
  (`con_weights_boot()`), and parametric bootstrap tests for `lavaan` models
  (`conTestD()`).

## Installation

Install the released version from CRAN:

```r
install.packages("restriktor", dependencies = TRUE)
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("LeonardV/restriktor")
```

`restriktor` requires R >= 4.0.0.

## Getting started

The example below uses the built-in `ZelazoKolb1972` data: the age (in months)
at which infants started walking, for four treatment groups. The theory states
that the active training group walks first, followed by the passive training
group, the control group, and the no-exercise group.

```r
library(restriktor)

# fit the unrestricted ANOVA model (cell means: -1 removes the intercept)
fit <- lm(Age ~ -1 + Group, data = ZelazoKolb1972)

# constraint syntax based on the factor-level names
constraints <- 'GroupActive < GroupPassive < GroupControl < GroupNo'

# restricted estimation
rfit <- restriktor(fit, constraints = constraints)
summary(rfit)

# informative hypothesis test (F-bar test)
iht(rfit)

# GORIC: compare the order-restricted hypothesis against its complement
goric(fit, hypotheses = list(H1 = constraints), comparison = "complement")
```

The constraint syntax also supports equalities (`==`), linear combinations
(e.g. `2*x1 + x2 > 0`), and matrix notation for full control. See
`?restriktor` and `?goric` for details.

## Documentation

- [restriktor.org](https://restriktor.org) — tutorials, examples, and
  background on informative hypothesis testing and the GORIC(A).
- Package vignettes — GORIC(A) tutorials, evidence-synthesis workflows, and
  guidelines for interpreting GORIC(A) output and benchmarks:
  `browseVignettes("restriktor")`.

## Citing restriktor

If you use `restriktor` in your work, please cite it. Run
`citation("restriktor")` in R for the preferred reference.

Key methodological references:

- Silvapulle, M. J., & Sen, P. K. (2005). *Constrained Statistical Inference:
  Order, Inequality, and Shape Constraints.* Wiley.
- Kuiper, R. M., Hoijtink, H., & Silvapulle, M. J. (2011). An Akaike-type
  information criterion for model selection under inequality constraints.
  *Biometrika, 98*(2), 495–501.
- Vanbrabant, L., Van de Schoot, R., & Rosseel, Y. (2015). Constrained
  statistical inference: sample-size tables for ANOVA and regression.
  *Frontiers in Psychology, 5*, 1565.

## Getting help and contributing

- Questions and bug reports: [GitHub issues](https://github.com/LeonardV/restriktor/issues)
  or info@restriktor.org.
- Contributions are welcome — please open an issue first to discuss what you
  would like to change.

## License

GPL (>= 2)
