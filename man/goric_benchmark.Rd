\name{benchmark_functions}
\alias{benchmark}
\alias{benchmark_means}
\alias{benchmark_asymp}
\alias{print.benchmark}
\alias{plot.benchmark}
\title{Benchmark Functions for GORIC(A) Analysis}
\description{
  The `benchmark` functions perform benchmarking for models using the 
  Generalized Order-Restricted Information Criterion (Approximation) (GORIC(A)). 
}
\usage{
  benchmark(object, model_type = c("asymp", "means"), \dots)
  
  benchmark_means(object, pop_es = NULL, ratio_pop_means = NULL, 
                  group_size = NULL, alt_group_size = NULL, 
                  quant = NULL, iter = NULL, 
                  control = list(), 
                  ncpus = 1, seed = NULL,
                  iter_adequacy_band = c(0.495, 0.505),
                  iter_stability_tol = 1,
                  iter_min = 500, iter_step = 100, iter_max = 2000,
                  hypo_rate_threshold = 1,
                  threshold_rlw = 1, priorICweights = NULL, \dots)

  benchmark_asymp(object, pop_est = NULL, sample_size = NULL, 
                  alt_sample_size = NULL, quant = NULL, iter = NULL, 
                  control = list(), 
                  ncpus = 1, seed = NULL,
                  iter_adequacy_band = c(0.495, 0.505),
                  iter_stability_tol = 1,
                  iter_min = 500, iter_step = 100, iter_max = 2000,
                  hypo_rate_threshold = 1,
                  threshold_rlw = 1, priorICweights = NULL, \dots)
  
  \method{print}{benchmark}(x, output_type = c("rgw", "gw", "lw", "rlw", "ld",
                                               "rgw_log", "rlw_log", "all"),
                            hypo_rate_threshold = NULL, threshold_rlw = NULL,
                            color = TRUE, percentiles = NULL, ...)
  
  \method{plot}{benchmark}(x, output_type = c("rgw", "rlw", "gw", "ld"), 
                           percentiles = NULL, x_lim = c(), log_scale = FALSE,
                           opacity = 0.50, nrow_grid = NULL, ncol_grid = 1, 
                           distr_grid = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{con_goric} (a GORIC(A) object from the \code{goric} function).}
  \item{model_type}{If "means", the benchmarks are based on the (adjusted) group
  means, else model_type = "asymp" (default). For "means" with a GORIC(A) object
  based on an \code{lm} object, the model may contain an intercept, multiple
  factors and/or covariates: the group means are then derived from the
  parameters (see \code{\link{group_means}}). See details for more information.}
  \item{x}{An object of class \code{benchmark} (i.e., of class
  \code{benchmark_means} or \code{benchmark_asymp}).}
  \item{pop_es}{A scalar or a vector of population Cohen's f (effect-size) values. 
  By default, it benchmarks ES = 0 (no-effect) and the observed Cohen's f.}
  \item{pop_est}{A 1 x k vector or an n x k matrix of population 
  estimates to benchmark. By default, estimates are determined based on replacing 
  the inequalities by equalities for the best hypothesis (no-effect) 
  and the observed estimates from the sample are used.}
  \item{ratio_pop_means}{A 1 x k vector denoting the relative difference
  between the k group means. Note that a ratio of \code{c(3, 2, 1)} gives the 
  same as \code{c(1, 0, -1)}, as the consecutive relative differences are 1 in 
  both ratios. By default, the relative differences from the data are used.
  The order of the groups is that of \code{group_means(model)} (in the output:
  \code{names(group_means_observed)}).}
  \item{group_size}{If the GORICA object is based on estimates and their covariance 
  matrix (instead of on a model/fit object), this should be a 1 x k vector 
  or a scalar to denote the group sizes. If a scalar is specified, it is assumed 
  that each group is of that size.}
  \item{alt_group_size}{An 1 x k vector or a scalar to denote alternative 
  group sizes, if you want to use sizes different from those in the data. This 
  can be used, for example, to see the values to which the GORIC(A) weights will 
  converge (and thus to see the maximum value of the weights). If a scalar is 
  specified, it is assumed that each group is of that size. By default, the group 
  sizes from the data are used.}
  \item{sample_size}{A scalar to denote the (total) sample sizes. Only used if 
  the GORIC object is based on estimates and their covariance matrix (instead 
  of on a model/fit object) or \code{alt_sample_size} is not \code{NULL}. }
  \item{alt_sample_size}{A scalar to denote an alternative sample size if you 
  want to use a different sample size from the one in the data. This can be used, 
  for example, to see the values to which the GORIC(A) weights will converge 
  (and thus to see the maximum value of the weights).}
  \item{quant}{Quantiles for benchmarking results. Defaults 5\%, 35\%, 50\%, 65\%, 95\%.}
  \item{iter}{The number of iterations (draws) for benchmarking. By default
  (\code{NULL}), the number of draws is set automatically: the benchmark starts with
  \code{iter_min} draws and adds \code{iter_step} draws at a time (keeping the earlier
  ones), until the percentile of the sample value within the benchmark distribution of
  the \code{"Observed"} population has stabilized (see \code{iter_stability_tol}) or
  \code{iter_max} draws are reached. If a number is given, exactly that many draws are
  used; a message is then given if that percentile has not (yet) stabilized.}
  \item{iter_min}{The number of draws to start with when \code{iter = NULL}. Defaults to
  \code{500}.}
  \item{iter_step}{The number of draws added per round when \code{iter = NULL}.
  Defaults to \code{100}.}
  \item{iter_max}{The maximum number of draws when \code{iter = NULL}. Defaults to
  \code{2000}.}
  \item{iter_stability_tol}{The tolerance, in percentage points, used to decide whether
  the percentile has stabilized: this is the case when, for both the GORIC(A) weights
  and the log-likelihood weights, the percentile of the sample value changed by less
  than \code{iter_stability_tol} after the last \code{iter_step} draws (when
  \code{iter = NULL}) or over the last 20\% of the draws (when \code{iter} is given).
  Defaults to \code{1}.}
  \item{iter_adequacy_band}{A numeric vector of length 2 with the bounds of the band
  around the median (0.5) used in an informational check of whether the sample value
  lies at the median of the benchmark distribution of the \code{"Observed"} population
  (returned as \code{median_bias_check_gw} and \code{median_bias_check_lw}). It is not
  used to decide on the number of draws. Defaults to \code{c(0.495, 0.505)}.}
  \item{hypo_rate_threshold}{A numeric value specifying the threshold for the 
  hypothesis rate. The function calculates the proportion of ratio-of-goric-weights 
  that exceeds this threshold. Defaults to 1 in \code{benchmark_means} and
  \code{benchmark_asymp}. In \code{print}, the default (\code{NULL}) uses the threshold
  stored in the object; another value recomputes the rate for that print call only,
  without rerunning the benchmark. See the note on the hypothesis rate.}
  \item{threshold_rlw}{A numeric value specifying the threshold for the rlw rate: the
  proportion of ratio-of-log-likelihood-weights that exceeds this threshold. Defaults
  to 1 in \code{benchmark_means} and \code{benchmark_asymp}. In \code{print}, the
  default (\code{NULL}) uses the threshold stored in the object; another value
  recomputes the rate for that print call only. See the note on the hypothesis rate.}
  \item{control}{A list of control parameters. For more information, see details \link{goric}.}
  \item{ncpus}{Number of CPUs to use for parallel processing. Defaults to \code{1}.
  See details for more information.}
  \item{seed}{A seed for random number generation.}
  \item{priorICweights}{Optional prior model probabilities used in the re-calculated and simulated 
  GORIC(A) values (and thus in the benchmark GORIC(A) weights and ratios of weights). 
  By default (\code{NULL}), the ones of \code{object} are used (i.e., \code{object$priorICweights}, 
  which are equal weights unless specified otherwise in \code{goric}), such that the benchmarks 
  are based on the same prior model probabilities as the sample values. 
  See \code{priorICweights} in \link{goric} for the required number of elements.}
  \item{output_type}{A character vector specifying the type of output to print 
  or plot. Options are \code{"all"}, \code{"gw"} (goric(a) weights), 
  \code{"rgw"} (ratio of goric(a) weights), \code{"rlw"} (ratio of log-likelihood values), 
  and \code{"ld"} (log-likelihood difference). Defaults to \code{"rgw"} for print 
  and \code{"rgw"} for plot.}
  \item{color}{If TRUE, the output will include ANSI color coding. Set \code{color = FALSE} 
  when using this function in R Markdown documents to avoid rendering issues with color codes.}
  \item{opacity}{The opacity of the density curves in \code{plot}. Values range from 0
  to 1, with lower values corresponding to more transparent colors. Defaults to
  \code{0.50}. (This argument used to be called \code{alpha}; \code{alpha} still works
  for now, with a deprecation warning, but may be removed in a future version.)}
  \item{nrow_grid}{An integer value representing the number of rows in the grid layout.}
  \item{ncol_grid}{An integer value representing the number of columns in the grid layout.}
  \item{distr_grid}{If TRUE, the facet_grid function is used to create a grid of separate plots for each effect-size (estimates).}
  \item{percentiles}{A numeric vector specifying the percentiles to be shown (by
  \code{print} and \code{plot}). By default (\code{NULL}) the percentiles are inherited
  from the quantiles used for benchmarking, see \code{quant}; other values are
  recomputed from the stored draws, without rerunning the benchmark.}
  \item{x_lim}{A numeric vector of length 2 specifying the x-axis limits. Defaults to \code{c()}.}
  \item{log_scale}{logical, If TRUE, The x-axis is transformed using a base-10 logarithmic scale.
        This transformation adjusts the way the data is visualized on the x-axis, but does not 
        alter the underlying data values themselves.}
  \item{\ldots}{See \link{goric}.}
}

\details{
  The function \code{benchmark_asymp} is named as such because it generates data from a 
  multivariate normal distribution with means equal to the population parameter 
  estimates and a covariance matrix derived from the original data. This is based 
  on the assumption that parameter estimates are asymptotically normally distributed. 
  This assumption is valid for many statistical models, including parameters from 
  a generalized linear model (GLM). In such models, as the sample size increases, 
  the distribution of the parameter estimates tends to a normal distribution, 
  allowing us to utilize the multivariate normal distribution for benchmarking.

  \code{benchmark_means} benchmarks the group means of a given GORIC(A) object 
  by evaluating various population effect sizes and comparing the observed 
  group means against these benchmarks.
  The effect size is Cohen's f of the group means \eqn{\mu_j} (with group
  sizes \eqn{N_j}): \eqn{f^2 = SSB / SSW = \sum_j N_j (\mu_j - \bar{\mu})^2 / \sum_j N_j s^2_j},
  with \eqn{\bar{\mu}} the weighted mean of the group means and \eqn{s^2_j}
  the residual variance (based on N, so that SSW is the residual sum of
  squares). If the GORIC(A) object is based on estimates and their covariance
  matrix, \eqn{s^2_j = N_j} times the variance of the j-th estimate, so then
  the covariance matrix should be based on N (not N - k).
  The population group means are \eqn{\bar{\mu} + c (\mu_j - \bar{\mu})}
  (or based on \code{ratio_pop_means}), with \eqn{c} such that their effect
  size equals \code{pop_es}. If the GORIC(A) object is based on an \code{lm}
  object, the group means are derived from the parameters via
  \code{\link{group_means}}, such that the model may contain an intercept,
  multiple factors (with or without interaction) and/or covariates (the group
  means are then adjusted means at 0 for the covariates). The data are then
  simulated for the parameters of the model (the ones in the hypotheses),
  where the population values of the parameters are such that they imply the
  population group means. If the GORIC(A) object is based on estimates and
  their covariance matrix, the estimates should be the group means (so, no
  intercept) and \code{group_size} should be specified.
  \code{benchmark_asymp} benchmarks the population estimates of a given 
  GORIC(A) object by evaluating various population estimates and comparing them 
  against the observed estimates.
  
  \code{print.benchmark} prints the results of benchmark analyses performed on 
  objects of class \code{benchmark}.
  
  \code{plot.benchmark} generates density plots for benchmark analyses of objects 
  of class \code{benchmark}. 
  
  The benchmark function leverages the \pkg{future} package for parallel processing, 
  allowing users to speed up computations by distributing tasks across multiple 
  cores or machines. If the user does not specify a parallelization plan using 
  \code{future::plan()}, the package will choose an appropriate strategy based 
  on the user's operating system. Specifically, on Windows, the package defaults 
  to using \code{multisession}, which creates separate R sessions for each 
  parallel task. On Unix-like systems (such as Linux and macOS), the package 
  defaults to \code{multicore}, which uses forked R processes to avoid the 
  overhead of setting up separate R sessions.
  
  The \code{plan()} must be specified before running the benchmark function, e.g.,
  \code{future::plan(future::multisession, workers = ncpus)}
  }


\section{Note on the hypothesis rate}{
  The hypothesis rate (\code{hypothesis_rate}) denotes, per comparison (row), how often
  the preferred hypothesis is supported (ratio of GORIC(A) weights > threshold) if this
  population were true, with this sample size; that is, the power to support the
  preferred hypothesis under this population. Likewise, the rlw rate
  (\code{rate_rlw}) denotes how often the ratio of log-likelihood weights exceeds its
  threshold under this population. For the \code{"No-effect"} population, the rates
  are not printed, but they are stored in the object.

  For the \code{"Observed"} population, the population values are the sample
  estimates. Since these vary from sample to sample, so does its rate. When the rate is
  close to 0 or 1, this hardly matters; when it is in between, the rate in another
  sample could be quite different (Preacher & Merkle, 2012).

  As a sensitivity check for the \code{"Observed"} population, the rates when the draws
  with a ratio of weights between threshold / 1.1 and threshold \eqn{\times}{x} 1.1 (near-ties) are
  discarded can be found in the output elements \code{hypothesis_rate_band} and
  \code{rate_rlw_band}: the number of draws with a ratio above threshold \eqn{\times}{x} 1.1 divided
  by the number of draws outside that band. Notably, this is not a corrected rate, but
  meant as a check: a clear difference with the rate based on all draws shows that the
  rate is driven by near-ties.

  The percentiles of the benchmark distributions are not affected by this, also for the
  \code{"Observed"} population: Preacher and Merkle found that percentile intervals of
  bootstrap distributions of (differences in) IC values perform well. For confidence
  intervals of the GORIC(A) values, their differences, weights, and ratios of weights,
  see \code{\link{goric_CI}}.
}


\references{
  Preacher, K. J., & Merkle, E. C. (2012). The problem of model selection uncertainty in 
  structural equation modeling. \emph{Psychological Methods, 17}(1), 1--14. 
  \doi{10.1037/a0026804}
}

\value{
  \code{benchmark_means} and \code{benchmark_asymp} return a list of 
  class \code{benchmark_means}, \code{benchmark}, and \code{list} or 
  \code{benchmark_asymp}, \code{benchmark}, and \code{list} containing the 
  results of the benchmark analysis.
  Among others, it contains \code{hypothesis_rate} and \code{rate_rlw} (per
  population, the rates at the thresholds \code{hypo_rate_threshold} and
  \code{threshold_rlw}), and \code{hypothesis_rate_band} and \code{rate_rlw_band}
  (per population, a matrix with, per comparison, \code{nr_draws_above}: the number
  of draws with a ratio of weights above threshold \eqn{\times}{x} 1.1; \code{nr_draws_outside}:
  the number of draws outside [threshold / 1.1, threshold \eqn{\times}{x} 1.1]; and \code{rate} =
  \code{nr_draws_above}/\code{nr_draws_outside}); see the note on the hypothesis rate.
  
  \code{print.benchmark} prints formatted benchmark analysis results to the
  console and invisibly returns the object, with the rates (including
  \code{hypothesis_rate_band} and \code{rate_rlw_band}) recomputed at the
  thresholds used in that print call.
  
  \code{plot.benchmark} returns a gtable object that can be displayed or further 
  customized using various functions from the gridExtra and grid packages. This 
  allows for flexible and detailed adjustments to the appearance and layout of the plot.
}

\author{Leonard Vanbrabant and Rebecca Kuiper}


\examples{
set.seed(1234)
# Generate data for 4 groups with different group sizes
group1 <- rnorm(10, mean = 5, sd = 0.1)
group2 <- rnorm(20, mean = 5.5, sd = 1)
group3 <- rnorm(30, mean = 6, sd = 0.5)
group4 <- rnorm(40, mean = 6.5, sd = 0.8)

# Combine data into a data frame
data <- data.frame(
  value = c(group1, group2, group3, group4),
  group = factor(rep(1:4, times = c(10,20,30,40)))
)

# Perform ANOVA
anova_result <- aov(value ~ -1 + group, data = data)

# model/hypothesis
h1 <- 'group1 < group2 < group3 < group4'
h2 <- 'group1 > group2 < group3 < group4'

# fit h1 and h2 model against the unconstrained model (i.e., failsafe to avoid
# selecting a weak hypothesis)
fit_goric <- goric(anova_result, hypotheses = list(H1 = h1, H2 = h2), 
                   comparison = "unconstrained", type = "goric")

# by default: ES = 0 \& ES = observed ES
# In practice you want more iterations: by default (iter = NULL), the number of
# iterations is set automatically (from 500 up to at most 2000).
\donttest{
# multisession supports windows machines
# future::plan(future::multisession, workers = 8)
# multicore supports unix machines
# future::plan(future::multicore, workers = 8)
benchmark_results_mean <- benchmark(fit_goric, iter = 10, model_type = "means")
print(benchmark_results_mean)
}

# by default the ratio of GORIC weights for the preferred hypothesis (here h1) is
# plotted against its competitors (i.e., h2 and the unconstrained). The
# readability of the plot can be improved by setting the x_lim option.
\donttest{plot(benchmark_results_mean, output_type = "rgw")}

# specify custom effect-sizes
\donttest{
benchmark_results_mean_es <- benchmark(fit_goric, iter = 10,
                                       pop_es = c(0, 0.1),
                                       model_type = "means")
print(benchmark_results_mean_es)  
}  

# Benchmark asymptotic estimates
\donttest{
fit_gorica <- goric(anova_result, hypotheses = list(h1=h1), 
                    comparison = "complement", type = "gorica")

# by default: no-effect \& estimates from the sample are used
benchmark_results_asymp <- benchmark(fit_gorica, sample_size = 30, iter = 5,
                                     model_type = "asymp")
print(benchmark_results_asymp)
}

\donttest{
# specify custom population estimates
my_pop_est <- rbind("no" = c(0,0,0,0), "observed"= coef(anova_result))

benchmark_results_asymp <- benchmark(fit_gorica, sample_size = 30, 
                                     iter = 5, pop_est = my_pop_est,
                                     model_type = "asymp")
print(benchmark_results_asymp)
plot(benchmark_results_asymp, x_lim = c(0, 75))
}
}
