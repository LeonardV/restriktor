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
                  hypo_rate_threshold = 1, threshold_rlw = 1, \dots)

  benchmark_asymp(object, pop_est = NULL, sample_size = NULL, 
                  alt_sample_size = NULL, quant = NULL, iter = NULL, 
                  control = list(), 
                  ncpus = 1, seed = NULL,
                  iter_adequacy_band = c(0.495, 0.505),
                  iter_stability_tol = 1,
                  iter_min = 500, iter_step = 100, iter_max = 2000,
                  hypo_rate_threshold = 1, threshold_rlw = 1, \dots)
  
  \method{print}{benchmark}(x, output_type = c("rgw", "gw", "lw", "rlw", "ld", 
                            "rgw_log", "rlw_log", "all"), 
                            hypo_rate_threshold = NULL, threshold_rlw = NULL, 
                            color = TRUE, percentiles = NULL, \dots)
  
  \method{plot}{benchmark}(x, output_type = c("rgw", "rlw", "gw", "ld"), 
                           percentiles = NULL, x_lim = c(), log_scale = FALSE,
                           opacity = 0.50, nrow_grid = NULL, ncol_grid = 1, 
                           distr_grid = FALSE, \dots)
}
\arguments{
  \item{object}{An object of class \code{con_goric} (a GORIC(A) object from the \code{goric} function).}
  \item{model_type}{If "means", the model parameters reflect (adjusted) means, else
  model_type = "asymp" (default). See details for more information about \code{asymp}.}
  \item{x}{An object of class \code{benchmark} (i.e., of class
  \code{benchmark_means} or \code{benchmark_asymp}).}
  \item{pop_es}{A scalar or a vector of population Cohen's f (effect-size) values. 
  By default, it benchmarks ES = 0 (no-effect) and the observed Cohen's f.}
  \item{pop_est}{A 1 x k vector or an n x k matrix of population 
  estimates to benchmark. By default, all estimates are set to zero (no-effect) 
  and the observed estimates from the sample are used.}
  \item{ratio_pop_means}{A 1 x k vector denoting the relative difference
  between the k group means. Note that a ratio of \code{c(3, 2, 1)} gives the 
  same as \code{c(1, 0, -1)}, as the consecutive relative differences are 1 in 
  both ratios. By default, the relative differences from the data are used.}
  \item{group_size}{A 1 x k vector or a scalar to denote the k group sizes. If 
  a scalar is specified, it is assumed that each group is of that size. If the 
  GORIC(A) object is based on a fitted (\code{lm}) model, the group sizes are 
  by default the (cell) counts of the factor variable(s) in that model; a 
  user-specified \code{group_size} overrides these (with a message if they 
  differ). If the GORIC(A) object is based on estimates and their covariance 
  matrix (instead of on a model/fit object), \code{group_size} must be specified.}
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
  \item{iter}{The number of iterations (bootstrap draws) per population for 
  benchmarking. Defaults to \code{NULL}, in which case the number of draws is 
  determined automatically: it starts at \code{iter_min} draws and adds 
  \code{iter_step} draws at a time until the percentile of the sample value 
  within the benchmark distribution under the 'Observed' population (for both 
  the GORIC(A) weights and the log-likelihood weights) has changed by less 
  than \code{iter_stability_tol} percentage points in each of two consecutive 
  rounds, or until \code{iter_max} draws are reached. Without an 'Observed' 
  population (i.e., a user-specified \code{pop_es} or \code{pop_est}), 
  \code{iter_min} draws are used. If a number is specified, exactly that many 
  draws are used; a message is then given if the percentile has not 
  stabilized over the last 20\% of the draws or if \code{iter < iter_min}.}
  \item{iter_adequacy_band}{A vector of two increasing proportions (between 0 
  and 1) used for an informational GORICA-based check of whether the sample 
  value lies at about the median (the 50th percentile) of its own benchmark 
  distribution under the 'Observed' population. The result is stored in the 
  output (\code{median_bias_check_gw} and \code{median_bias_check_lw}) but is 
  not used to determine \code{iter}. Defaults to \code{c(0.495, 0.505)}.}
  \item{iter_stability_tol}{A non-negative number (in percentage points): the 
  maximum change in the percentile of the sample value (see \code{iter}) for 
  it to be considered stable. Defaults to \code{1}. With 
  \code{iter_stability_tol = 0}, the percentile is never considered stable, so 
  that (for \code{iter = NULL}) \code{iter_max} draws are used.}
  \item{iter_min}{A whole number \code{>= 1}: the initial number of draws if 
  \code{iter = NULL}; for a user-specified \code{iter}, a message is given if 
  \code{iter} is smaller than \code{iter_min}. Defaults to \code{500}.}
  \item{iter_step}{A whole number \code{>= 1}: the number of draws added per 
  round if \code{iter = NULL}. Defaults to \code{100}.}
  \item{iter_max}{A whole number \code{>= 1}: the maximum number of draws if 
  \code{iter = NULL}. Defaults to \code{2000}.}
  \item{hypo_rate_threshold}{A numeric value specifying the threshold for the 
  hypothesis rate: the proportion of the ratio-of-GORIC(A)-weights draws (of the 
  preferred hypothesis versus each competing hypothesis) that exceeds this 
  threshold. In \code{benchmark_means} and \code{benchmark_asymp}, it is used 
  for the \code{hypothesis_rate} element of the output (defaults to \code{1}). 
  In \code{print.benchmark}, it defaults to \code{NULL}, meaning that the 
  threshold stored in the benchmark object is used; another value can be 
  specified to recompute the printed hypothesis rate without rerunning the 
  benchmark.}
  \item{threshold_rlw}{The same as \code{hypo_rate_threshold}, but for the ratio 
  of log-likelihood weights: the proportion of its draws that exceeds this 
  threshold (the \code{rate_rlw} element of the output). Defaults to \code{1} in 
  \code{benchmark_means} and \code{benchmark_asymp}, and to \code{NULL} (the 
  threshold stored in the benchmark object) in \code{print.benchmark}.}
  \item{control}{A list of control parameters. For more information, see details \link{goric}.}
  \item{ncpus}{Number of CPUs to use for parallel processing. Defaults to \code{1}.
  See details for more information.}
  \item{seed}{A seed for random number generation.}
  \item{output_type}{A character string specifying the type of output to print 
  or plot. Options are \code{"gw"} (GORIC(A) weights), \code{"lw"} 
  (log-likelihood weights; print only), \code{"rgw"} (ratio of GORIC(A) 
  weights), \code{"rlw"} (ratio of log-likelihood weights), \code{"ld"} 
  (log-likelihood difference), \code{"rgw_log"} and \code{"rlw_log"} (log of 
  the ratio of GORIC(A) weights and of the ratio of log-likelihood weights, 
  respectively; print only), and \code{"all"} (all of these; print only). 
  Defaults to \code{"rgw"} for both print and plot.}
  \item{color}{If TRUE, the output will include ANSI color coding. Set \code{color = FALSE} 
  when using this function in R Markdown documents to avoid rendering issues with color codes.}
  \item{opacity}{The opacity of the density plots. Values range from 0 to 1, 
  with lower values corresponding to more transparent colors. (This argument 
  used to be called \code{alpha}; \code{alpha} is still accepted, with a 
  deprecation warning.)}
  \item{nrow_grid}{An integer value representing the number of rows in the grid layout.}
  \item{ncol_grid}{An integer value representing the number of columns in the grid layout.}
  \item{distr_grid}{If TRUE, the facet_grid function is used to create a grid of separate plots for each effect-size (estimates).}
  \item{percentiles}{A numeric vector specifying the percentiles (as proportions, 
  e.g. \code{c(0.05, 0.95)}) to be shown. By default (\code{NULL}), the 
  percentiles used for benchmarking (see \code{quant}) are shown. Otherwise, 
  they are recomputed from the stored draws, without rerunning the benchmark.}
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
  \code{benchmark_asymp} benchmarks the population estimates of a given 
  GORIC(A) object by evaluating various population estimates and comparing them 
  against the observed estimates.
  
  \code{print.benchmark} prints the results of benchmark analyses performed on 
  objects of class \code{benchmark}.
  
  \code{plot.benchmark} generates density plots for benchmark analyses of objects 
  of class \code{benchmark}. 
  
  If \code{object} is based on an \code{lm} model with (continuous) covariates 
  (e.g., \code{lm(y ~ -1 + group + x)}), \code{benchmark_means} treats only the 
  coefficients of the factor variable(s) as group means: Cohen's f, the group 
  sizes and the population means (\code{pop_es}) refer to these group means 
  only, while the covariate coefficients are kept at their observed estimates 
  in each population. If \code{object} is based on estimates and their 
  covariance matrix, all estimates are treated as group means.
  
  The benchmark function leverages the \pkg{future} package for parallel processing, 
  allowing users to speed up computations by distributing tasks across multiple 
  cores or machines. If \code{ncpus > 1} and no parallel plan is set (i.e., the 
  current \code{future::plan()} is sequential), the function sets one up itself 
  for the duration of the call, with \code{ncpus} workers: on Windows, 
  \code{multisession} (separate R sessions), and on Unix-like systems (such as 
  Linux and macOS), \code{multicore} (forked R processes, which avoids the 
  overhead of setting up separate R sessions). The original plan is restored 
  afterwards. A (non-sequential) plan set by the user beforehand, e.g.,
  \code{future::plan(future::multisession, workers = 4)}, is used as is.
  }


\value{
  \code{benchmark_means} and \code{benchmark_asymp} return a list of 
  class \code{benchmark_means}, \code{benchmark}, and \code{list} or 
  \code{benchmark_asymp}, \code{benchmark}, and \code{list} containing the 
  results of the benchmark analysis.
  
  \code{print.benchmark} prints formatted benchmark analysis results to the 
  console and invisibly returns (a copy of) \code{x} including the printed tables.
  
  \code{plot.benchmark} returns an object of class \code{benchmark_plot}: a list 
  with the (\pkg{ggplot2}) plots (one per competing hypothesis) and the grid 
  layout, which is drawn (using \code{gridExtra::grid.arrange}) when printed. 
  The separate plots (element \code{plots}) can be further customized.
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
# In practice you want to increase the number of iterations (by default,
# iter = NULL, it is determined automatically).
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
