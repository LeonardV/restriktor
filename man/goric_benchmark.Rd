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
                  quant = NULL, iter = 1000, 
                  control = list(convergence_crit = 1e-03, 
                                 chunk_size = 1e4), 
                  ncpus = 1, cl = NULL, seed = NULL, \ldots)

  benchmark_asymp(object, pop_est = NULL, sample_size = NULL, 
                  alt_sample_size = NULL, quant = NULL, iter = 1000, 
                  control = list(convergence_crit = 1e-03, 
                                 chunk_size = 1e4), 
                  ncpus = 1, cl = NULL, seed = NULL, \ldots)
  
  \method{print}{benchmark}(x, output_type = c("rgw", "gw", "rlw", "ld", "all"), 
                            color = TRUE, ...)
  
  \method{plot}{benchmark}(x, output_type = c("rgw", "rlw", "gw", "ld"), 
                           percentiles = c(0.05, 0.95), x_lim = c(), 
                           alpha = 0.50, nrow_grid = NULL, ncol_grid = 1, 
                           distr_grid = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{con_goric} (a GORIC(A) object from the \code{goric} function).}
  \item{model_type}{If "means" (default), the model parameters relect (adjusted) means, else
  model_type = "asymp". See details for more information about \code{asymp}.}
  \item{x}{An object of class \code{benchmark} or \code{benchmark}.}
  \item{pop_es}{A scalar or a vector of population Cohen's f (effect-size) values. 
  By default, it benchmarks ES = 0 and the observed Cohen's f.}
  \item{pop_est}{A 1 x k vector or an n x k matrix of population 
  estimates to benchmark. By default, all estimates are set to zero and the 
  observed estimates from the sample are used.}
  \item{ratio_pop_means}{A 1 x k vector denoting the relative difference
  between the k group means. Note that a ratio of \code{c(3, 2, 1)} gives the 
  same as \code{c(1, 0, -1)}, as the consecutive relative differences are 1 in 
  both ratios. By default, the relative differences from the data are used.}
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
  \item{quant}{Quantiles for benchmarking results. Defaults 2.5\%, 5\%, 35\%, 50\%, 65\%, 95\%, and 97.5\%}
  \item{iter}{The number of iterations for benchmarking. Defaults to \code{1000}.}
  \item{control}{A list of control parameters including \code{convergence_crit} 
  and \code{chunk_size}. For more information, see \link{goric}.}
  \item{ncpus}{Number of CPUs to use for parallel processing. Defaults to \code{1}.}
  \item{cl}{A cluster object for parallel processing. If not supplied, a cluster 
  on the local machine is created for the duration of the \code{benchmark} call.}
  \item{seed}{A seed for random number generation.}
  \item{output_type}{A character vector specifying the type of output to print 
  or plot. Options are \code{"all"}, \code{"gw"} (goric(a) weights), 
  \code{"rgw"} (ratio of goric(a) weights), \code{"rlw"} (ratio of likelihood values), 
  and \code{"ld"} (likelihood difference). Defaults to \code{"all"} for print and \code{"gw"} for plot.}
  \item{color}{If TRUE, the output will include ANSI color coding. Set \code{color = FALSE} 
  when using this function in R Markdown documents to avoid rendering issues with color codes.}
  \item{alpha}{Alpha refers to the opacity of a geom. Values of alpha range from 0 to 1, with lower values corresponding to more transparent colors.}
  \item{nrow_grid}{An integer value representing the number of rows in the grid layout.}
  \item{ncol_grid}{An integer value representing the number of columns in the grid layout.}
  \item{distr_grid}{If TRUE, the facet_grid function is used to create a grid of separate plots for each effect-size (estimates).}
  \item{percentiles}{A numeric vector specifying the lower and upper percentiles. Defaults to \code{0.025, 0.975}.}
  \item{x_lim}{A numeric vector of length 2 specifying the x-axis limits. Defaults to \code{c()}.}
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
  }

\value{
  \code{benchmark_means} and \code{benchmark_asymp} return a list of 
  class \code{benchmark_means}, \code{benchmark}, and \code{list} or 
  \code{benchmark_asymp}, \code{benchmark}, and \code{list} containing the 
  results of the benchmark analysis.
  
  \code{print.benchmark} does not return a value. It prints formatted benchmark 
  analysis results to the console.
  
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
# In practice you want to increase the number of iterations (default = 1000).
\donttest{
benchmark_results_mean <- benchmark(fit_goric, iter = 10, model_type = "means")
print(benchmark_results_mean)
}

# by default the ratio of GORIC weights for the preferred hypothesis (here h1) is
# plotted against its competitors (i.e., h2 and the unconstrained). To improve
# the readability of the plot, the argument hypothesis_comparison can be used to
# focus on a specif competitor. Further readability can be achieved by setting
# the x_lim option. 
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
