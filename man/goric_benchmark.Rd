\name{benchmark_functions}
\alias{benchmark}
\alias{benchmark_means}
\alias{benchmark_asymp}
\alias{print.benchmark}
\alias{plot.benchmark}
\title{Benchmark Functions for GORIC(A) Analysis}
\description{
  The `benchmark_means` and `benchmark_asymp` functions perform benchmarking 
  for models using the Generalized Order-Restricted Information Criterion (GORIC). 
}
\usage{
  benchmark(object, model_type = NULL, \dots)
  
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
  
  \method{print}{benchmark}(x, output_type = c("all", "gw", "rgw", "rlw", "ld"), \ldots)
  
  \method{plot}{benchmark}(x, output_type = "gw", CI = 0.95, x_lim = c(), \ldots)
}
\arguments{
  \item{object}{An object of class \code{con_goric} (a GORIC object from the \code{goric} function).}
  \item{model_type}{If "means", the model parameters relect (adjusted) means, else
  model_type = "asymp".}
  \item{x}{An object of class \code{benchmark_means} or \code{benchmark_asymp}.}
  \item{pop_es}{A scalar or a vector of population Cohen's f (effect-size) values. 
  By default, it benchmarks ES = 0 and the observed Cohen's f.}
  \item{pop_est}{A 1 x k vector or an n x k matrix of population 
  estimates to benchmark. By default, all estimates are set to zero and the 
  estimates from the sample are used.}
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
  sizes from the data will be used.}
  \item{sample_size}{A scalar to denote the (total) sample sizes. Only used if 
  the GORIC object is based on estimates and their covariance matrix (instead 
  of on a model/fit object) and \code{alt_sample_size} is not \code{NULL}. }
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
  \item{CI}{A numeric value specifying the confidence interval level. Defaults to \code{0.95}.}
  \item{x_lim}{A numeric vector of length 2 specifying the x-axis limits. Defaults to \code{c()}.}
  \item{\ldots}{Not used for now.}
}
\details{
  \code{benchmark_means} benchmarks the group means of a given GORIC object 
  by evaluating various population effect sizes and comparing the observed 
  group means against these benchmarks.
  \code{benchmark_asymp} benchmarks the population estimates of a given 
  GORIC object by evaluating various population estimates and comparing them 
  against the observed estimates.
  
  \code{print.benchmark} prints the results of benchmark analyses performed on 
  objects of class \code{benchmark_means} or \code{benchmark_asymp}.
  
  \code{plot.benchmark} generates density plots for benchmark analyses of objects 
  of class \code{benchmark_means} or \code{benchmark_asymp}.
}
\value{
  \code{benchmark_means} and \code{benchmark_asymp} return a list of 
  class \code{benchmark_means}, \code{benchmark}, and \code{list} or 
  \code{benchmark_asymp}, \code{benchmark}, and \code{list} containing the 
  results of the benchmark analysis.
  
  \code{print.benchmark} does not return a value. It prints formatted benchmark 
  analysis results to the console.
  
  \code{plot.benchmark} returns a ggplot object that can be displayed or further 
  customized using ggplot2 functions.
}

\author{Leonard Vanbrabant and Rebecca Kuiper}


\examples{
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

# fit h1 model against its complement
fit_goric <- goric(anova_result, hypotheses = list(h1=h1), 
                   comparison = "complement", type = "goric")

# by default: ES = 0 \& ES = observed ES
# In practice you want to increase the number of iterations.
benchmark_results_mean <- benchmark(fit_goric, iter = 30, model_type = "means")

print(benchmark_results_mean)
plot(benchmark_results_mean, output_type = "gw")

# specify custom effect-sizes
benchmark_results_mean_es <- benchmark(fit_goric, iter = 30,
                                       pop_es = c(0, 0.1, 0.2),
                                       model_type = "means")
print(benchmark_results_mean_es)  
  
  
# Benchmark asymptotic estimates
fit_gorica <- goric(anova_result, hypotheses = list(h1=h1), 
                    comparison = "complement", type = "gorica")

# by default: no-effect \& estimates from the sample are used
benchmark_results_asymp <- benchmark(fit_gorica, sample_size = 100, iter = 30,
                                     model_type = "asymp")
print(benchmark_results_asymp)

# specify custom population estimates
my_pop_est <- rbind("no" = c(0,0,0,0), "medium" = c(4,5,6,7), "observed"= coef(anova_result))
benchmark_results_asymp <- benchmark(fit_gorica, sample_size = 100, 
                                     iter = 30, pop_est = my_pop_est,
                                     model_type = "asymp")
print(benchmark_results_asymp)
}
