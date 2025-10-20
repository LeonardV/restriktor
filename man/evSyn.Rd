\name{evSyn}
\alias{evSyn}
\alias{evsyn}
\alias{evSyn_est}
\alias{evSyn_LL}
\alias{evSyn_ICvalues}
\alias{evSyn_ICweights}
\alias{evSyn_ICratios}
\alias{evSyn_gorica}
\alias{evSyn_escalc}
\alias{print.evSyn}
\alias{print.summary.evSyn}
\alias{summary.evSyn}
\alias{plot.evSyn}

\title{GORIC(A) Evidence synthesis}
\description{GORIC(A) evidence synthesis aggregates the evidence for theory-based 
hypotheses from multiple studies that may use diverse designs to investigate the 
same central theory.}

\usage{
# I try to guess the input type.
evSyn(object, input_type = NULL, \ldots)

# input: Parameter estimates and covariance matrix
evSyn_est(object, \ldots, VCOV = list(), hypotheses = list(),
          type_ev = c("added", "equal", "average"), 
          comparison = c("unconstrained", "complement", "none"),
          hypo_names = c(), type = c("gorica", "goricac"), 
          order_studies = c("input_order", "ascending", "descending"),
          study_names = c(), study_sample_nobs = NULL)

# input: Log-likelihood and penalty values
evSyn_LL(object, \ldots, PT = list(), type_ev = c("added", "equal", "average"),
         hypo_names = c(), 
         order_studies = c("input_order", "ascending", "descending"),
         study_names = c())

# input: GORIC(A), ORIC, AIC values
evSyn_ICvalues(object, \ldots, type_ev = c("added", "average"), hypo_names = c(), 
               order_studies = c("input_order", "ascending", "descending"), 
               study_names = c())

# input: AIC or ORIC or GORIC or GORICA weights or (Bayesian) posterior 
# model probabilities
evSyn_ICweights(object, \ldots, priorWeights = NULL, hypo_names = c(),
                order_studies = c("input_order", "ascending", "descending"), 
                study_names = c())

# input: Ratio of AIC or ORIC or GORIC or GORICA weights or 
# (Bayesian) posterior model probabilities
evSyn_ICratios(object, \ldots, priorWeights = NULL, hypo_names = c(),
               order_studies = c("input_order", "ascending", "descending"), 
               study_names = c())

# input: result from the goric() function
# Note that this is a wrapper function for evSyn_LL.
evSyn_gorica(object, \dots, type_ev = c("added", "equal", "average"), 
             hypo_names = c(), 
             order_studies = c("input_order", "ascending", "descending"), 
             study_names = c())

# input: Result from the escalc() function from the metafor package. 
# Note that this is a wrapper function for evSyn_est.
evSyn_escalc(data, yi_col = "yi", vi_cols = "vi", cluster_col, outcome_col, \ldots)

\method{print}{evSyn}(x, digits = max(3, getOption("digits") - 4), \dots)

\method{summary}{evSyn}(object, \dots)

\method{print}{summary.evSyn}(x, digits = max(3, getOption("digits") - 4), \dots)  

\method{plot}{evSyn}(x, output_type = "gorica_weights", xlab = NULL, \dots)
}

\arguments{
  \item{object}{Currently, the following objects can be processed:
    \itemize{
      \item a list of vectors with (standardized) parameter estimates (the \code{VCOV} argument is required);   
      \item a list of vectors with log-likelihood values (the \code{PT} argument is required);
      \item a list of vectors with GORIC(A) weights;
      \item a list of vectors with ratio of GORIC(A) weights;
      \item a list of vectors with GORIC(A) values;
      \item a list of objects of class goric;
      \item a data.frame of class escalc from the \pkg{metafor} package;
    }
  }
  
  \item{input_type}{
  \code{character} Specifies the type of input provided to the function. Valid options are:
  \describe{
    \item{\code{"est_vcov"}}{Indicates that the input consists of estimates and covariance matrices. Invokes the \code{evSyn_est()} function.}
    \item{\code{"ll_pt"}}{Indicates that the input consists of log-likelihoods and penalty values. Invokes the \code{evSyn_LL()} function.}
    \item{\code{"icweights"}}{Indicates that the input consists of information criterion (IC) weights. Invokes the \code{evSyn_ICweights()} function.}
    \item{\code{"icratios"}}{Indicates that the input consists of IC weights ratios. Invokes the \code{evSyn_ICratios()} function.}
    \item{\code{"icvalues"}}{Indicates that the input consists of IC values. Invokes the \code{evSyn_ICvalues()} function.}
    \item{\code{"gorica"}}{Indicates that the input is of class goric from the \code{goric} function. Invokes the \code{evSyn_gorica()} function.}
    \item{\code{"escalc"}}{Indicates that the input is of class escalc from the \pkg{metafor} package. Invokes the \code{evSyn_escalc()} function.}
  }
  If \code{input_type} is \code{NULL}, the function attempts to infer the input type based on the structure of the \code{object} and other arguments.
}
  
  \item{VCOV}{a list of covariance matrices of the (standardized) parameter estimates 
  of interest.} 
  
  \item{PT}{a list of vectors with penalty values.}
  
  \item{x}{an object of class \code{evSyn}}

  \item{type_ev}{type of evidence-synthesis approach: Equal-evidence approach 
  (\code{type_ev = "equal"}), Added-evidence approach (\code{type_ev = "added"}), or
  Average-evidence approach (\code{type_ev = "average"}). See details for more 
  information.}

  \item{type}{If \code{"gorica"}, a proxi for the log-likihood is calculated by assuming that 
  the parameter estimates (of interest) are multivariate normally distributed. 
  Then, either the unbiased covariance matrix of the estimates (based on 'N'; e.g., for an lm object) is used
  or the biased one obtained via the vcov function (e.g., for a lavaan object). 
  In the latter case, a corresponding message will be printed. If \code{"goricac"}, 
  then the small-sample-size correction is applied.}

  \item{hypotheses}{When applying the same set of hypotheses to each study, 
  the syntax structure should be as follows: "hypotheses = list(H1, H2, \ldots)". 
  However, if a different set of hypotheses is applied to each study, the syntax 
  structure should be as follows: hypotheses = list(set1 = list(H11, H12), set2 = list(H21, H22)).
  See \code{\link{goric}} how to specify the hypotheses syntax or see the example 
  section below.}

  \item{comparison}{denotes the failsafe hypothesis; one can choose from "unconstrained", "complement", and "none".
  If "\code{none}", the model is only compared against the models provided by the user.
  By default, if one user-specified hypothesis, "\code{comparison = 'complement'}"; 
  if multiple user-specified hypotheses, "\code{comparison = 'unconstrained'}".
  Note that the complement can only be computed for one model/hypothesis at a time.} 
  
  \item{priorWeights}{vector that represents the prior belief for this model. By default,
  equal prior weights are used (i.e., 1/(#hypotheses)). Notably, in case the prior
  weights do not sum to 1, it will be rescaled such that it does; which implies that
  relative importance can be used and not per se weights.}
  
  \item{hypo_names}{character vector for labelling the hypotheses. By default the
  names are set to H1, H2, \ldots}

  \item{output_type}{a plot options. If "\code{gorica_weights}" (default), 
  a plot with the cumulative goric(a)
  weights and goric(a) weights per study is displayed. If "\code{ll_weights}",
  a plot with the cumulative log-likelihood weights and log-likelihood weights 
  per study is displayed.}

  \item{study_names}{
  A character vector specifying the names or labels of the individual studies.
  When provided, these names will be used throughout the output, including in 
  tables, summaries, and plots (in place of numeric study indices).
  If \code{NULL} (default), studies are automatically numbered 
  (1, 2, ..., K) following the order determined by \code{order_studies}.}

  \item{study_sample_nobs}{
  A numeric vector specifying the sample size (\code{n}) for each study.
  When provided, this information will be included in summaries and can be used 
  for weighted averaging or diagnostic purposes in the evidence synthesis. 
  If \code{NULL} (default), all studies are assumed to have equal weight 
  with respect to sample size.}
  
  \item{order_studies}{
  A character string that determines how studies are ordered in the output and plots. 
  If \code{"input_order"} (default), studies are displayed in the order they were 
  provided. If \code{"ascending"} or \code{"descending"}, studies are ordered based 
  on the ascending or descending cumulative support for the preferred hypothesis, 
  respectively.}
  
  \item{xlab}{a vector specifying custom labels for the x-axis (be aware that 
  study ordering may be changed based on \code{order_studies}). The length of 
  the vector must match the number of studies in the dataset. If not provided, 
  the x-axis labels default to a sequence from 1 to the number of studies 
  (possibly re-ordered based on \code{order_studies}).}

  \item{data}{an object of class "escalc" from the \pkg{metafor} package.}

  \item{cluster_col}{a character string specifying the column in \code{data} that 
  contains the cluster/study identifier. Observations within the same cluster/study should 
  share the same value in this column. 
  When not specified, the function looks for a column with one of the following names:
  "trial", "study", "author", "authors", "Trial", "Study", "Author", or "Authors".
  }
  
  \item{outcome_col}{a column name in \code{data} containing the outcome 
  identifiers/labels for each observation within a cluster/study. 
  The hypothesis/-es should be based on this labelling.
  By default (i.e., if 'outcome_col' is null), the function assumes the label is 
  "theta" 
  % ms juist via deze de labels voor de hypo opgeven!
  (one outcome variable); 
  when there is additionally a column called 'outcome', it will use that as the outcome identifier,
  employing its unique levels/values as the labels (more outcome variables).
  If 'outcome_col' is not null but user-specified, then the labeling should be based
  on the unique level(s)/value(s) in the corresponding data column.
  Notably, outcome_col is related to the labeling; 
  the column with estimates is referred to as 'yi_col' and is, 
  by default, set to "yi".
  For more information, see Example 5 below.
  }
  
  \item{yi_col}{a character string specifying the column in \code{data} that 
  contains the outcome values / estimates for each observation. 
  The default is \code{"yi"}.
  For more information, see Example 5 below.
  }
  
  \item{vi_cols}{a character vector specifying the columns in \code{data} that 
  contain the variance and covariance values for each observation. The default is, 
  in the case of one outcome variable, \code{"vi"};
  in the case of multiple (namely, O) variables, \code{c("v1i", \ldots, "vOi"}.
  For more information, see Example 5 below.
  }
  
  \item{digits}{the number of significant digits to use when printing.}

  \item{\ldots}{This depends on the class of the object.}
}

\details{
In the \strong{added-evidence} approach, evidence from each study or dataset is 
cumulatively aggregated. This means that for every new study, the log-likelihood 
and the penalty term are added to the cumulative totals. The strength of the 
aggregated evidence in this approach depends on the nature of the evidence itself. 
Simply having more studies doesn't necessarily mean stronger evidence if those 
studies provide weak or contradictory evidence.\cr 
Opt for this approach when you anticipate each new piece of evidence to provide an 
incremental contribution to the overall evidence, without the need to normalize or 
average across datasets. It's especially suitable when you believe that the 
aggregated evidence from multiple studies is stronger than if the data were combined 
into a single study.

The \strong{equal-evidence} approach aggregates the cumulative evidence in the same 
manner as the added-evidence approach. However, when calculating the GORICA, the 
cumulative evidence is divided by the number of studies. This ensures that the 
contribution from each study or dataset remains equal, regardless of the total count. 
Conceptually, aggregating evidence from multiple studies in this approach can be 
likened to obtaining evidence from a single larger study, similar to how a 
meta-analysis treats combined evidence.\cr
Choose this method when you want each study to contribute equally to the overall 
evidence, irrespective of the size or scope of each individual dataset. It's ideal 
for situations where you view the combined evidence from multiple studies as 
equivalent to that from a single, larger study.

The \strong{average-evidence} method can be conceptualized as a form of multiverse 
analysis. When faced with a single dataset, there are often numerous analytical 
choices available, such as handling missing data, selecting variables, or choosing 
statistical methods. Each choice can lead to a different analysis or model, creating 
a "multiverse" of possible outcomes.\cr
For each of these analyses, an "evidence" score can be calculated, indicating how 
well the model fits the data. Some models might offer a superior fit, while others 
might not align as closely with the data. The average-evidence method aggregates 
these scores, providing an average measure of fit across all considered models.
This approach offers an overarching perspective on the general trend across all 
analyses. If the average evidence suggests a good fit, it indicates that the 
majority of the chosen analyses align well with the data. This method is invaluable 
for assessing the robustness of results, ensuring that findings are not merely 
artifacts of a specific analytical choice but are consistent across various model 
specifications on the same dataset.\cr
Opt for the average-evidence approach when you wish to gauge the central tendency 
of evidence across multiple analytical choices. It's especially beneficial when 
aiming to determine the robustness of results across various model specifications 
applied to the same dataset.
}


\value{An object of class evSyn for which a print, summary and plot function is available. 
The output comprises, among other things, the cumulative and final evidence 
for the theory-based hypotheses.}


\author{Leonard Vanbrabant and Rebecca M. Kuiper}

\examples{
## By following these examples, you can appropriately specify hypotheses based on 
## your research questions and analytical framework.

# The hypotheses (i.e., constraints) have to be in a list. It is recommended to name
# each hypothesis in the list. Otherwise, the hypotheses are named accordingly 'H1', 'H2', \ldots.

# Example using text-based syntax (the labels x1, x2, and x3 are the names of 
# coef(model) or names(vector))
H1 <- '(x1, x2, x3) > 0'
H2 <- '(x1, x3) > 0; x2 < 0'
H3 <- 'x1 > 0; x2 < 0; x3 = 0'
hypotheses = list(hypo1 = H1, hypo2 = H2, hypo3 = H3)
# Here, it is assumed that each study uses the same hypothesis specification

# Example using mixed syntax: 
H1 <- 'x1 > x2 > x3 > 0'
hypotheses = list(Ha = H1, Hb = 'x1 = x2 > x3')
# Here, it is assumed that each study uses the same hypothesis specification

# Example using a different set of hypotheses for each study: 
# Set for Study 1
H11 <- 'group1 = group2 > group3' 
H12 <- 'group2 > group1 > group3'   
# Set for Study 2
H21 <- 'gr1 = gr2 > gr3'
H22 <- 'gr2 > gr1 > gr3'
#
# correct
hypotheses = list(set1 = list(H11, H12), set2 = list(H21, H22))
# note that the list names set1 and set2 are redundant and can be left out. 
# It is crucial to ensure that the hypotheses across each set are ordered in a similar manner.
#
# NOT correct
#hypotheses = list(set1 = list(H12, H11), set2 = list(H21, H22))


## Example 1: based on one estimate and its variance
# 4 studies

# estimates (beta) (from the 4 primary studies)
est_1 <- c(beta1 = 0.09)
est_2 <- c(beta1 = 0.14)
est_3 <- c(beta1 = 1.09)
est_4 <- c(beta1 = 1.781)
Param_studies <- list(est_1, est_2, est_3, est_4)

# standard error of the beta's (from the 4 primary studies)
vcov_est_1 <- matrix(c(0.029^2), nrow = 1)
vcov_est_2 <- matrix(c(0.054^2), nrow = 1)
vcov_est_3 <- matrix(c(0.093^2), nrow = 1)
vcov_est_4 <- matrix(c(0.179^2), nrow = 1)
CovMx_studies <- list(vcov_est_1, vcov_est_2, vcov_est_3, vcov_est_4)

# Set of hypotheses for each study
# Note: in this case the same for each study
H0   <- "beta1 = 0"
Hpos <- "beta1 > 0"
Hneg <- "beta1 < 0"
hypotheses <- list(H0 = H0, Hpos = Hpos, Hneg = Hneg)

# Since this covers the whole space / covers all theories, we do not need a safeguard-hypothesis:
comparison <- "none"

evS4_added <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                    hypotheses = hypotheses,
                    type_ev = "added", 
                    comparison = "none")
evS4_added
summary(evS4_added)
plot(evS4_added)

evS4_equal <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                    hypotheses = hypotheses,
                    type_ev = "equal", 
                    comparison = "none")
evS4_equal
summary(evS4_equal)
plot(evS4_equal)


## Example 2: based on (multiple) estimates and their covariance matrices
# 2 studies

# estimates and their covariance matrices (from the 2 primary studies)
est_1 <- c(1.88, 2.54, 0.02)
names(est_1) <- c("group1", "group2", "group3")
vcov_est_1 <- diag(c(0.2149074, 0.2149074, 0.1408014))

est_2 <- c(0.98, 0.02, 0.27)
names(est_2) <- c("gr1", "gr2", "gr3") 
vcov_est_2 <- diag(c(0.1382856, 0.1024337, 0.0987754))

# list of estimates 
object <- list(est_1, est_2)
# list of covariance matrices of the estimates
VCOV <- CovMx_studies <- list(vcov_est_1, vcov_est_2)

# Hypotheses (study-specific)
# names(est_1) # Specify restrictions using those names
H11 <- 'group1 = group2 > group3'
H12 <- 'group2 > group1 > group3'
# names(est_2) # Specify restrictions using those names
H21 <- 'gr1 = gr2 > gr3'
H22 <- 'gr2 > gr1 > gr3' 

# list of (competing) hypotheses
hypotheses <- list(H1 = list(H11, H12), H2 = list(H21, H22))

# Evidence synthesis (added-evidence)
evS2_added <- evSyn(object, VCOV = VCOV, hypotheses = hypotheses,
                    type_ev = "added", comparison = "unconstrained") 
evS2_added
plot(evS2_added)


## Example 3: based on (multiple) estimates and their covariance matrices
# 3 studies

# generate data
ratio <- c(1,1.1,1.2)
n <- c(30, 50, 100)

# Generate data1
n1 <- n[1]
x11 <- rnorm(n1)
x12 <- rnorm(n1)
x13 <- rnorm(n1)
data <- cbind(x11, x12, x13)
# Standardize data - since parameters for continuous variables will be compared
data1 <- as.data.frame(scale(data))
y1 <- ratio[1]*data1$x11 + ratio[2]*data1$x12 + ratio[3]*data1$x13 + rnorm(n1)
# Note: since there is one outcome, the outcome does not need to be standardized.
fit.lm1 <- lm(y1 ~ 1 + x11 + x12 + x13, data = data1)

# Generate data2
n2 <- n[2]
x21 <- rnorm(n2)
x22 <- rnorm(n2)
x23 <- rnorm(n2)
data <- cbind(x21, x22, x23)
data2 <- as.data.frame(scale(data))
y2 <- ratio[1]*data2$x21 + ratio[2]*data2$x22 + ratio[3]*data2$x23 + rnorm(n2)
fit.lm2 <- lm(y2 ~ 1 + x21 + x22 + x23, data = data2)

# Generate data3
n3 <- n[3]
x31 <- rnorm(n3)
x32 <- rnorm(n3)
x33 <- rnorm(n3)
data <- cbind(x31, x32, x33)
data3 <- as.data.frame(scale(data))
y3 <- ratio[1]*data3$x31 + ratio[2]*data3$x32 + ratio[3]*data3$x33 + rnorm(n3)
fit.lm3 <- lm(y3 ~ 1 + x31 + x32 + x33, data = data3)

# Extract estimates and their covariance matrix (per study)
est_1 <- coef(fit.lm1)
est_2 <- coef(fit.lm2)
est_3 <- coef(fit.lm3)
vcov_est_1 <- vcov(fit.lm1)
vcov_est_2 <- vcov(fit.lm2)
vcov_est_3 <- vcov(fit.lm3)

names(est_1) <- names(est_2) <- names(est_3) <- c("intercept", "x1", "x2", "x3")

# Parameter estimate values from the primary studies
Param_studies <- list(est_1, est_2, est_3)

# standard error of the beta's
CovMx_studies <- list(vcov_est_1, vcov_est_2, vcov_est_3)

# Set of hypotheses for each study. Note: in this case the same for each study
hypothesis <- 'x1 < x2 < x3'  

# In our theory, we compare estimates of continuous variables, so we standardized 
# the data beforehand to ensure comparability. In 'Param_studies' and 'CovMx_studies', 
# the intercept can be omitted without affecting the GORIC(A) weights, as there are 
# no restrictions on it. Since we have only one theory-based hypothesis, we will 
# utilize the more powerful complement of the hypothesis (Vanbrabant, Van Loey, Kuiper, 2019). 
# The complement represents the remaining 11 theories, while the unconstrained 
# scenario includes all 12 possible theories, including H1.

# Evidence synthesis (added-evidence)
evS3 <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
              hypotheses = list(H1 = hypothesis),
              type_ev = "added", 
              comparison = "complement") 
evS3
plot(evS3)


## Example 4: based on loglikelihood values and penalty values

# Make LL and PT values a list
LL <- as.list(data.frame(t(myLLs)))
penalty.values <- as.list(data.frame(t(myPTs)))

# Evidence synthesis
#
# Added-evidence
evS_LL_added <- evSyn(object = LL, PT = penalty.values, 
                      hypo_names = colnames(myLLs),
                      type_ev = "added")
evS_LL_added
#
# Equal-evidence
evS_LL_equal <- evSyn(object = LL, PT = penalty.values, 
                      hypo_names = colnames(myLLs),
                      type_ev = "equal")
evS_LL_equal


## Example 5: based on AIC, ORIC, GORIC(A) values

# Make GORIC values a list
goric.values <- as.list(data.frame(t(myGORICs)))

# Evidence synthesis (added-evidence)
evS_Gv <- evSyn(goric.values, 
                hypo_names = colnames(myGORICs))
evS_Gv


## Example 6: based on an escalc object

## Example 6a: one outcome, default labeling
if (requireNamespace("metadat", quietly = TRUE)) {
  # access data within metafor
  data_a <- metadat::dat.berkey1998[metadat::dat.berkey1998$outcome == "PD", ][, c(1:4, 6:7)]

  # yi gives the estimates and vi their variances.
  # So, here:
  # cluster_col = "trial" using the first variable with a cluster variable 
  # name (year and author are other possibilities)
  # outcome_col = NULL # No 'outcome' variable, so 'theta' is used as labeling
  # yi_col = "yi"
  # vi_cols = "vi"
  #  

  # The hypothesis should thus now be based on the default labelling theta:
  H1 <- "abs(theta) > 0.5" 
  evS_escalc_a <- evSyn(
    data_a, 
    hypotheses = list(H1 = H1), 
    type_ev = "equal",
    comparison = "complement"
  )
  evS_escalc_a
}


## Example 6b: one outcome, default labeling coming from 'outcome' variable
if (requireNamespace("metadat", quietly = TRUE)) {
  data_b <- metadat::dat.berkey1998[metadat::dat.berkey1998$outcome == "PD", ][, c(1:7)]
  # yi gives the estimates and vi their variances.
  # So, here:
  # cluster_col = "trial" # using the first variable with a cluster variable 
  # name (year and author are other possibilities)
  # outcome_col = "outcome" # was NULL, but found 'outcome' variable, so use that 
  # labeling, that is, "PD"
  # yi_col = "yi"
  # vi_cols = "vi"
  #  
  # The hypothesis should thus now be based on the labeling based on the 'outcome' variable:
  H1 <- "abs(PD) > 0.5" 
  evS_escalc_b <- evSyn(data_b, 
                        hypotheses = list(H1 = H1), 
                        type_ev = "equal",
                        comparison = "complement")
  evS_escalc_b
}


## Example 6c: 2 outcomes, default labeling coming from 'outcome' variable
if (requireNamespace("metadat", quietly = TRUE)) {
  data_c <- metadat::dat.berkey1998
  # yi gives the estimates (for both outcomes), vi their variances,
  # and v1i & v2i combined their covariance matrices.
  # So, here:
  # cluster_col = "trial" # using the first variable with a cluster variable 
  # name (year and author are other possibilities)
  # outcome_col = "outcome" # was NULL, but found 'outcome' variable, so use 
  # that labeling, that is, "PD"
  # yi_col = "yi"
  # vi_cols = c("v1i", "v2i")
  #  
  # The hypothesis should thus now be based on the labeling based on the 'outcome' variable:
  H1 <- "abs(AL) > 0.5, abs(PD) > 0.5" 
  evS_escalc_c <- evSyn(data_c, 
                        hypotheses = list(H1 = H1), 
                        type_ev = "equal",
                        comparison = "complement")
  evS_escalc_c
}

}


