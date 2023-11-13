\name{evSyn}
\alias{evSyn}
\alias{evsyn}
\alias{evSyn_est.list}
\alias{evSyn_ICweights.list}
\alias{evSyn_ICvalues.list}
\alias{evSyn_LL.list}
\alias{print.evSyn}
\alias{print.summary.evSyn}
\alias{summary.evSyn}
\alias{plot.evSyn}

\title{GORIC(A) Evidence synthesis}
\description{GORIC(A) evidence synthesis aggregates the evidence for theory-based 
hypotheses from multiple studies that may use diverse designs to investigate the 
same central theory.}

\usage{
evSyn(object, \ldots)

\method{evSyn_est}{list}(object, \ldots, VCOV = list(), hypotheses = list(),
      type = c("equal", "added", "average"), 
      comparison = c("unconstrained", "complement", "none"),
      hypo_names = c())

\method{evSyn_LL}{list}(object, \ldots, PT = list(), type = c("equal", "added", "average"),
      hypo_names = c())

\method{evSyn_ICvalues}{list}(object, \ldots, hypo_names = c())

\method{evSyn_ICweights}{list}(object, \ldots, priorWeights = NULL, hypo_names = c())

\method{print}{evSyn}(x, digits = max(3, getOption("digits") - 4), \dots)

\method{summary}{evSyn}(object, \dots)

\method{print}{summary.evSyn}(x, digits = max(3, getOption("digits") - 4), \dots)  

\method{plot}{evSyn}(x, \dots)
}

\arguments{
  \item{object}{Currently, the following objects can be processed:
    \itemize{
      \item a list of vectors with (standardized) parameter estimates (the \code{VCOV} argument is required);   
      \item a list of vectors with log-likelihood values (the \code{PT} argument is required);
      \item a list of vectors with GORIC(A) weights;
      \item a list of vectors with GORIC(A) values;
    }
  }
  \item{VCOV}{a list of covariance matrices of the (standardized) parameter estimates 
  of interest.} 
  
  \item{PT}{a list of vectors with penalty values.}
  
  \item{x}{an object of class \code{evSyn}}

  \item{type}{type of evidence-synthesis approach: Equal-evidence approach 
  (\code{type = "equal"}), Added-evidence approach (\code{type = "added"}), or
  Average-evidence approach (\code{type = "average"}). See details for more 
  information.}

  \item{hypotheses}{When applying the same set of hypotheses to each study, 
  the syntax structure should be as follows: "hypotheses = list(H1, H2, \ldots)". 
  However, if a different set of hypotheses is applied to each study, the syntax 
  structure should be as follows: hypotheses = list(set1 = list(H11, H12), set2 = list(H21, H22)).
  See \code{\link{goric}} how to specify the hypotheses syntax or see the example 
  section below.
}

  \item{comparison}{if "\code{unconstrained}" (default) the unconstrained model is 
  included in the set of models. If "\code{complement}" then the restricted object 
  is compared against its complement. Note that the complement can only be computed 
  for one model/hypothesis at a time (for now). If "\code{none}" the model is only 
  compared against the models provided by the user.} 
  
  \item{priorWeights}{vector that represents the prior belief for this model. By default,
  equal prior weights are used (i.e., 1/(#hypotheses)). Notably, in case the prior
  weights do not sum to 1, it will be rescaled such that it does; which implies that
  relative importance can be used and not per se weights.}
  
  \item{hypo_names}{character vector for labelling the hypotheses. By default the
  names are set to H1, H2, \ldots}

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


\author{Leonard Vanbrabant and Rebecca Kuiper}

\examples{
## By following these examples, you can appropriately specify hypotheses based on 
## your research questions and analytical framework.

# The hypotheses (i.e., constraints) have to be in a list. It is recommended to name
# each hypothesis in the list. Otherwise the hypotheses are named accordingly 'H1', 'H2', \ldots

# text-based syntax (the labels x1, x2, and x2 are the names of coef(model) or names(vector))
h1 <- '(x1, x2, x3) > 0'
h2 <- '(x1, x3) > 0; x2 = 0'
h3 <- 'x1 > 0; x2 < 0; x3 = 0'
hypotheses = list(hypo1 = h1, hypo2 = h2, hypo3 = h3)

# mixed syntax:  
hypotheses = list(Ha = h1, Hb = 'x1 = x2 > x3')

# the same set of hypotheses for each study:
# hypotheses = list(H1, H2, \ldots)

# a different set of hypotheses for each study: 
# note that the list names set1 and set2 are redundant and can be left out. 
# It is crucial to ensure that the hypotheses across each set are ordered in a similar manner.

H11 <- 'group1 = group2 > group3' 
H12 <- 'group2 > group1 > group3'   

H21 <- 'gr1 = gr2 > gr3'
H22 <- 'gr2 > gr1 > gr3'

# correct
hypotheses = list(set1 = list(H11, H12), set2 = list(H21, H22))

# NOT correct
hypotheses = list(set1 = list(H12, H11), set2 = list(H21, H22))

## Example 1 - 4 studies
est_1 <- c(beta1 = 0.09)
est_2 <- c(beta1 = 0.14)
est_3 <- c(beta1 = 1.09)
est_4 <- c(beta1 = 1.781)
Param_studies <- list(est_1, est_2, est_3, est_4)

# standard error of the beta's (from the primary studies)
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
                    type = "added", 
                    comparison = "none")
evS4_added
summary(evS4_added)

plot(evS4_added)

evS4_equal <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                    hypotheses = hypotheses,
                    type = "equal", 
                    comparison = "none")

evS4_equal
summary(evS4_equal)
plot(evS4_equal)

## Example 2 - 2 studies
est_1 <- c(1.88, 2.54, 0.02)
names(est_1) <- c("group1", "group2", "group3")
vcov_est_1 <- diag(c(0.2149074, 0.2149074, 0.1408014))

est_2 <- c(0.98, 0.02, 0.27)
names(est_2) <- c("gr1", "gr2", "gr3") 
vcov_est_2 <- diag(c(0.1382856, 0.1024337, 0.0987754))

# beta values from the analyses
object <- list(est_1, est_2)
# standard error of the beta's (from the S primary studies)
VCOV <- CovMx_studies <- list(vcov_est_1, vcov_est_2)

# names(est_1) # Specify restrictions using those names
H11 <- 'group1 = group2 > group3'
H12 <- 'group2 > group1 > group3'

# names(est_2) # Specify restrictions using those names
H21 <- 'gr1 = gr2 > gr3'
H22 <- 'gr2 > gr1 > gr3' 

# hypotheses
hypotheses <- list(H1 = list(H11, H12), H2 = list(H21, H22))

evS2_added <- evSyn(object, VCOV = VCOV, hypotheses = hypotheses,
                    type = "added", comparison = "unconstrained") 
evS2_added
plot(evS2_added)

## Example 3 - 3 studies

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

# Evidence synthesis
evS3 <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
              hypotheses = list(H1 = hypothesis),
              type = "added", 
              comparison = "complement") 
evS3
plot(evS3)

## Example 4 - loglikelihood values and penalty values
# make it a list
LL <- as.list(data.frame(t(myLLs)))
penalty.values <- as.list(data.frame(t(myPTs)))

evS_LL_added <- evSyn(object = LL, PT = penalty.values, type = "added")
evS_LL_equal <- evSyn(object = LL, PT = penalty.values, type = "equal")

evS_LL_added
evS_LL_equal


## Example 5 - AIC, ORIC, GORIC(A) values
goric.values <- as.list(data.frame(t(myGORICs)))

evS_Gv <- evSyn(goric.values)
evS_Gv
}


