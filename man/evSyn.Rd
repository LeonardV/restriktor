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
\alias{evSyn}
\alias{leave1studyout}
\alias{leave1studyout.evSyn}
\alias{leave1studyout.default}

\title{GORIC(A) Evidence synthesis}
\description{GORIC(A) evidence synthesis aggregates the evidence for theory-based 
hypotheses from multiple studies that may use diverse designs to investigate the 
same central theory.

The function \code{leave1studyout} performs a leave-one-study-out sensitivity (outlier) analysis 
on a fitted \code{evSyn} object: it re-runs the evidence synthesis S times, each 
time omitting one study, and reports which hypothesis is preferred in each 
reduced set.
}


\usage{
# I try to guess the input type.
evSyn(object, input_type = NULL, \ldots)

# input: Parameter estimates and covariance matrix
evSyn_est(object, \ldots, VCOV = list(), hypotheses = list(),
          type_ev = c("added", "equal", "average"), 
          comparison = c("unconstrained", "complement", "none"),
          hypo_names = c(), priorICweights = NULL, 
          type = c("gorica", "goricac"), 
          order_studies = c("input_order", "ascending", "descending"),
          study_names = c(), study_sample_nobs = NULL, study_weights = NULL,
          penalty_factor = 2, CI = FALSE, CI_level = 0.95, CI_B = 1000L)

# input: Log-likelihood and penalty values
evSyn_LL(object, \ldots, PT = list(), type_ev = c("added", "equal", "average"),
         hypo_names = c(), priorICweights = NULL, 
         type = c("goric", "goricc", "gorica", "goricac"),
         order_studies = c("input_order", "ascending", "descending"),
         study_names = c(), study_weights = NULL, penalty_factor = 2)

# input: GORIC(A), ORIC, AIC values
evSyn_ICvalues(object, \ldots, type_ev = c("added", "average"), 
               hypo_names = c(), priorICweights = NULL, 
               type = c("goric", "goricc", "gorica", "goricac"), 
               order_studies = c("input_order", "ascending", "descending"), 
               study_names = c(), study_weights = NULL)

# input: AIC or ORIC or GORIC or GORICA weights; or (Bayesian) posterior 
# model probabilities
evSyn_ICweights(object, \ldots, type_ev = c("added", "average"), 
                hypo_names = c(), priorICweights = NULL, 
                type = c("goric", "goricc", "gorica", "goricac"),
                order_studies = c("input_order", "ascending", "descending"), 
                study_names = c(), study_weights = NULL)

# input: Ratio of AIC or ORIC or GORIC or GORICA weights; or 
# Bayes factors
evSyn_ICratios(object, \ldots, type_ev = c("added", "average"), 
               hypo_names = c(), priorICweights = NULL, 
               type = c("goric", "goricc", "gorica", "goricac"),
               order_studies = c("input_order", "ascending", "descending"), 
               study_names = c(), study_weights = NULL)

# input: result from the goric() function
# Note that this is a wrapper function for evSyn_LL.
evSyn_gorica(object, \dots, type_ev = c("added", "equal", "average"), 
             hypo_names = c(), priorICweights = NULL, 
             order_studies = c("input_order", "ascending", "descending"), 
             study_names = c(), study_weights = NULL,
             CI = FALSE, CI_level = 0.95, CI_B = 1000L)

# input: Result from the escalc() function from the metafor package (or a 
# regular data.frame using the same labelling). 
# Note that this is a wrapper function for evSyn_est.
evSyn_escalc(data, yi_col = "yi", vi_cols = "vi", study_sample_nobs = "ni", 
             type = NULL, cluster_col = c("trial", "study", "author", "authors", 
             "Trial", "Study", "Author", "Authors"), outcome_col = NULL, \ldots)

\method{print}{evSyn}(x, digits = max(3, getOption("digits") - 4), \dots)

\method{summary}{evSyn}(object, probs = c(0, .25, .50, .75, 1), \dots)

\method{print}{summary.evSyn}(x, digits = max(3, getOption("digits") - 4), \dots)  

\method{plot}{evSyn}(x, output_type = "gorica_weights", xlab = NULL, 
xlab_unordered = NULL, angle_x = 30, \dots)

\method{leave1studyout}{evSyn}(object, \dots)
}

\arguments{
  \item{object}{Currently, the following objects can be processed:
    \itemize{
      \item a list of vectors with (standardized) parameter estimates (the \code{VCOV} argument is required);   
      \item a list of vectors with log-likelihood values (the \code{PT} argument is required);
      \item a list of vectors with GORIC(A) values;
      \item a list of vectors with GORIC(A) weights;
      \item a list of vectors with ratios of GORIC(A) weights (versus a reference hypothesis / central theory);
      \item a list of objects of class goric;
      \item a data.frame of class escalc from the \pkg{metafor} package or a regular data.frame using the same labelling;
    }
    In the case of the first or last option, \code{evSyn} will calculate the GORICA (incl. log-likelihood and penalty) for each study using \code{goric}.
  }
  
  \item{input_type}{
  \code{character} Specifies the type of input provided to the function. Valid options are:
  \describe{
    \item{\code{"est_vcov"}}{Indicates that the input consists of estimates and covariance matrices. Invokes the \code{evSyn_est()} function.}
    \item{\code{"ll_pt"}}{Indicates that the input consists of log-likelihoods and penalty values. Invokes the \code{evSyn_LL()} function.}
    \item{\code{"icvalues"}}{Indicates that the input consists of IC values. Invokes the \code{evSyn_ICvalues()} function.}
    \item{\code{"icweights"}}{Indicates that the input consists of information criterion (IC) weights. Invokes the \code{evSyn_ICweights()} function.}
    \item{\code{"icratios"}}{Indicates that the input consists of IC weights ratios (versus a reference hypothesis / central theory). Invokes the \code{evSyn_ICratios()} function.}
    \item{\code{"gorica"}}{Indicates that the input is of class goric from the \code{goric} function. Invokes the \code{evSyn_gorica()} function.}
    \item{\code{"escalc"}}{Indicates that the input is of class escalc from the \pkg{metafor} package (or a regular data.frame using the same labelling). Invokes the \code{evSyn_escalc()} function.}
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
  information. Notably, when the input is a list of i) GORIC(A) values, ii) GORIC(A) weights, 
  or iii) ratios of GORIC(A) weights, equal-evidence approach is not possible 
  (the added-evidence approach is then used instead).}

  \item{type}{If \code{"gorica"}, a proxy for the log-likihood is calculated by assuming that 
  the parameter estimates (of interest) are multivariate normally distributed. 
  %Then, either the unbiased covariance matrix of the estimates (based on 'N'; e.g., for an lm object) is used
  %or the biased one obtained via the vcov function (e.g., for a \code{lavaan} object). 
  %In the latter case, a corresponding message will be printed. 
  If \code{"goricac"}, then the small-sample-size correction is applied.}

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

  \item{hypo_names}{character vector for labelling the hypotheses. By default, the
  names are set to H1, H2, \ldots}

  \item{priorICweights}{vector that represents the prior belief for the specified hypotheses. 
  By default, equal prior weights are used (i.e., 1/(#hypotheses)). 
  Notably, in case the prior weights do not sum to 1, it will be rescaled such that it does; 
  the prior weights thus reflect relative a priori likeliness of the hypotheses.
  %% For now, like this; if code changes, then perhaps:
  %In the case of ratios of IC weights ('icratios'), the ratio of prior weights should be given.
  %In that case, they do not sum to 1; and 
  %the ratio of weights for the reference hypothesis equals 1.
  Notably, the priorICweights will be used in the cumulative and final results of evidence synthesis.
  They will not be used in ordering the studies with 'order_studies' = \code{"ascending"} or \code{"descending"}; 
  this will be based equal a priori likeliness of the hypotheses.
  Notably, with quite unequal priorICweights (or study_weights), 
  the hypothesis used for sorting/ordering can differ from the one that is preferred in the final weighted result. }
  
  \item{output_type}{a plot option. If "\code{gorica_weights}" (default), 
  a plot with the cumulative goric(a)
  weights and goric(a) weights per study is displayed. If "\code{ll_weights}",
  a plot with the cumulative log-likelihood weights and log-likelihood weights 
  per study is displayed.}
  
  \item{xlab_unordered}{a character vector specifying the original (unordered) 
  study labels. This argument should be used when studies have been re-ordered using
  \code{order_studies}. The labels will be reordered internally to match the 
  plotting order.}
  
  \item{angle_x}{numeric value specifying the angle (in degrees) of the x-axis 
  labels; by default, angle_x = 30. Useful when labels are long and overlap.}
  
  \item{study_names}{
  a character vector specifying the names or labels of the individual (unordered) studies.
  When provided, these names will be used throughout the output, including in 
  tables, summaries, and plots (in place of numeric study indices).
  If \code{NULL} (default), studies are automatically numbered 
  (1, 2, ..., S) following the order determined by \code{order_studies}.}

  \item{penalty_factor}{The penalty factor adjusts the penalty in the GORIC(A) formula
  (GORIC(A) = -2 x log-likelihood + penalty_factor x penalty); by default, 2. 
  It is used in all studies and in the evidence synthesis (see also \code{\link{goric}}). 
  For goric objects as input, the penalty factor used in these objects is used 
  (which must be the same for all objects).}

  \item{CI}{only for input consisting of (standardized) estimates and their covariance 
  matrices (input type \code{"est_vcov"} or \code{"escalc"}) or of goric objects 
  (input type \code{"gorica"}): if \code{TRUE}, confidence 
  intervals (CIs) are computed for the final (pooled) GORICA values, GORICA weights, and 
  ratios of GORICA weights, using \code{\link{evSyn_CI}} with the parametric bootstrap 
  (\code{method = "boot"}). They are printed as extra columns (lower and upper bound) for the 
  final weights in the print output, as extra rows in the final results of the summary, and, 
  for the ratios of weights, in the conclusion. Defaults to \code{FALSE}; the CIs can also be 
  computed afterwards with \code{evSyn_CI(object)}, see \code{\link{evSyn_CI}} for other options 
  (as \code{\link{goric_CI}} for \code{\link{goric}}). They are not available for 
  \code{type_ev = "average"} (since the analyses are then dependent) and not for the other input types 
  (which contain no information about the sampling variability). 
  The state of the random number generator is the same afterwards as without the CIs. 
  See also \code{\link{evSyn_forest}} for a forest-type plot of the study-specific and final 
  ratios of weights with their CIs.}

  \item{CI_level}{the confidence level of the CIs (default 0.95).}

  \item{CI_B}{the number of bootstrap draws (per study) used for the CIs (default 1000).}

  \item{study_sample_nobs}{a numeric vector specifying the sample size (\code{n}) 
  for each study; only needed if \code{type = "goricac"}. 
  In case the input is a escalc object (or a dataframe), it is inherited from the 
  'study_sample_nobs' column, which by default is labelled "ni".}
  
  \item{study_weights}{ 
  a numeric vector specifying the weight for each study; all weights must be larger than 0. 
  Only the relative sizes of the weights matter: the weights are not rescaled and 
  the results are the same for, for example, \code{c(1, 2, 1)}, \code{c(.25, .5, .25)}, 
  and \code{c(.75, 1.5, .75)}. 
  If the weights do not sum to 1 or S (the number of studies/analyses), 
  a message is given to make the user aware of this.
  It is in the same order as the studies/analyses in \code{object} 
  (i.e., the input order, also when \code{order_studies} changes the order of the studies in the output).
  If \code{NULL} (default), all studies are assumed to have equal weight. 
  Note that study-specific sample size is already taken into account via 'VCOV'.
  When provided, the study weights will be used to conduct weighted evidence synthesis.
  %Weighted evidence synthesis is possible when the GORIC(A) values can be calculated.
  %That is, it is not possible when the input consists of IC weights or IC weights ratios
  %(i.e., 'input_type' = \code{"icweights"} or \code{"icratios"});
  %then, the study weights cancel out.
  Notably, the study weights will be used in the cumulative and final results of evidence synthesis, 
  namely in calculating the (cumulative and final) IC weights (and ratios of IC weights): 
  in the added-evidence approach, these are based on weighted sums of the IC values; 
  in the equal-evidence approach, on a weighted sum of the log-likelihood values and a weighted average of the penalty values; 
  and, in the average-evidence approach, on weighted averages of the IC values.
  When the input consists of IC weights or ratios of IC weights, the study weights act as powers of these.
  The study weights are not used in the IC values themselves (which thus still equal 
  -2 x log-likelihood + penalty_factor x penalty), nor in the log-likelihood and penalty values, 
  nor in the log-likelihood weights. 
  Hence, with unequal study weights, the IC weights are not the IC weights of the reported IC values. 
  They will not be used in the ordering of the studies/analyses itself (which happens if 'order_studies' = \code{"ascending"} or \code{"descending"}). 
  The re-ordering will be based on the height of support for each study, thus assuming equal study weights.
  Notably, with quite unequal study_weights (or priorICweights), 
  the hypothesis used for sorting/ordering can differ from the one that is preferred in the final weighted result. 
  }

  \item{order_studies}{
  either a i) character string or ii) numeric vector:
  i) A character string that determines how studies are ordered in the output and plots. 
  If \code{"input_order"} (default), studies are displayed in the order they were 
  provided. If \code{"ascending"} or \code{"descending"}, studies are ordered based 
  on the ascending or descending support for the preferred hypothesis, 
  respectively. 
  In this case, the order will be determined using equal a priori likeliness of the hypotheses
  and equal study weights.
  ii) A numeric vector with length the number of studies (S). 
  This vector then gives the ordering of studies to be used, which one could base on, 
  say, year of the study.
  Notably, when order_studies equals "input_order", the ordering 1:S is used.}
  
  \item{xlab}{a vector specifying custom labels for the x-axis (be aware that 
  study ordering may be changed based on \code{order_studies}). The length of 
  the vector must match the number of studies in the dataset. If not provided, 
  the x-axis labels default to a sequence from 1 to the number of studies 
  (possibly re-ordered based on \code{order_studies}).}

  \item{data}{an object of class "escalc" from the \pkg{metafor} package or a regular data.frame using the same labelling.}

  \item{cluster_col}{a character string specifying the column in \code{data} that 
  contains the cluster/study identifier. Observations within the same cluster/study should 
  share the same value in this column. 
  When not specified, the function looks for a column with one of the following names:
  "trial", "study", "author", "authors", "Trial", "Study", "Author", or "Authors".
  }
  
  \item{outcome_col}{a column name in \code{data} containing the outcome 
  identifiers/labels for each observation within a cluster/study. 
  The hypothesis/-es should be based on this labelling.
  By default (i.e., if 'outcome_col' is null), the function assumes that the 
  parameter label used in the hypothesis is 
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
  in the case of multiple (namely, O) variables, \code{c("v1i", \ldots, "vOi")}.
  When specifying other column names for multiple outcomes, give one column per outcome 
  (otherwise, \code{"v1i"}, \ldots, \code{"vOi"} are used); 
  for studies with one outcome, the first element of \code{vi_cols} is used.
  For more information, see Example 5 below.
  }
  
  \item{digits}{the number of significant digits to use when printing.}

  \item{probs}{numeric vector of probabilities (in [0, 1]) for the percentiles 
  that are reported in the variability section of the summary (see Value); 
  by default, the minimum, first quartile, median, third quartile, and maximum. 
  For example, \code{probs = c(0, .025, .25, .50, .75, .975, 1)} adds a 95\% range; 
  notably, with a small number of studies, the 2.5\% and 97.5\% percentiles 
  lie (very) close to the minimum and maximum.}

  \item{\ldots}{This depends on the class of the object.}
}

\details{
In the \strong{added-evidence} approach, evidence for a central hypothesis/theory 
from each study or dataset is cumulatively aggregated. 
This means that for every new study, the log-likelihood and the penalty term are 
added to the cumulative totals. Stated differently, the (cumulative) IC weights 
are used as an a priori weigth for the next study.
\cr 
Opt for this approach when you believe that the 
aggregated evidence from multiple studies is stronger than if the data were combined 
into a single study (when that would be possible).

The \strong{equal-evidence} approach also aggregates the cumulative evidence 
central hypothesis/theory but as if the data were combined into a single study 
(when that would be possible).  
In that way, it is comparable to meta-analysis but without the need that the 
data can be combined. 
This means that for every new study, the log-likelihood 
part (reflecting the fit of the hypothesis to the data) is added to the 
cumulative total, while the penalty part (reflecting the size of the hypothesis) 
is averaged. Because of this, this approach can only be used when there is 
information about the log-likelihood and penalty (and not only the IC values or 
weights).
\cr
Choose this method when you view the combined evidence from multiple studies as 
equivalent to that from a single, larger study.

The \strong{average-evidence} method can be conceptualized as a form of multiverse 
analysis. When faced with a single dataset, there are often numerous analytical 
choices available, such as handling missing data, selecting variables, or choosing 
statistical methods. Each choice can lead to a different analysis or model, creating 
a "multiverse" of possible outcomes.
\cr
For each of these analyses, evidence for your hypothesis/-es of interest can be 
calculated using IC values or weights, indicating how likely your hypothesis is. 
In some analyses, a hypothesis might obtain high support, while in others it might not. 
The average-evidence method aggregates the evidence for a hypothesis from all analyses, 
providing an average measure of evidence across all considered analyses.
\cr
Opt for the average-evidence approach when you want to combine evidence from 
multiple analyses on one dataset.

\strong{Leave-one-out sensitivity analysis.}
\code{leave1studyout} repeats the evidence synthesis S times, each time dropping 
one of the S studies. For each reduced set, the aggregated IC values per hypothesis 
are computed using the same \code{type_ev} (and the same priorICweights and study_weights) 
as in the original \code{evSyn} call:
\itemize{
  \item \code{"added"}: column sums of the IC matrix with study \emph{s} removed.
  \item \code{"equal"}: the column means of the PT matrix minus column sums of the LL matrix with study \emph{s}
    removed.
  \item \code{"average"}: column means of the IC matrix with study \emph{s} removed.
}
As in \code{evSyn}, the study weights (if specified) are used in calculating the 
aggregated IC weights, but not in the reported aggregated IC values.
The hypothesis with the highest aggregated IC weight 
is identified as the preferred hypothesis for that reduced set 
(with equal priorICweights and study weights, this is the one with the lowest aggregated IC value). 
\cr
When the input of \code{evSyn} consists of IC weights (\code{"icweights"}), 
the IC weights are aggregated directly (as in \code{evSyn}); 
no IC values are available and thus reported then.
When the input consists of ratios of IC weights (\code{"icratios"}), 
the aggregated differences in IC values versus the reference hypothesis are reported 
(instead of IC values).
\cr
Notably, when leaving a study out, the combined sample size lowers. 
Therefore, the maximum log-likelihood goes down and thus the IC values go up (become less negative).
Since the sample sizes of the studies may differ per study, the combined sample size 
for each of the S evidence syntheses may differ.
Therefore, it is not fair to compare the rendered IC values across the S evidence syntheses. 
Although research is needed, one can best compare the rendered IC weights.
\cr
This analysis can be used to detect influential 
studies and to assess the robustness of the overall conclusion.
}


\value{An object of class evSyn for which a print, summary and plot function is available. 
The output comprises, among other things, the cumulative and final evidence 
for the theory-based hypotheses.
For confidence intervals for the final (pooled) results (when the input consists of 
(standardized) estimates and their covariance matrices or of goric objects), see the argument \code{CI} and 
\code{\link{evSyn_CI}}; for a forest-type plot of the study-specific and final ratios of 
weights with their CIs, see \code{\link{evSyn_forest}}.

The summary (with at least two studies/analyses) also contains measures of 
variability over the studies (or analyses), in the element \code{variability}. 
These are based on the study-specific IC weights (including the priorICweights) 
and compare the overall best hypothesis (i.e., the one with the highest final IC weight) 
with each other hypothesis:
\itemize{
  \item percentiles (see \code{probs}) of the study-specific IC weights per hypothesis 
    and of the study-specific ratios of IC weights (best versus each other hypothesis);
  \item the geometric mean of these ratios (cf. the gP-BF in Klaassen et al., 2018); 
    for the average-evidence approach, this equals the final ratio of IC weights;
  \item the evidence rate (ER; Klaassen et al., 2018): the proportion of studies 
    whose ratio favours the same hypothesis as the geometric mean 
    (i.e., > 1 when the geometric mean > 1). 
    An ER near 1 indicates homogeneous studies; Klaassen et al. regard ER >= .9 as sufficient;
  \item the stability rate (SR; Klaassen et al., 2018): the proportion of studies 
    whose ratio is more extreme than the geometric mean. 
    An SR around .5 (Klaassen et al. use .4 to .6) indicates that the ratios are evenly 
    spread around the geometric mean; below .5, a minority of studies with large ratios 
    inflates the geometric mean; above .5, a minority of studies with small ratios deflates it;
  \item the preference rate: the proportion of studies in which the overall best 
    hypothesis also has the highest study-specific IC weight.
}
With unequal study weights, the geometric mean, ER, SR, and preference rate are weighted 
(i.e., weighted geometric mean and weighted proportions); the percentiles are not.
Notably, when the best hypothesis overlaps with another one (e.g., the best hypothesis 
versus the unconstrained hypothesis when the estimates are in agreement with the best one), 
the ratio of IC weights has an upper bound (determined by the difference in penalties), 
so that many studies may obtain the same (maximum) ratio; ER and SR should then be 
interpreted with care.
When log-likelihood values are available, the number of studies in which the 
best and the other hypothesis have the same log-likelihood value (column \code{Ceiling}) is reported as well.

Based on these measures, the summary also gives a conclusion (element \code{variability$conclusion}), 
using the cut-off values suggested by Klaassen et al. (2018): ER >= .90 and .40 <= SR <= .60 
(which are subjective guidelines, not strict rules). As for the \code{goric} output, this is done in two steps:
first, the best hypothesis is compared with the failsafe hypothesis (unconstrained or complement), 
to check whether it is not weak (if the failsafe is the best, all order-restricted hypotheses are weak); 
second, if it is not weak, it is compared with the other order-restricted hypotheses in the set. 
For each comparison, the final ratio of IC weights is given, together with the ER and SR verdicts 
and, if applicable, the number of studies in which the ratio reached its maximum (ceiling).
If there is no failsafe hypothesis in the set (i.e., no hypothesis labelled 'Unconstrained' or 'Complement'), 
it is assumed that the hypotheses in the set together cover all possibilities (i.e., the full parameter space); 
if not, a failsafe hypothesis should be added to the set (preferably the default one). With fewer than 10 studies, a note is given that ER >= .90 then requires 
that all studies point in the same direction.

\code{leave1studyout.evSyn} returns a named list with the following components:
\describe{
  \item{\code{OverallGorica} (or \code{OverallGoric}, \code{OverallGoricc}, 
    \code{OverallGoricac})}{a numeric matrix of dimension S x H (studies by 
    hypotheses) containing the aggregated IC values when the corresponding study 
    is left out. Row names are of the form \code{"Leave '<study name>' out:"}. 
    Not available when the input of \code{evSyn} consists of IC weights; 
    when it consists of ratios of IC weights, \code{OverallICdiff} is given instead, containing 
    the aggregated differences in IC values versus the reference hypothesis.}
  \item{\code{OverallGoricaWeights} (or \code{OverallGoricWeights}, \code{OverallGoriccWeights}, 
    \code{OverallGoricacWeights})}{a numeric matrix of dimension S x H (studies by 
    hypotheses) containing the aggregated IC weights when the corresponding study 
    is left out. Row names are of the form \code{"Leave '<study name>' out:"}.}
  \item{\code{OverallPrefHypo}}{a character matrix of dimension S x 1 giving the 
    name of the preferred hypothesis (i.e., the one with the highest aggregated IC 
    weight) for each leave-one-out set.}
  \item{\code{type_ev}}{the evidence-synthesis approach used, inherited from the 
    input \code{evSyn} object.}
  \item{\code{type}}{the IC type used (\code{"goric"}, \code{"goricc"}, 
    \code{"gorica"}, or \code{"goricac"}), if available in the input object.}
}
}


\references{
  Klaassen, F., Zedelius, C. M., Veling, H., Aarts, H., and Hoijtink, H. (2018). 
  All for one or some for all? Evaluating informative hypotheses using multiple N = 1 studies. 
  \emph{Behavior Research Methods, 50}, 2276--2291. \doi{10.3758/s13428-017-0992-5}
}

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

# GORICA evidence synthesis (added-approach)
evS4_added <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                    hypotheses = hypotheses,
                    type_ev = "added", # default
                    comparison = "none")

# Output
# Final results:
evS4_added
# Study-specific, cumulative, and final results:
summary(evS4_added)
# evSyn plot:
plot(evS4_added) # Note that this renders a ggplot which you can change accordingly, as done below.


# GORICA evidence synthesis (equal-approach)
evS4_equal <- evSyn(object = Param_studies, VCOV = CovMx_studies, 
                    hypotheses = hypotheses,
                    type_ev = "equal", 
                    comparison = "none")
                    
# Output
# Final results:
evS4_equal
# Study-specific, cumulative, and final results:
summary(evS4_equal)
# evSyn plot:
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

# GORICA vidence synthesis (added-evidence)
evS2_added <- evSyn(object, VCOV = VCOV, hypotheses = hypotheses,
                    type_ev = "added", # default
                    comparison = "unconstrained") # default (if multiple hypotheses)

# Output
# Final results:
evS2_added
# Study-specific, cumulative, and final results:
summary(evS2_added)
# evSyn plot:
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
              type_ev = "added", # default
              comparison = "complement") # default (if one hypothesis)

# Output
# Final results:
evS3
# Study-specific, cumulative, and final results:
summary(evS3)
# evSyn plot:
plot(evS3)



## Example 4: based on loglikelihood values and penalty values

# Make LL and PT values a list
LL <- as.list(data.frame(t(myLLs)))
penalty.values <- as.list(data.frame(t(myPTs)))


# Evidence synthesis

# Added-evidence
evS_LL_added <- evSyn(object = LL, PT = penalty.values, 
                      hypo_names = colnames(myLLs),
                      type_ev = "added") # default
evS_LL_added
plot(evS_LL_added)

# Equal-evidence
evS_LL_equal <- evSyn(object = LL, PT = penalty.values, 
                      hypo_names = colnames(myLLs),
                      type_ev = "equal")
evS_LL_equal
plot(evS_LL_equal)



## Example 5: based on AIC, ORIC, GORIC(A) values

# Make GORIC values a list
goric.values <- as.list(data.frame(t(myGORICs)))

# Evidence synthesis (added-evidence)
evS_Gv <- evSyn(goric.values, 
                hypo_names = colnames(myGORICs))
evS_Gv
plot(evS_Gv)


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

  # The hypothesis should thus now be based on the default labelling 'theta':
  H1 <- "abs(theta) > 0.5" 
  
  # Evidence synthesis (equal-evidence)
  evS_escalc_a <- evSyn(
    data_a, 
    hypotheses = list(H1 = H1), 
    type_ev = "equal",
    comparison = "complement" # default (if one hypothesis)
  )
  
  # Output
  evS_escalc_a
  plot(evS_escalc_a)
  
  
  #Plot options
  #The evSyn plot is a ggplot and can be changed as such. 
  
  #As an example, one could change the x-axis colors based on the value(s) of another variable (in a way, looking at subgroups).
  #Next, we set the colors based on two subgroups of studies, those published before 1983 (dark-gray) and those published that year or later (black):
  evSynplot <- plot(evS_escalc_a)
  cols_subgroup <- ifelse(data_a$year < 1983, "darkgray", "black")
  library(ggplot2)
  evSynplot_subgroup <- evSynplot + theme(axis.text.x = element_text(colour = cols_subgroup))
  evSynplot_subgroup
  
  # As another example, one can overwrite the subtitle:
  evSynplot_subgroup_subtitle <- evSynplot_subgroup + labs(subtitle = 'For studies published < and >= 1983')
  evSynplot_subgroup_subtitle
  
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
  
  # The hypothesis should thus now be based on the labeling based on the 'outcome' variable:
  H1 <- "abs(PD) > 0.5" 
  
  # Evidence synthesis (equal-evidence)
  evS_escalc_b <- evSyn(data_b, 
                        hypotheses = list(H1 = H1), 
                        type_ev = "equal",
                        comparison = "complement") # default (if one hypothesis)
  
  # Output
  evS_escalc_b
  plot(evS_escalc_b)
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
  # that labeling, that is, "AL" and "PD"
  # yi_col = "yi"
  # vi_cols = c("v1i", "v2i")
   
  # The hypothesis should thus now be based on the labeling based on the 'outcome' variable:
  H1 <- "abs(AL) > 0.5, abs(PD) > 0.5" 
  
  # Evidence synthesis (equal-evidence)
  evS_escalc_c <- evSyn(data_c, 
                        hypotheses = list(H1 = H1), 
                        type_ev = "equal",
                        comparison = "complement") # default (if one hypothesis)
  
  # Output
  evS_escalc_c
  plot(evS_escalc_c)
  
  
  #Plot options
  #The evSyn plot is a ggplot and can be changed as such. 
  
  #As an example, one could change the x-axis colors based on the value(s) of another variable (in a way, looking at subgroups).
  #Next, we set the colors based on two subgroups of studies, those published before 1983 (dark-gray) and those published that year or later (black):
  evSynplot <- plot(evS_escalc_c)
  cols_subgroup <- ifelse(data_a$year < 1983, "darkgray", "black")
  library(ggplot2)
  evSynplot_subgroup <- evSynplot + theme(axis.text.x = element_text(colour = cols_subgroup))
  evSynplot_subgroup
  
  # As another example, one can overwrite the subtitle:
  evSynplot_subgroup_subtitle <- evSynplot_subgroup + labs(subtitle = 'For studies published < and >= 1983')
  evSynplot_subgroup_subtitle

}


## Example 6d: 2 outcomes, assigning study weigths
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
   
  # The hypothesis should thus now be based on the labeling based on the 'outcome' variable:
  H1 <- "abs(AL) > 0.5, abs(PD) > 0.5" 
  
  # Evidence synthesis (equal-evidence)
  evS_escalc_c <- evSyn(data_c, 
                        hypotheses = list(H1 = H1), 
                        type_ev = "equal",
                        comparison = "complement") # default (if one hypothesis)
  
  # Output
  evS_escalc_c
  plot(evS_escalc_c)
}
  

}


