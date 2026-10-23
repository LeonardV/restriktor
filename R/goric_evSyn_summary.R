# Helper: measures of variability over the studies (or analyses), based on the 
# study-specific IC weights (as reported, so without the prior IC weights (priorICweights), which are
# only used in the aggregated results).
# The overall best hypothesis is the one with the highest final (cumulative) IC 
# weight; it is compared with each other hypothesis via the ratio of IC weights.
#   - Percentiles (by default, min, Q1, median, Q3, max) of the study-specific 
#     IC weights (per hypothesis) and of the study-specific ratios of IC weights 
#     (best versus each other hypothesis). These are unweighted, i.e., they 
#     describe the spread over the studies as observed.
#   - Geometric mean of the study-specific ratios (cf. the gP-BF in Klaassen et al., 2018); 
#     with study weights, the weighted geometric mean.
#   - Evidence rate (ER; Klaassen et al., 2018): the (weighted) proportion of 
#     studies whose ratio favours the same hypothesis as the geometric mean 
#     (i.e., > 1 if the geometric mean > 1, and < 1 if it is < 1).
#   - Stability rate (SR; Klaassen et al., 2018): the (weighted) proportion of 
#     studies whose ratio is more extreme than the geometric mean (i.e., larger 
#     than it if the geometric mean > 1, smaller if it is < 1). Values around .5 
#     indicate that the ratios are evenly spread around the geometric mean; 
#     values below .5 that a minority of studies with large ratios inflate it; 
#     values above .5 that a minority of studies with small ratios deflate it.
#   - Preference rate: the (weighted) proportion of studies in which the overall 
#     best hypothesis also has the highest study-specific IC weight.
.evSyn_variability <- function(W, final_w, study_weights = NULL, 
                               probs = c(0, .25, .50, .75, 1), LL = NULL) {
  W <- as.matrix(W)
  S <- nrow(W)
  if (S < 2L || is.null(final_w)) return(NULL)
  hn <- colnames(W)
  if (is.null(hn)) hn <- paste0("H", seq_len(ncol(W)))
  if (is.null(study_weights)) study_weights <- rep(1, S)
  sw <- study_weights / sum(study_weights)
  
  best <- which.max(final_w)
  others <- setdiff(seq_len(ncol(W)), best)
  
  pct <- function(x) quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  pct_names <- paste0(format(100 * probs, trim = TRUE), "%")
  pct_names[probs == 0] <- "Min"
  pct_names[probs == 1] <- "Max"
  pct_names[probs == .5] <- "Median"
  
  # Percentiles of the study-specific IC weights (rows: hypotheses)
  weights_pct <- t(apply(W, 2, pct))
  if (length(probs) == 1L) weights_pct <- t(weights_pct)
  dimnames(weights_pct) <- list(hn, pct_names)
  
  # Study-specific ratios of IC weights: best versus each other hypothesis
  R <- W[, best] / W[, others, drop = FALSE]
  colnames(R) <- paste0(hn[best], " vs. ", hn[others])
  
  ratios_pct <- t(apply(R, 2, pct))
  if (length(probs) == 1L) ratios_pct <- t(ratios_pct)
  dimnames(ratios_pct) <- list(colnames(R), pct_names)
  
  rate_stats <- t(apply(R, 2, function(r) {
    ok <- !is.na(r)
    gm <- exp(sum(sw[ok] * log(r[ok])) / sum(sw[ok]))
    if (is.nan(gm)) gm <- NA_real_
    if (!is.na(gm) && gm >= 1) {
      ER <- sum(sw[ok] * (r[ok] > 1)) / sum(sw[ok])
      SR <- sum(sw[ok] * (r[ok] > gm)) / sum(sw[ok])
    } else if (!is.na(gm)) {
      ER <- sum(sw[ok] * (r[ok] < 1)) / sum(sw[ok])
      SR <- sum(sw[ok] * (r[ok] < gm)) / sum(sw[ok])
    } else {
      ER <- SR <- NA_real_
    }
    c(GeomMean = gm, ER = ER, SR = SR)
  }))
  if (length(others) == 1L && nrow(rate_stats) != 1L) rate_stats <- t(rate_stats)
  rownames(rate_stats) <- colnames(R)
  colnames(rate_stats) <- c("Geom. mean", "ER", "SR")
  
  # Ceiling: number of studies in which the ratio reached its maximum, that is, 
  # in which the best hypothesis and the other one have the same log-likelihood 
  # value (they overlap there); the ratio then only reflects the difference in 
  # penalties. Only available when log-likelihood values are available.
  if (!is.null(LL)) {
    LL <- as.matrix(LL)
    ceil <- vapply(others, function(k) {
      d <- abs(LL[, best] - LL[, k])
      sum(d <= 1e-8 * pmax(1, abs(LL[, best])), na.rm = TRUE)
    }, numeric(1))
    rate_stats <- cbind(rate_stats, Ceiling = ceil)
  }
  
  # Preference rate: (weighted) proportion of studies in which the overall best 
  # hypothesis also has the highest study-specific IC weight
  is_best_s <- apply(W, 1, function(w) isTRUE(w[best] >= max(w, na.rm = TRUE)))
  pref_rate <- sum(sw * is_best_s)
  
  list(best = hn[best], probs = probs,
       weights_percentiles = weights_pct,
       ratios_percentiles = ratios_pct,
       ratios_rates = rate_stats,
       preference_rate = pref_rate,
       weighted = length(unique(round(sw, 12))) > 1)
}

# Helper: conclusion (text) based on the measures of variability over the
# studies (or analyses), in combination with the cut-off values suggested by
# Klaassen et al. (2018): ER >= .90 and .40 <= SR <= .60.
# As for the goric output, this is done in two steps:
#   Step 1: compare the overall best hypothesis with the failsafe hypothesis
#           (unconstrained or complement), to check whether it is not weak.
#   Step 2: if the best hypothesis is not weak, compare it with the other
#           (order-restricted) hypotheses in the set.
# For each comparison, the final (synthesized) ratio of IC weights is given,
# together with how consistent (ER) and how evenly spread (SR) the
# study-specific ratios are; and, if applicable, how many studies reached the
# maximum (ceiling) of the ratio, which happens when the two hypotheses
# overlap in a study (i.e., have the same log-likelihood value).
.evSyn_variability_conclusion <- function(v, final_w, type_label = "GORICA",
                                          type_ev = "added", S,
                                          ER_cut = .90, SR_low = .40, SR_up = .60,
                                          CI = NULL) {
  if (is.null(v) || is.null(final_w)) return(NULL)
  eps <- sqrt(.Machine$double.eps) # tolerance in comparing ER and SR to the cut-off values

  unit    <- .evSyn_unit_label(type_ev, capital = FALSE)
  units   <- .evSyn_unit_label(type_ev, plural = TRUE, capital = FALSE)
  hn      <- names(final_w)
  best    <- v$best
  rates   <- v$ratios_rates
  others  <- hn[hn != best]
  fmt <- function(val) {
    if (is.na(val)) return("NA")
    if (abs(val) >= 1000 || (abs(val) <= 0.001 && val != 0)) sprintf("%.2e", val) else sprintf("%.2f", val)
  }
  pct <- function(p) paste0(round(100 * p), "%")

  # Text " (95%CI: [lower, upper])" for the final ratio of weights of 'm'
  # versus 'k', if CIs are available (evSyn(..., CI = TRUE); see evSyn_CI())
  ci_txt <- function(m, k) {
    if (is.null(CI)) return("")
    rd <- CI$ratio_weights
    j <- which(rd$model == m & rd$vs == k)
    if (length(j) != 1L || is.na(rd$lower[j]) || is.na(rd$upper[j])) return("")
    paste0(" (", format(100 * CI$level), "%CI: [", fmt(rd$lower[j]), ", ", fmt(rd$upper[j]), "])")
  }

  # Identify the failsafe hypothesis (by name)
  is_failsafe <- tolower(hn) %in% c("unconstrained", "complement", "hu", "hunc", "hc")
  failsafe <- if (any(is_failsafe)) hn[which(is_failsafe)[1]] else NULL

  # Text for one comparison (best versus 'other')
  pair_text <- function(other) {
    i  <- match(other, others)
    gm <- rates[i, "Geom. mean"]; ER <- rates[i, "ER"]; SR <- rates[i, "SR"]
    ceil_n <- if ("Ceiling" %in% colnames(rates)) rates[i, "Ceiling"] else NA
    fin <- final_w[best] / final_w[other]
    txt <- paste0("  * Final: ", sQuote(best), " has ", fmt(fin),
                  " times more support than ", sQuote(other), ci_txt(best, other), ".")
    if (is.na(gm) || is.na(ER)) {
      return(c(txt, "    (The variability over the studies could not be determined for this comparison.)"))
    }
    if (gm >= 1) {
      txt <- c(txt, paste0("    In ", pct(ER), " of the ", units, " (ER = ", fmt(ER), "), ",
                           sQuote(best), " has more support than ", sQuote(other), "."))
    } else {
      txt <- c(txt, paste0("    Notably, on average (geometric mean of the ratios = ", fmt(gm), "), ",
                           sQuote(other), " has more support than ", sQuote(best), " in the ", units,
                           "; this holds in ", pct(ER), " of the ", units, " (ER = ", fmt(ER), "). ",
                           "The final preference for ", sQuote(best), " is due to the prior IC weights (priorICweights) and/or the study weights."))
    }
    # ER verdict
    if (ER >= ER_cut - eps) {
      txt <- c(txt, paste0("    Since ER >= ", fmt(ER_cut), ", the ", units,
                           " are sufficiently homogeneous for this comparison."))
    } else {
      txt <- c(txt, paste0("    Since ER < ", fmt(ER_cut), ", the ", units,
                           " are not homogeneous for this comparison: ", pct(1 - ER), " of the ", units,
                           " point in the other direction. ",
                           "Hence, the final result may not hold for every ", unit,
                           "; inspect the ", unit, "-specific results (and possible moderators)."))
    }
    # Ceiling (overlap) and SR verdict
    if (!is.na(ceil_n) && ceil_n > 0) {
      bound_word <- if (tolower(other) %in% c("unconstrained", "hu", "hunc")) "ceiling" else "(possible) bound"
      txt <- c(txt, paste0("    In ", ceil_n, " of the ", S, " ", units, ", ", sQuote(best), " and ", sQuote(other),
                           " have the same fit (log-likelihood value), so that the ratio only reflects ",
                           "the difference in penalties (i.e., in complexity) there. ",
                           "This happens when the estimates are in agreement with both hypotheses. ",
                           if (tolower(other) %in% c("unconstrained", "hu", "hunc")) {
                             paste0("Since ", sQuote(best), " is contained in ", sQuote(other),
                                    ", this ratio is its maximum (ceiling): ",
                                    "it cannot become larger in these ", units, ", whatever the sample size.")
                           } else {
                             paste0("If one of the two hypotheses is contained in the other, this ratio is at its bound ",
                                    "(whatever the sample size): its maximum (ceiling) if ", sQuote(best),
                                    " is contained in ", sQuote(other), ", and its minimum if ", sQuote(other),
                                    " is contained in ", sQuote(best), ".")
                           }))
      if (ceil_n >= S / 2) {
        txt <- c(txt, paste0("    Because of this ", bound_word, ", SR (= ", fmt(SR), ") is not informative ",
                             "for this comparison."))
        return(txt)
      } else {
        txt <- c(txt, paste0("    Because of this ", bound_word, ", interpret SR (= ", fmt(SR), ") with care."))
      }
    }
    if (is.na(ceil_n) || ceil_n == 0 || ceil_n < S / 2) {
      if (SR >= SR_low - eps && SR <= SR_up + eps) {
        txt <- c(txt, paste0("    Since ", fmt(SR_low), " <= SR (= ", fmt(SR), ") <= ", fmt(SR_up),
                             ", the ", unit, "-specific ratios are evenly spread around their geometric mean."))
      } else if (SR < SR_low) {
        txt <- c(txt, paste0("    Since SR (= ", fmt(SR), ") < ", fmt(SR_low), ", a minority of ", units,
                             " with relatively large ratios inflates the geometric mean; ",
                             "the typical ", unit, " shows less support than the average suggests."))
      } else {
        txt <- c(txt, paste0("    Since SR (= ", fmt(SR), ") > ", fmt(SR_up), ", a minority of ", units,
                             " with relatively small ratios deflates the geometric mean; ",
                             "the typical ", unit, " shows more support than the average suggests."))
      }
    }
    txt
  }

  out <- character(0)

  if (!is.null(failsafe) && best == failsafe) {
    # The failsafe is the best: all (order-restricted) hypotheses are weak.
    out <- c(out, paste0("- The failsafe hypothesis ", sQuote(failsafe), " has the highest final ",
                         type_label, " weight. As a result, the order-restricted hypotheses are considered weak, ",
                         "and inspecting their relative support is not meaningful."),
             paste0("- The variability over the ", units, " for this conclusion (", sQuote(failsafe),
                    " versus each order-restricted hypothesis):"))
    for (o in others) out <- c(out, pair_text(o))
  } else if (!is.null(failsafe)) {
    # Step 1: best versus failsafe
    out <- c(out, paste0("- Step 1: The order-restricted hypothesis ", sQuote(best),
                         " is the best in the set, as it has the highest final ", type_label, " weight. ",
                         "Since it has a higher weight than the failsafe hypothesis ", sQuote(failsafe),
                         ", it is not considered weak:"))
    out <- c(out, pair_text(failsafe))
    i_fs <- match(failsafe, others)
    if (rates[i_fs, "ER"] < ER_cut - eps) {
      out <- c(out, paste0("    Notably, since ER < ", fmt(ER_cut), ", ", sQuote(best),
                           " is weak (i.e., not better than ", sQuote(failsafe), ") in some ", units,
                           "; keep this in mind in Step 2."))
    }
    # Step 2: best versus the other (order-restricted) hypotheses
    rest <- setdiff(others, failsafe)
    if (length(rest) > 0) {
      out <- c(out, paste0("- Step 2: We can now inspect the relative support for ", sQuote(best),
                           " against the other order-restricted hypotheses:"))
      for (o in rest) out <- c(out, pair_text(o))
    }
  } else {
    # No failsafe in the set
    out <- c(out, paste0("- The hypothesis ", sQuote(best), " is the best in the set, as it has the highest final ",
                         type_label, " weight. Since there is no failsafe hypothesis (unconstrained or complement) ",
                         "in the set, it is assumed that the hypotheses in the set together cover all possibilities ",
                         "(i.e., the full parameter space). If not, a failsafe hypothesis should be added to the set ",
                         "(preferably the default one), to protect against choosing the best out of a set of weak hypotheses."),
             paste0("- The relative support for ", sQuote(best), " against the other hypotheses:"))
    for (o in others) out <- c(out, pair_text(o))
  }

  out <- c(out, paste0("- In ", pct(v$preference_rate), " of the ", units, ", ", sQuote(best),
                       " also has the highest ", unit, "-specific ", type_label, " weight (preference rate)."))

  # Cut-off values and small number of studies
  out <- c(out, paste0("Note: The cut-off values ER >= ", fmt(ER_cut), " and ", fmt(SR_low), " <= SR <= ", fmt(SR_up),
                       " are those suggested by Klaassen et al. (2018); they are subjective guidelines, not strict rules."))
  if (S < 10) {
    out <- c(out, paste0("Note: With S = ", S, " ", units, ", ER can only take the values 0, 1/", S, ", ..., 1; ",
                         "so ER >= ", fmt(ER_cut), " requires that all ", units, " point in the same direction, ",
                         "and one ", unit, " changes ER (and SR) by ", fmt(1 / S), "."))
  }
  if (isTRUE(v$weighted)) {
    out <- c(out, "Note: ER, SR, and the preference rate are weighted by the study weights.")
  }
  if (!("Ceiling" %in% colnames(rates))) {
    out <- c(out, paste0("Note: For this type of input (no log-likelihood values), it cannot be determined whether ",
                         "ratios reached their maximum (ceiling); this happens when two hypotheses overlap in a ", unit,
                         " (e.g., the best hypothesis versus the unconstrained one when the estimates are in agreement ",
                         "with the best one)."))
  }
  out
}

# -------------------------------------------------------------------------
# Measures of variability over the studies (or analyses) and the corresponding 
# conclusion, for an evSyn object (see also evSyn(..., variability = TRUE)).
evSyn_variability <- function(object, probs = c(0, .25, .50, .75, 1)) {
  if (!inherits(object, "evSyn")) {
    stop("\nrestriktor ERROR: The object should be of class 'evSyn' ",
         "(i.e., output from the evSyn() function).", call. = FALSE)
  }
  if (!is.numeric(probs) || any(is.na(probs)) || any(probs < 0 | probs > 1)) {
    stop("\nrestriktor ERROR: The argument 'probs' should contain probabilities (in between 0 and 1).",
         call. = FALSE)
  }
  x <- object
  type_label <- if (!is.null(x[["type"]])) toupper(x[["type"]]) else "GORICA"
  S <- x[["n_studies"]]
  if (is.null(x[["GORICA_weight_m"]]) || is.null(S) || S < 2) {
    message("\nrestriktor Message: Measures of variability require at least two ",
            .evSyn_unit_label(x[["type_ev"]], plural = TRUE, capital = FALSE), ".")
    return(invisible(NULL))
  }
  # based on the study-specific IC weights; the study weights (in input order) 
  # are re-ordered to follow the (possibly re-ordered) studies.
  sw_var <- x[["study_weights"]]
  if (!is.null(sw_var) && !is.null(x[["order_studies"]])) sw_var <- sw_var[x[["order_studies"]]]
  final_w_var <- if (!is.null(x[["Cumulative_GORICA_weights"]])) x[["Cumulative_GORICA_weights"]]["Final", ] else NULL
  v <- .evSyn_variability(x[["GORICA_weight_m"]], final_w_var, study_weights = sw_var, probs = probs,
                          LL = x[["LL_m"]])
  if (is.null(v)) return(invisible(NULL))
  # Conclusion based on these measures (with the cut-off values of Klaassen et al., 2018)
  v$conclusion <- .evSyn_variability_conclusion(v, final_w_var, type_label = type_label,
                                                type_ev = x[["type_ev"]], S = S, CI = x[["CI"]])
  v$type_ev <- x[["type_ev"]]
  v$type_label <- type_label
  class(v) <- "evSyn_variability"
  v
}

print.evSyn_variability <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  .evSyn_print_variability(x, type_ev = x$type_ev, type_label = x$type_label, digits = digits)
  invisible(x)
}


summary.evSyn <- function(object, probs = c(0, .25, .50, .75, 1),
                          variability = isTRUE(object[["show_variability"]]), ...) {
  x <- object
  class(x) <- NULL
  
  # Determine label: default to "GORICA" unless type is explicitly given
  # Use [[ ]] to avoid partial matching (e.g., $type matching $type_ev)
  if (!is.null(x[["type"]])) {
    type_label <- switch(
      x$type,
      goric   = "GORIC",
      goricc  = "GORICC",
      gorica  = "GORICA",
      goricac = "GORICAC"
    )
  } else {
    type_label <- "GORICA"
  }
  
  ans <- list(
    type = x[["type"]],
    type_ev = x[["type_ev"]],
    #n_studies = nrow(x[["GORICA_weight_m"]][,, drop = FALSE]),
    n_studies = x[["n_studies"]], 
    order_studies = x[["order_studies"]],
    study_names = x[["study_names"]],
    study_weights = x[["study_weights"]],
    hypotheses = x[["hypotheses"]],
    GORICA_weight_m = x[["GORICA_weight_m"]],
    GORICA_m = x[["GORICA_m"]],
    LL_weights_m = x[["LL_weights_m"]],
    Cumulative_LL_weights = x[["Cumulative_LL_weights"]],
    Cumulative_LL = x[["Cumulative_LL"]],
    LL_m = x[["LL_m"]],
    PT_m = x[["PT_m"]],
    Cumulative_GORICA_weights = x[["Cumulative_GORICA_weights"]],
    Cumulative_GORICA = x[["Cumulative_GORICA"]],
    #
    # If icratios
    Href        = x[["Href"]],
    ICdiff_m    = x[["ICdiff_m"]], # diff in IC values versus reference hypo
    Cumulative_ICdiff = x[["Cumulative_ICdiff"]], # cum diff in IC values versus reference hypo
    Cumulative_ratioICweights = x[["Cumulative_ratioICweights"]], # cum ratios
    GwRatio_m   = x[["GwRatio_m"]] # ratio of IC weights!
  )
  
  if (!is.null(x[["PT_m"]])) {
    sequence <- .evSyn_sequence_labels(ans$n_studies, x[["type_ev"]])
    Cumulative_PT <- apply(x[["PT_m"]][,, drop = FALSE], 2, cumsum)  
    Cumulative_PT <- matrix(Cumulative_PT, nrow = nrow(x[["PT_m"]]), 
                            dimnames = list(sequence, colnames(x[["PT_m"]])))
    if (x[["type_ev"]] %in% c("equal", "average")) {
      Cumulative_PT <- Cumulative_PT / seq_len(ans$n_studies)
    } 
    ans$Cumulative_PT <- Cumulative_PT
  }
  
  final <- c()
  
  # Confidence intervals for the final results (if evSyn(..., CI = TRUE); see evSyn_CI()):
  # rows with the lower and upper bounds, below the final IC weights and values.
  ci_rows <- function(what, hyps) {
    ci <- x[["CI"]]
    if (is.null(ci)) return(NULL)
    tab <- if (what == "weights") ci$weights else ci$IC
    idx <- match(hyps, tab$model)
    if (anyNA(idx)) return(NULL)
    lab <- paste0(type_label, " ", what, ": ", format(100 * ci$level), "%CI ")
    out <- rbind(tab$lower[idx], tab$upper[idx])
    dimnames(out) <- list(paste0(lab, c("lower", "upper")), hyps)
    out
  }

  if (!is.null(x[["Cumulative_GORICA_weights"]])) {
    fcgw <- t(x[["Cumulative_GORICA_weights"]]["Final", ])
    rownames(fcgw) <- paste0(type_label, " weights")
    colnames(fcgw)  <- colnames(x[["Cumulative_GORICA_weights"]])
    # Prior IC weights (priorICweights; shown when unequal): row right before the IC weights
    prior_w <- x[["priorICweights"]]
    if (!is.null(prior_w) && length(prior_w) == ncol(fcgw) && length(unique(round(prior_w, 12))) > 1) {
      fpw <- matrix(prior_w, nrow = 1, dimnames = list(paste0("prior ", type_label, " weights"), colnames(fcgw)))
      final <- rbind(final, fpw)
    }
    final <- rbind(final, fcgw)
    final <- rbind(final, ci_rows("weights", colnames(fcgw)))
  }
  
  if (!is.null(x[["Cumulative_GORICA"]])) {
    fcgv <- t(x[["Cumulative_GORICA"]]["Final", ])
    rownames(fcgv) <- paste0(type_label, " values")
    final <- rbind(final, fcgv)
    final <- rbind(final, ci_rows("values", colnames(fcgv)))
  }
  
  if (!is.null(x[["Cumulative_LL_weights"]])) {
    fcgw <- t(x[["Cumulative_LL_weights"]]["Final", ])
    rownames(fcgw) <- "Log-likelihood weights"
    colnames(fcgw)  <- colnames(x[["Cumulative_LL_weights"]])
    final <- rbind(final, fcgw)
  }
  
  if (!is.null(x[["Cumulative_LL"]])) {
    fcllv <- t(x[["Cumulative_LL"]][ans$n_studies, ])
    rownames(fcllv) <- "Log-likelihood values"
    final <- rbind(final, fcllv)
  }
  
  if (!is.null(ans$Cumulative_PT)) { 
    fcptv <- t(ans$Cumulative_PT[ans$n_studies, ])
    rownames(fcptv) <- "Penalty term values"
    final <- rbind(final, fcptv)
  }
  
  # If icratios
  if (!is.null(x[["Cumulative_ratioICweights"]])) {
    f_cumICratio <- t(x[["Cumulative_ratioICweights"]]["Final", ])
    rownames(f_cumICratio) <- paste0("Ratio of ", type_label, " weights")
    #paste0("Ratio of ", type_label, " weights \n", "(versus reference hypothesis ", names(x[["Href"]]), ")")
    final <- rbind(final, f_cumICratio)
  }
  
  # If icratios
  if (!is.null(x[["Cumulative_ICdiff"]])) {
    f_cumICdiff <- t(x[["Cumulative_ICdiff"]]["Final", ])
    rownames(f_cumICdiff) <- paste0("Diff. in ", type_label, " values")
                #paste0("Difference in ", type_label, " values \n", "(versus reference hypothesis ", names(x[["Href"]]), ")")
    final <- rbind(final, f_cumICdiff)
  }
  
  ans$Final_Cumulative_results <- final
  ans$CI <- x[["CI"]]
  
  ans$priorICweights_study <- x[["priorICweights_study"]] # not active for now (see .evSyn_priorICweights_study)
  ans$Ratio_GORICA_weight_mu <- x[["ratio_GORICA_weight_mu"]]
  ans$Ratio_GORICA_weight_mc <- x[["ratio_GORICA_weight_mc"]]
  ans$Final_ratio_GORICA_weights <- x[["Final_ratio_GORICA_weights"]]
  ans$Final_ratio_LL_weights <- x[["Final_ratio_LL_weights"]]
  
  # Measures of variability over the studies (or analyses), incl. conclusion: 
  # only when requested (evSyn(..., variability = TRUE) or summary(..., variability = TRUE));
  # see evSyn_variability().
  ans$variability <- if (isTRUE(variability) && !is.null(x[["n_studies"]]) && x[["n_studies"]] > 1) {
    suppressMessages(evSyn_variability(object, probs = probs))
  } else NULL
  
  ans$messages <- list(mix_weights = x[["messages"]]$mix_weights)
  class(ans) <- "summary.evSyn"
  
  ans
  
  # TO DO x$messageAdded
  
}



## print.summary function
# Helper: note on the confidence intervals (CIs) of the final results
# (evSyn(..., CI = TRUE)), used in print.evSyn() and print.summary.evSyn().
.evSyn_CI_note <- function(ci, indent = "") {
  pct <- paste0(format(100 * ci$level), "%")
  meth <- switch(ci$method,
                 boot  = paste0("parametric bootstrap, B = ", ci$B,
                                " draws per study, independently over the studies"),
                 delta = "delta method",
                 mc    = "delta method, with Monte Carlo draws for the weights",
                 ci$method)
  if (isTRUE(ci$MI)) {
    meth <- "multiple imputation: Rubin's rules"
  } else if (identical(ci$type_ev, "average")) meth <- paste0(meth, "; dependent analyses, rho = ", format(ci$rho))
  if (identical(ci$dependence, "cluster")) meth <- paste0(meth, "; two-level synthesis, rho = ", format(ci$rho), " within clusters")
  if (identical(ci$dependence, "rve")) meth <- paste0(meth, "; robust variance estimation over ", ci$n_clusters, " clusters")
  paste0(indent, "(", pct, "CI: confidence intervals for the final (aggregated) ", toupper(ci$type),
         " values and weights, and, in the conclusion, for the final ratios of weights; computed with evSyn_CI() (",
         meth, "). Use evSyn_CI(object) for more details and other options; see ?evSyn_CI.)")
}

# Helper: print the measures of variability over the studies (or analyses) and the
# corresponding conclusion; used by print.summary.evSyn() and print.evSyn().
.evSyn_print_variability <- function(v, type_ev, type_label = "GORICA",
                                     digits = max(3, getOption("digits") - 4)) {
  if (is.null(v)) return(invisible(NULL))
  indentation <- "    "
  unit_cap <- .evSyn_unit_label(type_ev)
  unit_pl <- .evSyn_unit_label(type_ev, plural = TRUE, capital = FALSE)
  fmt_m <- function(m) apply(m, c(1, 2), function(val) format_numeric(val, digits = digits))
  pr_m <- function(m) {
    out <- capture.output(print(fmt_m(m), row.names = TRUE, right = TRUE, quote = FALSE))
    cat(paste0(indentation, out, "\n"), sep = "")
  }
  cat(paste0("\nVariability over the ", unit_pl, " (overall best hypothesis: ", v$best, "):\n"))
  cat(paste0("\n    ", unit_cap, "-specific ", type_label, " weights (percentiles):\n"))
  pr_m(v$weights_percentiles)
  cat("    ---\n")
  cat(paste0("\n    ", unit_cap, "-specific ratios of ", type_label, " weights (percentiles):\n"))
  pr_m(v$ratios_percentiles)
  cat("    ---\n")
  cat(paste0("\n    ", if (v$weighted) "Weighted geometric mean" else "Geometric mean", 
             ", evidence rate (ER), and stability rate (SR) of these ratios:\n"))
  pr_m(v$ratios_rates)
  if ("Ceiling" %in% colnames(v$ratios_rates)) {
    cat(paste0("\n    (Ceiling: number of ", unit_pl, " in which both hypotheses have the same log-likelihood value, ",
               "so that the ratio only reflects the difference in penalties.)\n"))
  }
  cat(paste0("\n    Preference rate (proportion of ", unit_pl, " in which ", v$best, 
             " has the highest ", type_label, " weight): ", 
             format_numeric(v$preference_rate, digits = digits), "\n"))
  cat(paste0("\n    Note: ER is the proportion of ", unit_pl, " whose ratio favours the same hypothesis ",
             "as the geometric mean (ER near 1: homogeneous ", unit_pl, "); \n",
             "    SR is the proportion of ", unit_pl, " whose ratio is more extreme than the geometric mean ",
             "(SR near .5: ratios evenly spread around it; \n",
             "    below .5: a few large ratios inflate it; above .5: a few small ratios deflate it); see Klaassen et al. (2018).\n"))
  if (v$weighted) {
    cat(paste0("    The geometric mean, ER, SR, and preference rate use the study weights; the percentiles do not.\n"))
  }
  cat("    ---\n")
  if (!is.null(v$conclusion)) {
    cat(paste0("\nConclusion (incl. variability over the ", unit_pl, "):\n"))
    # Wrap each line, keeping its indentation (bullets are continued below their text)
    for (ln in v$conclusion) {
      lead <- nchar(ln) - nchar(sub("^ +", "", ln))
      ext  <- lead + if (grepl("^ *(-|\\*) ", ln)) 2 else 0
      cat(strwrap(sub("^ +", "", ln), width = getOption("width", 80) - 2, 
                  indent = lead, exdent = ext), sep = "\n")
    }
  }
  invisible(NULL)
}

print.summary.evSyn <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  
  # added or equal approach
  type_ev <- x[["type_ev"]]
  # number of studies
  S <- x[["n_studies"]]
  
  # Determine label: default to "GORICA" unless type is explicitly given
  # Use [[ ]] to avoid partial matching (e.g., $type matching $type_ev)
  if (!is.null(x[["type"]])) {
    type_label <- switch(
      x$type,
      goric   = "GORIC",
      goricc  = "GORICC",
      gorica  = "GORICA",
      goricac = "GORICAC"
    )
  } else {
    type_label <- "GORICA"
  }
  
  cat(sprintf("restriktor (%s): ", packageDescription("restriktor", fields = "Version")))
  cat(paste(type_ev, "Evidence Synthesis results:\n"), sep = "")
  
  indentation <- "    "  # Four spaces for indentation
  
  # 'Study' or, when type_ev = "average" (multiple analyses of one study), 'Analysis'
  unit_cap <- .evSyn_unit_label(type_ev)
  
  cat(paste0("\n", unit_cap, "-specific results:\n"))
  
  if (!is.null(x[["GORICA_weight_m"]])) {
    cat(paste0("\n    ", type_label, " weights:\n"))
    formatted_gw <- apply(x[["GORICA_weight_m"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    input_gw <- cbind(1:S, formatted_gw)
    colnames(input_gw)[1] <- paste0("   ", unit_cap, " nr.")
    captured_output <- capture.output(print(input_gw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["GORICA_m"]])) {
    cat(paste0("\n    ", type_label, " values:\n"))
    formatted_gv <- apply(x[["GORICA_m"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    input_gv <- cbind(1:S, formatted_gv)
    colnames(input_gv)[1] <- paste0("   ", unit_cap, " nr.")
    captured_output <- capture.output(print(input_gv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["LL_weights_m"]])) {
    cat("\n    Log-likelihood weights:\n")  
    formatted_llw <- apply(x[["LL_weights_m"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    input_llw <- cbind(1:S, formatted_llw)
    colnames(input_llw)[1] <- paste0("   ", unit_cap, " nr.")
    captured_output <- capture.output(print(input_llw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["LL_m"]])) {
    cat("\n    Log-likelihood values:\n")  
    formatted_llv <- apply(x[["LL_m"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    input_llv <- cbind(1:S, formatted_llv)
    colnames(input_llv)[1] <- paste0("   ", unit_cap, " nr.")
    captured_output <- capture.output(print(input_llv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["PT_m"]])) {
    cat("\n    Penalty term values:\n")  
    formatted_ptv <- apply(x[["PT_m"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    input_ptv <- cbind(1:S, formatted_ptv)
    colnames(input_ptv)[1] <- paste0("   ", unit_cap, " nr.")
    captured_output <- capture.output(print(input_ptv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  # If icratios
  if (!is.null(x[["GwRatio_m"]])) {
    cat(paste0("\n    ", "Ratio of ", type_label, " weights",
               "\n    ", "(versus reference hypothesis ", names(x[["Href"]]), "):\n"))
    formatted_gw <- apply(x[["GwRatio_m"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    input_gw <- cbind(1:S, formatted_gw)
    colnames(input_gw)[1] <- paste0("   ", unit_cap, " nr.")
    captured_output <- capture.output(print(input_gw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  if (!is.null(x[["ICdiff_m"]])) {
    cat(paste0("\n    ", "Difference in ", type_label, " values",
               "\n    ", "(versus reference hypothesis ", names(x[["Href"]]), "):\n"))
    formatted_gv <- apply(x[["ICdiff_m"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    input_gv <- cbind(1:S, formatted_gv)
    colnames(input_gv)[1] <- paste0("   ", unit_cap, " nr.")
    captured_output <- capture.output(print(input_gv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }

  
  cat("\nCumulative results:\n")
  # The cumulative results are labelled by the position of the studies (i.e., 
  # 'Study nr.s 1-s'), which is the 'Study nr.' in the study-specific results 
  # above. If the studies were re-ordered (order_studies), state the order of 
  # the studies (by study name) such that it is clear which studies are included.
  ord <- x[["order_studies"]]
  if (!is.null(ord) && !identical(as.integer(ord), seq_len(S))) {
    st_names <- x[["study_names"]]
    if (is.null(st_names)) st_names <- as.character(ord)
    cat(paste0("    (", .evSyn_unit_label(type_ev, plural = TRUE), " are ordered as: ", 
               paste0(seq_len(S), " = ", st_names, collapse = ", "), ")\n"))
  }
  
  if (!is.null(x[["Cumulative_GORICA_weights"]])) {
    cat(paste0("\n    ", type_label, " weights:\n"))
    formatted_cgw <- apply(x[["Cumulative_GORICA_weights"]][1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cgw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["Cumulative_GORICA"]])) {
    cat(paste0("\n    ", type_label, " values:\n"))
    formatted_cgv <- apply(x[["Cumulative_GORICA"]][1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cgv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["Cumulative_LL_weights"]])) {
    cat("\n    Log-likelihood weights:\n")  
    formatted_cllw <- apply(x[["Cumulative_LL_weights"]][1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cllw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["LL_m"]])) {
    cat("\n    Log-likelihood values:\n")  
    formatted_cllv <- apply(x[["Cumulative_LL"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cllv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["PT_m"]])) {
    cat("\n    Penalty term values:\n")  
    formatted_cptv <- apply(x[["Cumulative_PT"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cptv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  # If icratios
  if (!is.null(x[["Cumulative_ratioICweights"]])) {
    cat(paste0("\n    ", "Ratio of ", type_label, " weights",
               "\n    ", "(versus reference hypothesis ", names(x[["Href"]]), "):\n"))
    formatted_crgw <- apply(x[["Cumulative_ratioICweights"]][1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_crgw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  # If icratios
  if (!is.null(x[["Cumulative_ICdiff"]])) {
    cat(paste0("\n    ", "Difference in ", type_label, " values",
               "\n    ", "(versus reference hypothesis ", names(x[["Href"]]), "):\n"))
    formatted_cgv <- apply(x[["Cumulative_ICdiff"]][1:S, , drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_cgv, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  
  cat("\nFinal results:\n")
  formatted_final <- apply(x[["Final_Cumulative_results"]][,,drop = FALSE], c(1,2), function(x) format_numeric(x, digits = digits))
  captured_output <- capture.output(print(formatted_final, row.names = TRUE, right = TRUE, quote = "FALSE"))
  adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
  cat(paste0(adjusted_output, "\n"), sep = "")
  if (any(grepl("^prior .* weights$", rownames(x[["Final_Cumulative_results"]])))) {
    cat(paste0("\n    (prior ", type_label, " weights: the priorICweights, which are used in the ", type_label,
               if (isTRUE(x[["priorICweights_study"]])) " weights (study-specific, cumulative, and final) and their ratios, but not in the " else
                           " weights (cumulative and final) and their ratios, but not in the study-specific ones (nor in the measures of variability based on these), the ",
               type_label, " values, the log-likelihood and penalty values, and the log-likelihood weights.)\n"))
  }
  if (!is.null(x[["CI"]])) {
    cat(paste0("\n", .evSyn_CI_note(x[["CI"]], indent = "    "), "\n"))
  }
  
  
  # The study weights are only used in calculating the IC weights; not in the 
  # IC values (which equal -2 x LL + penalty_factor x penalty), the log-likelihood 
  # (LL) and penalty values, and the LL weights.
  sw <- x[["study_weights"]]
  if (!is.null(sw) && length(unique(round(sw / sum(sw), 12))) > 1) {
    cat(paste0("\n    Note: The study weights are only used in calculating the ", type_label, 
               " weights; not in the ", type_label, " values (which equal -2 x LL + penalty_factor x penalty), \n",
               "    the log-likelihood (LL) and penalty values, and the LL weights. Hence, here, the ", 
               type_label, " weights are not based on the ", type_label, " values shown.\n"))
  }
  
  cat("\nFinal ratios:\n")
  
  if (!is.null(x[["Final_ratio_GORICA_weights"]])) {
    cat(paste0("\n    ", type_label, " weights:\n"))
    formatted_frgw <- apply(x[["Final_ratio_GORICA_weights"]], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_frgw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  if (!is.null(x[["Final_ratio_LL_weights"]])) {
    cat("\n    Log-likelihood weights:\n")  
    formatted_frllw <- apply(x[["Final_ratio_LL_weights"]], c(1,2), function(x) format_numeric(x, digits = digits))
    captured_output <- capture.output(print(formatted_frllw, row.names = TRUE, right = TRUE, quote = "FALSE"))
    adjusted_output <- gsub("^", indentation, captured_output, perl = TRUE)
    cat(paste0(adjusted_output, "\n"), sep = "")
    cat("    ---\n")
  }
  
  # Measures of variability over the studies (or analyses), incl. conclusion
  .evSyn_print_variability(x[["variability"]], type_ev = type_ev, type_label = type_label, digits = digits)
  
  cat("\n")
  message(x[["messages"]]$mix_weights)
  
  # # TO DO
  # if (!is.null(x[["messageAdded"]])) {
  #   cat("\n")
  #   message(x[["messageAdded"]])
  # }
  
}
