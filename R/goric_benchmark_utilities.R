# used in print.benchmark()
capitalize_first_letter <- function(input_string) {
  paste0(toupper(substring(input_string, 1, 1)), substring(input_string, 2))
}

# [CHANGE 2026-10 | audit] N7: new remove_self_row()/remove_self_col(): drop the preferred hypothesis' self-comparison by name (replaces remove_single_value_rows/col, which dropped by value)
# used in get_results_benchmark(): drop the preferred hypothesis'
# self-comparison (ratio = 1, log-ratio/difference = 0 by construction) BY
# NAME -- the row 'self_row' of a benchmark table, or the column 'self_col'
# of a matrix of draws. Previously rows/columns were dropped by VALUE (all
# values equal to 1 or 0), which was decided separately for the tables
# (Sample + quantiles) and the draws (every draw), and separately per
# population, so that a genuine comparison that happened to be (almost)
# constant (e.g. H1 vs. unconstrained when H1 is nearly always satisfied)
# could be dropped from one but not from the other -- leaving tables, draws,
# hypothesis_rate, rate_rlw and overlap misaligned.
remove_self_row <- function(data, self_row) {
  if (is.null(data)) return(data)
  data[rownames(data) != self_row, , drop = FALSE]
}

remove_self_col <- function(data, self_col) {
  if (is.null(data) || is.null(colnames(data))) return(data)
  data[, colnames(data) != self_col, drop = FALSE]
# [/CHANGE 2026-10]
}

# Function to filter columns based on exact matching of hypothesis_comparison
filter_columns <- function(data, hypothesis_comparison) {
  colnames_to_keep <- sapply(hypothesis_comparison, function(hypo) {
    grep(paste0("\\b", hypo, "\\b"), names(data), value = TRUE, fixed = FALSE)
  })
  colnames_to_keep <- unique(unlist(colnames_to_keep))
  colnames_to_keep <- intersect(names(data), colnames_to_keep) # Preserve original order
  return(data[, colnames_to_keep, drop = FALSE])
}

filter_vector <- function(data, hypothesis_comparison) {
  names_to_keep <- sapply(hypothesis_comparison, function(hypo) {
    grep(paste0("\\b", hypo, "\\b"), names(data), value = TRUE, fixed = FALSE)
  })
  names_to_keep <- unique(unlist(names_to_keep))
  names_to_keep <- intersect(names(data), names_to_keep) # Preserve original order
  return(data[names_to_keep])
}


# Functie om de power te berekenen voor een enkele density en critical value
calculate_power <- function(density_h1, critical_value) {
  sum(density_h1$y[density_h1$x > critical_value]) * mean(diff(density_h1$x))
}

# [CHANGE 2026-10 | Rebecca] new function compute_overlap(): overlapping coefficient of two benchmark distributions (kernel densities on a common grid)
# Overlapping coefficient (OVL) between two samples' distributions: the
# proportion of area shared by their kernel density estimates -- 1 = fully
# overlapping (identical) distributions, 0 = no overlap at all. This is the
# numeric counterpart of what overlaying, say, the 'Observed' and 'No-effect'
# population's benchmark draws in a density/histogram plot shows visually.
# Both densities are evaluated on the SAME x-grid (spanning the combined
# range of both samples, with a small margin so neither tail gets clipped)
# so their y-values are directly comparable point-by-point; Gaussian kernel
# with bw = "nrd0" matches the density estimate already used elsewhere in
# this file (see calculate_power() above). Like a plain density()-based plot,
# this does not respect the [0, 1] bounds that 'gw'/'lw' are constrained to
# -- a Gaussian kernel can put a little mass just outside that range, same
# as it would visually spill past the axis limits in a density plot.
# [CHANGE 2026-10 | audit] E5/O10: note on infinite draws; the overlap is NA (with a note) when a sample contains non-finite draws
#
# NOTE on infinite draws: a ratio of weights (rgw/rlw) can be Inf when the
# alternative's weight underflows to 0 (the ratios are formed on the log
# scale -- see compute_log_ratios() -- so this only happens when the
# log-ratio itself exceeds ~709, i.e. for truly overwhelming evidence). The
# quantile/percentile/rate computations in get_results_benchmark() handle
# Inf natively (an infinite draw simply counts as larger than any finite
# draw), so such draws are NOT dropped there. A kernel density, however,
# cannot be fitted on an infinite value. Such draws are NOT silently dropped
# here either: dropping them would give the overlap of the distributions
# CONDITIONAL on the draws being finite (e.g. 0.98 when only 10% of one
# sample is finite and happens to coincide with the other sample), which is
# not the overlap of the benchmark distributions. Instead, the overlap is
# NA (with an attribute "note" saying why) whenever either sample contains
# a non-finite draw; the log-ratio output_types (rgw_log/rlw_log) are always
# finite, so their overlap is always available and is the quantity to look
# at in such cases.
# [/CHANGE 2026-10]
compute_overlap <- function(draws1, draws2, n = 512) {
  # [CHANGE 2026-10 | audit] E5/O10: NA (with note) for non-finite draws instead of a crash on Inf
  if (any(!is.finite(draws1)) || any(!is.finite(draws2))) {
    return(structure(NA_real_, note = overlap_note_non_finite))
  }
  if (length(draws1) < 2 || length(draws2) < 2 ||
      sd(draws1) == 0 || sd(draws2) == 0) {
    # Can't fit a (non-degenerate) density on too few draws, or on draws
    # that are all identical (zero variance) -- not enough information for
    # an overlap coefficient to be meaningful.
    # [CHANGE 2026-10 | audit] E5/O10: NA with note (degenerate draws)
    return(structure(NA_real_, note = overlap_note_degenerate))
  }
  rng <- range(c(draws1, draws2))
  pad <- diff(rng) * 0.1
  from <- rng[1] - pad
  to   <- rng[2] + pad
  d1 <- tryCatch(
    density(draws1, kernel = "gaussian", bw = "nrd0", from = from, to = to, n = n),
    error = function(e) NULL
  )
  d2 <- tryCatch(
    density(draws2, kernel = "gaussian", bw = "nrd0", from = from, to = to, n = n),
    error = function(e) NULL
  )
  # [CHANGE 2026-10 | audit] E5/O10: NA with note when density() fails
  if (is.null(d1) || is.null(d2)) return(structure(NA_real_, note = overlap_note_degenerate))
  dx <- mean(diff(d1$x))
  min(sum(pmin(d1$y, d2$y)) * dx, 1) # cap at 1 for the (rare) numerical-integration overshoot
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] E5/O10: reasons for NA overlaps (notes) and new collect_overlap_notes()
# The reasons compute_overlap() can return NA, as stored in its "note"
# attribute and propagated (as attribute "notes" -- named by population, or
# by hypothesis column for the matrix-valued statistics) by
# compute_overlap_vs_observed()/compute_overlap_vs_observed_matrix() to
# print.benchmark(), where the NA is shown as "NA (<note>)" rather than as a
# number or a blank.
overlap_note_non_finite <- "non-finite draws; see rgw_log/rlw_log"
overlap_note_degenerate <- "too few or constant draws"

# Collects the "note" attributes of a list of compute_overlap() results into
# one named character vector (names = names(lst)), or NULL if none has one.
collect_overlap_notes <- function(lst) {
  notes <- vapply(lst, function(v) {
    nt <- attr(v, "note")
    if (is.null(nt)) NA_character_ else nt
  }, character(1))
  notes <- notes[!is.na(notes)]
  if (length(notes) == 0) NULL else notes
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | Rebecca] new function determine_overlap_reference(): reference population for the overlap column ("Observed", else the last population)
# Which population the "Overlap with ..." column is computed against.
# Prefers the "Observed" category when one is present (the usual case: the
# default 'No-effect' plus the user's actually-observed estimates). When
# there is no "Observed" category -- a custom pop_es/pop_est run, which
# never includes one -- falls back to the LAST population in
# 'pop_names' (i.e. the last pop_es/pop_est value supplied), so overlap is
# still shown rather than silently omitted. Returns NULL only when there are
# no populations at all.
determine_overlap_reference <- function(pop_names) {
  observed_name <- pop_names[grepl("= Observed$", pop_names)]
  if (length(observed_name) == 1) {
    return(observed_name)
  }
  if (length(pop_names) == 0) {
    return(NULL)
  }
  pop_names[length(pop_names)]
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | Rebecca] new function compute_overlap_vs_observed(): overlap of the reference population with each other population (gw/lw)
# Overlap of the reference population's (see determine_overlap_reference())
# benchmark distribution against EACH other population present -- by default
# just 'No-effect', but there can be more than one if the user supplied
# multiple pop_es/pop_est values. 'draws_combined' is a named list of draw
# vectors keyed by population (e.g. gw_combined/lw_combined from
# get_results_benchmark()), with names like "pop_es = No-effect"/
# "pop_es = Observed" (or "pop_est = ..." for benchmark_asymp). Returns a
# named numeric vector (named by the other population, plus the reference
# population's own entry fixed at 1 -- its overlap with itself), or NULL if
# there are no populations at all.
compute_overlap_vs_observed <- function(draws_combined) {
  reference_name <- determine_overlap_reference(names(draws_combined))
  if (is.null(reference_name)) {
    return(NULL)
  }
  reference_draws <- draws_combined[[reference_name]]
  other_names <- setdiff(names(draws_combined), reference_name)
  # [CHANGE 2026-10 | audit] E5/O10: keep the per-population overlap results to collect the notes of NA overlaps
  overlaps_lst <- lapply(other_names, function(nm) {
    compute_overlap(reference_draws, draws_combined[[nm]])
  })
  names(overlaps_lst) <- other_names
  overlaps <- vapply(overlaps_lst, as.numeric, numeric(1))
  # [/CHANGE 2026-10]
  names(overlaps) <- other_names
  overlaps[reference_name] <- 1
  # [CHANGE 2026-10 | audit] E5/O10: notes of NA overlaps as attribute
  # why an overlap is NA (if any), named by population -- see compute_overlap()
  attr(overlaps, "notes") <- collect_overlap_notes(overlaps_lst)
  overlaps
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | Rebecca] new function compute_overlap_vs_observed_matrix(): overlap per hypothesis column (rgw/rlw/ld)
# Same idea as compute_overlap_vs_observed(), but for the matrix-valued
# per-draw statistics ('rgw'/'rlw'/'ld'): each population's entry is a
# matrix with one column per alternative hypothesis being compared against
# the preferred one (e.g. column "H2"), rather than a single vector. Overlap
# is computed per shared column, so the result is a list (named by the other
# population, plus the reference population's own entry -- a same-length
# vector of 1s, one per column) of named numeric vectors (named by
# hypothesis column). Returns NULL under the same conditions as
# compute_overlap_vs_observed() (no populations at all), or when there are
# [CHANGE 2026-10 | audit] N7: comment (no columns when no successful draws)
# no columns to compare (e.g. no successful draws at all).
compute_overlap_vs_observed_matrix <- function(draws_combined) {
  reference_name <- determine_overlap_reference(names(draws_combined))
  if (is.null(reference_name)) {
    return(NULL)
  }
  reference_mat <- draws_combined[[reference_name]]
  other_names <- setdiff(names(draws_combined), reference_name)
  cols <- colnames(reference_mat)
  if (is.null(cols) || length(cols) == 0) {
    return(NULL)
  }
  overlaps <- lapply(other_names, function(nm) {
    other_mat <- draws_combined[[nm]]
    shared_cols <- intersect(cols, colnames(other_mat))
    # [CHANGE 2026-10 | audit] E5/O10: keep the per-column overlap results and their notes of NA overlaps
    ov_lst <- lapply(shared_cols, function(cn) {
      compute_overlap(reference_mat[, cn], other_mat[, cn])
    })
    names(ov_lst) <- shared_cols
    ov <- vapply(ov_lst, as.numeric, numeric(1))
    names(ov) <- shared_cols
    # why an overlap is NA (if any), named by hypothesis column -- see
    # compute_overlap()
    attr(ov, "notes") <- collect_overlap_notes(ov_lst)
    ov
    # [/CHANGE 2026-10]
  })
  names(overlaps) <- other_names
  self_overlap <- rep(1, length(cols))
  names(self_overlap) <- cols
  overlaps[[reference_name]] <- self_overlap
  overlaps
}
# [/CHANGE 2026-10]

# Function to calculate density
# calculate_density <- function(data, var, sample_value) {
#   dens <- density(data[[var]], kernel = "gaussian", na.rm = TRUE, bw = "nrd0")
#   data.frame(x = dens$x, y = dens$y, sample_value = sample_value)
# }


# Extract the parts within the parentheses
extract_in_parentheses <- function(name) {
  regmatches(name, regexpr("\\(.*\\)", name))
}

remove_spaces_in_parentheses <- function(string) {
  string <- gsub("\\(\\s+", "(", string)
  string <- gsub("\\s+\\)", ")", string)
  return(string)
}

construct_colnames <- function(list_name, colnames, pref_hypo_name) {
  remove_spaces_in_parentheses(paste0(list_name, " (", pref_hypo_name, " ", colnames, ")", sep = ""))
}

combine_matrices_cbind <- function(lst) {
  df_list <- lapply(lst, as.data.frame)
  combined_df <- do.call(cbind, df_list)
  return(combined_df)
}


combine_matrices_cbind <- function(lst) {
  # Vind het maximum aantal rijen in de lijst
  max_rows <- max(sapply(lst, nrow))
  
  # Zet de matrices om naar data frames en vul met NA waar nodig
  df_list <- lapply(lst, function(x) {
    as.data.frame(rbind(x, matrix(NA, nrow = max_rows - nrow(x), ncol = ncol(x))))
  })
  
  # Combineer de data frames kolomgewijs
  combined_df <- do.call(cbind, df_list)
  
  return(combined_df)
}


# model_type = "means" ----------------------------------------------------
detect_intercept <- function(model) {
  coefficients <- model$b.unrestr
  intercept_names <- c("(Intercept)", "Intercept", "(const)", "const",
                       "(Int)", "Int", "(Cons)", "Cons", "b0", "beta0")

  names_lower <- tolower(names(coefficients))
  intercept_names_lower <- tolower(intercept_names)

  detected_intercepts <- intercept_names[names_lower %in% intercept_names_lower]

  if (length(detected_intercepts) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# [CHANGE 2026-10 | audit] A12/E1: Cohen's f from a common residual variance sigma2 (was reconstructed from VCOV * (N - 1), which under-estimated sigma2)
# Compute Cohen's f based on the group means, the group sizes N and the
# (common) within-group / residual error variance sigma2:
#   f = sqrt( sum_g N_g (mu_g - mu)^2 / sum(N) ) / sigma,
# with mu the N-weighted grand mean (i.e. sqrt(SS_between / SS_within) with
# SS_within = sum(N) * sigma2). See residual_variance_from_vcov() for how
# sigma2 is obtained when no fitted model is available.
# (Previously sigma2 was reconstructed from the covariance matrix of the
# group means as VCOV * (N - 1); since Var(mean_g) = sigma2 / N_g, that
# under-estimated sigma2 by a factor (N_g - 1) / N_g, so that f was too
# large and generate_scaled_means() hit the wrong target f.)
compute_cohens_f <- function(group_means, N, sigma2) {
  total_mean <- sum(group_means * N) / sum(N)
  ss_between <- sum(N * (group_means - total_mean)^2)
  ss_within <- sum(N) * sigma2 # equates: summing sigma2 over i = 1 to N
# [/CHANGE 2026-10]
  # [CHANGE 2026-10 | Rebecca] Cohen's f as sqrt(SS_between / SS_within)
  cohens_f <- sqrt(ss_between/ss_within)

  return(cohens_f)
}


# [CHANGE 2026-10 | audit] A12: new residual_variance_from_vcov(): sigma2 = mean(N_g * VCOV[g, g]) for est + VCOV input, with a warning if the per-group values differ
# Residual (within-group) error variance sigma2 from the covariance matrix of
# the group means, for when no fitted model is available (input is est +
# VCOV; with a fitted lm model, benchmark_means() uses RSS/N = 
# deviance(fit)/nobs(fit) directly). ASSUMPTION: the group means are independent sample means of
# groups with a common error variance, so that Var(mean_g) = sigma2 / N_g,
# i.e. sigma2 = N_g * VCOV[g, g] for every group g; sigma2 is the average of
# these per-group values. If they differ by more than 'tol' (relative to
# their mean), the assumption is apparently violated (e.g. unequal error
# variances across groups, or estimates that are not plain group means) and
# a warning is given, since Cohen's f (and thus the population means for a
# given pop_es) is then only approximate.
residual_variance_from_vcov <- function(N, VCOV, tol = 0.1) {
  per_group <- N * diag(as.matrix(VCOV))
  sigma2 <- mean(per_group)
  rel_spread <- (max(per_group) - min(per_group)) / sigma2
  if (is.finite(rel_spread) && rel_spread > tol) {
    warning("\nrestriktor WARNING: The residual error variance (needed for Cohen's f) ",
            "is derived from the covariance matrix of the estimates as N_g * VCOV[g, g], ",
            "assuming independent group means with a common error variance. These ",
            "per-group values differ considerably (", paste(signif(per_group, 4), collapse = ", "),
            "), so this assumption seems violated; the (observed and population) ",
            "Cohen's f values are only approximate. Their average (", signif(sigma2, 4),
            ") is used.", call. = FALSE)
  }
  sigma2
}
# [/CHANGE 2026-10]


# Compute ratio data based on group_means
# compute_ratio_data <- function(group_means) {
#   # ratio_data <- rep(NA, ngroups)
#   # ratio_data[order(group_means) == 1] <- 1
#   # ratio_data[order(group_means) == 2] <- 2
#   # The choice of the smallest and the second smallest mean makes the scaling 
#   # more robust against changes in the other group means. Since these values 
#   # represent the lower bound of the data, the scale is less sensitive to the 
#   # spread of higher values.
#     
#   # For example:
#   # The value of 2.28 indicates that this particular group mean is 2.28 times 
#   # the scale factor d above the smallest mean. This means the group mean is 
#   # further from the smallest mean compared to the second smallest mean, and 
#   # helps in understanding the relative differences between the group means in 
#   # a normalized manner.
#   
#   # Aantal groepen
#   ngroups <- length(group_means)
#   # Lege vector voor ratio_data
#   ratio_data <- rep(NA, ngroups)
#   # Sorteer de indices van de waarden
#   sorted_indices <- order(group_means)
#   # Wijs 1 en 2 toe aan de kleinste en tweede kleinste waarden
#   ratio_data[sorted_indices[1]] <- 1
#   ratio_data[sorted_indices[2]] <- 2
#   # Bereken d: verschil tussen de tweede kleinste en de kleinste waarde
#   d <- group_means[sorted_indices[2]] - group_means[sorted_indices[1]]
#   # Bereken de ratio's voor de overige waarden
#   for (i in seq_len(ngroups)) {
#     if (!(i %in% sorted_indices[1:2])) {
#       ratio_data[i] <- 1 + (group_means[i] - group_means[sorted_indices[1]]) / d
#     }
#   }
#   return(ratio_data)
# }


# [CHANGE 2026-10 | audit] E1: generate_scaled_means() keeps the pattern of the means (centered deviations x positive factor; was division by min(), which flipped the ordering); clear error for all-equal means
# Population means with Cohen's f equal to 'target_f', keeping the PATTERN of
# 'group_means' (the observed group means, or a user-specified
# ratio_pop_means): the deviations from the (N-weighted) grand mean -- the
# quantity Cohen's f is based on, see compute_cohens_f() -- are multiplied by
# a single positive factor d, so that the ordering of the means, the sign of
# each deviation and the relative differences between the means are all
# preserved, whatever the sign of the means themselves. (Previously the means
# were divided by min(group_means) and rescaled by a positive factor: with a
# negative minimum this FLIPPED the ordering of the means, and with a zero
# minimum it gave NaN/Inf.) Since Cohen's f is shift-invariant and scales
# linearly with d, d = target_f / f(group_means) is available in closed form.
# The returned means are centered (N-weighted grand mean 0), so a pattern
# only matters up to a common shift: c(3, 2, 1) and c(1, 0, -1) give the same
# population means. The 'Observed' population in benchmark_means() uses the
# observed estimates as-is, rather than going through this function.
generate_scaled_means <- function(group_means, target_f, N, sigma2) {
  if (target_f == 0) {
    # If targeted effect size f is 0, then set all the means to 0.
    new_means <- rep(0, length(group_means))
  } else {
    grand_mean <- sum(group_means * N) / sum(N)
    deviations <- group_means - grand_mean
    f_pattern <- compute_cohens_f(group_means, N, sigma2)
    if (!is.finite(f_pattern) || f_pattern <= 0 ||
        isTRUE(all.equal(unname(deviations), rep(0, length(deviations))))) {
      stop("\nrestriktor ERROR: The group means (", paste(group_means, collapse = ", "),
           ") are all equal, so they cannot be scaled to the population effect size ",
           "(Cohen's f) ", target_f, ": the relative differences between the group ",
           "means are needed to do so. Specify a pattern of population means with ",
           "distinct values via the argument 'ratio_pop_means' (e.g., ratio_pop_means ",
           "= c(1, 2, 3)).", call. = FALSE)
    }
    d <- target_f / f_pattern
    new_means <- deviations * d
  }
  names(new_means) <- names(group_means)
# [/CHANGE 2026-10]

  return(new_means)
}

# Compute the population means based on the input parameters
# compute_population_means <- function(pop_es, ratio_pop_means, var_e, ngroups) {
#   means_pop_all <- matrix(NA, ncol = ngroups, nrow = length(pop_es))
#   nr_es <- length(pop_es)
#   for (teller_es in seq_len(nr_es)) {
#     #teller_es = 1
#     
#     # Determine mean values, with ratio of ratio.m
#     # Solve for x here
#     
#     # If all equal, then set population means to all 0
#     if (length(unique(ratio_pop_means)) == 1) {
#       means_pop <- rep(0, ngroups)
#     } else {
#       fun <- function(d) {
#         means_pop = ratio_pop_means * d 
#         (1/sqrt(var_e)) * sqrt((1/ngroups) * sum((means_pop - mean(means_pop))^2)) - pop_es[teller_es] #  AANPASSSEN NAAR NIEUWE FORMULE
#       }
#       d <- uniroot(fun, lower = 0, upper = 100)$root
#       # Construct means_pop
#       means_pop <- ratio_pop_means*d
#     }
#     means_pop_all[teller_es, ] <- means_pop
#   }
#   return(means_pop_all)  
# }

# [CHANGE 2026-10 | audit] removed: parallel_function_means() (dead code, B11/X8); benchmark_means() uses run_benchmark_simulation() with parallel_function_asymp()
# [CHANGE 2026-10 | audit] E5/E11: new compute_log_ratios()/extract_draw_results(): ratios of weights formed on the log scale (finite) and per-draw statistics for the workers
# Log of the ratio of the preferred hypothesis' GORIC(A) weight (rgw) and
# log-likelihood weight (rlw) to that of every hypothesis in a goric object,
# computed directly from the IC and log-likelihood values (and the prior IC
# weights) -- i.e., on the LOG scale -- rather than as log(ratio.gw): the
# weights themselves (and hence ratio.gw/ratio.lw from goric()) underflow to
# 0 / overflow to Inf for very large IC differences, whereas the log-ratio
# is finite whenever the IC values are finite. The ratio itself is exp() of
# this, so it can still be Inf (for a log-ratio above ~709) -- see
# compute_overlap() for how such draws are handled. Names match the
# columns of goric()$ratio.gw ("vs. <hypothesis>").
compute_log_ratios <- function(result, type, priorICweights, pref_hypo) {
  IC <- result[[type]]
  loglik <- result$loglik
  if (is.null(priorICweights)) priorICweights <- rep(1, length(IC))
  # weight_i = prior_i * exp(-IC_i/2) / sum(...), so
  # log(w_pref / w_j) = (IC_j - IC_pref)/2 + log(prior_pref) - log(prior_j);
  # likewise loglik-weight_i = exp(loglik_i) / sum(...).
  rgw_log <- 0.5 * (IC - IC[pref_hypo]) + log(priorICweights[pref_hypo]) - log(priorICweights)
  rlw_log <- loglik[pref_hypo] - loglik
  rgw_log[pref_hypo] <- 0 # self-comparison (ratio 1), also when prior = 0
  names(rgw_log) <- names(rlw_log) <- paste0("vs. ", result$model)
  list(rgw_log = rgw_log, rlw_log = rlw_log)
}

# The per-draw statistics returned by the workers below: gw/lw, the
# log-ratios (always finite, see compute_log_ratios()) and the ratios
# themselves (exp of the log-ratios).
extract_draw_results <- function(results_goric, pref_hypo) {
  log_ratios <- compute_log_ratios(results_goric$result, results_goric$type,
                                   results_goric$priorICweights, pref_hypo)
  ld <- results_goric$result$loglik[pref_hypo] - results_goric$result$loglik
  names(ld) <- names(log_ratios$rgw_log)
  list(
    gw  = results_goric$result[pref_hypo, 7], # goric(a) weight
    lw  = results_goric$result$loglik.weights[pref_hypo], # (unpenalized) log-likelihood weight
    rgw = exp(log_ratios$rgw_log), # ratio goric(a) weights
    rlw = exp(log_ratios$rlw_log), # ratio log-likelihood weights
    rgw_log = log_ratios$rgw_log,
    rlw_log = log_ratios$rlw_log,
    ld  = ld # loglik difference
  )
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] A9: new draw_failed(); worker parallel_function_asymp() muffles goric() warnings (draw kept, messages returned) and only drops a draw on an error; E2: same criterion/sample size, N3: priors of the object
# model_type = "asymp" ----------------------------------------------------

# A benchmark draw that failed: the worker (parallel_function_asymp())
# returns an empty list with attribute "error" (the error message) for a
# draw in which goric() gave an error; NULL is also treated as failed (for
# safety). Warnings inside a draw do NOT fail it -- see the worker.
draw_failed <- function(x) {
  is.null(x) || !is.null(attr(x, "error"))
}

# Worker for every benchmark draw (called from run_benchmark_simulation(),
# for both benchmark_means() and benchmark_asymp()): goric() on the i-th
# simulated estimate.
#
# NOTE on warnings and errors inside a draw: a warning given by goric()
# (e.g. a cosmetic one about the storage of ':=' definitions, or the
# non-convergence of the mix_weights bootstrap) does not invalidate the
# result -- goric() still returns a valid result -- so such warnings are
# muffled and the draw is KEPT; the warning messages are returned as
# attribute "warnings" and summarized (number of draws with a warning, first
# message) by run_benchmark_simulation(). Only an error invalidates a draw:
# the draw is then dropped (an empty list with attribute "error", see
# draw_failed()) and counted/reported likewise. (Previously every warning
# turned the draw into a dropped one -- so that with ':=' hypotheses, where
# every goric(est, VCOV) call warned, all draws were dropped and the
# benchmark crashed with "attempt to set an attribute on NULL" -- and only
# a per-draw console message remained as a trace.)
parallel_function_asymp <- function(i, est, VCOV, hypos, pref_hypo, comparison,
                                    type, control, mix_weights, penalty_factor,
                                    priorICweights = NULL, sample_nobs = NULL, ...) {
  warns <- character(0)
  results_goric <- tryCatch(
    withCallingHandlers(
      goric(est[i, ], VCOV = VCOV,
            hypotheses = hypos,
            comparison = comparison,
            # same criterion (gorica/goricac) and sample size as the
            # benchmarked object
            type = type,
            sample_nobs = sample_nobs,
            control = control,
            mix_weights = mix_weights,
            penalty_factor = penalty_factor,
            priorICweights = priorICweights,
            ...),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(list(), error = conditionMessage(e))
    }
  )

  if (draw_failed(results_goric)) {
    return(results_goric)
  }

  out <- extract_draw_results(results_goric, pref_hypo)
  if (length(warns) > 0) {
    attr(out, "warnings") <- warns
  }
  out
# [/CHANGE 2026-10]
}


# Define a function to extract and combine values from all elements in each pop_es list
extract_and_combine_values <- function(pop_es_list, value_name) {
  # [CHANGE 2026-10 | audit] A9: failed draws detected via draw_failed()
  failed <- vapply(pop_es_list, draw_failed, logical(1))
  empty_lists_count <- sum(failed)
  out <- do.call(rbind, lapply(pop_es_list[!failed], function(sub_list) sub_list[[value_name]]))
  attr(out, "empty_lists_count") <- empty_lists_count
  
  return(out)
}


# [CHANGE 2026-10 | audit] B31: new check_benchmark_weights(): clear error when the GORIC(A) weights are NA/NaN (e.g., goricac with infinite penalty)
# Stops with a clear message when the GORIC(A) weights of a goric object are
# NA/NaN (e.g. goricac with an infinite small-sample penalty when
# sample_nobs <= number of parameters + 1): such an object cannot be
# benchmarked (previously this surfaced much later as "attempt to select
# less than one element"). 'refit = TRUE' for the object refitted inside
# benchmark_means()/benchmark_asymp() (with alt_group_size/alt_sample_size
# and/or the goric -> gorica conversion).
check_benchmark_weights <- function(object, refit = FALSE) {
  w <- object$result[[paste0(object$type, ".weights")]]
  if (is.null(w)) w <- object$result[, ncol(object$result)]
  if (anyNA(w) || any(!is.finite(w))) {
    what <- if (refit) {
      paste0("the GORIC(A) object refitted for the benchmark (type '", object$type,
             "', sample size ", if (is.null(object$sample_nobs)) "NULL" else object$sample_nobs,
             ")")
    } else {
      paste0("the GORIC(A) object (type '", object$type, "')")
    }
    stop("\nrestriktor ERROR: The ", object$type, " weights of ", what, " are NA/NaN ",
         "(", paste(signif(w, 4), collapse = ", "), "), so it cannot be benchmarked. ",
         "This happens, for instance, for the goricac when the sample size is not larger ",
         "than the number of parameters + 1, so that its small-sample penalty is infinite; ",
         "check the goric object (and the sample size) first.", call. = FALSE)
  }
  invisible(TRUE)
}
# [/CHANGE 2026-10]


# [CHANGE 2026-10 | audit] B34: new unique_population_names(): duplicate population names made unique (was a crash)
# Population names (names of pop_es / rownames of pop_est) must be unique:
# the results are stored per population by name. Duplicates are made unique
# with make.unique(), with a message. (Previously duplicate names crashed
# with "invalid 'type' (list) of argument".)
unique_population_names <- function(rnames, what = "pop_es") {
  if (anyDuplicated(rnames)) {
    new_names <- make.unique(rnames)
    message("\nrestriktor Message: The population names of '", what, "' are not unique (",
            paste(rnames[duplicated(rnames)], collapse = ", "), "). They have been made ",
            "unique: ", paste(new_names, collapse = ", "), ".")
    rnames <- new_names
  }
  rnames
}
# [/CHANGE 2026-10]


## 
get_results_benchmark <- function(x, object, pref_hypo, pref_hypo_name,
                                  # [CHANGE 2026-10 | Rebecca] thresholds for hypothesis_rate and rate_rlw passed in
                                  quant, names_quant, nr.hypos,
                                  hypo_rate_threshold = 1,
                                  threshold_rlw = 1) {
  results <- x
  
    # Use lapply to apply the extract_and_combine_values function to each element in the results list
  gw_combined  <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "gw"))
  # [CHANGE 2026-10 | Rebecca] draws of the log-likelihood weights (lw)
  lw_combined  <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "lw"))
  rgw_combined <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "rgw"))
  rlw_combined <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "rlw"))
  ld_combined  <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "ld"))
  # [CHANGE 2026-10 | audit] E5/E11: draws of the log-ratios (always finite) and the Sample values of the (log-)ratios computed on the log scale
  # log(rgw)/log(rlw), as computed per draw on the log scale (see
  # compute_log_ratios()): always finite, also when rgw/rlw itself is Inf.
  rgw_log_combined <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "rgw_log"))
  rlw_log_combined <- lapply(results, function(pop_es_list) extract_and_combine_values(pop_es_list, "rlw_log"))
  # The 'Sample' values of the (log-)ratios, computed in the same way as the
  # draws (on the log scale; the ratio is exp() of the log-ratio).
  sample_log_ratios <- compute_log_ratios(object$result, object$type,
                                          object$priorICweights, pref_hypo)
  # [/CHANGE 2026-10]

  # [CHANGE 2026-10 | Rebecca] reference population for the "Pctl. median ref. pop." and overlap columns, determined once
  # Which population the "median ref. pop." columns below (and the "Overlap
  # with ..." column -- see overlap_reference_pop further down, which now
  # just reuses this) are computed against: "Observed" when present, else
  # the last pop_es/pop_est population supplied. Determined once here, since
  # it doesn't depend on gw/lw/rgw/etc. -- see determine_overlap_reference().
  reference_pop_name <- determine_overlap_reference(names(gw_combined))
  # [/CHANGE 2026-10]

  # Calculate CI_benchmarks_gw for each pop_es category
  CI_benchmarks_gw <- lapply(gw_combined, function(gw_values) {
    CI_benchmarks_gw <- matrix(c(object$result[pref_hypo, 7], quantile(gw_values, 
                                                                       quant, na.rm = TRUE)), 
                               nrow = 1)
    colnames(CI_benchmarks_gw) <- names_quant
    rownames(CI_benchmarks_gw) <- pref_hypo_name
    CI_benchmarks_gw

  })

  # [CHANGE 2026-10 | Rebecca] percentile of the Sample value (Pctl. Sample) and of the reference population's median (Pctl. median ref. pop.) for gw; new lw benchmarks with the same columns
  # Percentile of the observed goric(a) weight ('Sample' value) within its
  # own benchmark distribution, for each pop_es category (0-100 scale,
  # matching the others). Column/print header is "Pctl. Sample".
  pctl_Sample_gw <- lapply(gw_combined, function(gw_values) {
    Fn <- ecdf(gw_values)
    pctl_Sample_gw <- matrix(Fn(object$result[pref_hypo, 7]) * 100, nrow = 1)
    colnames(pctl_Sample_gw) <- "Pctl. Sample"
    rownames(pctl_Sample_gw) <- pref_hypo_name
    pctl_Sample_gw
  })

  # Percentile of the reference population's own median goric(a) weight,
  # located within *each* population's own benchmark distribution (0-100
  # scale): where the reference/last population's typical draw falls
  # relative to this population's distribution. For the reference
  # population's own row this is hardcoded to exactly 50 -- ecdf() of a
  # finite/even-n sample evaluated at its own median need not land on
  # exactly 50. Column/print header is "Pctl. median ref. pop.".
  median_gw_ref <- median(gw_combined[[reference_pop_name]], na.rm = TRUE)
  pctl_medianRefPop_gw <- lapply(names(gw_combined), function(nm) {
    value <- if (identical(nm, reference_pop_name)) {
      50
    } else {
      Fn <- ecdf(gw_combined[[nm]])
      Fn(median_gw_ref) * 100
    }
    out <- matrix(value, nrow = 1)
    colnames(out) <- "Pctl. median ref. pop."
    rownames(out) <- pref_hypo_name
    out
  })
  names(pctl_medianRefPop_gw) <- names(gw_combined)

  # Same as CI_benchmarks_gw/pctl_Sample_gw above, but for the (unpenalized)
  # log-likelihood weight -- used together with pctl_Sample_gw by the
  # iter-adequacy check (see goric_percentile_test()/check_iter_adequacy()),
  # since a user may ultimately be interested in either output_type.
  CI_benchmarks_lw <- lapply(lw_combined, function(lw_values) {
    CI_benchmarks_lw <- matrix(c(object$result$loglik.weights[pref_hypo], quantile(lw_values,
                                                                       quant, na.rm = TRUE)),
                               nrow = 1)
    colnames(CI_benchmarks_lw) <- names_quant
    rownames(CI_benchmarks_lw) <- pref_hypo_name
    CI_benchmarks_lw
  })

  pctl_Sample_lw <- lapply(lw_combined, function(lw_values) {
    Fn <- ecdf(lw_values)
    pctl_Sample_lw <- matrix(Fn(object$result$loglik.weights[pref_hypo]) * 100, nrow = 1)
    colnames(pctl_Sample_lw) <- "Pctl. Sample"
    rownames(pctl_Sample_lw) <- pref_hypo_name
    pctl_Sample_lw
  })

  median_lw_ref <- median(lw_combined[[reference_pop_name]], na.rm = TRUE)
  pctl_medianRefPop_lw <- lapply(names(lw_combined), function(nm) {
    value <- if (identical(nm, reference_pop_name)) {
      50
    } else {
      Fn <- ecdf(lw_combined[[nm]])
      Fn(median_lw_ref) * 100
    }
    out <- matrix(value, nrow = 1)
    colnames(out) <- "Pctl. median ref. pop."
    rownames(out) <- pref_hypo_name
    out
  })
  names(pctl_medianRefPop_lw) <- names(lw_combined)
  # [/CHANGE 2026-10]


  # Initialize matrices to store CI benchmarks for current pop_es category
  CI_benchmarks_rgw <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_rlw <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_rlw_ge1 <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  # [CHANGE 2026-10 | Rebecca] tables for the log-ratios rgw_log/rlw_log
  CI_benchmarks_rgw_log <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_rlw_log <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_ld <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))
  CI_benchmarks_ld_ge0 <- matrix(NA, nrow = nr.hypos, ncol = 1 + length(quant))

  # Fill the first column with sample values
  # [CHANGE 2026-10 | audit] E5/E11: Sample values of the ratios and log-ratios from compute_log_ratios() (exp of the log-ratio)
  CI_benchmarks_rgw[, 1] <- exp(sample_log_ratios$rgw_log)
  CI_benchmarks_rlw[, 1] <- exp(sample_log_ratios$rlw_log)
  CI_benchmarks_rgw_log[, 1] <- sample_log_ratios$rgw_log
  CI_benchmarks_rlw_log[, 1] <- sample_log_ratios$rlw_log
  # folded to >= 1: exp(|log-ratio|)
  CI_benchmarks_rlw_ge1[, 1] <- exp(abs(sample_log_ratios$rlw_log))
  # [/CHANGE 2026-10]
  CI_benchmarks_ld[, 1] <- object$result$loglik[pref_hypo] - object$result$loglik 
  CI_benchmarks_ld_ge0[, 1] <- abs(object$result$loglik[pref_hypo] - object$result$loglik) 
  
  CI_benchmarks_rgw_all <- list()
  CI_benchmarks_rlw_all <- list()
  CI_benchmarks_rlw_ge1_all <- list()
  # [CHANGE 2026-10 | Rebecca] lists for the rgw_log/rlw_log tables
  CI_benchmarks_rgw_log_all <- list()
  CI_benchmarks_rlw_log_all <- list()
  CI_benchmarks_ld_all <- list()
  CI_benchmarks_ld_ge0_all <- list()

  # [CHANGE 2026-10 | Rebecca] lists for the Pctl. Sample / Pctl. median ref. pop. results; reference population's medians per output type
  pctl_Sample_rgw_all <- list()
  pctl_Sample_rlw_all <- list()
  pctl_Sample_rlw_ge1_all <- list()
  pctl_Sample_rgw_log_all <- list()
  pctl_Sample_rlw_log_all <- list()
  pctl_Sample_ld_all <- list()
  pctl_Sample_ld_ge0_all <- list()

  medianRefPop_rgw_all <- list()
  medianRefPop_rlw_all <- list()
  medianRefPop_rlw_ge1_all <- list()
  medianRefPop_rgw_log_all <- list()
  medianRefPop_rlw_log_all <- list()
  medianRefPop_ld_all <- list()
  medianRefPop_ld_ge0_all <- list()

  # Reference population's own per-hypothesis median value, for each
  # output_type that the "Pctl. median ref. pop." column is computed for.
  # Fixed across all iterations of the loop below (it doesn't depend on
  # 'name'), so computed once here rather than per-population. rgw_log/
  # rlw_log need their own log() here since the rgw_log_combined/
  # rlw_log_combined lists aren't built until after this loop (see below).
  ref_median_rgw <- apply(rgw_combined[[reference_pop_name]], 2, median, na.rm = TRUE)
  ref_median_rlw <- apply(rlw_combined[[reference_pop_name]], 2, median, na.rm = TRUE)
  # [CHANGE 2026-10 | audit] E5: reference medians of the log-ratios from the log-scale draws
  ref_median_rgw_log <- apply(rgw_log_combined[[reference_pop_name]], 2, median, na.rm = TRUE)
  ref_median_rlw_log <- apply(rlw_log_combined[[reference_pop_name]], 2, median, na.rm = TRUE)
  ref_median_ld <- apply(ld_combined[[reference_pop_name]], 2, median, na.rm = TRUE)
  # rlw_ge1/ld_ge0 (see 'Prepare rlw_ge1 and ld_ge0 matrices' below) aren't
  # available as standalone combined lists the way rgw/rlw/ld are -- they're
  # folded per-population inside the loop below -- so the reference
  # population's own version is folded here too, just for this median.
  ref_rlw_ge1 <- rlw_combined[[reference_pop_name]]
  ref_rlw_ge1[ref_rlw_ge1 < 1] <- 1 / ref_rlw_ge1[ref_rlw_ge1 < 1]
  ref_median_rlw_ge1 <- apply(ref_rlw_ge1, 2, median, na.rm = TRUE)
  ref_median_ld_ge0 <- apply(abs(ld_combined[[reference_pop_name]]), 2, median, na.rm = TRUE)
  # [/CHANGE 2026-10]

  # Loop through each pop_es category to fill in the CI benchmark lists
  for (name in names(results)) {
    rgw_combined_values <- rgw_combined[[name]]
    rlw_combined_values  <- rlw_combined[[name]]
    ld_combined_values  <- ld_combined[[name]]
    
    # Prepare rlw_ge1 and ld_ge0 matrices
    rlw_ge1 <- rlw_combined_values
    rlw_ge1[rlw_combined_values < 1] <- 1 / rlw_combined_values[rlw_combined_values < 1]
    ld_ge0 <- abs(ld_combined_values)

    # [CHANGE 2026-10 | Rebecca] rgw_log/rlw_log: log of the ratios per draw (scale-invariant, symmetric around 0)
    # log(rgw)/log(rlw) -- a scale-invariant, symmetric-around-0 alternative
    # to rgw/rlw themselves: log(rgw)/log(rlw) equals the log-odds (logit) of
    # the corresponding pairwise-rescaled weight (e.g. lw_pref/(lw_pref+lw_k)),
    # since the shared (lw_pref+lw_k) normalizer cancels out of the ratio --
    # so it carries no new information beyond rgw/rlw, but is far better
    # suited to things like the 'Observed'-vs-'No Effect' overlap comparison:
    # unbounded and symmetric around 0 (rather than living on the heavily
    # right-skewed (0, Inf) scale that rgw/rlw do, where '1' isn't a natural
    # visual center), and not distorted by how much weight other hypotheses
    # in the set happen to be carrying (that dilution cancels out of the
    # ratio -- and so out of its log -- exactly, per draw). Self-comparison
    # (pref vs pref) is log(1) = 0 here, rather than the 1 that rgw/rlw use,
    # so it gets cleaned via the 0-baseline (like ld) rather than the
    # [CHANGE 2026-10 | audit] E5/E11: log-ratios taken from the log-scale draws (finite also when the ratio is Inf)
    # 1-baseline used for rgw/rlw below. Taken from the draws' own log-scale
    # computation (compute_log_ratios()) rather than as log() of the ratio,
    # so these are finite even when the ratio itself overflowed to Inf.
    rgw_log_combined_values <- rgw_log_combined[[name]]
    rlw_log_combined_values <- rlw_log_combined[[name]]
    # [/CHANGE 2026-10]
    # [/CHANGE 2026-10]
    
    # Loop through the hypotheses and calculate the quantiles
    for (j in seq_len(nr.hypos)) {
      CI_benchmarks_rgw[j, 2:(1 + length(quant))] <- quantile(rgw_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_rlw[j, 2:(1 + length(quant))] <- quantile(rlw_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_rlw_ge1[j, 2:(1 + length(quant))] <- quantile(rlw_ge1[, j], quant, na.rm = TRUE)
      # [CHANGE 2026-10 | Rebecca] quantiles of the log-ratios
      CI_benchmarks_rgw_log[j, 2:(1 + length(quant))] <- quantile(rgw_log_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_rlw_log[j, 2:(1 + length(quant))] <- quantile(rlw_log_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_ld[j, 2:(1 + length(quant))] <- quantile(ld_combined_values[, j], quant, na.rm = TRUE)
      CI_benchmarks_ld_ge0[j, 2:(1 + length(quant))] <- quantile(ld_ge0[, j], quant, na.rm = TRUE)
    }
    # [CHANGE 2026-10 | Rebecca] per hypothesis: percentile of the Sample value and of the reference population's median within this population's distribution (incl. rgw_log/rlw_log, rlw_ge1, ld_ge0); stored per population
    # Loop through the hypotheses and calculate the percentile of the sample
    # finding ('Sample' value; pctl_Sample_* below) as well as the percentile
    # of the reference population's own median value within this
    # population's own distribution ("Pctl. median ref. pop."; medianRefPop_*
    # below).
    percentile_rgw <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_rlw <- matrix(NA, nrow = nr.hypos, ncol = 1 )
    percentile_rlw_ge1 <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_rgw_log <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_rlw_log <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_ld <- matrix(NA, nrow = nr.hypos, ncol = 1)
    percentile_ld_ge0 <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_rgw <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_rlw <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_rlw_ge1 <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_rgw_log <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_rlw_log <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_ld <- matrix(NA, nrow = nr.hypos, ncol = 1)
    medianRefPop_ld_ge0 <- matrix(NA, nrow = nr.hypos, ncol = 1)
    # Hardcode this population's own row to exactly 50 when it IS the
    # reference population -- see the comment above ref_median_rgw etc.
    is_ref_pop <- identical(name, reference_pop_name)
    for (j in seq_len(nr.hypos)) {
      Fn <- ecdf(rgw_combined_values[, j])
      percentile_rgw[j, 1] <- Fn(CI_benchmarks_rgw[j, 1]) * 100
      medianRefPop_rgw[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_rgw[j]) * 100
      #
      Fn <- ecdf(rlw_combined_values[, j])
      percentile_rlw[j, 1] <- Fn(CI_benchmarks_rlw[j, 1]) * 100
      medianRefPop_rlw[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_rlw[j]) * 100
      #
      Fn <- ecdf(rlw_ge1[, j])
      percentile_rlw_ge1[j, 1] <- Fn(CI_benchmarks_rlw_ge1[j, 1]) * 100
      medianRefPop_rlw_ge1[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_rlw_ge1[j]) * 100
      #
      Fn <- ecdf(rgw_log_combined_values[, j])
      percentile_rgw_log[j, 1] <- Fn(CI_benchmarks_rgw_log[j, 1]) * 100
      medianRefPop_rgw_log[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_rgw_log[j]) * 100
      #
      Fn <- ecdf(rlw_log_combined_values[, j])
      percentile_rlw_log[j, 1] <- Fn(CI_benchmarks_rlw_log[j, 1]) * 100
      medianRefPop_rlw_log[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_rlw_log[j]) * 100
      #
      Fn <- ecdf(ld_combined_values[, j])
      percentile_ld[j, 1] <- Fn(CI_benchmarks_ld[j, 1]) * 100
      medianRefPop_ld[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_ld[j]) * 100
      #
      Fn <- ecdf(ld_ge0[, j])
      percentile_ld_ge0[j, 1] <- Fn(CI_benchmarks_ld_ge0[j, 1]) * 100
      medianRefPop_ld_ge0[j, 1] <- if (is_ref_pop) 50 else Fn(ref_median_ld_ge0[j]) * 100
    }

    # Label the percentiles (percentage of the benchmark distribution at or
    # below the sample value / reference-population median, i.e. on a 0-100
    # scale)
    percentile_names <- paste(pref_hypo_name, names(object$ratio.gw[pref_hypo, ]))
    rownames(percentile_rgw) <- rownames(percentile_rlw) <-
      rownames(percentile_rlw_ge1) <- rownames(percentile_rgw_log) <-
      rownames(percentile_rlw_log) <- rownames(percentile_ld) <-
      rownames(percentile_ld_ge0) <-
      rownames(medianRefPop_rgw) <- rownames(medianRefPop_rlw) <-
      rownames(medianRefPop_rlw_ge1) <-
      rownames(medianRefPop_rgw_log) <- rownames(medianRefPop_rlw_log) <-
      rownames(medianRefPop_ld) <- rownames(medianRefPop_ld_ge0) <- percentile_names
    # rlw_ge1/ld_ge0's "Sample" percentile columns now match the "Pctl.
    # Sample" naming used by every other output_type's Sample-percentile
    # column (was "percentile" -- an inconsistency in the original naming
    # that motivated renaming the percentile_rlw_ge1/percentile_absdifLL
    # output fields below in the first place).
    colnames(percentile_rgw) <- colnames(percentile_rlw) <-
      colnames(percentile_rlw_ge1) <-
      colnames(percentile_rgw_log) <- colnames(percentile_rlw_log) <-
      colnames(percentile_ld) <- colnames(percentile_ld_ge0) <- "Pctl. Sample"
    colnames(medianRefPop_rgw) <- colnames(medianRefPop_rlw) <-
      colnames(medianRefPop_rlw_ge1) <-
      colnames(medianRefPop_rgw_log) <- colnames(medianRefPop_rlw_log) <-
      colnames(medianRefPop_ld) <- colnames(medianRefPop_ld_ge0) <- "Pctl. median ref. pop."

    # Store this pop_es category's percentiles so they survive past this
    # loop iteration (mirrors the CI_benchmarks_*_all pattern below)
    pctl_Sample_rgw_all[[name]] <- percentile_rgw
    pctl_Sample_rlw_all[[name]] <- percentile_rlw
    pctl_Sample_rlw_ge1_all[[name]] <- percentile_rlw_ge1
    pctl_Sample_rgw_log_all[[name]] <- percentile_rgw_log
    pctl_Sample_rlw_log_all[[name]] <- percentile_rlw_log
    pctl_Sample_ld_all[[name]] <- percentile_ld
    pctl_Sample_ld_ge0_all[[name]] <- percentile_ld_ge0

    medianRefPop_rgw_all[[name]] <- medianRefPop_rgw
    medianRefPop_rlw_all[[name]] <- medianRefPop_rlw
    medianRefPop_rlw_ge1_all[[name]] <- medianRefPop_rlw_ge1
    medianRefPop_rgw_log_all[[name]] <- medianRefPop_rgw_log
    medianRefPop_rlw_log_all[[name]] <- medianRefPop_rlw_log
    medianRefPop_ld_all[[name]] <- medianRefPop_ld
    medianRefPop_ld_ge0_all[[name]] <- medianRefPop_ld_ge0
    # [/CHANGE 2026-10]

    # Set column names for the CI benchmarks
    colnames(CI_benchmarks_rgw) <- colnames(CI_benchmarks_rlw) <-
      # [CHANGE 2026-10 | Rebecca] column names incl. rgw_log/rlw_log tables
      colnames(CI_benchmarks_rlw_ge1) <- colnames(CI_benchmarks_rgw_log) <-
      colnames(CI_benchmarks_rlw_log) <- colnames(CI_benchmarks_ld) <-
      colnames(CI_benchmarks_ld_ge0) <- names_quant

    # Set row names for the CI benchmarks
    rownames(CI_benchmarks_rgw) <- rownames(CI_benchmarks_rlw) <-
      # [CHANGE 2026-10 | Rebecca] row names incl. rgw_log/rlw_log tables
      rownames(CI_benchmarks_rlw_ge1) <- rownames(CI_benchmarks_rgw_log) <-
      rownames(CI_benchmarks_rlw_log) <- rownames(CI_benchmarks_ld) <-
      rownames(CI_benchmarks_ld_ge0) <- paste(pref_hypo_name, names(object$ratio.gw[pref_hypo, ]))

    # Store CI benchmarks in lists
    CI_benchmarks_rgw_all[[name]] <- CI_benchmarks_rgw
    CI_benchmarks_rlw_all[[name]] <- CI_benchmarks_rlw
    CI_benchmarks_rlw_ge1_all[[name]] <- CI_benchmarks_rlw_ge1
    # [CHANGE 2026-10 | Rebecca] store the rgw_log/rlw_log tables
    CI_benchmarks_rgw_log_all[[name]] <- CI_benchmarks_rgw_log
    CI_benchmarks_rlw_log_all[[name]] <- CI_benchmarks_rlw_log
    CI_benchmarks_ld_all[[name]] <- CI_benchmarks_ld
    CI_benchmarks_ld_ge0_all[[name]] <- CI_benchmarks_ld_ge0
  }
  
  # [CHANGE 2026-10 | audit] N7: drop the self-comparison row by name (remove_self_row) from every table, in the same way for every output type
  # The preferred hypothesis' self-comparison (always 1 for the ratios, 0 for
  # the log-ratios and log-likelihood differences) is dropped -- by name, and
  # in exactly the same way for every output_type and population -- from
  # the benchmark tables (row 'self_row') and from the draws (column
  # 'self_col'), so that tables, draws, hypothesis_rate, rate_rlw and overlap
  # always refer to the same set of (alternative) hypotheses, in the same
  # order. See remove_self_row()/remove_self_col().
  self_col <- names(object$ratio.gw[pref_hypo, ])[pref_hypo]
  self_row <- paste(pref_hypo_name, self_col)

  CI_benchmarks_rgw_all_cleaned <- lapply(CI_benchmarks_rgw_all, function(pop_es_list) {
    remove_self_row(pop_es_list, self_row)
  })
  
  CI_benchmarks_rlw_all_cleaned <- lapply(CI_benchmarks_rlw_all, function(pop_es_list) {
    remove_self_row(pop_es_list, self_row)
  })
  
  CI_benchmarks_rlw_ge1_all_cleaned <- lapply(CI_benchmarks_rlw_ge1_all, function(pop_es_list) {
    remove_self_row(pop_es_list, self_row)
  })

  # [CHANGE 2026-10 | Rebecca] cleaned rgw_log table
  CI_benchmarks_rgw_log_all_cleaned <- lapply(CI_benchmarks_rgw_log_all, function(pop_es_list) {
    remove_self_row(pop_es_list, self_row)
  })

  # [CHANGE 2026-10 | Rebecca] cleaned rlw_log table
  CI_benchmarks_rlw_log_all_cleaned <- lapply(CI_benchmarks_rlw_log_all, function(pop_es_list) {
    remove_self_row(pop_es_list, self_row)
  })

  CI_benchmarks_ld_all_cleaned <- lapply(CI_benchmarks_ld_all, function(pop_es_list) {
    remove_self_row(pop_es_list, self_row)
  })
  
  CI_benchmarks_ld_ge0_all_cleaned <- lapply(CI_benchmarks_ld_ge0_all, function(pop_es_list) {
    remove_self_row(pop_es_list, self_row)
  })
  # [/CHANGE 2026-10]

  # [CHANGE 2026-10 | Rebecca] align the Pctl. Sample / Pctl. median ref. pop. matrices with the cleaned tables by row name (align_rows)
  # [CHANGE 2026-10 | audit] N7: comment (remove_self_row)
  # remove_self_row() drops the preferred hypothesis' self-comparison
  # row from the benchmarks_* matrices above (ratio/difference vs. itself is
  # always 1 or 0). The percentile_*_all matrices were never subject to that
  # filter (a percentile is essentially never exactly 1 or 0), so without this
  # they'd have one row more than their benchmarks_* counterpart, and cbind()-ing
  # them together (e.g. in print/summary) would silently misalign rows instead
  # of erroring. Subset by rowname to keep both sets of matrices in lockstep.
  align_rows <- function(percentile_list, cleaned_list) {
    Map(function(perc, bench) perc[rownames(bench), , drop = FALSE],
        percentile_list, cleaned_list)
  }
  pctl_Sample_rgw_all_cleaned <- align_rows(pctl_Sample_rgw_all, CI_benchmarks_rgw_all_cleaned)
  pctl_Sample_rlw_all_cleaned <- align_rows(pctl_Sample_rlw_all, CI_benchmarks_rlw_all_cleaned)
  pctl_Sample_rlw_ge1_all_cleaned <- align_rows(pctl_Sample_rlw_ge1_all, CI_benchmarks_rlw_ge1_all_cleaned)
  pctl_Sample_rgw_log_all_cleaned <- align_rows(pctl_Sample_rgw_log_all, CI_benchmarks_rgw_log_all_cleaned)
  pctl_Sample_rlw_log_all_cleaned <- align_rows(pctl_Sample_rlw_log_all, CI_benchmarks_rlw_log_all_cleaned)
  pctl_Sample_ld_all_cleaned <- align_rows(pctl_Sample_ld_all, CI_benchmarks_ld_all_cleaned)
  pctl_Sample_ld_ge0_all_cleaned <- align_rows(pctl_Sample_ld_ge0_all, CI_benchmarks_ld_ge0_all_cleaned)

  # Same alignment treatment for the new "Pctl. median ref. pop." matrices.
  pctl_medianRefPop_rgw_all_cleaned <- align_rows(medianRefPop_rgw_all, CI_benchmarks_rgw_all_cleaned)
  pctl_medianRefPop_rlw_all_cleaned <- align_rows(medianRefPop_rlw_all, CI_benchmarks_rlw_all_cleaned)
  pctl_medianRefPop_rlw_ge1_all_cleaned <- align_rows(medianRefPop_rlw_ge1_all, CI_benchmarks_rlw_ge1_all_cleaned)
  pctl_medianRefPop_rgw_log_all_cleaned <- align_rows(medianRefPop_rgw_log_all, CI_benchmarks_rgw_log_all_cleaned)
  pctl_medianRefPop_rlw_log_all_cleaned <- align_rows(medianRefPop_rlw_log_all, CI_benchmarks_rlw_log_all_cleaned)
  pctl_medianRefPop_ld_all_cleaned <- align_rows(medianRefPop_ld_all, CI_benchmarks_ld_all_cleaned)
  pctl_medianRefPop_ld_ge0_all_cleaned <- align_rows(medianRefPop_ld_ge0_all, CI_benchmarks_ld_ge0_all_cleaned)
  # [/CHANGE 2026-10]

  # [CHANGE 2026-10 | audit] N7: drop the self-comparison column by name (remove_self_col) from all draws
  rgw_combined <- lapply(rgw_combined, function(pop_es_list) {
    remove_self_col(pop_es_list, self_col)
  })

  rlw_combined <- lapply(rlw_combined, function(pop_es_list) {
    remove_self_col(pop_es_list, self_col)
  })

  # [CHANGE 2026-10 | Rebecca] rgw_log/rlw_log draws
  rgw_log_combined <- lapply(rgw_log_combined, function(pop_es_list) {
    remove_self_col(pop_es_list, self_col)
  })

  rlw_log_combined <- lapply(rlw_log_combined, function(pop_es_list) {
    remove_self_col(pop_es_list, self_col)
  })

  ld_combined <- lapply(ld_combined, function(pop_es_list) {
    remove_self_col(pop_es_list, self_col)
  })
  # [/CHANGE 2026-10]

  # [CHANGE 2026-10 | Rebecca] hypothesis_rate (rgw > hypo_rate_threshold) and rate_rlw (rlw > threshold_rlw) stored on the object
  # Rate at which each alternative hypothesis's ratio-GORIC(A)-weight
  # bootstrap draws exceed a threshold q (default 1, i.e. how often the
  # alternative hypothesis is preferred over the preferred hypothesis within
  # a given bootstrap draw) -- one rate per alternative hypothesis (self
  # column already dropped from rgw_combined above, same convention used
  # throughout). Computed here (at the 'hypo_rate_threshold' this function
  # was called with -- benchmark_means()/benchmark_asymp() expose this as
  # their own 'hypo_rate_threshold' argument, default 1) so it is a proper
  # output field on the returned benchmark object rather than only
  # appearing as a side effect of calling print.benchmark() -- previously
  # x$hypothesis_rate was computed fresh inside print.benchmark() on its
  # own local copy of x and was never saved back onto the object the user
  # holds. print.benchmark() still recomputes this itself (so its own
  # 'hypo_rate_threshold' argument keeps letting a user ask for a different
  # threshold at print time, without rerunning the -- often expensive --
  # bootstrap, as shown in the Guidelines vignette); that recomputation just
  # overwrites this value for the duration of that print call.
  #
  # Unlike print.benchmark() (which hides the "No-effect" category's
  # hypothesis rate from what gets printed/plotted, since there's no
  # alternative hypothesis to compare against under the null and it reads
  # as confusing there), this object-level field keeps the actually
  # computed value for every category, including "No-effect" -- it's just
  # not shown anywhere by default; a user who wants it can read it directly
  # off the object (x$hypothesis_rate[["pop_es = No-effect"]]).
  hypothesis_rate <- lapply(rgw_combined, calculate_hypothesis_rate, q = hypo_rate_threshold)

  # Same mechanism as hypothesis_rate above, but for the log-likelihood-weight
  # ratio (rlw) rather than the GORIC(A)-weight ratio (rgw) -- the rate at
  # which each alternative hypothesis's rlw bootstrap draws exceed a
  # threshold q (its own 'threshold_rlw', independent of
  # hypothesis_rate's 'hypo_rate_threshold', since the two ratios need not
  # be evaluated at the same cutoff). Computed on the raw (unfolded) rlw --
  # NOT on rlw_ge1 -- since folding to always-≥1 discards the direction
  # a threshold rate depends on (a folded ratio exceeding 1 is true for
  # almost every draw, telling you nothing). Also not needed for
  # rgw_log/rlw_log: since log() is strictly monotonic on a positive ratio,
  # rate(log(rgw) > log(q)) is identical to rate(rgw > q) -- no new
  # information, just a log-transformed threshold.
  #
  # Named 'rate_rlw' rather than 'hypothesis_rate_lw' (and its argument
  # 'threshold_rlw' rather than 'hypo_rate_threshold_lw'), unlike the rgw
  # version: rgw is literally the quantity GORIC(A) hypothesis selection is
  # based on, so "rate rgw exceeds q" reads as a hypothesis-support rate.
  # rlw doesn't carry that same interpretation, so calling it a "hypothesis
  # rate" would be misleading -- it's just the exceedance rate of rlw.
  rate_rlw <- lapply(rlw_combined, calculate_hypothesis_rate, q = threshold_rlw)
  # [/CHANGE 2026-10]

  # [CHANGE 2026-10 | Rebecca] overlap of the reference population with the other populations, per output type
  # Overlap (0-1 overlapping coefficient) between the reference population's
  # benchmark distribution and each other population's -- the numeric
  # counterpart of the overlap visible when plotting them together. The
  # reference population is "Observed" when present, else the last
  # pop_es/pop_est population supplied (see determine_overlap_reference());
  # it is the same population across gw/lw/rgw/rlw/rgw_log/rlw_log/ld since
  # they're all built from the same set of populations, so it only needs
  # computing once here and is exposed below as overlap_reference_pop for
  # print.R to build the "Overlap with <name>" column header.
  # Computed separately for gw, lw, rgw, rlw and ld (rather than assuming
  # gw and rgw -- or lw and rlw -- give the same overlap): rgw is only a
  # bijective transform of gw (and so shares its overlap) in the special
  # case of exactly 2 hypotheses; with more hypotheses, or for lw/rlw
  # (lw is not normalized to sum to 1 the way gw is), that need not hold.
  overlap_reference_pop <- reference_pop_name
  overlap_gw  <- compute_overlap_vs_observed(gw_combined)
  overlap_lw  <- compute_overlap_vs_observed(lw_combined)
  overlap_rgw <- compute_overlap_vs_observed_matrix(rgw_combined)
  overlap_rlw <- compute_overlap_vs_observed_matrix(rlw_combined)
  overlap_rgw_log <- compute_overlap_vs_observed_matrix(rgw_log_combined)
  overlap_rlw_log <- compute_overlap_vs_observed_matrix(rlw_log_combined)
  overlap_ld  <- compute_overlap_vs_observed_matrix(ld_combined)
  # [/CHANGE 2026-10]

  OUT <- list(
    benchmarks_gw = CI_benchmarks_gw,
    # [CHANGE 2026-10 | Rebecca] lw benchmarks
    benchmarks_lw = CI_benchmarks_lw,
    benchmarks_rgw = CI_benchmarks_rgw_all_cleaned,
    benchmarks_rlw = CI_benchmarks_rlw_all_cleaned,
    benchmarks_rlw_ge1 = CI_benchmarks_rlw_ge1_all_cleaned,
    # [CHANGE 2026-10 | Rebecca] rgw_log/rlw_log benchmarks
    benchmarks_rgw_log = CI_benchmarks_rgw_log_all_cleaned,
    benchmarks_rlw_log = CI_benchmarks_rlw_log_all_cleaned,
    benchmarks_difLL = CI_benchmarks_ld_all_cleaned,
    benchmarks_absdifLL = CI_benchmarks_ld_ge0_all_cleaned,
    # [CHANGE 2026-10 | Rebecca] new output fields: Pctl. Sample, Pctl. median ref. pop., overlap (with reference), hypothesis_rate/rate_rlw with thresholds
    pctl_Sample_gw  = pctl_Sample_gw,
    pctl_Sample_lw  = pctl_Sample_lw,
    pctl_Sample_rgw = pctl_Sample_rgw_all_cleaned,
    pctl_Sample_rlw = pctl_Sample_rlw_all_cleaned,
    pctl_Sample_rlw_ge1 = pctl_Sample_rlw_ge1_all_cleaned,
    pctl_Sample_rgw_log = pctl_Sample_rgw_log_all_cleaned,
    pctl_Sample_rlw_log = pctl_Sample_rlw_log_all_cleaned,
    pctl_Sample_difLL = pctl_Sample_ld_all_cleaned,
    pctl_Sample_absdifLL = pctl_Sample_ld_ge0_all_cleaned,
    pctl_medianRefPop_gw  = pctl_medianRefPop_gw,
    pctl_medianRefPop_lw  = pctl_medianRefPop_lw,
    pctl_medianRefPop_rgw = pctl_medianRefPop_rgw_all_cleaned,
    pctl_medianRefPop_rlw = pctl_medianRefPop_rlw_all_cleaned,
    pctl_medianRefPop_rlw_ge1 = pctl_medianRefPop_rlw_ge1_all_cleaned,
    pctl_medianRefPop_rgw_log = pctl_medianRefPop_rgw_log_all_cleaned,
    pctl_medianRefPop_rlw_log = pctl_medianRefPop_rlw_log_all_cleaned,
    pctl_medianRefPop_difLL = pctl_medianRefPop_ld_all_cleaned,
    pctl_medianRefPop_absdifLL = pctl_medianRefPop_ld_ge0_all_cleaned,
    overlap_gw = overlap_gw,
    overlap_lw = overlap_lw,
    overlap_rgw = overlap_rgw,
    overlap_rlw = overlap_rlw,
    overlap_rgw_log = overlap_rgw_log,
    overlap_rlw_log = overlap_rlw_log,
    overlap_ld = overlap_ld,
    overlap_reference_pop = overlap_reference_pop,
    hypothesis_rate = hypothesis_rate,
    hypo_rate_threshold = hypo_rate_threshold,
    rate_rlw = rate_rlw,
    threshold_rlw = threshold_rlw,
    # [/CHANGE 2026-10]
    combined_values = list(gw_combined = gw_combined,
                           # [CHANGE 2026-10 | Rebecca] lw draws
                           lw_combined = lw_combined,
                           rgw_combined = rgw_combined,
                           rlw_combined = rlw_combined,
                           # [CHANGE 2026-10 | Rebecca] rgw_log/rlw_log draws
                           rgw_log_combined = rgw_log_combined,
                           rlw_log_combined = rlw_log_combined,
                           ld_combined = ld_combined)
  )

  return(OUT)
}


# Calculate the error probability
calculate_error_probability <- function(object, hypos, pref_hypo, est, 
                                        VCOV, control, ...) {
  # Error probability based on complement of preferred hypothesis in data
  nr_hypos <- dim(object$result)[1]
  # [CHANGE 2026-10 | audit] E2: weights column selected by the criterion type (hard-coded 'gorica.weights' gave numeric(0) for goricc/goricac)
  # The weights column of a goric object is named after its criterion:
  # "goric.weights", "gorica.weights", "goricc.weights" or "goricac.weights"
  # (see goric()). It is therefore selected by the object's type rather than
  # hard-coded (previously 'gorica.weights' was used for every non-goric
  # type, which gave numeric(0) for goricc/goricac objects).
  weights_col <- function(fit) paste0(fit$type, ".weights")
  # [/CHANGE 2026-10]
  if (nr_hypos == 2 && object$comparison == "complement") {
    # TO DO also here re-run with GORICA, as we do for sample value as well?
    #       is ws al opgelost als we goric en gorica resultaten gelijk maken!!!
    #       Dus dan laten staan + re-run met gorica niet nodig dan ook!
    # [CHANGE 2026-10 | audit] E2: weights column by type
    error_prob <- 1 - object$result[[weights_col(object)]][pref_hypo]
  } else {
    if (pref_hypo == nr_hypos && object$comparison == "unconstrained") {
      error_prob <- "The unconstrained (i.e., the failsafe) containing all possible orderings is preferred."
    } else {
      H_pref <- hypos[[pref_hypo]]
      # [CHANGE 2026-10 | audit] error probability: same penalty_factor as the (refitted) object instead of the default; priors not applied to the two-model comparison
      # The preferred hypothesis is compared with its complement using the
      # same criterion (gorica/gorica(c)), sample size and penalty_factor as
      # the (refitted) benchmarked object (otherwise, the error probability
      # would be based on the default penalty_factor). The priorICweights of
      # the object are not used: they refer to the original set of
      # hypotheses and do not apply to this two-model comparison (preferred
      # hypothesis vs its complement), which uses equal prior weights.
      penalty_factor <- object$penalty_factor
      if (is.null(penalty_factor)) {
        penalty_factor <- 2
      }
      # [/CHANGE 2026-10]
      if (is.null(object$model.org)) {
        results_goric_pref <- goric(est, VCOV = VCOV,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    # [CHANGE 2026-10 | audit] E2: same criterion, sample size and penalty_factor as the benchmarked object
                                    type = object$type,
                                    sample_nobs = object$sample_nobs,
                                    penalty_factor = penalty_factor,
                                    control = control,
                                    ...)
      } else {
        fit_data <- object$model.org
        results_goric_pref <- goric(fit_data,
                                    hypotheses = list(H_pref = H_pref),
                                    comparison = "complement",
                                    type = object$type,
                                    # [CHANGE 2026-10 | audit] same penalty_factor as the benchmarked object
                                    penalty_factor = penalty_factor,
                                    control = control, 
                                    ...)
      }
      # [CHANGE 2026-10 | audit] E2: weights column by type
      error_prob <- results_goric_pref$result[[weights_col(results_goric_pref)]][2]
    }
  }
  return(error_prob)
}


# [CHANGE 2026-10 | Rebecca] new function check_iter_adequacy(): stability check of the Sample value's percentile for a fixed iter (messages), with informational median-bias check
# Diagnostic check: is 'iter' large enough for a stable benchmark?
#
# NOTE on what "stable" means here -- and what it deliberately does NOT mean.
# Under the "Observed" population (only present when pop_es/pop_est was left
# NULL, so it is auto-added), the benchmark draws are simulated centered
# exactly on the observed estimate. It is tempting to then expect the
# observed goric(a) weight to sit at the 50th percentile of its own benchmark
# distribution -- but that is NOT generally true, even with unlimited draws.
# gw/lw/rgw/rlw/ld are nonlinear functions of the (multivariate normal)
# bootstrap draws, computed through an order-restricted/inequality-constrained
# log-likelihood, which involves projecting each draw onto a constraint cone.
# That projection is nonlinear, and once more than one constraint is
# simultaneously relevant -- i.e. the observed estimates sit close to more
# than one boundary of the hypothesis at once -- the resulting distribution
# need not be symmetric around the plug-in value, so its median can genuinely
# differ from the observed value. This is the same mechanism behind
# chi-bar-square (mixture) distributions in order-restricted inference, and
# is analogous to why bootstrap bias-correction (the "BC" in BCa) exists at
# all: a nonlinear statistic's bootstrap median need not equal the original
# estimate. So a percentile far from 50 is not necessarily a sign that 'iter'
# is too low -- it can be a genuine feature of the benchmark distribution,
# especially when (some of) the hypotheses are close to being (an) equality
# constraint(s) -- and is not something to strive to eliminate by adjusting
# 'iter'. Because of this, "is the percentile near 50" is NOT what is
# checked (or acted on) here; see goric_percentile_test() for a diagnostic
# that DOES test that directly, kept purely informational for now (returned,
# not printed, and not used to decide anything).
#
# What IS checked here is whether the percentile has stabilized: whether
# growing the number of draws still meaningfully changes where the sample
# value falls in the benchmark distribution. If it does not, more draws are
# unlikely to change anything else about the benchmark either, regardless of
# whether that percentile happens to be near 50.
#
# This is checked using 'gw' rather than 'rgw'/'rlw' on purpose: gw is a
# single, bounded [0, 1] weight per draw, whereas rgw/rlw = gw_pref / gw_k is
# a ratio of two such weights (see calc_ICweights() in
# goric_calculate_IC_weights.R). When gw_k is small, that ratio amplifies
# small (essentially unavoidable, e.g. optimizer-precision-level)
# fluctuations in gw_k multiplicatively, so rgw/rlw are a noisier basis for
# this particular check than gw itself.
check_iter_adequacy <- function(benchmark_results, observed_name, iter,
                                band = c(0.495, 0.505),
                                iter_min = 500, iter_step = 100, iter_max = 2000,
                                stability_tol = 1,
                                control = list(), ...) {
  if (!observed_name %in% names(benchmark_results$pctl_Sample_gw)) {
    # No "Observed" population in this run (user supplied a custom pop_es/
    # pop_est without an "Observed" category) -- nothing to check.
    return(invisible(NULL))
  }

  sample_gw <- benchmark_results$benchmarks_gw[[observed_name]][1, 1]
  sample_lw <- benchmark_results$benchmarks_lw[[observed_name]][1, 1]
  gw_draws  <- benchmark_results$combined_values$gw_combined[[observed_name]]
  lw_draws  <- benchmark_results$combined_values$lw_combined[[observed_name]]
  gw_draws  <- gw_draws[!is.na(gw_draws)]
  lw_draws  <- lw_draws[!is.na(lw_draws)]

  # Informational only, for now: a formal GORICA-based test of whether the
  # sample value sits near the 50th percentile of its own ('Observed')
  # benchmark distribution -- see the note above this function for why that
  # need not hold even asymptotically. Returned on the object (see
  # median_bias_check_gw/median_bias_check_lw in benchmark_means()/
  # benchmark_asymp()) but not printed and not used below to decide anything.
  median_bias_check_gw <- goric_percentile_test(gw_draws, sample_gw, band = band, control = control, ...)
  median_bias_check_lw <- goric_percentile_test(lw_draws, sample_lw, band = band, control = control, ...)

  # Suggests running with 'iter' left at its default (iter = NULL) instead of
  # manually raising a too-low fixed 'iter' -- only sensible if the adaptive
  # procedure could actually try more draws than the user's fixed 'iter', i.e.
  # if that 'iter' is below iter_max; if they already used iter_max or more,
  # the automatic procedure wouldn't go any further either.
  suggest_default <- function() {
    if (iter >= iter_max) return("")
    paste0(
      " Instead of manually raising 'iter', you could also run this with 'iter' left at its ",
      "default (iter = NULL): draws are then added automatically, starting at iter_min = ",
      iter_min, " and increasing by iter_step = ", iter_step, " at a time, up to a maximum of ",
      "iter_max = ", iter_max, "."
    )
  }

  # Stability check: since 'iter' here is a single, fixed, non-growing run
  # (no rounds to compare across draw-by-draw, unlike the auto-'iter' case in
  # run_benchmark_simulation()), stability is instead assessed by comparing
  # the percentile computed from the first 80% of the draws against the
  # percentile computed from the full 'iter' draws -- i.e. did the last 20%
  # of draws still meaningfully move where the sample value falls in the
  # benchmark distribution?
  n <- length(gw_draws) # == length(lw_draws)
  n80 <- max(1, floor(0.8 * n))
  percentile_gw_80  <- 100 * mean(gw_draws[seq_len(n80)] <= sample_gw)
  percentile_lw_80  <- 100 * mean(lw_draws[seq_len(n80)] <= sample_lw)
  percentile_gw_full <- 100 * mean(gw_draws <= sample_gw)
  percentile_lw_full <- 100 * mean(lw_draws <= sample_lw)
  stable_gw <- abs(percentile_gw_full - percentile_gw_80) < stability_tol
  stable_lw <- abs(percentile_lw_full - percentile_lw_80) < stability_tol

  if (!(stable_gw && stable_lw)) {
    # The closing suggestion is tailored to which output_type(s) actually
    # have not (yet) stabilized at this 'iter': if only one of 'gw'/'lw' is
    # the problem, raising 'iter' only matters to someone who cares about
    # that one; someone only interested in the other (already-stable) one
    # doesn't need to.
    describe_type <- function(pct_80, pct_full) {
      paste0(
        "percentile went from ", sprintf("%.1f", pct_80), " (based on the first 80% of ",
        "'iter') to ", sprintf("%.1f", pct_full), " (based on the full 'iter')."
      )
    }
    closing <- if (!stable_gw && !stable_lw) {
      "Consider increasing 'iter' for a more stable benchmark."
    } else if (!stable_gw) {
      paste0(
        "Consider increasing 'iter' for a more stable benchmark when you are (also) interested ",
        "in GORIC(A) weights ('gw') and not (only) in log-likelihood weights ('lw')."
      )
    } else {
      paste0(
        "Consider increasing 'iter' for a more stable benchmark when you are (also) interested ",
        "in log-likelihood weights ('lw') and not (only) in GORIC(A) weights ('gw')."
      )
    }
    message(
      "\nrestriktor Message: For the user-specified 'iter' = ", iter, ", the percentile of the ",
      "value based on your data (called the 'Sample' value in the output) within the benchmark ",
      "distribution under the 'Observed' population has not (yet) stabilized: it changed by ",
      # [CHANGE 2026-10 | audit] N8: message states stability_tol; note when iter_stability_tol <= 0
      stability_tol, " percentage point(s) or more over the last 20% of the draws",
      if (stability_tol <= 0) {
        paste0(" (note that with iter_stability_tol = ", stability_tol, ", the percentile ",
               "can never be considered stable)")
      } else "",
      ".\n",
      # [/CHANGE 2026-10]
      "For output_type = 'gw', ", describe_type(percentile_gw_80, percentile_gw_full), "\n",
      "For output_type = 'lw', ", describe_type(percentile_lw_80, percentile_lw_full), "\n",
      closing, suggest_default()
    )
  } else if (iter < iter_min) {
    # This stability check itself has little power at a very low 'iter'
    # (comparing an 80%/20% split of, say, 10 draws is barely more than
    # comparing single draws) -- so a "looks stable" verdict can slip
    # through undetected. Flag any 'iter' below iter_min (the adaptive
    # procedure's own starting point) outright, regardless of what the
    # stability check itself concluded, rather than silently trusting a
    # check that may not have had the power to catch real instability.
    message(
      "\nrestriktor Message: The user-specified 'iter' = ", iter, " is below the recommended ",
      "starting point of iter_min = ", iter_min, " draws. In that case, there may be too little ",
      "power to detect whether the percentile of the 'Sample value' under the 'Observed' ",
      "population has actually stabilized.\n",
      "You may want to consider increasing 'iter'; either manually or by running the code with ",
      "'iter' left at its default (iter = NULL): draws are then added automatically, starting ",
      "at iter_min = ", iter_min, " and increasing by iter_step = ", iter_step, " at a time, up ",
      "to a maximum of iter_max = ", iter_max, "."
    )
  }
  invisible(list(median_bias_check_gw = median_bias_check_gw,
                median_bias_check_lw = median_bias_check_lw))
}
# [/CHANGE 2026-10]


# [CHANGE 2026-10 | Rebecca] new function goric_percentile_test(): GORICA-based test whether the Sample value is near the median of the draws (informational)
# Is 'sample_value' consistent with being the median of 'draws'? Rather than
# a binomial test against p = 0.5, this expresses "close to the median" as a
# small, fixed interval around 0.5 (band -- a region of practical
# equivalence, deliberately NOT scaled to the percentile estimate's own
# standard error, since doing so would just reconstruct a Wald interval that
# goric(a) would then re-analyze with that same standard error a second
# time) and lets goric(a) itself -- type = "gorica", using the percentile
# estimate's own sampling variance as VCOV -- weigh the evidence for "the
# true percentile lies in that band" against its complement. This keeps the
# whole diagnostic inside goric(a)'s own machinery rather than a separate
# hypothesis-testing framework. Returns the percentile (0-100) 'sample_value'
# falls at within 'draws', and the resulting gorica(a) weight for the
# "in-band" hypothesis (>= 0.5 is taken to mean converged, i.e. goric(a)
# favors "close to the median" over its complement).
#
# NOTE: despite the name, this is NOT a check for whether 'iter' is large
# enough -- see the long note above check_iter_adequacy() for why the
# 'Observed' population's median need not be at the 50th percentile at all
# (a real feature of the benchmark distribution near constraint boundaries,
# not a sign of too few draws). This function is kept around as an
# informational, GORICA-based diagnostic of that median-vs-sample-value gap
# specifically -- both check_iter_adequacy() and run_benchmark_simulation()
# call it and return its result on the benchmark object (median_bias_check_gw
# / median_bias_check_lw) for possible future use, but neither currently acts
# on it or prints it; the adequacy checks themselves are based on whether the
# percentile has stabilized with more draws instead (see both functions).
goric_percentile_test <- function(draws, sample_value, band = c(0.495, 0.505),
                                  control = list(), ...) {
  draws <- draws[!is.na(draws)]
  n <- length(draws)
  k <- sum(draws <= sample_value)
  phat <- k / n

  # Continuity-corrected proportion, used for the variance only, so that
  # phat = 0 or 1 (possible at small iter) doesn't collapse VCOV to exactly
  # zero.
  phat_v <- (k + 0.5) / (n + 1)
  VCOV <- matrix(phat_v * (1 - phat_v) / n)
  rownames(VCOV) <- colnames(VCOV) <- "p"
  est <- c(p = phat)

  H1 <- paste(band[1], "< p <", band[2])
  fit <- goric(est, VCOV = VCOV, hypotheses = list(H1 = H1),
              comparison = "complement", type = "gorica",
              control = control, ...)
  # [CHANGE 2026-10 | audit] E2: weights column by type
  gw_H1 <- fit$result[[paste0(fit$type, ".weights")]][fit$result$model == "H1"]

  list(percentile = 100 * phat, n = n, band = band,
      gw = gw_H1, converged = isTRUE(gw_H1 >= 0.5))
}
# [/CHANGE 2026-10]


# [CHANGE 2026-10 | Rebecca] description of run_benchmark_simulation() (adaptive iter)
# Run the pop_es/pop_est simulation loop used by benchmark_means()/
# benchmark_asymp(), growing the number of draws adaptively when the user
# leaves 'iter' unspecified (iter = NULL): start at iter_min draws and, if
# the "Observed" population's percentile (for output_type = 'gw' and/or
# 'lw') is still changing meaningfully round to round, add iter_step more
# draws -- WITHOUT discarding or redrawing the ones already computed --
# [CHANGE 2026-10 | audit] N8: stop after two consecutive stable rounds
# repeating until the percentile has stabilized (i.e., stayed within
# stability_tol for two consecutive rounds) or iter_max is reached.
# Growth is deliberately based on STABILITY of the percentile, not on
# whether it is close to 50 -- see the long note above check_iter_adequacy()
# for why the 'Observed' population's median need not be at the 50th
# percentile at all, so that is not something more draws can be expected to
# fix. If the user supplies a fixed numeric 'iter' instead, this runs
# exactly one round of that many draws (the original, non-adaptive
# behaviour); benchmark_means()/benchmark_asymp() then call
# check_iter_adequacy() themselves afterwards for that fixed-iter case (this
# function does not, to avoid messaging twice).
#
# Returns list(parallel_function_results = <as before, one element per
# pop_es/pop_est category>, iter = <final number of draws used>,
# median_bias_check_gw/median_bias_check_lw = <the last round's
# goric_percentile_test() result, informational only -- see the note above
# that function>).
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] N8: new validate_iter_args(): validation of iter, iter_min/step/max, iter_stability_tol and iter_adequacy_band
# Checks the iter-related arguments of benchmark_means()/benchmark_asymp()
# up front, so that invalid values give a clear error instead of an obscure
# failure (or an endless/empty loop) inside run_benchmark_simulation().
validate_iter_args <- function(iter, iter_min, iter_step, iter_max,
                               iter_stability_tol, iter_adequacy_band) {
  is_count <- function(v, min_value) {
    is.numeric(v) && length(v) == 1 && !is.na(v) && is.finite(v) &&
      v >= min_value && v == round(v)
  }
  if (!is.null(iter) && !is_count(iter, 1)) {
    stop("\nrestriktor ERROR: The argument 'iter' should be NULL (the default; then the ",
         "number of draws is set automatically) or a single whole number >= 1.", call. = FALSE)
  }
  if (!is_count(iter_min, 1)) {
    stop("\nrestriktor ERROR: The argument 'iter_min' should be a single whole number >= 1.",
         call. = FALSE)
  }
  if (!is_count(iter_step, 1)) {
    stop("\nrestriktor ERROR: The argument 'iter_step' should be a single whole number >= 1.",
         call. = FALSE)
  }
  if (!is_count(iter_max, 1)) {
    stop("\nrestriktor ERROR: The argument 'iter_max' should be a single whole number >= 1.",
         call. = FALSE)
  }
  if (!(is.numeric(iter_stability_tol) && length(iter_stability_tol) == 1 &&
        !is.na(iter_stability_tol) && iter_stability_tol >= 0)) {
    stop("\nrestriktor ERROR: The argument 'iter_stability_tol' should be a single number >= 0 ",
         "(in percentage points).", call. = FALSE)
  }
  if (!(is.numeric(iter_adequacy_band) && length(iter_adequacy_band) == 2 &&
        !anyNA(iter_adequacy_band) && all(iter_adequacy_band >= 0) &&
        all(iter_adequacy_band <= 1) && iter_adequacy_band[1] < iter_adequacy_band[2])) {
    stop("\nrestriktor ERROR: The argument 'iter_adequacy_band' should be a vector of two ",
         "increasing proportions between 0 and 1 (e.g., c(0.495, 0.505)).", call. = FALSE)
  }
  invisible(TRUE)
}
# [/CHANGE 2026-10]


# [CHANGE 2026-10 | Rebecca] new function run_benchmark_simulation(): simulation loop for benchmark_means()/benchmark_asymp() with adaptive iter (grow by iter_step until the percentile is stable), parallel draws and progress messages
run_benchmark_simulation <- function(nr_es, rnames, name_prefix, center_matrix,
                                     colnames_vec, VCOV, hypos, pref_hypo,
                                     comparison, control, mix_weights,
                                     penalty_factor, Heq, object, iter,
                                     es_labels = rnames,
                                     iter_min = 500, iter_step = 100,
                                     iter_max = 2000, band = c(0.495, 0.505),
                                     # [CHANGE 2026-10 | audit] E2/B24: arguments stability_tol, type, sample_nobs (same criterion/sample size as the object) and observed_label
                                     stability_tol = 1,
                                     # criterion (gorica/goricac) and sample
                                     # size used for every draw: the same as
                                     # those of the (refitted) benchmarked
                                     # 'object', so that the 'Sample' value
                                     # and the benchmark distribution are based on
                                     # the same criterion.
                                     type = object$type,
                                     sample_nobs = object$sample_nobs,
                                     # how the 'Observed' population is
                                     # referred to in the messages below
                                     observed_label = "the 'Observed' population",
                                     ...) {
                                     # [/CHANGE 2026-10]
  # [CHANGE 2026-10 | audit] E2: goricac requires the sample size
  if (type == "goricac" && is.null(sample_nobs)) {
    stop("\nrestriktor ERROR: The benchmark is based on the GORICAC (goricac), which ",
         "requires the sample size. Please specify it via the argument 'sample_size' ",
         "(benchmark_asymp) or 'group_size' (benchmark_means).", call. = FALSE)
  }
  # [/CHANGE 2026-10]

  auto_iter <- is.null(iter)
  sample_gw <- object$result[pref_hypo, 7]
  sample_lw <- object$result$loglik.weights[pref_hypo]
  # [CHANGE 2026-10 | audit] B34: comment: populations indexed by position (names made unique)
  # The populations are indexed by POSITION throughout (the names are made
  # unique by benchmark_means()/benchmark_asymp(), but position is what the
  # results are accumulated by).
  obs_pos   <- which(rnames == "Observed")

  parallel_function_results <- vector("list", nr_es)
  names(parallel_function_results) <- paste0(name_prefix, rnames)
  est_accum <- vector("list", nr_es)

  n_done <- 0L
  target <- if (auto_iter) min(iter_min, iter_max) else iter

  # Previous round's percentile per output_type, used to detect when it has
  # stopped moving (see 'stabilized' below) -- NA until a second round
  # exists to compare against.
  prev_percentile_gw <- NA_real_
  prev_percentile_lw <- NA_real_
  # [CHANGE 2026-10 | audit] N8: count consecutive stable rounds and the batch sizes of added draws
  # Number of consecutive rounds in which the percentile stayed within
  # stability_tol (for both 'gw' and 'lw'), and the sizes of the batches of
  # draws added in those rounds (the last batch can be smaller than
  # iter_step when it is capped by iter_max).
  n_stable_rounds <- 0L
  batch_sizes <- integer(0)
  # [/CHANGE 2026-10]

  # NULL unless/until computed inside the repeat loop below -- stays NULL for
  # a fixed, user-specified 'iter' (auto_iter = FALSE), since that path
  # breaks out of the loop before reaching the percentile-stability check
  # (benchmark_means()/benchmark_asymp() do that check themselves afterwards,
  # via check_iter_adequacy(), for the fixed-iter case).
  chk_gw <- NULL
  chk_lw <- NULL

  # [CHANGE 2026-10 | audit] A9: successful draws per population; B32: progress handler local to this call (global progressr handlers untouched)
  # Successful draws per population (see the worker: a draw only fails on an
  # error in goric(); warnings are muffled and the draw is kept).
  n_ok <- function(pos) {
    sum(!vapply(parallel_function_results[[pos]], draw_failed, logical(1)))
  }
  # Progress bar for this call only: the handler is passed to
  # with_progress() rather than set globally via progressr::handlers(), so
  # the user's own (global) progressr handlers are left untouched.
  progress_handler <- progressr::handler_txtprogressbar(char = ">")
  # [/CHANGE 2026-10]

  repeat {
    n_new <- target - n_done

    progressr::with_progress({
      p <- progressr::progressor(along = seq_len(n_new * nr_es))

      for (teller_es in seq_len(nr_es)) {
        cat("Calculating benchmark for", es_labels[teller_es],
            "-- draws", n_done + 1, "to", target, "\n")

        new_draws <- mvtnorm::rmvnorm(n = n_new, center_matrix[teller_es, ], sigma = VCOV)
        colnames(new_draws) <- colnames_vec
        est_accum[[teller_es]] <- rbind(est_accum[[teller_es]], new_draws)
        est_full <- est_accum[[teller_es]]

        # Wrapper function for future_lapply
        wrapper_function_asymp <- function(i) {
          p() # Update progress
          parallel_function_asymp(i,
                                  est = est_full, VCOV = VCOV,
                                  hypos = hypos, pref_hypo = pref_hypo,
                                  # [CHANGE 2026-10 | audit] E2: same criterion (gorica/goricac) and sample size for every draw
                                  comparison = comparison, type = type,
                                  sample_nobs = sample_nobs,
                                  control = control, mix_weights = mix_weights,
                                  penalty_factor = penalty_factor,
                                  # [CHANGE 2026-10 | audit] N3: prior weights of the goric object
                                  # same prior weights as the (refitted) goric object
                                  priorICweights = object$priorICweights,
                                  Heq = Heq, ...)
        }

        new_results <- future_lapply(
          seq(n_done + 1, target),
          wrapper_function_asymp,
          future.seed = TRUE # Ensures safe and reproducible random number generation
        )

        # [CHANGE 2026-10 | audit] A9: append the new draws (incl. failed ones, for the accounting below)
        parallel_function_results[[teller_es]] <- c(parallel_function_results[[teller_es]],
                                                    new_results)
      }
    # [CHANGE 2026-10 | audit] B32: local progress handler
    }, handlers = progress_handler)

    # [CHANGE 2026-10 | audit] N8: record the size of this batch of draws
    batch_sizes <- c(batch_sizes, as.integer(target - n_done))
    n_done <- target

    if (!auto_iter) break # fixed iter: exactly one round, done

    if (length(obs_pos) == 1) {
      # [CHANGE 2026-10 | audit] A9: failed draws (draw_failed) give NA
      gw_draws <- vapply(parallel_function_results[[obs_pos]], function(x) {
        if (draw_failed(x)) NA_real_ else as.numeric(x$gw)
      }, numeric(1))
      # [CHANGE 2026-10 | audit] A9: failed draws (draw_failed) give NA
      lw_draws <- vapply(parallel_function_results[[obs_pos]], function(x) {
        if (draw_failed(x)) NA_real_ else as.numeric(x$lw)
      }, numeric(1))
      gw_draws <- gw_draws[!is.na(gw_draws)]
      lw_draws <- lw_draws[!is.na(lw_draws)]
      # [CHANGE 2026-10 | audit] A9: no successful draw of the 'Observed' population yet: nothing to check (failures reported after the loop)
      if (length(gw_draws) == 0) {
        # every draw of the 'Observed' population failed so far: nothing to
        # check; the failures are reported (and, if they persist, turned
        # into an error) after the loop.
        stabilized <- TRUE
        chk_gw <- chk_lw <- NULL
      } else {
      # [/CHANGE 2026-10]
      # Informational only, for now -- see the note above goric_percentile_test().
      chk_gw <- goric_percentile_test(gw_draws, sample_gw, band = band, control = control)
      chk_lw <- goric_percentile_test(lw_draws, sample_lw, band = band, control = control)

      # Has the percentile stopped moving with more draws? This -- not
      # closeness to the 50th percentile -- is what growth is based on; see
      # the long note above check_iter_adequacy() for why the percentile
      # need not be near 50 at all, even once it has fully stabilized.
      # FALSE (not NA) on the very first round, since there is no earlier
      # [CHANGE 2026-10 | audit] N8: comment: the percentile has to be stable for two consecutive rounds
      # percentile yet to compare against. Since a percentile can drift by a
      # percentage point or two by chance alone, a single round-to-round
      # comparison could call "stabilized" too eagerly; therefore, the
      # percentile has to stay within stability_tol for TWO consecutive
      # rounds (i.e., over the last two batches of added draws) before the
      # growing stops.
      # [/CHANGE 2026-10]
      stable_gw <- isTRUE(!is.na(prev_percentile_gw) &&
                            abs(chk_gw$percentile - prev_percentile_gw) < stability_tol)
      stable_lw <- isTRUE(!is.na(prev_percentile_lw) &&
                            abs(chk_lw$percentile - prev_percentile_lw) < stability_tol)
      # [CHANGE 2026-10 | audit] N8: stop growing only after two consecutive stable rounds
      n_stable_rounds <- if (stable_gw && stable_lw) n_stable_rounds + 1L else 0L
      stabilized <- n_stable_rounds >= 2L
      prev_percentile_gw <- chk_gw$percentile
      prev_percentile_lw <- chk_lw$percentile
      # [CHANGE 2026-10 | audit] A9: end of the else-branch (successful draws available)
      }
    } else {
      # No "Observed" category (custom pop_es/pop_est) -- nothing to check
      # against, so stop growing after the first (iter_min-sized) batch.
      stabilized <- TRUE
      chk_gw <- chk_lw <- NULL
    }

    if (stabilized || n_done >= iter_max) {
      if (auto_iter && !is.null(chk_gw)) {
        # [CHANGE 2026-10 | audit] N8/A9: actual batch sizes and the number of successful draws for the messages
        # number of draws over which the last (one or two) comparison(s)
        # were made; the last batch can be smaller than iter_step if capped
        # by iter_max.
        n_last <- if (length(batch_sizes) > 1) batch_sizes[length(batch_sizes)] else NA
        n_last2 <- if (length(batch_sizes) > 2) sum(utils::tail(batch_sizes, 2)) else NA
        # the number of draws reported is the number of SUCCESSFUL draws
        # (failed draws, if any, are reported separately below)
        n_used <- n_ok(obs_pos)
        percentiles_txt <- paste0(
          "For output_type = 'gw', the percentile is ", sprintf("%.1f", chk_gw$percentile),
          "; for output_type = 'lw', the percentile is ", sprintf("%.1f", chk_lw$percentile),
          ".\n"
        )
        # [/CHANGE 2026-10]
        if (stabilized) {
          message(
            "\nrestriktor Message: 'iter' was not specified, so it was set automatically.\n",
            # [CHANGE 2026-10 | audit] A9: number of successful draws
            "Using iter = ", n_used, " draws, the percentile of the value based on your data ",
            "(called the 'Sample' value in the output) within the benchmark distribution under ",
            # [CHANGE 2026-10 | audit] N8/B24: message: stable over the last two rounds (observed_label)
            observed_label, " has stabilized: in each of the last two rounds of added ",
            "draws (", n_last2, " draws in total), it changed by less than ", stability_tol,
            " percentage point(s).\n",
            percentiles_txt,
            # [/CHANGE 2026-10]
            "So, no further draws were added. Note that this percentile is not necessarily ",
            "expected to be near 50 -- that need not indicate a problem, particularly when ",
            "(some of) the hypotheses are close to being (an) equality constraint(s)."
          )
        } else {
          # [CHANGE 2026-10 | audit] N8: reason why the percentile did not stabilize (tol <= 0, iter_min >= iter_max, or not stable)
          reason <- if (stability_tol <= 0) {
            paste0("since iter_stability_tol = ", stability_tol, ", the percentile can never be ",
                   "considered stable (it would have to change by less than ", stability_tol,
                   " percentage points).\n")
          } else if (is.na(n_last)) {
            paste0("since there were no further draws to compare the percentile against ",
                   "(iter_min >= iter_max), it could not be checked whether the percentile of ",
                   "the value based on your data (called the 'Sample' value in the output) within ",
                   "the benchmark distribution under ", observed_label, " has stabilized.\n")
          } else {
            paste0("since the percentile of the value based on your data (called the 'Sample' ",
                   "value in the output) within the benchmark distribution under ", observed_label,
                   " had not (yet) stabilized: it did not stay within ", stability_tol,
                   " percentage point(s) for two consecutive rounds of added draws (over the last ",
                   n_last, " draws, it changed by ",
                   if (n_stable_rounds > 0) "less than " else "", stability_tol,
                   " percentage point(s)", if (n_stable_rounds > 0) "" else " or more",
                   ").\n")
          }
          # [/CHANGE 2026-10]
          message(
            # [CHANGE 2026-10 | audit] N8/A9: message when iter_max is reached (successful draws, reason, percentiles)
            "\nrestriktor Message: 'iter' was not specified, so it was ",
            if (is.na(n_last)) "set to" else "increased automatically up to",
            " its maximum of iter = ", n_used, " draws (iter_max = ", iter_max, "), ",
            reason,
            percentiles_txt,
            "Consider re-running with a manually specified, larger 'iter' (or a larger ",
            "'iter_max') for a more stable benchmark."
            # [/CHANGE 2026-10]
          )
        }
      }
      break
    }

    target <- min(n_done + iter_step, iter_max)
  }

  # [CHANGE 2026-10 | audit] A9: accounting of failed (dropped) and warned (kept) draws per population with first messages; error if all draws of a population failed; iter = successful draws
  # Accounting of failed draws (goric() error: dropped) and draws with a
  # (muffled) warning (kept), per population -- see parallel_function_asymp().
  pop_names <- names(parallel_function_results)
  n_failed <- integer(nr_es)
  n_warned <- integer(nr_es)
  first_error <- rep(NA_character_, nr_es)
  first_warning <- rep(NA_character_, nr_es)
  for (teller_es in seq_len(nr_es)) {
    res <- parallel_function_results[[teller_es]]
    failed <- vapply(res, draw_failed, logical(1))
    warned <- vapply(res, function(x) !is.null(attr(x, "warnings")), logical(1))
    n_failed[teller_es] <- sum(failed)
    n_warned[teller_es] <- sum(warned)
    if (any(failed)) {
      first_error[teller_es] <- attr(res[[which(failed)[1]]], "error")
    }
    if (any(warned)) {
      first_warning[teller_es] <- attr(res[[which(warned)[1]]], "warnings")[1]
    }
  }
  names(n_failed) <- names(n_warned) <- names(first_error) <- names(first_warning) <- pop_names
  n_success <- as.integer(n_done) - n_failed
  names(n_success) <- pop_names
  if (any(n_success == 0)) {
    bad <- which(n_success == 0)[1]
    stop("\nrestriktor ERROR: All ", n_done, " benchmark draws for the population '",
         pop_names[bad], "' failed (goric() gave an error for every draw), so no ",
         "benchmark can be computed. The (first) error message was:\n  ",
         first_error[bad], call. = FALSE)
  }
  if (any(n_warned > 0)) {
    message("\nrestriktor Message: In ", sum(n_warned), " of the ", n_done * nr_es,
            " benchmark draws (",
            paste0(pop_names, ": ", n_warned, " of ", n_done, collapse = "; "),
            "), goric() gave a warning. These warnings were muffled and the draws were ",
            "kept (a warning does not invalidate the GORIC(A) result of a draw). The first ",
            "warning message was:\n  ", trimws(first_warning[!is.na(first_warning)][1]))
  }
  if (any(n_failed > 0)) {
    message("\nrestriktor Message: ", sum(n_failed), " of the ", n_done * nr_es,
            " benchmark draws (",
            paste0(pop_names, ": ", n_failed, " of ", n_done, collapse = "; "),
            ") failed (goric() gave an error) and were discarded; the benchmark is based on ",
            "the remaining draws (", paste0(pop_names, ": ", n_success, collapse = "; "),
            "). The first error message was:\n  ", trimws(first_error[!is.na(first_error)][1]))
  }
  # the number of successful draws: a single number when it is the same for
  # every population (the usual case: no failed draws), else one per population
  iter_out <- if (length(unique(n_success)) == 1) unname(n_success[1]) else n_success

  list(parallel_function_results = parallel_function_results, iter = iter_out,
      iter_requested = as.integer(n_done),
      n_failed_draws = n_failed, n_warned_draws = n_warned,
      draw_errors = first_error, draw_warnings = first_warning,
  # [/CHANGE 2026-10]
      median_bias_check_gw = chk_gw, median_bias_check_lw = chk_lw)
}
# [/CHANGE 2026-10]


calculate_hypothesis_rate <- function(x, q = 1) {
  return(colMeans(x > q))
}


# called by the benchmark.print() function
print_section <- function(header, content_printer, nchar, text_color, reset) {
  cat("\n")
  cat(strrep("=", nchar), "\n")
  cat(paste0(text_color, header, reset), "\n")
  cat(strrep("-", nchar), "\n")
  content_printer()
  #cat(strrep("-", nchar), "\n")
  cat("\n")
}


format_value <- function(value) {
  if (is.na(value)) {
    return("")
  }
  if (abs(value) >= 1000 || (abs(value) <= 0.001 && value != 0)) {
    return(sprintf("%.3e", value))
  } else {
    return(sprintf("%.3f", value))
  }
}

# [CHANGE 2026-10 | Rebecca] new function format_overlap_value(): self-overlap printed as "1"
# Used only for the "Overlap with ..." column: an overlap of exactly 1 there
# is always the reference population's own (self-)entry, fixed at 1 by
# construction rather than actually computed (see compute_overlap_vs_observed()/
# compute_overlap_vs_observed_matrix() in this file) -- printing "1.000" would
# suggest a computed value carrying that much precision, which is misleading.
# Every other value in this column (and every value in every other column)
# still goes through the normal format_value().
format_overlap_value <- function(value) {
  # [CHANGE 2026-10 | audit] E5/O10: an NA overlap is printed as "NA" (reason appended by print_rounded_es_value())
  if (is.na(value)) {
    # An overlap that could not be computed (see compute_overlap()) is shown
    # as "NA" (print_rounded_es_value() appends the reason, if known) rather
    # than as a blank, so it is not mistaken for 'not applicable'.
    return("NA")
  }
  if (value == 1) {
  # [/CHANGE 2026-10]
    return("1")
  }
  format_value(value)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | Rebecca] new function format_median_ref_pop_value(): reference population's own entry printed as "50"
# Used only for the "Pctl. median ref. pop." column: the reference
# population's own row there is hardcoded to exactly 50 by construction
# (see get_results_benchmark()'s pctl_medianRefPop_* computation, and its
# is_ref_pop/is_reference handling) rather than actually computed -- an
# ecdf() evaluated at a finite/even-n sample's own median need not land on
# exactly 50, so printing "50.000" there would suggest a computed value
# carrying that much precision, which is misleading. Every other value in
# this column (and every value in every other column) still goes through
# format_value().
format_median_ref_pop_value <- function(value) {
  if (!is.na(value) && value == 50) {
    return("50")
  }
  format_value(value)
}
# [/CHANGE 2026-10]


# [CHANGE 2026-10 | Rebecca] new function overlap_column(): builds the "Overlap with <reference>" column for the printed tables
# Builds the "Overlap with Observed" column added to each output_type's
# per-population benchmark table in print.benchmark(). 'overlap_source' is
# the relevant overlap_* field on the benchmark object: for gw/lw it is a
# plain vector named by population (a single overlap value per population);
# for rgw/rlw/ld it is a list named by population, each entry itself a
# vector named by alternative hypothesis (there can be more than one per
# population). 'n_rows' is the number of rows in the benchmark table this
# column is being attached to (nrow() of e.g. benchmarks_ratio_goric_weights
# [CHANGE 2026-10 | audit] N7: values aligned with the table rows by name (row_names, pref_hypo_name)
# [[pop_es]]). For rgw/rlw/ld, the values are aligned with the table's rows
# BY NAME when 'row_names' and 'pref_hypo_name' are supplied (see
# align_by_hypothesis()); otherwise by position. For gw/lw the single
# overlap value is repeated across the (single) row. Returns a column of NA
# when there's nothing to report -- no "Observed" category in this run, or
# (defensively) a length mismatch.
overlap_column <- function(overlap_source, pop_es, n_rows, row_names = NULL,
                           pref_hypo_name = NULL) {
# [/CHANGE 2026-10]
  na_col <- rep(NA_real_, n_rows)
  if (is.null(overlap_source) || !(pop_es %in% names(overlap_source))) {
    # Note: for the gw/lw case, overlap_source is a plain named vector, and
    # `[[` on that with a non-matching name errors ("subscript out of
    # bounds") rather than returning NULL the way it would for a list (the
    # rgw/rlw/ld case) -- so membership is checked explicitly first.
    return(na_col)
  }
  vals <- overlap_source[[pop_es]]
  if (is.null(vals) || length(vals) == 0) {
    return(na_col)
  }
  # [CHANGE 2026-10 | audit] E5/O10 + N7: carry the notes of NA overlaps along; align values and notes by hypothesis name
  # Why an overlap is NA (see compute_overlap()): for gw/lw a single note per
  # population (attribute "notes" on the vector overlap_source, which `[[`
  # above drops); for rgw/rlw/ld one per hypothesis column (attribute on the
  # population's own vector 'vals'). Carried along as a same-length character
  # vector (NA where there is nothing to note), so that
  # print_rounded_es_value() can print "NA (<note>)" in the overlap column.
  notes <- if (is.null(names(vals))) attr(overlap_source, "notes")[pop_es] else attr(vals, "notes")
  notes_full <- rep(NA_character_, length(vals))
  names(notes_full) <- names(vals)
  if (!is.null(notes) && !all(is.na(notes))) {
    if (is.null(names(vals))) {
      notes_full[] <- unname(notes)
    } else {
      notes_full[names(notes)] <- notes
    }
  }
  # rgw/rlw/ld case with the table's rownames available: align by name (see
  # align_by_hypothesis()).
  if (!is.null(names(vals)) && !is.null(row_names) && !is.null(pref_hypo_name)) {
    out <- align_by_hypothesis(vals, row_names, pref_hypo_name)
    attr(out, "notes") <- align_by_hypothesis(notes_full, row_names, pref_hypo_name)
    return(out)
  }
  # [/CHANGE 2026-10]
  vals <- unname(vals)
  # [CHANGE 2026-10 | audit] E5/O10: notes without names
  notes_full <- unname(notes_full)
  if (length(vals) == n_rows) {
    # [CHANGE 2026-10 | audit] E5/O10: return the notes as attribute
    return(structure(vals, notes = notes_full))
  }
  # gw/lw case: a single overlap value for the whole population, repeated
  # across every row (normally just one: the preferred hypothesis).
  if (length(vals) == 1) {
    # [CHANGE 2026-10 | audit] E5/O10: return the notes as attribute
    return(structure(rep(vals, n_rows), notes = rep(notes_full, n_rows)))
  }
  # Length mismatch that isn't the gw/lw broadcast case -- shouldn't
  # normally happen, but pad/truncate defensively rather than risk silently
  # mis-aligning a row with the wrong hypothesis's overlap value.
  out <- na_col
  # [CHANGE 2026-10 | audit] E5/O10: notes for the padded/truncated case
  out_notes <- rep(NA_character_, n_rows)
  keep <- seq_len(min(n_rows, length(vals)))
  out[keep] <- vals[keep]
  out_notes[keep] <- notes_full[keep]
  attr(out, "notes") <- out_notes
  # [/CHANGE 2026-10]
  out
}
# [/CHANGE 2026-10]


# [CHANGE 2026-10 | Rebecca] new function recompute_percentile_table(): recompute the Sample + percentile columns at print time for user-supplied percentiles
# Rebuilds the "Sample" + percentile% columns of an already-built benchmark
# table (e.g. x$benchmarks_goric_weights[[pop_es]]) at PRINT time, for a
# caller-supplied set of percentiles, instead of using the percentiles that
# happened to be requested via 'quant' back when the (often expensive,
# bootstrap-based) benchmark object was computed. Mirrors exactly what
# plot.benchmark()'s own 'percentiles' argument already does (see
# goric_benchmark_plot.R): both recompute quantile() fresh from the raw
# combined draws stored on the object (x$combined_values), so you can look at
# different percentiles without rerunning benchmark_means()/benchmark_asymp().
# 'existing_mat' is the current table (used only for its "Sample" column and
# its row names/count -- both stay unchanged); 'combined_data' is the
# matching raw-draws entry from x$combined_values for this population: a
# plain numeric vector for gw/lw, or a matrix (one column per row of
# 'existing_mat', in the same order -- gw/lw and rgw/rlw/rgw_log/rlw_log/ld
# align this way throughout the codebase, e.g. overlap_column() above relies
# on the same positional correspondence) for rgw/rlw/rgw_log/rlw_log/ld.
# [CHANGE 2026-10 | audit] N7: new align_by_hypothesis(): align per-hypothesis values (hypothesis_rate, rate_rlw, overlap) with the table rows by name
# Aligns a vector of per-alternative-hypothesis values (e.g. hypothesis_rate,
# rate_rlw or overlap; named by the draws' column names, like "vs. H2") with
# the rows of a benchmark table (rownames "<preferred hypothesis> vs. H2", see
# get_results_benchmark()) BY NAME. A single unnamed value (e.g. the NA used
# for the 'No-effect' population in print.benchmark()) is repeated for every
# row. Rows without a matching value get NA.
align_by_hypothesis <- function(vals, row_names, pref_hypo_name) {
  n_rows <- length(row_names)
  if (is.null(vals) || length(vals) == 0) {
    return(rep(NA_real_, n_rows))
  }
  if (is.null(names(vals))) {
    if (length(vals) == 1) return(rep(unname(vals), n_rows))
    if (length(vals) == n_rows) return(unname(vals))
    return(rep(NA_real_, n_rows))
  }
  unname(vals[match(row_names, paste(pref_hypo_name, names(vals)))])
}

# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] N7: argument pref_hypo_name (alignment by name)
recompute_percentile_table <- function(existing_mat, combined_data, percentiles,
                                       pref_hypo_name = NULL) {
  sample_col <- existing_mat[, 1, drop = FALSE]
  pct_names <- paste0(percentiles * 100, "%")
  if (is.null(dim(combined_data))) {
    # gw/lw case: a single row, one shared set of draws.
    q <- unname(quantile(combined_data, probs = percentiles, na.rm = TRUE))
    new_mat <- matrix(c(sample_col, q), nrow = nrow(existing_mat))
  } else {
    # rgw/rlw/rgw_log/rlw_log/ld case: one column of draws per row (hypothesis).
    # [CHANGE 2026-10 | audit] N7: draw columns matched to the table rows by name (else by position); robust matrix shape for a single percentile or no rows
    # Columns of draws are matched to the table's rows by name (see
    # align_by_hypothesis()) when possible, otherwise by position.
    col_idx <- seq_len(nrow(existing_mat))
    if (!is.null(pref_hypo_name) && !is.null(colnames(combined_data))) {
      col_idx <- match(rownames(existing_mat),
                       paste(pref_hypo_name, colnames(combined_data)))
    }
    q <- vapply(col_idx, function(j) {
      if (is.na(j) || j > ncol(combined_data)) return(rep(NA_real_, length(percentiles)))
      quantile(combined_data[, j], probs = percentiles, na.rm = TRUE, names = FALSE)
    }, numeric(length(percentiles)))
    # one row per table row, also for a single percentile or no rows at all
    q <- matrix(q, nrow = length(col_idx), ncol = length(percentiles), byrow = TRUE)
    # [/CHANGE 2026-10]
    new_mat <- cbind(sample_col, q)
  }
  colnames(new_mat) <- c("Sample", pct_names)
  rownames(new_mat) <- rownames(existing_mat)
  new_mat
}
# [/CHANGE 2026-10]


# [CHANGE 2026-10 | Rebecca] new function format_threshold_display(): threshold shown in the column header
# Displays a hypo_rate_threshold value for the "hypothesis_rate" column
# header below -- e.g. 1 -> "1", 1.5 -> "1.5" -- without the trailing
# ".000" sprintf("%.3f")-style formatting used elsewhere in this file would
# add for a plain round number.
format_threshold_display <- function(q) {
  if (is.null(q) || is.na(q)) {
    return("threshold")
  }
  as.character(q)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | Rebecca] new function header_lines_for_columns(): two-line (grouped) column headers for the printed tables
# Column-header display text for print_grouped_header_table() below, keyed
# off the (internal, formatter-matching) column names already set on the
# benchmark tables by print.benchmark() -- "Sample", a percentile like "5%",
# "Pctl. Sample", "Pctl. median ref. pop.", "hypothesis_rate", or an
# "Overlap with <population>" header. Returns, per column: 'line1' (blank
# for a column that instead shares a merged group label -- see 'group'
# below), 'line2', and 'group' (NA for a column with its own two-line
# header; otherwise the shared label text for the columns sharing a merged
# top cell, e.g. "Percentiles distributions" for the quantile columns).
# Layout per Rebecca's spec:
#   Sample                  -> "Sample" / "value" (own two-line header)
#   5%, 35%, ... (any N%)   -> grouped "Percentiles distributions", showing
#                              the percentile itself on line 2
#   Pctl. Sample            -> grouped "Percentiles (pctl.) comparisons",
#                              "Sample value" on line 2
#   Pctl. median ref. pop.  -> same group, "Median ref. pop." on line 2
#   hypothesis_rate         -> "Hypothesis rate" / "(rgw > <threshold>)"
#                              (own two-line header; 'hypo_rate_threshold'
#                              is the actual q this particular print call
#                              used -- see print.benchmark() -- shown here
#                              so it's visible at a glance why this printed
#                              rate can differ from x$hypothesis_rate
#                              (computed at x$hypo_rate_threshold, which may
#                              be a different value) or from another
#                              print() call using a different threshold).
#                              Called "Hypothesis rate" because rgw (the
#                              ratio-of-GORIC(A)-weights) is the quantity
#                              GORIC(A) hypothesis selection is actually
#                              based on, so "exceeds q" reads as "how often
#                              this alternative would be preferred".
#   rate_rlw                -> "Rate" / "(rlw > <threshold>)" -- same
#                              mechanism as hypothesis_rate (rate at which
#                              bootstrap draws exceed a threshold, own
#                              independent 'threshold_rlw'), applied to rlw
#                              (ratio-of-log-likelihood-weights) instead.
#                              Deliberately NOT labelled "Hypothesis rate"
#                              here -- rlw isn't the quantity GORIC(A)
#                              hypothesis selection is based on, so calling
#                              this a hypothesis rate would be misleading;
#                              "(rlw > <threshold>)" on its own already says
#                              exactly what's being counted.
#   Overlap with <pop>      -> "Overlap" / "with <pop>" (own two-line header,
#                              second line already dynamic via overlap_header)
# Anything unrecognized (shouldn't normally occur) falls back to a blank
# line 1 and the column's own name on line 2, so nothing is ever dropped.
header_lines_for_columns <- function(colnames_vec, hypo_rate_threshold = NULL,
                                     threshold_rlw = NULL) {
  n <- length(colnames_vec)
  line1 <- character(n)
  line2 <- character(n)
  group <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    cn <- colnames_vec[i]
    if (identical(cn, "Sample")) {
      line1[i] <- "Sample"
      line2[i] <- "value"
    } else if (grepl("^[0-9.]+%$", cn)) {
      group[i] <- "Percentiles distributions"
      line2[i] <- cn
    } else if (identical(cn, "Pctl. Sample")) {
      group[i] <- "Percentiles (pctl.) comparisons"
      line2[i] <- "Sample value"
    } else if (identical(cn, "Pctl. median ref. pop.")) {
      group[i] <- "Percentiles (pctl.) comparisons"
      line2[i] <- "Median ref. pop."
    } else if (identical(cn, "hypothesis_rate")) {
      line1[i] <- "Hypothesis rate"
      line2[i] <- paste0("(rgw > ", format_threshold_display(hypo_rate_threshold), ")")
    } else if (identical(cn, "rate_rlw")) {
      line1[i] <- "Rate"
      line2[i] <- paste0("(rlw > ", format_threshold_display(threshold_rlw), ")")
    } else if (grepl("^Overlap with ", cn)) {
      line1[i] <- "Overlap"
      line2[i] <- sub("^Overlap with ", "with ", cn)
    } else {
      line2[i] <- cn
    }
  }
  list(line1 = line1, line2 = line2, group = group)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | Rebecca] new function print_grouped_header_table(): prints a table with a two-line header with merged group labels
# Prints 'formatted_df' (a character matrix, already formatted -- see
# print_rounded_es_value() below) with a 2-line column header, where
# columns sharing the same non-NA 'group' (from header_lines_for_columns())
# get a single label centered across their combined width on header line 1
# -- a text "merged cell" -- instead of each getting a separate line-1 cell.
# This deliberately does NOT use R's own print(): base R print() has no
# notion of a header spanning two lines or multiple columns, and computes
# its column widths/line-wrapping independently of any second header line,
# so a naive two-line header would drift out of alignment with the data
# the moment the console width changes (print() may then wrap the table
# into column blocks, and nothing keeps a second header line in sync with
# that). Here, every column width -- and the width each merged group label
# needs -- is computed and padded by hand up front, so the two header lines
# and the data stay aligned regardless of getOption("width").
print_grouped_header_table <- function(formatted_df, rn, hypo_rate_threshold = NULL,
                                       threshold_rlw = NULL) {
  ncol_ <- ncol(formatted_df)
  colnames_vec <- colnames(formatted_df)
  hdr <- header_lines_for_columns(colnames_vec, hypo_rate_threshold = hypo_rate_threshold,
                                  threshold_rlw = threshold_rlw)

  col_w <- vapply(seq_len(ncol_), function(j) {
    max(nchar(hdr$line1[j]), nchar(hdr$line2[j]), nchar(formatted_df[, j]))
  }, integer(1))

  # If a group's shared label is wider than the columns it spans (plus the
  # single-space separator between them), widen that group's last column to
  # make room, rather than truncating or overflowing the label.
  group_ids <- unique(stats::na.omit(hdr$group))
  for (g in group_ids) {
    idx <- which(hdr$group == g)
    span_w <- sum(col_w[idx]) + (length(idx) - 1)
    if (nchar(g) > span_w) {
      col_w[idx[length(idx)]] <- col_w[idx[length(idx)]] + (nchar(g) - span_w)
    }
  }

  rn_w <- max(nchar(rn), 0)
  pad_left <- function(s, w) formatC(s, width = -w)
  pad_center <- function(s, w) {
    total <- w - nchar(s)
    left <- total %/% 2
    right <- total - left
    paste0(strrep(" ", left), s, strrep(" ", right))
  }

  # Header line 1, built in column order: a merged group's columns collapse
  # into one centered label spanning their combined width; an ungrouped
  # column shows its own (possibly blank) line1 text.
  line1_parts <- character(0)
  j <- 1
  while (j <= ncol_) {
    if (!is.na(hdr$group[j])) {
      idx <- which(hdr$group == hdr$group[j])
      span_w <- sum(col_w[idx]) + (length(idx) - 1)
      line1_parts <- c(line1_parts, pad_center(hdr$group[j], span_w))
      j <- max(idx) + 1
    } else {
      line1_parts <- c(line1_parts, pad_left(hdr$line1[j], col_w[j]))
      j <- j + 1
    }
  }
  line2_parts <- mapply(pad_left, hdr$line2, col_w)

  cat(pad_left("", rn_w), " ", paste(line1_parts, collapse = " "), "\n", sep = "")
  cat(pad_left("", rn_w), " ", paste(line2_parts, collapse = " "), "\n", sep = "")
  for (i in seq_len(nrow(formatted_df))) {
    cat(pad_left(rn[i], rn_w), " ",
        paste(mapply(pad_left, formatted_df[i, ], col_w), collapse = " "), "\n", sep = "")
  }
}
# [/CHANGE 2026-10]


# [CHANGE 2026-10 | Rebecca] print_rounded_es_value(): reference-population label, thresholds, per-column formatters (overlap, median ref. pop.), grouped two-line header
# called by the benchmark.print() function. is_reference: TRUE when 'pop_es'
# is the reference population that overlap/pctl_medianRefPop are computed
# against (see determine_overlap_reference()/x$overlap_reference) -- appends
# [CHANGE 2026-10 | audit] B24: configurable reference_label
# 'reference_label' to the printed population label, e.g. "Population
# effect-size = Observed (Reference population)".
print_rounded_es_value <- function(df, pop_es, model_type, text_color, reset,
                                   is_reference = FALSE, hypo_rate_threshold = NULL,
                                   # [CHANGE 2026-10 | audit] E5/B24: arguments threshold_rlw, overlap_notes and reference_label
                                   threshold_rlw = NULL,
                                   overlap_notes = attr(df, "overlap_notes"),
                                   reference_label = "(Reference population)") {
  if (model_type == "benchmark_asymp") {
    pop_es_value <- gsub("pop_est = ", "", pop_es)
    label <- "Population estimates"
  } else {
    pop_es_value <- gsub("pop_es = ", "", pop_es)
    label <- "Population effect-size"
  }
  if (is_reference) {
    # [CHANGE 2026-10 | audit] B24: reference_label appended
    pop_es_value <- paste0(pop_es_value, " ", reference_label)
  }
  cat(sprintf("%s = %s%s%s\n", label, text_color, pop_es_value, reset))

  #formatted_column <- sprintf("%.3f", df)
  # The "Overlap with ..." column (if present -- see overlap_column() /
  # print.benchmark()) is formatted with format_overlap_value() instead of
  # format_value(), so its self-overlap entries print as "1" rather than
  # "1.000" (that value isn't actually computed, it's fixed by construction --
  # see format_overlap_value()'s own comment). Likewise the "Pctl. median
  # ref. pop." column (if present) uses format_median_ref_pop_value() so the
  # reference population's own hardcoded-50 entry prints as "50" rather than
  # "50.000". Every other column keeps using format_value() as before.
  is_overlap_col <- grepl("^Overlap with ", colnames(df))
  is_median_ref_col <- colnames(df) == "Pctl. median ref. pop."
  formatted_cols <- lapply(seq_len(ncol(df)), function(j) {
    col_formatter <- if (is_overlap_col[j]) {
      format_overlap_value
    } else if (is_median_ref_col[j]) {
      format_median_ref_pop_value
    } else {
      format_value
    }
    vapply(df[, j], col_formatter, character(1))
  })
  formatted_df <- do.call(cbind, formatted_cols)
  rownames(formatted_df) <- rownames(df)
  colnames(formatted_df) <- colnames(df)
  # [CHANGE 2026-10 | audit] E5/O10: the reason of an NA overlap is appended ("NA (<note>)")
  # An NA overlap gets its reason appended, e.g. "NA (non-finite draws; see
  # rgw_log/rlw_log)" -- 'overlap_notes' is the "notes" attribute of
  # overlap_column()'s result (one entry per row, NA where nothing to note).
  if (any(is_overlap_col) && !is.null(overlap_notes) &&
      length(overlap_notes) == nrow(formatted_df)) {
    j <- which(is_overlap_col)[1]
    has_note <- !is.na(overlap_notes) & is.na(df[, j])
    formatted_df[has_note, j] <- paste0("NA (", overlap_notes[has_note], ")")
  }
  # [/CHANGE 2026-10]
  print_grouped_header_table(formatted_df, rownames(df), hypo_rate_threshold = hypo_rate_threshold,
                             threshold_rlw = threshold_rlw)
# [/CHANGE 2026-10]
  cat("\n")
}


print_formatted_matrix <- function(mat, text_color, reset) {
  row_width <- max(nchar(rownames(mat))) + 0  
  col_widths <- apply(mat, 2, function(col) max(nchar(as.character(col))))
  col_widths <- pmax(col_widths, nchar(colnames(mat))) + 2  
  
  cat("Population Estimates (PE):\n")
  #cat(sprintf("%-10s", ""))
  cat(sprintf(paste0("%-", row_width, "s"), "")) 
  for (k in 1:ncol(mat)) {
    cat(sprintf(paste0("%", col_widths[k], "s"), colnames(mat)[k]))
  }
  cat("\n")
  
  # Print rows
  for (i in 1:nrow(mat)) {
    cat(sprintf(paste0("%-", row_width, "s"), rownames(mat)[i]))
    for (j in 1:ncol(mat)) {
      cat(sprintf(paste0("%s%", col_widths[j], "s%s"), text_color, mat[i, j], reset))
    }
    cat("\n")
  }
}

# called by the benchmark_asymp() function
check_rhs_constants <- function(rhs_list) {
  constants_check <- lapply(rhs_list, function(element) {
    return(any(element != 0))
  })
  hypotheses_with_constants <- names(constants_check)[unlist(constants_check)]
  if (length(hypotheses_with_constants) > 0) {
    warning_message <- paste0("\nrestriktor WARNING: The following hypotheses contain constants",
                              " greater or less than 0: ", 
                              paste(hypotheses_with_constants, collapse = ", "),
                              ". The default population estimates are likely incorrect.",
                              " Consider providing custom population estimates via the",
                              " pop_est argument.")
    warning(warning_message, call. = FALSE)
  }
}

# restricted least squares
theta_restricted <- function(theta, V, R, rhs) {
  correction <- V %*% t(R) %*% solve(R %*% V %*% t(R)) %*% (R %*% theta - rhs)
  as.vector(theta - correction)
}
