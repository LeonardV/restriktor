## input options
# 1. est + vcov
# 2. LL + PT
# 3. IC values (AIC, ORIC, GORIC, GORICA)
# 4. IC weights or (Bayesian) posterior model probs. 
# 5. Output from gorica function
# 6. Output from escalc (metafor package)

# In case of an equal-evidence approach, aggregating evidence from, say, 5 studies 
# with n=100 observations is the same as obtaining evidence from 1 study 
# (as if it was possible) with n=500 observations (like meta-analysis does).
# In the added-evidence approach, the aggregated evidence from, says, 5 studies 
# is stronger than as if the data were combined (as if that was possible).

# evSyn_est       <- function(object, ...) UseMethod("evSyn_est")
# evSyn_LL        <- function(object, ...) UseMethod("evSyn_LL")
# evSyn_ICvalues  <- function(object, ...) UseMethod("evSyn_ICvalues")
# [CHANGE 2026-10 | Rebecca] commented-out generic moved (order of the routes)
# evSyn_ICweights <- function(object, ...) UseMethod("evSyn_ICweights")
# evSyn_ICratios  <- function(object, ...) UseMethod("evSyn_ICratios")
# evSyn_escalc    <- function(object, ...) UseMethod("evSyn_escalc")
# -------------------------------------------------------------------------
## est (list + vec) + cov (list + mat)
## LL (list + vec) + PT (list + vec)
## IC weights (list + vec) rowSums = 1
## IC values (list + vec)

## object = est, VCOV = cov
## object = LL, PT = PT
## object = IC weights + rowsum check
## object = IC values
## object = Ratio IC weights


# [CHANGE 2026-10 | Rebecca] TO DO note on wording for type_ev = "average"
# TO DO when type_ev = "average" then we should not say 'studies' but 'analyses'.
#       So, messages and labels should be changed then (names of arguments should not).


# -------------------------------------------------------------------------
# [CHANGE 2026-10 | audit] B22: .validate_order_studies replaced by .evSyn_order_studies (study names accepted, clear error); A7/B14: new helper .evSyn_check_list_input (equal numbers of hypotheses per study, no NA/NaN)
# Helper: validate and process the 'order_studies' argument.
# Accepts a character string ("input_order", "ascending", "descending"; the
# default vector of choices is taken as "input_order"), a numeric vector
# specifying a custom study order (a permutation of 1:S), or a character
# vector with the study names in the requested order (a permutation of
# 'study_names'; mapped to the positions). Returns the character string or
# the integer order; anything else gives a clear error.
.evSyn_order_studies <- function(order_studies, S, study_names = NULL) {
  choices <- c("input_order", "ascending", "descending")
  if (is.null(order_studies)) {
    return("input_order")
  }
  if (is.numeric(order_studies)) {
    if (length(order_studies) != S) {
      stop("\nrestriktor ERROR: When 'order_studies' is a numeric vector, ",
           "its length (now, ", length(order_studies), ") must equal the number of studies (i.e., ", S, ").",
           call. = FALSE)
    }
    if (anyNA(order_studies) || !setequal(order_studies, 1:S)) {
      stop("\nrestriktor ERROR: When 'order_studies' is a numeric vector, ",
           "it must be a permutation of 1:", S, ".",
           call. = FALSE)
    }
    return(as.integer(order_studies))
  }
  if (is.character(order_studies) && !anyNA(order_studies)) {
    if (identical(order_studies, choices)) {
      return("input_order")
    }
    if (length(order_studies) == 1L && order_studies %in% choices) {
      return(order_studies)
    }
    # The study names are matched exactly (before a possible abbreviation of
    # the choices is considered), such that a study name that happens to be
    # an abbreviation of one of the choices is never taken as that choice.
    if (!is.null(study_names) && length(order_studies) == S &&
        !anyDuplicated(order_studies) &&
        setequal(order_studies, as.character(study_names))) {
      return(match(order_studies, as.character(study_names)))
    }
    if (length(order_studies) == 1L) {
      m <- pmatch(order_studies, choices)
      if (!is.na(m)) {
        return(choices[m])
      }
    }
  }
  stop("\nrestriktor ERROR: The argument 'order_studies' must be one of 'input_order', ",
       "'ascending', 'descending', a numeric vector with a permutation of 1:", S,
       ", or a character vector with a permutation of the study names in 'study_names'",
       if (is.null(study_names)) " (which is not specified)",
       ". Now, it is: ", paste(deparse(order_studies, nlines = 1), collapse = ""), ".",
       call. = FALSE)
}

# Helper: check the elements of a list input ('object', or 'PT'): each element
# must be a numeric vector without NA/NaN (infinite values are allowed) and
# all elements must have the same length, i.e., the number of hypotheses must
# be identical across studies (no silent recycling in rbind()). When 'ref'
# (another such list, e.g., 'object' when checking 'PT') is given, the lengths
# must also match those of 'ref', per study.
.evSyn_check_list_input <- function(x, name = "object", ref = NULL, ref_name = "object") {
  if (!is.list(x) || length(x) == 0 || !all(vapply(x, is.numeric, logical(1)))) {
    stop("\nrestriktor ERROR: The argument '", name, "' must be a list of numeric vectors ",
         "(one vector for each study).", call. = FALSE)
  }
  has_na <- vapply(x, anyNA, logical(1))
  if (any(has_na)) {
    stop("\nrestriktor ERROR: The argument '", name, "' contains NA or NaN values ",
         "(study ", paste(which(has_na), collapse = ", "), "). ",
         "Please remove or replace these values (infinite values are allowed).",
         call. = FALSE)
  }
  len <- lengths(x)
  if (length(unique(len)) > 1L) {
    stop("\nrestriktor ERROR: The number of hypotheses must be identical across studies. ",
         "The elements of '", name, "' have lengths (", paste(len, collapse = ", "), ").",
         call. = FALSE)
  }
  if (!is.null(ref)) {
    if (length(x) != length(ref)) {
      stop("\nrestriktor ERROR: The number of elements in '", name, "' (", length(x),
           ") must equal the number of elements in '", ref_name, "' (", length(ref),
           "), i.e., the number of studies.", call. = FALSE)
    }
    if (!all(len == lengths(ref))) {
      stop("\nrestriktor ERROR: The number of values in '", name, "' must match the ",
           "number of values in '", ref_name, "' for each study. Found lengths (",
           paste(len, collapse = ", "), ") versus (", paste(lengths(ref), collapse = ", "), ").",
           call. = FALSE)
    }
  }
  invisible(x)
}

# Helper: match a named weight vector (priorICweights, study_weights) to the
# names of the hypotheses or studies. When the weights carry names, these
# must be exactly the names in 'ref' (possibly in another order; the weights
# are then re-ordered accordingly) or, when 'alias' is given (e.g., the
# hypothesis names of the input vectors, which correspond one-to-one to the
# labels in 'ref'), exactly the names in 'alias'. An error is given when the
# names do not match. Returns the (re-ordered) weights without names; unnamed
# weights are returned as is (they are matched by position).
# Note: 'ref' may be numeric (e.g., study names that are years); the weights
# are indexed by the names as character strings (not by position).
.evSyn_match_weight_names <- function(w, ref, name = "priorICweights", what = "hypotheses",
                                      alias = NULL) {
  nms <- names(w)
  if (is.null(nms)) {
    return(w)
  }
  if (is.null(ref) || length(w) != length(ref)) {
    # the length is checked (and reported) by check_weights()
    return(unname(w))
  }
  ref <- as.character(ref)
  if (!is.null(alias) && (length(alias) != length(ref) || anyDuplicated(alias))) {
    alias <- NULL
  }
  if (!is.null(alias)) {
    alias <- as.character(alias)
  }
  # the failsafe labels 'Complement'/'Unconstrained' may be given in any case
  # (goric() uses 'complement'/'unconstrained')
  norm_fs <- function(x) {
    is_fs <- tolower(x) %in% c("complement", "unconstrained")
    x[is_fs] <- tolower(x[is_fs])
    x
  }
  nms_n <- norm_fs(nms)
  ref_n <- norm_fs(ref)
  if (anyNA(nms) || any(nms == "") || anyDuplicated(nms_n)) {
    ok <- FALSE
  } else if (setequal(nms_n, ref_n)) {
    w <- w[match(ref_n, nms_n)]
    names(w) <- ref
    ok <- TRUE
  } else if (!is.null(alias) && setequal(nms, alias)) {
    # the names of the input are used (in any order); map them to the labels
    w <- w[alias]
    names(w) <- ref
    ok <- TRUE
  } else {
    ok <- FALSE
  }
  if (!ok) {
    stop("\nrestriktor ERROR: The names of '", name, "' (", paste(nms, collapse = ", "),
         ") do not match the names of the ", what, " (", paste(ref, collapse = ", "),
         if (!is.null(alias) && !identical(alias, ref)) {
           paste0("; or the names of the input, i.e., ", paste(alias, collapse = ", "))
         }, "). ",
         "Please use exactly these names (in any order) or an unnamed vector ",
         "(which is matched by position).", call. = FALSE)
  }
  unname(w[ref])
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B15: new helper .evSyn_priorWeights_compat (deprecated 'priorWeights' -> 'priorICweights' with a warning)
# -------------------------------------------------------------------------
# Helper: backwards compatibility for the deprecated argument 'priorWeights'
# (renamed to 'priorICweights'). 'dots' is the list of arguments passed via '...'.
.evSyn_priorWeights_compat <- function(priorICweights, dots) {
  priorWeights <- dots[["priorWeights"]]
  if (is.null(priorWeights)) {
    return(priorICweights)
  }
  if (!is.null(priorICweights)) {
    stop("\nrestriktor ERROR: both 'priorICweights' and the deprecated argument 'priorWeights' are specified; ",
         "please use 'priorICweights' only.", call. = FALSE)
  }
  warning("restriktor WARNING: argument 'priorWeights' is deprecated; use 'priorICweights'",
          call. = FALSE)
  priorWeights
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] new helper .evSyn_cum_weighted: study-weighted cumulative sums; a study with weight 0 contributes exactly 0 (no 0 * Inf)
# Helper: cumulative (study-)weighted sums of the rows of M.
# Row s equals sum_{i <= s} w_i * M[i, ] / mean(w_1, ..., w_s), with w =
# study_weights_S (which sum to S). This equals the weighting used for the
# cumulative IC values; with equal study weights it reduces to cumsum.
# A study with weight 0 contributes nothing: it is not counted, so row s then 
# equals row s-1 (the mean is taken over the studies with a positive weight, 
# see .evSyn_n_pos). If the weights of studies 1, ..., s are all 0, there is 
# no evidence yet and row s is set to 0 (so that the corresponding IC weights 
# equal the prior IC weights). The contribution of a study with weight 0 is 
# set to exactly 0 (and not computed as 0 * value), such that infinite values 
# (e.g., an IC difference of Inf when an IC weight is 0) of such a study do 
# not result in NaN (0 * Inf) in the cumulative results.
.evSyn_cum_weighted <- function(M, study_weights_S) {
  M <- M[, , drop = FALSE]
  S <- nrow(M)
  Mw <- M * study_weights_S
  Mw[study_weights_S == 0, ] <- 0
  out <- apply(Mw, 2, cumsum)
  out <- matrix(out, nrow = S, dimnames = dimnames(M))
  denom <- cumsum(study_weights_S) / pmax(.evSyn_n_pos(study_weights_S), 1)
  out <- out / denom
  out[denom == 0, ] <- 0
  out
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] new helper .evSyn_n_pos: cumulative number of positively weighted studies
# Helper: cumulative number of studies with a positive study weight (equals 
# 1, ..., S when all study weights are positive). Used instead of the number
# of studies when averaging, such that a study with weight 0 is not counted.
.evSyn_n_pos <- function(study_weights_S) {
  cumsum(study_weights_S > 0)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B21: new helper .evSyn_check_study_weights (central validation via check_weights, names matched to the studies, zero weights allowed)
# Helper: validate 'study_weights' (length S, non-negative, finite, at least
# one positive; zero weights are allowed). Returns a list with the study
# weights (summing to 1 or to S, as before) and the study weights summing to S.
# When the study weights carry names, these must be the study names (in
# 'study_names', or "1", ..., "S" when not specified; in any order) and the
# weights are matched to the studies by name (see .evSyn_match_weight_names).
.evSyn_check_study_weights <- function(study_weights, S, study_names = NULL) {
  if (is.null(study_weights)) {
    study_weights <- rep(1/S, S)
  } else {
    if (is.null(study_names) || length(study_names) != S) {
      study_names <- as.character(seq_len(S))
    }
    study_weights <- .evSyn_match_weight_names(study_weights, study_names,
                                               name = "study_weights",
                                               what = "studies (see 'study_names')")
    study_weights <- check_weights(study_weights, name = "study_weights",
                                   length_expected = S,
                                   what = "one for each study",
                                   rescale = FALSE)
    if (!isTRUE(all.equal(sum(study_weights), 1)) &&
        !isTRUE(all.equal(sum(study_weights), S))) {
      study_weights <- study_weights / sum(study_weights)
      message("\nrestriktor Message: The argument 'study_weights' should add up to 1 or to S = ", S, ". It has been rescaled accordingly.")
    }
  }
  list(study_weights = study_weights,
       study_weights_S = S * study_weights / sum(study_weights)) # Now, they sum up to S
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B2/B21: new helper .evSyn_check_priorICweights (central validation, names matched to the hypotheses or the input names)
# Helper: validate 'priorICweights' (length NrHypos_incl, non-negative,
# finite, at least one positive; rescaled such that they sum to 1).
# When the prior weights carry names, these must be the hypothesis names
# (i.e., the column names of the output, in any order) and the prior weights
# are matched to the hypotheses by name (see .evSyn_match_weight_names).
# When the input vectors (or hypothesis sets) carry hypothesis names and these
# differ from the labels in 'hnames' (because 'hypo_names' is given), the
# names of the input ('input_names', in the order of 'hnames') are accepted
# as well.
.evSyn_check_priorICweights <- function(priorICweights, NrHypos_incl, hnames = NULL,
                                        input_names = NULL) {
  if (is.null(priorICweights)) {
    return(rep(1/NrHypos_incl, NrHypos_incl))
  }
  priorICweights <- .evSyn_match_weight_names(priorICweights, hnames,
                                              name = "priorICweights",
                                              what = "hypotheses (i.e., the column names of the output)",
                                              alias = input_names)
  if (is.numeric(priorICweights) && !anyNA(priorICweights) &&
      all(is.finite(priorICweights)) &&
      !isTRUE(all.equal(sum(priorICweights), 1))) {
    message("\nrestriktor Message: The argument 'priorICweights' should add up to 1. It has been rescaled accordingly.")
  }
  check_weights(priorICweights, name = "priorICweights",
                length_expected = NrHypos_incl,
                what = "one for each hypothesis including a possible failsafe hypothesis",
                rescale = TRUE)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] new helper .evSyn_rescale_study_weights: study weights as powers of IC weights (added: sum to #positive studies; average: sum to 1)
# Helper: rescale the study weights of studies 1, ..., s as used when combining
# IC weights (i.e., when the IC weights are 'powers'; see evSyn_ICweights): 
# - added:   summing to the number of positively weighted studies (as if IC 
#            values were summed), 
# - average: summing to 1 (as if IC values were averaged).
# If all study weights are 0, all rescaled weights are 0 (no evidence yet).
# This is used both for the cumulative/final results and for determining the 
# preferred hypothesis when ordering the studies, such that these agree.
.evSyn_rescale_study_weights <- function(study_weights_S, type_ev) {
  if (sum(study_weights_S) == 0) {
    return(rep(0, length(study_weights_S)))
  }
  w <- study_weights_S / sum(study_weights_S)
  if (type_ev == "average") {
    w
  } else {
    sum(study_weights_S > 0) * w
  }
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] new helper .evSyn_log_prod_weights: products of IC weights on the log scale (no underflow)
# Helper: (study-)weighted sum of log IC weights over studies, i.e., the log of
# prod_i W[i, ]^expo[i]. Computed on the log scale to avoid underflow when
# many studies are combined. A study with exponent 0 contributes nothing
# (also when one of its weights is 0, i.e., log(0) = -Inf).
.evSyn_log_prod_weights <- function(logW, expo) {
  logW <- logW[, , drop = FALSE]
  terms <- logW * expo
  terms[expo == 0, ] <- 0
  colSums(terms)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B16: new helper .evSyn_cum_IC: cumulative IC values (added/equal/average) with study weights and penalty_factor
# Helper: cumulative IC values based on the log-likelihood and penalty values 
# (possibly weighted using the study weights, see .evSyn_cum_weighted), for 
# the added-, equal-, and average-evidence approach. Row s corresponds to 
# studies 1, ..., s.
.evSyn_cum_IC <- function(LL_m, PT_m, study_weights_S, type_ev, penalty_factor = 2) {
  S <- nrow(LL_m)
  cumLL <- .evSyn_cum_weighted(LL_m, study_weights_S)
  cumPT <- .evSyn_cum_weighted(PT_m, study_weights_S)
  # number of (positively weighted) studies so far; used for averaging
  s <- pmax(.evSyn_n_pos(study_weights_S), 1)
  switch(type_ev,
         added   = -2 * cumLL + penalty_factor * cumPT,
         equal   = -2 * cumLL + penalty_factor * cumPT / s,
         average = (-2 * cumLL + penalty_factor * cumPT) / s)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B14: new helper .evSyn_IC_weights_rows: (prior-weighted) IC weights per row via log-sum-exp (ic_weights_log)
# Helper: (prior-weighted) IC weights for each row of a matrix with IC values,
# computed on the log scale (see ic_weights_log).
.evSyn_IC_weights_rows <- function(IC_m, priorICweights = NULL) {
  IC_m <- IC_m[, , drop = FALSE]
  out <- IC_m
  for (s in seq_len(nrow(IC_m))) {
    out[s, ] <- ic_weights_log(IC_m[s, ], priorICweights)
  }
  out
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] new helper .evSyn_pref_hypo: preferred hypothesis = highest final prior-weighted IC weight
# Helper: overall preferred hypothesis, i.e., the hypothesis with the highest
# final (prior-weighted) IC weight. 'IC' is the vector with final IC values
# (or IC differences); the first hypothesis is taken in case of ties.
.evSyn_pref_hypo <- function(IC, priorICweights) {
  w <- ic_weights_log(IC, priorICweights)
  which.max(w)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] new helper .evSyn_align_goric_hypos: goric objects aligned by hypothesis name (and text), not by position
# Helper: align a list of goric objects to the hypothesis set of the first 
# object. The hypotheses are identified by their names (x$result$model, i.e., 
# the user-specified names or H1, H2, ... when unnamed, plus the possible 
# failsafe hypothesis 'unconstrained' or 'complement'). When the hypotheses of
# a study are given in another order, its results are re-ordered to the order 
# of study 1. If the names differ across studies, an error is given. 
# Additionally, when the hypothesis text (x$hypotheses_usr) is available, it 
# is checked whether a hypothesis with the same name has the same text: 
# - for unnamed hypotheses (default names H1, H2, ...) differing text is 
#   ambiguous and an error is given (name the hypotheses in goric() instead),
# - for user-named hypotheses a warning is given (the text may differ, e.g.,
#   because of study-specific parameter names).
.evSyn_align_goric_hypos <- function(object) {
  hypo_text <- function(x) {
    h <- x$hypotheses_usr
    if (is.null(h) || !is.list(h) || !all(vapply(h, is.character, logical(1)))) {
      return(NULL)
    }
    gsub("[[:space:]]+", "", vapply(h, function(z) paste(z, collapse = ";"), character(1)))
  }
  hypo_named <- function(x) {
    h <- x$hypotheses_usr
    !is.null(h) && !is.null(names(h)) && all(names(h) != "")
  }
  nms <- lapply(object, function(x) as.character(x$result$model))
  ref <- nms[[1]]
  if (anyDuplicated(ref)) {
    stop("\nrestriktor ERROR: The hypothesis names of a goric object must be unique. ",
         "Found: ", paste(sQuote(ref), collapse = ", "), ".", call. = FALSE)
  }
  txt_ref <- hypo_text(object[[1]])
  for (s in seq_along(object)) {
    if (length(nms[[s]]) != length(ref) || !setequal(nms[[s]], ref)) {
      stop("\nrestriktor ERROR: The hypotheses (names) must be identical across the goric objects. ",
           "Study 1 uses (", paste(ref, collapse = ", "), "), while study ", s, 
           " uses (", paste(nms[[s]], collapse = ", "), "). ",
           "Please use the same (named) hypotheses in goric() for all studies.",
           call. = FALSE)
    }
    idx <- match(ref, nms[[s]])
    if (!identical(idx, seq_along(ref))) {
      message("\nrestriktor Message: The hypotheses of study ", s, " are given in a ", 
              "different order than those of study 1. The hypotheses are matched by name.")
      object[[s]]$result <- object[[s]]$result[idx, , drop = FALSE]
      rownames(object[[s]]$result) <- NULL
    }
    # Check the hypothesis text, if available (user-specified hypotheses only, 
    # i.e., excluding the failsafe hypothesis)
    txt_s <- hypo_text(object[[s]])
    if (s > 1 && !is.null(txt_ref) && !is.null(txt_s) && 
        length(txt_ref) == length(txt_s)) {
      idx_usr <- idx[idx <= length(txt_s)]
      txt_s <- txt_s[idx_usr]
      differs <- txt_s != txt_ref
      if (any(differs)) {
        hn <- ref[seq_along(txt_ref)][differs]
        msg <- paste0("The hypothesis ", paste(sQuote(hn), collapse = ", "), 
                      " differs between study 1 (", paste(sQuote(txt_ref[differs]), collapse = ", "), 
                      ") and study ", s, " (", paste(sQuote(txt_s[differs]), collapse = ", "), ").")
        if (hypo_named(object[[1]]) && hypo_named(object[[s]])) {
          warning("\nrestriktor WARNING: ", msg, 
                  " Since the hypotheses are named, they are assumed to represent the same theory.",
                  call. = FALSE)
        } else {
          stop("\nrestriktor ERROR: ", msg, 
               " The hypotheses are unnamed (H1, H2, ...), so it is ambiguous which hypotheses ",
               "should be combined. Please name the hypotheses in goric() identically across studies.",
               call. = FALSE)
        }
      }
    }
  }
  object
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B13: new helper .evSyn_check_input_names: named input vectors checked and re-ordered by name across studies; hypo_names mismatch/permutation warned
# Helper: check the hypothesis names of the input vectors (object), if any.
# If all studies carry names, these must denote the same set of hypotheses;
# when the order differs from study 1, the vectors are re-ordered to the order
# of study 1. If only some studies carry names, the positions are used (with a
# warning). Returns the (possibly re-ordered) object.
.evSyn_check_input_names <- function(object, hypo_names = NULL) {
  nms <- lapply(object, names)
  has_names <- vapply(nms, function(x) !is.null(x) && all(!is.na(x)) && all(x != ""), logical(1))
  if (!any(has_names)) {
    return(object)
  }
  if (!all(has_names)) {
    warning("\nrestriktor WARNING: The hypothesis names are missing for some of the studies. ",
            "The hypotheses are matched across studies by position.",
            call. = FALSE)
    return(object)
  }
  ref <- nms[[1]]
  if (anyDuplicated(ref)) {
    stop("\nrestriktor ERROR: The hypothesis names of a study must be unique. ",
         "Found: ", paste(sQuote(ref), collapse = ", "), ".", call. = FALSE)
  }
  for (s in seq_along(object)) {
    if (!setequal(nms[[s]], ref) || length(nms[[s]]) != length(ref)) {
      stop("\nrestriktor ERROR: The hypothesis names must be identical across studies. ",
           "Study 1 uses (", paste(ref, collapse = ", "), "), while study ", s,
           " uses (", paste(nms[[s]], collapse = ", "), ").",
           call. = FALSE)
    }
    if (!identical(nms[[s]], ref)) {
      message("\nrestriktor Message: The hypotheses of study ", s, " are given in a ",
              "different order than those of study 1. The hypotheses are matched by name.")
      object[[s]] <- object[[s]][ref]
    }
  }
  if (!is.null(hypo_names) && length(hypo_names) == length(ref)) {
    if (!setequal(hypo_names, ref)) {
      warning("\nrestriktor WARNING: The names in 'hypo_names' (", paste(hypo_names, collapse = ", "),
              ") differ from the hypothesis names of the input (", paste(ref, collapse = ", "),
              "). The names in 'hypo_names' are used as labels, in the order of the input.",
              call. = FALSE)
    } else if (!identical(as.character(hypo_names), ref)) {
      .evSyn_warn_hypo_names_permuted(hypo_names, ref)
    }
  }
  object
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] new helper .evSyn_input_hypo_names: input names accepted as alias for a named priorICweights
# Helper: the hypothesis names of the input vectors (after
# .evSyn_check_input_names, i.e., identical across studies), or NULL when
# (some of) the input vectors are unnamed. These names are accepted as an
# alias of the labels (e.g., for a named 'priorICweights').
.evSyn_input_hypo_names <- function(object) {
  nms <- lapply(object, names)
  has_names <- vapply(nms, function(x) !is.null(x) && all(!is.na(x)) && all(x != ""), logical(1))
  if (!all(has_names) || anyDuplicated(nms[[1]])) {
    return(NULL)
  }
  nms[[1]]
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B13: new helper .evSyn_warn_hypo_names_permuted
# Helper: warning when 'hypo_names' consists of the hypothesis names of the
# input, but in another order. The names in 'hypo_names' are labels that are
# applied in the order of the input (the hypotheses are not re-ordered), so
# this is most likely not what the user intends.
.evSyn_warn_hypo_names_permuted <- function(hypo_names, ref) {
  warning("\nrestriktor WARNING: The names in 'hypo_names' (", paste(hypo_names, collapse = ", "),
          ") are the hypothesis names of the input (", paste(ref, collapse = ", "),
          ") in another order. The names in 'hypo_names' are used as labels, in the ",
          "order of the input; the hypotheses are NOT re-ordered. If the hypotheses ",
          "should be labelled by their own names, use 'hypo_names' in the order of the ",
          "input (or leave 'hypo_names' unspecified).",
          call. = FALSE)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] new helper .evSyn_align_second_input: PT aligned to the (re-ordered) LL values per study by name
# Helper: align a second, parallel input (e.g., the penalty values 'PT' that 
# go with the log-likelihood values in 'object') to the (possibly re-ordered, 
# see .evSyn_check_input_names) 'object', per study. For each study, the names
# of PT[[s]] must denote the same hypotheses as the names of object[[s]]: when
# the order differs, PT[[s]] is re-ordered to the order of object[[s]]; when
# the sets of names differ, or when only one of the two is named, an error is
# given. Unnamed vectors (in both) are matched by position.
.evSyn_align_second_input <- function(object, second, name = "PT", 
                                      object_name = "object") {
  is_named <- function(x) {
    nms <- names(x)
    !is.null(nms) && all(!is.na(nms)) && all(nms != "")
  }
  for (s in seq_along(object)) {
    if (s > length(second)) {
      break
    }
    obj_named <- is_named(object[[s]])
    sec_named <- is_named(second[[s]])
    if (!obj_named && !sec_named) {
      next
    }
    if (obj_named != sec_named) {
      stop("\nrestriktor ERROR: For study ", s, ", the hypothesis names are given for '",
           if (obj_named) object_name else name, "' but not for '", 
           if (obj_named) name else object_name, "'. ",
           "Please name the elements of both (identically) or of neither.",
           call. = FALSE)
    }
    ref <- names(object[[s]])
    nms <- names(second[[s]])
    if (length(nms) != length(ref) || !setequal(nms, ref) || anyDuplicated(nms)) {
      stop("\nrestriktor ERROR: For study ", s, ", the hypothesis names of '", name, 
           "' (", paste(nms, collapse = ", "), ") must be identical to those of '", 
           object_name, "' (", paste(ref, collapse = ", "), ").",
           call. = FALSE)
    }
    if (!identical(nms, ref)) {
      second[[s]] <- second[[s]][ref]
    }
  }
  second
}
# [/CHANGE 2026-10]

# -------------------------------------------------------------------------
evSyn <- function(object, input_type = NULL, ...) {

  args <- list(...)

  # [CHANGE 2026-10 | audit] B15: deprecated 'priorWeights' mapped to 'priorICweights' in the dispatcher
  # Backwards compatibility: 'priorWeights' is renamed to 'priorICweights'.
  if (!is.null(args[["priorWeights"]])) {
    args$priorICweights <- .evSyn_priorWeights_compat(args[["priorICweights"]], args)
    args$priorWeights <- NULL
  }
  # [/CHANGE 2026-10]

  VCOV <- args$VCOV
  PT   <- args$PT
  
  type_ev <- args$type_ev
  
  call_sub <- function(fun, args, object) {
    do.call(fun, c(list(object), args))
  }
  
  isGoric <- if (is.list(object)) {
    vapply(object, function(x) inherits(x, "con_goric"), logical(1))
  } else {
    FALSE
  }
  
  # Handle input_type explicitly if provided
  if (!is.null(input_type)) {
    it <- tolower(input_type)
    if (it == "est_vcov") {
      return(call_sub(evSyn_est, args, object))
    } else if (it == "ll_pt") {
      return(call_sub(evSyn_LL, args, object))
    } else if (it == "icweights") {
      return(call_sub(evSyn_ICweights, args, object))
    } else if (it == "icratios") {
      return(call_sub(evSyn_ICratios, args, object))
    } else if (it == "icvalues") {
      return(call_sub(evSyn_ICvalues, args, object))
    } else if (it %in% c("goric", "gorica", "goricc", "goricac")) {
      return(call_sub(evSyn_gorica, args, object))
    } else if (it == "escalc") {
      return(call_sub(evSyn_escalc, args, object))
    } else {
      stop(paste0("\nrestriktor ERROR: Unknown input_type ", sQuote(input_type), "."))
    }
  }
  
  if (all(isGoric)) {
    return(call_sub(evSyn_gorica, args, object))
  }
  
  if (any(inherits(object, c("escalc", "data.frame")))) {
    return(call_sub(evSyn_escalc, args, object))
  } 
  
  # [CHANGE 2026-10 | audit] empty list refused; all elements must be numeric (was any())
  if (!is.list(object) || length(object) == 0 ||
      !all(vapply(object, is.numeric, logical(1)))) {
    stop("\nrestriktor ERROR: object must be a list of numeric vectors.", call. = FALSE)
  }
  
  if (!is.null(VCOV) && !is.null(PT)) {
    stop("\nrestriktor ERROR: both VCOV and PT are found, which confuses me.", call. = FALSE)
  }
  
  if (!is.null(VCOV)) {
    # [CHANGE 2026-10 | audit] comment
    # estimates (the number of estimates may differ across studies)
    return(call_sub(evSyn_est, args, object))
  # [CHANGE 2026-10 | audit] A7/B14/A8: list input validated (.evSyn_check_list_input, also PT); input type inferred by .evSyn_detect_input_type with a message (replaces inline sum-to-1 / last-element-1 detection and per-route 'equal' messages)
  }
  
  # From here on, the input is a list of vectors with log-likelihood values,
  # IC values, IC weights, or ratios of IC weights: one value for each
  # hypothesis, the same hypotheses in each study, no NA/NaN.
  .evSyn_check_list_input(object, name = "object")

  if (!is.null(PT)) {
    .evSyn_check_list_input(PT, name = "PT", ref = object, ref_name = "object")
    return(call_sub(evSyn_LL, args, object))
  }

  # Infer the input type from the values (IC weights, ratios of IC weights,
  # or IC values; see .evSyn_detect_input_type) and say so. Note: when the
  # equal-evidence approach is requested for these input types, the route
  # functions fall back to the added-evidence approach (with a message).
  detected <- .evSyn_detect_input_type(object)
  message("\nrestriktor Message: ", detected$msg,
          " Specify 'input_type' (i.e., 'icvalues', 'icweights', or 'icratios') ",
          "to override the inferred input type.")
  fun <- switch(detected$type,
                icweights = evSyn_ICweights,
                icratios  = evSyn_ICratios,
                icvalues  = evSyn_ICvalues)
  return(call_sub(fun, args, object))
}
  # [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] A8/B20: new helper .evSyn_detect_input_type (ratios only if all > 0 and a common 1-position; ambiguous input requires input_type)
# Helper: infer the input type of a list of numeric vectors (one for each
# study), when 'input_type' is not specified:
# - 'icweights': all values lie between 0 and 1 and the values of each study
#                sum to 1. When the values of each study sum to 1 only
#                approximately (within 1e-3; e.g., rounded weights), the input
#                is ambiguous and an error is given (specify 'input_type').
# - 'icratios':  all values are positive and there is a hypothesis (position)
#                with a value of exactly 1 in each study (the reference
#                hypothesis). When each study contains a value of 1, but the
#                values are not all positive or the 1 is not at the same
#                position in each study, the input is ambiguous (IC values
#                may equal 1 as well) and an error is given.
# - 'icvalues':  otherwise.
# Returns a list with the inferred type and a message explaining the choice.
.evSyn_detect_input_type <- function(object) {
  tol_sum <- 1e-3
  in01 <- all(vapply(object, function(x) all(x >= 0 & x <= 1), logical(1)))
  sums <- vapply(object, sum, numeric(1))
  if (in01 && all(abs(sums - 1) <= sqrt(.Machine$double.eps))) {
    return(list(type = "icweights",
                msg = paste0("The input is treated as IC weights (input_type = 'icweights'), ",
                             "since the values of each study lie between 0 and 1 and sum to 1.")))
  }
  if (in01 && all(abs(sums - 1) <= tol_sum)) {
    stop("\nrestriktor ERROR: The input type is ambiguous: the values of each study lie ",
         "between 0 and 1 and sum approximately (but not exactly) to 1, as rounded IC ",
         "weights do. Please specify 'input_type': 'icweights' (after rescaling the ",
         "values of each study such that they sum to 1), 'icvalues', or 'icratios'.",
         call. = FALSE)
  }
  has_one <- vapply(object, function(x) any(x == 1), logical(1))
  if (all(has_one)) {
    all_pos <- all(vapply(object, function(x) all(x > 0), logical(1)))
    common  <- Reduce(intersect, lapply(object, function(x) which(x == 1)))
    if (all_pos && length(common) > 0) {
      return(list(type = "icratios",
                  msg = paste0("The input is treated as ratios of IC weights (input_type = 'icratios'), ",
                               "since all values are positive and hypothesis ", common[1],
                               " has a value of exactly 1 in each study (the reference hypothesis). ",
                               "Note that IC values may equal 1 as well; if the input consists of ",
                               "IC values, specify input_type = 'icvalues'.")))
    }
    stop("\nrestriktor ERROR: The input type is ambiguous: each study contains a value of ",
         "exactly 1 (as ratios of IC weights have for the reference hypothesis), but ",
         if (!all_pos) "not all values are positive" else
           "the value of 1 is not at the same position (hypothesis) in each study",
         ". Please specify 'input_type': 'icratios' (ratios of IC weights) or 'icvalues' (IC values).",
         call. = FALSE)
  }
  list(type = "icvalues",
       msg = "The input is treated as IC values (input_type = 'icvalues').")
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B23: new helper .evSyn_type_ev_IC ('equal' falls back to 'added' with a message instead of a match.arg error)
# Helper: for input consisting of IC values, IC weights, or ratios of IC
# weights, the equal-evidence approach is not possible (there are no separate
# log-likelihood and penalty values); the added-evidence approach is used
# instead (with a message), as documented. Returns the matched 'type_ev'.
.evSyn_type_ev_IC <- function(type_ev, what) {
  # also an abbreviation of "equal" (e.g., "eq"), as accepted by match.arg()
  is_equal <- is.character(type_ev) && length(type_ev) == 1L && !is.na(type_ev) &&
    identical(pmatch(type_ev, c("added", "equal", "average")), 2L)
  if (is_equal) {
    message("\nrestriktor Message: When the input consists of ", what,
            ", the equal-evidence approach is not possible. The added-evidence approach is used instead.")
    type_ev <- "added"
  }
  match.arg(type_ev, c("added", "average"))
}
# [/CHANGE 2026-10]


# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on the (standardized) parameter estimates and 
# the covariance matrix
evSyn_est <- function(object, ..., VCOV = list(), hypotheses = list(),
                      type_ev = c("added", "equal", "average"), 
                      comparison = c("unconstrained", "complement", "none"),
                      # [CHANGE 2026-10 | Rebecca] priorICweights argument
                      hypo_names = c(), priorICweights = NULL,
                      type = c("gorica", "goricac"),
                      order_studies = c("input_order", "ascending", "descending"),
                      study_names = c(),
                      # [CHANGE 2026-10 | Rebecca] study_weights argument
                      study_sample_nobs = NULL,
                      study_weights = NULL) {
  
  # [CHANGE 2026-10 | audit] B15: priorWeights compat (removed from dots passed to goric()); Heq: a hypothesis named 'Heq' is dropped and regenerated by goric()
  # Backwards compatibility: 'priorWeights' is renamed to 'priorICweights'.
  # The deprecated argument is removed from the arguments that are passed on
  # to goric() (which does not know it).
  dots <- list(...)
  priorICweights <- .evSyn_priorWeights_compat(priorICweights, dots)
  dots$priorWeights <- NULL
  
  # 'Heq' (passed on to goric()): the equality-restricted version of the
  # (single) order-restricted hypothesis is added to the set by goric() when
  # it is compared to its complement; the result then consists of Heq, the
  # hypothesis, and its complement. As in goric(), a hypothesis named 'Heq'
  # (e.g., from a benchmark) is removed from the set(s) and regenerated.
  Heq <- isTRUE(dots[["Heq"]])
  if (Heq) {
    drop_Heq <- function(h) {
      if (is.list(h) && !is.null(names(h))) h[names(h) != "Heq"] else h
    }
    if (is.list(hypotheses) && length(hypotheses) > 0 &&
        all(vapply(hypotheses, is.list, logical(1)))) {
      hypotheses <- lapply(hypotheses, drop_Heq)
    } else {
      hypotheses <- drop_Heq(hypotheses)
    }
  }
  # [/CHANGE 2026-10]
  
  # [CHANGE 2026-10 | audit] regression fix: missing(comparison) evaluated before 'comparison' is assigned
  # Note: missing() must be evaluated before any assignment to 'comparison'.
  comparison_missing <- missing(comparison)
  if (comparison_missing) {
    if (length(hypotheses) == 1) {
      comparison <- "complement"
    } else {
      comparison <- "unconstrained"
      # [CHANGE 2026-10 | Rebecca] note on the complement default for lists of single hypotheses
      # Note: in the case of a list in a list,
      # and only one hypothesis in the sub-lists, then:
      # comparison <- "complement"
      # which is done below.
      # [/CHANGE 2026-10]
    }
  }
  
  comparison <- match.arg(comparison)
  #type <- match.arg(type) 
  # I want to all for c("goric", "goricc", "gorica", "goricac"), which
  # may be overwritten next, so these should not be shown as the options
  
  if (missing(type)) { 
    type <- "gorica" 
    type_missing <- TRUE
  } else if (!is.null(type) && type %in% c("goric", "goricc", "gorica", "goricac")) {
    type_missing <- FALSE
  } else if (is.null(type)) {
    type <- "gorica"
    type_missing <- TRUE
  } else {
    message(paste0("\nrestriktor Message: The value for the argument type (i.e., '", type, "') is not valid. \n", 
            "Since the input is a list of estimates, the GORICA will be used."))
    type <- "gorica"
    type_missing <- NULL
  }
  #
  if (type == "goric") {
    message("\nrestriktor Message: Since the input is a list of estimates, the GORICA will be used, not the GORIC.")
    type = "gorica"
  } else if (type == "goricc") {
    message("\nrestriktor Message: Since the input is a list of estimates, the GORICAC will be used, not the GORICC.")
    type = "goricac"
  } 
  #
  #if (missing(study_sample_nobs)) { study_sample_nobs <- NULL } 
  if (type == "goricac" && is.null(study_sample_nobs)) {
    stop("\nrestriktor ERROR: To compute the GORICAC, the argument 'study_sample_nobs' is required. ",
         "Please provide a numeric vector with the sample sizes of all primary studies.",
         call. = FALSE)
  }
  type <- match.arg(type)
  
  if (missing(type_ev)) 
    type_ev <- "added"
  type_ev <- match.arg(type_ev)
  
  # number of primary studies
  S <- length(object)
  V <- length(VCOV)
  
  if (S != V) {
    stop("\nrestriktor ERROR: The number of elements in the 'object' list (i.e., the number of (standardized) estimates) must match the number of elements in the 'VCOV' list.",
         call. = FALSE)
  }
  
  if (type == "goricac" && length(study_sample_nobs) == 1) {
    study_sample_nobs <- rep(study_sample_nobs, S)
    message("\nrestriktor Message: The argument 'study_sample_nobs' contains a single value; all primary studies are assumed to have the same sample size.")
  } else if (type == "goricac" && length(study_sample_nobs) != S) {
    stop("\nrestriktor ERROR: The argument 'study_sample_nobs' must be a numeric vector containing S = ", S, " values (one for each study). \n",
         "Alternatively, 'study_sample_nobs' can be a scalar if all studies have the same sample size.",
         call. = FALSE)
  }
  
  # [CHANGE 2026-10 | audit] B21/B22: study_weights and order_studies validated via the helpers (replaces the match.arg block)
  # Check the study weights (zero weights are allowed; see .evSyn_cum_weighted)
  study_weights <- .evSyn_check_study_weights(study_weights, S, study_names)
  study_weights_S <- study_weights$study_weights_S # Now, they sum up to S
  study_weights <- study_weights$study_weights

  # Check the order of the studies (character string, permutation of 1:S, or
  # permutation of the study names)
  order_studies <- .evSyn_order_studies(order_studies, S, study_names)
  # [/CHANGE 2026-10]
  
  # Ensure hypotheses are nested
  if (!all(vapply(hypotheses, is.list, logical(1)))) {
    hypotheses <- rep(list(hypotheses), S)
  }
  
  # check if VCOV and hypotheses are both a non-empty list
  if ( !is.list(VCOV) && length(VCOV) == 0 ) {
    stop("\nrestriktor ERROR: The argument 'VCOV' must be a list of covariance matrices corresponding to the (standardized) parameter estimates of interest.",
         call. = FALSE)
  } 
  # TO DO check of square matrices in list, kan evt door nu fout in goric() te laten gebeuren, maar
  #       dan is het meegeven van study nr wel fijn!
  
  if ( !is.list(hypotheses) && length(hypotheses) == 0 ) {
    stop("\nrestriktor ERROR: hypotheses must be a list.", call. = FALSE)  
  } 
  
  
  # check if the matrices are all symmetrical
  #VCOV_isSym <- vapply(VCOV, isSymmetric, logical(1), check.attributes = FALSE)
  #if (!all(VCOV_isSym)) {
  #  stop(sprintf("\nrestriktor ERROR: the %sth covariance matrix in VCOV is not symmetric.", which(!VCOV_isSym)), call. = FALSE)  
  #}
  
  # number of hypotheses must be equal for each study. In each study a set of 
  # shared theories (i.e., hypotheses) are compared.
  len_H <- vapply(hypotheses, length, integer(1))
  NrHypos <- unique(len_H)
  
  if (length(unique(len_H)) > 1) {
    stop("\nrestriktor ERROR: The number of hypotheses must be identical across all studies.",
         call. = FALSE)
  }
  
  if (length(object) != length(len_H)) {
    stop("\nrestriktor ERROR: The number of hypothesis sets (", length(len_H), 
         ") does not match the number of studies (", S, ").",
         call. = FALSE)
  }
  
  # [CHANGE 2026-10 | Rebecca] comparison = 'complement' when every sub-list holds one hypothesis (fix: only when comparison was not specified)
  # Note: in the case of a list in a list, 
  # and only one hypothesis in the sub-lists, then:
  # comparison <- "complement"
  # BUT only when it was not set to something in the first place!
  # [CHANGE 2026-10 | audit] regression fix: use comparison_missing (missing() on the local copy was always FALSE)
  if (comparison_missing && all(len_H == 1)) {
    comparison <- "complement"
  # [CHANGE 2026-10 | audit] regression fix: default 'unconstrained' otherwise
  } else if (comparison_missing) {
    comparison <- "unconstrained"
  }
  # [/CHANGE 2026-10]
  
  complement_check <- all(len_H == 1)
  if (comparison == "complement") {
    #if ((sameHypo && !comp_check_same) | (!sameHypo && !comp_check_diff)) {
    if (!complement_check) {  
      warning("\nrestriktor WARNING: Only one order-restricted hypothesis is currently supported when comparison = 'complement'. ",
              "The comparison type has been set to 'unconstrained' instead.",
              "Notably, the relative support between informative hypotheses is independent from this choice.",
              call. = FALSE)
      comparison <- "unconstrained"
    }
  }
  
  NrHypos_incl <- NrHypos + 1
  if (comparison == "none") {
    NrHypos_incl <- NrHypos
  }
  # [CHANGE 2026-10 | audit] A6: Heq only for one hypothesis vs its complement; per-study hypothesis sets aligned by name (error on different names, warning when only some are named) instead of by position
  # Heq is only valid when a single order-restricted hypothesis is compared
  # to its complement (as in goric(), with one warning instead of one per study)
  if (Heq && !(comparison == "complement" && NrHypos == 1)) {
    warning("\nrestriktor WARNING: The 'Heq' argument is ignored. ",
            "The 'Heq' option is only valid when a single order-restricted hypothesis ",
            "is compared to its complement (comparison = 'complement').",
            call. = FALSE)
    Heq <- FALSE
    dots$Heq <- NULL
  }
  if (Heq) {
    NrHypos_incl <- NrHypos_incl + 1 # nl, also Heq itself
  }
  
  # Align the per-study hypothesis sets by name: when the hypotheses of all
  # studies are named, each set must contain the same names as the set of
  # study 1; the hypotheses are matched by name, so they may be given in
  # another order (they are then re-ordered to the order of study 1), but
  # different names give an error. When the hypotheses of some (but not all)
  # studies are named, the hypotheses are matched by position (with a warning).
  list_hypo_names <- lapply(hypotheses, names)
  set_named <- vapply(list_hypo_names, function(x) {
    !is.null(x) && all(!is.na(x)) && all(x != "")
  }, logical(1))
  if (all(set_named)) {
    ref_names <- list_hypo_names[[1]]
    if (anyDuplicated(ref_names)) {
      stop("\nrestriktor ERROR: The hypothesis names within a hypothesis set must be unique. ",
           "Found: ", paste(sQuote(ref_names), collapse = ", "), ".", call. = FALSE)
    }
    for (s in seq_len(S)) {
      if (anyDuplicated(list_hypo_names[[s]]) || !setequal(list_hypo_names[[s]], ref_names)) {
        stop("\nrestriktor ERROR: The hypothesis names must be identical across the hypothesis ",
             "sets of all studies (the hypotheses are matched by name). Study 1 uses (",
             paste(ref_names, collapse = ", "), "), while study ", s, " uses (",
             paste(list_hypo_names[[s]], collapse = ", "), "). Please use the same names ",
             "in each set (or leave all hypotheses unnamed, in which case they are matched by position).",
             call. = FALSE)
      }
      if (!identical(list_hypo_names[[s]], ref_names)) {
        message("\nrestriktor Message: The hypotheses of study ", s, " are given in a ",
                "different order than those of study 1. The hypotheses are matched by name.")
        hypotheses[[s]] <- hypotheses[[s]][ref_names]
      }
    }
  } else {
    ref_names <- NULL
    if (any(set_named)) {
      warning("\nrestriktor WARNING: The hypothesis names are missing for some of the studies ",
              "(or some hypotheses within a set are unnamed). The hypotheses are matched ",
              "across studies by position and named 'H1', 'H2', ... (unless 'hypo_names' is given).",
              call. = FALSE)
    }
  }

  if (is.null(hypo_names)) {
    element_hypo_names <- ref_names
  # [/CHANGE 2026-10]
  } else {
    # [CHANGE 2026-10 | Rebecca] hypo_names validated (length and type)
    if (length(hypo_names) != NrHypos) {
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos, " names, \n",
           "namely one for each specified hypothesis. It now consists of ", length(hypo_names), ".",
           call. = FALSE)
    }
    if (!all(is.character(hypo_names))) {
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos, " names. \n",
           "Now, (some of) the elements are not characters.",
           call. = FALSE)
    }
    # [/CHANGE 2026-10]
    # [CHANGE 2026-10 | audit] B13: warning when hypo_names is a permutation of the hypothesis names
    if (!is.null(ref_names) && setequal(hypo_names, ref_names) &&
        !identical(as.character(hypo_names), ref_names)) {
      .evSyn_warn_hypo_names_permuted(hypo_names, ref_names)
    }
    # [/CHANGE 2026-10]
    element_hypo_names <- hypo_names
  }
  
  if (NrHypos == 1 && comparison == "complement") {
    if (!is.null(element_hypo_names)) {
      element_hypo_names <- c(element_hypo_names, "Complement")
      hnames_idx <- element_hypo_names != ""
    } else {
      element_hypo_names <- vector("character", 2L)
      hnames_idx <- element_hypo_names != ""
    }
    hnames <- c("H1", "Complement")
    hnames_idx <- element_hypo_names != ""
    element_hypo_names[!hnames_idx] <- hnames[!hnames_idx]
    hnames <- element_hypo_names
    
    hypotheses <- lapply(hypotheses, function(h) {
      names(h)[1:(length(hnames) - 1L)] <- hnames[-max(length(hnames))]  
      return(h)
    })
    # [CHANGE 2026-10 | audit] Heq: goric() adds the equality-restricted hypothesis first
    if (Heq) {
      # goric() adds the equality-restricted hypothesis (named 'Heq') first
      hnames <- c("Heq", hnames)
    }
    # [/CHANGE 2026-10]
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = 1)
  } else if (comparison == "none") {
    if (!is.null(element_hypo_names)) {
      hnames_idx <- element_hypo_names != ""
    } else {
      element_hypo_names <- vector("character", unique(len_H))
      hnames_idx <- element_hypo_names != ""
    }
    hnames <- c(paste0("H", 1:NrHypos))
    element_hypo_names[!hnames_idx] <- hnames[!hnames_idx]
    hnames <- element_hypo_names
    
    hypotheses <- lapply(hypotheses, function(h) {
      names(h)[seq_len(length(hnames))] <- hnames 
      return(h)
    })
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  } else {
    if (!is.null(element_hypo_names)) {
      element_hypo_names <- c(element_hypo_names, "Unconstrained")
      hnames_idx <- element_hypo_names != ""
    } else {
      element_hypo_names <- vector("character", unique(len_H) + 1L)
      hnames_idx <- element_hypo_names != ""
    }
    hnames <- c(paste0("H", 1:NrHypos), "Unconstrained")
    hnames_idx <- element_hypo_names != ""
    element_hypo_names[!hnames_idx] <- hnames[!hnames_idx]
    hnames <- element_hypo_names
    
    hypotheses <- lapply(hypotheses, function(h) {
      names(h)[1:(length(hnames) - 1L)] <- hnames[-max(length(hnames))]
      return(h)
    })
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  }
  
  # [CHANGE 2026-10 | audit] B2/B21: priorICweights validated and matched by name (hypothesis names accepted as alias)
  # Check the prior IC weights (one for each hypothesis, incl. the failsafe
  # hypothesis; matched by name when named; rescaled to sum to 1). When the
  # hypotheses are named, these names (with 'Heq' and the failsafe hypothesis
  # at their positions) are accepted as well, also when 'hypo_names' is given.
  input_names <- NULL
  if (!is.null(ref_names)) {
    input_names <- hnames
    input_names[Heq + seq_len(NrHypos)] <- ref_names
  }
  priorICweights <- .evSyn_check_priorICweights(priorICweights, NrHypos_incl, hnames,
                                                input_names = input_names)
  # [/CHANGE 2026-10]
  
  LL_m <- LL_weights_m <- GORICA_m <- GORICA_weight_m <- PT <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  colnames(LL_m) <- colnames(LL_weights_m) <- colnames(GORICA_m) <- colnames(GORICA_weight_m) <- colnames(PT) <- hnames
  # rownames are set after determining the order of the studies
  #
  study_sample_nobs <- unlist(study_sample_nobs) # when it comes from escalc, then it is a list
  for (s in 1:S) {
    # [CHANGE 2026-10 | audit] goric() called via do.call with dots (without priorWeights); number of returned models checked (Heq)
    # Note: the remaining arguments (dots, i.e., '...' without the deprecated
    # 'priorWeights') are passed on to goric().
    res_goric <- do.call(goric, c(list(object[[s]], VCOV = VCOV[[s]],
                                       hypotheses = hypotheses[[s]],
                                       type = type, comparison = comparison,
                                       sample_nobs = study_sample_nobs[s]),
                                  dots))
    
    # the number of models in the result of goric() must match the number
    # of columns (e.g., goric() ignores 'Heq' when the hypothesis contains
    # no inequality restrictions)
    if (nrow(res_goric$result) != NrHypos_incl) {
      stop("\nrestriktor ERROR: For study ", s, ", goric() returned ", nrow(res_goric$result),
           " models (", paste(res_goric$result$model, collapse = ", "), "), while ",
           NrHypos_incl, " were expected (", paste(hnames, collapse = ", "), ").",
           if (Heq) " The 'Heq' option cannot be used for this set of hypotheses.",
           call. = FALSE)
    }
    # [/CHANGE 2026-10]
    
    if (comparison == "unconstrained") {
      ratio.weight_mu[s, ] <- res_goric$ratio.gw[, NrHypos_incl]
    } else if (comparison == "complement") {
      # [CHANGE 2026-10 | audit] Heq: the order-restricted hypothesis is row 1 + Heq
      # the (single) order-restricted hypothesis versus its complement
      ratio.weight_mu[s, ] <- res_goric$ratio.gw[1 + Heq, NrHypos_incl]
    } 
    
    LL_m[s, ] <- res_goric$result$loglik
    LL_weights_m[s, ] <- res_goric$result$loglik.weights
    GORICA_m[s, ] <- res_goric$result[[type]] #res_goric$result$gorica
    # [CHANGE 2026-10 | audit] B17: open question (prior in study-specific weights)
    # TO DO: open question (Leonard/Rebecca): should study-specific weights include priorICweights in all routes? Currently est/gorica route does not, ICvalues/ICweights routes do.
    GORICA_weight_m[s, ] <- res_goric$result[[paste0(type, ".weights")]]
    PT[s, ] <- res_goric$result$penalty
  }
  # [CHANGE 2026-10 | audit] penalty_factor taken from the goric() results
  # penalty factor used in goric() (i.e., IC = -2 * LL + penalty_factor * PT)
  penalty_factor <- res_goric$penalty_factor
  if (is.null(penalty_factor)) {
    penalty_factor <- 2
  }
  # [/CHANGE 2026-10]
  
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  if (is.numeric(order_studies)) {
    # User-specified numeric order vector
    # [CHANGE 2026-10 | audit] order_studies already validated; drop = FALSE (single study)
    orderStudies <- order_studies
    LL_m <- LL_m[orderStudies, , drop = FALSE]
    LL_weights_m <- LL_weights_m[orderStudies, , drop = FALSE]
    GORICA_m <- GORICA_m[orderStudies, , drop = FALSE]
    GORICA_weight_m <- GORICA_weight_m[orderStudies, , drop = FALSE]
    PT <- PT[orderStudies, , drop = FALSE]
    # [/CHANGE 2026-10]
  } else if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    # [CHANGE 2026-10 | audit] preferred hypothesis from the final prior-weighted, study-weighted IC weights (same normalisation as the final synthesis)
    # That is, the hypothesis with the highest final (prior-weighted) IC weight.
    OverallGoric <- .evSyn_cum_IC(LL_m, PT, study_weights_S, type_ev, penalty_factor)[S, ]
    OverallPrefHypo <- .evSyn_pref_hypo(OverallGoric, priorICweights)
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    orderStudies <- order(GORICA_weight_m[, OverallPrefHypo], decreasing = decreasing)
    #
    # [CHANGE 2026-10 | audit] drop = FALSE (single study)
    LL_m <- LL_m[orderStudies, , drop = FALSE]
    LL_weights_m <- LL_weights_m[orderStudies, , drop = FALSE]
    GORICA_m <- GORICA_m[orderStudies, , drop = FALSE]
    GORICA_weight_m <- GORICA_weight_m[orderStudies, , drop = FALSE]
    PT <- PT[orderStudies, , drop = FALSE]
    # [/CHANGE 2026-10]
  }
  
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    if (length(study_names) != S) {
      stop("restriktor ERROR: Length of study_names must match the number of studies.", call. = FALSE)
    }
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(LL_m) <- rownames(LL_weights_m) <- rownames(GORICA_m) <- rownames(GORICA_weight_m) <- rownames(PT) <- study_names
  # [CHANGE 2026-10 | audit] ratio.weight_mu re-ordered with the studies
  ratio.weight_mu <- ratio.weight_mu[orderStudies, , drop = FALSE]
  rownames(ratio.weight_mu) <- study_names
  # [CHANGE 2026-10 | audit] study weights and study_sample_nobs re-ordered with the studies
  # Re-order the study-specific settings as well, such that they match the 
  # (re-ordered) studies and results do not depend on the input order.
  study_weights <- study_weights[orderStudies]
  study_weights_S <- study_weights_S[orderStudies]
  if (!is.null(study_sample_nobs)) {
    study_sample_nobs <- study_sample_nobs[orderStudies]
  }
  # [/CHANGE 2026-10]
  
  CumulativeLLWeights <- CumulativeGoricaWeights <- CumulativeGorica <- matrix(NA, nrow = S+1, ncol = NrHypos_incl)
  colnames(CumulativeLLWeights) <- colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- hnames
  sequence <- paste0("Study nr.s 1-", 1:S, "   ")
  sequence[1] <- "Study nr.  1   "
  rownames(CumulativeLLWeights) <- rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
  #
  # [CHANGE 2026-10 | audit] B16: cumulative IC values/weights via .evSyn_cum_IC and log-sum-exp (study weights, penalty_factor, priorICweights); cumulative LL study-weighted and averaged for 'average' (replaces the per-type loops)
  # Cumulative IC values (possibly weighted using study weights):
  # - added:   sum of LL values and sum of PT values,
  # - equal:   sum of LL values and average of PT values,
  # - average: average of LL values and average of PT values.
  # Cumulative IC weights: based on the cumulative IC values, taking into 
  # account possible prior hypothesis weights (computed on the log scale).
  CumulativeGorica[1:S, ] <- .evSyn_cum_IC(LL_m, PT, study_weights_S, type_ev, penalty_factor)
  CumulativeGoricaWeights[1:S, ] <- .evSyn_IC_weights_rows(CumulativeGorica[1:S, , drop = FALSE], priorICweights)
  
  # cumulative log-likelihood values (possibly weighted using study weights,
  # in the same way as for the cumulative IC values: summed, or averaged in
  # the average-evidence approach); the log-likelihood weights are based on
  # these, so that they are consistent with the cumulative IC values.
  # Note: priorICweights are not used for the log-likelihood weights, since 
  #       these are not IC weights (and there is no penalty term).
  Cumulative_LL <- .evSyn_cum_weighted(LL_m, study_weights_S)
  if (type_ev == "average") {
    # average-evidence approach: average of the log-likelihood values (over
    # the positively weighted studies so far), like the IC and penalty values
    Cumulative_LL <- Cumulative_LL / pmax(.evSyn_n_pos(study_weights_S), 1)
  }
  Cumulative_LL <- matrix(Cumulative_LL, nrow = nrow(LL_m), 
                          dimnames = list(sequence, colnames(LL_m)))
  # [/CHANGE 2026-10]
  
  # cumulative log_likelihood weights
  for (l in 1:S) {
    # [CHANGE 2026-10 | audit] cumulative LL weights from the (weighted) Cumulative_LL
    CumulativeLL <- -2 * Cumulative_LL[l, ]
    minLL <- min(CumulativeLL)
    CumulativeLLWeights[l, ] <- exp(-0.5*(CumulativeLL-minLL)) / sum(exp(-0.5*(CumulativeLL-minLL)))
  }
  # final cumulative log-likelihood value
  Cumulative_LL_final <- -2*Cumulative_LL[S, , drop = FALSE]
  rownames(Cumulative_LL_final) <- "Final"
  minLL <- min(Cumulative_LL_final)
  # cumulative log-likelihood weights
  Final.LL.weights <- exp(-0.5*(Cumulative_LL_final-minLL)) / sum(exp(-0.5*(Cumulative_LL_final-minLL)))
  Final.LL.weights <- Final.LL.weights[,, drop = TRUE]
  Final.ratio.LL.weights <- Final.LL.weights %*% t(1/Final.LL.weights)
  rownames(Final.ratio.LL.weights) <- hnames
  
  
  # add final row
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  CumulativeLLWeights[(S+1), ] <- CumulativeLLWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  # [CHANGE 2026-10 | Rebecca] diagonal of the ratio matrix set to 1 (zero weights)
  diag(Final.ratio.GORICA.weights) <- 1 # If a weight is zero, then you get Inf and NaN; this way you get Inf and 1.
  
  rownames(Final.ratio.GORICA.weights) <- hnames
  
  
  # Output
  if (NrHypos == 1 && comparison == "complement") {
    # [CHANGE 2026-10 | audit] Heq: name of the order-restricted hypothesis
    colnames(ratio.weight_mu) <- c(paste0(hnames[1 + Heq], " vs. ", "Complement"))
    colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
  } else if (comparison == "none") {
    #colnames(ratio.weight_mu) <- c(paste0(hnames[1], " vs. ", "Complement"))
    colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
  } else { 
    # unconstrained
    colnames(ratio.weight_mu) <- c(paste0(colnames(CumulativeGorica), " vs. ", "Unconstrained"))
    colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. ", colnames(CumulativeGorica)))
  }
  
  out <- list(type = type,
              type_ev = type_ev,
              hypotheses = hypotheses,
              # [CHANGE 2026-10 | Rebecca] priorICweights in the output
              priorICweights = priorICweights,
              n_studies = S,
              order_studies = orderStudies,
              study_names = study_names,
              # [CHANGE 2026-10 | Rebecca] study_weights in the output
              study_weights = study_weights,
              study_sample_nobs = study_sample_nobs,
              # [CHANGE 2026-10 | audit] penalty_factor in the output
              penalty_factor = penalty_factor,
              PT_m = PT,
              GORICA_weight_m = GORICA_weight_m, 
              LL_weights_m = LL_weights_m,
              GORICA_m = GORICA_m, 
              LL_m = LL_m, 
              Cumulative_GORICA = CumulativeGorica, 
              Cumulative_LL = Cumulative_LL,
              Cumulative_GORICA_weights = CumulativeGoricaWeights,
              Cumulative_LL_weights = CumulativeLLWeights,                  
              ratio_GORICA_weight_mu = ratio.weight_mu, 
              Final_ratio_GORICA_weights = Final.ratio.GORICA.weights,
              Final_ratio_LL_weights = Final.ratio.LL.weights
  )
  # TO DO welke volgorde en dan in alle functies zo gelijk mogelijk maken ook
  
  class(out) <- c("evSyn_est", "evSyn")
  
  return(out)
}


# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on log likelihood and penalty values
evSyn_LL <- function(object, ..., PT = list(), 
                     type_ev = c("added", "equal", "average"),
                     # [CHANGE 2026-10 | Rebecca] priorICweights argument
                     hypo_names = c(), priorICweights = NULL,
                     type = c("goric", "goricc", "gorica", "goricac"),
                     order_studies = c("input_order", "ascending", "descending"),
                     study_names = c(),
                     # [CHANGE 2026-10 | audit] study_weights and penalty_factor arguments
                     study_weights = NULL,
                     penalty_factor = 2) {
  
  if (missing(type_ev)) 
    type_ev <- "added"
  type_ev <- match.arg(type_ev)
  
  if (missing(type)) 
    type <- "gorica"
  type <- match.arg(type)
  
  # [CHANGE 2026-10 | audit] B15 priorWeights compat; A7/B14 input checks; penalty_factor validated; B13 names aligned across studies and PT aligned to LL by name
  # Backwards compatibility: 'priorWeights' is renamed to 'priorICweights'.
  priorICweights <- .evSyn_priorWeights_compat(priorICweights, list(...))

  # check the log-likelihood values (object) and the penalty values (PT):
  # lists of numeric vectors without NA/NaN, with one value for each
  # hypothesis, the same hypotheses in each study
  if (length(PT) == 0) {
    stop("\nrestriktor ERROR: PT must be a list of penalty values (one vector for each study).",
         call. = FALSE)
  }
  .evSyn_check_list_input(object, name = "object")
  .evSyn_check_list_input(PT, name = "PT", ref = object, ref_name = "object")
  
  # penalty factor: IC = -2 * LL + penalty_factor * PT
  if (!is.numeric(penalty_factor) || length(penalty_factor) != 1L || 
      is.na(penalty_factor) || penalty_factor < 0) {
    stop("\nrestriktor ERROR: The argument 'penalty_factor' must be a single non-negative number.",
         call. = FALSE)
  }
  
  # If the input vectors carry hypothesis names, these must denote the same 
  # hypotheses across studies (matched by name); otherwise matched by position.
  object <- .evSyn_check_input_names(object, hypo_names)
  # the hypothesis names of the input (if any; accepted for a named 'priorICweights')
  input_names <- .evSyn_input_hypo_names(object)
  # The penalty values must match the (possibly re-ordered) log-likelihood 
  # values per study: matched by name when named, otherwise by position.
  PT <- .evSyn_align_second_input(object, PT, name = "PT", object_name = "object")
  # [/CHANGE 2026-10]
  
  LL_m <- object
  S <- length(LL_m)
  NrHypos <- length(LL_m[[1]]) - 1
  # [CHANGE 2026-10 | Rebecca] NrHypos_incl
  NrHypos_incl <- NrHypos + 1
  
  # [CHANGE 2026-10 | audit] B22: order_studies validated via .evSyn_order_studies
  # Check the order of the studies (character string, permutation of 1:S, or
  # permutation of the study names)
  order_studies <- .evSyn_order_studies(order_studies, S, study_names)
  
  if (is.null(hypo_names)) {
    # [CHANGE 2026-10 | Rebecca] hnames via NrHypos_incl
    hnames <- paste0("H", 1:NrHypos_incl)
  } else {
    # [CHANGE 2026-10 | Rebecca] hypo_names validated (length and type)
    if (length(hypo_names) != NrHypos_incl) {
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos_incl, " names, \n",
           "namely one for each specified hypothesis. It now consists of ", length(hypo_names), ".",
           call. = FALSE)
    }
    if (!all(is.character(hypo_names))) {
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos_incl, " names. \n",
           "Now, (some of) the elements are not characters.",
           call. = FALSE)
    }
    # [/CHANGE 2026-10]
    hnames <- hypo_names
  }
  
  # [CHANGE 2026-10 | audit] B2/B21: priorICweights validated and matched by name
  # Check the prior IC weights (one for each hypothesis; matched by name when
  # named; rescaled to sum to 1)
  priorICweights <- .evSyn_check_priorICweights(priorICweights, NrHypos_incl, hnames,
                                                input_names = input_names)
  # [/CHANGE 2026-10]
  
  # [CHANGE 2026-10 | audit] B21: study_weights validated
  # Check the study weights (zero weights are allowed; see .evSyn_cum_weighted)
  study_weights <- .evSyn_check_study_weights(study_weights, S, study_names)
  study_weights_S <- study_weights$study_weights_S # Now, they sum up to S
  study_weights <- study_weights$study_weights
  # [/CHANGE 2026-10]
  
  # [CHANGE 2026-10 | audit] unname before rbind; penalty_factor in the IC values
  LL_m <- do.call(rbind, lapply(LL_m, unname))
  PT <- do.call(rbind, lapply(PT, unname))
  IC <- -2 * LL_m + penalty_factor * PT
  #
  # [CHANGE 2026-10 | audit] B14/B17: study-specific IC weights via log-sum-exp (prior-free; open question noted)
  # TO DO: open question (Leonard/Rebecca): should study-specific weights include priorICweights in all routes? Currently est/gorica route does not, ICvalues/ICweights routes do.
  GORICA_weight_m <- .evSyn_IC_weights_rows(IC)
  
  
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  if (is.numeric(order_studies)) {
    # User-specified numeric order vector
    # [CHANGE 2026-10 | audit] order_studies already validated; drop = FALSE (single study)
    orderStudies <- order_studies
    LL_m <- LL_m[orderStudies, , drop = FALSE]
    PT <- PT[orderStudies, , drop = FALSE]
    IC <- IC[orderStudies, , drop = FALSE]
    GORICA_weight_m <- GORICA_weight_m[orderStudies, , drop = FALSE]
    # [/CHANGE 2026-10]
  } else if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    # [CHANGE 2026-10 | audit] preferred hypothesis from the final prior-weighted, study-weighted IC weights
    # That is, the hypothesis with the highest final (prior-weighted) IC weight.
    OverallGoric <- .evSyn_cum_IC(LL_m, PT, study_weights_S, type_ev, penalty_factor)[S, ]
    OverallPrefHypo <- .evSyn_pref_hypo(OverallGoric, priorICweights)
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    orderStudies <- order(GORICA_weight_m[, OverallPrefHypo], decreasing = decreasing)
    #
    # [CHANGE 2026-10 | audit] drop = FALSE (single study)
    LL_m <- LL_m[orderStudies, , drop = FALSE]
    PT <- PT[orderStudies, , drop = FALSE]
    IC <- IC[orderStudies, , drop = FALSE]
    GORICA_weight_m <- GORICA_weight_m[orderStudies, , drop = FALSE]
    # [/CHANGE 2026-10]
  }
  
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    # [CHANGE 2026-10 | Rebecca] study_names length checked
    if (length(study_names) != S) {
      stop("\nrestriktor ERROR: The argument 'study_names' should consist of ", S, " names, \n",
           "namely one for each study. It now consists of ", length(study_names), ".",
           call. = FALSE)
    }
    # [/CHANGE 2026-10]
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(LL_m) <- rownames(PT) <- rownames(IC) <- rownames(GORICA_weight_m) <- study_names
  # Set colnames
  colnames(LL_m) <- colnames(PT) <- colnames(IC) <- colnames(GORICA_weight_m) <- hnames
  # [CHANGE 2026-10 | audit] study weights re-ordered with the studies
  # Re-order the study weights as well, such that they match the (re-ordered) 
  # studies and results do not depend on the input order.
  study_weights <- study_weights[orderStudies]
  study_weights_S <- study_weights_S[orderStudies]
  # [/CHANGE 2026-10]
  
  sequence <- paste0("Study nr.s 1-", 1:S, "   ")
  sequence[1] <- "Study nr.  1   "
  #
  # [CHANGE 2026-10 | audit] B16: cumulative LL study-weighted and averaged for 'average'
  # cumulative log-likelihood values (possibly weighted using study weights,
  # in the same way as for the cumulative IC values)
  # Note: priorICweights are not used for the log-likelihood weights, since 
  #       these are not IC weights (and there is no penalty term).
  Cumulative_LL <- .evSyn_cum_weighted(LL_m, study_weights_S)
  if (type_ev == "average") {
    # average-evidence approach: average of the log-likelihood values (over
    # the positively weighted studies so far), like the IC and penalty values
    Cumulative_LL <- Cumulative_LL / pmax(.evSyn_n_pos(study_weights_S), 1)
  }
  # [/CHANGE 2026-10]
  Cumulative_LL <- matrix(Cumulative_LL, nrow = nrow(LL_m), 
                          dimnames = list(sequence, colnames(LL_m)))
  # final cumulative log-likelihood value
  Cumulative_LL_final <- -2*Cumulative_LL[S, , drop = FALSE]
  minLL <- min(Cumulative_LL_final)
  # cumulative log-likelihood weights
  Final.LL.weights <- exp(-0.5*(Cumulative_LL_final-minLL)) / sum(exp(-0.5*(Cumulative_LL_final-minLL)))
  Final.LL.weights <- Final.LL.weights[,, drop = TRUE]
  Final.ratio.LL.weights <- Final.LL.weights %*% t(1/Final.LL.weights) 
  
  sumLL <- 0
  LL_weights_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  CumulativeLLWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  rownames(LL_weights_m) <- study_names
  rownames(CumulativeLLWeights) <- c(sequence, "Final")
  colnames(LL_weights_m) <- colnames(CumulativeLLWeights) <- hnames
  for (l in 1:S) {
    LL <- -2*LL_m[l, ]
    delta_LL <- LL - min(LL)
    LL_weights_m[l, ] <- exp(-0.5 * delta_LL) / sum(exp(-0.5 * delta_LL))
    #
    # [CHANGE 2026-10 | audit] cumulative LL weights from the (weighted) Cumulative_LL
    CumulativeLL <- -2 * Cumulative_LL[l, ]
    minLL <- min(CumulativeLL)
    CumulativeLLWeights[l, ] <- exp(-0.5*(CumulativeLL-minLL)) / sum(exp(-0.5*(CumulativeLL-minLL)))
  }
  
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- hnames
  # [CHANGE 2026-10 | audit] B16: cumulative IC values/weights via .evSyn_cum_IC and log-sum-exp (replaces the per-type loops)
  # Cumulative IC values (possibly weighted using study weights; see .evSyn_cum_IC) 
  # and cumulative IC weights (taking into account possible prior hypothesis 
  # weights; computed on the log scale).
  CumulativeGorica[1:S, ] <- .evSyn_cum_IC(LL_m, PT, study_weights_S, type_ev, penalty_factor)
  CumulativeGoricaWeights[1:S, ] <- .evSyn_IC_weights_rows(CumulativeGorica[1:S, , drop = FALSE], priorICweights)
  # [/CHANGE 2026-10]
  
  # fill in the final row  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  CumulativeLLWeights[(S+1), ] <- CumulativeLLWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  # [CHANGE 2026-10 | Rebecca] diagonal of the ratio matrix set to 1 (zero weights)
  diag(Final.ratio.GORICA.weights) <- 1 # If a weight is zero, then you get Inf and NaN; this way you get Inf and 1.
  
  rownames(Final.ratio.LL.weights) <- rownames(Final.ratio.GORICA.weights) <- hnames
  colnames(Final.ratio.LL.weights) <- colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hnames)
  
  out <- list(type = type,
    type_ev = type_ev,
    #hypotheses = hypo_names,
    # [CHANGE 2026-10 | Rebecca] priorICweights in the output
    priorICweights = priorICweights,
    n_studies = S,
    order_studies = orderStudies,
    study_names = study_names,
    # [CHANGE 2026-10 | Rebecca] study_weights in the output
    study_weights = study_weights,
    #study_sample_nobs = study_sample_nobs,
    # [CHANGE 2026-10 | audit] penalty_factor in the output
    penalty_factor = penalty_factor,
    PT_m = PT, 
    GORICA_weight_m = GORICA_weight_m,
    LL_weights_m = LL_weights_m,
    GORICA_m = IC, 
    LL_m = LL_m, 
    Cumulative_GORICA_weights = CumulativeGoricaWeights,
    Cumulative_LL_weights = CumulativeLLWeights,
    Cumulative_GORICA = CumulativeGorica, 
    Cumulative_LL = Cumulative_LL,
    Final_ratio_GORICA_weights = Final.ratio.GORICA.weights,
    Final_ratio_LL_weights = Final.ratio.LL.weights
  )
  
  class(out) <- c("evSyn_LL", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA values
evSyn_ICvalues <- function(object, ..., type_ev = c("added", "average"), 
                           # [CHANGE 2026-10 | Rebecca] priorICweights argument
                           hypo_names = c(), priorICweights = NULL,
                           type = c("goric", "goricc", "gorica", "goricac"),
                           order_studies = c("input_order", "ascending", "descending"),
                           # [CHANGE 2026-10 | Rebecca] study_weights argument
                           study_names = c(),
                           study_weights = NULL) {
  
  if (missing(type_ev)) 
    type_ev <- "added"
  # [CHANGE 2026-10 | audit] B23: 'equal' falls back to 'added' with a message
  type_ev <- .evSyn_type_ev_IC(type_ev, "IC values")
  
  if (missing(type)) 
    type <- "gorica"
  type <- match.arg(type)

  # [CHANGE 2026-10 | audit] B15 priorWeights compat; A7/B14 input checks; B13 names aligned across studies
  # Backwards compatibility: 'priorWeights' is renamed to 'priorICweights'.
  priorICweights <- .evSyn_priorWeights_compat(priorICweights, list(...))

  # check the input: a list of numeric vectors without NA/NaN, with one value
  # for each hypothesis, the same hypotheses in each study
  .evSyn_check_list_input(object, name = "object")
  
  # If the input vectors carry hypothesis names, these must denote the same 
  # hypotheses across studies (matched by name); otherwise matched by position.
  object <- .evSyn_check_input_names(object, hypo_names)
  # the hypothesis names of the input (if any; accepted for a named 'priorICweights')
  input_names <- .evSyn_input_hypo_names(object)
  # [/CHANGE 2026-10]
  
  IC <- object
  S  <- length(IC)
  NrHypos <- length(IC[[1]]) - 1
  # [CHANGE 2026-10 | Rebecca] TO DO notes; NrHypos_incl
  # TO DO waarom -1 (op meerdere plekken), was vast ergens voor nodig.... 
  # TO DO wat ik kan bedenken maar we nu nog niets mee doen:
  # We assume that the last one is the failsafe Hunc
  # If not, users can specify the hypotheses names as well (in 'hypo_names').
  NrHypos_incl <- NrHypos + 1
  # [/CHANGE 2026-10]
  GORICA_weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  
  if (is.null(hypo_names)) {
    # [CHANGE 2026-10 | Rebecca] hnames via NrHypos_incl; hypo_names validated (length and type)
    #hnames <- paste0("H", 1:NrHypos)
    #hnames <- c(hnames, "unconstrained")
    hnames <- paste0("H", 1:NrHypos_incl)
  } else {
    if (length(hypo_names) != NrHypos_incl) {
      # [CHANGE 2026-10 | audit] NrHypos_incl in the message
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos_incl, " names, \n",
           "namely one for each specified hypothesis. It now consists of ", length(hypo_names), ".",
           call. = FALSE)
    }
    if (!all(is.character(hypo_names))) {
      # [CHANGE 2026-10 | audit] NrHypos_incl in the message
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos_incl, " names. \n",
           "Now, (some of) the elements are not characters.",
           call. = FALSE)
    }
      # [/CHANGE 2026-10]
    hnames <- hypo_names
  }
  
  # [CHANGE 2026-10 | audit] B2/B21: priorICweights validated and matched by name
  # Check the prior IC weights (one for each hypothesis; matched by name when
  # named; rescaled to sum to 1)
  priorICweights <- .evSyn_check_priorICweights(priorICweights, NrHypos_incl, hnames,
                                                input_names = input_names)
  # [/CHANGE 2026-10]
  
  # [CHANGE 2026-10 | audit] B21: study_weights validated
  # Check the study weights (zero weights are allowed; see .evSyn_cum_weighted)
  study_weights <- .evSyn_check_study_weights(study_weights, S, study_names)
  study_weights_S <- study_weights$study_weights_S # Now, they sum up to S
  study_weights <- study_weights$study_weights
  # [/CHANGE 2026-10]
  
  # [CHANGE 2026-10 | audit] B22: order_studies validated via .evSyn_order_studies (replaces the match.arg block)
  # Check the order of the studies (character string, permutation of 1:S, or
  # permutation of the study names)
  order_studies <- .evSyn_order_studies(order_studies, S, study_names)
  
  # [CHANGE 2026-10 | audit] unname before rbind
  IC <- do.call(rbind, lapply(IC, unname))
  #
  # [CHANGE 2026-10 | audit] B14: study-specific IC weights via log-sum-exp (replaces the loop)
  GORICA_weight_m <- .evSyn_IC_weights_rows(IC)
  
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  if (is.numeric(order_studies)) {
    # User-specified numeric order vector
    # [CHANGE 2026-10 | audit] order_studies already validated; drop = FALSE (single study)
    orderStudies <- order_studies
    IC <- IC[orderStudies, , drop = FALSE]
    GORICA_weight_m <- GORICA_weight_m[orderStudies, , drop = FALSE]
  } else if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    # [CHANGE 2026-10 | audit] preferred hypothesis from the final study-weighted IC values
    # That is, the hypothesis with the highest final (prior-weighted) IC weight.
    OverallGoric <- .evSyn_cum_weighted(IC, study_weights_S)[S, ]
    if (type_ev == "average") { 
      # average-evidence approach
      # [CHANGE 2026-10 | audit] average over the positively weighted studies
      OverallGoric <- OverallGoric / sum(study_weights_S > 0)
    } else {
      # type_ev == "added" (or when "equal", because then it is overruled to be "added")
      type_ev = "added"
    }
    # [CHANGE 2026-10 | audit] preferred hypothesis = highest prior-weighted IC weight
    OverallPrefHypo <- .evSyn_pref_hypo(OverallGoric, priorICweights)
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    orderStudies <- order(GORICA_weight_m[, OverallPrefHypo], decreasing = decreasing)
    #
    # [CHANGE 2026-10 | audit] drop = FALSE (single study)
    IC <- IC[orderStudies, , drop = FALSE]
    GORICA_weight_m <- GORICA_weight_m[orderStudies, , drop = FALSE]
  }
  
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  #
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(IC) <- colnames(GORICA_weight_m) <- hnames
  #
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    # [CHANGE 2026-10 | Rebecca] study_names length checked
    if (length(study_names) != S) {
      stop("\nrestriktor ERROR: The argument 'study_names' should consist of ", S, " names, \n",
           "namely one for each study. It now consists of ", length(study_names), ".",
           call. = FALSE)
    }
    # [/CHANGE 2026-10]
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(IC) <- rownames(GORICA_weight_m) <- study_names
  # [CHANGE 2026-10 | audit] study weights re-ordered with the studies
  # Re-order the study weights as well, such that they match the (re-ordered) 
  # studies and results do not depend on the input order.
  study_weights <- study_weights[orderStudies]
  study_weights_S <- study_weights_S[orderStudies]
  # [/CHANGE 2026-10]
  sequence <- paste0("Study nr.s 1-", 1:S, "   ")
  sequence[1] <- "Study nr.  1   "
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(sequence, "Final")
  #
  # [CHANGE 2026-10 | audit] cumulative IC values via .evSyn_cum_weighted (study weights, zero weights) and IC weights via log-sum-exp (replaces the per-type loops)
  # Cumulative IC values (possibly weighted using study weights):
  # - added:   sum of IC values (also when "equal", because then it is overruled to be "added"),
  # - average: average of IC values.
  # Cumulative IC weights: based on the cumulative IC values, taking into 
  # account possible prior hypothesis weights (computed on the log scale).
  CumulativeGorica[1:S, ] <- .evSyn_cum_weighted(IC, study_weights_S)
  if (type_ev == "average") { 
    CumulativeGorica[1:S, ] <- CumulativeGorica[1:S, , drop = FALSE] / pmax(.evSyn_n_pos(study_weights_S), 1)
  }
  CumulativeGoricaWeights[1:S, ] <- .evSyn_IC_weights_rows(CumulativeGorica[1:S, , drop = FALSE], priorICweights)
  # [/CHANGE 2026-10]
  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  # [CHANGE 2026-10 | Rebecca] diagonal of the ratio matrix set to 1 (zero weights)
  diag(Final.ratio.GORICA.weights) <- 1 # If a weight is zero, then you get Inf and NaN; this way you get Inf and 1.
  
  rownames(Final.ratio.GORICA.weights) <- hnames
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hnames)
  
  # [CHANGE 2026-10 | Rebecca] priorICweights applied to the study-specific IC weights
  # Use priorICweights (i.e., a priori likeliness for each hypotheses)
  # [CHANGE 2026-10 | audit] sweep() over the columns (no recycling across the wrong dimension); B17 open question noted
  # Note: sweep() multiplies each column (hypothesis) with its own prior weight.
  # TO DO: open question (Leonard/Rebecca): should study-specific weights include priorICweights in all routes? Currently est/gorica route does not, ICvalues/ICweights routes do.
  GORICA_weight_m <- sweep(GORICA_weight_m, 2, priorICweights, "*")
  GORICA_weight_m <- GORICA_weight_m / rowSums(GORICA_weight_m)
  # [/CHANGE 2026-10]
  
  out <- list(type             = type,
    type_ev           = type_ev,
    #hypotheses       = hypo_names,
    # [CHANGE 2026-10 | Rebecca] priorICweights in the output
    priorICweights = priorICweights,
    n_studies         = S,
    order_studies     = orderStudies,
    study_names       = study_names,
    # [CHANGE 2026-10 | Rebecca] study_weights in the output
    study_weights = study_weights,
    #study_sample_nobs = study_sample_nobs,
    GORICA_m          = IC, 
    GORICA_weight_m   = GORICA_weight_m,
    Cumulative_GORICA = CumulativeGorica, 
    Cumulative_GORICA_weights  = CumulativeGoricaWeights,
    Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  
  # if (!is.null(messageAdded)) {
  #   out <- append(out, messageAdded)
  #   # TO DO dit is ws nodig als message meegeven, doe ook bij andere varianten dan nog
  # }
  
  class(out) <- c("evSyn_ICvalues", "evSyn")
  
  return(out)
  
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA weights or 
# (Bayesian) posterior model probabilities
evSyn_ICweights <- function(object, ..., type_ev = c("added", "average"), 
                            # [CHANGE 2026-10 | Rebecca] priorICweights argument (replaces priorWeights)
                            hypo_names = c(), priorICweights = NULL, 
                            type = c("goric", "goricc", "gorica", "goricac"),
                            order_studies = c("input_order", "ascending", "descending"),
                            # [CHANGE 2026-10 | Rebecca] study_weights argument
                            study_names = c(),
                            study_weights = NULL) {
  
  if (missing(type_ev)) 
    type_ev <- "added"
  # [CHANGE 2026-10 | audit] B23: 'equal' falls back to 'added' with a message
  type_ev <- .evSyn_type_ev_IC(type_ev, "IC weights")
  
  if (missing(type)) 
    type <- "gorica"
  type <- match.arg(type)
  
  # [CHANGE 2026-10 | audit] B15 priorWeights compat; A7/B14 input checks; B13 names aligned across studies
  # Backwards compatibility: 'priorWeights' is renamed to 'priorICweights'.
  priorICweights <- .evSyn_priorWeights_compat(priorICweights, list(...))

  # check the input: a list of numeric vectors without NA/NaN, with one value
  # for each hypothesis, the same hypotheses in each study
  .evSyn_check_list_input(object, name = "object")
  
  # If the input vectors carry hypothesis names, these must denote the same 
  # hypotheses across studies (matched by name); otherwise matched by position.
  object <- .evSyn_check_input_names(object, hypo_names)
  # the hypothesis names of the input (if any; accepted for a named 'priorICweights')
  input_names <- .evSyn_input_hypo_names(object)
  # [/CHANGE 2026-10]
  
  Weights <- object
  # [CHANGE 2026-10 | Rebecca] input checked to be IC weights (between 0 and 1, summing to 1)
  # Check whether weights between 0 and 1 (and sum to 1)
  min0 <- all(abs(vapply(object, min, numeric(1)) >= 0)) 
  max1 <- all(abs(vapply(object, max, numeric(1)) <= 1)) 
  sum1 <- all(abs(vapply(object, sum, numeric(1)) - 1) <= sqrt(.Machine$double.eps))
  if (min0 + max1 + sum1 != 3) {
    text <- paste0("\nrestriktor ERROR: Please check the input. The function expects IC weights, but:. \n",
                   "ICweights are values between 0 and 1 and (per study) sum to 1. \n")
    # [CHANGE 2026-10 | audit] fix: negated condition (message was given when the check passed)
    if (!min0) {
      text <- paste0(text,
                     "Not all values are >= 0. \n")
    }
    # [CHANGE 2026-10 | audit] fix: negated condition
    if (!max1) {
      text <- paste0(text,
                     "Not all values are <= 1. \n")
    }
    # [CHANGE 2026-10 | audit] fix: negated condition (was max1)
    if (!sum1) {
      text <- paste0(text,
                     "For one or more studies, the values do not sum to 1. \n")
    }
    stop(text, call. = FALSE)
  }
    # [/CHANGE 2026-10]
  
  S <- length(Weights)
  # [CHANGE 2026-10 | audit] unname before rbind
  Weights <- do.call(rbind, lapply(Weights, unname)) 
  NrHypos <- ncol(Weights)
  # [CHANGE 2026-10 | Rebecca] NrHypos_incl (prior-weight code replaced by the central check below)
  NrHypos_incl <- NrHypos
  
  if (is.null(hypo_names)) {
    # [CHANGE 2026-10 | Rebecca] hypo_names validated (length and type)
    hypo_names <- paste0("H", 1:NrHypos) 
  } else {
    if (length(hypo_names) != NrHypos) {
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos, " names, \n",
           "namely one for each specified hypothesis. It now consists of ", length(hypo_names), ".",
           call. = FALSE)
    }
    if (!all(is.character(hypo_names))) {
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos, " names. \n",
           "Now, (some of) the elements are not characters.",
           call. = FALSE)
    }
    # [/CHANGE 2026-10]
  }
  
  # [CHANGE 2026-10 | audit] B2/B21: priorICweights validated and matched by name (replaces the match.arg block for order_studies)
  # Check the prior IC weights (one for each hypothesis; matched by name when
  # named; rescaled to sum to 1)
  priorICweights <- .evSyn_check_priorICweights(priorICweights, NrHypos_incl, hypo_names,
                                                input_names = input_names)
  # [/CHANGE 2026-10]
  
  
  # [CHANGE 2026-10 | audit] B21: study_weights validated
  # Check the study weights (zero weights are allowed; see .evSyn_cum_weighted)
  study_weights <- .evSyn_check_study_weights(study_weights, S, study_names)
  study_weights_S <- study_weights$study_weights_S # Now, they sum up to S
  study_weights <- study_weights$study_weights
  # [/CHANGE 2026-10]
  
  
  # [CHANGE 2026-10 | audit] B22: order_studies validated via .evSyn_order_studies
  # Check the order of the studies (character string, permutation of 1:S, or
  # permutation of the study names)
  order_studies <- .evSyn_order_studies(order_studies, S, study_names)
  
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  if (is.numeric(order_studies)) {
    # User-specified numeric order vector
    # [CHANGE 2026-10 | audit] order_studies already validated; drop = FALSE (single study)
    orderStudies <- order_studies
    Weights <- Weights[orderStudies, , drop = FALSE]
  } else if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    if (type_ev == "average") { 
      # average-evidence approach
      # [CHANGE 2026-10 | Rebecca] notes on combining IC weights with study and prior weights (average)
      # IC: Average of IC values, possibly weighted using study weights.
      # ICweights: based on IC values, but take into account possible prior hypothesis weights.
      # Notably, not possible to calculate IC values. Thus:
      # ICweights: product of IC weights, where the study and prior weights are now 'powers'.
      #            Btw Here the study weights sum to 1, because of taking the average IC values.
      #OverallGoric <- # cannot be determined now.
      #OverallPrefHypo <- which(OverallGoric == max(OverallGoric))
      # [/CHANGE 2026-10]
      # [CHANGE 2026-10 | audit] preferred hypothesis on the log scale with rescaled study weights and prior weights (same as the final results; fixes undefined OverallGoric)
      # Computed on the log scale (-2 * log of the weighted product of IC weights
      # is a difference in IC values), including possible prior hypothesis weights.
      # The study weights are rescaled in the same way as for the final results
      # below (see .evSyn_rescale_study_weights), such that the preferred 
      # hypothesis equals the finally preferred one.
      OverallICdiff <- -2 * .evSyn_log_prod_weights(log(Weights), .evSyn_rescale_study_weights(study_weights_S, type_ev))
      OverallPrefHypo <- .evSyn_pref_hypo(OverallICdiff, priorICweights)
      # [/CHANGE 2026-10]
    } else {
      # type_ev == "added" (or when "equal", because then it is overruled to be "added")
      # [CHANGE 2026-10 | Rebecca] notes on combining IC weights with study and prior weights (added)
      # IC: Sum of IC values, possibly weighted using study weights.
      # ICweights: based on IC values, but take into account possible prior hypothesis weights.
      # Notably, not possible to calculate IC values. Thus:
      # ICweights: product of IC weights, where the study and prior weights are now 'powers'.
      #            Btw Here the study weights sum to S, because of taking the average IC values.
      # [/CHANGE 2026-10]
      type_ev = "added"
      # [CHANGE 2026-10 | Rebecca] commented-out code
      #OverallGoric <- # cannot be determined now.
      #OverallPrefHypo <- which(OverallGoric == max(OverallGoric))
      # [CHANGE 2026-10 | audit] preferred hypothesis on the log scale with rescaled study weights and prior weights (fixes undefined OverallGoric)
      # Computed on the log scale (-2 * log of the weighted product of IC weights
      # is a difference in IC values), including possible prior hypothesis weights.
      # The study weights are rescaled in the same way as for the final results
      # below (see .evSyn_rescale_study_weights), such that the preferred 
      # hypothesis equals the finally preferred one.
      OverallICdiff <- -2 * .evSyn_log_prod_weights(log(Weights), .evSyn_rescale_study_weights(study_weights_S, type_ev))
      OverallPrefHypo <- .evSyn_pref_hypo(OverallICdiff, priorICweights)
      # [/CHANGE 2026-10]
    }
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    # [CHANGE 2026-10 | audit] fix: order on Weights (GORICA_weight_m does not exist in this route)
    orderStudies <- order(Weights[, OverallPrefHypo[1]], decreasing = decreasing)
    #
    # [CHANGE 2026-10 | audit] drop = FALSE (single study)
    Weights <- Weights[orderStudies, , drop = FALSE]
  }
  #
  # Set colnames
  colnames(Weights) <- hypo_names
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    # [CHANGE 2026-10 | Rebecca] study_names length checked
    if (length(study_names) != S) {
      stop("\nrestriktor ERROR: The argument 'study_names' should consist of ", S, " names, \n",
           "namely one for each study. It now consists of ", length(study_names), ".",
           call. = FALSE)
    }
    # [/CHANGE 2026-10]
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(Weights) <- study_names
  # [CHANGE 2026-10 | audit] study weights re-ordered with the studies
  # Re-order the study weights as well, such that they match the (re-ordered) 
  # studies and results do not depend on the input order.
  study_weights <- study_weights[orderStudies]
  study_weights_S <- study_weights_S[orderStudies]
  # [/CHANGE 2026-10]
  
  CumulativeWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos))
  colnames(CumulativeWeights) <- hypo_names
  sequence <- paste0("Study nr.s 1-", 1:S, "   ")
  sequence[1] <- "Study nr.  1   "
  rownames(CumulativeWeights) <- c(sequence, "Final")
  #
  # [CHANGE 2026-10 | audit] cumulative IC weights on the log scale: product of IC weights with rescaled study weights as powers, times priorICweights (replaces the sequential product; no underflow, zero study weights)
  # Cumulative IC weights: product of the IC weights of studies 1 to s, where 
  # the (rescaled) study weights are used as powers, times the prior hypothesis 
  # weights; and then normalized. This is computed on the log scale (-2 * log 
  # of the product is a difference in IC values) to avoid underflow.
  # If the study weights of studies 1 to s are all zero, there is no evidence 
  # yet and the cumulative IC weights equal the prior IC weights.
  logWeights <- log(Weights)
  for (s in seq_len(S)) {
    # Rescaled study weights of studies 1 to s (see .evSyn_rescale_study_weights):
    # - average: as if there were IC values, then average IC values;
    #            therefore, use study weights which sum to 1,
    # - added (or "equal", because then it is overruled to be "added"):
    #            as if there were IC values, then sum IC values; therefore, 
    #            use study weights which sum to the number of (positively 
    #            weighted) studies so far (not to 1).
    stW <- .evSyn_rescale_study_weights(study_weights_S[1:s], type_ev)
    CumICdiff <- -2 * .evSyn_log_prod_weights(logWeights[1:s, , drop = FALSE], stW)
    CumulativeWeights[s, ] <- ic_weights_log(CumICdiff, priorICweights)
  # [/CHANGE 2026-10]
  }
  CumulativeWeights[(S+1), ] <- CumulativeWeights[S, ]
  
  Final.weights <- CumulativeWeights[S, ]
  Final.ratio.GORICA.weights <- Final.weights %*% t(1/Final.weights)
  # [CHANGE 2026-10 | Rebecca] diagonal of the ratio matrix set to 1 (zero weights)
  diag(Final.ratio.GORICA.weights) <- 1 # If a weight is zero, then you get Inf and NaN; this way you get Inf and 1.
  
  rownames(Final.ratio.GORICA.weights) <- hypo_names
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hypo_names)
  
  # [CHANGE 2026-10 | Rebecca] priorICweights applied to the study-specific IC weights
  # Use priorICweights (i.e., a priori likeliness for each hypotheses)
  # [CHANGE 2026-10 | audit] sweep() over the columns (no recycling across the wrong dimension); B17 open question noted
  # Note: sweep() multiplies each column (hypothesis) with its own prior weight.
  # TO DO: open question (Leonard/Rebecca): should study-specific weights include priorICweights in all routes? Currently est/gorica route does not, ICvalues/ICweights routes do.
  Weights <- sweep(Weights, 2, priorICweights, "*")
  Weights <- Weights / rowSums(Weights)
  # [/CHANGE 2026-10]
  
  out <- list(type             = type,
    type_ev           = type_ev,
    #hypotheses       = hypo_names,
    # [CHANGE 2026-10 | Rebecca] priorICweights in the output
    priorICweights = priorICweights,
    n_studies         = S,
    order_studies     = orderStudies,
    study_names       = study_names,
    # [CHANGE 2026-10 | Rebecca] study_weights in the output
    study_weights = study_weights, #rep(1/S, S),
    #study_sample_nobs = study_sample_nobs,
    GORICA_weight_m            = Weights,
    # [CHANGE 2026-10 | audit] logW_m (log IC weights as input) stored for leave1studyout()
    # log of the IC weights as input (per study; without prior IC weights and
    # study weights): -2 * logW_m are differences in IC values (up to a
    # study-specific constant), which is the information needed to re-do the
    # synthesis on the log scale (e.g., in leave1studyout()) without the loss
    # of information (underflow) that recovering it from the normalised
    # (prior-weighted) weights would give.
    logW_m                     = logWeights,
    # [/CHANGE 2026-10]
    Cumulative_GORICA_weights  = CumulativeWeights,
    Final_ratio_GORICA_weights = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn_ICweights", "evSyn")
  
  return(out)
}


# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on the ratio of AIC or ORIC or GORIC or GORICA 
# weights or (Bayesian) posterior model probabilities
evSyn_ICratios <- function(object, ..., type_ev = c("added", "average"), 
                           # [CHANGE 2026-10 | Rebecca] priorICweights argument (replaces priorWeights)
                           hypo_names = c(), priorICweights = NULL, 
                           type = c("goric", "goricc", "gorica", "goricac"),
                           order_studies = c("input_order", "ascending", "descending"),
                           # [CHANGE 2026-10 | Rebecca] study_weights argument
                           study_names = c(),
                           study_weights = NULL) {
  
  if (missing(type_ev)) 
    type_ev <- "added"
  # [CHANGE 2026-10 | audit] B23: 'equal' falls back to 'added' with a message
  type_ev <- .evSyn_type_ev_IC(type_ev, "ratios of IC weights")
  
  if (missing(type)) 
    type <- "gorica"
  type <- match.arg(type)
  
  # [CHANGE 2026-10 | audit] B15 priorWeights compat; A7/B14 input checks; B13 names aligned; Href: common reference hypothesis searched, otherwise ratios rescaled to a common Href (message)
  # Backwards compatibility: 'priorWeights' is renamed to 'priorICweights'.
  priorICweights <- .evSyn_priorWeights_compat(priorICweights, list(...))

  # check the input: a list of numeric vectors without NA/NaN, with one value
  # for each hypothesis, the same hypotheses in each study
  .evSyn_check_list_input(object, name = "object")
  
  # If the input vectors carry hypothesis names, these must denote the same 
  # hypotheses across studies (matched by name); otherwise matched by position.
  object <- .evSyn_check_input_names(object, hypo_names)
  # the hypothesis names of the input (if any; accepted for a named 'priorICweights')
  input_names <- .evSyn_input_hypo_names(object)
  
  # Determine reference hypothesis -- use in output headers
  # Which hypothesis is the best, for each study
  #
  # Note that because of a previous step, each study will have at least one 1
  # So, for each study, determine for which hypothesis/-es there is a one:
  refHypo_s <- lapply(object, function(x){which(x == 1)})
  # Search for a hypothesis that has a ratio of 1 in all studies 
  # (i.e., a common reference hypothesis); take the first candidate.
  Href <- NA_integer_
  for (h in refHypo_s[[1]]) {
    if (all(vapply(refHypo_s, function(x){any(x == h)}, logical(1)))) {
      Href <- h
      break
    }
  }
  # If did not find a Href for which all studies have a one,
  # transform the input such that all have the same Href.
  # Note: this must be done before the ratios (Weights) are determined below.
  if (is.na(Href)) {
    #Href <- Mode(as.data.frame(refHypo_s)) # not base function
    Href <- if (length(refHypo_s[[1]]) > 0) refHypo_s[[1]][1] else 1L
    object <- lapply(object, function(x){x / x[Href]})
    message("\nrestriktor Message: Not all studies used the same reference hypothesis. \n",
            "Now, evidence synthesis is done using Href = ", Href, " as the reference hypothesis. \n",
            "Note that the choice for Href does not affect the ratio of final IC weights.")
  }
  Href <- unname(Href)
  
  # [/CHANGE 2026-10]
  # [CHANGE 2026-10 | Rebecca] Weights are ratios of IC weights
  Weights <- object # Now, ratio of weights 
  S <- length(Weights)
  # [CHANGE 2026-10 | audit] unname before rbind
  Weights <- do.call(rbind, lapply(Weights, unname))
  NrHypos <- ncol(Weights)
  # [CHANGE 2026-10 | Rebecca] NrHypos_incl (prior-weight code replaced by the central check below)
  NrHypos_incl <- NrHypos
  
  if (is.null(hypo_names)) {
    # [CHANGE 2026-10 | Rebecca] hypo_names validated (length and type)
    hypo_names <- paste0("H", 1:NrHypos) 
  } else {
    if (length(hypo_names) != NrHypos) {
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos, " names, \n",
           "namely one for each specified hypothesis. It now consists of ", length(hypo_names), ".",
           call. = FALSE)
    }
    if (!all(is.character(hypo_names))) {
      stop("\nrestriktor ERROR: The argument 'hypo_names' should consist of ", NrHypos, " names. \n",
           "Now, (some of) the elements are not characters.",
           call. = FALSE)
    }
    # [/CHANGE 2026-10]
  }
  
  # [CHANGE 2026-10 | audit] name of the reference hypothesis (used in the output headers)
  # Name of the reference hypothesis (used in output headers)
  names(Href) <- hypo_names[Href]
  
  
  # [CHANGE 2026-10 | audit] B2/B21: priorICweights validated and matched by name
  # Check the prior IC weights (one for each hypothesis; matched by name when
  # named; rescaled to sum to 1)
  priorICweights <- .evSyn_check_priorICweights(priorICweights, NrHypos_incl, hypo_names,
                                                input_names = input_names)
  # [/CHANGE 2026-10]
  # [CHANGE 2026-10 | Rebecca] commented-out code: priorICweights as ratios
  # # If using ratios:
  # # Note that priorICweights is now also a ratio of hypotheses weights.
  # NrHypos_incl <- NrHypos
  # if (is.null(priorICweights)) {
  #   priorICweights <- rep(1, (NrHypos_incl))
  #   # Note that priorICweights is now also a ratio of hypotheses weights.
  #   # Now all hypotheses equally likely a priori, so ratios of 1.
  # }
  # # It should have the same reference hypothesis as the ICratios do.
  # # If, then: Check is done above.
  # # Could perhaps also use weights instead of ratios...
  # # Check if length is number of hypotheses in the set
  # if (length(priorICweights) != NrHypos_incl) {
  #   stop("\nrestriktor ERROR: The argument 'priorICweights' should consist of ", NrHypos_incl, " elements, \n",
  #        "namely one for each hypothesis including a possible failsafe hypothesis. \n", 
  #        "It now consists of ", length(priorICweights), ".\n",
  #        "Note that it should contain ratios of a priori hypotheses weights, \n",
  #        "since the input consists of ratios of IC ratios; \n", 
  #        "where all should use the same reference hypothesis (leading to a ratio of 1).",
  #        call. = FALSE)
  # }
  # [/CHANGE 2026-10]
  
  # [CHANGE 2026-10 | audit] B21: study_weights validated
  # Check the study weights (zero weights are allowed; see .evSyn_cum_weighted)
  study_weights <- .evSyn_check_study_weights(study_weights, S, study_names)
  study_weights_S <- study_weights$study_weights_S # Now, they sum up to S
  study_weights <- study_weights$study_weights
  # [/CHANGE 2026-10]
  
  # [CHANGE 2026-10 | audit] B22: order_studies validated via .evSyn_order_studies
  # Check the order of the studies (character string, permutation of 1:S, or
  # permutation of the study names)
  order_studies <- .evSyn_order_studies(order_studies, S, study_names)
  
  orderStudies <- 1:S
  # Check if order of studies should be changed.
  if (is.numeric(order_studies)) {
    # User-specified numeric order vector
    # [CHANGE 2026-10 | audit] order_studies already validated; drop = FALSE (single study)
    orderStudies <- order_studies
    Weights <- Weights[orderStudies, , drop = FALSE]
  } else if (order_studies %in% c("ascending", "descending")) {
    # Order needs to be changed based on the overall preferred hypothesis.
    # Determine what the overall preferred hypothesis is.
    # [CHANGE 2026-10 | audit] preferred hypothesis from the study-weighted cumulative IC differences (log scale), as in the final results
    # That is, the hypothesis with the highest final (prior-weighted) IC weight,
    # computed on the log scale (-2 * log of the weighted product of the ratios
    # is a difference in IC values).
    # The same computation as for the final cumulative IC differences below,
    # such that the preferred hypothesis equals the finally preferred one.
    OverallICdiff <- .evSyn_cum_weighted(-2 * log(Weights), study_weights_S)[S, ]
    # [/CHANGE 2026-10]
    if (type_ev == "average") { 
      # average-evidence approach
      # [CHANGE 2026-10 | audit] average over the positively weighted studies
      OverallICdiff <- OverallICdiff / sum(study_weights_S > 0)
    } else {
      # type_ev == "added" (or when "equal", because then it is overruled to be "added")
      type_ev = "added"
    }
    # [CHANGE 2026-10 | audit] preferred hypothesis = highest prior-weighted IC weight (fixes undefined OverallGoric)
    OverallPrefHypo <- .evSyn_pref_hypo(OverallICdiff, priorICweights)
    if (order_studies == "descending") {
      decreasing = TRUE
    } else {
      decreasing = FALSE
    }
    # [CHANGE 2026-10 | audit] fix: order on the normalised ratios (GORICA_weight_m does not exist in this route)
    # Order based on the study-specific IC weights (i.e., normalized ratios)
    orderStudies <- order((Weights / rowSums(Weights))[, OverallPrefHypo[1]], decreasing = decreasing)
    #
    # [CHANGE 2026-10 | audit] drop = FALSE (single study)
    Weights <- Weights[orderStudies, , drop = FALSE]
  }
  #
  # Set colnames
  colnames(Weights) <- hypo_names
  # Set rownames (after determining the order of the studies)
  if (is.null(study_names)) {
    # If no suggested study_names, then make them
    #paste0("Study ", orderStudies)
    study_names <- as.character(orderStudies)
  } else {
    # If suggested study_names:
    # Check if length correct
    # [CHANGE 2026-10 | Rebecca] study_names length checked
    if (length(study_names) != S) {
      stop("\nrestriktor ERROR: The argument 'study_names' should consist of ", S, " names, \n",
           "namely one for each study. It now consists of ", length(study_names), ".",
           call. = FALSE)
    }
    # [/CHANGE 2026-10]
    #
    # Re-order
    study_names <- study_names[orderStudies]
  }
  rownames(Weights) <- study_names
  # [CHANGE 2026-10 | audit] study weights re-ordered with the studies
  # Re-order the study weights as well, such that they match the (re-ordered) 
  # studies and results do not depend on the input order.
  study_weights <- study_weights[orderStudies]
  study_weights_S <- study_weights_S[orderStudies]
  # [/CHANGE 2026-10]
  
  
  sequence <- paste0("Study nr.s 1-", 1:S, "   ")
  sequence[1] <- "Study nr.  1   "
  #
  # [CHANGE 2026-10 | Rebecca] ICratios output matrices: study-specific weights, cumulative ratios and cumulative IC differences
  studyspecWeights <- matrix(NA, nrow = (S), ncol = (NrHypos))
  colnames(studyspecWeights) <- hypo_names
  # [CHANGE 2026-10 | audit] rownames = study_names
  rownames(studyspecWeights) <- study_names
  #
  CumulativeRatios <- matrix(NA, nrow = (S+1), ncol = (NrHypos))
  CumulativeWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos))
  CumulativeICdiff <- matrix(NA, nrow = (S+1), ncol = (NrHypos))
  colnames(CumulativeRatios) <- colnames(CumulativeWeights) <- colnames(CumulativeICdiff) <- hypo_names
  rownames(CumulativeRatios) <- rownames(CumulativeWeights) <- rownames(CumulativeICdiff) <- c(sequence, "Final")
  #
  CumulativeRatioWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos))
  colnames(CumulativeRatioWeights) <- hypo_names
  rownames(CumulativeRatioWeights) <- c(sequence, "Final")
  #
  # [/CHANGE 2026-10]
  # [CHANGE 2026-10 | audit] IC differences from the ratios; study-specific and cumulative IC weights via log-sum-exp with priorICweights; study-weighted cumulative IC differences; cumulative ratios vs Href (replaces the sequential product)
  # The difference in IC values (vs reference hypothesis) can be determined 
  # based on the ratios; the study-specific IC weights take into account 
  # possible prior hypothesis weights.
  IC_diff <- -2 * log(Weights)
  studyspecWeights[1:S, ] <- .evSyn_IC_weights_rows(IC_diff, priorICweights)
  # Cumulative IC differences (possibly weighted using study weights):
  # - added:   sum of IC differences (also when "equal", because then it is overruled to be "added"),
  # - average: average of IC differences.
  # Cumulative IC weights: based on the cumulative IC differences, taking into 
  # account possible prior hypothesis weights (computed on the log scale).
  # Note: The priorICweights should here perhaps be ratios (vs ref. hypo) as well...
  CumulativeICdiff[1:S, ] <- .evSyn_cum_weighted(IC_diff, study_weights_S)
  if (type_ev == "average") { 
    CumulativeICdiff[1:S, ] <- CumulativeICdiff[1:S, , drop = FALSE] / pmax(.evSyn_n_pos(study_weights_S), 1)
  }
  CumulativeWeights[1:S, ] <- .evSyn_IC_weights_rows(CumulativeICdiff[1:S, , drop = FALSE], priorICweights)
  # Ratio GORIC(A) weights (vs reference hypothesis)
  CumulativeRatioWeights[1:S, ] <- CumulativeWeights[1:S, , drop = FALSE] / CumulativeWeights[1:S, Href]
  # [/CHANGE 2026-10]
  
  # [CHANGE 2026-10 | Rebecca] final row of the cumulative IC differences
  # add final row
  CumulativeICdiff[(S+1), ] <- CumulativeICdiff[S, ] 
  CumulativeWeights[(S+1), ] <- CumulativeWeights[S, ]
  # [CHANGE 2026-10 | Rebecca] final row of the cumulative ratios
  #
  CumulativeRatioWeights[(S+1), ] <- CumulativeRatioWeights[S, ]
  
  Final.weights <- CumulativeWeights[S, ]
  # [CHANGE 2026-10 | Rebecca] ratios all vs all; diagonal set to 1 (zero weights)
  # The Final.weights are the ratios of weights vs the reference hypothesis
  Final.ratio.GORICA.weights <- Final.weights %*% t(1/Final.weights) 
  diag(Final.ratio.GORICA.weights) <- 1 # If a weight is zero, then you get Inf and NaN; this way you get Inf and 1.
  
  rownames(Final.ratio.GORICA.weights) <- hypo_names
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hypo_names)
  
  out <- list(type             = type,
    type_ev           = type_ev,
    # [CHANGE 2026-10 | Rebecca] priorICweights and Href in the output
    ##hypotheses       = hypo_names,
    priorICweights = priorICweights,
    Href = Href, 
    n_studies         = S,
    order_studies     = orderStudies,
    study_names       = study_names,
    # [CHANGE 2026-10 | Rebecca] ICratios output elements (study_weights, ICdiff_m, GwRatio_m, cumulative ratios/differences)
    study_weights = study_weights, #rep(1/S, S),
    ##study_sample_nobs = study_sample_nobs,
    #GORICA_m          = IC_diff, # diff in IC values versus reference hypo
    GORICA_weight_m   = studyspecWeights, # IC weights!
    #Cumulative_GORICA = CumulativeICdiff, # cum diff in IC values versus reference hypo
    ICdiff_m    = IC_diff, # diff in IC values versus reference hypo
    GwRatio_m   = Weights, # ratio of IC weights!
    Cumulative_ratioICweights  = CumulativeRatioWeights, # cum ratios!
    Cumulative_ICdiff = CumulativeICdiff, # cum diff in IC values versus reference hypo
    Cumulative_GORICA_weights  = CumulativeWeights, # This is cum GORICA weights (so, not cum ratios!)
    Final_ratio_GORICA_weights = Final.ratio.GORICA.weights) # These are the ratios, and all versus all (as usual)
    # [/CHANGE 2026-10]
  
  class(out) <- c("evSyn_ICratios", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# list with goric objects
evSyn_gorica <- function(object, ..., type_ev = c("added", "equal", "average"), 
                         # [CHANGE 2026-10 | Rebecca] priorICweights argument
                         hypo_names = c(), priorICweights = NULL,
                         order_studies = c("input_order", "ascending", "descending"),
                         # [CHANGE 2026-10 | Rebecca] study_weights argument
                         study_names = c(),
                         study_weights = NULL) {
  
  if (missing(type_ev)) 
    type_ev <- "added"
  type_ev <- match.arg(type_ev)

  # [CHANGE 2026-10 | audit] B15: priorWeights compat (removed from dots passed to evSyn_LL())
  # Backwards compatibility: 'priorWeights' is renamed to 'priorICweights'.
  # The deprecated argument is removed from the arguments passed on to evSyn_LL().
  dots <- list(...)
  priorICweights <- .evSyn_priorWeights_compat(priorICweights, dots)
  dots$priorWeights <- NULL
  # [/CHANGE 2026-10]

  # Check if all objects are of type "con_goric"
  # [CHANGE 2026-10 | audit] empty or non-list input refused
  if (!is.list(object) || length(object) == 0 ||
      !all(vapply(object, function(x) inherits(x, "con_goric"), logical(1)))) {
    stop("\nrestriktor ERROR: the object must be a list with fitted objects from the goric() function", 
         call. = FALSE)
  }
  # Check if all objects have the same type (e.g., gorica, goric, goricc, goricac)
  object_types <- vapply(object, function(x) x$type, character(1))
  if (length(unique(object_types)) > 1) {
    stop("\nrestriktor ERROR: All goric objects must be of the same type. Found types: ",
         paste(sQuote(unique(object_types)), collapse = ", "), ".",
         call. = FALSE)
  }
  # [CHANGE 2026-10 | audit] goric objects must share comparison and penalty_factor (penalty_factor via '...' must match); A6: hypotheses aligned by name via .evSyn_align_goric_hypos
  # Check if all objects use the same comparison (unconstrained, complement, none)
  object_comparison <- vapply(object, function(x) {
    if (is.null(x$comparison)) NA_character_ else as.character(x$comparison)
  }, character(1))
  if (length(unique(object_comparison)) > 1) {
    stop("\nrestriktor ERROR: All goric objects must use the same comparison ",
         "(i.e., 'unconstrained', 'complement', or 'none'). Found: ",
         paste(sQuote(unique(object_comparison)), collapse = ", "), ".",
         call. = FALSE)
  }
  # Check if all objects use the same penalty factor (IC = -2 * LL + penalty_factor * PT)
  object_pf <- vapply(object, function(x) {
    if (is.null(x$penalty_factor)) 2 else as.numeric(x$penalty_factor)
  }, numeric(1))
  if (length(unique(object_pf)) > 1) {
    stop("\nrestriktor ERROR: All goric objects must use the same 'penalty_factor'. Found: ",
         paste(unique(object_pf), collapse = ", "), ".",
         call. = FALSE)
  }
  penalty_factor <- object_pf[1]
  # A 'penalty_factor' passed via '...' must equal the one of the goric objects
  # (it is not passed on to evSyn_LL() twice)
  if (!is.null(dots[["penalty_factor"]])) {
    if (!is.numeric(dots[["penalty_factor"]]) || length(dots[["penalty_factor"]]) != 1L ||
        !isTRUE(all.equal(as.numeric(dots[["penalty_factor"]]), penalty_factor))) {
      stop("\nrestriktor ERROR: The argument 'penalty_factor' (now, ",
           paste(deparse(dots[["penalty_factor"]], nlines = 1), collapse = ""),
           ") is taken from the goric objects (i.e., ", penalty_factor,
           ") and cannot be changed in evSyn(). Please refit the goric objects ",
           "with the requested 'penalty_factor' or do not specify it in evSyn().",
           call. = FALSE)
    }
    dots$penalty_factor <- NULL
  }
  
  # Identify the hypotheses by name (and, when available, by the hypothesis 
  # text), such that the studies are aligned to the hypothesis set of study 1:
  # the log-likelihood and penalty values are NOT taken by position.
  object <- .evSyn_align_goric_hypos(object)
  # [/CHANGE 2026-10]
  
  # TO DO als small sample, dan ook sample_nobs nodig of kan het zonder?
  
  # [CHANGE 2026-10 | audit] hypothesis names taken from study 1; B13: permuted hypo_names warned; B21: named priorICweights matched (alias: goric names)
  # Hypothesis names (incl. the possible failsafe hypothesis) from study 1
  model_names <- as.character(object[[1]]$result$model)
  if (is.null(hypo_names)) {
    hypo_names <- model_names
  } else if (length(hypo_names) == length(model_names) &&
             setequal(hypo_names, model_names) &&
             !identical(as.character(hypo_names), model_names)) {
    # 'hypo_names' are labels, applied in the order of the hypotheses of the
    # goric objects (the hypotheses are not re-ordered)
    .evSyn_warn_hypo_names_permuted(hypo_names, model_names)
  }
  # A named 'priorICweights' may use the labels in 'hypo_names' or the
  # hypothesis names of the goric objects (matched here; the length is
  # checked by evSyn_LL())
  if (!is.null(priorICweights) && !is.null(names(priorICweights)) &&
      length(hypo_names) == length(model_names)) {
    priorICweights <- .evSyn_match_weight_names(priorICweights, hypo_names,
                                                name = "priorICweights",
                                                what = "hypotheses (i.e., the column names of the output)",
                                                alias = model_names)
  }
  # [/CHANGE 2026-10]
  
  # Create a list for the evSyn_LL.list function
  conList <- list(
    type = object[[1]]$type,
    object = lapply(object, function(x) x$result$loglik),
    PT = lapply(object, function(x) x$result$penalty),
    type_ev = type_ev,
    hypo_names = hypo_names,
    # [CHANGE 2026-10 | Rebecca] priorICweights passed to evSyn_LL()
    priorICweights = priorICweights,
    order_studies = order_studies,
    study_names = study_names,
    # [CHANGE 2026-10 | audit] study_weights and penalty_factor passed to evSyn_LL()
    study_weights = study_weights,
    penalty_factor = penalty_factor
  )
  
  # Call the evSyn_LL.list function and return the result
  # [CHANGE 2026-10 | audit] dots without the deprecated priorWeights
  result <- do.call(evSyn_LL, append(conList, dots))
  # Add the type from the goric objects (evSyn_LL does not carry type)
  result$type <- object[[1]]$type
  class(result) <- c(class(result), "evSyn_gorica")
  
  return(result)
}


evSyn_escalc <- function(data, yi_col = "yi", vi_cols = "vi", 
                         study_sample_nobs = "ni",
                         type = NULL,
                         cluster_col = c("trial", "study", "author", "authors", "Trial", "Study", "Author", "Authors"),
                         outcome_col = NULL, ...) {
  
  if (is.null(outcome_col)) {
    message("\nrestriktor Message: By default (i.e., if 'outcome_col' is null), ", 
            "the function assumes that the parameter label used in the hypothesis ",
            "is 'theta' (one outcome variable).")
  }
  
  results <- extract_est_vcov_outcomes(
    data = data,
    outcome_col = outcome_col,
    yi_col = yi_col,
    vi_cols = vi_cols,
    cluster_col = cluster_col,
    type = type,
    study_sample_nobs = study_sample_nobs
  )
  
  # Access the parameter estimates and vcov blocks
  yi_list <- results$yi_list
  vcov_blocks <- results$vcov_blocks
  ni_list <- results$ni_list
  
  conList <- list(
    object = yi_list,
    VCOV = vcov_blocks,
    type = type,
    study_sample_nobs = ni_list
  )
  
  # Call the evSyn_LL.list function and return the result
  result <- do.call(evSyn_est, append(conList, list(...)))
  class(result) <- c(class(result), "evSyn_escalc")
  
  return(result)
  
}
