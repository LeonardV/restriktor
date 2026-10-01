# -------------------------------------------------------------------------
# goric_add_hc.R
#
# Implements the 'add_Hc' option of goric()/goric.default(): in addition to
# evaluating the user-specified order-restricted hypotheses against each
# other (i.e., comparison = "none", without an implicit unconstrained
# safeguard), also evaluate one of those hypotheses -- by default the
# first one, or the one named/numbered via 'add_Hc' -- against its own
# complement, and fold that complement in as one extra row of the results,
# replacing the usual implicit "unconstrained" safeguard row.
#
# This is only a sensible/interpretable set-up when every other specified
# hypothesis is nested within (a subset of) the reference hypothesis: the
# complement of the reference hypothesis is then also the complement of
# the whole set, so {H_1, ..., H_k, Complement of H_ref} is a legitimate
# partition-like comparison (much like {H_1, Hc} for a single hypothesis),
# rather than mixing in a Hc that partially overlaps one or more of the
# H_i. check_hypotheses_subset() below performs a best-effort check of
# this nesting condition and warns (but does not error/block) when it
# cannot be confirmed -- see '?goric' for this documented prerequisite.
#
# The resulting object keeps class 'con_goric', with the same $result /
# $ratio.gw / $ratio.lw / $ratio.pw / $objectList / $hypotheses_usr shape
# as an ordinary goric() output (just with the final row being "Complement
# of <hypothesis>" instead of "unconstrained"), so the existing
# print.con_goric() / summary.con_goric() / coef.con_goric() methods -- and
# benchmark() -- work on it unchanged.
# -------------------------------------------------------------------------

# Resolves 'add_Hc' (a number or a hypothesis name) against the fitted
# hypothesis list, and returns its index and name.
resolve_add_Hc <- function(add_Hc, objectnames) {
  if (is.numeric(add_Hc)) {
    ref_idx <- as.integer(add_Hc)
  } else {
    ref_idx <- which(objectnames == add_Hc)
    if (length(ref_idx) == 0) {
      stop(paste0(
        "\nrestriktor ERROR: The argument add_Hc = '", add_Hc, "' does not match ",
        "any of the hypothesis names (", paste(objectnames, collapse = ", "), ")."
      ), call. = FALSE)
    }
  }
  list(idx = ref_idx, name = objectnames[ref_idx])
}


# Best-effort check of whether each non-reference hypothesis is a subset
# of (nested within) the reference hypothesis add_Hc, i.e., whether every
# parameter vector satisfying hypothesis H_i also satisfies the reference
# hypothesis H_ref.
#
# Implemented as a linear feasibility test via quadprog's solve.QP
# (already a restriktor dependency), using a trivial positive-definite
# objective purely to probe feasibility of a linear system: H_i is NOT a
# subset of H_ref if and only if there exists a point that satisfies H_i's
# constraints while violating at least one of H_ref's constraint rows. We
# test each such "H_i and not-(row k of H_ref)" system for feasibility
# (probing just past the boundary by 'eps'); if all are infeasible, H_i is
# confirmed a subset (up to that numerical tolerance).
#
# Returns a character vector of hypothesis names that could NOT be
# confirmed as a subset of the reference hypothesis (empty if all were
# confirmed).
check_hypotheses_subset <- function(objectList, ref_idx, eps = 1e-6) {

  qp_feasible <- function(Amat, bvec, meq, p) {
    if (nrow(Amat) == 0) return(TRUE)  # no constraints -> whole space, feasible
    res <- tryCatch(
      solve.QP(Dmat = diag(p), dvec = rep(0, p),
              Amat = t(Amat), bvec = bvec, meq = meq),
      error = function(e) NULL
    )
    !is.null(res)
  }

  ref <- objectList[[ref_idx]]
  Amat_ref <- as.matrix(ref$constraints)
  bvec_ref <- ref$rhs
  meq_ref  <- ref$neq
  p <- ncol(Amat_ref)

  not_confirmed <- character(0)

  for (i in seq_along(objectList)) {
    if (i == ref_idx) next

    hi <- objectList[[i]]
    Amat_i <- as.matrix(hi$constraints)
    bvec_i <- hi$rhs
    meq_i  <- hi$neq

    violates_ref <- FALSE
    n_ref <- nrow(Amat_ref)

    for (k in seq_len(n_ref)) {
      row_k   <- Amat_ref[k, , drop = FALSE]
      rhs_k   <- bvec_ref[k]
      is_eq_k <- k <= meq_ref

      # one or two candidate "violates row k of H_ref" systems to probe:
      # an inequality row can only be violated from one side (<), an
      # equality row can be violated from either side (> or <)
      extra_rows <- if (is_eq_k) list(row_k, -row_k) else list(-row_k)
      extra_rhs  <- if (is_eq_k) list(rhs_k + eps, -rhs_k + eps) else list(-rhs_k + eps)

      for (j in seq_along(extra_rows)) {
        # H_i's own equalities first, then the candidate violation row,
        # then H_i's own inequalities (quadprog requires equalities first)
        Amat_test <- rbind(
          if (meq_i > 0) Amat_i[seq_len(meq_i), , drop = FALSE] else matrix(numeric(0), 0, p),
          extra_rows[[j]],
          if (meq_i < nrow(Amat_i)) Amat_i[(meq_i + 1):nrow(Amat_i), , drop = FALSE] else matrix(numeric(0), 0, p)
        )
        bvec_test <- c(
          if (meq_i > 0) bvec_i[seq_len(meq_i)] else numeric(0),
          extra_rhs[[j]],
          if (meq_i < nrow(Amat_i)) bvec_i[(meq_i + 1):nrow(Amat_i)] else numeric(0)
        )

        if (qp_feasible(Amat_test, bvec_test, meq_i, p)) {
          violates_ref <- TRUE
          break
        }
      }
      if (violates_ref) break
    }

    if (violates_ref) {
      not_confirmed <- c(not_confirmed, names(objectList)[i])
    }
  }

  not_confirmed
}


# Computes and folds in the "complement of add_Hc" row for goric.default()
# when the add_Hc argument is used. 'ans' is the already-built con_goric
# result for the user-specified hypotheses (comparison = "unconstrained"
# or "none"); this drops the implicit unconstrained row (if present),
# appends a "Complement of <hypothesis>" row (from running add_Hc's
# hypothesis against its complement, the same way goric() ordinarily does
# for a single hypothesis), and recomputes the weights / weight-ratios
# over this combined set.
goric_add_complement <- function(ans, add_Hc, hypotheses, objectnames,
                                 object, type, VCOV, sample_nobs,
                                 penalty_factor, priorICweights,
                                 control, debug, ldots, check_add_Hc = TRUE) {

  num_hypotheses <- length(hypotheses)
  ref <- resolve_add_Hc(add_Hc, objectnames)
  ref_idx  <- ref$idx
  ref_name <- ref$name

  # best-effort nesting check -- warn, don't block, if inconclusive.
  # Skipped (check_add_Hc = FALSE) for goric() calls made internally by
  # benchmark()/benchmark_asymp()/benchmark_means() to re-derive or
  # simulate from an *already* add_Hc-checked object: the check depends
  # only on the hypotheses' constraint structure, not on the estimates, so
  # re-running it on every bootstrap draw is redundant -- and, for the
  # numeric-estimate-vector input mode those internal calls use, has been
  # observed to raise a spurious/inconsistent extra warning that isn't
  # meaningful here, which would otherwise discard that draw entirely (see
  # goric_benchmark_utilities.R's parallel_function_asymp()).
  not_confirmed <- if (!isTRUE(check_add_Hc)) {
    character(0)
  } else {
    tryCatch(
      check_hypotheses_subset(ans$objectList, ref_idx),
      error = function(e) {
        warning(paste0(
          "\nrestriktor WARNING: could not verify whether the specified hypotheses ",
          "are subsets of ", sQuote(ref_name), " (the add_Hc reference hypothesis); ",
          "proceeding without this check. Error: ", conditionMessage(e)
        ), call. = FALSE)
        character(0)
      }
    )
  }
  if (length(not_confirmed) > 0) {
    warning(paste0(
      "\nrestriktor WARNING: the hypothes", ifelse(length(not_confirmed) == 1, "is ", "es "),
      paste(sQuote(not_confirmed), collapse = ", "),
      ifelse(length(not_confirmed) == 1, " is", " are"),
      " not confirmed to be a subset of the reference hypothesis ", sQuote(ref_name),
      " (add_Hc). The 'add_Hc' option -- comparing the specified hypotheses together ",
      "with the complement of ", sQuote(ref_name), " -- is only meaningful when every ",
      "specified hypothesis is nested within (a subset of) the reference hypothesis; ",
      "see '?goric'."
    ), call. = FALSE)
  }

  # run the reference hypothesis against its own complement, the same way
  # goric() ordinarily handles a single hypothesis
  ldots_c <- ldots[setdiff(names(ldots), c("control", "comparison", "hypotheses"))]
  complement_args <- c(
    list(object = object, hypotheses = hypotheses[ref_idx], comparison = "complement",
        type = type, VCOV = VCOV, sample_nobs = sample_nobs,
        penalty_factor = penalty_factor, control = control, debug = debug),
    ldots_c
  )
  complement_run <- do.call(goric.default, complement_args)

  compl_row <- complement_run$result[complement_run$result$model == "complement", , drop = FALSE]
  nameCompl <- paste0("Complement of ", ref_name)
  compl_row$model <- nameCompl

  # drop any implicit safeguard row (unconstrained) and keep just the
  # num_hypotheses user-specified rows, then append the complement row
  main_rows <- ans$result[seq_len(num_hypotheses), c("model", "loglik", "penalty", type)]
  compl_row <- compl_row[, c("model", "loglik", "penalty", type)]
  names(compl_row) <- names(main_rows)
  df <- rbind(main_rows, compl_row)
  rownames(df) <- NULL

  NrHypos_incl <- num_hypotheses + 1
  if (is.null(priorICweights)) {
    priorICweights <- rep(1 / NrHypos_incl, NrHypos_incl)
  } else if (length(priorICweights) != NrHypos_incl) {
    stop(paste0(
      "\nrestriktor ERROR: When using add_Hc, the argument 'priorICweights' should ",
      "consist of ", NrHypos_incl, " elements, namely one for each of the ",
      num_hypotheses, " specified hypotheses and one for the complement of ",
      sQuote(ref_name), ". It now consists of ", length(priorICweights), " elements."
    ), call. = FALSE)
  } else if (sum(priorICweights) != 1) {
    priorICweights <- priorICweights / sum(priorICweights)
    message("\nrestriktor Message: The argument 'priorICweights' should add up to 1. It has been rescaled accordingly.")
  }

  model_comparison_metrics <- calculate_model_comparison_metrics(df, priorICweights)

  df$loglik.weights  <- model_comparison_metrics$loglik_weights
  df$penalty.weights <- model_comparison_metrics$penalty_weights
  df$goric.weights   <- model_comparison_metrics$goric_weights
  names(df)[4] <- type
  names(df)[7] <- paste0(type, ".weights")
  rownames(df) <- NULL

  ans$result   <- df
  ans$ratio.gw <- model_comparison_metrics$goric_rw
  ans$ratio.pw <- model_comparison_metrics$penalty_rw
  ans$ratio.lw <- model_comparison_metrics$loglik_rw

  best_hypo <- which.max(df[, 7])
  ans$best_hypo <- best_hypo
  rownames(ans$ratio.gw)[best_hypo] <- paste0(rownames(ans$ratio.gw)[best_hypo], " (best)")
  rownames(ans$ratio.pw)[best_hypo] <- paste0(rownames(ans$ratio.pw)[best_hypo], " (best)")
  rownames(ans$ratio.lw)[best_hypo] <- paste0(rownames(ans$ratio.lw)[best_hypo], " (best)")

  # append the complement's restricted-coefficient row (already computed
  # by the complement sub-call, including any def.function() columns)
  compl_coef_row <- complement_run$ormle$b.restr["complement", , drop = FALSE]
  main_coefs <- ans$ormle$b.restr[objectnames, , drop = FALSE]
  all_cols <- union(colnames(main_coefs), colnames(compl_coef_row))
  pad <- function(m) {
    missing_cols <- setdiff(all_cols, colnames(m))
    if (length(missing_cols) > 0) {
      extra <- matrix(NA_real_, nrow = nrow(m), ncol = length(missing_cols),
                      dimnames = list(rownames(m), missing_cols))
      m <- cbind(m, extra)
    }
    m[, all_cols, drop = FALSE]
  }
  coefs <- rbind(pad(main_coefs), pad(compl_coef_row))
  rownames(coefs) <- c(objectnames, nameCompl)
  ans$ormle$b.restr <- coefs

  ans$comparison     <- "none"
  ans$priorICweights <- priorICweights
  ans$add_Hc         <- ref_name

  ans
}
