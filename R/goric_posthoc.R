# Post-hoc analysis: the best hypothesis versus its complement ---------------
#
# After a goric(a) analysis with a set of hypotheses (e.g., several
# order-restricted hypotheses plus the unconstrained safeguard), the
# hypothesis with the highest GORIC(A) weight is the best in that set. This
# function evaluates that best hypothesis also against its own complement
# (by re-running goric() with only that hypothesis and comparison =
# "complement"). This gives:
#  - the ratio of GORIC(A) weights of the best hypothesis versus its
#    complement, and
#  - the conditional error probability, that is, the probability (given
#    the data) that choosing the best hypothesis is incorrect (i.e., the
#    GORIC(A) weight of the complement).
# Since the hypothesis is selected based on the data, this is a post-hoc
# (exploratory) analysis, which can be used as information for future
# research.
#
# The post-hoc analysis is performed (i.e., is applicable) when the best
# hypothesis is one of the specified hypotheses; not when the best is the
# unconstrained hypothesis, a complement, or Heq, and not when the output
# already contains the (single) hypothesis versus its complement
# (comparison = "complement").
#
# It can be called directly (goric_posthoc(goric_object)) or via goric()'s
# argument posthoc, in which case the result is stored in
# goric_object$posthoc and shown by print() and summary(). By default
# (posthoc = TRUE), goric() performs it only when there is more than one
# specified hypothesis and one of these is the best, and otherwise skips it
# silently; when the user explicitly specifies posthoc = TRUE, the reason
# for not performing it is reported. posthoc = FALSE switches it off (as
# done in all internal goric() calls, e.g., in benchmark() and evSyn()).

goric_posthoc <- function(object, priorICweights = NULL, ...) {

  if (!inherits(object, "con_goric")) {
    stop(paste0("\nrestriktor ERROR: the object should be of class 'con_goric' ",
                "(i.e., output from the goric() function). It now belongs to class(es): ",
                paste(class(object), collapse = ", "), "."), call. = FALSE)
  }
  if (!is.null(object$rerun_args)) {
    rerun_args <- object$rerun_args
  } else {
    stop(paste0("\nrestriktor ERROR: the goric object does not contain the information ",
                "needed for the post-hoc analysis (element 'rerun_args'). Please re-run ",
                "goric() with the current version of restriktor."), call. = FALSE)
  }

  type       <- object$type
  result     <- object$result
  comparison <- object$comparison
  Heq        <- isTRUE(object$Heq)

  # best hypothesis, based on the IC weights (as in the print output)
  best_idx  <- which.max(result[, 7])
  best_name <- as.character(result$model[best_idx])

  # names of the user-specified hypotheses (Heq excluded)
  hypo_names <- setdiff(object$objectNames, "Heq")

  out <- list(applicable = FALSE, reason = NULL, best = best_name, type = type)
  class(out) <- "goric_posthoc"

  # When is the post-hoc analysis (not) applicable?
  if (!(best_name %in% hypo_names)) {
    if (best_name == "unconstrained") {
      out$reason <- paste0("The unconstrained hypothesis is the best in the set. That is, ",
                           "none of the specified hypotheses is supported (enough); hence, ",
                           "a post-hoc analysis of the best hypothesis versus its complement ",
                           "is not meaningful.")
    } else if (best_name == "Heq") {
      out$reason <- paste0("The equality-restricted hypothesis Heq is the best in the set; ",
                           "hence, a post-hoc analysis of the best (order-restricted) hypothesis ",
                           "versus its complement is not performed.")
    } else {
      out$reason <- paste0("The best hypothesis in the set (", sQuote(best_name), ") is a ",
                           "complement; hence, a post-hoc analysis of the best (specified) ",
                           "hypothesis versus its complement is not performed.")
    }
    return(out)
  }
  if (comparison == "complement") {
    out$reason <- paste0("The output already contains the evaluation of the hypothesis ",
                         "versus its complement; hence, no post-hoc analysis is needed.")
    return(out)
  }

  # (equality-restricted only?) If so, its complement is the unconstrained
  # hypothesis (as in goric()).
  best_pos <- match(best_name, object$objectNames)
  ceq_only <- isTRUE(NROW(object$constraints[[best_pos]]) == object$neq[[best_pos]])

  # prior IC weights for the best hypothesis and its complement
  priorICweights_rescaled <- FALSE
  if (!is.null(priorICweights)) {
    check_priorICweights_values(priorICweights)
    if (length(priorICweights) != 2L) {
      stop(paste0("\nrestriktor ERROR: The argument 'priorICweights' should consist of 2 ",
                  "elements, namely one for the best hypothesis and one for its complement. ",
                  "It now consists of ", length(priorICweights), " elements."), call. = FALSE)
    }
    priorICweights <- priorICweights / sum(priorICweights)
  } else {
    pr <- object$priorICweights
    if (!is.null(pr) && length(pr) == nrow(result) && any(abs(pr - pr[1]) > 1e-12)) {
      # unequal prior weights in the set: the best hypothesis keeps its own
      # prior weight, its complement gets the rest
      priorICweights <- c(pr[best_idx], 1 - pr[best_idx])
      priorICweights_rescaled <- TRUE
    } else {
      priorICweights <- c(0.5, 0.5)
    }
  }

  # re-run goric() with the best hypothesis versus its complement
  hypotheses_best <- rerun_args$hypotheses[best_pos - as.integer(Heq)]
  names(hypotheses_best) <- best_name
  if (is.null(rerun_args$object)) {
    obj <- object$model.org
    class(obj) <- setdiff(class(obj), "goric")
  } else {
    obj <- rerun_args$object
  }
  ldots <- rerun_args$ldots
  ldots <- ldots[setdiff(names(ldots), c("control", "comparison", "hypotheses", "posthoc",
                                         "CI", "add_Hc", "Heq", "priorICweights"))]
  args <- c(list(object = obj, hypotheses = hypotheses_best,
                 comparison = if (ceq_only) "unconstrained" else "complement",
                 type = type, VCOV = rerun_args$VCOV, sample_nobs = rerun_args$sample_nobs,
                 penalty_factor = object$penalty_factor, Heq = FALSE,
                 priorICweights = priorICweights,
                 nr_nuis_param = rerun_args$nr_nuis_param,
                 posthoc = FALSE, CI = FALSE,
                 control = rerun_args$control, debug = FALSE),
            ldots)
  # (messages are suppressed: they would repeat those of the original analysis,
  # e.g., that the complement of an equality-restricted hypothesis is the
  # unconstrained hypothesis)
  res <- suppressMessages(do.call(goric.default, args))

  out$applicable     <- TRUE
  out$hypothesis     <- if (!is.null(object$hypotheses_usr)) object$hypotheses_usr[[best_name]] else NULL
  out$comparison     <- res$comparison # "complement" (or "unconstrained" if equality-restricted only)
  out$result         <- res$result
  out$ratio.gw       <- unname(res$ratio.gw[1, 2])
  out$error_prob     <- unname(res$result[2, 7])
  out$priorICweights <- priorICweights
  out$priorICweights_rescaled <- priorICweights_rescaled
  out$goric          <- res

  out
}


print.goric_posthoc <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  IC_lab <- toupper(x$type)
  cat("\nPost-hoc analysis (information for future research): the best hypothesis versus its complement\n")
  if (!isTRUE(x$applicable)) {
    cat(paste0(x$reason, "\n"))
    return(invisible(x))
  }

  dig <- paste0("%6.", digits, "f")
  df <- x$result[, 1:7]
  df[, -1] <- lapply(df[, -1], function(v) sprintf(dig, v))
  print(df, print.gap = 2, quote = FALSE, right = TRUE, row.names = FALSE)

  compl_lab <- if (identical(x$comparison, "unconstrained")) {
    "the unconstrained hypothesis (which is the complement of a hypothesis with only equality restrictions)"
  } else {
    "its complement"
  }
  ratio <- x$ratio.gw
  ratio_txt <- if (ratio > 10000) sprintf("%.2e", ratio) else sprintf("%.2f", ratio)
  if (ratio >= 1) {
    cat(paste0("The best hypothesis in the set, ", sQuote(x$best), ", has ", ratio_txt,
               " times more support than ", compl_lab, ".\n"))
  } else {
    inv <- 1 / ratio
    inv_txt <- if (inv > 10000) sprintf("%.2e", inv) else sprintf("%.2f", inv)
    cat(paste0("Although ", sQuote(x$best), " is the best hypothesis in the set, ",
               compl_lab, " has ", inv_txt, " times more support than ", sQuote(x$best), ".\n"))
  }
  cat(paste0("The conditional error probability, that is, the probability (given the data) ",
             "that choosing ", sQuote(x$best), " (over ", compl_lab, ") is incorrect, is ",
             sprintf("%.3f", x$error_prob), ".\n"))
  cat(paste0("Note: The hypothesis is selected based on the data; hence, this is a post-hoc ",
             "(exploratory) result, which can be used to inform future research, ",
             "but not as a confirmatory evaluation.\n"))
  if (isTRUE(x$priorICweights_rescaled)) {
    cat(paste0("Note: The prior ", IC_lab, " weights (priorICweights) used here are ",
               paste(sprintf("%.3f", x$priorICweights), collapse = " and "),
               ", namely the prior weight of ", sQuote(x$best),
               " in the original set and 1 minus that value for ", compl_lab, ".\n"))
  }
  invisible(x)
}
