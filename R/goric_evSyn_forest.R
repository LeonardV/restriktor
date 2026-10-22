# =============================================================================
# Forest-type plot for a GORICA evidence synthesis based on (standardized)
# estimates and their covariance matrices (input type 'est_vcov' or 'escalc'):
# per study the ratio of GORICA weights (best hypothesis versus another one)
# with its confidence interval (CI; from goric_CI()), and at the bottom the
# final (pooled) ratio with its CI (from evSyn_CI()).
# =============================================================================

evSyn_forest <- function(object, vs = NULL, best = NULL, level = 0.95,
                         method_study = c("invert", "boot", "delta"),
                         method_pooled = c("boot", "delta", "mc"),
                         B = 1000L, seed = NULL, ...) {

  if (!inherits(object, "evSyn") || !inherits(object, c("evSyn_est", "evSyn_gorica")) ||
      is.null(object$goric_studies)) {
    stop("\nrestriktor ERROR: The forest plot is only available for an evSyn object based on ",
         "(standardized) estimates and their covariance matrices (input type 'est_vcov' or ",
         "'escalc') or on goric objects (input type 'gorica'); for the other input types, there is no information about the sampling variability.",
         call. = FALSE)
  }
  method_study  <- match.arg(method_study)
  method_pooled <- match.arg(method_pooled)
  if (!is.null(seed)) set.seed(seed)

  type_label <- toupper(object$type)
  type_ev <- object$type_ev
  unit_cap <- .evSyn_unit_label(type_ev)
  units    <- .evSyn_unit_label(type_ev, plural = TRUE, capital = FALSE)
  goric_list <- object$goric_studies            # input order of the studies
  S <- length(goric_list)
  hn <- colnames(object$LL_m)                   # hypothesis names (as in evSyn)
  M <- length(hn)
  if (M < 2L) {
    stop("\nrestriktor ERROR: At least two hypotheses are needed for ratios of weights.", call. = FALSE)
  }
  prior <- object$priorICweights
  if (is.null(prior) || length(prior) != M) prior <- rep(1 / M, M)
  prior <- prior / sum(prior)
  final_w <- object$Cumulative_GORICA_weights["Final", hn]

  # Best hypothesis (default: highest final weight) and the comparison hypotheses
  if (is.null(best)) best <- hn[which.max(final_w)]
  if (!best %in% hn) {
    stop("\nrestriktor ERROR: 'best' should be one of: ", paste(sQuote(hn), collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(vs)) vs <- setdiff(hn, best)
  if (!all(vs %in% hn) || best %in% vs) {
    stop("\nrestriktor ERROR: 'vs' should consist of hypotheses (other than 'best') from: ",
         paste(sQuote(hn), collapse = ", "), ".", call. = FALSE)
  }
  b <- match(best, hn)

  # Study order and labels as in the evSyn object (possibly re-ordered)
  ord <- object$order_studies
  if (is.null(ord)) ord <- seq_len(S)
  st_names <- object$study_names
  if (is.null(st_names)) st_names <- as.character(ord)
  sw <- object$study_weights
  if (is.null(sw)) sw <- rep(1 / S, S)
  sw <- sw / sum(sw)

  # Study-specific CIs (goric_CI); the study-specific goric objects are based on
  # equal prior weights, so the ratios (and their bounds) are multiplied by the
  # ratio of the priorICweights, as for the study-specific weights in evSyn.
  study_rows <- list()
  for (pos in seq_len(S)) {
    s <- ord[pos]
    g <- goric_list[[s]]
    mg <- as.character(g$result$model)       # labels in the goric object (same order)
    ci <- suppressMessages(suppressWarnings(goric_CI(g, level = level, method = method_study, B = B)))
    rd <- ci$ratio_weights
    for (k in vs) {
      kk <- match(k, hn)
      j <- which(rd$model == mg[b] & rd$vs == mg[kk])
      # ratio of the evSyn priors, divided by that of the priors used in the goric object
      pg <- g$priorICweights
      if (is.null(pg) || length(pg) != M) pg <- rep(1 / M, M)
      pr <- (prior[b] / prior[kk]) / (pg[b] / pg[kk])
      study_rows[[length(study_rows) + 1L]] <- data.frame(
        comparison = paste0(best, " vs. ", k), label = st_names[pos], kind = "study",
        ratio = rd$ratio[j] * pr, lower = rd$lower[j] * pr, upper = rd$upper[j] * pr,
        same_fit = isTRUE(rd$same_ormle[j]), study_weight = sw[s],
        stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, study_rows)

  # Pooled CIs (evSyn_CI; re-used from evSyn(..., CI = TRUE) when possible)
  pooled_note <- NULL
  ci_p <- NULL
  if (type_ev == "average") {
    pooled_note <- paste0("No CI for the final ratio, since the ", units,
                          " (multiple analyses of the same data) are dependent.")
  } else if (!is.null(object$CI) && isTRUE(all.equal(object$CI$level, level)) &&
             object$CI$method == method_pooled) {
    ci_p <- object$CI
  } else {
    ci_p <- evSyn_CI(object, level = level, method = method_pooled, B = B)
  }
  pooled_rows <- lapply(vs, function(k) {
    fin <- final_w[best] / final_w[k]
    lo <- up <- NA_real_
    if (!is.null(ci_p)) {
      j <- which(ci_p$ratio_weights$model == best & ci_p$ratio_weights$vs == k)
      lo <- ci_p$ratio_weights$lower[j]; up <- ci_p$ratio_weights$upper[j]
    }
    data.frame(comparison = paste0(best, " vs. ", k), label = "Final (pooled)", kind = "pooled",
               ratio = unname(fin), lower = lo, upper = up, same_fit = FALSE, study_weight = NA_real_,
               stringsAsFactors = FALSE)
  })
  df <- rbind(df, do.call(rbind, pooled_rows))
  df$comparison <- factor(df$comparison, levels = paste0(best, " vs. ", vs))
  # first study on top, the pooled result at the bottom
  df$label <- factor(df$label, levels = rev(c(st_names, "Final (pooled)")))
  df$point_type <- ifelse(df$kind == "pooled", "Final (pooled)",
                          ifelse(df$same_fit, paste0(unit_cap, ": same fit (ratio at a bound if nested)"),
                                 unit_cap))

  weighted <- length(unique(round(sw, 12))) > 1
  pct <- paste0(format(100 * level), "%")
  cap <- paste0("Horizontal lines: ", pct, " CIs (", units, ": goric_CI(), method = '", method_study,
                if (is.null(ci_p)) "" else paste0("'; final: evSyn_CI(), method = '", ci_p$method), "').\n",
                "Dashed line: ratio = 1 (equal support).")
  if (weighted) cap <- paste0(cap, "\nPoint size: study weight.")
  if (length(unique(prior)) > 1) cap <- paste0(cap, "\nThe ratios include the ratio of the priorICweights.")
  if (!is.null(pooled_note)) cap <- paste0(cap, "\n", pooled_note)

  shape_vals <- c(16, 1, 18)
  names(shape_vals) <- c(unit_cap, paste0(unit_cap, ": same fit (ratio at a bound if nested)"), "Final (pooled)")
  shape_vals <- shape_vals[names(shape_vals) %in% unique(df$point_type)]

  p <- ggplot(df, aes(y = .data[["label"]], x = .data[["ratio"]])) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "gray50") +
    geom_linerange(aes(xmin = .data[["lower"]], xmax = .data[["upper"]],
                       colour = .data[["kind"]]), linewidth = 0.8, na.rm = TRUE) +
    {
      if (weighted) {
        geom_point(aes(shape = .data[["point_type"]], colour = .data[["kind"]],
                       size = ifelse(.data[["kind"]] == "pooled", NA, .data[["study_weight"]])),
                   na.rm = TRUE)
      } else {
        geom_point(aes(shape = .data[["point_type"]], colour = .data[["kind"]]), size = 3)
      }
    } +
    {
      if (weighted) geom_point(data = df[df$kind == "pooled", , drop = FALSE],
                               aes(shape = .data[["point_type"]], colour = .data[["kind"]]), size = 5)
    } +
    scale_x_log10(labels = function(b) vapply(b, function(v) {
      if (is.na(v)) return(NA_character_)
      format(v, scientific = abs(log10(v)) > 5, drop0trailing = TRUE, trim = TRUE)
    }, character(1))) +
    scale_shape_manual(values = shape_vals) +
    scale_colour_manual(values = c(study = "black", pooled = "#D95F02"), guide = "none") +
    { if (weighted) scale_size_continuous(range = c(1.5, 5), guide = "none") } +
    facet_wrap(~ comparison, scales = "free_x") +
    labs(x = paste0("Ratio of ", type_label, " weights (log scale)"), y = NULL, shape = NULL,
         title = paste0("Ratio of ", type_label, " weights per ",
                        .evSyn_unit_label(type_ev, capital = FALSE), " and final (pooled) ratio"),
         subtitle = paste(type_ev, "Evidence Synthesis results"),
         caption = cap) +
    theme(
      legend.position = "bottom",
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      legend.text = element_text(size = 11),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 13),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      strip.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0, size = 9)
    )

  attr(p, "forest_data") <- df
  p
}
