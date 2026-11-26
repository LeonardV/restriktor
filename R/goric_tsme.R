# Engine for GORIC(A) on time-series models with measurement error ----------
#
# Used by goric.Arima() (goric_arima.R) and goric.arme() (arme.R). A model is
# described by a 'model interface' (mi), a list with:
#   theta_hat : named vector of the free base parameters (the ML estimates)
#   V_base    : covariance matrix of theta_hat
#   map       : function(theta) -> named vector of all parameters that can be
#               used in hypotheses (base parameters and derived ones)
#   loglik    : function(theta) -> log-likelihood (-Inf if invalid), or NULL
#               when the likelihood is not available (then no GORIC)
#   starts    : list of extra starting values for the constrained ML
#   simulate  : function(theta) -> a new data set from the model, or NULL
#   refit     : function(data) -> a new fitted object (same class), or NULL
#   nobs      : number of observations
#   me        : whether measurement-error parameters are available
#   me_names  : names of the measurement-error parameters
#   label     : short description of the model (for printing)
#   defaults  : function(set, delta) -> list(hypotheses, comparison), or NULL
#   id_check  : function() -> NULL or a message on weak identification
#
# goric_tsme() returns the goric() output (class con_goric) with an extra
# element 'tsme' and class 'goric_tsme' added.


goric_tsme <- function(mi, fit, ..., hypotheses = NULL, comparison = NULL,
                       type = "gorica", default_hypotheses = "relevance",
                       delta = 0.1, sensitivity = FALSE,
                       delta_values = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50),
                       VCOV_method = c("delta", "bootstrap"),
                       B = 1000L, seed = NULL, boot_trim = 0.025,
                       debug = FALSE, class_name = "goric_tsme") {

  VCOV_method <- match.arg(VCOV_method)

  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta) || delta <= 0) {
    stop("\nrestriktor ERROR: the argument 'delta' must be a single number > 0.", call. = FALSE)
  }
  if (!is.logical(sensitivity) || length(sensitivity) != 1L || is.na(sensitivity)) {
    stop("\nrestriktor ERROR: the argument 'sensitivity' must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(sensitivity) &&
      (!is.numeric(delta_values) || anyNA(delta_values) || any(delta_values <= 0))) {
    stop("\nrestriktor ERROR: 'delta_values' must be a vector of numbers > 0.", call. = FALSE)
  }

  # type: GORIC (likelihood-based) or GORICA ------------------------------------
  # default: the GORIC when the likelihood can be evaluated, otherwise the GORICA
  if (is.null(type)) {
    if (!is.null(mi$loglik)) {
      type <- "goric"
    } else {
      type <- "gorica"
      message(paste0("\nrestriktor Message: the GORIC (the default) is not available, since ",
                     if (!is.null(mi$no_lik_reason)) mi$no_lik_reason else
                       "the likelihood cannot be evaluated", "; the GORICA is used instead."))
    }
  }
  type <- tsme_check_type(type, has_lik = !is.null(mi$loglik))
  use_lik <- type %in% c("goric", "goricc")
  type_a  <- if (type == "goricc") "goricac" else if (type == "goric") "gorica" else type

  # estimates and delta-method covariance matrix -------------------------------
  est_all <- mi$map(mi$theta_hat)
  J <- tsme_jacobian(mi$map, mi$theta_hat)
  V_delta <- J %*% mi$V_base %*% t(J)
  dimnames(V_delta) <- list(names(est_all), names(est_all))

  # hypotheses -------------------------------------------------------------------
  user_hypotheses <- !is.null(hypotheses)
  if (!user_hypotheses) {
    if (is.null(mi$defaults) || !mi$me) {
      stop(paste0("\nrestriktor ERROR: default hypotheses are only available for models ",
                  "with measurement error (", mi$me_models, "). Please specify the argument ",
                  "'hypotheses' (available parameters: ",
                  paste(names(est_all), collapse = ", "), ")."), call. = FALSE)
    }
    def <- mi$defaults(default_hypotheses)
    hypotheses <- def$hypotheses
    if (is.null(comparison)) comparison <- def$comparison
  }
  if (!is.list(hypotheses)) {
    stop(paste("\nrestriktor ERROR: The 'hypotheses' argument is not a list.",
               "Please provide a list, for example, hypotheses = list(H1 = 'pme < delta')."),
         call. = FALSE)
  }
  if (is.null(names(hypotheses))) names(hypotheses) <- paste0("H", seq_along(hypotheses))

  uses_delta <- any(vapply(hypotheses, function(h) any(grepl("\\bdelta\\b", h)), logical(1)))
  hypotheses_used <- tsme_insert_delta(hypotheses, delta)

  mentions <- function(p) any(vapply(hypotheses_used, function(h)
    any(grepl(paste0("(^|[^[:alnum:]._])", gsub(".", "\\.", p, fixed = TRUE),
                     "($|[^[:alnum:]._])"), h)), logical(1)))
  if (!mi$me && length(mi$me_names) &&
      any(vapply(mi$me_names, mentions, logical(1)))) {
    stop(paste0("\nrestriktor ERROR: the measurement-error parameters (",
                paste(mi$me_names, collapse = ", "), ") are only available for ",
                mi$me_models, "."), call. = FALSE)
  }
  involved <- vapply(names(est_all), mentions, logical(1))
  if (!any(involved)) {
    stop(paste0("\nrestriktor ERROR: No parameters involved in the hypotheses. ",
                "Available parameters: ", paste(names(est_all), collapse = ", "), "."),
         call. = FALSE)
  }

  # covariance matrix -------------------------------------------------------------
  boot_info <- NULL
  if (VCOV_method == "bootstrap") {
    boot <- tsme_boot_vcov(mi, B = B, seed = seed, trim = boot_trim)
    VCOV_all <- boot$VCOV
    boot_info <- boot[c("B", "B_ok", "draws")]
    message(paste0("\nrestriktor Message: The covariance matrix of the estimates was ",
                   "obtained via a parametric bootstrap (", boot$B_ok, " of ", B,
                   " bootstrap samples could be used",
                   if (boot_trim > 0) paste0("; each parameter winsorized at its ", 100 * boot_trim,
                                             "% and ", 100 * (1 - boot_trim), "% quantiles") else "",
                   ")."))
  } else {
    VCOV_all <- V_delta
    message(paste0("\nrestriktor Message: The covariance matrix of the estimates was ",
                   "obtained from the ", mi$vcov_source,
                   if (mi$me) ", with the delta method for the derived parameters" else "",
                   ". It is the large-sample (asymptotic) covariance matrix."))
  }

  b_inv <- est_all[involved]
  V_inv <- VCOV_all[involved, involved, drop = FALSE]
  tsme_check_vcov(V_inv, names(b_inv))

  id_note <- NULL
  if (mi$me && any(involved[names(involved) %in% mi$me_names]) && !is.null(mi$id_check)) {
    id_note <- mi$id_check()
    if (!is.null(id_note)) message(id_note)
  }

  # GORICA (also used for the penalties of the GORIC) ------------------------------
  extraArgs <- list(...)
  if (use_lik && isTRUE(extraArgs$CI)) {
    message("\nrestriktor Message: CI = TRUE is not available for the GORIC of this object; it is ignored.")
    extraArgs$CI <- NULL
  }
  objectList <- c(list(object = b_inv, VCOV = V_inv, sample_nobs = mi$nobs,
                       hypotheses = hypotheses_used, comparison = comparison,
                       type = type_a, debug = debug),
                  extraArgs)
  attr(objectList$object, "class_org") <- class(fit)[1]
  objectList <- add_nr_nuis_param(objectList, type_a, function() {
    max(0, length(mi$theta_hat) - sum(involved))
  })
  if (use_lik) objectList$posthoc <- FALSE

  res <- do.call(goric.default, objectList)

  # GORIC: replace the fit part by the order-restricted maximum log-likelihood
  lik_info <- NULL
  if (use_lik) {
    psi_fun <- function(theta) mi$map(theta)[names(b_inv)]
    lik_info <- tsme_lik_fits(mi, res, psi_fun, b_inv)
    res <- tsme_rebuild(res, lik_info$loglik, lik_info$ormle, type = type)
    # post-hoc analysis (best hypothesis versus its complement), likelihood-based
    posthoc_user <- if (!is.null(extraArgs$posthoc)) extraArgs$posthoc else TRUE
    if (isTRUE(posthoc_user)) {
      res$posthoc <- tsme_posthoc(mi, res, objectList, psi_fun, b_inv, lik_info, type)
    }
  }

  # sensitivity analysis for delta ---------------------------------------------------
  sens <- NULL
  if (isTRUE(sensitivity)) {
    if (!uses_delta) {
      message(paste0("\nrestriktor Message: the sensitivity analysis is skipped, since none ",
                     "of the hypotheses contains 'delta'."))
    } else {
      sens <- tsme_sensitivity(mi, objectList, hypotheses, delta_values, use_lik, type, b_inv)
    }
  }

  # store and return --------------------------------------------------------------------
  se_all <- sqrt(pmax(diag(VCOV_all), 0))
  res$tsme <- list(
    estimates   = data.frame(estimate = est_all, se = se_all, involved = involved,
                             row.names = names(est_all)),
    VCOV        = VCOV_all,
    VCOV_method = VCOV_method,
    bootstrap   = boot_info,
    me          = mi$me,
    label       = mi$label,
    hypotheses  = hypotheses,
    hypotheses_used = hypotheses_used,
    comparison  = comparison,
    delta       = if (uses_delta) delta else NULL,
    default_hypotheses = if (!user_hypotheses) default_hypotheses else NULL,
    identification = id_note,
    sensitivity = sens,
    loglik_unrestricted = if (use_lik) lik_info$loglik_unc else NULL,
    # what is needed to redo the analysis (e.g., in benchmark_tsme())
    call_args   = list(type = type, VCOV_method = VCOV_method, B = B, boot_trim = boot_trim,
                       extra = extraArgs),
    fit         = fit,
    mi          = mi
  )
  class(res) <- unique(c(class_name, "goric_tsme", class(res)))
  res
}


# type ------------------------------------------------------------------------------
tsme_check_type <- function(type, has_lik) {
  if (!type %in% c("goric", "goricc", "gorica", "goricac")) {
    message("\nrestriktor Message: type = '", type, "' is not available; the GORICA is used.")
    return("gorica")
  }
  if (type %in% c("goric", "goricc") && !has_lik) {
    message(paste0("\nrestriktor Message: the GORIC (type = '", type, "') is not available ",
                   "for this model (the likelihood cannot be evaluated; see the help file); ",
                   "the GORICA is used."))
    return(if (type == "goricc") "goricac" else "gorica")
  }
  if (type == "goricc") {
    message("\nrestriktor Message: type = 'goricc' is not available for this object; the GORIC is used.")
    return("goric")
  }
  type
}


# Replace the word 'delta' in the hypotheses by its value
tsme_insert_delta <- function(hypotheses, delta) {
  val <- format(delta, scientific = FALSE, digits = 15)
  lapply(hypotheses, function(h) gsub("\\bdelta\\b", val, h))
}


# Central-difference Jacobian (avoids a dependency on numDeriv)
tsme_jacobian <- function(f, x) {
  f0 <- f(x)
  J  <- matrix(0, length(f0), length(x), dimnames = list(names(f0), names(x)))
  for (j in seq_along(x)) {
    h  <- 1e-6 * max(1, abs(x[j]))
    xp <- x; xm <- x
    xp[j] <- x[j] + h; xm[j] <- x[j] - h
    J[, j] <- (f(xp) - f(xm)) / (2 * h)
  }
  J
}


tsme_check_vcov <- function(V, nms) {
  if (any(!is.finite(V))) {
    stop(paste0("\nrestriktor ERROR: the covariance matrix of the parameters involved in ",
                "the hypotheses contains non-finite values. This happens when the model is ",
                "(nearly) not identified (e.g., an autoregressive parameter close to 0). ",
                "Consider VCOV_method = 'bootstrap'."), call. = FALSE)
  }
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-12 * max(abs(ev))) {
    stop(paste0("\nrestriktor ERROR: the covariance matrix of the parameters involved in ",
                "the hypotheses (", paste(nms, collapse = ", "), ") is singular. ",
                "This happens when the hypotheses combine parameters that are functions of ",
                "each other (e.g., lambda and pme), or when the model is (nearly) not ",
                "identified."), call. = FALSE)
  }
  invisible(TRUE)
}


# Parametric bootstrap of the covariance matrix ----------------------------------------
tsme_boot_vcov <- function(mi, B = 1000L, seed = NULL, trim = 0.025) {
  if (is.null(mi$simulate) || is.null(mi$refit)) {
    stop("\nrestriktor ERROR: VCOV_method = 'bootstrap' is not available for this model.",
         call. = FALSE)
  }
  if (!is.numeric(trim) || length(trim) != 1L || is.na(trim) || trim < 0 || trim >= 0.5) {
    stop("\nrestriktor ERROR: 'boot_trim' must be a number in [0, 0.5).", call. = FALSE)
  }
  est0 <- mi$map(mi$theta_hat)
  tsme_with_seed(seed, {
    draws <- lapply(seq_len(B), function(b) {
      d   <- tryCatch(mi$simulate(mi$theta_hat), error = function(e) NULL)
      if (is.null(d)) return(NULL)
      fit <- tryCatch(suppressWarnings(suppressMessages(mi$refit(d))), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      e <- tryCatch(mi$refit_estimates(fit), error = function(e) NULL)
      if (is.null(e) || length(e) != length(est0) || any(!is.finite(e))) return(NULL)
      e[names(est0)]
    })
  })
  draws <- do.call(rbind, Filter(Negate(is.null), draws))
  B_ok  <- if (is.null(draws)) 0L else nrow(draws)
  if (B_ok < max(20L, 0.5 * B)) {
    stop(paste0("\nrestriktor ERROR: only ", B_ok, " of ", B, " bootstrap samples ",
                "could be used; the covariance matrix cannot be obtained this way."),
         call. = FALSE)
  }
  # Winsorize each parameter: near non-identification, a few bootstrap
  # estimates are extreme and would otherwise dominate the covariance matrix.
  draws_w <- draws
  if (trim > 0) {
    for (j in seq_len(ncol(draws_w))) {
      q <- stats::quantile(draws_w[, j], c(trim, 1 - trim), names = FALSE)
      draws_w[, j] <- pmin(pmax(draws_w[, j], q[1]), q[2])
    }
  }
  list(VCOV = stats::cov(draws_w), B = B, B_ok = B_ok, draws = draws)
}


# Evaluate 'expr' with a seed, without changing the user's random-number stream
tsme_with_seed <- function(seed, expr) {
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  expr
}


# Order-restricted maximum likelihood ------------------------------------------------------
#
# Maximizes loglik(theta) subject to Amat %*% psi(theta) >= bvec, of which the
# first meq rows are equalities (psi = the parameters involved in the
# hypotheses, possibly nonlinear functions of theta). A quadratic-penalty
# method is used with an increasing penalty, started from several points.

tsme_constrained_ml <- function(loglik, theta_start_list, psi_fun, Amat, bvec, meq,
                                n_scale = 1) {
  Amat <- as.matrix(Amat)
  viol <- function(theta) {
    v <- drop(Amat %*% psi_fun(theta)) - bvec
    c(if (meq > 0) v[seq_len(meq)], pmin(0, v[setdiff(seq_along(v), seq_len(meq))]))
  }
  obj <- function(theta, mu) {
    ll <- loglik(theta)
    if (!is.finite(ll)) return(1e10)
    p <- viol(theta)
    if (any(!is.finite(p))) return(1e10)
    -ll + mu * sum(p^2)
  }
  best <- NULL
  for (th0 in theta_start_list) {
    th <- th0
    if (!is.finite(loglik(th))) next
    for (k in 0:8) {
      mu <- n_scale * 10^k
      o <- tryCatch(stats::optim(th, obj, mu = mu, method = "BFGS",
                                 control = list(maxit = 500, reltol = 1e-12)),
                    error = function(e) NULL)
      if (is.null(o) || o$value >= 1e10) {
        o <- tryCatch(stats::optim(th, obj, mu = mu, method = "Nelder-Mead",
                                   control = list(maxit = 4000, reltol = 1e-12)),
                      error = function(e) NULL)
      }
      if (is.null(o) || o$value >= 1e10) break
      th <- o$par
      if (max(abs(viol(th))) < 1e-7 * (1 + max(abs(bvec)))) break
    }
    ll <- loglik(th)
    v  <- max(abs(viol(th)))
    if (is.finite(ll) && v < 1e-4 * (1 + max(abs(bvec))) &&
        (is.null(best) || ll > best$loglik)) {
      best <- list(theta = th, loglik = ll, violation = v)
    }
  }
  best
}


# Log-likelihoods of all hypotheses in a goric result ----------------------------------------
tsme_lik_fits <- function(mi, res, psi_fun, b_inv) {
  n_scale <- max(1, mi$nobs)

  # unconstrained ML (also refines the estimates)
  o <- stats::optim(mi$theta_hat, function(th) {
    ll <- mi$loglik(th); if (!is.finite(ll)) 1e10 else -ll
  }, method = "BFGS", control = list(maxit = 1000, reltol = 1e-12))
  theta_unc <- if (o$value < -mi$loglik(mi$theta_hat)) o$par else mi$theta_hat
  ll_unc <- mi$loglik(theta_unc)
  starts <- c(list(theta_unc), mi$starts)
  psi_unc <- psi_fun(theta_unc)

  models <- res$result$model
  hyp_names <- names(res$constraints)
  ll <- setNames(rep(NA_real_, length(models)), models)
  ormle <- matrix(NA_real_, length(models), length(b_inv),
                  dimnames = list(models, names(b_inv)))
  fail <- character(0)

  fit_h <- function(A, b, meq) {
    v <- drop(A %*% psi_unc) - b
    feasible <- (meq == 0 || all(abs(v[seq_len(meq)]) < 1e-10)) &&
                all(v[setdiff(seq_along(v), seq_len(meq))] >= 0)
    if (feasible) return(list(theta = theta_unc, loglik = ll_unc))
    tsme_constrained_ml(mi$loglik, starts, psi_fun, A, b, meq, n_scale)
  }

  for (h in hyp_names) {
    A <- res$constraints[[h]]; b <- res$rhs[[h]]; meq <- res$neq[[h]]
    f <- fit_h(A, b, meq)
    if (is.null(f)) { fail <- c(fail, h); next }
    ll[h] <- f$loglik
    ormle[h, ] <- psi_fun(f$theta)
  }

  if ("unconstrained" %in% models) {
    ll["unconstrained"] <- ll_unc
    ormle["unconstrained", ] <- psi_unc
  }
  if ("complement" %in% models) {
    h  <- hyp_names[1]
    A  <- res$constraints[[h]]; b <- res$rhs[[h]]; meq <- res$neq[[h]]
    cm <- tsme_complement_fit(mi, A, b, meq, psi_unc, theta_unc, ll_unc, starts,
                              psi_fun, n_scale)
    ll["complement"] <- cm$loglik
    ormle["complement", ] <- psi_fun(cm$theta)
  }
  if (length(fail)) {
    warning(paste0("\nrestriktor WARNING: the order-restricted maximum likelihood could not ",
                   "be obtained for: ", paste(fail, collapse = ", "), "."), call. = FALSE)
  }
  list(loglik = ll, ormle = ormle, loglik_unc = ll_unc, theta_unc = theta_unc,
       starts = starts)
}


# Complement of a hypothesis: if the unconstrained estimate is in the
# hypothesis (and it has no equalities), the maximum over the complement lies
# on its boundary: the best fit with one of the restrictions as an equality.
tsme_complement_fit <- function(mi, A, b, meq, psi_unc, theta_unc, ll_unc, starts,
                                psi_fun, n_scale) {
  v <- drop(A %*% psi_unc) - b
  in_h <- meq == 0 && all(v >= 0)
  if (!in_h) return(list(theta = theta_unc, loglik = ll_unc))
  fits <- lapply(seq_len(nrow(A)), function(k)
    tsme_constrained_ml(mi$loglik, starts, psi_fun, A[k, , drop = FALSE], b[k], 1L, n_scale))
  fits <- Filter(Negate(is.null), fits)
  fits[[which.max(vapply(fits, `[[`, numeric(1), "loglik"))]]
}


# Rebuild a goric result with likelihood-based log-likelihoods -----------------------------------
tsme_rebuild <- function(res, ll, ormle, type = "goric") {
  r  <- res$result
  PT <- r$penalty
  pf <- if (!is.null(res$penalty_factor)) res$penalty_factor else 2
  ic <- -2 * ll + pf * PT
  prior <- if (!is.null(res$priorICweights)) res$priorICweights else rep(1 / length(ll), length(ll))
  w_of <- function(x) { x <- x - max(x); e <- exp(x); e / sum(e) }
  gw <- w_of(-ic / 2 + log(prior))
  lw <- w_of(ll)
  pw <- w_of(-PT)
  out <- data.frame(model = r$model, loglik = unname(ll), penalty = PT, goric = unname(ic),
                    loglik.weights = lw, penalty.weights = pw, goric.weights = gw,
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  res$result <- out
  ratio <- function(w) { m <- outer(w, w, "/"); dimnames(m) <- list(r$model, paste0("vs. ", r$model)); m }
  res$ratio.gw <- ratio(gw)
  res$ratio.lw <- ratio(lw)
  res$ratio.pw <- ratio(pw)
  res$best_hypo <- which.max(gw)
  if (!is.null(res$ormle)) res$ormle$b.restr <- ormle
  res$type <- type
  class(res) <- c("con_goric")
  res
}


# Post-hoc analysis for the GORIC: the best hypothesis versus its complement
tsme_posthoc <- function(mi, res, objectList, psi_fun, b_inv, lik_info, type) {
  best <- res$result$model[res$best_hypo]
  not_applicable <- function(reason) structure(
    list(applicable = FALSE, reason = reason, best = best, type = type), class = "goric_posthoc")
  if (best %in% c("unconstrained", "complement")) {
    return(not_applicable(paste0("The post-hoc analysis is not applicable, since the best ",
                                 "hypothesis is the ", best, " hypothesis.")))
  }
  if (identical(res$comparison, "complement")) {
    return(not_applicable("The output already contains the hypothesis versus its complement."))
  }
  ol <- objectList
  ol$hypotheses <- objectList$hypotheses[best]
  ol$comparison <- "complement"
  ol$posthoc <- FALSE
  g <- tryCatch(suppressMessages(do.call(goric.default, ol)), error = function(e) NULL)
  if (is.null(g) || !"complement" %in% g$result$model) {
    return(not_applicable("The post-hoc analysis could not be performed for this hypothesis."))
  }
  A <- g$constraints[[best]]; b <- g$rhs[[best]]; meq <- g$neq[[best]]
  psi_unc <- psi_fun(lik_info$theta_unc)
  cm <- tsme_complement_fit(mi, A, b, meq, psi_unc, lik_info$theta_unc, lik_info$loglik_unc,
                            lik_info$starts, psi_fun, max(1, mi$nobs))
  ll <- c(lik_info$loglik[[best]], cm$loglik)
  orm <- rbind(res$ormle$b.restr[best, ], psi_fun(cm$theta))
  dimnames(orm) <- list(c(best, "complement"), names(b_inv))
  g2 <- tsme_rebuild(g, ll, orm, type = type)
  structure(list(applicable = TRUE, reason = NULL, best = best, type = type,
                 hypothesis = objectList$hypotheses[[best]],
                 comparison = "complement", result = g2$result,
                 ratio.gw = g2$result$goric.weights[1] / g2$result$goric.weights[2],
                 error_prob = g2$result$goric.weights[2],
                 priorICweights = c(0.5, 0.5), priorICweights_rescaled = FALSE,
                 goric = g2),
            class = "goric_posthoc")
}


# Sensitivity analysis: the weights for several values of delta -----------------------------------
tsme_sensitivity <- function(mi, objectList, hypotheses, delta_values, use_lik, type, b_inv) {
  psi_fun <- function(theta) mi$map(theta)[names(b_inv)]
  rows <- lapply(delta_values, function(d) {
    ol <- objectList
    ol$hypotheses <- tsme_insert_delta(hypotheses, d)
    ol$posthoc <- FALSE
    r <- tryCatch(suppressMessages(do.call(goric.default, ol)), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    if (use_lik) {
      li <- tryCatch(suppressWarnings(tsme_lik_fits(mi, r, psi_fun, b_inv)), error = function(e) NULL)
      if (is.null(li)) return(NULL)
      r <- tsme_rebuild(r, li$loglik, li$ormle, type = type)
    }
    wcol <- grep("^goric(a)?(c)?\\.weights$", names(r$result), value = TRUE)[1]
    w <- setNames(r$result[[wcol]], r$result$model)
    data.frame(delta = d, t(w), check.names = FALSE)
  })
  sens <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(sens)) return(NULL)
  wcols <- setdiff(names(sens), "delta")
  sens$preferred <- wcols[max.col(as.matrix(sens[, wcols, drop = FALSE]), ties.method = "first")]
  rownames(sens) <- NULL
  sens
}


# Printing -------------------------------------------------------------------------------------------
print.goric_tsme <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  a <- x$tsme
  cat("\n", a$label, "\n", sep = "")
  cat("\nEstimates (", if (a$VCOV_method == "bootstrap") "bootstrap" else "delta-method",
      " standard errors); * = involved in the hypotheses:\n", sep = "")
  tab <- a$estimates
  out <- data.frame(estimate = format(round(tab$estimate, digits), nsmall = digits),
                    se = format(round(tab$se, digits), nsmall = digits),
                    " " = ifelse(tab$involved, "*", ""), check.names = FALSE,
                    row.names = rownames(tab))
  print(out, quote = FALSE)
  if (a$me) {
    cat("  (r = ME variance; s2 = innovation variance of the signal; sx2 = variance of the\n",
        "  signal; gamma0 = variance of y; lambda = reliability; pme = r / gamma0)\n", sep = "")
  }
  if (!is.null(a$delta)) cat("\n'delta' in the hypotheses is set to ", a$delta, ".\n", sep = "")
  if (!is.null(a$identification)) {
    cat("\nNote:", sub("^\nrestriktor Message: ", "", a$identification), "\n")
  }
  if (x$type %in% c("goric", "goricc")) {
    cat("\nNote: the log-likelihood values are the order-restricted maximum log-likelihoods ",
        "of the model; the penalties are those of the GORICA (based on the covariance ",
        "matrix of the parameters in the hypotheses).\n", sep = "")
  }

  NextMethod()

  if (!is.null(a$sensitivity)) {
    cat("\nSensitivity analysis: weights for several values of delta\n")
    s <- a$sensitivity
    num <- vapply(s, is.numeric, logical(1)) & names(s) != "delta"
    s[num] <- lapply(s[num], function(v) format(round(v, 3), nsmall = 3))
    print(s, row.names = FALSE, quote = FALSE)
  }
  invisible(x)
}


# Benchmark based on a parametric bootstrap from the model -----------------------------------------
#
# Simulates data sets from the fitted model under a population restriction
# (e.g., "pme = 0": no ME), refits the model and redoes the analysis. The
# observed weights are compared with this bootstrap distribution, which, unlike
# benchmark(), does not assume that the estimates are normally distributed.

benchmark_tsme <- function(object, pop = NULL, B = 500L, seed = NULL,
                           probs = c(0.05, 0.35, 0.50, 0.65, 0.95)) {
  if (!inherits(object, "goric_tsme")) {
    stop("\nrestriktor ERROR: the object must be the output of goric() for an arima or arme object.",
         call. = FALSE)
  }
  a  <- object$tsme
  mi <- a$mi
  if (is.null(mi$simulate) || is.null(mi$refit) || is.null(mi$loglik)) {
    stop("\nrestriktor ERROR: benchmark_tsme() is not available for this model.", call. = FALSE)
  }
  if (is.null(pop)) {
    if (!a$me) stop("\nrestriktor ERROR: please specify 'pop', e.g., pop = 'ar1 = 0.5'.", call. = FALSE)
    pop <- paste(mi$me_pop_default, collapse = "; ")
  }

  # population parameters: the ML estimates under the restriction(s) in 'pop'
  est_all <- mi$map(mi$theta_hat)
  J <- tsme_jacobian(mi$map, mi$theta_hat)
  V_all <- J %*% mi$V_base %*% t(J)
  dimnames(V_all) <- list(names(est_all), names(est_all))
  pop_par <- names(est_all)[vapply(names(est_all), function(p)
    grepl(paste0("(^|[^[:alnum:]._])", gsub(".", "\\.", p, fixed = TRUE),
                 "($|[^[:alnum:]._])"), pop), logical(1))]
  if (!length(pop_par)) stop("\nrestriktor ERROR: no parameters found in 'pop'.", call. = FALSE)
  g0 <- suppressMessages(goric.default(est_all[pop_par],
                                       VCOV = V_all[pop_par, pop_par, drop = FALSE],
                                       hypotheses = list(Hpop = pop), comparison = "none",
                                       type = "gorica", posthoc = FALSE))
  psi_fun <- function(theta) mi$map(theta)[pop_par]
  f0 <- tsme_constrained_ml(mi$loglik, c(list(mi$theta_hat), mi$starts), psi_fun,
                            g0$constraints$Hpop, g0$rhs$Hpop, g0$neq$Hpop, max(1, mi$nobs))
  if (is.null(f0)) stop("\nrestriktor ERROR: the model could not be fitted under 'pop'.", call. = FALSE)
  theta_pop <- f0$theta

  # rerun the analysis on bootstrap samples
  args <- c(list(hypotheses = a$hypotheses, comparison = a$comparison,
                 type = a$call_args$type, delta = if (is.null(a$delta)) 0.1 else a$delta,
                 VCOV_method = "delta", posthoc = FALSE),
            a$call_args$extra[setdiff(names(a$call_args$extra), "posthoc")])
  obs_w <- object$result[[grep("^goric(a)?(c)?\\.weights$", names(object$result))[1]]]
  names(obs_w) <- object$result$model

  W <- tsme_with_seed(seed, {
    lapply(seq_len(B), function(b) {
      d   <- tryCatch(mi$simulate(theta_pop), error = function(e) NULL)
      if (is.null(d)) return(NULL)
      fit <- tryCatch(suppressWarnings(suppressMessages(mi$refit(d))), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      r <- tryCatch(suppressWarnings(suppressMessages(do.call(goric, c(list(fit), args)))),
                    error = function(e) NULL)
      if (is.null(r)) return(NULL)
      w <- r$result[[grep("^goric(a)?(c)?\\.weights$", names(r$result))[1]]]
      names(w) <- r$result$model
      w[names(obs_w)]
    })
  })
  W <- do.call(rbind, Filter(Negate(is.null), W))
  if (is.null(W) || nrow(W) < 10) stop("\nrestriktor ERROR: too few usable bootstrap samples.", call. = FALSE)

  pct <- vapply(names(obs_w), function(h) 100 * mean(W[, h] <= obs_w[[h]]), numeric(1))
  quant <- t(apply(W, 2, stats::quantile, probs = probs, na.rm = TRUE))
  structure(list(pop = pop, pop_estimates = mi$map(theta_pop)[union(pop_par, names(est_all))],
                 B = B, B_ok = nrow(W), observed = obs_w, percentile = pct,
                 quantiles = quant, draws = W, type = object$type),
            class = "benchmark_tsme")
}


print.benchmark_tsme <- function(x, digits = 3, ...) {
  lab <- if (x$type %in% c("goric", "goricc")) "GORIC" else "GORICA"
  cat("\nParametric-bootstrap benchmark of the ", lab, " weights\n", sep = "")
  cat("Population: the model fitted under '", x$pop, "' (", x$B_ok, " of ", x$B,
      " bootstrap samples used)\n\n", sep = "")
  tab <- data.frame(observed = round(x$observed, digits),
                    percentile = round(x$percentile, 1),
                    round(x$quantiles, digits), check.names = FALSE)
  print(tab)
  cat("\n'percentile' = percentage of bootstrap samples (under the population) with a\n",
      "weight smaller than or equal to the observed weight. E.g., a percentile of 97 for\n",
      "the hypothesis with ME means that such a high weight is rare (3%) without ME.\n", sep = "")
  invisible(x)
}
