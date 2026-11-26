# object of class Arima (stats::arima) --------------------------------------
#
# goric.Arima() evaluates hypotheses about the parameters of an (AR)(I)MA
# model fitted with stats::arima(), using the GORICA or (for non-seasonal
# models without differencing) the GORIC.
#
# - Coefficients can be used by their arima names: ar1, ar2, ..., ma1, ...,
#   sar1, ..., sma1, ..., intercept and the names of the xreg columns; and
#   'sigma2' (the innovation variance).
# - For an ARMA(1,1) model (order = c(1, 0, 1), no seasonal part), which is
#   equivalent to an AR(1) model with measurement error (ME),
#     x_t = phi x_{t-1} + e_t,  e_t ~ N(0, s2),
#     y_t = mu + x_t + u_t,     u_t ~ N(0, r),
#   the following parameters are derived (back-transformation; the sign
#   convention of arima, (1 + ma1 L) a_t, is used):
#     r      = -ma1 sigma2 / ar1                       (ME variance)
#     s2     = sigma2 (ar1 + ma1)(1 + ar1 ma1) / ar1   (innovation variance signal)
#     sx2    = s2 / (1 - ar1^2)                        (variance signal)
#     gamma0 = sx2 + r                                 (variance y)
#     lambda = sx2 / gamma0                            (reliability)
#     pme    = r / gamma0 = 1 - lambda                 (proportion ME variance)
#   These are 'unconstrained' estimates: r (and pme) can be negative, which
#   is needed to evaluate, e.g., the complement of r > 0.
# - The word 'delta' in a hypothesis is replaced by the value of the argument
#   'delta' (e.g., "pme < delta"); with sensitivity = TRUE, the hypotheses are
#   evaluated for all values in 'delta_values'.
# - The general machinery (covariance matrix, GORIC via order-restricted
#   maximum likelihood, sensitivity analysis, bootstrap benchmark) is in
#   goric_tsme.R.


goric.Arima <- function(object, ..., hypotheses = NULL,
                        comparison = NULL,
                        type = "gorica",
                        default_hypotheses = c("relevance", "minimum", "zero", "sign"),
                        delta = 0.1,
                        sensitivity = FALSE,
                        delta_values = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50),
                        VCOV_method = c("delta", "bootstrap"),
                        B = 1000L, seed = NULL, boot_trim = 0.025,
                        y = NULL, xreg = NULL,
                        sample_nobs = NULL,
                        debug = FALSE) {

  if (!inherits(object, "Arima")) {
    stop("\nrestriktor ERROR: the object must be of class 'Arima' (from stats::arima()).",
         call. = FALSE)
  }
  default_hypotheses <- match.arg(default_hypotheses)
  VCOV_method <- match.arg(VCOV_method)

  # the data (needed for the GORIC, the bootstrap and benchmark_tsme()) are
  # taken from the call of arima(), unless given via 'y' (and 'xreg')
  # (type = NULL: the GORIC when the data are found, otherwise the GORICA)
  need_data <- (!is.null(type) && type %in% c("goric", "goricc")) || VCOV_method == "bootstrap"
  mi <- arima_interface(object, y = y, xreg = xreg, envir = parent.frame(),
                        need_data = need_data)

  goric_tsme(mi, fit = object, ..., hypotheses = hypotheses, comparison = comparison,
             type = type, default_hypotheses = default_hypotheses, delta = delta,
             sensitivity = sensitivity, delta_values = delta_values,
             VCOV_method = VCOV_method, B = B, seed = seed, boot_trim = boot_trim,
             debug = debug, class_name = "goric_arima")
}


# Names of the derived measurement-error parameters
arima_me_names <- function() c("r", "s2", "sx2", "gamma0", "lambda", "pme")


# Back-transformation ARMA(1,1) -> AR(1) + ME; par = c(ar1, ma1, sigma2)
arima_me_transform <- function(par) {
  phi <- par[1]; theta <- par[2]; sa2 <- par[3]
  r      <- -theta * sa2 / phi
  s2     <- sa2 * (phi + theta) * (1 + theta * phi) / phi
  sx2    <- s2 / (1 - phi^2)
  gamma0 <- sa2 * (1 + 2 * phi * theta + theta^2) / (1 - phi^2)
  c(r = r, s2 = s2, sx2 = sx2, gamma0 = gamma0,
    lambda = sx2 / gamma0, pme = r / gamma0)
}


# Default hypothesis sets for one measurement-error variance (r, pme)
tsme_defaults_univariate <- function(set) {
  switch(set,
    # Is the ME negligible, relevant, or (implausibly) negative?
    relevance = list(
      hypotheses = list(H_negl = "pme > -delta; pme < delta",
                        H_ME   = "pme > delta",
                        H_neg  = "pme < -delta"),
      comparison = "none"),
    # Is the ME at least relevant (versus not)?
    minimum = list(
      hypotheses = list(H_ME = "pme > delta"),
      comparison = "complement"),
    # Is the ME variance zero, positive or negative?
    zero = list(
      hypotheses = list(H0 = "r = 0", Hpos = "r > 0", Hneg = "r < 0"),
      comparison = "none"),
    # Is the ME variance positive (versus negative)?
    sign = list(
      hypotheses = list(Hpos = "r > 0"),
      comparison = "complement")
  )
}


# Estimates and covariance matrix of an Arima object
arima_estimates <- function(object) {
  arma  <- object$arma   # c(p, q, P, Q, period, d, D)
  order <- c(p = arma[1], d = arma[6], q = arma[2],
             P = arma[3], D = arma[7], Q = arma[4], period = arma[5])

  cf   <- stats::coef(object)
  mask <- if (!is.null(object$mask)) object$mask else rep(TRUE, length(cf))
  V <- object$var.coef

  # names that cannot be used in a hypothesis (e.g., an xreg named after the
  # expression used in the call) are made syntactically valid
  nm_valid <- make.names(names(cf), unique = TRUE)
  renamed <- NULL
  if (!identical(nm_valid, names(cf))) {
    changed <- nm_valid != names(cf)
    renamed <- paste0("'", names(cf)[changed], "' -> '", nm_valid[changed], "'")
    if (!is.null(V)) {
      idx <- match(rownames(V), names(cf))
      dimnames(V) <- list(nm_valid[idx], nm_valid[idx])
    }
    names(cf) <- nm_valid
  }
  free <- names(cf)[mask]
  if (is.null(V) || length(free) == 0L) {
    stop("\nrestriktor ERROR: the arima object contains no (free) coefficients with a covariance matrix.",
         call. = FALSE)
  }
  V <- V[free, free, drop = FALSE]

  n   <- object$nobs
  sa2 <- object$sigma2
  base_est <- c(cf[free], sigma2 = sa2)
  base_V   <- matrix(0, length(base_est), length(base_est),
                     dimnames = list(names(base_est), names(base_est)))
  base_V[free, free] <- V
  base_V["sigma2", "sigma2"] <- 2 * sa2^2 / n

  # ARMA(1,1), no seasonal part, no differencing, ar1 and ma1 free
  me <- all(arma[1:4] == c(1, 1, 0, 0)) && arma[6] == 0 && arma[7] == 0 &&
        all(c("ar1", "ma1") %in% free)

  map <- if (me) {
    function(theta) c(theta, arima_me_transform(unname(theta[c("ar1", "ma1", "sigma2")])))
  } else {
    function(theta) theta
  }
  list(estimate = map(base_est), me = me, order = order, base_est = base_est,
       base_V = base_V, map = map, coef_all = cf, mask = mask, renamed = renamed)
}


# Model interface (see goric_tsme.R) for an Arima object
arima_interface <- function(object, y = NULL, xreg = NULL, envir = parent.frame(),
                            need_data = FALSE) {
  est  <- arima_estimates(object)
  if (!is.null(est$renamed)) {
    message(paste0("\nrestriktor Message: the following coefficient names are changed to ",
                   "be usable in the hypotheses: ", paste(est$renamed, collapse = ", "),
                   ". (Tip: for xreg, use a named matrix or data frame, or an object with ",
                   "the name you want, e.g., xreg = year.)"))
  }
  arma <- object$arma
  p <- arma[1]; q <- arma[2]; d <- arma[6]
  simple <- arma[3] == 0 && arma[4] == 0 && arma[6] == 0 && arma[7] == 0 && all(est$mask)
  cf <- est$coef_all
  arma_names <- c(if (p) paste0("ar", seq_len(p)), if (q) paste0("ma", seq_len(q)))
  include_mean <- "intercept" %in% names(cf)
  xreg_names <- setdiff(names(cf), c(arma_names, "intercept"))
  m <- length(xreg_names)
  method <- tryCatch(if (!is.null(object$call$method)) eval(object$call$method, envir) else "CSS-ML",
                     error = function(e) "CSS-ML")

  # data ---------------------------------------------------------------------------
  yy <- y; X <- xreg
  if (simple) {
    if (is.null(yy)) yy <- tryCatch(eval(object$call$x, envir), error = function(e) NULL)
    if (m > 0 && is.null(X)) X <- tryCatch(eval(object$call$xreg, envir), error = function(e) NULL)
    if (!is.null(yy)) yy <- as.numeric(yy)
    if (m > 0 && !is.null(X)) X <- as.matrix(X)
    ok <- !is.null(yy) && (m == 0 || (!is.null(X) && ncol(X) == m && nrow(X) == length(yy)))
    if (ok && anyNA(yy)) ok <- FALSE
    if (!ok) yy <- NULL
  }
  if (need_data && (!simple || is.null(yy))) {
    message(paste0("\nrestriktor Message: ",
                   if (!simple) paste0("the likelihood (GORIC), bootstrap and benchmark are only ",
                                       "available for non-seasonal arima models without ",
                                       "differencing and without fixed coefficients.")
                   else paste0("the data could not be retrieved from the call of arima(); ",
                               "please provide them via the arguments 'y' (and 'xreg').")))
  }

  # log-likelihood (exact, Gaussian), theta = (free coefs, sigma2) ---------------------
  loglik <- NULL
  if (simple && !is.null(yy)) {
    nn <- length(yy)
    loglik <- function(theta) {
      theta <- unname(theta)
      phi <- theta[seq_len(p)]; th <- theta[p + seq_len(q)]
      k <- p + q
      mu <- if (include_mean) theta[k + 1] else 0
      beta <- theta[k + include_mean + seq_len(m)]
      s2a <- theta[length(theta)]
      if (!is.finite(s2a) || s2a <= 0) return(-Inf)
      if (p && any(Mod(polyroot(c(1, -phi))) <= 1 + 1e-8)) return(-Inf)
      if (q && any(Mod(polyroot(c(1, th))) < 1)) return(-Inf)
      w <- yy - mu
      if (m) w <- w - drop(X %*% beta)
      mod <- stats::makeARIMA(phi, th, Delta = numeric(0))
      kl <- tryCatch(stats::KalmanLike(w, mod, nit = 0L, update = FALSE), error = function(e) NULL)
      if (is.null(kl) || !is.finite(kl$Lik) || !is.finite(kl$s2) || kl$s2 <= 0) return(-Inf)
      sumlog <- nn * (2 * kl$Lik - log(kl$s2))
      ssq <- kl$s2 * nn
      -0.5 * (nn * log(2 * pi * s2a) + sumlog + ssq / s2a)
    }
    # sanity check: the log-likelihood at the estimates should equal that of arima()
    if (method != "CSS" && is.finite(object$loglik) &&
        abs(loglik(est$base_est) - object$loglik) > 1e-3 * (1 + abs(object$loglik))) {
      message(paste0("\nrestriktor Message: the log-likelihood could not be reproduced from the ",
                     "data that were found (", round(loglik(est$base_est), 3), " versus ",
                     round(object$loglik, 3), " in the arima object); please provide the ",
                     "data via the arguments 'y' (and 'xreg'). The GORIC and bootstrap are ",
                     "not available now."))
      loglik <- NULL; yy <- NULL
    }
  }

  # simulation and refitting (bootstrap, benchmark_tsme()) -------------------------------
  simulate <- refit <- NULL
  if (!is.null(loglik)) {
    nn <- length(yy)
    simulate <- function(theta) {
      theta <- unname(theta)
      phi <- theta[seq_len(p)]; th <- theta[p + seq_len(q)]; k <- p + q
      mu <- if (include_mean) theta[k + 1] else 0
      beta <- theta[k + include_mean + seq_len(m)]
      s2a <- theta[length(theta)]
      e <- if (p + q > 0) {
        as.numeric(stats::arima.sim(model = list(ar = phi, ma = th), n = nn, sd = sqrt(s2a)))
      } else stats::rnorm(nn, sd = sqrt(s2a))
      e + mu + (if (m) drop(X %*% beta) else 0)
    }
    refit <- function(ynew) {
      f <- stats::arima(ynew, order = c(p, 0, q), include.mean = include_mean,
                        xreg = if (m) X else NULL, method = method)
      if (f$code != 0) NULL else f
    }
  }

  starts <- list()
  th0 <- est$base_est
  if (est$me) {
    s1 <- th0; s1["ma1"] <- 0                          # AR(1), no ME
    s2 <- th0; s2[c("ar1", "ma1")] <- 0                # white noise
    s3 <- th0; s3["ar1"] <- 0.5; s3["ma1"] <- -0.3     # AR(1) with ME
    starts <- list(s1, s2, s3)
  }

  no_lik_reason <- if (!simple) {
    "the likelihood is only available for non-seasonal arima models without differencing and without fixed coefficients"
  } else if (is.null(loglik)) {
    "the data could not be retrieved from the call of arima() (provide them via the argument 'y', and 'xreg' if used)"
  } else NULL

  list(theta_hat = est$base_est, V_base = est$base_V, map = est$map,
       loglik = loglik, no_lik_reason = no_lik_reason, starts = starts, simulate = simulate, refit = refit,
       refit_estimates = function(f) arima_estimates(f)$estimate,
       nobs = object$nobs, me = est$me, me_names = arima_me_names(),
       me_models = "an ARMA(1,1) model (arima(y, order = c(1, 0, 1)))",
       me_pop_default = "pme = 0",
       label = sprintf("ARIMA(%d,%d,%d) model fitted with stats::arima()%s",
                       p, d, q, if (est$me) " = AR(1) model with measurement error (ME)" else ""),
       vcov_source = "'var.coef' of the arima object (with Var(sigma2) = 2 sigma2^2 / n)",
       defaults = if (est$me) tsme_defaults_univariate else NULL,
       id_check = if (est$me) function() arima_identification_check(est) else NULL)
}


# Is the AR(1) + ME model (weakly) identified?
arima_identification_check <- function(est) {
  phi <- est$base_est[["ar1"]]; theta <- est$base_est[["ma1"]]
  V <- est$base_V
  z_phi <- phi / sqrt(V["ar1", "ar1"])
  z_sum <- (phi + theta) / sqrt(V["ar1", "ar1"] + V["ma1", "ma1"] + 2 * V["ar1", "ma1"])
  tsme_id_message(phi, z_phi, phi + theta, z_sum, "ar1", "ar1 + ma1")
}


tsme_id_message <- function(phi, z_phi, sum = NULL, z_sum = NULL,
                            lab_phi = "phi", lab_sum = NULL) {
  msgs <- character(0)
  if (is.finite(z_phi) && abs(z_phi) < 2) {
    msgs <- c(msgs, sprintf("%s (= %.3f) is not clearly different from 0 (estimate / SE = %.2f)",
                            lab_phi, phi, z_phi))
  }
  if (!is.null(z_sum) && is.finite(z_sum) && abs(z_sum) < 2) {
    msgs <- c(msgs, sprintf(paste0("%s (= %.3f) is not clearly different from 0 ",
                                   "(estimate / SE = %.2f), i.e., the series is close to ",
                                   "white noise"), lab_sum, sum, z_sum))
  }
  if (!length(msgs)) return(NULL)
  paste0("\nrestriktor Message: the measurement-error parameters are weakly identified: ",
         paste(msgs, collapse = "; and "), ". Their estimates are then not ",
         "(approximately) normally distributed and the GORICA can be overconfident. ",
         "Consider type = 'goric' (likelihood-based) and/or benchmark_tsme(), and interpret ",
         "the results with care (more time points, or a shorter time interval between ",
         "measurements, give more information).")
}
