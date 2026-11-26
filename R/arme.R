# AR(p) and VAR(1) models with measurement error ---------------------------
#
# arme() fits, by exact Gaussian maximum likelihood (Kalman filter),
#   univariate:   x_t = phi_1 x_{t-1} + ... + phi_p x_{t-p} + e_t,  e_t ~ N(0, s2)
#                 y_t = mu + x_t + u_t,                             u_t ~ N(0, r)
#   multivariate: x_t = Phi x_{t-1} + e_t,  e_t ~ N(0, Sigma)
#                 y_t = mu + x_t + u_t,     u_t ~ N(0, diag(r_1, ..., r_k))
# The measurement-error (ME) variances r are not restricted to be positive
# (only the implied covariance matrix of y must be positive definite), such
# that the complement of 'r > 0' can be evaluated with goric().
#
# goric.arme() evaluates hypotheses about the parameters (incl. the derived
# reliabilities 'lambda' and proportions of ME variance 'pme'); see
# goric_tsme.R for the general machinery.


arme <- function(y, p = 1L, start = NULL, control = list()) {
  cl <- match.call()
  if (is.data.frame(y)) y <- as.matrix(y)
  multi <- is.matrix(y) && ncol(y) > 1
  if (multi) {
    if (p != 1L) stop("\nrestriktor ERROR: for multivariate data, only p = 1 (VAR(1)) is available.",
                      call. = FALSE)
    Y <- y
    if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
    if (anyNA(Y)) stop("\nrestriktor ERROR: missing values are not allowed.", call. = FALSE)
    spec <- arme_spec_var(ncol(Y))
    data <- Y
  } else {
    yy <- as.numeric(y)
    if (anyNA(yy)) stop("\nrestriktor ERROR: missing values are not allowed.", call. = FALSE)
    if (!is.numeric(p) || length(p) != 1L || p < 1) {
      stop("\nrestriktor ERROR: 'p' must be a positive integer.", call. = FALSE)
    }
    spec <- arme_spec_uni(as.integer(p))
    data <- yy
  }
  ll <- function(theta) spec$loglik(theta, data)

  starts <- if (!is.null(start)) list(start) else spec$starts(data)
  fits <- lapply(starts, function(th0) {
    if (!is.finite(ll(th0))) return(NULL)
    o <- tryCatch(stats::optim(th0, function(th) { v <- ll(th); if (is.finite(v)) -v else 1e10 },
                               method = "BFGS",
                               control = utils::modifyList(list(maxit = 2000, reltol = 1e-12), control)),
                  error = function(e) NULL)
    if (is.null(o) || o$value >= 1e10) return(NULL)
    o
  })
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("\nrestriktor ERROR: the model could not be fitted.", call. = FALSE)
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  theta <- best$par
  names(theta) <- spec$names

  H <- tryCatch(stats::optimHess(theta, function(th) -ll(th)), error = function(e) NULL)
  V <- if (!is.null(H)) tryCatch(solve(H), error = function(e) NULL) else NULL
  if (is.null(V) || any(!is.finite(V)) ||
      min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    warning(paste0("\nrestriktor WARNING: the covariance matrix of the estimates could not be ",
                   "obtained (the model may not be identified, e.g., an autoregressive ",
                   "parameter close to 0)."), call. = FALSE)
    V <- matrix(NA_real_, length(theta), length(theta))
  }
  dimnames(V) <- list(spec$names, spec$names)

  structure(list(coefficients = theta, vcov = V, loglik = -best$value,
                 nobs = if (multi) nrow(data) else length(data),
                 derived = spec$derived(theta), p = p, k = spec$k,
                 data = data, spec = spec, convergence = best$convergence, call = cl),
            class = "arme")
}


# Specification: univariate AR(p) + ME ------------------------------------------------------
arme_spec_uni <- function(p) {
  phin <- if (p == 1) "phi" else paste0("phi", seq_len(p))
  nms  <- c("mu", phin, "s2", "r")
  unpack <- function(theta) {
    theta <- unname(theta)
    list(mu = theta[1], phi = theta[1 + seq_len(p)], s2 = theta[p + 2], r = theta[p + 3])
  }
  state <- function(phi, s2) {
    Tm <- matrix(0, p, p); Tm[1, ] <- phi
    if (p > 1) Tm[cbind(2:p, 1:(p - 1))] <- 1
    V <- matrix(0, p, p); V[1, 1] <- s2
    P0 <- matrix(solve(diag(p^2) - kronecker(Tm, Tm), c(V)), p)
    list(T = Tm, V = V, P0 = P0)
  }
  valid <- function(u) is.finite(u$s2) && u$s2 > 0 &&
    all(Mod(polyroot(c(1, -u$phi))) > 1 + 1e-8)
  loglik <- function(theta, y) {
    u <- unpack(theta)
    if (!valid(u)) return(-Inf)
    st <- state(u$phi, u$s2)
    # plausibility bound pme >= -1 (signal variance at most twice Var(y))
    if (u$r < -st$P0[1, 1] / 2) return(-Inf)
    mod <- list(phi = u$phi, theta = numeric(0), Delta = numeric(0),
                Z = c(1, rep(0, p - 1)), a = rep(0, p), P = st$P0, T = st$T,
                V = st$V, h = u$r, Pn = st$P0)
    n <- length(y)
    kl <- tryCatch(suppressWarnings(stats::KalmanLike(y - u$mu, mod, nit = 0L, update = FALSE)),
                   error = function(e) NULL)
    if (is.null(kl) || !is.finite(kl$Lik) || !is.finite(kl$s2) || kl$s2 <= 0) return(-Inf)
    sumlog <- n * (2 * kl$Lik - log(kl$s2))
    -0.5 * (n * log(2 * pi) + sumlog + kl$s2 * n)
  }
  acov_x <- function(u, n) {      # autocovariances of the signal, lags 0..n-1
    st <- state(u$phi, u$s2)
    g <- numeric(n); M <- st$P0
    for (h in seq_len(n)) { g[h] <- M[1, 1]; M <- st$T %*% M }
    g
  }
  derived <- function(theta) {
    u <- unpack(theta)
    sx2 <- if (valid(u)) state(u$phi, u$s2)$P0[1, 1] else NA_real_
    g0 <- sx2 + u$r
    c(sx2 = sx2, gamma0 = g0, lambda = sx2 / g0, pme = u$r / g0)
  }
  starts <- function(y) {
    v <- stats::var(y)
    fy <- tryCatch(stats::ar.yw(y, aic = FALSE, order.max = p, demean = TRUE),
                   error = function(e) NULL)
    a  <- if (!is.null(fy)) as.numeric(fy$ar) else rep(0.1, p)
    vp <- if (!is.null(fy)) fy$var.pred else 0.9 * v
    lapply(c(0, 0.2, 0.5), function(f) {
      ph <- a / (1 - f)
      if (any(Mod(polyroot(c(1, -ph))) <= 1.01)) ph <- a
      th <- c(mean(y), ph, max(1e-4 * v, (1 - f) * vp), f * v)
      names(th) <- nms
      th
    })
  }
  simulate <- function(theta, n) {
    u <- unpack(theta)
    g <- acov_x(u, n); g[1] <- g[1] + u$r
    L <- chol(stats::toeplitz(g))
    u$mu + drop(crossprod(L, stats::rnorm(n)))
  }
  id_check <- function(fit) {
    if (p != 1) return(NULL)
    se <- sqrt(fit$vcov["phi", "phi"])
    tsme_id_message(fit$coefficients[["phi"]], fit$coefficients[["phi"]] / se)
  }
  list(names = nms, k = 1L, p = p, loglik = loglik, derived = derived, starts = starts,
       simulate = simulate, id_check = id_check,
       me_names = c("r", "sx2", "gamma0", "lambda", "pme"),
       pop_default = "pme = 0",
       defaults = tsme_defaults_univariate,
       label = sprintf("AR(%d) model with measurement error (ME), fitted with arme()", p))
}


# Specification: VAR(1) + ME -----------------------------------------------------------------
arme_spec_var <- function(k) {
  idx <- seq_len(k)
  phin <- paste0("Phi", rep(idx, each = k), rep(idx, times = k))   # Phi_ij: row i, column j
  lt <- which(lower.tri(diag(k), diag = TRUE), arr.ind = TRUE)
  lt <- lt[order(lt[, 2], lt[, 1]), , drop = FALSE]
  sign <- paste0("Sigma", lt[, 1], lt[, 2])
  nms <- c(paste0("mu", idx), phin, sign, paste0("r", idx))
  np <- c(k, k * k, nrow(lt), k)
  unpack <- function(theta) {
    theta <- unname(theta)
    mu <- theta[seq_len(k)]
    Phi <- matrix(theta[k + seq_len(k * k)], k, k, byrow = TRUE)
    Q <- matrix(0, k, k)
    Q[lt] <- theta[k + k * k + seq_len(nrow(lt))]
    Q[lt[, 2:1, drop = FALSE]] <- theta[k + k * k + seq_len(nrow(lt))]
    r <- theta[sum(np[1:3]) + seq_len(k)]
    list(mu = mu, Phi = Phi, Q = Q, r = r)
  }
  lyap <- function(Phi, Q) matrix(solve(diag(k^2) - kronecker(Phi, Phi), c(Q)), k)
  valid <- function(u) {
    all(is.finite(c(u$Phi, u$Q, u$r))) &&
      all(Mod(eigen(u$Phi, only.values = TRUE)$values) < 1 - 1e-8) &&
      !is.null(tryCatch(chol(u$Q), error = function(e) NULL))
  }
  loglik <- function(theta, Y) {
    u <- unpack(theta)
    if (!valid(u)) return(-Inf)
    n <- nrow(Y)
    Yc <- sweep(Y, 2, u$mu)
    P <- lyap(u$Phi, u$Q)
    # plausibility bound pme_i >= -1 (signal variance at most twice Var(y_i))
    if (any(u$r < -diag(P) / 2)) return(-Inf)
    Rm <- diag(u$r, k)
    a <- rep(0, k); ll <- 0; c0 <- k * log(2 * pi)
    tPhi <- t(u$Phi)
    steady <- FALSE
    for (t in seq_len(n)) {
      if (!steady) {
        cF <- tryCatch(chol(P + Rm), error = function(e) NULL)
        if (is.null(cF)) return(-Inf)
        Finv <- chol2inv(cF)
        logdet <- 2 * sum(log(diag(cF)))
        K <- P %*% Finv
        Pn <- u$Phi %*% (P - K %*% P) %*% tPhi + u$Q
        if (max(abs(Pn - P)) < 1e-12 * (1 + max(abs(P)))) steady <- TRUE
        P <- Pn
      }
      v <- Yc[t, ] - a
      ll <- ll - 0.5 * (c0 + logdet + sum(v * drop(Finv %*% v)))
      a <- drop(u$Phi %*% (a + drop(K %*% v)))
    }
    ll
  }
  derived <- function(theta) {
    u <- unpack(theta)
    sx2 <- if (valid(u)) diag(lyap(u$Phi, u$Q)) else rep(NA_real_, k)
    g0 <- sx2 + u$r
    out <- c(sx2, g0, sx2 / g0, u$r / g0)
    names(out) <- c(paste0("sx2_", idx), paste0("gamma0_", idx),
                    paste0("lambda", idx), paste0("pme", idx))
    out
  }
  starts <- function(Y) {
    n <- nrow(Y)
    Yc <- sweep(Y, 2, colMeans(Y))
    X <- Yc[-n, , drop = FALSE]; Z <- Yc[-1, , drop = FALSE]
    Phi0 <- t(solve(crossprod(X), crossprod(X, Z)))
    ev <- max(Mod(eigen(Phi0, only.values = TRUE)$values))
    if (ev >= 0.95) Phi0 <- Phi0 * 0.9 / ev
    res <- Z - X %*% t(Phi0)
    Q0 <- stats::cov(res)
    v <- diag(stats::cov(Y))
    lapply(c(0, 0.2), function(f) {
      Phi <- Phi0
      if (f > 0) {
        Phi <- Phi0 / (1 - f)
        ev <- max(Mod(eigen(Phi, only.values = TRUE)$values))
        if (ev >= 0.95) Phi <- Phi * 0.9 / ev
      }
      Q <- Q0 * (1 - f)
      th <- c(colMeans(Y), as.vector(t(Phi)), Q[lt], f * v)
      names(th) <- nms
      th
    })
  }
  simulate <- function(theta, n) {
    u <- unpack(theta)
    Sx <- lyap(u$Phi, u$Q)
    S <- matrix(0, n * k, n * k)
    G <- Sx                                  # Cov(x_i, x_j) = Sx Phi'^(j - i), j >= i
    for (h in 0:(n - 1)) {
      for (i in seq_len(n - h)) {
        ri <- (i - 1) * k + idx; rj <- (i + h - 1) * k + idx
        S[ri, rj] <- G
        if (h > 0) S[rj, ri] <- t(G)
      }
      G <- G %*% t(u$Phi)
      if (max(abs(G)) < 1e-14) break
    }
    S <- S + diag(rep(u$r, n))
    L <- chol(S)
    z <- drop(crossprod(L, stats::rnorm(n * k)))
    sweep(matrix(z, n, k, byrow = TRUE), 2, -u$mu)
  }
  defaults <- function(set) tsme_defaults_var(k, set)
  list(names = nms, k = k, p = 1L, loglik = loglik, derived = derived, starts = starts,
       simulate = simulate, id_check = function(fit) NULL,
       me_names = c(paste0("r", idx), paste0("sx2_", idx), paste0("gamma0_", idx),
                     paste0("lambda", idx), paste0("pme", idx)),
       pop_default = paste0("pme", idx, " = 0"),
       defaults = defaults,
       label = sprintf("VAR(1) model for %d variables with measurement error (ME), fitted with arme()", k))
}


# Default hypothesis sets for k measurement-error variances
tsme_defaults_var <- function(k, set) {
  idx <- seq_len(k)
  negl <- function(i) sprintf("pme%d > -delta; pme%d < delta", i, i)
  mep  <- function(i) sprintf("pme%d > delta", i)
  switch(set,
    relevance = {
      pats <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), k)))
      hyps <- list()
      for (j in seq_len(nrow(pats))) {
        pr <- pats[j, ]
        nm <- if (!any(pr)) "H_negl" else paste0("H_ME", paste(idx[pr], collapse = ""))
        hyps[[nm]] <- paste(ifelse(pr, mep(idx), negl(idx)), collapse = "; ")
      }
      list(hypotheses = hyps, comparison = "unconstrained")
    },
    minimum = list(hypotheses = list(H_ME = paste(mep(idx), collapse = "; ")),
                   comparison = "complement"),
    zero = list(hypotheses = list(H0 = paste(sprintf("r%d = 0", idx), collapse = "; "),
                                  Hpos = paste(sprintf("r%d > 0", idx), collapse = "; ")),
                comparison = "unconstrained"),
    sign = list(hypotheses = list(Hpos = paste(sprintf("r%d > 0", idx), collapse = "; ")),
                comparison = "complement")
  )
}


# Model interface for goric_tsme()
arme_interface <- function(object) {
  spec <- object$spec
  data <- object$data
  map <- function(theta) {
    names(theta) <- spec$names
    c(theta, spec$derived(theta))
  }
  n <- object$nobs
  list(theta_hat = object$coefficients, V_base = object$vcov, map = map,
       loglik = function(theta) spec$loglik(theta, data),
       starts = spec$starts(data),
       simulate = function(theta) spec$simulate(theta, n),
       refit = function(d) arme(d, p = object$p, start = object$coefficients),
       refit_estimates = function(f) map(f$coefficients),
       nobs = n, me = TRUE, me_names = spec$me_names,
       me_models = "arme() fits", me_pop_default = spec$pop_default,
       label = spec$label,
       vcov_source = "inverse of the (numerical) Hessian of the log-likelihood",
       defaults = spec$defaults,
       id_check = function() spec$id_check(object))
}


goric.arme <- function(object, ..., hypotheses = NULL,
                       comparison = NULL,
                       type = "gorica",
                       default_hypotheses = c("relevance", "minimum", "zero", "sign"),
                       delta = 0.1,
                       sensitivity = FALSE,
                       delta_values = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50),
                       VCOV_method = c("delta", "bootstrap"),
                       B = 1000L, seed = NULL, boot_trim = 0.025,
                       sample_nobs = NULL,
                       debug = FALSE) {
  default_hypotheses <- match.arg(default_hypotheses)
  VCOV_method <- match.arg(VCOV_method)
  if (any(!is.finite(object$vcov)) && VCOV_method == "delta") {
    stop(paste0("\nrestriktor ERROR: the covariance matrix of the estimates is not available; ",
                "consider VCOV_method = 'bootstrap' (the penalties of the GORIC also need it)."), call. = FALSE)
  }
  mi <- arme_interface(object)
  goric_tsme(mi, fit = object, ..., hypotheses = hypotheses, comparison = comparison,
             type = type, default_hypotheses = default_hypotheses, delta = delta,
             sensitivity = sensitivity, delta_values = delta_values,
             VCOV_method = VCOV_method, B = B, seed = seed, boot_trim = boot_trim,
             debug = debug, class_name = "goric_arme")
}


# Methods for arme objects -------------------------------------------------------------------------
coef.arme   <- function(object, ...) object$coefficients
vcov.arme   <- function(object, ...) object$vcov
logLik.arme <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients), nobs = object$nobs,
            class = "logLik")
}
nobs.arme <- function(object, ...) object$nobs

print.arme <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("\n", x$spec$label, "\n\n", sep = "")
  se <- sqrt(pmax(diag(x$vcov), 0))
  tab <- cbind(estimate = x$coefficients, se = se)
  print(round(tab, digits))
  cat("\nDerived: \n")
  print(round(x$derived, digits))
  cat("  (sx2 = variance of the signal; gamma0 = variance of y; lambda = reliability;\n",
      "  pme = proportion of ME variance = r / gamma0)\n", sep = "")
  cat("\nlog-likelihood:", format(round(x$loglik, 3)), " n =", x$nobs, "\n")
  invisible(x)
}
