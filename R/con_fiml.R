## Author: Leonard Vanbrabant and Claude
## Full information maximum likelihood (FIML) for linear models with
## missing data. Self-contained implementation (no norm/norm2/lavaan needed):
##  - EM algorithm for the saturated multivariate normal model,
##  - observed-data log-likelihood and analytic score,
##  - observed information via finite differences of the score,
##  - delta method to obtain the ACOV of the regression coefficients,
##  - ECM algorithm for *order-restricted* FIML estimation (restriktor).
##
## The joint vector z = (y, x1, ..., xk, aux1, ..., auxm) is assumed
## multivariate normal and the missingness MAR (Rubin, 1976). The linear
## model y = b0 + b'x + e is a just-identified reparameterization of the
## saturated model, so the unrestricted FIML estimates follow directly from
## the saturated EM estimates. Auxiliary variables are included in the
## saturated part and improve the missing-data handling without changing
## the target model (saturated-correlates approach; Graham, 2003).


# data preparation --------------------------------------------------------

# recover the original data (including missing values) from a fitted lm
# object and return the numeric data matrix z = [y, X (no intercept), aux]
fiml_lm_data <- function(object, auxiliary = c()) {
  if (!is.null(weights(object))) {
    stop("restriktor ERROR: missing = \"fiml\" is not available for weighted regression.",
         call. = FALSE)
  }
  if (!is.null(object$offset)) {
    stop("restriktor ERROR: missing = \"fiml\" is not available for models with an offset.",
         call. = FALSE)
  }
  if (attr(terms(object), "intercept") != 1L) {
    stop("restriktor ERROR: missing = \"fiml\" requires a model with an intercept.",
         call. = FALSE)
  }

  form <- formula(object)
  env  <- environment(form)
  if (is.null(env)) {
    env <- parent.frame()
  }
  datacall <- getCall(object)$data

  if (!is.null(datacall)) {
    data <- try(eval(datacall, envir = env), silent = TRUE)
    if (inherits(data, "try-error") || is.null(data)) {
      stop(paste("restriktor ERROR: the original data could not be recovered from the model call.",
                 "Please make sure the data set is available in the environment in which",
                 "the model was fitted."), call. = FALSE)
    }
    data <- as.data.frame(data)
  } else {
    # variables live in the formula environment
    vnames <- unique(c(all.vars(terms(object)), auxiliary))
    data <- try(as.data.frame(mget(vnames, envir = as.environment(env))), silent = TRUE)
    if (inherits(data, "try-error")) {
      stop(paste("restriktor ERROR: the original data could not be recovered.",
                 "Please fit the model using the 'data' argument."), call. = FALSE)
    }
  }

  # apply a subset argument the same way lm() did
  subsetcall <- getCall(object)$subset
  if (!is.null(subsetcall)) {
    subs <- try(eval(subsetcall, envir = data, enclos = env), silent = TRUE)
    if (inherits(subs, "try-error")) {
      stop(paste("restriktor ERROR: the 'subset' argument of the fitted model",
                 "could not be evaluated on the recovered data."), call. = FALSE)
    }
    if (is.logical(subs)) {
      subs <- which(subs)
    }
    data <- data[subs, , drop = FALSE]
  }

  # convert character variables to factors, as lm() does internally
  for (j in seq_along(data)) {
    if (is.character(data[[j]])) {
      data[[j]] <- factor(data[[j]])
    }
  }

  # model frame with missing values preserved
  mf <- model.frame(delete.response(terms(object)), data = data, na.action = na.pass)
  y  <- eval(form[[2L]], envir = data, enclos = env)
  if (!is.numeric(y) || !is.null(dim(y))) {
    stop("restriktor ERROR: missing = \"fiml\" requires a single numeric response variable.",
         call. = FALSE)
  }

  # dummy coding for factors: build the model matrix row-wise so that rows
  # with missing values are preserved (model.matrix drops or breaks on them)
  X <- fiml_model_matrix(delete.response(terms(object)), mf)
  # remove intercept column
  int_idx <- which(colnames(X) == "(Intercept)")
  X <- X[, -int_idx, drop = FALSE]

  # check that the coefficient names match the fitted model
  cnames <- names(coef(object))[-1L]
  if (!identical(as.character(colnames(X)), as.character(cnames))) {
    stop("restriktor ERROR: the FIML model matrix could not be matched to the fitted model.",
         call. = FALSE)
  }

  # auxiliary variables
  A <- NULL
  if (length(auxiliary) > 0L) {
    missing_aux <- setdiff(auxiliary, colnames(data))
    if (length(missing_aux) > 0L) {
      stop("restriktor ERROR: auxiliary variable(s) not found in the data: ",
           paste(sQuote(missing_aux), collapse = ", "), call. = FALSE)
    }
    A <- data[, auxiliary, drop = FALSE]
    is_num <- vapply(A, is.numeric, logical(1))
    if (!all(is_num)) {
      stop(paste("restriktor ERROR: auxiliary variables must be numeric.",
                 "Please recode factor variables into numeric (dummy) variables:",
                 paste(sQuote(auxiliary[!is_num]), collapse = ", ")), call. = FALSE)
    }
    A <- as.matrix(A)
    dup <- intersect(colnames(A), c(deparse(form[[2L]]), colnames(X)))
    if (length(dup) > 0L) {
      stop("restriktor ERROR: auxiliary variable(s) also appear in the model: ",
           paste(sQuote(dup), collapse = ", "), call. = FALSE)
    }
  }

  Z <- cbind(y, X, A)
  colnames(Z)[1L] <- deparse(form[[2L]])

  # remove rows for which all variables are missing
  all_miss <- rowSums(!is.na(Z)) == 0L
  if (any(all_miss)) {
    Z <- Z[!all_miss, , drop = FALSE]
  }
  if (anyNA(Z[, 1L]) && all(is.na(Z[, 1L]))) {
    stop("restriktor ERROR: the response variable is completely missing.", call. = FALSE)
  }

  list(Z = Z, nx = ncol(X), naux = if (is.null(A)) 0L else ncol(A),
       coef_names = names(coef(object)))
}


# build a model matrix in which rows with missing values are kept as NA.
# for rows with complete predictor information the result equals
# model.matrix(); incomplete rows get NA in every column that involves a
# missing variable.
fiml_model_matrix <- function(trms, mf) {
  n <- nrow(mf)
  # complete rows: standard model.matrix
  cc <- complete.cases(mf)
  if (!any(cc)) {
    stop(paste("restriktor ERROR: no complete cases available to construct the model matrix.",
               "FIML for models with factor variables requires at least some complete cases."),
         call. = FALSE)
  }
  mf_cc <- mf[cc, , drop = FALSE]
  # set factor contrasts on the complete part
  X_cc <- model.matrix(trms, mf_cc)
  X <- matrix(NA_real_, nrow = n, ncol = ncol(X_cc),
              dimnames = list(rownames(mf), colnames(X_cc)))
  X[cc, ] <- X_cc

  # incomplete rows: impute a dummy value per variable, build the model
  # matrix, and afterwards set every column that involves a missing
  # variable to NA
  ic <- which(!cc)
  if (length(ic) > 0L) {
    mf_ic <- mf[ic, , drop = FALSE]
    filled <- mf_ic
    for (j in seq_len(ncol(filled))) {
      col <- filled[[j]]
      nas <- is.na(col)
      if (!any(nas)) next
      if (is.factor(col)) {
        col[nas] <- levels(col)[1L]
      } else {
        col[nas] <- 0
      }
      filled[[j]] <- col
    }
    # ensure identical factor levels/contrasts as the complete part
    for (j in seq_len(ncol(filled))) {
      if (is.factor(filled[[j]])) {
        filled[[j]] <- factor(filled[[j]], levels = levels(mf[[j]]),
                              ordered = is.ordered(mf[[j]]))
      }
    }
    X_ic <- model.matrix(trms, filled)
    # map each model-matrix column to the model-frame variables it involves
    fac <- attr(trms, "factors")
    asgn <- attr(X_ic, "assign")
    for (k in seq_len(ncol(X_ic))) {
      if (asgn[k] == 0L) next # intercept
      vars_k <- rownames(fac)[fac[, asgn[k]] > 0L]
      vars_k <- intersect(vars_k, colnames(mf_ic))
      if (length(vars_k) == 0L) next
      miss_k <- rowSums(is.na(mf_ic[, vars_k, drop = FALSE])) > 0L
      X_ic[miss_k, k] <- NA_real_
    }
    X[ic, ] <- X_ic
  }

  X
}


# missing data patterns ---------------------------------------------------

fiml_patterns <- function(Z) {
  obs <- !is.na(Z)
  pat_id <- apply(obs, 1L, function(r) paste(as.integer(r), collapse = ""))
  idx_list <- split(seq_len(nrow(Z)), pat_id)

  pats <- lapply(idx_list, function(idx) {
    o <- which(obs[idx[1L], ])
    list(o    = o,
         n    = length(idx),
         Z    = Z[idx, o, drop = FALSE],
         idx  = idx)
  })
  unname(pats)
}


# observed-data log-likelihood for the saturated MVN model ----------------

fiml_sat_loglik <- function(mu, Sigma, pats) {
  ll <- 0
  for (pat in pats) {
    o <- pat$o
    So <- Sigma[o, o, drop = FALSE]
    cholSo <- tryCatch(chol(So), error = function(e) NULL)
    if (is.null(cholSo)) return(-Inf)
    logdet <- 2 * sum(log(diag(cholSo)))
    cent <- sweep(pat$Z, 2L, mu[o])
    q <- backsolve(cholSo, t(cent), transpose = TRUE)
    ll <- ll - 0.5 * (pat$n * (length(o) * log(2 * pi) + logdet) + sum(q * q))
  }
  ll
}


# analytic score of the observed-data log-likelihood wrt
# theta = (mu, vech(Sigma)) (lower-triangular column-major order)
fiml_sat_score <- function(mu, Sigma, pats) {
  p <- length(mu)
  g_mu <- numeric(p)
  G <- matrix(0, p, p)
  for (pat in pats) {
    o <- pat$o
    Soinv <- chol2inv(chol(Sigma[o, o, drop = FALSE]))
    cent <- sweep(pat$Z, 2L, mu[o])
    g_mu[o] <- g_mu[o] + Soinv %*% colSums(cent)
    Sc <- crossprod(cent)
    G[o, o] <- G[o, o] + 0.5 * (Soinv %*% Sc %*% Soinv - pat$n * Soinv)
  }
  # gradient wrt the unique elements of Sigma: off-diagonal elements appear
  # twice in the symmetric matrix
  g_sigma <- numeric(p * (p + 1) / 2)
  pos <- 1L
  for (j in seq_len(p)) {
    for (i in j:p) {
      g_sigma[pos] <- if (i == j) G[i, i] else 2 * G[i, j]
      pos <- pos + 1L
    }
  }
  c(g_mu, g_sigma)
}


# vech helpers (lower-triangular, column-major)
fiml_vech <- function(S) S[lower.tri(S, diag = TRUE)]
fiml_vech_inv <- function(v, p) {
  S <- matrix(0, p, p)
  S[lower.tri(S, diag = TRUE)] <- v
  S[upper.tri(S)] <- t(S)[upper.tri(S)]
  S
}


# E-step: expected complete-data sufficient statistics T1 = E(sum z) and
# T2 = E(sum z z') given the observed data and (mu, Sigma)
fiml_estep <- function(mu, Sigma, pats, p) {
  T1 <- numeric(p)
  T2 <- matrix(0, p, p)
  for (pat in pats) {
    o <- pat$o
    m <- setdiff(seq_len(p), o)
    Zo <- pat$Z
    if (length(m) == 0L) {
      T1 <- T1 + colSums(Zo)
      T2 <- T2 + crossprod(Zo)
    } else {
      Soinv <- chol2inv(chol(Sigma[o, o, drop = FALSE]))
      B <- Sigma[m, o, drop = FALSE] %*% Soinv
      Cm <- Sigma[m, m, drop = FALSE] - B %*% Sigma[o, m, drop = FALSE]
      Em <- matrix(mu[m], nrow = pat$n, ncol = length(m), byrow = TRUE) +
        sweep(Zo, 2L, mu[o]) %*% t(B)
      T1[o] <- T1[o] + colSums(Zo)
      T1[m] <- T1[m] + colSums(Em)
      T2[o, o] <- T2[o, o] + crossprod(Zo)
      T2om <- crossprod(Zo, Em)
      T2[o, m] <- T2[o, m] + T2om
      T2[m, o] <- T2[m, o] + t(T2om)
      T2[m, m] <- T2[m, m] + crossprod(Em) + pat$n * Cm
    }
  }
  T2 <- (T2 + t(T2)) / 2
  list(T1 = T1, T2 = T2)
}


# EM algorithm for the saturated MVN model --------------------------------

fiml_em_sat <- function(Z, max.iter = 5000L, tol = 1e-09, verbose = FALSE) {
  p <- ncol(Z)
  n <- nrow(Z)
  pats <- fiml_patterns(Z)

  # starting values: available-case means and covariances
  mu <- colMeans(Z, na.rm = TRUE)
  Sigma <- suppressWarnings(cov(Z, use = "pairwise.complete.obs"))
  Sigma[is.na(Sigma)] <- 0
  diag_fallback <- apply(Z, 2L, var, na.rm = TRUE)
  diag_fallback[is.na(diag_fallback) | diag_fallback <= 0] <- 1
  diag(Sigma)[is.na(diag(Sigma)) | diag(Sigma) <= 0] <-
    diag_fallback[is.na(diag(Sigma)) | diag(Sigma) <= 0]
  # bend to positive definiteness if needed
  ev <- eigen(Sigma, symmetric = TRUE)
  if (min(ev$values) < 1e-06 * max(ev$values)) {
    ev$values <- pmax(ev$values, 1e-06 * max(ev$values))
    Sigma <- ev$vectors %*% (ev$values * t(ev$vectors))
  }

  converged <- FALSE
  iter <- 0L
  theta_old <- c(mu, fiml_vech(Sigma))
  for (iter in seq_len(max.iter)) {
    ss <- fiml_estep(mu, Sigma, pats, p)
    mu <- ss$T1 / n
    Sigma <- ss$T2 / n - tcrossprod(mu)
    Sigma <- (Sigma + t(Sigma)) / 2
    theta_new <- c(mu, fiml_vech(Sigma))
    if (verbose) {
      cat("EM iteration", iter, ": loglik =", fiml_sat_loglik(mu, Sigma, pats), "\n")
    }
    # convergence on the (relative) parameter change; a log-likelihood
    # criterion stops too early because EM converges linearly
    if (max(abs(theta_new - theta_old) / (abs(theta_old) + 1)) < tol) {
      converged <- TRUE
      break
    }
    theta_old <- theta_new
  }
  ll_old <- fiml_sat_loglik(mu, Sigma, pats)
  if (!converged) {
    warning("restriktor WARNING: the EM algorithm for the saturated model did not converge.",
            " You may increase the number of iterations via control = list(max.iter = ...).",
            call. = FALSE)
  }

  names(mu) <- colnames(Z)
  dimnames(Sigma) <- list(colnames(Z), colnames(Z))

  list(mu = mu, Sigma = Sigma, loglik = ll_old, n = n, pats = pats,
       iter = iter, converged = converged)
}


# observed information matrix of theta = (mu, vech(Sigma)) via central
# finite differences of the analytic score
fiml_sat_information <- function(mu, Sigma, pats) {
  p <- length(mu)
  theta <- c(mu, fiml_vech(Sigma))
  q <- length(theta)
  H <- matrix(0, q, q)
  for (k in seq_len(q)) {
    h <- 1e-05 * (abs(theta[k]) + 1e-02)
    th_p <- theta; th_p[k] <- th_p[k] + h
    th_m <- theta; th_m[k] <- th_m[k] - h
    s_p <- fiml_sat_score(th_p[1:p], fiml_vech_inv(th_p[-(1:p)], p), pats)
    s_m <- fiml_sat_score(th_m[1:p], fiml_vech_inv(th_m[-(1:p)], p), pats)
    H[, k] <- (s_p - s_m) / (2 * h)
  }
  H <- (H + t(H)) / 2
  -H
}


# mapping from theta = (mu, vech(Sigma)) to the regression coefficients
# (intercept, slopes) of y (= first variable) on x (variables 2:(nx+1)).
# auxiliary variables (if any) occupy the remaining positions and only
# affect the estimates through the saturated (mu, Sigma).
fiml_lm_coef_map <- function(theta, p, nx) {
  mu <- theta[1:p]
  Sigma <- fiml_vech_inv(theta[-(1:p)], p)
  if (nx == 0L) {
    return(mu[1L])
  }
  xi <- 2:(nx + 1L)
  beta <- solve(Sigma[xi, xi, drop = FALSE], Sigma[xi, 1L])
  b0 <- mu[1L] - sum(beta * mu[xi])
  c(b0, beta)
}


# residual variance implied by theta
fiml_lm_sigma2 <- function(mu, Sigma, nx) {
  if (nx == 0L) {
    return(Sigma[1L, 1L])
  }
  xi <- 2:(nx + 1L)
  beta <- solve(Sigma[xi, xi, drop = FALSE], Sigma[xi, 1L])
  as.numeric(Sigma[1L, 1L] - sum(beta * Sigma[xi, 1L]))
}


# unrestricted FIML fit for a linear model --------------------------------

con_fiml_lm <- function(object, auxiliary = c(), control = list()) {

  max.iter <- if (is.null(control$max.iter)) 1000L else control$max.iter
  # backwards compatibility with the former emControl argument
  if (!is.null(control$iter.max)) max.iter <- control$iter.max
  tol <- if (is.null(control$tol)) 1e-09 else control$tol
  if (!is.null(control$criterion)) tol <- control$criterion
  verbose <- isTRUE(control$verbose)

  prep <- fiml_lm_data(object, auxiliary = auxiliary)
  Z  <- prep$Z
  p  <- ncol(Z)
  nx <- prep$nx

  # stage 1: saturated EM
  em <- fiml_em_sat(Z, max.iter = max.iter, tol = tol, verbose = verbose)

  # FIML estimates of the regression coefficients (the linear model is a
  # just-identified reparameterization of the saturated model)
  theta <- c(em$mu, fiml_vech(em$Sigma))
  est <- fiml_lm_coef_map(theta, p, nx)
  names(est) <- prep$coef_names
  sigma2 <- fiml_lm_sigma2(em$mu, em$Sigma, nx)

  # observed information and ACOV of theta
  info <- fiml_sat_information(em$mu, em$Sigma, em$pats)
  acov_theta <- tryCatch(chol2inv(chol(info)), error = function(e) {
    MASS::ginv(info)
  })

  # delta method: numeric Jacobian of the coefficient mapping
  q <- length(theta)
  J <- matrix(0, length(est), q)
  for (k in seq_len(q)) {
    h <- 1e-06 * (abs(theta[k]) + 1e-02)
    th_p <- theta; th_p[k] <- th_p[k] + h
    th_m <- theta; th_m[k] <- th_m[k] - h
    J[, k] <- (fiml_lm_coef_map(th_p, p, nx) - fiml_lm_coef_map(th_m, p, nx)) / (2 * h)
  }
  VCOV <- J %*% acov_theta %*% t(J)
  VCOV <- (VCOV + t(VCOV)) / 2
  dimnames(VCOV) <- list(prep$coef_names, prep$coef_names)

  list(est = est, VCOV = VCOV, sigma2 = sigma2, loglik = em$loglik,
       N = em$n, mu = em$mu, Sigma = em$Sigma, pats = em$pats, Z = Z,
       nx = nx, naux = prep$naux, iter = em$iter, converged = em$converged,
       control = list(max.iter = max.iter, tol = tol))
}


# order-restricted FIML via ECM -------------------------------------------
#
# the complete-data log-likelihood of z = (y, x, a) factorizes into
#   f(x; mu_x, Sigma_xx) * f(y | x; b, s2) * f(a | y, x; c0, C, Psi),
# with variation-independent parameter blocks. Hence EM with three exact
# conditional maximization steps (ECM) is a genuine EM algorithm and the
# observed-data log-likelihood increases monotonically. The inequality
# constraints Amat %*% b >= bvec (first meq rows equalities) only affect
# the f(y | x) block, whose M-step is a quadratic program.

fiml_ecm_restricted <- function(fiml, Amat, bvec, meq,
                                max.iter = 5000L, tol = 1e-09, verbose = FALSE) {
  Z <- fiml$Z
  p <- ncol(Z)
  n <- nrow(Z)
  nx <- fiml$nx
  naux <- fiml$naux
  pats <- fiml$pats

  xi <- if (nx > 0L) 2:(nx + 1L) else integer(0)
  wi <- 1:(nx + 1L)                       # (y, x)
  ai <- if (naux > 0L) (nx + 2L):p else integer(0)

  # starting values: unrestricted FIML solution
  mu <- fiml$mu
  Sigma <- fiml$Sigma

  build_implied <- function(b, s2, mux, Sxx, c0, C, Psi) {
    mu_new <- numeric(p)
    Sigma_new <- matrix(0, p, p)
    # (y, x) block
    if (nx > 0L) {
      beta <- b[-1L]
      mu_new[xi] <- mux
      mu_new[1L] <- b[1L] + sum(beta * mux)
      Sigma_new[xi, xi] <- Sxx
      Syx <- as.numeric(Sxx %*% beta)
      Sigma_new[1L, xi] <- Syx
      Sigma_new[xi, 1L] <- Syx
      Sigma_new[1L, 1L] <- sum(beta * Syx) + s2
    } else {
      mu_new[1L] <- b[1L]
      Sigma_new[1L, 1L] <- s2
    }
    # auxiliary block: a = c0 + C w + e,  e ~ N(0, Psi)
    if (naux > 0L) {
      Sww <- Sigma_new[wi, wi, drop = FALSE]
      mu_new[ai] <- c0 + as.numeric(C %*% mu_new[wi])
      Saw <- C %*% Sww
      Sigma_new[ai, wi] <- Saw
      Sigma_new[wi, ai] <- t(Saw)
      Sigma_new[ai, ai] <- Saw %*% t(C) + Psi
    }
    Sigma_new <- (Sigma_new + t(Sigma_new)) / 2
    list(mu = mu_new, Sigma = Sigma_new)
  }

  converged <- FALSE
  b <- NULL; s2 <- NULL
  theta_old <- c(mu, fiml_vech(Sigma))

  for (iter in seq_len(max.iter)) {
    # E-step
    ss <- fiml_estep(mu, Sigma, pats, p)
    T1 <- ss$T1
    T2 <- ss$T2

    # CM-step 1: marginal x-distribution
    if (nx > 0L) {
      mux <- T1[xi] / n
      Sxx <- T2[xi, xi, drop = FALSE] / n - tcrossprod(mux)
      Sxx <- (Sxx + t(Sxx)) / 2
    } else {
      mux <- numeric(0)
      Sxx <- matrix(0, 0L, 0L)
    }

    # CM-step 2: restricted regression of y on (1, x) via quadratic
    # programming; expected complete-data SSE is quadratic in b
    M <- rbind(cbind(n, t(T1[xi])),
               cbind(T1[xi], T2[xi, xi, drop = FALSE]))
    d <- c(T1[1L], T2[xi, 1L])
    if (all(Amat == 0L)) {
      b <- solve(M, d)
    } else {
      qp <- quadprog::solve.QP(Dmat = M, dvec = d, Amat = t(Amat),
                               bvec = bvec, meq = meq)
      b <- qp$solution
    }
    s2 <- as.numeric((T2[1L, 1L] - 2 * sum(b * d) + t(b) %*% M %*% b) / n)
    if (s2 < .Machine$double.eps) {
      stop("restriktor ERROR: the restricted FIML residual variance became zero.",
           call. = FALSE)
    }

    # CM-step 3: saturated regression of the auxiliary variables on (1, y, x)
    c0 <- NULL; C <- NULL; Psi <- NULL
    if (naux > 0L) {
      Mw <- rbind(cbind(n, t(T1[wi])),
                  cbind(T1[wi], T2[wi, wi, drop = FALSE]))
      Dw <- rbind(t(T1[ai]), T2[wi, ai, drop = FALSE])
      G <- solve(Mw, Dw)                  # (1 + nw) x naux
      c0 <- G[1L, ]
      C <- t(G[-1L, , drop = FALSE])      # naux x nw
      Psi <- (T2[ai, ai, drop = FALSE] - t(Dw) %*% G) / n
      Psi <- (Psi + t(Psi)) / 2
    }

    imp <- build_implied(b, s2, mux, Sxx, c0, C, Psi)
    mu <- imp$mu
    Sigma <- imp$Sigma

    if (verbose) {
      cat("ECM iteration", iter, ": loglik =", fiml_sat_loglik(mu, Sigma, pats), "\n")
    }
    theta_new <- c(mu, fiml_vech(Sigma))
    if (max(abs(theta_new - theta_old) / (abs(theta_old) + 1)) < tol) {
      converged <- TRUE
      break
    }
    theta_old <- theta_new
  }
  ll_old <- fiml_sat_loglik(mu, Sigma, pats)
  if (!converged) {
    warning("restriktor WARNING: the restricted FIML (ECM) algorithm did not converge.",
            " You may increase the number of iterations via control = list(max.iter = ...).",
            call. = FALSE)
  }

  names(b) <- names(fiml$est)

  list(b.restr = b, sigma2 = s2, loglik = ll_old, mu = mu, Sigma = Sigma,
       iter = iter, converged = converged)
}


# model-implied R-squared: 1 - s2 / var(y)
fiml_R2 <- function(b, s2, Sigma, nx) {
  if (nx == 0L) {
    return(0)
  }
  xi <- 2:(nx + 1L)
  beta <- b[-1L]
  expl <- as.numeric(t(beta) %*% Sigma[xi, xi, drop = FALSE] %*% beta)
  expl / (expl + s2)
}
