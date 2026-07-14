## ═══════════════════════════════════════════════════════════════════════════
## HOIJ_bootstrapD_experiment.R  (v1)
##
## Doel: kwantificeren hoe goed de tweede-orde infinitesimale jackknife
## (HOIJ-2, Giordano et al. 2019) de refit-lus van bootstrapD() kan
## vervangen bij order-constrained toetsing in SEM (restriktor-setting),
## voor (1) de ENKELE bootstrap en (2) de DUBBELE bootstrap (inner-only).
##
## Kernidee (zie discussie): de D-statistiek stroomafwaarts van
## (theta.hat.boot, I.hat.boot) is al een tweede-orde-object (twee QP's
## op de kwadratische benadering van de log-likelihood). HOIJ vervangt
## uitsluitend de dure stap: de ML-refit per trekking. De QP's — en
## daarmee de niet-gladde active-set-logica — blijven per trekking EXACT.
##
## Vergelijkingsprotocol: per trekking worden identieke bootstrapindices
## aan de exacte route (refit) en de HOIJ-route (expansie) gevoerd, zodat
## D_exact en D_hoij trekking-voor-trekking vergelijkbaar zijn. Alle
## Monte-Carlo-variatie valt daarmee uit de vergelijking weg; wat
## overblijft is zuivere approximatiefout.
##
## Ontwerpbeslissingen (changelog-conventie; nieuw bestand, prefix [D]):
## [D1] HOIJ-kernfuncties (compute_loglik_casewise, compute_all_J,
##      make_grad_F, calibrate_alpha, compute_T_tensor_grad) zijn
##      LETTERLIJK overgenomen uit HOIJ_benchmarkstudie.R (v5) zodat de
##      schaalconventies (Scores = -s_i/N, info.observed = J_som/N,
##      T = d3(-loglik/N)) identiek geverifieerd blijven.
## [D2] Informatiematrix per trekking: de exacte route gebruikt
##      lavTech(fit, "information.observed") (unit-schaal). NB: dit wijkt
##      af van bootstrapD(), dat lavInspect(h1, "information") (expected)
##      gebruikt. Reden: de HOIJ-surrogaat benadert de GEOBSERVEERDE
##      informatie ((1/N) sum_i w_i J_i); zou de exacte route expected
##      gebruiken, dan vermengt de vergelijking HOIJ-fout met het
##      observed-vs-expected-verschil, dat niets met HOIJ te maken heeft.
##      De D-statistiek zelf is onder beide keuzes gedefinieerd; voor de
##      approximatievraag is observed-vs-observed de zuivere vergelijking.
## [D3] HOIJ-informatiesurrogaat: I_w(theta0) = H_obs + J_dw, met
##      J_dw = (Delta_w %*% J_i)/N (exact: herweging van casewise
##      krommingen), optioneel + T·dtheta als eerste-orde correctie voor
##      de evaluatie in theta_hoij i.p.v. theta0 (INFO_T_CORRECTIE). De
##      T-correctie gebruikt de full-data T als benadering van de
##      herwogen T; de fout daarvan is O(N^-1/2)·O(dtheta) = hogere orde.
## [D4] Demping als in benchmark [C2]: ||d2|| begrensd op
##      KAPPA_DAMP*||d1||; ongedempte variant (hoij2_raw) en zuiver
##      eerste-orde (ij) lopen mee als referentie, conform [C3].
## [D5] Geen stille fallbacks (conform [C15]): mislukte prep, QP of
##      convergentie geeft NA + teller, nooit heimelijk een andere
##      methode.
## [D6] Type-B "unconstrained fit": de nul-restrictie-QP uit bootstrapD
##      heeft als oplossing exact theta zelf; de waarde wordt hier
##      analytisch berekend (-0.5*theta'I theta), semantisch identiek.
## [D7] Constraints worden 1x geparset (lavaan:::lav_constraints_parse)
##      op de originele h1-fit; dit veronderstelt LINEAIRE restricties
##      (constante Jacobiaan), zoals in de order-constrained setting
##      standaard is. Voor niet-lineaire restricties zou de JAC per
##      trekking herberekend moeten worden — buiten scope.
## [D8] Bollen-Stine wordt gereproduceerd als datatransformatie + gewone
##      multinomiale resampling van rijen (identiek aan bootstrapD).
##      Daardoor valt de BS-variant exact binnen het gewichtsvector-
##      raamwerk van Giordano et al.: de HOIJ-expansie wordt verankerd in
##      EEN extra fit van h1 op de getransformeerde volledige data.
##      De parametrische variant valt buiten het raamwerk (verse
##      mvrnorm-data, geen herweging) en wordt hier niet meegenomen.
## [D9] Dubbele bootstrap: alleen de INNER lus wordt door HOIJ vervangen
##      (de verdedigbare stap). De outer lus blijft exact; per outer
##      trekking is 1 extra ankerfit op de her-getransformeerde data
##      nodig. Kosten: van ~R_OUT*(2 + R_IN) fits naar ~R_OUT*3 fits.
##      D.original van de inner lus komt in BEIDE routes uit dezelfde
##      exacte fit.boot.h1, zodat het plugin-p-verschil zuiver de
##      inner-lus-approximatie meet.
## [D10] Niet-PD informatiematrices (solve.QP vereist PD) krijgen een
##      minimale ridge (eigenwaarde-verschuiving); toepassing wordt in
##      BEIDE routes geteld en gerapporteerd (ridge_exact/ridge_hoij).
## [D11] Rekenkosten worden primair in AANTAL MODELFITS gerapporteerd
##      (algoritmische complexiteit), wandkloktijd secundair; de
##      afgeleide-kosten van de prep (casewise J's, T-tensor) worden
##      apart geteld als gradient-/loglik-evaluaties.
## [D12] FDB (fast double bootstrap) is met 1 inner trekking per outer
##      trekking al goedkoop; de winst zit vrijwel volledig bij
##      double.bootstrap = "standard". FDB is hier niet geimplementeerd.
## [D13] Experiment 3 (grens-nabij): per dataset wordt de ruis (x, e_m,
##      e_y) gefixeerd en de effectgrootte e (a = b = e) via uniroot
##      gekalibreerd zodat D_orig(e) op de asymptotische chi-bar-kwadraat-
##      kritieke waarde bij alpha ligt. De mengselgewichten (1/4, 1/2,
##      1/4) veronderstellen orthogonale restricties en dienen UITSLUITEND
##      als kalibratiedoel; de gerealiseerde p_exact (die rond alpha
##      scattert) wordt gerapporteerd en draagt de vergelijking. De
##      beslissingsvraag: hoe vaak wijkt verwerp(p_hoij < alpha) af van
##      verwerp(p_exact < alpha) op identieke trekkingen — d.w.z. hoe
##      vaak duwt zuivere approximatiefout de p over de grens.
## [D14] Experiment 3 gebruikt uitsluitend bollen.stine (D* is daar een
##      nulverdeling, dus p_exact ~ alpha is daar de toets-relevante
##      situatie) en R_BOOT3 = 999 (gangbare toetsresolutie; het
##      discretisatie-effect op p is dan 1/999).
## [D15] BUGFIX na eerste pilotrun: bs_transform kreeg de datamatrix in
##      data.frame-kolomvolgorde (x,m,y) terwijl icov en Sigma.hat in
##      lavaans interne variabelenvolgorde (m,y,x) staan; de transformatie
##      legde H0 daardoor NIET op. Signatuur: geinfleerde D*-verdeling bij
##      effect-DGP's (exp1 zwak/BS: q95 ~ 38-77 i.p.v. ~4-8) en p_exact
##      ver van alpha in exp3 ondanks geslaagde D_orig-kalibratie; bij
##      nul-DGP's ~ correct omdat de transform daar ~ de identiteit is.
##      Fix: X komt nu uit h0_fit@Data@X (interne volgorde) en de functie
##      verifieert met stopifnot dat cov(Xt) = Sigma.hat(h0). NB: de
##      exact-vs-HOIJ-vergelijkingen uit de eerste run blijven intern
##      geldig (beide routes zagen identieke data en indices); de
##      procedureniveau-p's onder bollen.stine uit die run zijn dat niet.
##
## Output: ./hoij_bootstrapD_output/  (RDS met per-trekking-vergelijking
##         + samenvattende tabellen op de console)
## Geverifieerd op: R 4.3+, lavaan 0.6-17, quadprog 1.5-8.
## ═══════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(lavaan)
  library(quadprog)
  library(MASS)
})

## ─────────────────────────────────────────────────────────────────────────────
## 0. CONFIG
## ─────────────────────────────────────────────────────────────────────────────

SEED_BASIS <- 20260706

## Experiment 1: enkele bootstrap (ordinary + bollen.stine)
RUN_EXP1    <- TRUE
SCENARIOS   <- c("nul", "zwak")   # DGP onder H0 en onder zwak alternatief
K_DATASETS  <- 3                  # datasets per scenario (pilot: 3 | serieus: >= 50)
N_OBS       <- 200                # steekproefomvang
R_BOOT      <- 400                # bootstraptrekkingen (pilot: 400 | serieus: >= 1000)

## Experiment 2: dubbele bootstrap "standard", inner-only HOIJ [D9]
RUN_EXP2    <- TRUE
K_DB        <- 2                  # datasets (pilot: 2 | serieus: >= 20)
R_OUTER     <- 30                 # outer trekkingen (pilot: 30 | serieus: >= 500)
R_INNER     <- 100                # inner trekkingen (pilot: 100 | serieus: >= 500)
SCENARIO_DB <- "zwak"

## Experiment 3: grens-nabij (beslissingsovereenstemming bij alpha) [D13/D14]
RUN_EXP3    <- TRUE
K3          <- 12                 # datasets (pilot: 12 | serieus: >= 100)
R_BOOT3     <- 999                # trekkingen (toetsresolutie 1/999)

## HOIJ-instellingen (conventies benchmark)
KAPPA_DAMP       <- 0.5           # [C2]/[D4]
INFO_T_CORRECTIE <- TRUE          # [D3]
ALPHA_SPREAD_TOL <- 0.1           # [C10]
ALPHA_NIVEAU     <- 0.05          # nominaal niveau (voor adj.alpha in exp2)
TYPE_TOETS       <- "A"           # "A": H0 gelijkheden vs H1 ongelijkheden

out_dir <- file.path(getwd(), "hoij_bootstrapD_output")
if (!dir.exists(out_dir)) dir.create(out_dir)

## ─────────────────────────────────────────────────────────────────────────────
## 1. MODEL, CONSTRAINTS EN DGP
##    Klassieke one-sided mediationtoets: H0: a = b = 0 vs H1: a > 0, b > 0
## ─────────────────────────────────────────────────────────────────────────────

model_h1 <- '
  m ~ a*x
  y ~ b*m + c*x
'
model_h0 <- '
  m ~ 0*x
  y ~ 0*m + c*x
'
constraints_txt <- '
a > 0
b > 0
'

## DGP-parameters per scenario
dgp_pars <- list(
  nul  = c(a = 0.0, b = 0.0, cp = 0.2),
  zwak = c(a = 0.3, b = 0.3, cp = 0.2)
)

gen_data <- function(N, scen) {
  p  <- dgp_pars[[scen]]
  x  <- rnorm(N)
  m  <- p["a"] * x + rnorm(N)
  y  <- p["b"] * m + p["cp"] * x + rnorm(N)
  data.frame(x = x, m = m, y = y)
}

## Fit-helper met globale fitteller [D11]
FIT_TELLER <- new.env(); FIT_TELLER$n <- 0L
fit_sem <- function(model, data) {
  FIT_TELLER$n <- FIT_TELLER$n + 1L
  fit <- tryCatch(
    sem(model, data = data, fixed.x = FALSE, estimator = "ML",
        se = "none", test = "none", baseline = FALSE, h1 = FALSE,
        warn = FALSE),
    error = function(e) NULL)
  if (is.null(fit) || !lavInspect(fit, "converged")) return(NULL)
  fit
}
reset_teller <- function() { n <- FIT_TELLER$n; FIT_TELLER$n <- 0L; n }

## ─────────────────────────────────────────────────────────────────────────────
## 2. HOIJ-KERNFUNCTIES — letterlijk uit HOIJ_benchmarkstudie.R [D1]
## ─────────────────────────────────────────────────────────────────────────────

compute_loglik_casewise <- function(fit, theta) {
  X <- fit@Data@X[[1]]
  N <- nrow(X); p <- ncol(X)
  GLIST <- lavaan:::lav_model_x2GLIST(fit@Model, x = theta)
  implied <- lavaan:::lav_model_implied(fit@Model, GLIST = GLIST)
  Sigma <- implied$cov[[1]]
  mu <- implied$mean[[1]]
  if (is.null(mu) || length(mu) == 0) mu <- colMeans(X)
  Sigma_inv <- tryCatch(solve(Sigma), error = function(e) MASS::ginv(Sigma))
  log_det <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  const <- -0.5 * (p * log(2 * pi) + log_det)
  X_centered <- sweep(X, 2, mu, "-")
  quad_form <- rowSums((X_centered %*% Sigma_inv) * X_centered)
  as.numeric(const - 0.5 * quad_form)
}

compute_all_J <- function(fit, theta0, delta = 1e-5) {
  D <- length(theta0); N <- nrow(fit@Data@X[[1]])
  J_array <- array(0, dim = c(N, D, D))
  for (k in 1:D) {
    for (l in k:D) {
      if (k == l) {
        tp <- theta0; tp[k] <- tp[k] + delta
        tm <- theta0; tm[k] <- tm[k] - delta
        ll_p <- compute_loglik_casewise(fit, tp)
        ll_0 <- compute_loglik_casewise(fit, theta0)
        ll_m <- compute_loglik_casewise(fit, tm)
        J_array[, k, k] <- -((ll_p - 2*ll_0 + ll_m) / (delta^2))
      } else {
        tpp <- theta0; tpp[k] <- tpp[k]+delta; tpp[l] <- tpp[l]+delta
        tpm <- theta0; tpm[k] <- tpm[k]+delta; tpm[l] <- tpm[l]-delta
        tmp_ <- theta0; tmp_[k] <- tmp_[k]-delta; tmp_[l] <- tmp_[l]+delta
        tmm <- theta0; tmm[k] <- tmm[k]-delta; tmm[l] <- tmm[l]-delta
        J_array[, k, l] <- -((compute_loglik_casewise(fit, tpp) -
                                compute_loglik_casewise(fit, tpm) -
                                compute_loglik_casewise(fit, tmp_) +
                                compute_loglik_casewise(fit, tmm)) / (4*delta^2))
        J_array[, l, k] <- J_array[, k, l]
      }
    }
  }
  J_array
}

make_grad_F <- function(fit) {
  lavmodel       <- fit@Model
  lavsamplestats <- fit@SampleStats
  lavdata        <- fit@Data
  lavcache       <- fit@Cache
  function(theta) {
    GLIST <- lavaan:::lav_model_x2GLIST(lavmodel, x = theta)
    as.numeric(lavaan:::lav_model_gradient(
      lavmodel       = lavmodel,
      GLIST          = GLIST,
      lavsamplestats = lavsamplestats,
      lavdata        = lavdata,
      lavcache       = lavcache))
  }
}

calibrate_alpha <- function(grad_F, theta0, H_observed, h = 1e-5) {
  D_loc <- length(theta0)
  H_grad <- matrix(NA_real_, D_loc, D_loc)
  for (k in 1:D_loc) {
    tp <- theta0; tp[k] <- tp[k] + h
    tm <- theta0; tm[k] <- tm[k] - h
    H_grad[, k] <- (grad_F(tp) - grad_F(tm)) / (2 * h)
  }
  H_grad <- (H_grad + t(H_grad)) / 2
  idx    <- abs(H_grad) > 1e-6 * max(abs(H_grad))
  ratio  <- as.numeric(H_observed)[idx] / as.numeric(H_grad)[idx]
  alpha  <- median(ratio)
  spread <- max(abs(ratio / alpha - 1))
  list(alpha = alpha, spread = spread)
}

compute_T_tensor_grad <- function(grad_F, theta, alpha, h = 1e-4) {
  D_loc <- length(theta)
  T_arr <- array(0, dim = c(D_loc, D_loc, D_loc))
  g0    <- grad_F(theta)
  for (l in 1:D_loc) {
    for (m in l:D_loc) {
      if (l == m) {
        tp <- theta; tp[l] <- tp[l] + h
        tm <- theta; tm[l] <- tm[l] - h
        col <- (grad_F(tp) - 2 * g0 + grad_F(tm)) / h^2
      } else {
        tpp <- theta; tpp[l] <- tpp[l] + h; tpp[m] <- tpp[m] + h
        tpm <- theta; tpm[l] <- tpm[l] + h; tpm[m] <- tpm[m] - h
        tmp_ <- theta; tmp_[l] <- tmp_[l] - h; tmp_[m] <- tmp_[m] + h
        tmm <- theta; tmm[l] <- tmm[l] - h; tmm[m] <- tmm[m] - h
        col <- (grad_F(tpp) - grad_F(tpm) - grad_F(tmp_) + grad_F(tmm)) / (4*h^2)
      }
      T_arr[, l, m] <- alpha * col
      T_arr[, m, l] <- alpha * col
    }
  }
  T_arr <- (T_arr +
              aperm(T_arr, c(2, 1, 3)) + aperm(T_arr, c(3, 2, 1)) +
              aperm(T_arr, c(1, 3, 2)) + aperm(T_arr, c(2, 3, 1)) +
              aperm(T_arr, c(3, 1, 2))) / 6
  T_arr
}

## ─────────────────────────────────────────────────────────────────────────────
## 3. HOIJ-PREP EN PER-TREKKING-EVALUATIE
## ─────────────────────────────────────────────────────────────────────────────

## Alle trekking-onafhankelijke grootheden, 1x per ankerfit.
## Retourneert NULL met attribuut "reden" bij falen [D5].
hoij_prepare <- function(fit) {
  faal <- function(r) { out <- NULL; attr(out, "reden") <- r; out }
  theta0 <- tryCatch(coef(fit, type = "free"), error = function(e) NULL)
  if (is.null(theta0)) return(faal("coef"))
  D <- length(theta0); N <- nrow(fit@Data@X[[1]])
  Scores <- tryCatch(lavScores(fit, scaling = TRUE), error = function(e) NULL)
  if (is.null(Scores)) return(faal("lavScores"))
  H_obs <- tryCatch(lavTech(fit, "information.observed"), error = function(e) NULL)
  if (is.null(H_obs)) return(faal("H_obs"))
  H.inv <- tryCatch(solve(H_obs), error = function(e) NULL)
  if (is.null(H.inv)) return(faal("H_inv"))
  grad_F <- make_grad_F(fit)
  cal <- tryCatch(calibrate_alpha(grad_F, theta0, H_obs), error = function(e) NULL)
  if (is.null(cal) || !is.finite(cal$alpha) || cal$spread > ALPHA_SPREAD_TOL)
    return(faal("alpha_kalibratie"))
  T_arr <- tryCatch(compute_T_tensor_grad(grad_F, theta0, cal$alpha),
                    error = function(e) NULL)
  if (is.null(T_arr)) return(faal("T_tensor"))
  J_all <- tryCatch(compute_all_J(fit, theta0, delta = 1e-5),
                    error = function(e) NULL)
  if (is.null(J_all)) return(faal("J_casewise"))
  list(theta0 = theta0, D = D, N = N,
       Scores = Scores, H_obs = H_obs, H.inv = H.inv,
       Tmat  = matrix(T_arr, nrow = D),          # D x D^2 (voor K-hat-term)
       Tmat2 = matrix(T_arr, nrow = D * D),      # D^2 x D (voor info-correctie [D3])
       J_all_2d = matrix(J_all, nrow = N, ncol = D * D))
}

## PD-waarborg voor solve.QP [D10]
as_pd <- function(M, eps = 1e-8) {
  M <- (M + t(M)) / 2
  ev_min <- min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
  if (is.finite(ev_min) && ev_min > eps) list(M = M, ridged = FALSE)
  else list(M = M + diag(eps - ev_min + 1e-12, nrow(M)), ridged = TRUE)
}

## D-statistiek: gedeeld door exacte en HOIJ-route, zodat de QP-laag
## per definitie identiek is. Amat: rijen = restricties, gelijkheden eerst.
D_stat <- function(theta, I_unit, N, Amat, rhs, n_eq, type = TYPE_TOETS) {
  pd <- as_pd(I_unit)
  I_unit <- pd$M
  dvec <- drop(I_unit %*% theta)
  out.ineq <- tryCatch(
    solve.QP(Dmat = I_unit, dvec = dvec, Amat = t(Amat), bvec = rhs, meq = n_eq),
    error = function(e) NULL)
  if (is.null(out.ineq)) return(list(D = NA_real_, ridged = pd$ridged))
  if (type == "A") {
    out.eq <- tryCatch(
      solve.QP(Dmat = I_unit, dvec = dvec, Amat = t(Amat), bvec = rhs,
               meq = nrow(Amat)),
      error = function(e) NULL)
    if (is.null(out.eq)) return(list(D = NA_real_, ridged = pd$ridged))
    D <- N * 2 * (out.eq$value - out.ineq$value)
  } else {
    val.un <- -0.5 * drop(t(theta) %*% I_unit %*% theta)   # [D6]
    D <- N * 2 * (out.ineq$value - val.un)
  }
  list(D = D, ridged = pd$ridged)
}

## HOIJ-evaluatie van D voor een matrix bootstrapindices IDX (R x N).
## Retourneert per trekking D voor ij / hoij2 / hoij2_raw + diagnostiek.
hoij_D_draws <- function(prep, IDX, Amat, rhs, n_eq) {
  R <- nrow(IDX); N <- prep$N; D <- prep$D
  W <- t(apply(IDX, 1, function(idx) tabulate(idx, nbins = N) - 1L))
  C_mat <- (W %*% prep$Scores) %*% prep$H.inv       # rijen: c_vec per trekking
  JW_2d <- (W %*% prep$J_all_2d) / N                # rijen: J_dw per trekking
  HT    <- prep$H.inv %*% prep$Tmat
  res <- data.frame(D_ij = rep(NA_real_, R), D_hoij2 = NA_real_,
                    D_hoij2_raw = NA_real_, s_damp = NA_real_,
                    ridge = FALSE)
  for (i in seq_len(R)) {
    c_vec  <- C_mat[i, ]
    J_dw_i <- matrix(JW_2d[i, ], D, D)
    Bc     <- drop(prep$H.inv %*% J_dw_i %*% c_vec)
    kron_cc <- as.vector(tcrossprod(c_vec))
    Ac     <- 0.5 * drop(HT %*% kron_cc)
    d1     <- -c_vec
    d2     <- Bc - Ac
    n1 <- sqrt(sum(d1^2)); n2 <- sqrt(sum(d2^2))
    s  <- if (n2 > 0) min(1, KAPPA_DAMP * n1 / n2) else 1
    res$s_damp[i] <- s
    th_ij   <- prep$theta0 + d1
    th_h    <- prep$theta0 + d1 + s * d2
    th_raw  <- prep$theta0 + d1 + d2
    ## Informatiesurrogaat [D3]: I_w(theta0) = H_obs + J_dw (herweging
    ## van casewise krommingen), optioneel + T-correctie naar theta_hoij.
    I_w0 <- prep$H_obs + J_dw_i
    info_at <- function(th) {
      if (!INFO_T_CORRECTIE) return(I_w0)
      dth <- th - prep$theta0
      I_w0 + matrix(prep$Tmat2 %*% dth, D, D)
    }
    o_ij  <- D_stat(th_ij,  I_w0,          N, Amat, rhs, n_eq)
    o_h   <- D_stat(th_h,   info_at(th_h), N, Amat, rhs, n_eq)
    o_raw <- D_stat(th_raw, info_at(th_raw), N, Amat, rhs, n_eq)
    res$D_ij[i]        <- o_ij$D
    res$D_hoij2[i]     <- o_h$D
    res$D_hoij2_raw[i] <- o_raw$D
    res$ridge[i]       <- o_ij$ridged || o_h$ridged || o_raw$ridged
  }
  res
}

## Exacte D uit een lavaan-fit (gedeelde D_stat-laag) [D2]
exact_D_from_fit <- function(fit, Amat, rhs, n_eq) {
  th <- coef(fit, type = "free")
  I  <- tryCatch(lavTech(fit, "information.observed"), error = function(e) NULL)
  if (is.null(I)) return(list(D = NA_real_, ridged = NA))
  D_stat(th, I, N = nrow(fit@Data@X[[1]]), Amat, rhs, n_eq)
}

## ─────────────────────────────────────────────────────────────────────────────
## 4. BOLLEN-STINE-TRANSFORMATIE [D8] + CONSTRAINTS [D7]
## ─────────────────────────────────────────────────────────────────────────────

## [D15] X komt uit h0_fit@Data@X (lavaans INTERNE variabelenvolgorde,
## consistent met icov en Sigma.hat); de stopifnot verifieert dat de
## getransformeerde data de H0-geimpliceerde covariantie exact
## reproduceren — de check die de volgorde-bug direct had gevangen.
bs_transform <- function(h0_fit) {
  X <- h0_fit@Data@X[[1]]
  implied   <- lavaan:::lav_model_implied(h0_fit@Model)
  Sigma.hat <- implied$cov[[1]]
  sigma.sqrt <- lav_matrix_symmetric_sqrt(Sigma.hat)
  S.inv.sqrt <- lav_matrix_symmetric_sqrt(h0_fit@SampleStats@icov[[1]])
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Xt <- Xc %*% S.inv.sqrt %*% sigma.sqrt
  if (h0_fit@Model@meanstructure) {
    Mu.hat <- implied$mean[[1]]
    Xt <- scale(Xt, center = -1 * Mu.hat, scale = FALSE)
  }
  N <- nrow(Xt)
  cov_N <- crossprod(scale(Xt, center = TRUE, scale = FALSE)) / N
  stopifnot("bs_transform: cov(Xt) reproduceert Sigma.hat(h0) niet" =
              max(abs(cov_N - Sigma.hat)) < 1e-8 * max(1, max(abs(Sigma.hat))))
  colnames(Xt) <- lavNames(h0_fit, "ov")
  Xt
}

parse_constraints <- function(h1_fit, constraints) {
  pt <- as.list(parTable(h1_fit))
  conInfo <- lavaan:::lav_constraints_parse(pt, constraints = constraints,
                                            theta = coef(h1_fit))
  Amat <- rbind(conInfo$ceq.JAC, conInfo$cin.JAC)
  rhs  <- c(conInfo$ceq.rhs, conInfo$cin.rhs)
  list(Amat = Amat, rhs = rhs, n_eq = NROW(conInfo$ceq.JAC))
}

## ─────────────────────────────────────────────────────────────────────────────
## 5. PRE-FLIGHT ZELFTEST op het mediationanker (conventie [C11])
## ─────────────────────────────────────────────────────────────────────────────

selftest_mediatie <- function() {
  cat("── Zelftest op mediationmodel (lavaan ",
      as.character(packageVersion("lavaan")), ") ──\n", sep = "")
  set.seed(SEED_BASIS)
  dat <- gen_data(300, "zwak")
  fit <- fit_sem(model_h1, dat)
  stopifnot("ankerfit convergeert niet" = !is.null(fit))
  th0 <- coef(fit, type = "free"); D <- length(th0); N <- nrow(dat)
  h <- 1e-6
  S_num <- sapply(1:D, function(k) {
    tp <- th0; tp[k] <- tp[k] + h; tm <- th0; tm[k] <- tm[k] - h
    (compute_loglik_casewise(fit, tp) - compute_loglik_casewise(fit, tm)) / (2*h)
  })
  Sc <- lavScores(fit, scaling = TRUE)
  ## NB: anders dan de Holzinger-zelftest in de benchmark heeft dit
  ## padmodel veel (numeriek) nul-elementen; ratio's alleen nemen waar
  ## de noemer boven de ruisdrempel ligt (zelfde filter als
  ## calibrate_alpha), anders vervuilt 0/0-ruis de mediaan.
  ix_s <- abs(S_num) > 1e-6 * max(abs(S_num))
  r_sc <- median(as.numeric(Sc)[ix_s] / as.numeric(S_num)[ix_s]) * N
  cat(sprintf("  (a) lavScores-schaal : N*ratio = %+.6f (verwacht -1)\n", r_sc))
  stopifnot("lavScores-schaal wijkt af" = abs(r_sc + 1) < 0.01)
  J_sum <- apply(compute_all_J(fit, th0), c(2, 3), sum)
  H_obs <- lavTech(fit, "information.observed")
  ix_H <- abs(J_sum) > 1e-6 * max(abs(J_sum))
  r_H <- median(as.numeric(H_obs)[ix_H] / as.numeric(J_sum)[ix_H]) * N
  cat(sprintf("  (b) info.observed    : N*ratio = %+.6f (verwacht +1)\n", r_H))
  stopifnot("informatie-identiteit wijkt af" = abs(r_H - 1) < 0.01)
  ## (c) 1-trekking-consistentie: gewichten = 1_N moeten D exact reproduceren
  prep <- hoij_prepare(fit)
  stopifnot("hoij_prepare faalt in zelftest" = !is.null(prep))
  con <- parse_constraints(fit, constraints_txt)
  IDX1 <- matrix(rep(1:N, 1), nrow = 1)              # elke case 1x: Delta_w = 0
  r1 <- hoij_D_draws(prep, IDX1, con$Amat, con$rhs, con$n_eq)
  D0 <- exact_D_from_fit(fit, con$Amat, con$rhs, con$n_eq)$D
  cat(sprintf("  (c) w = 1_N          : D_hoij = %.8f, D_exact = %.8f\n",
              r1$D_hoij2, D0))
  stopifnot("HOIJ reproduceert D niet bij w = 1_N" =
              abs(r1$D_hoij2 - D0) < 1e-6 * max(1, abs(D0)))
  ## (d) [D15] BS-transformatie legt H0 op: a en b moeten 0 zijn op de
  ## getransformeerde data (de check die de volgorde-bug had gevangen)
  h0_st <- fit_sem(model_h0, dat)
  stopifnot("h0-fit in zelftest convergeert niet" = !is.null(h0_st))
  f_tr <- fit_sem(model_h1, as.data.frame(bs_transform(h0_st)))
  ab_tr <- coef(f_tr)[c("a", "b")]
  cat(sprintf("  (d) BS legt H0 op    : |a|,|b| na transform = %.2e, %.2e\n",
              abs(ab_tr[1]), abs(ab_tr[2])))
  stopifnot("BS-transformatie legt H0 niet op" = all(abs(ab_tr) < 1e-6))
  cat("  Zelftest GESLAAGD.\n\n")
  invisible(TRUE)
}

## ─────────────────────────────────────────────────────────────────────────────
## 6. EXPERIMENT 1 — ENKELE BOOTSTRAP (exact vs HOIJ, identieke indices)
## ─────────────────────────────────────────────────────────────────────────────

run_exp1_dataset <- function(scen, k, boot_type) {
  set.seed(SEED_BASIS + 1000 * k + 17 * nchar(scen) + 3 * (boot_type == "bollen.stine"))
  dat <- gen_data(N_OBS, scen)
  N <- nrow(dat)
  
  ## Originele fits (zoals bootstrapD): h1 op ruwe data, h0 alleen voor BS
  reset_teller()
  h1 <- fit_sem(model_h1, dat)
  if (is.null(h1)) return(NULL)
  con <- parse_constraints(h1, constraints_txt)
  D_orig <- exact_D_from_fit(h1, con$Amat, con$rhs, con$n_eq)$D
  
  ## Databron voor de bootstrap: ruwe data (ordinary) of BS-getransformeerd
  if (boot_type == "bollen.stine") {
    h0 <- fit_sem(model_h0, dat)
    if (is.null(h0)) return(NULL)
    Xb <- bs_transform(h0)                          # [D15]
  } else {
    Xb <- as.matrix(dat)
  }
  fits_setup <- reset_teller()
  
  ## HOIJ-anker: 1 fit van h1 op de bootstrapdatabron [D8]
  reset_teller()
  anker <- fit_sem(model_h1, as.data.frame(Xb))
  if (is.null(anker)) return(NULL)
  prep <- hoij_prepare(anker)
  fits_hoij <- reset_teller()
  if (is.null(prep)) {
    cat(sprintf("    [%s/k=%d/%s] hoij_prepare FAALT: %s\n",
                scen, k, boot_type, attr(prep, "reden")))
    return(NULL)
  }
  
  ## Gedeelde bootstrapindices
  IDX <- t(replicate(R_BOOT, sample.int(N, N, replace = TRUE)))
  
  ## Exacte route: refit per trekking
  reset_teller()
  D_ex <- rep(NA_real_, R_BOOT)
  for (b in seq_len(R_BOOT)) {
    fb <- fit_sem(model_h1, as.data.frame(Xb[IDX[b, ], , drop = FALSE]))
    if (!is.null(fb)) D_ex[b] <- exact_D_from_fit(fb, con$Amat, con$rhs, con$n_eq)$D
  }
  fits_exact <- reset_teller()
  
  ## HOIJ-route: expansie per trekking (0 fits)
  hd <- hoij_D_draws(prep, IDX, con$Amat, con$rhs, con$n_eq)
  
  ok <- is.finite(D_ex) & is.finite(hd$D_hoij2)
  pv <- function(Dv) mean(Dv[is.finite(Dv)] > D_orig)
  list(
    scen = scen, k = k, boot_type = boot_type, D_orig = D_orig,
    per_draw = data.frame(D_exact = D_ex, hd),
    samenvatting = data.frame(
      scen = scen, k = k, boot = boot_type,
      n_ok_exact  = sum(is.finite(D_ex)),
      n_ok_hoij   = sum(is.finite(hd$D_hoij2)),
      cor_hoij    = if (sum(ok) > 2) cor(D_ex[ok], hd$D_hoij2[ok]) else NA,
      medae_hoij  = median(abs(D_ex[ok] - hd$D_hoij2[ok])),
      maxae_hoij  = max(abs(D_ex[ok] - hd$D_hoij2[ok])),
      medae_ij    = median(abs(D_ex - hd$D_ij), na.rm = TRUE),
      medae_raw   = median(abs(D_ex - hd$D_hoij2_raw), na.rm = TRUE),
      q95_exact   = quantile(D_ex, 0.95, na.rm = TRUE),
      q95_hoij    = quantile(hd$D_hoij2, 0.95, na.rm = TRUE),
      p_exact     = pv(D_ex),
      p_hoij      = pv(hd$D_hoij2),
      p_ij        = pv(hd$D_ij),
      frac_gedempt = mean(hd$s_damp < 1, na.rm = TRUE),
      ridge       = mean(hd$ridge),
      fits_exact  = fits_exact,
      fits_hoij   = fits_hoij + fits_setup   # anker (+ h0 bij BS) + h1-origineel
    )
  )
}

## ─────────────────────────────────────────────────────────────────────────────
## 7. EXPERIMENT 2 — DUBBELE BOOTSTRAP "standard", INNER-ONLY HOIJ [D9]
## ─────────────────────────────────────────────────────────────────────────────

run_exp2_dataset <- function(k) {
  set.seed(SEED_BASIS + 90000 + 1000 * k)
  dat <- gen_data(N_OBS, SCENARIO_DB)
  N <- nrow(dat)
  
  h1 <- fit_sem(model_h1, dat); h0 <- fit_sem(model_h0, dat)
  if (is.null(h1) || is.null(h0)) return(NULL)
  con <- parse_constraints(h1, constraints_txt)
  D_orig <- exact_D_from_fit(h1, con$Amat, con$rhs, con$n_eq)$D
  Xb <- bs_transform(h0)                         # outer BS-databron [D15]
  reset_teller()
  
  D_outer     <- rep(NA_real_, R_OUTER)
  plugin_ex   <- rep(NA_real_, R_OUTER)
  plugin_hoij <- rep(NA_real_, R_OUTER)
  fits_ex_pad <- 0L; fits_h_pad <- 0L
  prep_faal   <- 0L
  
  for (r in seq_len(R_OUTER)) {
    idx_out <- sample.int(N, N, replace = TRUE)
    X_out <- Xb[idx_out, , drop = FALSE]
    d_out <- as.data.frame(X_out)
    
    ## Outer fits — EXACT in beide routes (gedeeld)
    reset_teller()
    fb_h1 <- fit_sem(model_h1, d_out)
    fb_h0 <- fit_sem(model_h0, d_out)
    f_shared <- reset_teller()
    fits_ex_pad <- fits_ex_pad + f_shared
    fits_h_pad  <- fits_h_pad  + f_shared
    if (is.null(fb_h1) || is.null(fb_h0)) next
    
    D_outer[r] <- exact_D_from_fit(fb_h1, con$Amat, con$rhs, con$n_eq)$D
    D_orig_in  <- D_outer[r]      # inner-D.original: uit fit.boot.h1, gedeeld [D9]
    
    ## Inner BS-databron: hertransformatie t.o.v. fit.boot.h0
    X_in <- bs_transform(fb_h0)                    # [D15]
    IDX_in <- t(replicate(R_INNER, sample.int(N, N, replace = TRUE)))
    
    ## (i) EXACTE inner lus
    reset_teller()
    D_in_ex <- rep(NA_real_, R_INNER)
    for (b in seq_len(R_INNER)) {
      fin <- fit_sem(model_h1, as.data.frame(X_in[IDX_in[b, ], , drop = FALSE]))
      if (!is.null(fin))
        D_in_ex[b] <- exact_D_from_fit(fin, con$Amat, con$rhs, con$n_eq)$D
    }
    fits_ex_pad <- fits_ex_pad + reset_teller()
    plugin_ex[r] <- mean(D_in_ex[is.finite(D_in_ex)] > D_orig_in)
    
    ## (ii) HOIJ inner lus: 1 ankerfit op X_in, daarna expansie
    reset_teller()
    anker <- fit_sem(model_h1, as.data.frame(X_in))
    fits_h_pad <- fits_h_pad + reset_teller()
    if (is.null(anker)) { prep_faal <- prep_faal + 1L; next }
    prep <- hoij_prepare(anker)
    if (is.null(prep)) { prep_faal <- prep_faal + 1L; next }
    hd <- hoij_D_draws(prep, IDX_in, con$Amat, con$rhs, con$n_eq)
    plugin_hoij[r] <- mean(hd$D_hoij2[is.finite(hd$D_hoij2)] > D_orig_in)
  }
  
  p_outer <- mean(D_outer[is.finite(D_outer)] > D_orig)
  adj_p <- function(pl) mean(pl[is.finite(pl)] < p_outer)
  adj_a <- function(pl) unname(quantile(pl, ALPHA_NIVEAU, na.rm = TRUE))
  ok <- is.finite(plugin_ex) & is.finite(plugin_hoij)
  
  list(
    k = k,
    per_outer = data.frame(D_outer, plugin_ex, plugin_hoij),
    samenvatting = data.frame(
      k = k, p_outer = p_outer,
      cor_plugin  = if (sum(ok) > 2) cor(plugin_ex[ok], plugin_hoij[ok]) else NA,
      mae_plugin  = mean(abs(plugin_ex[ok] - plugin_hoij[ok])),
      max_plugin  = max(abs(plugin_ex[ok] - plugin_hoij[ok])),
      adj_p_exact = adj_p(plugin_ex),  adj_p_hoij = adj_p(plugin_hoij),
      adj_a_exact = adj_a(plugin_ex),  adj_a_hoij = adj_a(plugin_hoij),
      prep_faal   = prep_faal,
      fits_exact  = fits_ex_pad, fits_hoij = fits_h_pad
    )
  )
}

## ─────────────────────────────────────────────────────────────────────────────
## 7b. EXPERIMENT 3 — GRENS-NABIJ: beslissingsovereenstemming bij alpha [D13]
## ─────────────────────────────────────────────────────────────────────────────

## Asymptotische chi-bar-kwadraat-kritieke waarde (kalibratiedoel; [D13])
chibar_krit <- function(alpha) {
  uniroot(function(c) 0.5 * (1 - pchisq(c, 1)) + 0.25 * (1 - pchisq(c, 2)) - alpha,
          c(1e-3, 30))$root
}

run_exp3_dataset <- function(k) {
  set.seed(SEED_BASIS + 500000 + 1000 * k)
  ## Ruis fixeren zodat D_orig(e) een deterministische functie van e is
  x  <- rnorm(N_OBS); em <- rnorm(N_OBS); ey <- rnorm(N_OBS)
  make_dat <- function(e) {
    m <- e * x + em
    data.frame(x = x, m = m, y = e * m + 0.2 * x + ey)
  }
  c_krit <- chibar_krit(ALPHA_NIVEAU)
  
  ## Kalibratie: D_orig(e) = c_krit via uniroot; con 1x parsen op e=0.3-fit
  reset_teller()
  fit_ref <- fit_sem(model_h1, make_dat(0.3))
  if (is.null(fit_ref)) return(NULL)
  con <- parse_constraints(fit_ref, constraints_txt)
  D_van_e <- function(e) {
    f <- fit_sem(model_h1, make_dat(e))
    if (is.null(f)) return(NA_real_)
    exact_D_from_fit(f, con$Amat, con$rhs, con$n_eq)$D
  }
  h_lo <- D_van_e(0) - c_krit
  h_hi <- D_van_e(0.8) - c_krit
  if (!is.finite(h_lo) || !is.finite(h_hi)) return(NULL)
  if (h_lo >= 0) { e_star <- 0 }                    # ruis alleen al boven c_krit
  else if (h_hi <= 0) { e_star <- 0.8 }             # bracket te smal (rapporteren)
  else e_star <- uniroot(function(e) D_van_e(e) - c_krit,
                         c(0, 0.8), tol = 1e-3)$root
  fits_kal <- reset_teller()
  
  dat <- make_dat(e_star)
  reset_teller()
  h1 <- fit_sem(model_h1, dat); h0 <- fit_sem(model_h0, dat)
  if (is.null(h1) || is.null(h0)) return(NULL)
  D_orig <- exact_D_from_fit(h1, con$Amat, con$rhs, con$n_eq)$D
  Xb <- bs_transform(h0)                            # BS: D* is nulverdeling [D14][D15]
  anker <- fit_sem(model_h1, as.data.frame(Xb))
  if (is.null(anker)) return(NULL)
  prep <- hoij_prepare(anker)
  fits_setup <- reset_teller()
  if (is.null(prep)) return(NULL)
  
  IDX <- t(replicate(R_BOOT3, sample.int(N_OBS, N_OBS, replace = TRUE)))
  reset_teller()
  D_ex <- rep(NA_real_, R_BOOT3)
  for (b in seq_len(R_BOOT3)) {
    fb <- fit_sem(model_h1, as.data.frame(Xb[IDX[b, ], , drop = FALSE]))
    if (!is.null(fb)) D_ex[b] <- exact_D_from_fit(fb, con$Amat, con$rhs, con$n_eq)$D
  }
  fits_exact <- reset_teller()
  hd <- hoij_D_draws(prep, IDX, con$Amat, con$rhs, con$n_eq)
  
  pv <- function(Dv) mean(Dv[is.finite(Dv)] > D_orig)
  p_ex <- pv(D_ex); p_h <- pv(hd$D_hoij2); p_ij <- pv(hd$D_ij)
  data.frame(
    k = k, e_kal = round(e_star, 4), D_orig = D_orig,
    p_exact = p_ex, p_hoij = p_h, p_ij = p_ij,
    verwerp_exact = p_ex < ALPHA_NIVEAU,
    verwerp_hoij  = p_h  < ALPHA_NIVEAU,
    flip = (p_ex < ALPHA_NIVEAU) != (p_h < ALPHA_NIVEAU),
    q95_exact = quantile(D_ex, 0.95, na.rm = TRUE),
    q95_hoij  = quantile(hd$D_hoij2, 0.95, na.rm = TRUE),
    frac_gedempt = mean(hd$s_damp < 1, na.rm = TRUE),
    ridge = mean(hd$ridge),
    fits_exact = fits_exact,
    fits_hoij  = fits_setup + fits_kal   # kalibratie telt in BEIDE routes; hier
    # conservatief volledig bij HOIJ geboekt
  )
}

## ─────────────────────────────────────────────────────────────────────────────
## 8. UITVOERING
## ─────────────────────────────────────────────────────────────────────────────

selftest_mediatie()

if (RUN_EXP1) {
  cat("═══ EXPERIMENT 1: enkele bootstrap (R =", R_BOOT, ", N =", N_OBS, ") ═══\n")
  res1 <- list(); i <- 0L
  t1 <- Sys.time()
  for (scen in SCENARIOS) for (bt in c("ordinary", "bollen.stine"))
    for (k in seq_len(K_DATASETS)) {
      r <- run_exp1_dataset(scen, k, bt)
      if (!is.null(r)) { i <- i + 1L; res1[[i]] <- r }
    }
  cat(sprintf("  klaar in %.1f min\n", as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  
  samen1 <- do.call(rbind, lapply(res1, `[[`, "samenvatting"))
  rownames(samen1) <- NULL
  num <- vapply(samen1, is.numeric, TRUE)
  samen1[num] <- lapply(samen1[num], function(x) round(x, 4))
  cat("\nPer dataset (D_exact vs D_hoij2 op identieke trekkingen):\n")
  print(samen1, row.names = FALSE)
  
  cat("\nGeaggregeerd per bootstraptype x scenario:\n")
  agg1 <- aggregate(cbind(cor_hoij, medae_hoij, medae_ij, p_exact, p_hoij,
                          frac_gedempt) ~ boot + scen, data = samen1, FUN = mean)
  agg1[-(1:2)] <- lapply(agg1[-(1:2)], round, 4)
  print(agg1, row.names = FALSE)
  
  saveRDS(res1, file.path(out_dir, "exp1_enkele_bootstrap.rds"))
}

if (RUN_EXP2) {
  cat("\n═══ EXPERIMENT 2: dubbele bootstrap (R_out =", R_OUTER,
      ", R_in =", R_INNER, ") ═══\n")
  res2 <- list(); i <- 0L
  t2 <- Sys.time()
  for (k in seq_len(K_DB)) {
    r <- run_exp2_dataset(k)
    if (!is.null(r)) { i <- i + 1L; res2[[i]] <- r }
  }
  cat(sprintf("  klaar in %.1f min\n", as.numeric(difftime(Sys.time(), t2, units = "mins"))))
  samen2 <- do.call(rbind, lapply(res2, `[[`, "samenvatting"))
  rownames(samen2) <- NULL
  num <- vapply(samen2, is.numeric, TRUE)
  samen2[num] <- lapply(samen2[num], function(x) round(x, 4))
  cat("\nPer dataset (plugin-p exact vs HOIJ; outer lus gedeeld):\n")
  print(samen2, row.names = FALSE)
  saveRDS(res2, file.path(out_dir, "exp2_dubbele_bootstrap.rds"))
}

if (RUN_EXP3) {
  cat("\n═══ EXPERIMENT 3: grens-nabij, beslissing bij alpha =", ALPHA_NIVEAU,
      "(R =", R_BOOT3, ", bollen.stine) ═══\n")
  cat(sprintf("  kalibratiedoel D_orig = chi-bar-kritiek = %.4f [D13]\n",
              chibar_krit(ALPHA_NIVEAU)))
  res3 <- list(); i <- 0L
  t3 <- Sys.time()
  for (k in seq_len(K3)) {
    r <- run_exp3_dataset(k)
    if (!is.null(r)) { i <- i + 1L; res3[[i]] <- r }
  }
  cat(sprintf("  klaar in %.1f min\n", as.numeric(difftime(Sys.time(), t3, units = "mins"))))
  samen3 <- do.call(rbind, res3)
  rownames(samen3) <- NULL
  num <- vapply(samen3, is.numeric, TRUE)
  samen3[num] <- lapply(samen3[num], function(x) round(x, 4))
  cat("\nPer dataset (p_exact gekalibreerd rond alpha; identieke trekkingen):\n")
  print(samen3, row.names = FALSE)
  cat(sprintf("\nBeslissingsovereenstemming: %d/%d datasets gelijk (flips: %d)\n",
              sum(!samen3$flip), nrow(samen3), sum(samen3$flip)))
  cat(sprintf("Getekend p-verschil (hoij - exact): mediaan %+.4f, bereik [%+.4f, %+.4f]\n",
              median(samen3$p_hoij - samen3$p_exact),
              min(samen3$p_hoij - samen3$p_exact),
              max(samen3$p_hoij - samen3$p_exact)))
  cat(sprintf("NB: p-resolutie = 1/%d = %.4f; flips zijn verwacht wanneer\n",
              R_BOOT3, 1 / R_BOOT3))
  cat("|p_exact - alpha| kleiner is dan de systematische staartfout (~0.01-0.02).\n")
  saveRDS(res3, file.path(out_dir, "exp3_grens_nabij.rds"))
}

cat("\nNB [D11]: fits_hoij telt alle modelfits van de HOIJ-route (anker +\n")
cat("gedeelde outer fits); de afgeleide-prep (casewise J's, T-tensor) kost\n")
cat("daarnaast O(D^2) gradient-/loglik-evaluaties per anker, geen fits.\n")
cat("Resultaten opgeslagen in:", out_dir, "\n")
