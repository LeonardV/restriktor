coef.con_goric <- function(object, ...)  {
  return(object$ormle$b.restr)
}

coef.gorica_est <- function(object, ...)  {
  # [CHANGE 2026-10 | audit] A13: coef.gorica_est() appends the defined parameters (:=)
  b <- object$b.restr
  # defined parameters (':=') are appended, as coef.restriktor() does
  if (!is.null(object$parTable$op) && any(object$parTable$op == ":=") &&
      is.function(object$CON$def.function)) {
    b <- c(b, object$CON$def.function(b))
  }
  b
  # [/CHANGE 2026-10]
}

coef_named_vector <- function(x, VCOV = NULL, ...)  {
  # TO DO wat als vcov niet bestaat of nog niet matrix (maar bijv dpo)?
  # ms zelfs suppressWarnings(vcov(...)) gebruiken dan
  if (!is.vector(coef(x))) {
    # Some fit objects (like mlm) render matrices (>1 x >1 matrices).
    # So, if not a vector, make it one, with names(!) coming from vcov.
    est <- as.vector(coef(x))
    if (is.null(VCOV)) {
      names(est) <- rownames(vcov(x))
      # TO DO dit werkt niet, er is niet alleen : maar ook _ en laatste mag ws niet (runt iig niet)
    } else { 
      # TO DO eigenlijk nog checken of rownames bestaan
      names(est) <- rownames(VCOV)
    }
    # [CHANGE 2026-10 | audit] B5: message about the mlm coefficient names moved to goric.lm() (message_mlm_coef_names())
    # the message about these names is given once in goric.lm() (see
    # message_mlm_coef_names()), irrespective of the comparison
  } else {
    est <- coef(x)
  }
  return(est)
}

# [CHANGE 2026-10 | audit] B5: new function message_mlm_coef_names(): message with the mlm coefficient names to use in hypotheses
# Message (once, from goric.lm()) about the coefficient names of an mlm object:
# the names of vcov() with ':' and '(' ')' replaced by '.', e.g., 'Age.GroupNo'
# and 'Age..Intercept.' for the intercept of response 'Age'.
message_mlm_coef_names <- function(object) {
  labs <- gsub("[:()]", ".", rownames(vcov(object)))
  message("\nrestriktor Message: The coefficient matrix of the mlm object has been ",
          "converted into a vector. The coefficient names are the row names of ",
          "the covariance matrix (vcov) with ':', '(' and ')' replaced by '.', ",
          "e.g., 'y1..Intercept.' for the intercept of response 'y1'. ",
          "Use these names when specifying hypotheses: ",
          paste(sQuote(labs), collapse = ", "), ".")
  invisible(labs)
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] M1: new function loglik_unrestr(): unrestricted log-likelihood, also for mlm
# log-likelihood of the unrestricted model. logLik.lm() does not support
# multiple responses, so for mlm objects con_loglik_lm() is used.
loglik_unrestr <- function(model.org)  {
  if (inherits(model.org, "mlm")) {
    con_loglik_lm(model.org)
  } else {
    logLik(model.org)
  }
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] M3: new function coef_vec_mlm(): mlm coefficient matrix as a named vector in vcov() order
# For mlm objects the coefficients form a p x ny matrix. Convert it to a
# vector in the order of vcov(), with the names of vcov() (e.g., 'Age:GroupNo').
# Other objects are returned as is.
coef_vec_mlm <- function(b, model.org)  {
  if (inherits(model.org, "mlm") && is.matrix(b)) {
    b <- structure(as.vector(b), names = rownames(vcov(model.org)))
  }
  b
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] E6/N2/B2: new function check_weights(): central validation of priorICweights/study_weights
# Central validation of weight vectors (priorICweights, study_weights):
# finite, non-negative, positive sum and (optionally) a given length.
# Returns the weights rescaled to sum 1 (if rescale = TRUE). Zero weights are
# allowed (unless allow_zero = FALSE); the calling code must handle them
# (e.g., on the log scale, see ic_weights_log()).
check_weights <- function(w, name = "priorICweights", length_expected = NULL,
                          what = NULL, rescale = TRUE, allow_zero = TRUE) {
  if (is.null(w)) {
    return(NULL)
  }
  # a 1 x k (or k x 1) matrix is accepted as a vector
  if (is.numeric(w) && !is.null(dim(w))) {
    w <- as.vector(w)
  }
  if (!is.null(length_expected) && length(w) != length_expected) {
    stop("\nrestriktor ERROR: The argument '", name, "' should consist of ",
         length_expected, ngettext(length_expected, " element", " elements"),
         if (!is.null(what)) paste0(", namely ", what),
         ". It now consists of ", length(w), 
         ngettext(length(w), " element.", " elements."), call. = FALSE)
  }
  if (!is.numeric(w) || anyNA(w) || any(!is.finite(w))) {
    stop("\nrestriktor ERROR: The argument '", name, "' should be a numeric vector ",
         "with finite values (no NA, NaN, Inf).", call. = FALSE)
  }
  if (any(w < 0)) {
    stop("\nrestriktor ERROR: The argument '", name, "' should only contain ",
         "non-negative values.", call. = FALSE)
  }
  if (!allow_zero && any(w == 0)) {
    stop("\nrestriktor ERROR: The argument '", name, "' should only contain ",
         "positive values.", call. = FALSE)
  }
  if (sum(w) <= 0) {
    stop("\nrestriktor ERROR: The argument '", name, "' should contain at least ",
         "one positive value.", call. = FALSE)
  }
  if (rescale && !isTRUE(all.equal(sum(w), 1))) {
    w <- w / sum(w)
  }
  w
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] E7/E11/B14: new function ic_weights_log(): IC weights via log-sum-exp (prior 0 -> weight 0; -Inf/NaN IC handled)
# Numerically stable (prior-weighted) IC weights:
#   w_i = prior_i * exp(-IC_i / 2) / sum_j prior_j * exp(-IC_j / 2)
# computed via log-sum-exp. A zero prior gives a zero weight (not NaN),
# and very large IC differences do not underflow to 0/0.
# An IC of -Inf (infinite support) gives weight 1 for that model (shared
# equally among several models with IC = -Inf) and 0 for the others; a NaN/NA
# IC gives NaN weights.
ic_weights_log <- function(IC, prior = NULL) {
  if (is.null(prior)) {
    prior <- rep(1, length(IC))
  }
  prior <- unname(prior)
  lw <- -IC / 2 + log(prior)              # log(0) = -Inf for a zero prior
  lw[prior == 0] <- -Inf
  if (anyNA(lw)) {
    return(rep(NaN, length(IC)))
  }
  if (any(lw == Inf)) {
    w <- as.numeric(lw == Inf)
    return(w / sum(w))
  }
  m <- max(lw)
  if (!is.finite(m)) {
    # all models have IC = Inf (or prior 0): no information at all
    return(rep(NaN, length(IC)))
  }
  w <- exp(lw - m)
  w / sum(w)
}
# [/CHANGE 2026-10]

check_sample_nobs <- function(sample_nobs, ...)  {
  if (length(sample_nobs) > 1) { 
    # Probably group sizes, then sample size is sum of group sizes
    sample_nobs <- sum(sample_nobs)
    message(paste0(
    "\nrestriktor Message: The argument 'sample_nobs' contains more than one value. ",
    "It is assumed that these represent group sizes; their total (", sample_nobs, 
    ") is used as the overall sample size."
    ))
  }
  return(sample_nobs)
}

check_N_with_sample_nobs <- function(N, sample_nobs, ...)  {
  # Check on N
  if (!is.null(sample_nobs) && sample_nobs != N) {
    message(paste0(
    # [CHANGE 2026-10 | Rebecca] message wording
    "\nrestriktor Message: The (specified) 'sample_nobs' (or its sum = ", sample_nobs, 
    ") differs from the sample size derived from the fitted model (", N, "). ",
    "The model-based value is used instead."
    ))
  }
  return(N)
}

VCOV.unbiased <- function(model.org, sample_nobs = NULL, ...)  {
  sample_nobs <- check_sample_nobs(sample_nobs)
  N <- NULL 
  if (!is.na(model.org$df.residual) && !is.null(model.org$df.residual) && !is.null(model.org$rank)) {
    # Use cov.mx based on N not N-k, such that output goric and gorica are the same
    # Btw if est & VCOV are used instead of fitted object, then their gorica results differ...
    N_min_k <- model.org$df.residual
    N <- N_min_k + model.org$rank
    VCOV <- vcov(model.org) * N_min_k / N
  } else if (!is.null(model.org$x) && !is.null(model.org$rank)) { 
    # Note: In rlm object model.org$df.residual is NA
    # Use cov.mx based on N not N-k, such that output goric and gorica are the same
    # Btw if est & VCOV are used instead of fitted object, then their gorica results differ...
    N <- dim(model.org$x)[1] 
    N_min_k <- N - model.org$rank
    VCOV <- vcov(model.org)*N_min_k/N
  } else {
    VCOV <- vcov(model.org)
    # TO DO Voor als dpoMatrix (kan dat hierboven ook gebeuren?), dan ws:
    #as.matrix(suppressWarnings(vcov(object))) # Behoud dit zijn namen ook?
    message(
    "\nrestriktor Message: The covariance matrix of the estimates was obtained via ",
    "'vcov()'. This represents the biased (restricted) sample covariance matrix, ",
    "not the unbiased version based on the full sample size ('N')."
    )
# TO DO als pdf van Rmd file maak, dan loopt dit door....
    # is dan ms toch een Rmd instelling....
    }
  # Check on N
  #if(!is.null(N) && sample_nobs != N) {
  if (!is.null(N) && !is.null(sample_nobs) && sample_nobs != N) {
    message(paste0(
    # [CHANGE 2026-10 | Rebecca] message wording
    "\nrestriktor Message: The (specified) 'sample_nobs' (or its sum = ", sample_nobs, 
    ") differs from the sample size determined from the fitted model (", N, "). ",
    "The unbiased covariance matrix is computed using the model-based value."
    ))
  }
  
  return(VCOV)
}

message.VCOV <- function(...)  {
  message(
  "\nrestriktor Message: The covariance matrix of the estimates was obtained via ",
  "'vcov()'. This is the biased (restricted) sample covariance matrix, ",
  "not the unbiased version based on the total sample size ('N')."
  )
} 

message.VCOVvb <- function(...)  {
  message(
  "\nrestriktor Message: The covariance matrix of the estimates was obtained via ",
  "the 'vb' argument from the metafor package."
  )
}

check.type <- function(type, class, ...)  {
  # [CHANGE 2026-10 | audit] type is case-insensitive (as in goric.default())
  # case-insensitive, like goric.default() (otherwise e.g. "GORICAC" would
  # silently be turned into "gorica")
  type <- tolower(type)
  # [/CHANGE 2026-10]
  if (type == "goric") {
    message("\nrestriktor Message: object of class ", class, " is only supported for",
            "type = 'gorica(c)'. The GORICA will be used, not the the GORIC.")
    type = "gorica"
  } else if (type == "goricc") {
    message("\nrestriktor Message: object of class ", class, " is only supported for", 
            "type = 'gorica(c)'. The GORICAC will be used, not the the GORICC.")
    type = "goricac"
  } else if (!c(type %in% c("gorica", "goricac"))) {
    message("\nrestriktor Message: object of class ", class, " is only supported for",
            "type = 'gorica(c)'. The GORICA will be used.")
    type = "gorica"
  } 
  return(type)
}

# [CHANGE 2026-10 | audit] E7: new function weight_ratio_matrix(): ratio matrix of weights (0 and Inf handled)
# ratio of weights w_i / w_j. A zero weight gives 0 (row) or Inf (column);
# 0/0 (NaN) only occurs for two zero weights, the diagonal is set to 1.
weight_ratio_matrix <- function(w, modelnames) {
  rw <- outer(w, w, "/")
  diag(rw) <- 1
  rownames(rw) <- modelnames
  colnames(rw) <- paste0("vs. ", modelnames)
  rw
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B9/E7/B1: calculate_model_comparison_metrics(): priorICweights and type arguments; IC column by exact name (no partial matching); all weights via ic_weights_log(); Heq/unconstrained rows by exact name
calculate_model_comparison_metrics <- function(x, priorICweights, type = NULL) {
  modelnames <- as.character(x$model)
  # the IC column is named after the type ('goric', 'gorica', ...); it is
  # selected explicitly (x$goric would only partially match 'gorica')
  if (is.null(type)) {
    type <- intersect(c("goric", "goricc", "gorica", "goricac"), names(x))[1]
  }
  IC <- x[[type]]
  # All weights are computed on the log scale (log-sum-exp, see
  # ic_weights_log()): a zero prior gives a zero weight (not NaN) and large
  # IC differences do not result in 0/0.
  ## Log-likelihood
  loglik_weights = ic_weights_log(-2 * x$loglik)
  loglik_rw = weight_ratio_matrix(loglik_weights, modelnames)
  
  ## penalty
  penalty_weights = ic_weights_log(2 * x$penalty)
  penalty_rw = weight_ratio_matrix(penalty_weights, modelnames)
  
  ## goric
  goric_weights = ic_weights_log(IC, priorICweights)
  goric_rw = weight_ratio_matrix(goric_weights, modelnames)
  
  # if user specified hypotheses is >= 2 and comparison = unconstrained
  # add extra column with goric weights excluding unconstrained model.
  # The special rows are identified by their exact (reserved) names; user
  # hypotheses cannot carry these names (see goric.default()).
  mn_unc_idx <- which(modelnames == "unconstrained")
  if (length(modelnames) > 2 && length(mn_unc_idx) == 1L && 
      which.max(goric_weights) != mn_unc_idx) {
    goric_weights_without_unc = ic_weights_log(IC[-mn_unc_idx], 
                                               priorICweights[-mn_unc_idx])
    goric_weights_without_unc <- append(goric_weights_without_unc, NA, 
                                        after = mn_unc_idx - 1L)
  } else { goric_weights_without_unc <- NULL }
  
  mn_heq_idx <- which(modelnames == "Heq")
  if (length(modelnames) > 2 && length(mn_heq_idx) == 1L && 
      which.max(goric_weights) != mn_heq_idx) {
    goric_weights_without_heq = ic_weights_log(IC[-mn_heq_idx], 
                                               priorICweights[-mn_heq_idx])
    goric_weights_without_heq <- append(goric_weights_without_heq, NA, 
                                        after = mn_heq_idx - 1L)
                                        # [/CHANGE 2026-10]
  } else { goric_weights_without_heq <- NULL }
  
  out <- list(loglik_weights = loglik_weights, 
              penalty_weights = penalty_weights,
              goric_weights = goric_weights,
              goric_weights_without_unc = goric_weights_without_unc,
              goric_weights_without_heq = goric_weights_without_heq,
              loglik_rw = loglik_rw,
              penalty_rw = penalty_rw,
              # [CHANGE 2026-10 | Rebecca] priorICweights returned
              goric_rw = goric_rw,
              priorICweights = priorICweights)
  
  return(out)
}

# compute penalty term, where range restrictions are treated as ceq.
PT_Amat_meq <- function(Amat, meq) {
  PT_meq  <- meq
  # check for linear dependence
  RREF <- GaussianElimination(t(Amat)) # qr(Amat)$rank
  # remove linear dependent rows
  PT_Amat <- Amat[RREF$pivot, , drop = FALSE] 
  
  if (nrow(Amat) > 1) {
    # check for range restrictions, e.g., -1 < beta < 1
    idx_range_restrictions <- detect_range_restrictions(Amat)
    # [CHANGE 2026-10 | audit] B7: design note on range restrictions and the PT
    # range restrictions are treated as equalities for computing PT (goric).
    # Note (design): a pair of opposite rows is a range irrespective of the
    # bounds, so 'x1 > 1; x1 < 1' obtains the PT of the equality x1 = 1, and
    # an inactive (non-binding) bound added to an inequality (e.g., 'x1 > 0;
    # x1 < 100') changes the PT from that of one inequality to that of one
    # equality. Only linearly independent rows are kept above.
    # [/CHANGE 2026-10]
    n_range_restrictions <- nrow(idx_range_restrictions)
    PT_meq <- meq + n_range_restrictions
    # reorder PT_Amat: ceq first, ciq second, needed for QP.solve()
    meq_order_idx <- RREF$pivot %in% c(idx_range_restrictions)
    PT_Amat <- rbind(PT_Amat[meq_order_idx, ], PT_Amat[!meq_order_idx, ])
  }
  
  return(list(PT_meq = PT_meq, PT_Amat = PT_Amat, RREF = RREF))
}


# Create a function to sort the elements in each string
sort_combination <- function(combination) {
  split_combination <- strsplit(combination, " vs. ")[[1]]
  #sorted_combination <- sort(split_combination)
  # Check if the combination includes "complement"
  if ("complement" %in% split_combination) {
    # Find the other element that is not "complement"
    other_element <- split_combination[split_combination != "complement"]
    sorted_combination <- c(other_element, "complement")
  } else if ("unconstrained" %in% split_combination) {
    # Find the other element that is not "unconstrained"
    other_element <- split_combination[split_combination != "unconstrained"]
    sorted_combination <- c(other_element, "unconstrained")
  } else {
    sorted_combination <- sort(split_combination)
  }
  paste(sorted_combination, collapse = " vs. ")
}


calculate_weight_bar <- function(Amat, meq, VCOV, mix_weights, seed, control,
                                 verbose, ...) {
  wt.bar <- NA
  if (nrow(Amat) == meq) {
    # equality constraints only
    wt.bar <- rep(0L, ncol(VCOV) + 1)
    wt.bar.idx <- ncol(VCOV) - qr(Amat)$rank + 1
    wt.bar[wt.bar.idx] <- 1
  } else if (all(c(Amat) == 0)) { 
    # unrestricted case
    wt.bar <- c(rep(0L, ncol(VCOV)), 1)
  } else if (mix_weights == "boot") { 
    # compute chi-square-bar weights based on Monte Carlo simulation
    wt.bar <- con_weights_boot(VCOV = VCOV,
                               Amat = Amat, 
                               meq  = meq, 
                               R    = ifelse(is.null(control$mix_weights_bootstrap_limit), 
                                             1e5L, control$mix_weights_bootstrap_limit),
                               seed = seed,
                               convergence_crit = ifelse(is.null(control$convergence_crit), 
                                                         5e-03, control$convergence_crit),
                               chunk_size = ifelse(is.null(control$chunk_size), 
                                                   5000L, control$chunk_size),
                               verbose = verbose, ...)
    attr(wt.bar, "mix_weights_bootstrap_limit") <- control$mix_weights_bootstrap_limit 
  } else if (mix_weights == "pmvnorm" && meq < nrow(Amat)) {
    # compute chi-square-bar weights based on pmvnorm
    wt.bar <- rev(con_weights(Amat %*% VCOV %*% t(Amat), meq = meq, 
                              tolerance = ifelse(is.null(control$tolerance), 1e-15, control$tolerance), 
                              ridge_constant = ifelse(is.null(control$ridge_constant), 1e-05, control$ridge_constant), 
                              ...))
    
    # Check if wt.bar contains NaN values
    if (any(is.nan(wt.bar))) {
      mix_weights <- "boot"
      wt.bar <- con_weights_boot(VCOV = VCOV,
                                 Amat = Amat, 
                                 meq  = meq, 
                                 R    = ifelse(is.null(control$mix_weights_bootstrap_limit), 
                                               1e5L, control$mix_weights_bootstrap_limit),
                                 seed = seed,
                                 convergence_crit = ifelse(is.null(control$convergence_crit), 
                                                           5e-03, control$convergence_crit),
                                 chunk_size = ifelse(is.null(control$chunk_size), 
                                                     5000L, control$chunk_size),
                                 verbose = verbose, ...)
      attr(wt.bar, "mix_weights_bootstrap_limit") <- control$mix_weights_bootstrap_limit   
    } 
  }
  return(wt.bar)
}



# Functie om rijen te filteren uit de parameter tabel
# extract_constraints <- function(parameter_table, hypotheses) {
#   # 1. Verwijder alle spaties in hypotheses
#   clean_hypotheses <- gsub("\\s+", "", hypotheses)
#   
#   # 2. Definieer de model- en constraint-operators
#   model_operators <- c("=~", "<~", "~*~", "~~", "~", "\\|", "%")
#   constraint_operators <- c("<", ">", "=", "==", ":=")
#   
#   # 3. Splits de hypotheses in afzonderlijke constraints
#   subconstraints <- unlist(strsplit(clean_hypotheses, ",|;|&|\\n"))
#   
#   # Functie om model-onderdelen te parseren
#   parse_model <- function(expression) {
#     for (op in model_operators) {
#       if (grepl(op, expression, fixed = TRUE)) {
#         parts <- unlist(strsplit(expression, op, fixed = TRUE))
#         if (length(parts) == 2) {
#           return(list(lhs = parts[1], op = op, rhs = parts[2]))
#         }
#       }
#     }
#     return(NULL)
#   }
#   
#   # Parse alle subconstraints
#   parsed_constraints <- lapply(subconstraints, function(constraint) {
#     for (c_op in constraint_operators) {
#       if (grepl(c_op, constraint, fixed = TRUE)) {
#         terms <- unlist(strsplit(constraint, c_op, fixed = TRUE))
#         if (length(terms) == 2) {
#           return(list(lhs = parse_model(terms[1]), 
#                       c_op = c_op, 
#                       rhs = parse_model(terms[2])))
#         }
#       }
#     }
#     return(list(lhs = parse_model(constraint), c_op = NULL, rhs = NULL))
#   })
#   
#   # Filter parameter_table voor elke constraint
#   results <- lapply(parsed_constraints, function(pc) {
#     if (!is.null(pc$lhs) && !is.null(pc$rhs)) {
#       # Filter met beide zijden van de constraint
#       parameter_table[
#         (parameter_table$lhs == pc$lhs$lhs & parameter_table$op == pc$lhs$op & parameter_table$rhs == pc$lhs$rhs) |
#           (parameter_table$lhs == pc$rhs$lhs & parameter_table$op == pc$rhs$op & parameter_table$rhs == pc$rhs$rhs), ]
#     } else if (!is.null(pc$lhs)) {
#       # Filter alleen met lhs (indien rhs ontbreekt)
#       parameter_table[
#         parameter_table$lhs == pc$lhs$lhs & parameter_table$op == pc$lhs$op & parameter_table$rhs == pc$lhs$rhs, ]
#     } else {
#       NULL
#     }
#   })
#   
#   # Combineer alle resultaten
#   result <- do.call(rbind, results)
#   
#   return(result)
# }
# 
# 
# # Voorbeeldgebruik
# hypotheses <- "dem60=~y2>dem65=~y6, dem60=~y3>dem65=~y7,dem60=~y4>dem65=~y8"
# 
# model1 <- '
#     A =~ Ab + Al + Af + An + Ar + Ac 
#     B =~ Bb + Bl + Bf + Bn + Br + Bc 
# '
# # Use the lavaan sem function to execute the confirmatory factor analysis
# fit1 <- sem(model1, data = sesamesim, std.lv = TRUE)
# 
# parameter_table <- parameterTable(fit1)
# hypotheses1 <-
#   " A=~Ab > .6 & A=~Al > .6 & A=~Af > .6 & A=~An > .6 & A=~Ar > .6 & A=~Ac >.6 & 
# B=~Bb > .6 & B=~Bl > .6 & B=~Bf > .6 & B=~Bn > .6 & B=~Br > .6 & B=~Bc >.6"
# 
# 
# model2 <- '
#     A  =~ Ab + Al + Af + An + Ar + Ac 
#     B =~ Bb + Bl + Bf + Bn + Br + Bc
# 
#     A ~ B + age + peabody
# '
# fit2 <- sem(model2, data = sesamesim, std.lv = TRUE)
# hypotheses2 <- "A~B > A~peabody = A~age = 0; 
#                A~B > A~peabody > A~age = 0; 
# A~B > A~peabody > A~age > 0"
# parameter_table <- parameterTable(fit2)
# extract_constraints(parameter_table, hypotheses2)



# [CHANGE 2026-10 | audit] B3/B4: new function goric_heq_constraints(): Heq of a hypothesis (redundant/implied inequalities removed, ranges refused, NULL without inequalities)
# Construct the equality-restricted hypothesis (Heq) belonging to an 
# order-restricted hypothesis by replacing '<' and '>' by '='. Before that,
# inequality restrictions that are identical to another one (only the one 
# with the largest rhs is kept, e.g., x1 > 0.2 is removed when x1 > 0.5 is 
# also specified) or that are implied by an equality restriction are removed,
# since the equality versions of those would conflict (x1 = 0.5 and x1 = 0.2).
# Only exact duplicates of the Jacobian rows are recognised: an inequality
# that is implied by a combination of others, or a scalar multiple of another
# one, is not. In that case Heq is inconsistent and the error is caught when
# fitting Heq (see fit_hypothesis()).
# Returns NULL when no inequality restriction remains (Heq would be identical
# to the hypothesis itself), and gives an error for a range restriction
# (0 < x1 < 1), for which Heq is not defined. If the syntax cannot be parsed,
# the plain substitution is returned.
goric_heq_constraints <- function(object, hypothesis) {
  heq_error <- function(...) {
    stop(structure(class = c("restriktor_heq_error", "error", "condition"),
                   list(message = paste0(...), call = NULL)))
  }
  
  # a constraint matrix (list(constraints =, rhs =, neq =)): all rows become
  # equality restrictions, after the same redundancy/range checks
  if (is.list(hypothesis)) {
    return(goric_heq_constraints_matrix(hypothesis, heq_error))
  }
  if (!is.character(hypothesis)) {
    stop("\nrestriktor ERROR: Heq = TRUE requires a character hypothesis or a ",
         "list with a constraint matrix (constraints, rhs, neq).", call. = FALSE)
  }
  Hceq_default <- gsub("<|>", "=", hypothesis)
  
  Hceq <- tryCatch({
    singles <- unlist(lapply(hypothesis, function(h) {
      syntax <- process_constraint_syntax(clean_constraints(h))
      unlist(lapply(lapply(syntax, process_abs_and_expand), 
                    expand_compound_constraints))
    }), use.names = FALSE)
    singles <- unique(singles[singles != ""])
    
    is_def  <- grepl(":=", singles, fixed = TRUE)
    is_ineq <- grepl("[<>]", singles) & !is_def
    is_eq   <- grepl("=", singles, fixed = TRUE) & !is_def & !is_ineq
    
    if (sum(is_ineq) < 1L) {
      # no inequality restrictions: Heq would equal the hypothesis itself
      NULL
    } else if (sum(is_ineq) < 2L && !any(is_eq)) {
      Hceq_default
    } else {
      # Amat row and rhs of each single restriction
      get_row <- function(s) {
        defs <- singles[is_def]
        cc <- con_constraints(object, constraints = paste(c(defs, s), collapse = "\n"))
        if (NROW(cc$Amat) != 1L) {
          stop("not a single linear restriction")
        }
        list(key     = paste(c(cc$Amat) + 0, collapse = "|"), 
             key_neg = paste(-c(cc$Amat) + 0, collapse = "|"),
             rhs     = cc$bvec)
      }
      idx_ineq <- which(is_ineq)
      rows_ineq <- lapply(singles[idx_ineq], get_row)
      keys_eq <- unlist(lapply(singles[is_eq], function(s) {
        r <- get_row(s)
        c(r$key, r$key_neg)
      }))
      keys_ineq <- vapply(rows_ineq, `[[`, character(1), "key")
      keys_ineq_neg <- vapply(rows_ineq, `[[`, character(1), "key_neg")
      rhs_ineq  <- vapply(rows_ineq, `[[`, numeric(1), "rhs")
      
      # range restrictions (a lower and an upper bound on the same linear
      # combination, e.g., 0.2 < x1 < 0.5): Heq is not defined
      is_range <- keys_ineq_neg %in% keys_ineq
      if (any(is_range)) {
        # the same bound on both sides (x1 > 1; x1 < 1) is an equality
        same_bound <- vapply(seq_along(keys_ineq), function(i) {
          j <- which(keys_ineq == keys_ineq_neg[i])
          length(j) > 0L && isTRUE(all.equal(rhs_ineq[i], -rhs_ineq[j[1]]))
        }, logical(1))
        if (any(is_range & !same_bound)) {
          heq_error("\nrestriktor ERROR: The equality-restricted hypothesis ",
                    "(Heq) cannot be formed: the hypothesis contains a range ",
                    "restriction (", 
                    paste(singles[idx_ineq][is_range & !same_bound], collapse = "; "),
                    "), for which there is no equality version. ",
                    "Use Heq = FALSE or compare the hypothesis to an ",
                    "equality-restricted hypothesis of your own.")
        }
      }
      
      drop <- rep(FALSE, length(singles))
      # inequalities implied by an equality restriction
      drop[idx_ineq[keys_ineq %in% keys_eq]] <- TRUE
      # identical inequality rows: keep the one with the largest rhs
      ord <- order(-rhs_ineq)
      dupl <- duplicated(keys_ineq[ord])
      drop[idx_ineq[ord][dupl]] <- TRUE
      
      if (!any(is_ineq & !drop)) {
        # all inequalities are implied by the equality restrictions
        NULL
      } else if (!any(drop)) {
        Hceq_default
      } else {
        gsub("<|>", "=", paste(singles[!drop], collapse = "; "))
      }
    }
  }, error = function(e) {
    if (inherits(e, "restriktor_heq_error")) {
      stop(e)
    }
    Hceq_default
  })
  
  Hceq
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B4: new function goric_heq_constraints_matrix(): Heq for constraint-matrix hypotheses
# Heq for a hypothesis given as a constraint matrix (list(constraints =, 
# rhs =, neq =)): the inequality rows become equality rows. As for the 
# character version, inequality rows identical to an equality row or to 
# another inequality row (the one with the largest rhs is kept) are removed
# first; a range restriction (a row and its negative) gives an error.
# Returns NULL when no inequality row remains.
goric_heq_constraints_matrix <- function(hypothesis, heq_error) {
  names(hypothesis) <- tolower(names(hypothesis))
  Amat <- hypothesis$constraints
  if (is.null(Amat)) {
    stop("\nrestriktor ERROR: The list objects must be named 'constraints', ",
         "'rhs' and 'neq'.", call. = FALSE)
  }
  if (!is.matrix(Amat)) {
    Amat <- rbind(Amat)
  }
  bvec <- if (is.null(hypothesis$rhs)) rep(0, nrow(Amat)) else hypothesis$rhs
  meq  <- if (is.null(hypothesis$neq)) 0L else as.integer(hypothesis$neq)
  if (length(bvec) != nrow(Amat) || meq < 0L || meq > nrow(Amat)) {
    stop("\nrestriktor ERROR: 'rhs' must have one element per row of ",
         "'constraints' and 'neq' cannot exceed the number of rows.", call. = FALSE)
  }
  
  key     <- apply(Amat, 1, function(r) paste(r + 0, collapse = "|"))
  key_neg <- apply(Amat, 1, function(r) paste(-r + 0, collapse = "|"))
  is_eq   <- seq_len(nrow(Amat)) <= meq
  idx_ineq <- which(!is_eq)
  if (length(idx_ineq) == 0L) {
    return(NULL)
  }
  keys_ineq     <- key[idx_ineq]
  keys_ineq_neg <- key_neg[idx_ineq]
  rhs_ineq      <- bvec[idx_ineq]
  keys_eq       <- c(key[is_eq], key_neg[is_eq])
  
  # range restrictions (a row and its negative, e.g., x1 > 0.2 and -x1 > -0.5)
  is_range <- keys_ineq_neg %in% keys_ineq
  if (any(is_range)) {
    same_bound <- vapply(seq_along(keys_ineq), function(i) {
      j <- which(keys_ineq == keys_ineq_neg[i])
      length(j) > 0L && isTRUE(all.equal(rhs_ineq[i], -rhs_ineq[j[1]]))
    }, logical(1))
    if (any(is_range & !same_bound)) {
      heq_error("\nrestriktor ERROR: The equality-restricted hypothesis ",
                "(Heq) cannot be formed: the constraint matrix contains a range ",
                "restriction (row(s) ", 
                paste(idx_ineq[is_range & !same_bound], collapse = ", "),
                "), for which there is no equality version. ",
                "Use Heq = FALSE or compare the hypothesis to an ",
                "equality-restricted hypothesis of your own.")
    }
  }
  
  drop <- rep(FALSE, nrow(Amat))
  drop[idx_ineq[keys_ineq %in% keys_eq]] <- TRUE
  ord <- order(-rhs_ineq)
  drop[idx_ineq[ord][duplicated(keys_ineq[ord])]] <- TRUE
  if (!any(!is_eq & !drop)) {
    return(NULL)
  }
  list(constraints = Amat[!drop, , drop = FALSE], rhs = bvec[!drop], 
       neq = sum(!drop))
}
# [/CHANGE 2026-10]

# [CHANGE 2026-10 | audit] B4: new function fit_hypothesis(): fit one hypothesis, clear error when Heq is infeasible
# Fit one hypothesis (restriktor() or con_gorica_est()). For the generated
# Heq hypothesis, an error about an inconsistent set of equality restrictions
# (e.g., 'x1 > x2 > 0.1; x1 > 0.05' -> 'x1 = x2 = 0.1; x1 = 0.05') is
# reported with a clear message instead of the raw quadprog/restriktor error;
# any other error is re-thrown unchanged.
fit_hypothesis <- function(name, fun, args) {
  if (!identical(name, "Heq")) {
    return(do.call(fun, args))
  }
  tryCatch(do.call(fun, args), error = function(e) {
    msg <- trimws(gsub("\\s+", " ", conditionMessage(e)))
    if (!grepl("constraints are inconsistent|constraints are conflicting|cannot exceed the number of constraints", msg)) {
      stop(e)
    }
    heq_txt <- if (is.character(args$constraints)) {
      paste(args$constraints, collapse = "; ")
    } else {
      "all rows of the constraint matrix as equalities"
    }
    stop("\nrestriktor ERROR: The equality-restricted hypothesis (Heq) cannot ",
         "be formed for this hypothesis: replacing its inequality restrictions ",
         "by equalities (", heq_txt, ") gives ",
         "a set of restrictions that cannot be fitted (", msg, "). ",
         "This happens when an inequality is implied by (a combination of) ",
         "the other restrictions. Use Heq = FALSE or remove the redundant ",
         "restriction(s).", call. = FALSE)
  })
}
# [/CHANGE 2026-10]
