coef.con_goric <- function(object, ...)  {
  return(object$ormle$b.restr)
}

coef.gorica_est <- function(object, ...)  {
  return(object$b.restr)
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
    message("\nrestriktor Message: The coefficients from the fitted model have been converted into a vector. ",
            "The coefficient names are taken from the row names of the covariance matrix (vcov). ",
            "Use these names when specifying hypotheses. ",
            "Note: replace any ':' characters with '.' in hypothesis labels."
            )
    # TO DO mlm: ws nog zeggen dat intercept dan wordt: DV..Intercept.
    # TO DO mlm: Werkt dit wel voor mlm, bij mij volgens mij niet....
  } else {
    est <- coef(x)
  }
  return(est)
}

# log-likelihood of the unrestricted model. logLik.lm() does not support
# multiple responses, so for mlm objects con_loglik_lm() is used.
loglik_unrestr <- function(model.org)  {
  if (inherits(model.org, "mlm")) {
    con_loglik_lm(model.org)
  } else {
    logLik(model.org)
  }
}

# For mlm objects the coefficients form a p x ny matrix. Convert it to a
# vector in the order of vcov(), with the names of vcov() (e.g., 'Age:GroupNo').
# Other objects are returned as is.
coef_vec_mlm <- function(b, model.org)  {
  if (inherits(model.org, "mlm") && is.matrix(b)) {
    b <- structure(as.vector(b), names = rownames(vcov(model.org)))
  }
  b
}

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
  if (!is.null(length_expected) && length(w) != length_expected) {
    stop("\nrestriktor ERROR: The argument '", name, "' should consist of ",
         length_expected, " elements",
         if (!is.null(what)) paste0(", namely ", what),
         ". It now consists of ", length(w), " elements.", call. = FALSE)
  }
  if (rescale && !isTRUE(all.equal(sum(w), 1))) {
    w <- w / sum(w)
  }
  w
}

# Numerically stable (prior-weighted) IC weights:
#   w_i = prior_i * exp(-IC_i / 2) / sum_j prior_j * exp(-IC_j / 2)
# computed via log-sum-exp. A zero prior gives a zero weight (not NaN),
# and very large IC differences do not underflow to 0/0.
ic_weights_log <- function(IC, prior = NULL) {
  if (is.null(prior)) {
    prior <- rep(1, length(IC))
  }
  lw <- -IC / 2 + log(prior)              # log(0) = -Inf for a zero prior
  lw[prior == 0] <- -Inf
  m <- max(lw)
  if (!is.finite(m)) {
    return(rep(NaN, length(IC)))
  }
  w <- exp(lw - m)
  w / sum(w)
}

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

# ratio of weights w_i / w_j. A zero weight gives 0 (row) or Inf (column);
# 0/0 (NaN) only occurs for two zero weights, the diagonal is set to 1.
weight_ratio_matrix <- function(w, modelnames) {
  rw <- outer(w, w, "/")
  diag(rw) <- 1
  rownames(rw) <- modelnames
  colnames(rw) <- paste0("vs. ", modelnames)
  rw
}

calculate_model_comparison_metrics <- function(x, priorICweights) {
  modelnames <- as.character(x$model)
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
  goric_weights = ic_weights_log(x$goric, priorICweights)
  goric_rw = weight_ratio_matrix(goric_weights, modelnames)
  
  # if user specified hypotheses is >= 2 and comparison = unconstrained
  # add extra column with goric weights excluding unconstrained model.
  mn_unc_idx <- grep("unconstrained", modelnames)
  if (length(modelnames) > 2 && length(mn_unc_idx) > 0 && 
      which.max(goric_weights) != mn_unc_idx) {
    goric_weights_without_unc = ic_weights_log(x$goric[-mn_unc_idx], 
                                               priorICweights[-mn_unc_idx])
    goric_weights_without_unc <- c(goric_weights_without_unc, NA)
  } else { goric_weights_without_unc <- NULL }
  
  mn_heq_idx <- grep("Heq", modelnames)
  if (length(modelnames) > 2 && length(mn_heq_idx) > 0 && 
      which.max(goric_weights) != mn_heq_idx) {
    goric_weights_without_heq = ic_weights_log(x$goric[-mn_heq_idx], 
                                               priorICweights[-mn_heq_idx])
    goric_weights_without_heq <- c(NA, goric_weights_without_heq)
  } else { goric_weights_without_heq <- NULL }
  
  out <- list(loglik_weights = loglik_weights, 
              penalty_weights = penalty_weights,
              goric_weights = goric_weights,
              goric_weights_without_unc = goric_weights_without_unc,
              goric_weights_without_heq = goric_weights_without_heq,
              loglik_rw = loglik_rw,
              penalty_rw = penalty_rw,
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
    # range restrictions are treated as equalities for computing PT (goric)
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



# Construct the equality-restricted hypothesis (Heq) belonging to an 
# order-restricted hypothesis by replacing '<' and '>' by '='. Redundant 
# inequality restrictions (e.g., x1 > 0.2 when x1 > 0.5 is also specified, or 
# an inequality implied by an equality restriction) are removed first. Otherwise, 
# Heq would contain conflicting equality restrictions (x1 = 0.5 and x1 = 0.2).
# If this does not work out, the plain substitution is returned.
goric_heq_constraints <- function(object, hypothesis) {
  Hceq_default <- gsub("<|>", "=", hypothesis)
  if (!is.character(hypothesis)) {
    return(Hceq_default)
  }
  
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
    
    if (sum(is_ineq) < 1L || (sum(is_ineq) < 2L && !any(is_eq))) {
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
      rhs_ineq  <- vapply(rows_ineq, `[[`, numeric(1), "rhs")
      
      drop <- rep(FALSE, length(singles))
      # inequalities implied by an equality restriction
      drop[idx_ineq[keys_ineq %in% keys_eq]] <- TRUE
      # identical inequality rows: keep the one with the largest rhs
      ord <- order(-rhs_ineq)
      dupl <- duplicated(keys_ineq[ord])
      drop[idx_ineq[ord][dupl]] <- TRUE
      
      if (!any(drop)) {
        Hceq_default
      } else {
        gsub("<|>", "=", paste(singles[!drop], collapse = "; "))
      }
    }
  }, error = function(e) Hceq_default)
  
  Hceq
}
