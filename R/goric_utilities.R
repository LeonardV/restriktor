coef.con_goric <- function(object, ...)  {
  return(object$ormle$b.restr)
}

coef.gorica_est <- function(object, ...)  {
  return(object$b.restr)
}

# TO DO check iig messages
coef.named.vector <- function(x, VCOV = NULL, ...)  {
  # TO DO wat als vcov niet bestaat of nog niet matrix (maar bijv dpo)?
  # ms zelfs suppressWarnings(vcov(...)) gebruiken dan
  if (!is.vector(coef(x))) {
    # Some fit objects (like mlm) render matrices (>1 x >1 matrices).
    # So, if not a vector, make it one, with names(!) coming from vcov.
    est <- as.vector(coef(x))
    if (is.null(VCOV)) {
      names(est) <- rownames(vcov(x))
      # TO DO dit werkt niet, er is niet alleen : maar ook _ en laatste mag ws niet (runt iig niet)
      message("\nRestriktor Message: The estimates from the fit object -- i.e, coef(fit) -- are vectorized, ", 
                    "and their names are based on the rownames used in vcov(). ",
                    "Hence, these names should be used in the hypotheses specification ",
                    "(where ':' should be replaced by '.'.")
    } else { 
      # TO DO eigenlijk nog checken of rownames bestaan
      names(est) <- rownames(VCOV)
      message("\nRestriktor Message: The estimates from the fit object -- i.e, coef(fit) -- are vectorized, ", 
              "and their names are based on the rownames from their covariance matrix. ",
              "Hence, these names should be used in the hypotheses specification ",
              "(where ':' should be replaced by '.'.")
    }
  } else {
    est <- coef(x)
  }
  return(est)
}

VCOV.unbiased <- function(model.org, ...)  {
  # TO DO Controleer of goede check en message 
  
  if (!is.na(model.org$df.residual) && !is.null(model.org$df.residual) && !is.null(model.org$rank)) {
    # Use cov.mx based on N not N-k, such that output goric and gorica are the same
    # Btw if est & VCOV are used instead of fitted object, then their gorica results differ...
    N_min_k <- model.org$df.residual
    N <- N_min_k + model.org$rank
    VCOV <- vcov(model.org)*N_min_k/N
  } else if (!is.null(model.org$x) && !is.null(model.org$rank)) { 
    # Note: In rlm object model.org$df.residual is NA
    # Use cov.mx based on N not N-k, such that output goric and gorica are the same
    # Btw if est & VCOV are used instead of fitted object, then their gorica results differ...
    N <- dim(model.org$x)[1] 
    N_min_k <- N - model.org$rank
    VCOV <- vcov(model.org)*N_min_k/N
  } else {
    VCOV <- vcov(model.org)
    # TO DO Voor als dpo (kan dat hierboven ook?), dan ws:
    #as.matrix(suppressWarnings(vcov(object))) # Behoud idt zijn namen ook?
    message("Note: The covariance matrix of the estimates is obtained via the vcov function; therefore, 
                it is the (biased) restricted sample covariance matrix and not the unbiased sample covariance matrix (based on 'N').")
  }
  
  return(VCOV)
}
  

# TO DO check message
message.VCOV <- function(...)  {
  message("\nRestriktor Message: The covariance matrix of the estimates is obtained via the vcov function; ", 
                       "therefore, it is the (biased) restricted sample covariance matrix and not the unbiased sample covariance matrix (based on 'N').")
} 

# TO DO check message
message.VCOVvb <- function(...)  {
  message("\nRestriktor Message: The covariance matrix of the estimates is obtained via the 'vb' argument from the metaforpackage.")
}

  

calculate_model_comparison_metrics <- function(x) {
  modelnames <- as.character(x$model)
  ## Log-likelihood
  LL = -2 * x$loglik
  delta_LL = LL - min(LL)
  loglik_weights = exp(0.5 * -delta_LL) / sum(exp(0.5 * -delta_LL))
  loglik_rw = loglik_weights %*% t(1/loglik_weights)
  diag(loglik_rw) = 1
  
  ## penalty
  penalty_weights = exp(-x$penalty) / sum(exp(-x$penalty))
  penalty_rw = penalty_weights %*% t(1/penalty_weights)
  diag(penalty_rw) = 1
  
  ## goric
  delta_goric = x$goric - min(x$goric)
  goric_weights = exp(0.5 * -delta_goric) / sum(exp(0.5 * -delta_goric))
  goric_rw = goric_weights %*% t(1/goric_weights)
  diag(goric_rw) = 1
  
  # if user specified hypotheses is >= 2 and comparison = unconstrained
  # add extra column with goric weights excluding unconstrained model.
  mn_unc_idx <- grep("unconstrained", modelnames)
  if (length(modelnames) > 2 && length(mn_unc_idx) > 0 && 
      which.max(goric_weights) != mn_unc_idx) {
    delta_goric = x$goric[-mn_unc_idx] - min(x$goric[mn_unc_idx])
    goric_weights_without_unc = exp(0.5 * -delta_goric) / 
      sum(exp(0.5 * -delta_goric))
    goric_weights_without_unc <- c(goric_weights_without_unc, NA)
  } else { goric_weights_without_unc <- NULL }
  
  mn_heq_idx <- grep("Heq", modelnames)
  if (length(modelnames) > 2 && length(mn_heq_idx) > 0 && 
      which.max(goric_weights) != mn_heq_idx) {
    delta_goric = x$goric[-mn_heq_idx] - min(x$goric[mn_heq_idx])
    goric_weights_without_heq = exp(0.5 * -delta_goric) / 
      sum(exp(0.5 * -delta_goric))
    goric_weights_without_heq <- c(NA, goric_weights_without_heq)
  } else { goric_weights_without_heq <- NULL }
  
  rownames(goric_rw) = modelnames
  rownames(penalty_rw) = modelnames
  rownames(loglik_rw) = modelnames
  colnames(goric_rw) = paste0("vs. ", modelnames)
  colnames(penalty_rw) = paste0("vs. ", modelnames)
  colnames(loglik_rw) = paste0("vs. ", modelnames)
  
  out <- list(loglik_weights = loglik_weights, 
              penalty_weights = penalty_weights,
              goric_weights = goric_weights,
              goric_weights_without_unc = goric_weights_without_unc,
              goric_weights_without_heq = goric_weights_without_heq,
              loglik_rw = loglik_rw,
              penalty_rw = penalty_rw,
              goric_rw = goric_rw)
  
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
                                                         1e-03, control$convergence_crit),
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
                                                           1e-03, control$convergence_crit),
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

