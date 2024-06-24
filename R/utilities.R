## utility functions
coef.restriktor <- function(object, ...)  {
  
  b.def <- c()
  b.restr <- object$b.restr
  
  if (any(object$parTable$op == ":=")) {
    b.def <- object$CON$def.function(object$b.restr)
  }
  
  if (inherits(object, "conMLM")) {
    OUT <- rbind(b.restr, b.def)
  } else {
    OUT <- c(b.restr, b.def)
  }
  
  return(OUT)
}


coef.con_goric <- function(object, ...)  {
  return(object$ormle$b.restr)
}

coef.gorica_est <- function(object, ...)  {
  return(object$b.restr)
}


logLik.restriktor <- function(object, ...) {
  return(object$loglik)
}


model.matrix.restriktor <- function(object, ...) {
  return(model.matrix(object$model.org))
}


tukeyChi <- function(x, c = 4.685061, deriv = 0, ...) {
  u <- x / c
  out <- abs(x) > c
  if (deriv == 0) { # rho function
    r <- 1 - (1 - u^2)^3
    r[out] <- 1
  } else if (deriv == 1) { # rho' = psi function
    r <- 6 * x * (1 - u^2)^2 / c^2
    r[out] <- 0
  } else if (deriv == 2) { # rho'' 
    r <- 6 * (1 - u^2) * (1 - 5 * u^2) / c^2
    r[out] <- 0
  } else {
    stop("deriv must be in {0,1,2}")
  }
  return(r)
}


# code taken from robustbase package.
# addapted by LV (3-12-2017).
robWeights <- function(w, eps = 0.1/length(w), eps1 = 0.001, ...) {
  stopifnot(is.numeric(w))
  cat("Robustness weights:", "\n")
  cat0 <- function(...) cat("", ...)
  n <- length(w)
  if (n <= 10) 
    print(w, digits = 5, ...)
  else {
    n1 <- sum(w1 <- abs(w - 1) < eps1)
    n0 <- sum(w0 <- abs(w) < eps)
    if (any(w0 & w1)) 
      warning("weights should not be both close to 0 and close to 1!\n", 
              "You should use different 'eps' and/or 'eps1'")
    if (n0 > 0 || n1 > 0) {
      if (n0 > 0) {
        formE <- function(e) formatC(e, digits = max(2, 5 - 3), width = 1)
        i0 <- which(w0)
        maxw <- max(w[w0])
        c3 <- paste0("with |weight| ", if (maxw == 0) 
          "= 0"
          else paste("<=", formE(maxw)), " ( < ", formE(eps), 
          ");")
        cat0(if (n0 > 1) {
          cc <- sprintf("%d observations c(%s)", n0, 
                        strwrap(paste(i0, collapse = ",")))
          c2 <- " are outliers"
          paste0(cc, if (nchar(cc) + nchar(c2) + nchar(c3) > 
                         getOption("width")) 
            "\n\t", c2)
        }
        else sprintf("observation %d is an outlier", 
                     i0), c3, "\n")
      }
      if (n1 > 0) 
        cat0(ngettext(n1, "one weight is", sprintf("%s%d weights are", 
                                                   if (n1 == n) 
                                                     "All "
                                                   else "", n1)), "~= 1.")
      n.rem <- n - n0 - n1
      if (n.rem <= 0) {
        if (n1 > 0) 
          cat("\n")
        return(invisible())
      }
    }
  }
}




# function taken from 'bain' package 
expand_compound_constraints <- function(hyp) {
  equality_operators <- gregexpr("[=<>]", hyp)[[1]]
  if(length(equality_operators) > 1){
    string_positions <- c(0, equality_operators, nchar(hyp)+1)
    return(sapply(1:(length(string_positions)-2), function(pos) {
      substring(hyp, (string_positions[pos]+1), (string_positions[pos+2]-1))
    }))
  } else {
    return(hyp)
  }
}

# function taken from 'bain' package 
expand_parentheses <- function(hyp) {
  parenth_locations <- gregexpr("[\\(\\)]", hyp)[[1]]
  if (!parenth_locations[1] == -1 && !grepl("abs\\(.*\\)", hyp) ) {
    if (length(parenth_locations) %% 2 > 0) stop("Not all opening parentheses are matched by a closing parenthesis, or vice versa.")
    expanded_contents <- strsplit(substring(hyp, (parenth_locations[1]+1), (parenth_locations[2]-1)), ",")[[1]]
    if (length(parenth_locations) == 2){
      return(paste0(substring(hyp, 1, (parenth_locations[1]-1)), expanded_contents, substring(hyp, (parenth_locations[2]+1), nchar(hyp))))
    } else {
      return(apply(
        expand.grid(expanded_contents, expand_parentheses(substring(hyp, (parenth_locations[2]+1), nchar(hyp)))),
        1, paste, collapse = ""))
    }
  } else {
    return(hyp)
  }
}

# used in print.benchmark()
capitalize_first_letter <- function(input_string) {
  paste0(toupper(substring(input_string, 1, 1)), substring(input_string, 2))
}

# used in get_results_benchmark_means()
remove_single_value_rows <- function(data, value) {
  rows_to_keep <- apply(data, 1, function(row) !all(row == value))
  data[rows_to_keep, , drop = FALSE]
}  


detect_intercept <- function(model) {
  coefficients <- model$b.unrestr
  intercept_names <- c("(Intercept)", "Intercept", "(const)", "const", 
                       "(Int)", "Int", "(Cons)", "Cons", "b0", "beta0")
  
  names_lower <- tolower(names(coefficients))
  intercept_names_lower <- tolower(intercept_names)
  
  detected_intercepts <- intercept_names[names_lower %in% intercept_names_lower]

  if (length(detected_intercepts) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


format_numeric <- function(x, digits = 3) {
  if (abs(x) <= 1e-8) {
    format(0, nsmall = digits)
  } else if (abs(x) >= 1e3 || abs(x) <= 1e-3) {
    format(x, scientific = TRUE, digits = digits)
  } else {
    format(round(x, digits), nsmall = digits) 
  }
}


calculate_model_comparison_metrics <- function(x) {
  modelnames <- as.character(x$model)
  ## Log-likelihood
  LL = -2 * x$loglik
  delta_LL = LL - min(LL)
  loglik_weights = exp(0.5 * -delta_LL) / sum(exp(0.5 * -delta_LL))
  loglik_rw  <- loglik_weights %*% t(1/loglik_weights)
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
  
  rownames(goric_rw)   = modelnames
  rownames(penalty_rw) = modelnames
  rownames(loglik_rw)  = modelnames
  colnames(goric_rw)   = paste0("vs. ", modelnames)
  colnames(penalty_rw) = paste0("vs. ", modelnames)
  colnames(loglik_rw)  = paste0("vs. ", modelnames)
  
  out <- list(loglik_weights  = loglik_weights, 
              penalty_weights = penalty_weights,
              goric_weights   = goric_weights,
              loglik_rw       = loglik_rw,
              penalty_rw      = penalty_rw,
              goric_rw        = goric_rw)
  
  return(out)
}

# this function is called from the goric_benchmark_anova() function
parallel_function_means <- function(i, N, var_e, means_pop, 
                                    hypos, pref_hypo, object, ngroups, sample, 
                                    control, form_model_org, ...) {  
  # Sample residuals
  epsilon <- rnorm(sum(N), sd = sqrt(var_e))
  
  # original model formula 
  if (length(form_model_org) > 0) {
    model <- form_model_org
    lhs <- all.vars(model)[1]
    sample[[lhs]] <- as.matrix(sample[, 2:(1 + ngroups)]) %*% matrix(means_pop, 
                                                                         nrow = ngroups) + epsilon
    df_boot <- data.frame(lhs = sample[[lhs]], sample[, 2:(1 + ngroups)])
    colnames(df_boot)[1] <- lhs
    
    has_intercept <- attr(terms(model), "intercept") == 1
    rhs <- as.character(attr(terms(model), "term.labels"))
    
    # Create the RHS with all other variables and optionally the intercept
    if (has_intercept) {
      new_rhs <- "."
    } else {
      new_rhs <- "-1 + ."
    }
    
    # Create the new formula
    new_model <- as.formula(paste(lhs, "~", new_rhs))
  } else {
    new_model <- y ~ 0 + .
    # Generate data
    sample$y <- as.matrix(sample[, 2:(1 + ngroups)]) %*% matrix(means_pop, 
                                                                nrow = ngroups) + epsilon
    df_boot <- data.frame(y = sample$y, sample[, 2:(1 + ngroups)])
  }
  
  
  # Obtain fit
  fit_boot <- lm(new_model, data = df_boot)  
  # GORICA or GORICA depending on what is done in data
  results_goric <- goric(fit_boot,
                         hypotheses = hypos,
                         comparison = object$comparison,
                         type = object$type,
                         control = control, 
                         ...)
  
  # Return the relevant results
  list(
    #test  = attr(results.goric$objectList[[results.goric$objectNames]]$wt.bar, "mvtnorm"),
    gw  = results_goric$result[pref_hypo, 7], # goric(a) weight
    rgw = results_goric$ratio.gw[pref_hypo, ], # ratio goric(a) weights
    rlw = results_goric$ratio.lw[pref_hypo, ], # ratio likelihood weights
    ld  = (results_goric$result$loglik[pref_hypo] - results_goric$result$loglik)
  )
}


# this function is called from the benchmark_asymp() function
parallel_function_asymp <- function(i, est, VCOV, hypos, pref_hypo, comparison,
                                    type, control, ...) {  
  
  results_goric <- goric(est[i, ], VCOV = VCOV,
                         hypotheses = hypos,
                         comparison = comparison,
                         type = type,
                         control = control, 
                         ...)
  
  # Return the relevant results
  list(
    #test  = attr(results.goric$objectList[[results.goric$objectNames]]$wt.bar, "mvtnorm"),
    gw  = results_goric$result[pref_hypo, 7], # goric(a) weight
    rgw = results_goric$ratio.gw[pref_hypo, ], # ratio goric(a) weights
    rlw = results_goric$ratio.lw[pref_hypo, ], # ratio likelihood weights
    ld  = (results_goric$result$loglik[pref_hypo] - results_goric$result$loglik)
  )
}



# Function to identify list and corresponding messages
identify_messages <- function(x) {
  messages_info <- list()
  hypo_messages <- names(x$objectList)
  for (object_name in hypo_messages) {
    if (length(x$objectList[[object_name]]$messages) > 0) {
      messages <- names(x$objectList[[object_name]]$messages)
      messages_info[[object_name]] <- messages
    } else {
      messages_info[[object_name]] <- "No messages"
    }
  }
  return(messages_info)
}


detect_range_restrictions <- function(Amat) {
  n <- nrow(Amat)
  range_restrictions <- matrix(0, ncol = 2, nrow = n)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (all(Amat[i, ] == -Amat[j, ])) {
        range_restrictions[i, ] <- c(i, j)
      }
    }
  }
  
  range_restrictions <- range_restrictions[range_restrictions[, 1] != 0, , drop = FALSE]
  return(range_restrictions)
}



PT_Amat_meq <- function(Amat, meq) {
  # check for range-restrictions and treat them as equalities
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


# correct misspecified constraints of format e.g., x1 < 1 & x1 < 2.
# x1 < 2 is removed since it is redundant. It has no impact on the LPs, but
# since the redundant matrix is not full row-rank the slower boot method is used. 
remove_redundant_constraints <- function(constraints, rhs) {
  df <- data.frame(constraints, rhs)
  df <- df[order(df$rhs, decreasing = TRUE),]  
  unique_constraints <- !duplicated(df[, -ncol(df)])
  df_reduced <- df[unique_constraints,]
  rhs <- df_reduced$rhs
  row.names(df_reduced) <- NULL
  colnames(df_reduced) <- NULL
  list(constraints = as.matrix(df_reduced[, -ncol(df_reduced)]), rhs = rhs) 
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
    wt.bar <- rev(con_weights(Amat %*% VCOV %*% t(Amat), meq = meq))
    
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
