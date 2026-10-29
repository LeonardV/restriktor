# Group (cell) means for lm objects ------------------------------------------
#
# group_means() expresses the group means of a fitted lm in terms of its
# regression coefficients, so that order-restricted hypotheses about group
# means can be specified (in goric) also when the model contains an intercept
# and/or more than one grouping variable (factor).
#
# - Grouping variables (factors) are: factor/ordered, character and logical
#   variables, and numeric variables taking on exactly two distinct values
#   (e.g., 0/1 dummies). All other (continuous) predictors are not used to
#   define the groups: their columns are set to 0, so the group means are the
#   (adjusted) means at a value of 0 for the continuous predictors (center
#   them to obtain the adjusted means at their average).
# - Each group mean is a linear combination of the coefficients, computed
#   from the model matrix (so it works for any coding/contrasts, with or
#   without intercept, and with or without interactions).
# - The result can be printed (to see the notation) and the syntax can be
#   copied into a hypothesis; goric.lm also adds the needed definitions
#   automatically when a hypothesis refers to these group-mean names.


group_means <- function(object, labels = c("number", "abbrev", "level"),
                        prefix = "mu", sep = "_", digits = 10) {

  labels <- match.arg(labels)

  if (!inherits(object, "lm")) {
    stop("\nrestriktor ERROR: group_means() is currently only available for ",
         "objects of class lm (incl. mlm, aov and glm).", call. = FALSE)
  }
  is_mlm <- inherits(object, "mlm")

  mf  <- model.frame(object)
  trm <- terms(object)
  resp_idx <- attr(trm, "response")
  vars <- as.character(attr(trm, "variables"))[-1]
  pred_vars <- if (resp_idx > 0L) vars[-resp_idx] else vars

  # which predictors are grouping variables (factors)?
  is_group <- vapply(pred_vars, function(v) {
    x <- mf[[v]]
    if (is.factor(x) || is.character(x) || is.logical(x)) return(TRUE)
    if (is.numeric(x) && is.null(dim(x))) {
      return(length(unique(x[!is.na(x)])) == 2L)
    }
    FALSE
  }, logical(1))
  group_vars <- pred_vars[is_group]
  cont_vars  <- pred_vars[!is_group]

  # warn if a factor seems to be coded by hand with numeric columns
  check_hand_coded_factors(mf, group_vars, cont_vars)

  if (length(group_vars) == 0L) {
    stop("\nrestriktor ERROR: no grouping variables (factors, dummies or ",
         "other variables with two values) found in the model.", call. = FALSE)
  }

  # levels per grouping variable
  lev_list <- lapply(group_vars, function(v) {
    x <- mf[[v]]
    if (is.factor(x))         levels(droplevels(x))
    else if (is.logical(x))   c(FALSE, TRUE)
    else                      sort(unique(x[!is.na(x)]))
  })
  names(lev_list) <- group_vars

  # all combinations; the first factor varies slowest (mu_11, mu_12, ...)
  grid <- expand.grid(rev(lapply(lev_list, seq_along)), KEEP.OUT.ATTRS = FALSE)
  grid <- grid[, rev(seq_len(ncol(grid))), drop = FALSE]
  names(grid) <- group_vars
  n_cells <- nrow(grid)

  # new model frame: grouping variables at the grid values, everything else 0
  newmf <- mf[rep(1L, n_cells), , drop = FALSE]
  rownames(newmf) <- NULL
  for (v in names(newmf)) {
    if (v %in% group_vars) {
      lv <- lev_list[[v]]
      val <- lv[grid[[v]]]
      x <- mf[[v]]
      if (is.factor(x)) {
        newmf[[v]] <- factor(val, levels = levels(x), ordered = is.ordered(x))
      } else if (is.character(x)) {
        newmf[[v]] <- as.character(val)
      } else {
        newmf[[v]] <- val
      }
    } else {
      x <- mf[[v]]
      if (is.numeric(x)) {
        newmf[[v]] <- if (is.null(dim(x))) rep(0, n_cells) else
          matrix(0, n_cells, ncol(x), dimnames = list(NULL, colnames(x)))
      }
    }
  }
  attr(newmf, "terms") <- trm
  X <- model.matrix(trm, newmf, contrasts.arg = object$contrasts)

  b <- coef(object)
  b_names <- if (is_mlm) rownames(b) else names(b)
  if (!identical(colnames(X), b_names)) {
    stop("\nrestriktor ERROR: could not match the model matrix to the ",
         "coefficients.", call. = FALSE)
  }
  aliased <- if (is_mlm) apply(is.na(b), 1, any) else is.na(b)
  if (any(aliased)) {
    used <- colSums(abs(X[, aliased, drop = FALSE])) > 0
    if (any(used)) {
      warning("\nrestriktor WARNING: the coefficient(s) ",
              paste(sQuote(b_names[aliased][used]), collapse = ", "),
              " are not estimable (NA); they are left out of the group means.",
              call. = FALSE)
    }
    X <- X[, !aliased, drop = FALSE]
  }

  # coefficient names as used in restriktor's hypothesis syntax
  par_names <- gsub("[:()]", ".", colnames(X))

  expr <- apply(X, 1, function(w) lincomb_to_text(w, par_names, digits))

  # names of the group means
  sub_list <- lapply(group_vars, function(v) {
    lv <- as.character(lev_list[[v]])
    switch(labels,
           number = as.character(seq_along(lv)),
           abbrev = if (all(!is.na(suppressWarnings(as.numeric(lv))))) {
                      clean_mu_label(lv)   # numeric levels (e.g., 0.5): as is
                    } else {
                      clean_mu_label(abbreviate(lv, minlength = 1L,
                                                named = FALSE))
                    },
           level  = clean_mu_label(lv))
  })
  names(sub_list) <- group_vars

  subs <- vapply(seq_len(n_cells), function(i) {
    s <- vapply(group_vars, function(v) sub_list[[v]][grid[i, v]], "")
    if (labels == "number" && all(lengths(lev_list) < 10L)) {
      paste(s, collapse = "")              # e.g. mu_12
    } else {
      paste(s, collapse = sep)             # e.g. mu_OJ_0.5
    }
  }, "")
  mu_names <- paste0(prefix, sep, subs)
  if (anyDuplicated(mu_names)) {
    stop("\nrestriktor ERROR: the group-mean labels are not unique; ",
         "use another 'labels' option.", call. = FALSE)
  }

  tab <- data.frame(name = mu_names,
                    lapply(group_vars, function(v) lev_list[[v]][grid[[v]]]),
                    definition = unname(expr),
                    stringsAsFactors = FALSE, check.names = FALSE)
  names(tab) <- c("name", group_vars, "definition")
  for (v in group_vars) {
    if (is.factor(mf[[v]])) tab[[v]] <- as.character(tab[[v]])
  }

  # number of observations per group (cell)
  key_data <- do.call(paste, c(lapply(group_vars, function(v) {
    as.character(mf[[v]])
  }), sep = "\r"))
  key_grid <- do.call(paste, c(lapply(group_vars, function(v) {
    as.character(lev_list[[v]][grid[[v]]])
  }), sep = "\r"))
  tab$n <- as.vector(table(factor(key_data, levels = key_grid)))

  syntax <- paste0(mu_names, " := ", expr)
  names(syntax) <- mu_names

  # weights of the (non-aliased) coefficients for each group mean, i.e.,
  # group means = L %*% coef(object)[!is.na(coef(object))]
  L <- X
  rownames(L) <- mu_names

  # mlm: group means per response variable, named <response>.<mu name>
  # (e.g., 'len.mu_12'), in line with the parameter names (e.g., 'len.dose1')
  if (is_mlm) {
    resp <- gsub("[:()]", ".", colnames(b))
    q <- length(resp)
    mu_names_r <- as.vector(t(outer(resp, mu_names, paste, sep = ".")))
    expr_r <- unlist(lapply(resp, function(r) {
      apply(X, 1, function(w) {
        lincomb_to_text(w, paste0(r, ".", par_names), digits)
      })
    }), use.names = FALSE)
    tab <- do.call(rbind, lapply(seq_len(q), function(k) {
      tk <- tab
      tk$name <- paste0(resp[k], ".", mu_names)
      tk$definition <- expr_r[(k - 1L) * n_cells + seq_len(n_cells)]
      cbind(response = colnames(b)[k], tk, stringsAsFactors = FALSE)
    }))
    tab <- tab[, c("name", "response", setdiff(names(tab), c("name", "response")))]
    rownames(tab) <- NULL
    syntax <- paste0(mu_names_r, " := ", expr_r)
    names(syntax) <- mu_names_r
    L <- kronecker(diag(q), X)
    dimnames(L) <- list(mu_names_r,
                        as.vector(t(outer(resp, par_names, paste, sep = "."))))
  }

  out <- list(table = tab, syntax = syntax, L = L,
              group_vars = group_vars, continuous_vars = cont_vars,
              intercept = attr(trm, "intercept") == 1L)
  class(out) <- "goric_group_means"
  out
}


print.goric_group_means <- function(x, ...) {
  cat("\nrestriktor (", packageDescription("restriktor", fields = "Version"),
      "): group means in terms of the model parameters\n\n", sep = "")
  cat("Grouping variable(s):", paste(x$group_vars, collapse = ", "), "\n")
  if (length(x$continuous_vars) > 0L) {
    cat("Continuous predictor(s):", paste(x$continuous_vars, collapse = ", "),
        "(set to 0; the group means are adjusted means at 0)\n")
  }
  cat("\n")
  tab <- x$table[, c(setdiff(names(x$table), c("definition", "n")), "n",
                     "definition")]
  print(tab, row.names = FALSE, right = FALSE)
  cat("\nThe group-mean names (e.g., ", sQuote(x$table$name[1]),
      ") can be used directly in the goric hypotheses;\n",
      "the definitions are added automatically. ",
      "Syntax to copy (if desired):\n\n", sep = "")
  cat(paste(x$syntax, collapse = "\n"), "\n")
  invisible(x)
}


# helpers ----------------------------------------------------------------------

# coefficient weights -> text like "-.Intercept. + 2*x - dose1"
lincomb_to_text <- function(w, par_names, digits = 10) {
  w <- round(w, digits)
  nz <- which(w != 0)
  if (length(nz) == 0L) return("0")
  terms_txt <- character(length(nz))
  for (k in seq_along(nz)) {
    j <- nz[k]
    a <- abs(w[j])
    coef_txt <- if (a == 1) par_names[j] else
      paste0(format(a, digits = digits, trim = TRUE), "*", par_names[j])
    sgn <- if (w[j] < 0) "-" else "+"
    terms_txt[k] <- if (k == 1L) {
      if (sgn == "-") paste0("-", coef_txt) else coef_txt
    } else {
      paste0(" ", sgn, " ", coef_txt)
    }
  }
  paste(terms_txt, collapse = "")
}

clean_mu_label <- function(x) {
  x <- gsub("[^A-Za-z0-9._]", ".", x)
  x
}

# used in goric.lm: if the (character) hypotheses refer to group means
# (in any of the label notations of group_means()), add their definitions
expand_group_means_lm <- function(object, hypotheses, prefix = "mu") {
  if (!is.list(hypotheses)) return(hypotheses)
  is_chr <- vapply(hypotheses, is.character, logical(1))
  if (!any(is_chr)) return(hypotheses)
  txt <- paste(unlist(hypotheses[is_chr]), collapse = "\n")
  if (!grepl(prefix, txt, fixed = TRUE)) return(hypotheses)
  syn <- character(0)
  for (lab in c("number", "abbrev", "level")) {
    # (warnings, e.g., on factors coded by hand, are shown once: for the
    # first notation only)
    gm <- tryCatch(if (lab == "number") {
                     group_means(object, labels = lab, prefix = prefix)
                   } else {
                     suppressWarnings(group_means(object, labels = lab,
                                                  prefix = prefix))
                   },
                   error = function(e) NULL)
    if (!is.null(gm)) {
      syn <- c(syn, gm$syntax[setdiff(names(gm$syntax), names(syn))])
    }
  }
  if (length(syn) == 0L) return(hypotheses)
  # do not overrule actual parameter names
  par_names <- if (inherits(object, "mlm")) {
    names(mlm_coef_vector(coef(object)))
  } else {
    gsub("[:()]", ".", names(coef(object)))
  }
  syn <- syn[!names(syn) %in% par_names]
  add_group_means_to_hypotheses(hypotheses, list(syntax = syn))
}

# add the definitions of the group means used in (character) hypotheses
add_group_means_to_hypotheses <- function(hypotheses, gm) {
  mu_names <- names(gm$syntax)
  used_any <- character(0)
  out <- lapply(hypotheses, function(h) {
    if (!is.character(h)) return(h)
    h_all <- paste(h, collapse = "\n")
    # names that the user defines him/herself (via :=) are left alone
    own_def <- regmatches(h_all, gregexpr("[A-Za-z0-9._]+(?=\\s*:=)",
                                          h_all, perl = TRUE))[[1]]
    tokens <- regmatches(h_all, gregexpr("[A-Za-z0-9._]+", h_all))[[1]]
    used <- intersect(mu_names, setdiff(unique(tokens), own_def))
    if (length(used) == 0L) return(h)
    used_any <<- union(used_any, used)
    paste(c(gm$syntax[used], h_all), collapse = "\n")
  })
  attr(out, "group_means_used") <- used_any
  out
}

# mlm: coefficient matrix (parameters x responses) -> named vector, in the
# order of vcov() (response by response) and with the names used in the
# hypotheses (e.g., 'len.dose0.5', 'len..Intercept.'); other input is
# returned unchanged.
mlm_coef_vector <- function(b) {
  if (!is.matrix(b) || is.null(colnames(b)) || is.null(rownames(b))) return(b)
  out <- as.vector(b)
  names(out) <- gsub("[:()]", ".",
                     paste0(rep(colnames(b), each = nrow(b)), ":",
                            rep(rownames(b), times = ncol(b))))
  out
}


# Warn when a factor seems to be coded by hand with several numeric columns
# instead of being one factor:
# 1. effects coding by hand: numeric predictors taking (only) the values
#    -1, 0, 1; these are treated as continuous (set to 0), so the groups of
#    that factor are missed;
# 2. dummy coding by hand: 2-valued numeric predictors that never take their
#    higher value at the same time (mutually exclusive); then all
#    combinations are formed, including non-existing groups (n = 0).
check_hand_coded_factors <- function(mf, group_vars, cont_vars) {
  advice <- paste0(
    "Use one factor for these columns instead (e.g., lm(y ~ f), with f a factor). ",
    "For the group means, it does not matter whether dummy or effects coding is ",
    "used; if you also want effects-coded estimates, additionally specify the ",
    "contrasts, e.g., lm(y ~ f, contrasts = list(f = \"contr.sum\")).")
  # 1. effects coding by hand
  eff <- cont_vars[vapply(cont_vars, function(v) {
    x <- mf[[v]]
    is.numeric(x) && is.null(dim(x)) &&
      setequal(unique(x[!is.na(x)]), c(-1, 0, 1))
  }, logical(1))]
  if (length(eff) > 0L) {
    warning("\nrestriktor WARNING: the variable(s) ",
            paste(sQuote(eff), collapse = ", "),
            " only take on the values -1, 0, and 1, so they seem to be ",
            "effects-coded columns of a factor coded by hand. They are now ",
            "treated as continuous variables (set to 0), so the groups of that ",
            "factor are not distinguished. ", advice, call. = FALSE)
  }
  # 2. dummy coding by hand: mutually exclusive 2-valued numeric variables
  num2 <- group_vars[vapply(group_vars, function(v) {
    x <- mf[[v]]
    is.numeric(x) && is.null(dim(x))
  }, logical(1))]
  if (length(num2) >= 2L) {
    hi <- lapply(num2, function(v) mf[[v]] == max(mf[[v]], na.rm = TRUE))
    names(hi) <- num2
    excl <- character(0)
    for (a in seq_along(num2)) {
      for (b in seq_along(num2)) {
        if (a < b && !any(hi[[a]] & hi[[b]], na.rm = TRUE)) {
          excl <- union(excl, num2[c(a, b)])
        }
      }
    }
    if (length(excl) > 0L) {
      warning("\nrestriktor WARNING: the variable(s) ",
              paste(sQuote(excl), collapse = ", "),
              " never take on their higher value at the same time, so they ",
              "seem to be dummy-coded columns of a factor coded by hand. Then, ",
              "also non-existing groups (with n = 0) are formed. ", advice,
              call. = FALSE)
    }
  }
  invisible(NULL)
}
