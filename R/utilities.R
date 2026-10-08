# Functie om haakjesbalans te valideren
validate_parentheses <- function(hyp) {
  open_count <- 0
  for (char in strsplit(hyp, NULL)[[1]]) {
    if (char == "(") {
      open_count <- open_count + 1
    } else if (char == ")") {
      open_count <- open_count - 1
    }
    if (open_count < 0) {
      return(FALSE) # Meer sluitende dan open haakjes
    }
  }
  return(open_count == 0) # Gelijk aantal open en sluitende haakjes
}

# Functie om constraints op te schonen
clean_constraints <- function(constraints) {
  constraints <- gsub("[#!].*(?=\n|$)", "", constraints, perl = TRUE) # Verwijder opmerkingen
  constraints <- gsub("[;|&]", "\n", constraints, perl = TRUE) # Normaliseer delimiters naar \n
  constraints <- gsub("[ \t]+", "", constraints, perl = TRUE) # Verwijder extra spaties
  constraints <- gsub("\n{2,}", "\n", constraints, perl = TRUE) # Verwijder extra nieuwe regels
  
  # Alleen komma's buiten haakjes vervangen
  result <- ""
  open_parens <- 0
  
  for (char in strsplit(constraints, NULL)[[1]]) {
    if (char == "(") {
      open_parens <- open_parens + 1
    } else if (char == ")") {
      open_parens <- open_parens - 1
    }
    
    if (char == "," && open_parens == 0) {
      result <- paste0(result, "\n")
    } else {
      result <- paste0(result, char)
    }
  }
  return(result)
}

# Functie om haakjesinhoud uit te breiden
expand_parentheses <- function(hyp) {
  if (!validate_parentheses(hyp)) {
    stop("Not all opening parentheses are matched by a closing parenthesis, or vice versa.")
  }
  
  # Vind de locaties van haakjes
  parenth_locations <- gregexpr("[\\(\\)]", hyp)[[1]]
  
  if (parenth_locations[1] != -1) {
    # Haal de inhoud tussen de eerste set haakjes op
    inner_content <- substring(hyp, parenth_locations[1] + 1, parenth_locations[2] - 1)
    
    # Controleer of er komma's in de inhoud staan
    # Geen komma = rekenkundige groepering -> niet uitbreiden
    if (!grepl(",", inner_content)) {
      return(hyp)
    }
    
    # Komma aanwezig = bain-stijl uitbreiding
    expanded_contents <- strsplit(inner_content, ",\\s*")[[1]]
    
    expanded_strings <- sapply(expanded_contents, function(content) {
      paste0(
        substring(hyp, 1, parenth_locations[1] - 1), 
        content, 
        substring(hyp, parenth_locations[2] + 1, nchar(hyp))
      )
    }, USE.NAMES = FALSE)
    
    if (any(grepl("\\(", expanded_strings))) {
      return(unlist(lapply(expanded_strings, expand_parentheses)))
    } else {
      return(expanded_strings)
    }
  } else {
    return(hyp)
  }
}


# Functie om abs() en haakjes correct te verwerken
process_abs_and_expand <- function(hyp) {
  abs_matches <- gregexpr("abs\\((?:[^()]*+|\\((?:[^()]*+|\\([^()]*\\))*\\))*\\)", hyp, perl = TRUE)[[1]]
  
  if (abs_matches[1] != -1) {
    abs_substrings <- regmatches(hyp, list(abs_matches))[[1]]
    placeholder <- paste0("PLACEHOLDER_", seq_along(abs_substrings))
    for (i in seq_along(abs_substrings)) {
      hyp <- sub(abs_substrings[i], placeholder[i], hyp, fixed = TRUE)
    }
  } else {
    abs_substrings <- NULL
    placeholder <- NULL
  }
  
  # Sorteer placeholders en corresponderende substrings op lengte
  # Als een placeholder al een deelstring van een andere placeholder bevat 
  # (bijvoorbeeld PLACEHOLDER_10 bevat PLACEHOLDER_1), kan gsub verkeerd werken 
  # en een onbedoeld resultaat opleveren. Dit fixen we hier
  order <- order(nchar(placeholder), decreasing = TRUE)
  placeholder <- placeholder[order]
  abs_substrings <- abs_substrings[order]
  
  res <- expand_parentheses(hyp)
  
  if (!is.null(abs_substrings)) {
    for (i in seq_along(placeholder)) {
      res <- gsub(placeholder[i], abs_substrings[i], res, fixed = TRUE)
    }
  }
  return(res)
}

# Functie om constraints te verwerken
process_constraint_syntax <- function(constraint.syntax) {
  syntax <- strsplit(constraint.syntax, split = "\n", perl = TRUE)[[1]]
  syntax <- syntax[syntax != ""] # Verwijder lege strings
  syntax <- gsub("==", "=", syntax, perl = TRUE) # Normaliseer gelijkheid
  syntax <- gsub("<=", "<", syntax, perl = TRUE) # Normaliseer <= naar <
  syntax <- gsub(">=", ">", syntax, perl = TRUE) # Normaliseer >= naar >
  
  return(syntax)
}


# function taken from 'bain' package, but adapted by LV (27-11-2024) 
expand_compound_constraints <- function(hyp_list) {
  lapply(hyp_list, function(hyp) {
    equality_operators <- gregexpr("[=<>]", hyp)[[1]]
    if (length(equality_operators) > 1 && equality_operators[1] != -1) {
      string_positions <- c(0, equality_operators, nchar(hyp) + 1)
      sapply(1:(length(string_positions) - 2), function(pos) {
        substring(hyp, string_positions[pos] + 1, string_positions[pos + 2] - 1)
      })
    } else {
      hyp
    }
  })
}


# coef.restriktor <- function(object, ...)  {
#   
#   b.def <- c()
#   b.restr <- object$b.restr
#   
#   if (any(object$parTable$op == ":=")) {
#     b.def <- object$CON$def.function(object$b.restr)
#   }
#   
#   if (inherits(object, "conMLM")) {
#     OUT <- rbind(b.restr, b.def)
#   } else {
#     OUT <- c(b.restr, b.def)
#   }
#   return(OUT)
# }


coef.restriktor <- function(object, ..., which = c("restr", "unrestr"))  {
  which <- match.arg(which)
  
  # kies basis coefs
  b_name <- if (which == "restr") "b.restr" else "b.unrestr"
  if (is.null(object[[b_name]])) {
    stop(sprintf("Object heeft geen '%s'.", b_name), call. = FALSE)
  }
  b_base <- object[[b_name]]
  
  # definities (:=) toepassen op dezelfde basisvector
  b_def <- NULL
  has_defs <- !is.null(object$parTable$op) && any(object$parTable$op == ":=")
  if (has_defs) {
    if (is.null(object$CON$def.function) || !is.function(object$CON$def.function)) {
      stop("Definitie-parameters (':=') aanwezig, maar object$CON$def.function ontbreekt.", call. = FALSE)
    }
    b_def <- object$CON$def.function(b_base)
  }
  
  if (inherits(object, "conMLM")) {
    OUT <- rbind(b_base, b_def)
  } else {
    OUT <- c(b_base, b_def)
  }
  
  OUT
}



list_to_df_rows <- function(x) {
  stopifnot(is.list(x), length(x) > 0L)
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))

  mat <- do.call(rbind, lapply(x, function(v) {
    out <- setNames(rep(NA_real_, length(all_names)), all_names)
    out[names(v)] <- as.numeric(v)
    out
  }))

  as.data.frame(mat, check.names = FALSE)
}


# Appends one named numeric row ('new_row') to an existing data.frame of
# coefficients ('df', as produced by list_to_df_rows() above), aligning by
# column NAME instead of relying on raw rbind()'s positional matching.
#
# goric.R's "complement"/"unconstrained" result-assembly (coefs <- rbind(
# coefs, one_vec_full)) builds 'one_vec_full' as a hypothesis's own
# coefficients plus, when present, its defined (':=') parameter value(s) --
# but whether coef(<hypothesis>, which = "restr") itself already includes
# that defined-parameter value (and so is already reflected in 'df's
# columns) can differ depending on how the hypothesis's 'object' was
# derived upstream (e.g. a plain numeric-estimate-vector + VCOV input, the
# mode goric_benchmark.R's internal/simulated goric() calls use, vs. a
# fitted model object). When the two disagree, 'df' and 'one_vec_full' end
# up with different lengths, and raw rbind() recycles positionally instead
# of erroring -- which both warns ("number of columns of result, ... is
# not a multiple of vector length ...") *and* silently misaligns the
# recycled values into the wrong columns. Aligning by name and NA-filling
# any column one side doesn't have avoids both: the added/missing
# parameter column is just NA for the rows that don't carry it.
rbind_named_row <- function(df, new_row) {
  all_names <- union(colnames(df), names(new_row))

  df_full <- as.data.frame(
    matrix(NA_real_, nrow = nrow(df), ncol = length(all_names),
           dimnames = list(rownames(df), all_names)),
    check.names = FALSE
  )
  if (nrow(df) > 0 && ncol(df) > 0) {
    df_full[, colnames(df)] <- df
  }

  new_row_full <- setNames(rep(NA_real_, length(all_names)), all_names)
  new_row_full[names(new_row)] <- as.numeric(new_row)

  out <- rbind(df_full, new_row_full)
  rownames(out) <- NULL
  out
}

# get_defs <- function(obj, b) {
#   f <- obj$CON$def.function
#   if (!is.function(f)) return(numeric(0))
#   
#   bd <- body(f)
#   if (is.null(bd) || identical(bd, quote(NULL))) return(numeric(0))
#   
#   res <- tryCatch(
#     {
#       if (length(formals(f)) == 0L) f() else f(b)
#     },
#     error = function(e) NULL
#   )
#   
#   if (is.null(res) || length(res) == 0L) numeric(0) else res
# }


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





format_numeric <- function(x, digits = 3) {
  if (abs(x) <= 1e-8) {
    format(0, nsmall = digits)
  } else if (abs(x) >= 1e3 || abs(x) <= 1e-3) {
    format(x, scientific = TRUE, digits = digits)
  } else {
    format(round(x, digits), nsmall = digits) 
  }
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

# correct mis-specified constraints of format e.g., x1 < 1 & x1 < 2.
# x1 < 2 is removed since it is redundant. It has no impact on the LPs, but
# since the redundant matrix is not full row-rank the slower boot method is used.
#
# Rewritten to operate on plain matrices/vectors instead of routing through
# a data.frame (build one, sort it, duplicated() it, tear it back down into
# a matrix). This function is called once per hypothesis on every single
# con_constraints() call -- including, since benchmark()'s Monte Carlo loop
# calls goric() once per simulated draw, thousands of times per benchmark()
# run -- for what is typically a tiny (1-5 row) constraint matrix, so the
# data.frame construction/teardown overhead was consistently one of the
# largest single line-level costs in a benchmark() profile (around a fifth
# of total self time), well ahead of any of the actual numerical work
# (constrained optimization, GORICA penalty calculation). Behavior is
# unchanged: same sort-by-rhs-descending (order() is a stable sort either
# way, so tied rows keep their original relative order exactly as before),
# same "duplicate constraint ROWS keep only the one with the largest rhs"
# dedup rule (duplicated() on a matrix is already row-wise -- see
# duplicated.matrix/duplicated.array -- the same as duplicated() on a
# data.frame), and the same final meq (count of surviving rows that were
# originally flagged as equalities). Verified byte-for-byte equivalent to
# the previous implementation (modulo one inert attribute-representation
# detail -- see below) across several thousand randomized inputs spanning
# every combination of row count, column count, meq, and duplicate/tie
# pattern exercised here.
#
# The only difference from the previous implementation is one that no code
# anywhere could actually observe: the old data.frame-based code happened
# to attach a dimnames = list(NULL, NULL) attribute to its output matrix
# when more than one column survived, but NO dimnames attribute at all
# when exactly one column survived (an incidental side effect of
# data.frame's `[` auto-simplifying a single remaining column to a plain
# vector before as.matrix() was applied) -- i.e. its own output was
# already inconsistent on this point, depending on ncol. This version
# always attaches dimnames = list(NULL, NULL), for both cases alike.
# Either form means exactly "no row/column names": rownames()/colnames()
# return NULL for both, and every downstream use here only ever reads the
# constraint matrix's values (via indexing/matrix algebra), never its
# dimnames directly.
remove_redundant_constraints <- function(constraints, rhs, meq) {
  ord <- order(rhs, decreasing = TRUE)
  constraints_s <- constraints[ord, , drop = FALSE]
  rhs_s <- rhs[ord]

  eq_s <- rep(0, length(rhs))
  if (meq > 0) {
    eq_s[seq_len(meq)] <- 1
  }
  eq_s <- eq_s[ord]

  keep <- !duplicated(constraints_s) # row-wise; see duplicated.matrix
  out_constraints <- constraints_s[keep, , drop = FALSE]
  dimnames(out_constraints) <- list(NULL, NULL)

  list(constraints = out_constraints,
       rhs = unname(rhs_s[keep]),
       meq = sum(eq_s[keep]))
}
