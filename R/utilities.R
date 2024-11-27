## utility functions
clean_constraints <- function(constraints) {
  constraints <- gsub("[#!].*(?=\n)", "", constraints, perl = TRUE) # Remove comments
  constraints <- gsub("[;|&]", "\n", constraints, perl = TRUE) # Standardize delimiters to \n
  constraints <- gsub("[ \t]+", "", constraints, perl = TRUE) # Remove extra whitespace
  constraints <- gsub("\n{2,}", "\n", constraints, perl = TRUE) # Remove extra newlines
  
  # Adjust parentheses formatting
  if (length(gregexpr("[\\(\\)]", constraints)[[1]]) %% 2 == 0) {
    constraints <- gsub("\\),", ")\n", constraints, perl = TRUE)
  } else {
    constraints <- gsub(",", "\n", constraints, perl = TRUE)
  }
  
  return(constraints)
}

# Function to process constraint syntax
process_constraint_syntax <- function(constraint.syntax) {
  syntax <- strsplit(constraint.syntax, split = "\n", perl = TRUE)[[1]]
  syntax <- syntax[syntax != ""] # Remove empty strings
  syntax <- gsub("==", "=", syntax, perl = TRUE) # Normalize equality
  syntax <- gsub("<=", "<", syntax, perl = TRUE) # Normalize less-than-equal
  syntax <- gsub(">=", ">", syntax, perl = TRUE) # Normalize greater-than-equal
  return(syntax)
}


# # function taken from 'bain' package 
# expand_parentheses <- function(hyp) {
#   parenth_locations <- gregexpr("[\\(\\)]", hyp)[[1]]
#   if (!parenth_locations[1] == -1 && !grepl("abs\\(.*\\)", hyp)) {
#     if (length(parenth_locations) %% 2 > 0) {
#       stop("Not all opening parentheses are matched by a closing parenthesis, or vice versa.")
#     }
#     expanded_contents <- strsplit(substring(hyp, (parenth_locations[1]+1), (parenth_locations[2]-1)), ",")[[1]]
#     if (length(parenth_locations) == 2){
#       return(paste0(substring(hyp, 1, (parenth_locations[1]-1)), expanded_contents, substring(hyp, (parenth_locations[2]+1), nchar(hyp))))
#     } else {
#       res <- apply(expand.grid(expanded_contents, expand_parentheses(substring(hyp, (parenth_locations[2]+1), nchar(hyp)))), 1, paste, collapse = "")
#       return(res)
#     }
#   } else {
#     return(hyp)
#   }
# }

# function taken from 'bain' package, but adapted by LV (27-11-2024)
expand_parentheses <- function(hyp) {
  # Vind de locaties van haakjes
  parenth_locations <- gregexpr("[\\(\\)]", hyp)[[1]]
  
  # Als er haakjes zijn
  if (parenth_locations[1] != -1 && !grepl("abs\\(.*\\)", hyp)) {
    # Controleer of de haakjes correct sluiten
    if (length(parenth_locations) %% 2 != 0) {
      stop("Not all opening parentheses are matched by a closing parenthesis, or vice versa.")
    }
    
    # Haal de inhoud binnen de eerste set haakjes
    expanded_contents <- strsplit(
      substring(hyp, parenth_locations[1] + 1, parenth_locations[2] - 1), 
      ",\\s*"
    )[[1]]
    
    # Maak een lijst van uitgebreide strings
    expanded_strings <- sapply(expanded_contents, function(content) {
      paste0(
        substring(hyp, 1, parenth_locations[1] - 1), 
        content, 
        substring(hyp, parenth_locations[2] + 1, nchar(hyp))
      )
    }, USE.NAMES = FALSE)
    
    # Als er nog meer haakjes zijn, verwerk ze recursief
    if (any(grepl("\\(", expanded_strings))) {
      return(unlist(lapply(expanded_strings, expand_parentheses)))
    } else {
      return(expanded_strings)
    }
  } else {
    # Geen haakjes aanwezig
    return(hyp)
  }
}


# function taken from 'bain' package 
# expand_compound_constraints <- function(hyp) {
#   equality_operators <- gregexpr("[=<>]", hyp)[[1]]
#   if(length(equality_operators) > 1){
#     string_positions <- c(0, equality_operators, nchar(hyp)+1)
#     res <- sapply(1:(length(string_positions)-2), function(pos) {
#       substring(hyp, (string_positions[pos]+1), (string_positions[pos+2]-1))
#     })
#     return(res)
#   } else {
#     return(hyp)
#   }
# }

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
