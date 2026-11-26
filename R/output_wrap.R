# Line breaks in output text -------------------------------------------------
#
# Much of the output of restriktor (e.g., conclusions, notes, messages, and
# warnings) is running text. Without line breaks, only part of it is visible
# in, for example, a pdf made from an Rmd file. Therefore, the text that
# restriktor writes via message(), warning(), stop(), and cat() is wrapped at
# getOption("restriktor.linewidth", 76) characters. Use
# options(restriktor.linewidth = Inf) to switch this off.
#
# Only lines of running text are wrapped; lines that look like (part of) a
# table (i.e., with two or more consecutive spaces after the indentation,
# used for aligning columns), lines with ANSI color codes, and lines without
# spaces are left as they are.
#
# The functions below mask base::message(), base::warning(), base::stop(),
# and base::cat() within the restriktor namespace only (so, only for the code
# of restriktor itself, not for the user or for other packages).

.rk_linewidth <- function() {
  w <- getOption("restriktor.linewidth", 76)
  if (!is.numeric(w) || length(w) != 1L || is.na(w) || w < 20) w <- 76
  w
}

# Wrap the long lines of running text in 'txt' (a single string, possibly
# with line breaks).
.rk_wrap <- function(txt, width = .rk_linewidth()) {
  if (!is.finite(width) || length(txt) != 1L || is.na(txt) || nchar(txt) <= width) {
    return(txt)
  }
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  # strsplit() drops a trailing empty line: restore the trailing line break(s)
  n_trailing <- nchar(txt) - nchar(sub("\n+$", "", txt))
  wrapped <- vapply(lines, function(line) {
    if (nchar(line) <= width) return(line)
    lead <- nchar(line) - nchar(sub("^ +", "", line))
    body <- substring(line, lead + 1L)
    # table-like lines, lines with color codes, and lines without spaces: as is
    if (grepl("  ", body, fixed = TRUE) || grepl("\033", line, fixed = TRUE) ||
        !grepl(" ", body, fixed = TRUE)) {
      return(line)
    }
    # continuation lines of a bullet ('- ', '* ', or '1. ') are indented
    # such that they line up with the text of the bullet
    bullet <- regmatches(body, regexpr("^([-*]|[0-9]+\\.) ", body))
    ex <- lead + if (length(bullet) == 1L) nchar(bullet) else 0L
    paste(strwrap(body, width = width, indent = lead, exdent = ex), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
  out <- paste(wrapped, collapse = "\n")
  if (n_trailing > 0L && !grepl("\n$", out)) {
    out <- paste0(out, strrep("\n", n_trailing))
  }
  out
}

message <- function(..., domain = NULL, appendLF = TRUE) {
  args <- list(...)
  if (length(args) == 1L && inherits(args[[1L]], "condition")) {
    return(base::message(args[[1L]]))
  }
  msg <- .makeMessage(..., domain = domain)
  base::message(.rk_wrap(msg), appendLF = appendLF)
}

warning <- function(..., call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
                    domain = NULL) {
  args <- list(...)
  if (length(args) == 1L && inherits(args[[1L]], "condition")) {
    return(base::warning(args[[1L]]))
  }
  msg <- .rk_wrap(.makeMessage(..., domain = domain))
  # the call of the function that gave the warning (not of this function)
  cl <- if (isTRUE(call.)) sys.call(-1L) else NULL
  if (isTRUE(immediate.)) {
    old <- options(warn = 1L)
    on.exit(options(old))
  }
  invisible(base::warning(simpleWarning(msg, call = cl)))
}

stop <- function(..., call. = TRUE, domain = NULL) {
  args <- list(...)
  if (length(args) == 1L && inherits(args[[1L]], "condition")) {
    base::stop(args[[1L]])
  }
  msg <- .rk_wrap(.makeMessage(..., domain = domain))
  # the call of the function that gave the error (not of this function)
  cl <- if (isTRUE(call.)) sys.call(-1L) else NULL
  base::stop(simpleError(msg, call = cl))
}

cat <- function(..., file = "", sep = " ", fill = FALSE, labels = NULL,
                append = FALSE) {
  args <- list(...)
  # only output to the console (or a sink), and only plain (character or
  # numeric) arguments are wrapped; otherwise, base::cat() as is
  plain <- all(vapply(args, function(a) is.null(a) || (is.atomic(a) && !is.factor(a)),
                      logical(1)))
  if (!identical(file, "") || !isFALSE(fill) || !is.null(labels) || !plain ||
      !is.finite(.rk_linewidth())) {
    return(base::cat(..., file = file, sep = sep, fill = fill, labels = labels,
                     append = append))
  }
  args <- Filter(function(a) length(a) > 0L, args)
  if (length(args) == 0L) return(invisible(NULL))
  pieces <- unlist(lapply(args, as.character), use.names = FALSE)
  # as base::cat(): the elements are separated by 'sep' (recycled)
  sep <- rep_len(sep, length(pieces))
  txt <- paste0(paste0(pieces[-length(pieces)], sep[-length(pieces)], collapse = ""),
                pieces[length(pieces)])
  base::cat(.rk_wrap(txt), sep = "")
}
