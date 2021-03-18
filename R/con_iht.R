iht <- function(...) {
  ldots <- list(...)
  if ("model"%in% names(ldots) | is.character(ldots[[1]][1])) {
    conTestD(...)
  } else {
    conTest(...)
  }
}

