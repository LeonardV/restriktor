# utility functions

coef.conLM <- function(object, ...)  {
  object$b.constr
}

model.matrix.conLM <- function(object, ...) {
  model.matrix(object$model.org)
}


# is this method fool proof?
dfEq_correction <- function(object, ...) {
  partable <- object
  # if parameters are constraint (==) to a number, df should be adjusted.
  opEq.idx <- !grepl("[^==]", partable$op)
  partable.df <- as.data.frame(partable)
  idx <- partable.df[opEq.idx,]
  # check if lhs contains any numbers
  lhs.idx <- grepl("(^[[:digit:]][.][[:digit:]]$)|(^[0-9]$)", as.vector(idx$lhs))
  # check if rhs contains any numbers
  rhs.idx <- grepl("(^[[:digit:]][.][[:digit:]]$)|(^[0-9]$)", as.vector(idx$rhs))
  # check if perhaps both sides contain numbers. This is probably not possible,
  # but they should not be included.
  check.idx <- length(which(rowSums(cbind(lhs.idx, rhs.idx)) == 2))
  # remove any numbers before checking duplicate variable names. 
  char.idx <- !grepl("(^[[:digit:]][.][[:digit:]]$)|(^[0-9]$)", 
                     c(as.vector(idx$lhs), as.vector(idx$rhs)))
  char <- c(as.vector(idx$lhs), as.vector(idx$rhs))[char.idx]
  diff.idx <- length(char) - length(unique(char))

  pEq.corr <- sum(rhs.idx) + sum(lhs.idx) - check.idx + diff.idx
  
  out <- pEq.corr
  attr(out, "char") <- char

  out
}