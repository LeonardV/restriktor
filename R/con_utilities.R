# utility functions

coef.conLM <- function(object, ...)  {
  object$b.constr
}

model.matrix.conLM <- function(object, ...) {
  model.matrix(object$model.org)
}




# # # is this method fool proof?
# conEqualZero <- function(object,
#                          bvec.idx = "(^[[:digit:]][.][[:digit:]]$)|(^[0-9]$)",
#                          ...) {
#   parTable <- object
#   opEq.idx <- !grepl("[^==]", parTable$op)
#   parTable.df <- as.data.frame(parTable)
#   idx <- parTable.df[opEq.idx,]
#   # check if lhs contains any numbers
#   lhs.idx <- grepl(bvec.idx, as.vector(idx$lhs))
#   # check if rhs contains any numbers
#   rhs.idx <- grepl(bvec.idx, as.vector(idx$rhs))
#   # remove bvec.idx before checking duplicate variable names.
#   char.idx <- !grepl(bvec.idx, c(as.vector(idx$lhs), as.vector(idx$rhs)))
#   char <- c(as.vector(lhs.idx), as.vector(rhs.idx))[char.idx]
#   char <- unique(char)
#   #diff.idx <- length(char) - length(unique(char))
#   #pEq.corr <- sum(rhs.idx) + sum(lhs.idx) - check.idx + diff.idx
# 
#   out <- char
#   
#     return(out)
# }