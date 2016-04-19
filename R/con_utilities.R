# utility functions
coef.conLM <- function(object, ...)  {
  object$b.constr
}

coef.conRLM <- function(object, ...)  {
  object$b.constr
}

model.matrix.conLM <- function(object, ...) {
  model.matrix(object$model.org)
}

model.matrix.conRLM <- function(object, ...) {
  model.matrix(object$model.org)
}
