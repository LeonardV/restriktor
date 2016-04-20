# utility functions
coef.conLM <- function(object, ...)  {
  object$b.constr
}


model.matrix.conLM <- function(object, ...) {
  model.matrix(object$model.org)
}

