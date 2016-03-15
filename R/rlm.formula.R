# adjusted rlm.formula function for inequality constraints from MASS package
rlm.formula <- function(formula, data, weights, ..., subset, na.action,
                           method = c("M", "MM", "model.frame"),
                           wt.method = c("inv.var", "case"),
                           model = TRUE, x.ret = TRUE, y.ret = FALSE, 
                           contrasts = NULL) {
  mf <- match.call(expand.dots = FALSE)
  mf$method <- mf$wt.method <- mf$model <- mf$x.ret <- mf$y.ret <- 
    mf$contrasts <- mf$Amat <- mf$bvec <- mf$meq <- mf$... <- NULL
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  method <- match.arg(method)
  wt.method <- match.arg(wt.method)
  if(method == "model.frame") return(mf)
  mt <- attr(mf, "terms")
  y <- model.response(mf)
  offset <- model.offset(mf)
  if(!is.null(offset)) y <- y - offset
  x <- model.matrix(mt, mf, contrasts)
  xvars <- as.character(attr(mt, "variables"))[-1L]
  if ((yvar <- attr(mt, "response")) > 0L)
    xvars <- xvars[-yvar]
  xlev <- if (length(xvars) > 0L) {
    xlev <- lapply(mf[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  weights <- model.weights(mf)
  if(!length(weights)) weights <- rep(1, nrow(x))
  fit <- rlm_fit(x, y, weights, method = method,
                    wt.method = wt.method, ...)
  fit$terms <- mt
  ## fix up call to refer to the generic, but leave arg name as `formula'
  cl <- match.call()
  cl[[1L]] <- as.name("rlm")
  fit$call <- cl
  fit$contrasts <- attr(x, "contrasts")
  fit$xlevels <- .getXlevels(mt, mf)
  fit$na.action <- attr(mf, "na.action")
  if(model) fit$model <- mf
  if(!x.ret) fit$x <- NULL
  if(y.ret) fit$y <- y
  fit
}
