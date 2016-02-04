# acknowledgement: code taken from ic.infer package
conLM.default <- function(model, constraints, se = "default", 
                          bvec = NULL, meq = 0L, control = NULL,
                          tol = sqrt(.Machine$double.eps), df.error = NULL, 
                          debug = FALSE, ...) {

    ## check model
    if (!(is.matrix(model))) {
        stop("ERROR: model must be of class lm or a covariance matrix.")
    }  
    else if (!(nrow(model)==ncol(model))) {
        stop("ERROR: If it is not a linear model, model must be a quadratic matrix.")
    }  
    else if (!(all(eigen(model,TRUE,only.values=TRUE)$values>0))) {
        stop("ERROR: matrix model must be positive definite.")
    }
    g <- nrow(model)-1
    if (is.null(df.error)) {
        stop("ERROR: df.error is required, when working from a covariance matrix.")
    }
    if (!(df.error > 2)) {
        stop("ERROR: df.error must be at least 2.")
    }

    Amat <- constraints
    ## preliminary calculations
    V <- solve(model[2:(g+1),2:(g+1)]) / df.error
    if (is.null(colnames(V))) {
      colnames(V) <- paste("X",1:g,sep="")
    }
    b <- solve(model[2:(g+1),2:(g+1)],model[2:(g+1),1])
    var.y <- model[1,1]
    s2 <- c(var.y - model[1,2:(g+1)]%*%b)
    b <- as.vector(b)
    names(b) <- colnames(V)
    orig.R2 <- model[1,2:(g+1)]%*%b/var.y
    ## check inputs
    if (is.vector(Amat)) {
        Amat <- matrix(Amat, 1, length(Amat))
    }
    if (!is.matrix(Amat)) {
        stop("Amat must be a matrix.")
    }
    if (is.null(bvec)) {
        bvec <- rep(0, nrow(Amat))
    }
    if (!is.vector(bvec)) {
        stop("bvec must be a vector.")
    }
    if (!nrow(Amat) == length(bvec)) {
        stop("mismatch between number of rows in Amat and elements in bvec")
    }
    ## inequality constraints only, all fulfilled
    if (all(Amat %*% c(b) - bvec >= 0 * bvec) & meq == 0) {
        out <- list(b.unconstr = b, b.constr = b,
                    R2 = orig.R2, residuals = NULL, fitted.values = NULL,
                    weights = NULL, orig.R2 = orig.R2,
                    df.error = df.error, s2 = s2, Sigma = s2*V,
                    origmodel = NULL, Amat = Amat, bvec = bvec, iact = NULL,
                    meq = meq, bootout = NULL)
    }
    else {
        ## equality constraints involved or some inequality constraints violated
        ## calculate restricted estimate
        out <- con_my_solve_QP_lm(Dmat = solve(V), dvec = solve(V, b),
                                  Amat = t(Amat), bvec = bvec, meq = meq)
        names(out$solution) <- names(b)
        out$solution[abs(out$solution) < tol] <- 0
        ## initialize output list
        out <- list(b.constr = out$solution, b.unconstr = b,
                    R2 = NULL, residuals = NULL, fitted.values = NULL,
                    weights = NULL, orig.R2 = orig.R2,
                    df.error = df.error, s2 = s2, Sigma = s2*V,
                    origmodel = NULL, Amat = Amat, bvec = bvec, iact = out$iact,
                    meq = meq, bootout = NULL)
        ### R2
        out$R2 <- model[1,2:(g+1)]%*%t(t(out$b.constr))/var.y
    }
    class(out) <- "conLM"
    out
}
