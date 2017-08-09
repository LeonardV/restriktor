# build a bare-bones parameter table from a fitted object
con_partable <- function(object, est = FALSE, label = FALSE,
                         as.data.frame. = FALSE) {

    #  we first check the class of object
    if (!(class(object)[1] %in% c("lm", "rlm", "glm", "mlm"))) {
       stop("It only works for lm, mlm, rlm and glm")
     }

    objectTerms <- terms(object)

    responseIndex <- attr(objectTerms, "response")
    varNames <- as.character(attr(objectTerms, "variables"))[-1]
    responseName <- varNames[responseIndex]

    predCoef  <- coef(object)
    if (class(object)[1] == "mlm") {
      predNames <- colnames(predCoef)  
    } else {
      predNames <- names(predCoef)
    }
    
    lhs <- rep(responseName, length(predNames))
     op <- rep("~", length(predNames))
    rhs <- predNames

    # intercept?
    if (attr(objectTerms, "intercept")) {
        int.idx <- which(rhs == "(Intercept)")
        op[int.idx] <- "~1"
        rhs[int.idx] <- ""
    }

    # always add residual variance?
    #lhs <- c(lhs, responseName)
    # op <- c(op, "~~")
    #rhs <- c(rhs, responseName)

    # construct minimal partable
    partable <- list(lhs = lhs, op = op, rhs = rhs)

    # include 'est' column?
    if (est) {
        #partable$est <- c(as.numeric(predCoef),
        #                  sum(resid(object)^2) / object$df.residual)
        partable$est <- as.numeric(predCoef)
    }

    # include 'label' column?
    if (label) {
        # partable$label <- c(predNames, responseName)
        partable$label <- predNames

        # convert all ':' to '.'
        partable$label <- gsub("[:()]", ".", partable$label)
    }

    # convert to data.frame?
    if (as.data.frame.) {
        partable <- as.data.frame(partable, stringsAsFactors = FALSE)
    }
    
    partable
}
