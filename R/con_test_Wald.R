con_test_Wald <- function(VCOV, JAC, theta.r) {

    # remove redundant rows from JAC *and* theta.r
    npar <- ncol(JAC)

    JAC.aug <- cbind(JAC, theta.r)
    Q <- qr(t(JAC.aug))
    JAC.full <- t(qr.X(Q, ncol = Q$rank))
    JAC <- JAC.full[,seq_len(npar),drop = FALSE]
    theta.r <- as.numeric(JAC.full[,(npar + 1L)])

    # restricted vcov
    info.r  <- JAC %*% VCOV %*% t(JAC)

    # Wald test statistic
    Wald <- as.numeric(t(theta.r) %*% solve( info.r ) %*% theta.r)

    # df
    Wald.df <- nrow(JAC)

    # p-value based on chisq
    Wald.pvalue <- 1 - pchisq(Wald, df=Wald.df)

    OUT <- list(test = "Wald",
                stat = Wald,
                df = Wald.df,
                p.value = Wald.pvalue)

    OUT
}
