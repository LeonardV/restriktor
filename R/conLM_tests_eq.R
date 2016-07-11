con_test_Wald <- function(Sigma, JAC, theta.r) {

    # remove redundant rows from JAC *and* theta.r
    npar <- ncol(JAC)

    JAC.aug <- cbind(JAC, theta.r)
    Q <- qr(t(JAC.aug))
    JAC.full <- t(qr.X(Q, ncol = Q$rank))
    JAC <- JAC.full[,seq_len(npar),drop = FALSE]
    theta.r <- as.numeric(JAC.full[,(npar + 1L)])

    # restricted vcov
    info.r  <- JAC %*% Sigma %*% t(JAC)

    # Wald test statistic
    Wald <- as.numeric(t(theta.r) %*% solve(info.r) %*% theta.r)

    # df
    Wald.df <- nrow(JAC)

    # p-value based on chisq
    Wald.pvalue <- 1 - pchisq(Wald, df = Wald.df)

    OUT <- list(test = "Wald",
                Ts = Wald,
                df = Wald.df,
                pvalue = Wald.pvalue)

    OUT
}



# con_test_score <- function(I, JAC, d0.r) {
#   
#   # score test statistic
#   score <- as.numeric(d0.r %*% solve(I) %*% d0.r)
#   #score <- as.numeric(t(d0.r) %*% solve(info.r) %*% d0.r) 
#   
#   # df
#   score.df <- nrow(JAC)
#   
#   # p-value based on chisq
#   score.pvalue <- 1 - pchisq(score, df = score.df)
#   
#   OUT <- list(test = "score",
#               Ts = score,
#               df = score.df,
#               pvalue = score.pvalue)
#   
#   OUT
# }
