# library(mvtnorm)
# library(restriktor)
# n <- 100
# p <- 4
# betas <- c(0.1,0.2,0.3,0.4,0.5)
# set.seed(3013073)
# X <- cbind(rmvnorm(n, mean = rep(0,p), sigma = diag(p)), rbinom(n,1,0.5))
# colnames(X) <- c("x1","x2","x3","x4","f1")
# z <- X %*% betas        
# y <- z + rnorm(n)
# DATA <- data.frame(y, X)
# 
# # intercept model
# model1 <- y ~  1 + x1 + x2 + x3 + x4
# 
# ############################ lm #################################
# linmod1 <- lm(model1, data = DATA)
# constraints1 <- 'x1 > 0; x2 > 0; x3 > 0'
# 
# fit1.goric <- goric(linmod1, hypotheses = list(constraints1), 
#                     comparison = "complement")
# 
# goric_obj <- fit1.goric
# pop.est <- rep(0, length(coef(goric_obj)))
# 
# # object = goric(a)
# # object = est + vcov
# 
# object <- fit1.goric #list(fit1.goric, fit2.goric)
# 
# 
# benchmarks_goric(fit1.goric, iters = 100)


# -------------------------------------------------------------------------
#library(benchmarks)

# -------------------------------------------------------------------------
# benchmarks_goric(x, iters = 10, verbose = TRUE)
out <- benchmarks_goric(fit1.goric, iters = 100, verbose = TRUE, seed = 1)
out

out <- benchmarks_goric(fit1.goric, pop.est = rbind(rep(0, 5)), iters = 100, 
                        verbose = TRUE, seed = 1)
out

pop_est <- rbind(rep(0, 5), c(0,.1,.2,.3,.4))
out <- benchmarks_goric(fit1.goric, pop.est = pop_est, iters = 500, 
                        verbose = TRUE)
out

pop_est <- rbind(rep(0, 5), c(0,1,2,3,4))
out <- benchmarks_goric(fit1.goric, pop.est = pop_est, iters = 50, 
                        verbose = TRUE)
out

#head(out$ratio_gorica_weights)
hist(out$gw_pref_hypo, 50)
abline(v=fit1.goric$result$goric.weights[2])

# out$ratio_gorica_weights
# out
# benchmarks_goric(fit2.goric, iters = 1, verbose = TRUE, sample.nobs = 50,
#                  mix.weights = "boot", mix.bootstrap = 99999, parallel = "multicore",
#                  ncpus = 1)


# TODO
# add n to output goric
# add parallel processing
# if only 1 popo.est specified, repeat it internaly for all params.
# add plots

benchmarks_goric <- function(object, pop.est = NULL, sample.nobs = NULL, 
                             quantiles = NULL, 
                             iters = 1000, seed = 90210, verbose = TRUE, ...) {
  
  x <- object
  # text-based hypothesis: 'x1 > x2 > x3'
  hypos <- x$hypotheses_usr
  # number of hypotheses incl. complement or unconstrained if available.
  nr_hypos <- length(x$result$model)
  # goric weights
  # note: cannot use colname because it changes goric(a).weights
  gw <- x$result[, 7]
  hypo_names <- x$result$model
  names(gw) <- hypo_names #Map(`names<-`, gw, hypo_names)
  # which goric weight is highest
  pref_hypo <- which.max(gw)
  # preferred goric weight 
  gw_pref_hypo_org <- x$result[pref_hypo, 7, drop = FALSE]
  rownames(gw_pref_hypo_org) <- names(pref_hypo)
  # ratio goric weight for preferred hypothesis
  rgw_pref_hypo_org <- x$ratio.gw[pref_hypo, -pref_hypo, drop = FALSE] 
  rownames(rgw_pref_hypo_org) <- names(pref_hypo)
  # number of parameters
  n_coef <- length(coef(x))
  
  # get unconstrained parameter estimates
  if (is.null(pop.est)) {
    # may be replace with null-hypothesis = 0
    pop_est <- lapply(x$objectList, function(y) y$b.unrestr)[[1]]
    names(pop_est ) <- NULL
    pop_est <- t(as.matrix(unlist(pop_est)))
  } else if (is.vector(pop.est)) {
    pop_est <- matrix(pop.est, nrow = 1, ncol = n_coef)
  } else {
    pop_est <- pop.est
  }
  
  nr_es <- nrow(pop_est)
  
  # unconstrained VCOV
  VCOV <- x$Sigma
  b.unrestr <- x$b.unrestr
  vnames <- names(b.unrestr)
  
  # sample size
  samplesize <- x$objectList[[1]]$df.residual 
  
  if (is.null(samplesize) && is.null(sample.nobs)) {
    stop("Restriktor ERROR: specify the sample-size using the 'sample.nobs' argument.", call. = FALSE)
  } 
  
  if (is.null(samplesize)) {
    samplesize <- sample.nobs  
  }
  
  if (!is.null(sample.nobs)) {
    VCOV <- VCOV * samplesize / sample.nobs
    samplesize <- sample.nobs
  }
  
  set.seed(seed)
  #quant <- c(.05, .35, .50, .65, .95)
  quant <- c(.01, .05, .25, .50, .75, .95, .99)
  names_quant <- c("1%", "5%", "25%", "50%", "75%", "95%", "99%")
  
  pvalue <- CI_benchmarks_rgw_all <- CI_benchmarks_all <- list()
  for (teller.es in 1:nr_es) {
    name <- paste0("nr.pop.est = ", teller.es)
    est <- MASS:::mvrnorm(n = iters, mu = pop_est[teller.es, ], Sigma = VCOV)
    est <- matrix(est, nrow = iters, ncol = length(pop_est[teller.es, ]))
    colnames(est) <- vnames
    
    # goric weights for the prefered hypothesis, i.e., max goric weight
    gw_pref_hypo  <- rep(NA, iters)
    rgw_pref_hypo <- matrix(NA, nrow = iters, ncol = nr_hypos-1)
    colnames(rgw_pref_hypo) <- paste(names(pref_hypo), names(x$ratio.gw[pref_hypo, -pref_hypo]))
    
    for (i in 1:iters) {
      if (verbose) cat("...iteration:", i, "\n")
      pop.est.CI <- est[i, ]
      results_gorica <- goric(pop.est.CI, VCOV = VCOV, hypotheses = hypos, 
                              comparison = x$comparison, type = "gorica")
      
      gw_pref_hypo[i] <- results_gorica$result[pref_hypo, 7]
      rgw_pref_hypo[i, ] <- results_gorica$ratio.gw[pref_hypo, -pref_hypo, drop = FALSE]
    }
    
    # how many samples are lager than the observed preferred goric weight
    pvalue[[name]] <- sum(gw_pref_hypo > gw_pref_hypo_org[1,1]) / iters
    
    CI_benchmarks_goric <- matrix(quantile(gw_pref_hypo, quant), nrow = 1) 
    
    colnames(CI_benchmarks_goric) <- names_quant
    rownames(CI_benchmarks_goric) <- names(pref_hypo)
    
    # -1 = we dont want the pref_hypo in the output
    CI_benchmarks_rgw <- matrix(NA, nrow = nr_hypos-1, ncol = length(quant))
    # select the best hypothesis from the example data
    #CI_benchmarks_rgw[, 1] <- x$ratio.gw[pref_hypo, -pref_hypo] 
    CI_benchmarks_rgw[, 1:length(quant)] <- t(apply(rgw_pref_hypo, 2, function(x) quantile(x, quant)))
      colnames(CI_benchmarks_rgw) <- names_quant
      rownames(CI_benchmarks_rgw) <- paste(names(pref_hypo), colnames(x$ratio.gw[pref_hypo, -pref_hypo, drop = FALSE]))
    
    CI_benchmarks_all[[name]] <- CI_benchmarks_goric
    CI_benchmarks_rgw_all[[name]] <- CI_benchmarks_rgw
  }
  
  # Error probability based on complement of preferred hypothesis in data
  # TO DO zie ANOVA fie
  if (nr_hypos == 2 && x$comparison == "complement") {
    error_prob <- 1 - x$result[pref_hypo, 7]
  } else {
    if (pref_hypo == nr_hypos) {
      error_prob <- "The unconstrained hypothesis is preferred..."
    } else {
      if ("gorica_est" %in% class(x$objectList[[1]])) {
        object_goric <- b.unrestr
        VCOV2 <- VCOV
      } else {
        H_pref <- hypos[[pref_hypo]]
        object_goric <- x$model.org
        VCOV2 <- NULL
      }
      
      results_goric_pref <- goric(object_goric, VCOV = VCOV2, 
                                  hypotheses = list(H_pref = H_pref), 
                                  comparison = "complement", type = x$type, ...)
      
      error_prob <- results_goric_pref$result[2, 7]
    }
  }
  
  
  out <- list(gw_pref_hypo = gw_pref_hypo,
              rgw_pref_hypo = rgw_pref_hypo,
              iters = iters,
              n_coef = n_coef, N = samplesize,
              pop_estim = pop_est, pop_VCOV = VCOV,
              pref_hypo = pref_hypo,
              gw_pref_hypo_org = gw_pref_hypo_org,
              rgw_pref_hypo_org = rgw_pref_hypo_org,
              error_prob_pref_hypo = error_prob,
              pvalue = pvalue,
              benchmarks_gw  = CI_benchmarks_all,
              benchmarks_gwr = CI_benchmarks_rgw_all)
  
  
  class(out) <- c("con_goric_benchmarks")
  
  out
  
} 

