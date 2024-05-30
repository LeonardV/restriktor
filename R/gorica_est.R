con_gorica_est <- function(object, constraints = NULL, VCOV = NULL,
                           rhs = NULL, neq = 0L, mix_weights = "pmvnorm", 
                           seed = NULL, control = list(), verbose = FALSE, 
                           debug = FALSE, ...) {
  
  if (is.null(VCOV)) {
    stop("Restriktor ERROR: variance-covariance matrix VCOV must be provided.")
  }
  # check method to compute chi-square-bar weights
  if (!(mix_weights %in% c("pmvnorm", "boot", "none"))) {
    stop("Restriktor ERROR: ", sQuote(mix_weights), " method unknow. Choose from \"pmvnorm\", \"boot\", or \"none\"")
  }
  
  # timing
  start.time0 <- start.time <- proc.time()[3]; timing <- list()
  # store call
  #mc <- match.call()
  # rename for internal use
  Amat <- constraints
  bvec <- rhs 
  meq  <- neq
  
  b.unrestr <- object
  b.unrestr[abs(b.unrestr) < ifelse(is.null(control$tol), 
                                    sqrt(.Machine$double.eps), 
                                    control$tol)] <- 0L
  Sigma <- VCOV
  # number of parameters
  p <- length(b.unrestr)
  # unrestricted log-likelihood
  ll.unrestr <- dmvnorm(rep(0, p), sigma = Sigma, log = TRUE)
  
  if (debug) {
    print(list(loglik.unc = ll.unrestr))
  }
  
  timing$preparation <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # deal with constraints
  if (!is.null(constraints)) {
    restr.OUT <- con_constraints(object, 
                                 VCOV        = Sigma,
                                 est         = b.unrestr,
                                 constraints = Amat, 
                                 bvec        = bvec, 
                                 meq         = meq, 
                                 mix_weights = mix_weights,
                                 se          = "none",
                                 debug       = debug)  
    # a list with useful information about the restriktions.}
    CON <- restr.OUT$CON
    # a parameter table with information about the observed variables in the object 
    # and the imposed restriktions.}
    parTable <- restr.OUT$parTable
    # constraints matrix
    Amat <- restr.OUT$Amat
    # rhs 
    bvec <- restr.OUT$bvec
    # neq
    meq <- restr.OUT$meq
  } else if (is.null(constraints)) { 
    # no constraints specified - needed for GORIC to include unconstrained model
    CON <- NULL
    parTable <- NULL
    Amat <- rbind(rep(0L, p))
    bvec <- rep(0L, nrow(Amat))
    meq  <- 0L
  } 
  
  # if only new parameters are defined and no constraints
  if (length(Amat) == 0L) {
    Amat <- rbind(rep(0L, p))
    bvec <- rep(0L, nrow(Amat))
    meq  <- 0L
  }
  
  ## create list for warning messages
  messages <- list()
  
  timing$constraints <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  if (ncol(Amat) != length(b.unrestr)) {
    stop("Restriktor ERROR: length coefficients and the number of",
         "\n       columns constraints-matrix must be identical")
  }
  
  if (!(nrow(Amat) == length(bvec))) {
    stop("nrow(Amat) != length(bvec)")
  }
  
  start.time <- proc.time()[3]
  
  # check if the constraints are not in line with the data, else skip optimization
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) && meq == 0) {
    b.restr  <- b.unrestr
    
    OUT <- list(CON         = CON,
                #call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.unrestr,
                loglik      = ll.unrestr, 
                Sigma       = Sigma,
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                wt.bar      = NULL,
                iact        = 0L, 
                control     = control)  
  } else {
    # compute constrained estimates using quadprog
    out.solver <- con_solver_gorica(est  = b.unrestr, 
                                    VCOV = Sigma, 
                                    Amat = Amat, 
                                    bvec = bvec, 
                                    meq  = meq)
    b.restr <- out.solver$solution
    names(b.restr) <- names(b.unrestr)
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    
    timing$optim <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    ll.restr <- dmvnorm(c(b.unrestr - b.restr), sigma = Sigma, log = TRUE)
    
    OUT <- list(CON         = CON,
                #call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.restr,
                loglik      = ll.restr, 
                Sigma       = Sigma,
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                wt.bar      = NULL,
                iact        = out.solver$iact, 
                control     = control)
  }

  
  Amat_meq_PT <- PT_Amat_meq(Amat, meq)
  RREF <- Amat_meq_PT$RREF
  PT_Amat <- Amat_meq_PT$PT_Amat
  PT_meq  <- Amat_meq_PT$PT_meq
  OUT$PT_meq  <- PT_meq
  OUT$PT_Amat <- PT_Amat
  
  if (mix_weights == "pmvnorm") {
    if (RREF$rank < nrow(PT_Amat) && RREF$rank != 0L) {
      messages$mix_weights_rank <- paste(
        "Restriktor message: Since the constraint matrix is not full row-rank, the level probabilities", 
        "are calculated using mix_weights = \"boot\" (the default is mix_weights = \"pmvnorm\").",
        "For more information see ?restriktor.\n"
      )
      mix_weights <- "boot"
    }
  } 

  ## determine level probabilities
  wt.bar <- calculate_weight_bar(Amat = PT_Amat, meq = PT_meq, VCOV = Sigma, 
                                   mix_weights = mix_weights, seed = seed, 
                                   control = control, verbose = verbose, ...)
  attr(wt.bar, "method") <- mix_weights
  OUT$wt.bar <- wt.bar
  
  if (debug) {
    print(list(mix_weights = wt.bar))
  }
  
  timing$mix_weights <- (proc.time()[3] - start.time)
  OUT$messages <- messages
  OUT$timing$total <- (proc.time()[3] - start.time0)
  
  class(OUT) <- c("gorica_est")
  
  OUT
}



con_gorica_est_lav <- function(x, standardized = FALSE, ...) {
  ## create empty list
  out <- list()
  ## number of groups 
  #num_groups <- lavInspect(x, what = "ngroups")
  ## get parameter table
  unstandardized_parTable <- parTable(x)
  unstandardized_parTable <- unstandardized_parTable[unstandardized_parTable[, "plabel"] != "", ]
  standardized_parTable   <- standardizedSolution(x, ci = FALSE, zstat = FALSE, se = FALSE)$est.std
  
  ## combine unstandardized and standardized parameter estimates  
  parameter_table <- cbind(unstandardized_parTable, est.std = standardized_parTable)
  
  ## Only user-specified labels
  parameter_table <- parameter_table[parameter_table$label != "" & parameter_table$free != 0L, ]
  ## remove any duplicate labels
  parameter_table <- parameter_table[!duplicated(parameter_table$label), ]
  ## use (un)standardized parameter estimates
  out$estimate <- 
    if (standardized) {
      parameter_table$est.std
    } else { 
      parameter_table$est
    }
  names(out$estimate) <- parameter_table$label
  ## extract (un)standardized VCOV
  out$VCOV <- 
    if (standardized) {
      lavInspect(x, "vcov.std.all")
    } else {
      lavInspect(x, "vcov")
    }
  ## remove not used columns of VCOV
  out$VCOV <- out$VCOV[parameter_table$label, parameter_table$label, drop = FALSE]
  
  out
}
