con_gorica_est <- function(object, constraints = NULL, VCOV = NULL,
                           B = 999, rhs = NULL, neq = 0L, mix.weights = "pmvnorm", 
                           mix.bootstrap = 99999L, parallel = "no", ncpus = 1L, cl = NULL, 
                           seed = NULL, control = list(), verbose = FALSE, 
                           debug = FALSE, ...) {
  
  # check class
  if (!(class(object)[1] == "numeric")) {
    stop("Restriktor ERROR: object must be of class numeric.")
  }
  
  if (is.null(VCOV)) {
    stop("Restriktor ERROR: variance-covariance matrix VCOV must be provided.")
  }
  
  # check method to compute chi-square-bar weights
  if (!(mix.weights %in% c("pmvnorm", "boot", "none"))) {
    stop("Restriktor ERROR: ", sQuote(mix.weights), " method unknown. Choose from \"pmvnorm\", \"boot\", or \"none\"")
  }
  
  # timing
  start.time0 <- start.time <- proc.time()[3]; timing <- list()
  # store call
  mc <- match.call()
  # rename for internal use
  Amat <- constraints
  bvec <- rhs 
  meq  <- neq
  
  b.unrestr <- object
  b.unrestr[abs(b.unrestr) < ifelse(is.null(control$tol), 
                                    sqrt(.Machine$double.eps), 
                                    control$tol)] <- 0L
  # ML unconstrained MSE
  Sigma <- VCOV
  # number of parameters
  p <- length(b.unrestr)
  ll.unrestr <- dmvnorm(rep(0, p), sigma = Sigma, log = TRUE)
  
  if (debug) {
    print(list(loglik.unc = ll.unrestr))
  }
  
  timing$preparation <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # deal with constraints
  if (!is.null(constraints)) {
    restr.OUT <- con_constraints(object, 
                                 constraints = Amat, 
                                 bvec        = bvec, 
                                 meq         = meq, 
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
  
  # compute the reduced row-echelon form of the constraints matrix
  rAmat <- GaussianElimination(t(Amat))
  if (mix.weights == "pmvnorm") {
    if (rAmat$rank < nrow(Amat) && rAmat$rank != 0L) {
      stop(paste("Restriktor ERROR: The constraint matrix must have full row-rank", 
                 "\n  (choose e.g. rows", paste(rAmat$pivot, collapse = " "),")."))
    }
  } 
  
  timing$constraints <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  if (ncol(Amat) != length(b.unrestr)) {
    stop("Restriktor ERROR: length coefficients and the number of",
         "\n       columns constraints-matrix must be identical")
  }
  
  if (!(nrow(Amat) == length(bvec))) {
    stop("nrow(Amat) != length(bvec)")
  }
  
  # compute chi-square-bar weights
  if (mix.weights != "none") {
    if (nrow(Amat) == meq) {
      # equality constraints only
      wt.bar <- rep(0L, ncol(Sigma) + 1)
      wt.bar.idx <- ncol(Sigma) - meq + 1
      wt.bar[wt.bar.idx] <- 1
    } else if (all(c(Amat) == 0)) { 
      # unrestricted case
      wt.bar <- c(rep(0L, p), 1)
    } else if (mix.weights == "boot") { 
      # compute chi-square-bar weights based on Monte Carlo simulation
      wt.bar <- con_weights_boot(VCOV     = Sigma,
                                 Amat     = Amat, 
                                 meq      = meq, 
                                 R        = mix.bootstrap,
                                 parallel = parallel,
                                 ncpus    = ncpus,
                                 cl       = cl,
                                 seed     = seed,
                                 verbose  = verbose)
      attr(wt.bar, "mix.bootstrap") <- mix.bootstrap 
    } else if (mix.weights == "pmvnorm" && (meq < nrow(Amat))) {
      # compute chi-square-bar weights based on pmvnorm
      wt.bar <- rev(con_weights(Amat %*% Sigma %*% t(Amat), meq = meq))
    } 
  } else {
    wt.bar <- NA
  }
  attr(wt.bar, "method") <- mix.weights
  
  if (debug) {
    print(list(mix.weights = wt.bar))
  }
  
  timing$mix.weights <- (proc.time()[3] - start.time)
  start.time <- proc.time()[3]
  
  # check if the constraints are not in line with the data, else skip optimization
  if (all(Amat %*% c(b.unrestr) - bvec >= 0 * bvec) & meq == 0) {
    b.restr  <- b.unrestr
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.unrestr,
                loglik      = ll.unrestr, 
                Sigma       = Sigma,
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                wt.bar      = wt.bar,
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
    b.restr[abs(b.restr) < ifelse(is.null(control$tol), 
                                  sqrt(.Machine$double.eps), 
                                  control$tol)] <- 0L
    
    timing$optim <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    
    ll.restr <- dmvnorm(c(b.unrestr - b.restr), sigma = Sigma, log = TRUE)
    
    OUT <- list(CON         = CON,
                call        = mc,
                timing      = timing,
                parTable    = parTable,
                b.unrestr   = b.unrestr,
                b.restr     = b.restr,
                loglik      = ll.restr, 
                Sigma       = Sigma,
                constraints = Amat, 
                rhs         = bvec, 
                neq         = meq, 
                wt.bar      = wt.bar,
                iact        = out.solver$iact, 
                control     = control)
  }
  
  
  OUT$timing$total <- (proc.time()[3] - start.time0)
  
  class(OUT) <- c("gorica_est")
  
  OUT
}