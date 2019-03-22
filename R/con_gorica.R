gorica_ll <- function(est, VCOV, Amat, bvec, meq, ...) {
  
  out.qp <- con_solver_gorica(est  = est, 
                              VCOV = VCOV, 
                              Amat = Amat, 
                              bvec = bvec, 
                              meq  = meq)
  
  b.restr   <- out.qp$solution
  b.unrestr <- out.qp$unconstrained.solution
  
  ll <- dmvnorm(c(b.unrestr - b.restr), sigma = VCOV, log = TRUE)
  
  out <- list(ll = ll, b.restr = b.restr, b.unrestr = b.unrestr, 
              Amat = Amat, bvec = bvec, meq = meq)
  
  out
}