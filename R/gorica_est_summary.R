summary.gorica_est <- function(object, type, ...) {
  z <- object
  
  ldots <- list(...)

  if (is.null(ldots$penalty_factor)) {
    ldots$penalty_factor <- 2
  }

  if (!inherits(z, "gorica_est")) {
    stop("object of class ", sQuote(class(z)), " is not supported.")
  }
  
  Amat    <- z$PT_Amat
  meq     <- z$PT_meq
  wt.bar  <- z$wt.bar
  
  ans     <- z$model.org
  
  # number of nuisance parameters (d; parameters of the model that are not in
  # the estimate vector) and of response variables (q > 1 for estimates of an
  # mlm), see penalty_goric() (goric_penalty.R). Only the GORICAC depends on
  # these: the GORICA penalty is the GORIC penalty minus d.
  n_nuis <- ldots$nr_nuis_param
  n_resp <- ldots$n_resp
  if (is.null(n_resp)) n_resp <- 1L
  if (is.null(n_nuis)) {
    if (type == "goricac") {
      stop("\nrestriktor ERROR: the argument 'nr_nuis_param' (number of nuisance ",
           "parameters) is needed when type = 'goricac'.", call. = FALSE)
    }
    n_nuis <- n_resp * (n_resp + 1) / 2
  }
  ## compute goric
  if (!(attr(wt.bar, "method") == "none")) {  
    ## REF: Kuiper, R.M.; Hoijtink, H.J.A.; Silvapulle, M. J. (2012) 
    ## Journal of statistical planning and inference, volume 142, pp. 2454 - 2463
    if (type %in% c("goric", "gorica")) {
      PT <- penalty_goric(Amat        = Amat, 
                          meq         = meq, 
                          LP          = wt.bar, 
                          correction  = FALSE, 
                          sample.nobs = NULL,
                          n_resp      = n_resp,
                          n_nuis      = n_nuis)
      if (type == "gorica") {
        PT <- PT - n_nuis
      }
    } else if (type %in% c("goricc", "goricac")) {
      PT <- penalty_goric(Amat        = Amat, 
                          meq         = meq, 
                          LP          = wt.bar, 
                          correction  = TRUE, 
                          sample.nobs = ldots$sample.nobs,
                          n_resp      = n_resp,
                          n_nuis      = n_nuis)
      if (type == "goricac") {
        PT <- PT - n_nuis
      }
    } else {
      stop("restriktor ERROR: unable to compute penalty term value.")  
    }
    
    ans$goric <- -2*z$loglik + ldots$penalty_factor*PT #-2*(z$loglik - PT)
    attr(ans$goric, "penalty") <- PT
    attr(ans$goric, "type") <- type
    attr(ans$goric, "loglik")  <- z$loglik 
  }
  
  ans$messages <- z$messages

  class(ans) <- c("summary.gorica_est")
  
  
  ans
}
