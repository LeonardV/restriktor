## input options
# 1. est + vcov
# 2. LL + PT
# 3. IC values (AIC, ORIC, GORIC, GORICA)
# 4. IC weights or (Bayesian) posterior model probs. 

# In case of an equal-evidence approach, aggregating evidence from, say, 5 studies 
# with n=100 observations is the same as obtaining evidence from 1 study 
# (as if it was possible) with n=500 observations (like meta-analysis does).
# In the added-evidence approach, the aggregated evidence from, says, 5 studies 
# is stronger than as if the data were combined (as if that was possible).


#' @param TypeEv The type of evidence-synthesis approach: Equal-evidence approach (0) or Added-evidence approach (1).
#' @param S The number of (primary) studies. That is, the results (evidence) of S studies will be aggregated.
#' @param Param_studies List of S 'named' vectors with the k_s (standardized) parameter estimates of interest of Study s. Thus, there are S items in the list and each item is a 'named' vector with k_s elements: the k_s number of parameter estimates ratioevant for that study. In case each study has the same number of parameters (k) which denote the same (in terms of hypothesis specification), Param_studies can be an S x k 'named' matrix. Note: The names of the vectors (or the column names of the S x 'k' matrix) with estimates should be used in the hypothesis specification.
#' @param CovMx_studies List of the S covariance matrices of the (standardized) parameter estimates of interest (of size k_s x k_s). In case number of parameters are the same, it can also be a S*k_s x k_s matrix. Note: The columns (and rows) do not need to be named.
#' @param Hypo_studies A vector of strings containing the NrHypos theory-based hypotheses. If SameHypo = 0, then there should be S specifications of the NrHypos theory-based hypotheses, that is S times NrHypos strings.
#' @param comparison Indicator of which safeguard-hypothesis should be used: "unconstrained" (default; i.e., all possible theories including the one specified), "none" (only advised when set of hyptheses cover all theories), or (only when 'NrHypos = 1') "complement" (i.e., the remaining theories).
#'


# TODO
# print functies



evSyn <- function(object, ...) { UseMethod("evSyn") }


evSyn <- function(object, ...) {
  
  arguments <- list(...)
  if (length(arguments)) {
    pnames <- c("VCOV", "PT", "constraints", "type", "comparison", "hypo_names")
    pm <- pmatch(names(arguments), pnames, nomatch = 0L)
    if (any(pm == 0L)) { 
      pm.idx <- which(pm == 0L)
      stop("Restriktor Error: ", names(arguments[pm.idx]), " invalid argument(s).")
    }
  }
  
  ## est + VCOV
  ## check if list contains a vector/matrix with the the parameter estimates
  if (is.list(object) & is.list(arguments$VCOV)) {
    obj_class <- unique(unlist(lapply(object, class)))
    
    # if matrix, then make it numeric vectors
    if ("matrix" %in% obj_class) {
      vnames <- lapply(object, function(x) colnames(x))
      object <- lapply(object, function(x) as.numeric(x))
      # colnames are not inherited, so these are fixed here.
      for (i in 1:length(object)) {
        names(object[[i]]) <- vnames[[i]]
      }
    } 
    evSyn.est(object, ...)
  } else if (inherits(object, c("matrix", "data.frame")) & 
             !is.null(arguments$PT)) {
    # LL + PT
    evSyn.LL(object, ...)
  } else if (inherits(object, c("matrix", "data.frame")) & 
             all(rowSums(object) - 1 > .Machine$double.eps) &
             is.null(arguments$PT) & 
             is.null(arguments$VCOV)) { 
    # IC values
    evSyn.ICvalues(object, ...)
  } else if (inherits(object, "matrix") & is.null(arguments$PT) & 
             all(rowSums(object) - 1 <= .Machine$double.eps) &
             is.null(arguments$VCOV)) { 
    # IC penalty weights
    evSyn.ICweights(object, ...)
  } else {
    stop("restriktor Error: I don't know how to handle the input.")
  }
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on the (standard) parameter estimates and the covariance matrix
evSyn.est <- function(object, VCOV = NULL, constraints = NULL,
                      type = c("equal", "added"), 
                      comparison = c("unconstrained", "complement", "none"), ...) {
  
  if (!is.list(object) && is.null(VCOV)) {
    stop("Restriktor Error: object must be a list with (standardized) parameter estimates for each study.")
  }
  if (length(object) != length(VCOV)) {
    stop("Restriktor Error: object must have the same length as VCOV.")
  }
  
  # number of primary studies
  S <- length(object)
  V <- length(VCOV)
  
  if (!is.list(VCOV)) {
    stop("Restriktor Error: VCOV must be a list with the variance-covariance matrices for each study.")
  } else {
    # check if the matrices are all symmetrical
    VCOV_isSym <- sapply(VCOV, isSymmetric)
    if (!all(VCOV_isSym)) {
      sprintf("Restriktor Error: the %sth covariance matrix in VCOV is not symmetric.", which(!VCOV_isSym))
      stop(call. = FALSE)  
    }
  }
  
  if (S != V) {
    stop("Restriktor Error: the number of studies does not match the number of VCOV.")
  }
  if (!is.list(object)) {
    stop("Restriktor Error: object must be list with (standardized) parameter estimates for each study.", call. = FALSE)
  }
  # if object is a list of class matrix, make it a vector
  if (any(sapply(object, is.matrix))) {
    object <- lapply(object, function(x) unlist(as.data.frame(x)))
  }
  
  if (!all(sapply(object, is.vector))) {
    stop("Restriktor Error: object must be a list with named numeric vectors", call. = FALSE)
  }
  
  # if only 1 hypothesis is provided, that hypothesis will be applied to all studies
  # This also means that the vector with parameter estimates must have equal names between the studies.
  # check if the user specified the same hypothesis in each study. 
  # number of hypotheses (this should be equal in each study, see check below)
  NrHypos <- unique(sapply(constraints, length))
  len_H <- length(unique(sapply(constraints, function(x) gsub("\\s", "", x))))
  if (len_H == 1L | NrHypos == 1) {
    SameHypo <- TRUE  
    Nrhypos <- 1L
  } else {
    SameHypo <- FALSE  
  }
    
  # number of constraints must be equal for each studie. In each study a set of shared theories (i.e., constraints) are compared.
  if (length(NrHypos) > 1L) {
    stop("Restriktor Error: the number of constraints must be equal in each study.", call. = FALSE)
  }
  if (comparison == "complement" & NrHypos > 1L) {
    warning("Restriktor Warning: if comparison = 'complement', only one order-restricted hypothesis\n",
            " is allowed (for now). Therefore, comparison is set to 'unconstrained'.",
            call. = FALSE)
  }
  
  # if NrHypos == 1, this means that we apply the same hypotheses to each study
  if (NrHypos == 1) {
    NrHypos <- length(constraints)
    NrHypos_incl <- NrHypos + 1
    sameHypo <- TRUE
    if (comparison == "none"){
      NrHypos_incl <- NrHypos
    }
  } else {
    NrHypos_incl <- NrHypos + 1
    if (comparison == "none"){
      NrHypos_incl <- NrHypos
    }
    sameHypo <- FALSE
  }
  
  LL <- PT <- weight_m <- GORICA_m <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
    rownames(LL) <- rownames(PT) <- paste0("study_", 1:S)

  CumulativeGoricaWeights <- CumulativeGorica <- matrix(NA, nrow = S+1, ncol = NrHypos_incl)
    rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(paste0("study ", 1:S), "final")
    
  if (NrHypos == 1 & comparison == "complement") {
    hnames <- c("H1", "Hc")
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = 1)
  } else if (comparison == "none"){
    hnames <- c(paste0("H", 1:NrHypos))
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  } else {
    hnames <- c(paste0("H", 1:NrHypos), "Hu")
    ratio.weight_mu <- matrix(data = NA, nrow = S, ncol = NrHypos_incl)
  }
    
  colnames(GORICA_m) <- colnames(weight_m) <- colnames(LL) <- colnames(PT) <- colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- hnames
  rownames(GORICA_m) <- rownames(ratio.weight_mu) <- rownames(weight_m) <- paste0("study_", 1:S)
  

  for(s in 1:S) {
    if (sameHypo) {
      res_goric <- goric(object[[s]], VCOV = VCOV[[s]], 
                         constraints = constraints,
                         type = 'gorica', comparison = comparison)
    } else {
      res_goric <- goric(object[[s]], VCOV = VCOV[[s]], 
                         constraints = constraints[[s]],
                         type = 'gorica', comparison = comparison)
    }
    
    if (comparison == "unconstrained") {
      ratio.weight_mu[s,] <- res_goric$ratio.gw[, NrHypos_incl]
    } else if (comparison == "complement") {
      ratio.weight_mu[s,] <- res_goric$ratio.gw[1, NrHypos_incl]
    } 
    
  
    LL[s, ] <- res_goric$result$loglik
    PT[s, ] <- res_goric$result$penalty
    GORICA_m[s, ] <- res_goric$result$gorica
    weight_m[s, ] <- res_goric$result$gorica.weights
  }
  
  sumPT <- sumLL <- 0
  if (type == "added") { 
    # added-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
    EvSyn_approach <- "Added-evidence approach"
  } else { 
    # equal-evidence approach
    for(s in 1:S) {
      sumLL <- sumLL + LL[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s,] <- -2 * sumLL + 2 * sumPT/s
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s,] <- exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric)))
    }
    EvSyn_approach <- "Equal-evidence approach"
  }
  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  
  rownames(Final.ratio.GORICA.weights) <- hnames
  
  hypotheses <- res_goric$constraints_usr
  names(hypotheses) <- hnames
  
  # Output
  if (NrHypos == 1 & comparison == "complement") {
    colnames(ratio.weight_mu) <- c("H1 vs. Hc1")
    colnames(Final.ratio.GORICA.weights) <- c("vs. H1", "vs. Hc")
    
    out <- list(GORICA_m = GORICA_m, GORICA.weight_m = weight_m, 
                ratio.GORICA.weight_mc = ratio.weight_mu, LL_m = LL, PT_m = PT,
                EvSyn_approach = EvSyn_approach, Cumulative.GORICA = CumulativeGorica, 
                Cumulative.GORICA.weights = CumulativeGoricaWeights,
                Final.ratio.GORICA.weights = Final.ratio.GORICA.weights,
                hypotheses = constraints_usr)
  } else if (comparison == "none") {
    colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. H", 1:NrHypos))
    
    out <- list(GORICA_m = GORICA_m, GORICA.weight_m = weight_m, LL_m = LL, PT_m = PT,
                EvSyn_approach = EvSyn_approach, Cumulative.GORICA = CumulativeGorica, 
                Cumulative.GORICA.weights = CumulativeGoricaWeights,
                Final.ratio.GORICA.weights = Final.ratio.GORICA.weights,
                hypotheses = hypotheses)
  } else { 
    # unconstrained
    colnames(ratio.weight_mu) <- c(paste0("H", 1:NrHypos, " vs. Unc."), "Unc. vs. Unc.")
    colnames(Final.ratio.GORICA.weights) <- c(paste0("vs. H", 1:NrHypos), "vs. Hu")
    
    out <- list(GORICA_m        = GORICA_m, 
                GORICA.weight_m = weight_m, 
                ratio.GORICA.weight_mu = ratio.weight_mu, 
                LL_m = LL, 
                PT_m = PT,
                EvSyn_approach    = EvSyn_approach, 
                Cumulative.GORICA = CumulativeGorica, 
                Cumulative.GORICA.weights  = CumulativeGoricaWeights,
                Final.ratio.GORICA.weights = Final.ratio.GORICA.weights)
  }
  
  class(out) <- c("evSyn.est", "evSyn")
  
  return(out)
}




# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on log likelihood and penalty values
evSyn.LL <- function(object, PT = NULL, type = c("added", "equal"),
                     hypo_names = NULL) {
  
  if (is.null(PT)) {
    stop("restriktor Error: no penalty term values found.")
  }
  
  if (inherits(object, "data.frame")) {
    object <- as.matrix(object)  
  }
  if (inherits(PT, "data.frame")) {
    PT <- as.matrix(PT)  
  }
  
  if (!inherits(object, c("matrix", "data.frame"))) {
    stop("Restriktor Error: object must be a data.frame or matrix with the log-likelihood values.")
  }
  
  if (!inherits(PT, c("matrix", "data.frame"))) {
    stop("Restriktor Error: PT must be a data.frame or matrix with the penalty term values.")
  }
  
  S <- nrow(object)
  NrHypos <- ncol(object) - 1
  if (is.null(hypo_names)) {
    hnames <- paste0("H", 1:(NrHypos + 1))
  } else {
    hnames <- hypo_names
  }
  
  weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(LL) <- colnames(PT) <- colnames(weight_m) <- hnames
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(paste0("Study_", 1:S), "Final")
  rownames(LL) <- rownames(PT) <- rownames(weight_m) <- paste0("Study_", 1:S)
  
  
  sumLL <- 0
  sumPT <- 0
  IC <- -2 * LL + 2 * PT
  if (type == "added") { 
    # added-ev approach
    for(s in 1:S) {
      minIC <- min(IC[s, ])
      weight_m[s, ] <- as.matrix(exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC))))
      
      sumLL <- sumLL + LL[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- as.matrix(-2 * sumLL + 2 * sumPT)
      
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- as.matrix(exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric))))
    }
    EvSyn_approach <- "Added-evidence approach"
  } else { 
    # equal-ev approach
    for(s in 1:S){
      sumLL <- sumLL + LL[s, ]
      sumPT <- sumPT + PT[s, ]
      CumulativeGorica[s, ] <- as.matrix(-2 * sumLL + 2 * sumPT/s)
      
      minGoric <- min(CumulativeGorica[s, ])
      CumulativeGoricaWeights[s, ] <- as.matrix(exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric))))
    }
    EvSyn_approach <- "Equal-evidence approach"
  }
  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S,]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S,]
  Final.rel.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  rownames(Final.rel.GORICA.weights) <- hnames
  colnames(Final.rel.GORICA.weights) <- paste0("vs. ", hnames)
  
  out <- list(LL_m = LL, 
              PT_m = PT, 
              GORICA_m = IC, 
              GORICA.weight_m   = weight_m,
              EvSyn_approach    = EvSyn_approach, 
              Cumulative.GORICA = CumulativeGorica, 
              Cumulative.GORICA.weights = CumulativeGoricaWeights,
              Final.rel.GORICA.weights  = Final.rel.GORICA.weights)
  
  
  class(out) <- c("evSyn.LL", "evSyn")
  
  return(out)
  
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA values
evSyn.ICvalues <- function(object, hypo_names = NULL) {
  
  if (inherits(object, "data.frame")) {
    object <- as.matrix(object)  
  }
  
  IC <- object
  S  <- nrow(IC)
  
  NrHypos <- ncol(IC) - 1
  weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  
  if (is.null(hypo_names)) {
    hnames <- paste0("H", 1:(NrHypos+1))
  } else {
    hnames <- hypo_names
  }
  
  weight_m <- matrix(NA, nrow = S, ncol = (NrHypos + 1))
  CumulativeGorica <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  CumulativeGoricaWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  colnames(CumulativeGorica) <- colnames(CumulativeGoricaWeights) <- colnames(IC) <- colnames(weight_m) <- hnames
  rownames(CumulativeGorica) <- rownames(CumulativeGoricaWeights) <- c(paste0("Study", 1:S), "Final")
  rownames(IC) <- rownames(weight_m) <- paste0("Study", 1:S)
  
  sumIC <- 0
  for(s in 1:S){
    minIC <- min(IC[s, ])
    weight_m[s, ] <- matrix(exp(-0.5*(IC[s, ]-minIC)) / sum(exp(-0.5*(IC[s, ]-minIC))))
  
    sumIC <- sumIC + IC[s, ]
    CumulativeGorica[s, ] <- sumIC
    minGoric <- min(CumulativeGorica[s, ])
    CumulativeGoricaWeights[s, ] <- matrix(exp(-0.5*(CumulativeGorica[s, ]-minGoric)) / sum(exp(-0.5*(CumulativeGorica[s, ]-minGoric))))
  }
  
  EvSyn_approach <- "Added-evidence approach (which is the only option when the input consists of IC's)"
  
  CumulativeGorica[(S+1), ] <- CumulativeGorica[S, ]
  CumulativeGoricaWeights[(S+1), ] <- CumulativeGoricaWeights[S, ]
  
  Final.GORICA.weights <- CumulativeGoricaWeights[S, ]
  Final.ratio.GORICA.weights <- Final.GORICA.weights %*% t(1/Final.GORICA.weights)
  
  rownames(Final.ratio.GORICA.weights) <- hnames
  colnames(Final.ratio.GORICA.weights) <- paste0("vs. ", hnames)
  
  out <- list(GORICA_m          = IC, 
              GORICA.weight_m   = weight_m,
              EvSyn_approach    = EvSyn_approach, 
              Cumulative.GORICA = CumulativeGorica, 
              Cumulative.GORICA.weights  = CumulativeGoricaWeights,
              Final.ratio.GORICA.weights = Final.ratio.GORICA.weights)
  
  class(out) <- c("evSyn.ICvalues", "evSyn")
  
  return(out)
}



# -------------------------------------------------------------------------
# GORIC(A) evidence synthesis based on AIC or ORIC or GORIC or GORICA weights or (Bayesian) posterior model probabilities
evSyn.ICweights <- function(object, PriorWeights = NULL, hypo_names = NULL) {
  
  if (inherits(object, "data.frame")) {
    object <- as.matrix(object)  
  }
  
  Weights <- object
  NrHypos <- ncol(Weights) - 1
  S <- nrow(object)
  
  if (is.null(PriorWeights)) {
    PriorWeights <- rep(1/(NrHypos + 1), (NrHypos + 1))
  }
  # To make it sum to 1 (if it not already did)
  PriorWeights <- PriorWeights / sum(PriorWeights) 
  
  if (is.null(hypo_names)) {
    hypo_names <- paste0("H", 1:(NrHypos+1))
  }
  
  CumulativeWeights <- matrix(NA, nrow = (S+1), ncol = (NrHypos + 1))
  colnames(CumulativeWeights) <- hypo_names
  rownames(CumulativeWeights) <- c(paste0("Study", 1:S), "Final")
  
  CumulativeWeights[1, ] <- PriorWeights * Weights[1, ] / sum( PriorWeights * Weights[1, ] )
  
  for(s in 2:S) {
    CumulativeWeights[s, ] <- CumulativeWeights[(s-1), ] * Weights[s, ] / sum( CumulativeWeights[(s-1), ] * Weights[s, ] )
  }
  EvSyn_approach <- "Added-evidence approach (which is the only option when the input consists of Weights)"
  
  CumulativeWeights[(S+1), ] <- CumulativeWeights[S, ]
  
  Final.weights <- CumulativeWeights[S, ]
  Final.rel.weights <- Final.weights %*% t(1/Final.weights)
  rownames(Final.rel.weights) <- hypo_names
  colnames(Final.rel.weights) <- paste0("vs. ", hypo_names)
  
  out <- list(weight_m           = Weights,
              EvSyn_approach     = EvSyn_approach,
              Cumulative.weights = CumulativeWeights,
              Final.rel.weights  = Final.rel.weights)
  
  class(out) <- c("evSyn.ICweights", "evSyn")
  
  return(out)
  
}



## plot function
plot.evSyn <- function(x, ...) {
  namesH <- colnames(x$GORICA_m)
  NrHypos_incl <- ncol(x$GORICA_m)
  S <- nrow(x$GORICA_m)
  Name_studies <- as.factor(1:S)
  
  weight_m <- x$GORICA.weight_m
  CumulativeGoricaWeights <- x$Cumulative.GORICA.weights
  
  # Create data frame for per study weights
  per_study_df <- data.frame(study = rep(Name_studies, NrHypos_incl),
                             weight = c(weight_m))
  per_study_df$weight_type <- "per study"
  
  # Create data frame for cumulative weights
  cumulative_df <- data.frame(study = rep(Name_studies, NrHypos_incl),
                              weight = c(CumulativeGoricaWeights[1:S, ]))
  cumulative_df$weight_type <- "cumulative"
  
  # Combine data frames
  plot_data <- rbind(per_study_df, cumulative_df)
  plot_data$variable <- rep(rep(namesH, each = S), times = 2)
  
  # Create plot
  ggplot(plot_data, aes(x = .data[['study']], 
                        y = .data[['weight']], 
                        shape    = .data[['weight_type']], 
                        linetype = .data[['weight_type']], 
                        color    = .data[['variable']])) +
    geom_point(size = 3) +
    geom_line(data = plot_data[plot_data[['weight_type']] == "cumulative", ], 
              aes(group = .data[['variable']]), linewidth = 1) +
    #scale_shape_manual(values = c(8, 16)) +
    #scale_linetype_manual(values = c("solid", "dashed")) +
    #scale_color_manual(values = c("#D81B60", "#1E88E5")) +
    theme(
      plot.margin = unit(c(1,1,1,1), "cm"),
      legend.position = "bottom",
      legend.margin = margin(t = 10, r = 0, b = 3, l = 0),
      #legend.title = element_text(size = 18),
      legend.key = element_blank(),
      legend.text = element_text(size = 12),
      axis.text.x  = element_text(size = 12), axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, vjust = -3), axis.title.y = element_text(size = 14, vjust = 5),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    scale_y_continuous(limits = c(0, 1), n.breaks = 10) +
    scale_x_discrete(expand = c(0, 0.05)) +
    labs(x = "Studies", y = "GORIC(A) weights",
         title = "GORIC(A) weights per study and cumulative",
         shape = "", color = "", linetype = "") +
    guides(shape = guide_legend(override.aes = list(linetype = 0)))
}


