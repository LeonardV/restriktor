


con_constraints <- function(model, VCOV, est, constraints, bvec = NULL, meq = 0L, 
                            debug = FALSE, auto_bound=FALSE, ...) {
  
##-------- build a bare-bones parameter table for this model--------------------
  # if model is a numeric vector
  if ("numeric" %in% class(model)) {
    parTable <- con_partable_est(model, est = TRUE, label = TRUE)
    parTable_org <- parTable
  } else {
    # if model is a fitted unrestricted object
    parTable <- con_partable(model, est = TRUE, label = TRUE)  
    parTable_org <- parTable
  }
  constraints <- unlist(constraints)

 
#-----when constraints are character type--------------------------
  
 
  if (is.character(constraints)) {
    
    #We split complex restrictions 
    sepc <- lengths(regmatches(constraints, gregexpr(":=", constraints)))
    semicol <- lengths(regmatches(constraints, gregexpr(";", constraints)))
    #operators check
    op1 <- lengths(regmatches(constraints, gregexpr("=~", constraints)))
    op2 <- lengths(regmatches(constraints, gregexpr("~", constraints)))
    op3 <- lengths(regmatches(constraints, gregexpr("~~", constraints)))
    if(sum(op1,op2,op3)>0){
      stop("You used incorrect operator. While defining the hypothesis make sure you use the correct specification model. ?restriktor for details on how to specify the constraint syntax or check the website, 
           https://restriktor.org/tutorial/syntax.html. ")
    } #note that this webside needs an update 
    
    if(all(grepl("[><]{2,}", constraints))==TRUE){
      stop("The message should not be seen if the code works fine! Unless you put <> or >< in between parameters which is not very smart.")
    }
    
    if(is.character(constraints) && sepc==0){
      hyp1<-list()
      for(i in 1:length(constraints)){
        semicolon_nr <- lengths(regmatches(constraints[i], gregexpr(";", constraints[i])))
        if(semicolon_nr >= 1){
          constraints[i]<-strsplit(constraints[i], split = ";")
          constraints[[i]]<-gsub("==","=",constraints[[i]])
          constraints[[i]]<-gsub("<=","<",constraints[[i]])
          constraints[[i]]<-gsub(">=",">",constraints[[i]])
          
          
        }else{
          constraints[i]<-constraints[i]
          constraints[[i]]<-gsub("==","=",constraints[[i]])
          constraints[[i]]<-gsub("<=","<",constraints[[i]])
          constraints[[i]]<-gsub(">=",">",constraints[[i]])
        }
        new<-unlist(constraints[i])
        if(all(grepl("[><=]{2,}", new))==FALSE){
          new_vect<-list()
          for(k in 1:length(new)){
            new_vect[[k]]<-expand_compound_constraints(new[k])
          }
          new_vect<-unlist(new_vect)
          hyp1[[i]]<-implode(new_vect,sep = ";")
          hyp1[[i]]<-gsub("=", "==",hyp1[[i]] )
        }else if(all(grepl("[><=]{2,}", new))==TRUE){
          stop("The message should not be seen if the code works fine! Unless you put <> or >< in between parameters which is not very smart.")
          
        }
      }
      constraints<-unlist(hyp1)
      
    }else if(is.character(constraints) && sepc>0 && semicol==0){
      hyp1<-list()
      for(i in 1:length(constraints)){
        constraints[i]<-strsplit(constraints[i], split = "\n")
        constraints[[i]]<-gsub("==","=",constraints[[i]])
        constraints[[i]]<-gsub("<=","<",constraints[[i]])
        constraints[[i]]<-gsub(">=",">",constraints[[i]])
        constraints[[i]]<-c(constraints[[i]],"")
        new<-unlist(constraints[i])
        if(all(grepl("[><=]{2,}", new))==FALSE){
          new_vect<-list()
          for(k in 1:length(new)){
            new_vect[[k]]<-expand_compound_constraints(new[k])
            sep_n <- lengths(regmatches(new_vect[[k]], gregexpr(":=", new_vect[[k]])))
            for (m in length(new_vect[[k]])) {
              if(sep_n[m]==0){
                new_vect[[k]][m]<-gsub("=", "==",new_vect[[k]][m])
                
              }
            }
          }
          new_vect<-unlist(new_vect)
          hyp1[[i]]<-implode(new_vect,sep = "\n")
        }else if(all(grepl("[><=]{2,}", new))==TRUE){
          stop("The message should not be seen if the code works fine! Unless you put <> or >< in between parameters which is not very smart.")
        }
        
      }
      constraints<-unlist(hyp1)
      
    }else if(is.character(constraints) && sepc>0 && semicol!=0){
      hyp1<-list()
      for(i in 1:length(constraints)){
        constraints[i]<-strsplit(constraints[i], split = ";")
        constraints[[i]]<-strsplit(constraints[[i]], split = "\n")
        constraints[[i]]<-unlist(constraints[[i]])
        constraints[[i]]<-gsub("==","=",constraints[[i]])
        constraints[[i]]<-gsub("<=","<",constraints[[i]])
        constraints[[i]]<-gsub(">=",">",constraints[[i]])
        constraints[[i]]<-c(constraints[[i]],"")
        new<-unlist(constraints[i])
        if(all(grepl("[><=]{2,}", new))==FALSE){
          new_vect<-list()
          for(k in 1:length(new)){
            new_vect[[k]]<-expand_compound_constraints(new[k])
            sep_n <- lengths(regmatches(new_vect[[k]], gregexpr(":=", new_vect[[k]])))
            for (m in length(new_vect[[k]])) {
              if(sep_n[m]==0){
                new_vect[[k]][m]<-gsub("=", "==",new_vect[[k]][m])
              }
            }
          }
          new_vect<-unlist(new_vect)
          hyp1[[i]]<-implode(new_vect,sep = "\n")
        }else if(all(grepl("[><=]{2,}", new))==TRUE){
          stop("The message should not be seen if the code works fine! Unless you put <> or >< in between parameters which is not very smart.")
        }
        
      }
      constraints<-unlist(hyp1)
      
      
    }
    
    # parse the constraints 
    CON <- lav_constraints_parse(constraints = constraints,
                                 partable    = parTable,
                                 debug       = debug,
                                 theta       = parTable$est)
     
    FLAT <- lavParseModelString(constraints)
    CON_FLAT <- attr(FLAT, "constraints")
    LIST <- list()
    lhs <- unlist(lapply(CON_FLAT, "[[", "lhs"))
    op  <- unlist(lapply(CON_FLAT, "[[", "op"))
    rhs <- unlist(lapply(CON_FLAT, "[[", "rhs"))
    LIST$lhs <- lhs
    LIST$op  <- op
    LIST$rhs <- c(LIST$rhs, rhs) 
    #i need to change these 
    parTable$lhs <- c(parTable$lhs, LIST$lhs)
    parTable$op <- c(parTable$op, LIST$op)     # each "==" need to change to ">" 
    parTable$rhs <- c(parTable$rhs, LIST$rhs)  # rhs represents 
    parTable$label <- c(parTable$label, rep("", length(lhs)))

    # equality constraints
    meq  <- nrow(con_constraints_ceq_amat(model, constraints = constraints))    
    # right-hand-side
    bvec <- con_constraints_rhs_bvec(model, constraints = constraints)
    # inequality constraints
    Amat <- con_constraints_con_amat(model, constraints = constraints)
    
    if(meq!=0 && auto_bound==TRUE){
      #change parameters 
 #     sign<-parTable$op 
 #     a<-match("==",sign)
 #     while (!is.na(a)) {
 #       l<-a[1]
 #       part1<-sign[1:(l-1)]
 #       part2<-sign[-(1:l)]
 #       midd<-c(">",">")
 #       sign<-c(part1,midd,part2)
  #      a<-match("==",sign)
  #    } 
  #    parTable$op<- sign
      
      Amat_eq<-matrix(Amat[1:meq,], nrow=meq)
      bound<-0.75
      prox<-SE<-c()
      Amat_prox_2<-matrix(NA,ncol = ncol(Amat))
      rep<-matrix()
      for(i in 1:nrow(Amat_eq)){
        #Calculate the Standard errors used to generate the equality bounds 
        SE[i]<-sum(Amat_eq[i,]%*%t(Amat_eq[i,])*VCOV)
        rep<-Amat_eq[rep(i, times = 2),]
        Amat_prox<-rep
        Amat_prox_2<-rbind(Amat_prox_2,Amat_prox)
        prox[i]<-(-1)*SE[i]*bound
      }
      Amat_prox<-Amat_prox_2[-1,]
      
      ans <- seq(2,nrow(Amat_prox),2)
      for(k in ans){
        Amat_prox[k,]<-(-1)*Amat_prox[k,]
      }
      #vector of the boundry values 
      new<-c()
      new1<-c()
      bvec_eq<-bvec[1:meq]
      #repeat the values dor the rhs twice 
      for(n in 1:meq){ 
        prox1<-prox[rep(n, time=2)]
        bvec_eq1<-bvec_eq[rep(n, time=2)]
        new1<-c(new1,bvec_eq1)
        new<-c(new,prox1)
      }
      rhs_eq<-new1+new
      bvec1<-c(rhs_eq,bvec[-(1:meq)])
      Amat<-rbind(Amat_prox,Amat[-(1:meq),])
      meq<-0
      bvec<-bvec1
    }
    
    # check for not supported constraint syntax or an empty matrix
    nsc_lhs.idx <- sum(grepl("<|>|=", parTable$lhs))
    nsc_rhs.idx <- sum(grepl("<|>|=", parTable$rhs))
 ### Now this message should not be the case anymore but I leaved it in case it might contradict something else and as a message that something went wrong while reparating compound hypothesis   
    if (all(Amat == 0) | nsc_lhs.idx > 0 | nsc_rhs.idx > 0) {
      stop("Restriktor ERROR: Sorry, but I have no idea how to deal with your constraint syntax. \n",
           "See ?restriktor for details on how to specify the constraint syntax or check the website \n", 
            "https://restriktor.org/tutorial/syntax.html. \n\n",  
           "Hint: constraints have to be specified pairwise, e.g., x1 < x2; x2 == x3; x1 > 2", sep = "",
          call. = FALSE
      )
    }
    
    # In case of abs() the contraints may incorretly be considered as non-linear. 
    # Here, we remove the abs() from the constraint function which is redundant 
    # for determining if the constraints are linear. 
    
    # check if any abs() functie exists in string. 
    if (any(grepl("abs\\(.*\\)", c(LIST$lhs, LIST$rhs)))) {
      LIST2 <- LIST
      
      # reomve abs( and ) from string
      LIST2$lhs <- gsub("abs\\(|\\)", "", LIST2$lhs)
      LIST2$rhs <- gsub("abs\\(|\\)", "", LIST2$rhs)
      
      parTable_org$free <- seq_len(length(parTable_org$lhs))
      cin.function <- lav_partable_constraints_ciq(partable = parTable_org, con = LIST2)
      ceq.function <- lav_partable_constraints_ceq(partable = parTable_org, con = LIST2)
      
      CON$cin.nonlinear.idx <- con_constraints_nonlinear_idx(func = cin.function, 
                                                             npar = length(parTable_org$est))
      CON$ceq.nonlinear.idx <- con_constraints_nonlinear_idx(func = ceq.function, 
                                                             npar = length(parTable_org$est))
    }
    
    
    CON$constraints <- constraints
    
#-------if constraints are not of type character and are not NULL--------------
  } else if (!is.character(constraints) && !is.null(constraints)) {
    if (is.vector(constraints) ) {
      constraints <- rbind(constraints)
    }
    CON <- NULL
    Amat <- constraints
    bvec <- if (is.null(bvec)) { rep(0L, nrow(Amat)) } else { bvec }
    meq  <- if (is.null(meq)) { 0L } else { meq }
    
    
    #Calculate the Standard errors used to generate the equality bounds 
    
    if(meq!=0 && auto_bound==TRUE){
      Amat_eq<-matrix(Amat[1:meq,], nrow=meq)
      bound<-0.75
      prox<-SE<-c()
      Amat_prox_2<-matrix(NA,ncol = ncol(Amat))
      rep<-matrix()
      for(i in 1:nrow(Amat_eq)){
        SE[i]<-sum(Amat_eq[i,]%*%t(Amat_eq[i,])*VCOV)
        rep<-Amat_eq[rep(i, times = 2),]
        Amat_prox<-rep
        Amat_prox_2<-rbind(Amat_prox_2,Amat_prox)
        prox[i]<-(-1)*SE[i]*bound
      }
      Amat_prox<-Amat_prox_2[-1,]
      
      ans <- seq(2,nrow(Amat_prox),2)
      for(k in ans){
        Amat_prox[k,]<-(-1)*Amat_prox[k,]
      }
      #vector of the boundry values 
      new<-c()
      new1<-c()
      bvec_eq<-bvec[1:meq]
      #repeat the values dor the rhs twice 
      for(n in 1:meq){ 
        prox1<-prox[rep(n, time=2)]
        bvec_eq1<-bvec_eq[rep(n, time=2)]
        new1<-c(new1,bvec_eq1)
        new<-c(new,prox1)
      }
      rhs_eq<-new1+new
      bvec1<-c(rhs_eq,bvec[-(1:meq)])
      Amat<-rbind(Amat_prox,Amat[-(1:meq),])
      meq<-0
      bvec<-bvec1
    }
    
 
    
  } else { 
    stop("no restriktions were specified.") 
  }
  
  if (!(nrow(Amat) == length(bvec))) {
    warning("restriktor WARNING: The number of constraints does not match 
                    the \'rhs\' (nrow(Amat) != length(rhs)).")
  }
  
  if (meq > nrow(Amat)) { 
    stop("restriktor ERROR: The maximum number of equality constraints = ", nrow(Amat), "\n")
  }
  
  if (length(CON$ceq.nonlinear.idx) > 0L || length(CON$cin.nonlinear.idx) > 0L) {
    stop("restriktor ERROR: can not handle (yet) nonlinear (in)equality restriktions")
  }
  
  if (debug && is.character(constraints)) {
    print(as.data.frame(parTable, stringsAsFactors = FALSE))
    print(CON)
  }
  
  
  # rAmat <- GaussianElimination(t(Amat))

  ## still to catch 
  #H1 <- 'x1 < 4; x1 > 1' # range restrictie
  #H1 <- 'x1 < 1; x1 > 1' # equality
  #H1 <- 'x1 > 3; x1 > 4' # 
  #H1 <- 'x1 > -1; x1 > 4'#
  
  # if (mix.weights == "pmvnorm") {
  #   if (rAmat$rank < nrow(Amat) && rAmat$rank != 0L) {
  #     ## check for inconsistent constraints: quadprog gives an error if constraints
  #     ## are inconsistent
  #     # consistent.check <- con_solver_gorica(est  = est, 
  #     #                                       VCOV = VCOV, 
  #     #                                       Amat = Amat, 
  #     #                                       bvec = bvec, 
  #     #                                       meq  = meq)
  #     
  #     ## remove any linear dependent rows from the constraint matrix. Amat
  #     ## must be of full row rank.
  #     # remove any zero vectors
  #     allZero.idx <- rowSums(abs(Amat)) == 0
  #     Amat <- Amat[!allZero.idx, , drop = FALSE]
  #     bvec <- bvec[!allZero.idx]
  #     # rank Amat
  #     rank <- qr(Amat)$rank 
  #     # singular value decomposition
  #     s <- svd(Amat)
  #     # continue untill Amat is of full-row rank
  #     while (rank != length(s$d)) {
  #       # check which singular values are zero
  #       zero.idx <- which(zapsmall(s$d) <= 1e-16)
  #       # remove linear dependent rows and reconstruct the constraint matrix
  #       Amat <- s$u[-zero.idx, ] %*% diag(s$d) %*% t(s$v)
  #       # zapping small ones to zero
  #       Amat <- zapsmall(Amat)
  #       bvec <- bvec[-zero.idx]
  #       s <- svd(Amat)
  #       if (debug) {
  #         cat("rank = ", rank, " ... non-zero det. = ", length(s$d), "\n")
  #       }
  #     }
  #   }
  # } else if (rAmat$rank < nrow(Amat) &&
  #            !(se %in% c("none", "boot.model.based", "boot.standard")) &&
  #            rAmat$rank != 0L) {
  #   warning(paste("Restriktor Warning: No standard errors could be computed.
  #                     The constraint matrix must be full row-rank.
  #                     Try to set se = \"none\", \"boot.model.based\" or \"boot.standard\".")) 
  # }
  
  OUT <- list(CON      = CON, 
              parTable = parTable,
              Amat     = Amat,
              bvec     = bvec, 
              meq      = meq)
  
  OUT
}

# It creates a matrix of equality constraints
con_constraints_ceq_amat <- function(object, constraints = NULL) {

  # build a bare-bones parameter table for this object
  if ("numeric" %in% class(object)) {
    lavpartable <- con_partable_est(object, est = TRUE, label = TRUE)
  } else {
  # if object is a fitted unrestricted object
    lavpartable <- con_partable(object, est = TRUE, label = TRUE)  
  }

  # parse the constraints
  CON <- lav_constraints_parse(constraints = constraints,
                               partable    = lavpartable,
                               theta       = lavpartable$est)

  CON$ceq.JAC
}

# It creates a matrix of inequality terms 
con_constraints_con_amat <- function(object, constraints = NULL) {
  
  # build a bare-bones parameter table for this object
  if ("numeric" %in% class(object)) {
    lavpartable <- con_partable_est(object, est = TRUE, label = TRUE)
  } else {
    # if object is a fitted unrestricted object
    lavpartable <- con_partable(object, est = TRUE, label = TRUE)  
  }
  # parse the constraints
  CON <- lav_constraints_parse(constraints = constraints,
                               partable    = lavpartable, 
                               theta       = lavpartable$est)

  rbind(CON$ceq.JAC, CON$cin.JAC)
}



con_constraints_rhs_bvec <- function(object, constraints = NULL) {

  # build a bare-bones parameter table for this object
  #lavpartable <- con_partable(object, est = TRUE, label = TRUE)

  # build a bare-bones parameter table for this object
  if ("numeric" %in% class(object)) {
    lavpartable <- con_partable_est(object, est = TRUE, label = TRUE)
  } else {
    # if object is a fitted unrestricted object
    lavpartable <- con_partable(object, est = TRUE, label = TRUE)  
  }
  
  
  # parse the constraints
  CON <- lav_constraints_parse(constraints = constraints,
                               partable    = lavpartable,
                               theta       = lavpartable$est)
  
  #returns numeric(0) if any con_constraints_ceq_amat or con_constraints_con_amat is empty and the rhs value for non empty matrix starting from the exuality 

  c(CON$ceq.rhs, CON$cin.rhs)
}






