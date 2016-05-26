# Code written by John Fox
# Last modified 9 July 2007
# Found on https://stat.ethz.ch/pipermail/r-help/2007-September/140021.html
# Slightly adapted by LV
GaussianElimination <- function(A, B, tol = sqrt(.Machine$double.eps), 
                                verbose = FALSE, fractions = FALSE){
  # A: coefficient matrix
  # B: right-hand side vector or matrix
  # tol: tolerance for checking for 0 pivot
  # verbose: if TRUE, print intermediate steps
  # fractions: try to express nonintegers as rational numbers
  # If B is absent returns the reduced row-echelon form of A.
  # If B is present, reduces A to RREF carrying B along.
  if (fractions) {
    mass <- require(MASS)
    if (!mass) {
      stop("fractions=TRUE needs MASS package")
    }
  }
  if ((!is.matrix(A)) || (!is.numeric(A))) {
    stop("argument must be a numeric matrix")
  }
  n <- nrow(A)
  m <- ncol(A)
  if (!missing(B)) {
    B <- as.matrix(B)
    if (!(nrow(B) == nrow(A)) || !is.numeric(B)) {
      stop("argument must be numeric and must match the number of row of A")
    }
    A <- cbind(A, B)
  }
  i <- j <- 1
  pivot.new <- c()
  while (i <= n && j <= m) {
    while (j <= m) {
      currentColumn <- A[,j]
      currentColumn[1:n < i] <- 0
      # find maximum pivot in current column at or below current row
      which <- which.max(abs(currentColumn))
      pivot <- currentColumn[which]
      if (abs(pivot) <= tol) { # check for 0 pivot
        j <- j + 1
        next
      }     
      if (which > i) {
        A[c(i, which),] <- A[c(which, i),]  # exchange rows
      }
      A[i,] <- A[i,] / pivot            # pivot
      row <- A[i,]
      A <- A - outer(A[,j], row)      # sweep
      A[i,] <- row                    # restore current row
      pivot.new <- c(pivot.new, j)
      if (verbose) {
        if (fractions) {
          print(fractions(A))
        } else {
          print(round(A, round(abs(log(tol,10)))))
        }
      }
      j <- j + 1
      break
    }
    i <- i + 1
  }
  # 0 rows to bottom
  zeros <- which(apply(A[,1:m], 1, function(x) max(abs(x)) <= tol))
  if (length(zeros) > 0){
    zeroRows <- A[zeros,]
    A <- A[-zeros,]
    A <- rbind(A, zeroRows)
  }
  rownames(A) <- NULL
  
  OUT <- list(A     = round(A, round(abs(log(tol, 10)))), 
              pivot = pivot.new, 
              rank  = length(pivot.new))
  
  OUT
}

