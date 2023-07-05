# is_closed_convex_cone <- function(A, b, x) {
#   convex <- closed <- FALSE
#   # is it convex
#   eigenvalues <- eigen(t(A) %*% A)$values
#   if (all(eigenvalues >= 0)) {
#     if(all((A %*% x) - b >= 0)) {
#       convex <- TRUE
#     }
#   } 
#   # is it closed
#   num_vars <- ncol(A)
#   rhs <- c(b, rep(0, num_vars))
#   constr_mat <- rbind(A, diag(num_vars))
#   res <- lp("min", objective.in = rep(1, num_vars+1), const.mat = constr_mat, 
#             const.dir = ">=", const.rhs = rhs, compute.sens = FALSE)
#   if (res$status == 0) {
#     closed <- TRUE
#   }
#   
#   out <- list(convex = convex, closed = closed)
#   return(out)
# }

