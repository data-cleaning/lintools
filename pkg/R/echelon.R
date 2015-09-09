#' Reduced row echelon form
#' 
#' Transform the equalities in a system of linear (in)equations or reduced row echelon form (RRE)
#' 
#' 
#' @section Details:
#' 
#' The parameters \code{A}, \code{b} and \code{neq} describe a system of the form \code{Ax<=b}, where
#' the first \code{neq} rows are equalities. The equalities are transformed to RRE form
#' where rows containing only zeros are deleted. 
#' 
#' 
#' @param A \code{[numeric]} matrix
#' @param b \code{[numeric]} vector
#' @param neq \code{[numeric]} The first \code{neq} rows of \code{A}, \code{b} are treated as equations.
#' @param tol \code{[numeric]} Values of magnitude less than \code{tol} are considered zero (for the purpose of handling
#' machine rounding errors).
#' 
#' @return 
#' A list with the following components:
#' \itemize{
#'   \item{\code{A}: the \code{A} matrix with equalities transformed to RRE form.}
#'   \item{\code{b}: the constant vector corresponding to \code{A}}
#'   \item{\code{neq}: the number of equalities in the resulting system.}
#' }
#' 
#' @export
echelon <- function(A, b, neq=nrow(A), tol=1e-8){
    check_sys(A=A,b=b,neq=neq,tol=tol)
  
    Ab <- cbind(A,b)
    
    ineq <- Ab[neq + seq_len(nrow(A)-neq),,drop=FALSE]
    Ab <- Ab[seq_len(neq),,drop=FALSE]
    
    k <- min(ncol(Ab),nrow(Ab))
    I <- seq_len(nrow(Ab))
    for ( i in 1:k ){
        I1 <- which(I >= i)
        ip <- I1[which.max(abs(Ab[I1,i]))]
        p <- Ab[ip,]
        if ( abs(p[i]) < tol ) next
        if ( ip > i ) Ab[c(ip,i),] <- Ab[c(i,ip),]
        Ab[-i,] <- Ab[-i,] - outer(Ab[-i,i],p/p[i])
    }
    
    d <- diag(Ab)
    id <- abs(d) > tol
    Ab[id,] <- Ab[id,]/d[id]
    I0 <- rowSums(abs(Ab) < tol) == ncol(Ab)
    Ab <- rbind(Ab[!I0,,drop=FALSE],Ab[I0,,drop=FALSE])
    i <- rowSums(abs(Ab) > tol) > 0 
    Ab <- rbind(Ab[i,,drop=FALSE],ineq)
    list(
      A = as.matrix(Ab[,1:ncol(A),drop=FALSE])
      , b = Ab[,ncol(A)+1,drop=TRUE]
      , neq = sum(i)
    )
}


