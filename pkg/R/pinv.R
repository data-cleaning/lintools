
#' Moore-Penrose pseudoinverse
#'
#' Compute the pseudoinverse of a matrix using the
#' SVD-construction 
#'
#' @param A [numeric] matrix
#' @param eps [numeric] tolerance for determining zero singular values
#'
#' @export
pinv <- function(A, eps=1e-8){
  L <- svd(A)
  d <- L$d
  i <- abs(d) > eps
  d[i] <- 1/d[i]
  L$v %*% diag(d) %*% t(L$u)
}















  
