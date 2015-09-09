
#' Find independent blocks of equations.
#'
#' @param A \code{[numeric]} Matrix
#' @param tol \code{[numeric]} Coefficients with absolute value \code{< tol} are treated as zero.
#'
#' @return A \code{list} containing \code{numeric} vectors, each vector indexing an independent
#' block of rows in the system \code{Ax <= b}.
#' 
#'
block_index <- function(A, tol=1e-8){
  
  
  block <- function(B){
    x1 <- FALSE
    x <- B[1,]
    while (sum(x1 != x)){
      x1 <- x
      b <- sapply( 1:nrow(B)
                   , function(i){
                     any(B[i,] & x)
                   }
      )
      x <- colSums(B[b,,drop=FALSE]) > 0 #this is another way of "or"ring all found rows
    }
    b
  }
  
  D <- abs(cbind(A)) > tol
  orignames <- row.names(D)
  row.names(D) <- 1:nrow(D)
  
  #remove empty rows
  b <- rowSums(D) == 0
  D <- D[!b,,drop=FALSE]
  
  # create a list which will contain the blocks
  blocks <- list()
  L <- 1
  
  # detect and remove blocks until no blocks are left
  while (nrow(D) > 0){
    
    # find block
    b <- block(D)
    
    # store the original row numbers of the detected block
    blocks[[L]] <- as.integer(row.names(D)[b])
    L <- L + 1
    
    # remove the detected block
    D <- D[!b,,drop=FALSE]
  }
  lapply(blocks,function(b) {names(b)<-orignames[b]; b})
}


