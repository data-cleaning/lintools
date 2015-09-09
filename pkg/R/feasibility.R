
# find straigtforward contradictions of the form 0 <= b or 0 == b 
has_contradiction <- function(A,b,tol){
  if (nrow(A)==0) return(FALSE)
  if (ncol(A)==0 & abs(b) > tol) return(TRUE)
  any(rowSums(abs(A)>tol) == 0 & (abs(b)>tol))
}

# find straightforward redundancies of the form 0 <= 0 or 0 == 0
is_tautology <- function(A,b,tol){
  rowSums(abs(A) > tol) == 0 & abs(b) < tol
}

shrink <- function(A, b, neq, tol){
  i <- !is_tautology(A,b,tol)
  j <- colSums(abs(A)>tol) == 0 
  list(
      A   = A[i,!j,drop=FALSE]
    , b   = b[i]
    , neq = sum(which(i) < neq)
  )
}

is_feasible <- function(A, b, neq=nrow(A), tol=1e-8){
  if (nrow(A)==0) return(TRUE)
  if ( has_contradiction(A,b,tol) ) return(FALSE)
  
  bi <- block_index(A,tol = tol)
  # sort so smaller blocks are treated first:
  bi <- bi[order(sapply(bi,length))]
  
  blocks <- lapply(bi, function(i) shrink(A=A[i,,drop=FALSE], b=b[i], neq = sum(i <= neq), tol=tol))
  
  for ( sys in blocks ){
    var <- 1    
    L <- eliminate(sys$A, sys$b, neq = sys$neq, variable=var)
    L <- shrink(L$A, L$b, neq=L$neq, tol=tol)
    feasible <- is_feasible(A=L$A, b=L$b, neq=L$neq, tol=tol) 
    if (!feasible) break
  }
  feasible
}







