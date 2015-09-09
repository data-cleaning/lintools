
# find straigtforward contradictions of the form 0 <= b or 0 == b 
has_contradiction <- function(A,b,tol){
  if (nrow(A)==0) return(FALSE)
  any(rowSums(abs(A)>tol) == 0 & (abs(b)>tol))
}

# find straightforward redundancies of the form 0 <= 0 or 0 == 0
is_tautology <- function(A,b,tol){
  rowSums(abs(A) > tol) == 0 & abs(b) < tol
}

shrink <- function(A, b, neq, variable=0, tol){
  i <- !is_tautology(A,b,tol)
  list(
      A   = A[i,-variable,drop=FALSE]
    , b   = b[i]
    , neq = sum(which(i) < neq)
  )
}

is_feasible <- function(A, b, neq=nrow(A), tol=1e-8){
  bi <- block_index(A,tol = tol)
  blocks <- lapply(bi,function(i) list(A=A[i,,drop=FALSE],b=b[i], neq = sum(i < neq) ))
  !any(sapply(blocks, function(sys){
    var <- 1    
    L <- eliminate(sys$A, sys$b, neq = sys$neq,variable=var)
    L <- shrink(L$A, L$b, neq=L$neq, variable=var, tol=tol)
    has_contradiction(A=L$A, b=L$b, tol=tol) || !is_feasible(A=L$A, b=L$b, neq=L$neq, tol=tol) 
  }))
}










