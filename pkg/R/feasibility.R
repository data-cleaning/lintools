
# find straigtforward contradictions of the form 0 <= b or 0 == b 
has_contradiction <- function(A,b,eps){
  if (nrow(A)==0) return(FALSE)
  if (ncol(A)==0 & all(abs(b) > eps)) return(TRUE)
  any(rowSums(abs(A)>eps) == 0 & (abs(b)>eps))
}

# find straightforward redundancies of the form 0 <= 0 or 0 == 0
is_tauepsogy <- function(A,b,eps){
  rowSums(abs(A) > eps) == 0 & abs(b) < eps
}


#' Check feasibility of a system of linear (in)equations
#'
#' @param A [\code{numeric}] matrix
#' @param b [\code{numeric}] vector
#' @param neq [\code{numeric}] The first \code{neq} rows in \code{A} and
#'   \code{b} are treated as linear equalities. 
#' @param nleq [\code{numeric}] The \code{nleq} rows after \code{neq} are treated as
#'   inequations of the form \code{a.x<=b}. All remaining rows are treated as strict inequations
#'   of the form \code{a.x<b}.
#' @param eps [\code{numeric}] Absolute values \code{< eps} are treated as zero.
#' @param method [\code{character}] At the moment, only the 'elimination' method is implemented.
#'
#'
#' @export
is_feasible <- function(A, b, neq=nrow(A), nleq=0, eps=1e-8, method="elimination"){
  if (nrow(A)==0) return(TRUE)
  if ( has_contradiction(A,b,eps) ) return(FALSE)
  
  bi <- block_index(A,eps = eps)
  # sort so smaller blocks are treated first:
  bi <- bi[order(sapply(bi,length))]
 
  ## TODO: all sorts of optimizations, including:
  # - use H-matrix from elimination step
  # - check singularity of A'A of equality section
  # - figure out a good variable elimination order
   
  for ( ii in bi ){
    block <- compact(
      A[ii,,drop=FALSE]
      , b=b[ii]
      , neq=sum(ii<=neq)
      , nleq= sum(ii>neq & ii <= nleq)
      )
    L <- eliminate(block$A, block$b, neq = block$neq, variable=1)
    L <- compact(L$A, L$b, neq=L$neq, nleq=L$nleq, eps=eps)
    feasible <- is_feasible(A=L$A, b=L$b, neq=L$neq, eps=eps) 
    if (!feasible) break
  }
  feasible
}







