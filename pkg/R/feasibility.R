
# find straigtforward contradictions of the form 0 <= b or 0 == b 
has_contradiction <- function(A,b, neq, nleq, eps){
  if (nrow(A)==0) return(FALSE)
  if (ncol(A)==0 & all(abs(b) > eps)) return(TRUE)

  ieq <- seq_len(neq)
  leq <- neq + seq_len(nleq)
  lt <- neq + nleq + seq_len(nrow(A)-neq-nleq)
  # nonzero entries
  AI <- rowSums(abs(A) > eps) == 0
  any( AI[ieq] & abs(b[ieq]) > eps) ||
    any( AI[leq] & b[leq] <= -eps) ||
    any( AI[lt] & b[lt] <= 0 )
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
#' @export
#' 
#' @examples 
#' # x + y == 0
#' # x > 0
#' # y > 0
#' A <- matrix(c(1,1,1,0,0,1),byrow=TRUE,nrow=3)
#' b <- rep(0,3)
#' is_feasible(A=A,b=b,neq=1,nleq=0)
#' 
is_feasible <- function(A, b, neq=nrow(A), nleq=0, eps=1e-8, method="elimination"){
  if (nrow(A)==0) return(TRUE)
  if ( has_contradiction(A=A,b=b,neq=neq,nleq=nleq,eps=eps) ) return(FALSE)
  
  bi <- block_index(A,eps = eps)
  # sort so smaller blocks are treated first:
  bi <- bi[order(sapply(bi,length))]
 
  ## TODO: all sorts of optimizations, including:
  # - use H-matrix from elimination step
  # - check singularity of A'A of equality section
  # - figure out a good variable elimination order
   
  for ( ii in bi ){
    L <- compact(
      A[ii,,drop=FALSE]
      , b=b[ii]
      , neq=sum(ii<=neq)
      , nleq= sum(ii>neq & ii <= nleq)
      , remove_columns=FALSE
      )
    L <- eliminate(L$A, L$b, neq = L$neq, variable=1)
    feasible <- is_feasible(A=L$A, b=L$b, neq=L$neq,nleq=L$nleq, eps=eps) 
    if (!feasible) break
  }
  feasible
}







