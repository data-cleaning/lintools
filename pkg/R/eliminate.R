
#' Eliminate a variable from a set of edit rules
#' 
#' Eliminating a variable amounts to deriving all (non-redundant) linear
#' (in)equations not containing that variable. Geometrically, it can be seen as
#' a projection of the solution space (records obeying all edits) along the
#' eliminated variable's axis. 
#' 
#' @section Details:
#' For equalities, Gaussian elimination is used. If inequalities are involved,
#' Fourier-Motzkin elimination is used. In principle, FM-elimination can
#' generate a large number of redundant inequations, especially when applied
#' recursively. Redundancies can be recognized by recording how new inequations
#' have been derived from the original set. This is stored in the \code{H} matrix
#' when multiple variables are to be eliminated (Kohler, 1967).
#' 
#' 
#'
#' @param A \code{[numeric]} Matrix 
#' @param b \code{[numeric]} vector
#' @param neq [\code{numeric}] The first \code{neq} rows in \code{A} and
#'   \code{b} are treated as linear equalities. The others as Linear
#'   inequalities of the form \eqn{Ax<=b}.
#' @param variable \code{[numeric|logical|character]} Index in columns of \code{A}, representing the variable to eliminate.
#' @param H \code{[numeric]} (optional) Matrix indicating how linear inequalities have been derived. 
#' @param h \code{[numeric]} (optional) number indicating how many variables have been eliminated from the original system
#' using Fourier-Motzkin elimination.
#' 
#'   
#' @export
#'
#' 
#' @return A \code{list} with the folowing components
#' \itemize{
#'   \item{\code{A}: the \code{A} corresponding to the system with variables eliminated.}
#'   \item{\code{b}: the constant vector corresponding to the resulting system}
#'   \item{\code{neq}: the number of equations}
#'   \item{\code{H}: The memory matrix storing how each row was derived}
#'   \item{\code{h}: The number of variables eliminated from the original system.}
#' }
#' 
#'
#' @references
#' D.A. Kohler (1967) Projections of convex polyhedral sets, Operational Research
#' Center Report , ORC 67-29, University of California, Berkely.
#' 
#' H.P. Williams (1986) Fourier's method of linear programming and its dual,
#' The American Mathematical Monthly 93, 681-695
#' @export
eliminate <- function(A, b, neq , variable, H=NULL, h=0){

    Ab <- cbind(A,b)
    if (is.character(variable)){
      var <- match(variable, colnames(A))[1]
    } else {
      var <- variable
    }
    
    
    
    ops <- rep("<=",nrow(A))
    ops[seq_len(neq)] <- "=="
   
    coefs <- Ab[,var]
    I <- coefs != 0
    
    eq <- I & seq_len(nrow(A)) <= neq
    
    upper <- which(!eq & coefs > 0)
    lower <- which(!eq & coefs < 0)
    
    coefs[!eq] <- abs(coefs[!eq])
    
    eq <- which(eq);

    # elimination possible?
    if ( (length(upper) > 0 && length(lower) > 0) ||
         (length(eq) >= 1 && (length(upper) > 0 || length(lower) > 0)) ||
         (length(eq) >= 2) ){
       h <- h+1
    } else {
      # return rows and columns where 'var' does not occur
      ii <- A[,var] == 0
      return(list(
        A = Ab[ii,-ncol(Ab), drop=FALSE]
        , b = b[ii]
        , neq = sum(ii<=neq)
        , H = H
        , h = h
      ))
    } 

    if ( is.null(H) ){
        H <- matrix(FALSE,nrow=nrow(A),ncol=nrow(A))
        diag(H) <- TRUE 
        colnames(A) <- rownames(A)
    }
    

    #normalize matrix, every row has coefficient 1 or -1 for var
    Ab[I,] <- Ab[I,] / coefs[I]

    # eqs and ineqs w/coef>0 ==> ineqs w/coef<0
    equpper <- c(eq, upper)
    I1 <- rep(equpper,each=length(lower))
    I2 <- rep(lower, times=length(equpper))
    ml <- Ab[I1,,drop=FALSE] + Ab[I2,,drop=FALSE]
    ol <- ifelse(ops[I1] != "<", ops[I2], ops[I1])
    dl <- H[I1,,drop=FALSE] | H[I2,,drop=FALSE]

    # eqs ==> ineqs w/coef>0
    I1 <- rep(eq,each=length(upper))
    I2 <- rep(upper, times=length(eq))
    mu <- Ab[I2,,drop=FALSE] - Ab[I1,,drop=FALSE]
    ou <- ops[I2]
    du <- H[I1,,drop=FALSE] | H[I2,,drop=FALSE]

    # eqs ==> eqs
    me <- Ab[logical(0),,drop=FALSE]
    de <- H[logical(0),,drop=FALSE]
    if ( length(eq)>1){
        me <- t(t(Ab[eq[-1],,drop=FALSE]) - Ab[eq[1],])
        de <- t(t(H[eq[-1],,drop=FALSE]) | H[eq[1],])       
    } 
    oe <- rep("==",nrow(me))

    Ab <- rbind(ml,mu,me,Ab[!I,,drop=FALSE])
    H <- rbind(dl,du,de,H[!I,,drop=FALSE])
    o <- c(ol,ou,oe,ops[!I])
    redundant <- rowSums(H) > h + 1 #| isObviouslyRedundant.matrix(E=m, operators=o)
    
    Ab <- Ab[!redundant,,drop=FALSE]
    H <- H[!redundant,,drop=FALSE]

    L <- normalize(
        A = Ab[,-ncol(Ab),drop=FALSE]
      , b = Ab[,ncol(Ab),drop=TRUE]
      , operators = o[!redundant]
    ) 
    list(A = L$A, b = L$b, neq=L$neq, H=H[L$order,,drop=FALSE], h=h)
}



#' Bring a system of (in)equalities in normal form
#' 
#' @section Details:
#' A set of equations is in normal form when the first \code{neq} entries are 
#' linear equalities and all inequalities are weak inequalities of the form
#' \code{a.x <= b}. To bring a system in normal form, the following steps are 
#' performed:
#' \itemize{
#'   \item{Equations of the form \code{a.x >= b} are replaced with \code{-a.x <= -b}}
#'   \item{Equations of the form \code{a.x > b} are replaced with \code{-a.x <= -b-unit}}
#'   \item{Equations of the form \code{a.x < b} are replaced with \code{a.x <= b-unit}}
#'   \item{Equalities are sorted on top}
#' }
#' 
#' 
#' 
#' 
#' @param A \code{[numeric]} Matrix
#' @param b \code{[numeric]} vector
#' @param operators \code{[character]} operators in \code{{<,<=,==,>=,>}}.
#' @param unit \code{[numeric]} (positive) Your unit of measurement. This is used to
#' replace strict inequations of the form \code{a.x < b} with \code{a.x <= b-unit}.
#' Typically, \code{unit} is related to the units in which your data 
#' is measured.  
#' 
#' @return A \code{list} with the folowing components
#' \itemize{
#'   \item{\code{A}: the \code{A} corresponding to the normalized sytem.}
#'   \item{\code{b}: the constant vector corresponding to the normalized system}
#'   \item{\code{neq}: the number of equations}
#'   \item{\code{order}: the index vector used to permute the original rows of \code{A}.}
#' }
#' 
#' 
#' 
#' @export
normalize <- function(A, b, operators, unit=1 ){
  
  geq <- operators  == ">="
  gt <- operators == ">"
  A[geq|gt,] <- -A[geq|gt,,drop=FALSE]
  b[geq] <- -b[geq]
  b[gt] <- -b[gt] - unit
  operators[geq|gt] <- "<="
 
  lt <- operators == "<"
  b[lt] <- b[lt] - unit
  operators[lt] <- "<="
  
   
  ii <- order(operators,decreasing=TRUE)
    
  list(A = A[ii,,drop=FALSE], b=b[ii], neq=sum(operators=="=="), order=ii)
}




