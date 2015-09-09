
all_finite <- function(x){
  stopifnot(is.numeric(x))
  storage.mode(x) <- "double"
  .Call("all_finite_double",x)
}


check_sys <- function(A, b, neq, x ,tol){
  if (!is.numeric(A)){
    stop("'A' is not a numeric matrix")
  }
  if (!is.numeric(b)){
    stop("'b' is not a numeric vector")
  }
  if (!missing(x) && !is.numeric(x)){
    stop("'x' is not a numeric vector")
  }
  if ( !is.matrix(A)){
    stop("'A' is not of class 'matrix'")
  }
  if (!missing(tol) && tol < 0){
    stop("Negative value for 'tol' specified")
  }
  if (!all_finite(A)){
    stop("Non-finite values detected in 'A'")
  }
  if (!all_finite(b)){
    stop("Non-finite values detected in 'b'")
  }
  if (!missing(x) && !all_finite(x)){
    stop("Non-finite values detected in 'x'")
  }
  if (!missing(tol) && !all_finite(tol)){
    stop("Non-finite value specified for 'tol'")
  }
  if (nrow(A) != length(b)){
    stop("Number of rows of 'A' not equal to length of 'b'")
  }
  if (!missing(x) && ncol(A) != length(x)){
    stop("Number of columns of 'A' not equal to length of 'x'")
  }
  if (!missing(neq) && neq > nrow(A)){
    stop("number of equations 'neq' larger than number of rows in matrix 'A'")
  }
  if (!missing(neq) && neq > length(b)){
    stop("number of equations 'neq' larger than length of 'b'")
  }
  
}

