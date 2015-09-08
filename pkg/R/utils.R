
all_finite <- function(x){
  stopifnot(is.numeric(x))
  storage.mode(x) <- "double"
  .Call("all_finite_double",x)
}

