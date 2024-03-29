

## compactify linear systems

## column removal, no x
  L <- compact(
    A = matrix(c(1,0),nrow=1)
    ,b = 1
    ,neq=1
    ,nleq=0
  )
  expect_equivalent(L$A,matrix(1))
  expect_equal(L$b,1)
  expect_equal(L$neq,1)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,TRUE))
  

## neq is computed correctly in case of duplicate rows
  # -x <= -1
  #  x <=  1
  #  x <=  1
  # conclusion: x == 1.
  L <- compact(
      A    = matrix(c(-1,1,1), nrow=3)
    , b    = matrix(c(-1,1,1), nrow=3)
    , neq  = 0
    , nleq = 3
  )
  expect_equal(L$neq,1)


## column removal, with x
  L <- compact(
    A = matrix(c(1,0),nrow=1)
    , x = c(2,8)
    ,b = 1
    ,neq=1
    ,nleq=0
  )
  expect_equivalent(L$A,matrix(1))
  expect_equal(L$b,1)
  expect_equal(L$x,2)
  expect_equal(L$neq,1)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,TRUE))


## row removal
  # x + y == 1
  # 0 + 0 == 0
  L <- compact(
    A = matrix(c(1,1,0,0),nrow=2,byrow=TRUE)
    , b = c(1,0)
    , neq=2
  )
  expect_equivalent(L$A,matrix(c(1,1),nrow=1))
  expect_equal(L$b,1)
  expect_equal(L$neq,1)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,FALSE))

## row removal, case w/nonzero b
  L <- compact(
    A = matrix(c(1,1,0,0),nrow=2,byrow=TRUE)
    , b = c(1,2)
    , neq=2
  )
  expect_equivalent(L$A,matrix(c(1,1,0,0),nrow=2,byrow=TRUE))
  expect_equal(L$b,c(1,2))
  expect_equal(L$neq,2)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,FALSE))


## Combine inequalities, simple case
  # x <= 0
  # x >= 0
  L <- compact(
    A=matrix(c(1,-1))
    , b=c(0,0)
    , neq=0
    , nleq=2
    )
 expect_equivalent(L$A, matrix(1,nrow=1))
 expect_equal(L$b,0)
 expect_equal(L$neq,1)
 expect_equal(L$nleq,0)
 expect_equal(L$cols_removed,FALSE)

## Combine inequalities, simple case with non-zero b
  # x <= 3
  # x >= 3
  L <- compact(
    A=matrix(c(1,-1))
    , b=c(3,-3)
    , neq=0
    , nleq=2
    )
 expect_equivalent(L$A, matrix(1,nrow=1))
 expect_equal(L$b,3)
 expect_equal(L$neq,1)
 expect_equal(L$nleq,0)
 expect_equal(L$cols_removed,FALSE)


## Combine inequalities, simple case that includes equalities
  # x + y == 1
  # x <= 3
  # x >= 3
  L <- compact(
    A=matrix(c(1,1, 1,0, -1,0),nrow=3,byrow=TRUE)
    , b=c(1,3,-3)
    , neq=1
    , nleq=2
    )
 expect_equivalent(L$A, matrix(c(1,1,1,0),nrow=2,byrow=TRUE))
 expect_equal(L$b,c(1,3))
 expect_equal(L$neq,2)
 expect_equal(L$nleq,0)
 expect_equal(L$cols_removed,c(FALSE,FALSE))


##  combined inequalities and row removal
  #  x + y ==  1
  #  x + 0 <=  1
  # -x + 0 <= -1
  #  0 + 0 <=  0
  L <- compact(
    A = matrix(c(1,1, 1,0, -1,0, 0,0),nrow=4,byrow=TRUE)
    , b = c(1,1,-1,0)
    , neq = 1
    , nleq = 3
  )
  expect_equivalent(L$A, matrix(c(1,1,1,0),nrow=2,byrow=TRUE))
  expect_equal(L$neq,2)
  expect_equal(L$nleq,0)
  expect_equal(L$cols_removed,c(FALSE,FALSE))

## earlier bugs
  
  A <- matrix(c(1,1,-1,0,0,-1),byrow=TRUE,nrow=3)
  b <- rep(0,3)
  expect_equal(compact(A=A,b=b,neq=1,nleq=2,eps=1e-8)$neq, 1)
  
  A <- structure(c(0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 
        -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        , .Dim = c(6L, 6L))
  b <- c(1130, 0, 75, 0, 0, 18915)
  expect_equal(compact(A,b,neq=2,nleq=4)$neq,1)
  

## implied equations in presence of equalities
  # a + b - c == 0
  #         c <= 1
  #        -c <= -1
  A <- matrix(c(1,1,-1, 0,0,1, 0,0,-1),nrow=3,byrow = TRUE)
  b <- c(0,1,-1)
  L <- compact(A,b,neq=1,nleq=2)
  expect_equal(L$neq,2)
  # one redundant equality remains
  expect_equal(L$A, matrix(c(1,1,-1,0,0,1),nrow=2,byrow=TRUE))
  

## implied equations, nonzero coefficients.
  #   x <=  0
  #  -x >=  0
  #   y <=  1
  #  -y <= -1
  A <- matrix(c(1,0, -1,0, 0,1, 0,-1),nrow=4,byrow=TRUE)
  b <- c(0,0,1,-1) 
  L <- compact(A,b,neq=0, nleq=4)  
  expect_equal(L$neq,2)
  expect_equal(L$A, diag(2))
  expect_equal(L$b, c(0,1))

## implied equalities, remaining inequation, and remaining strict inequation
 #  x+z == -2
 #  x   <=  1
 # -x   <= -1
 #  y   <=  8
 #  z   <   0
 A <- matrix(c(1,0,1, 1,0, 0, -1,0,0, 0,1,0, 0,0,1), nrow=5,byrow=TRUE)
 b <- c(-2, 1, -1, 8, 0)
 L <- compact(A,b, neq=1,nleq=3) 
 expect_equal(L$neq,2)
 expect_equal(L$nleq,1)
 expect_equal(L$A, A[c(1,2,4,5),,drop=FALSE])
 expect_equal(L$b, b[c(1,2,4,5)])


## implicit duplicates
#   x <=  1
#  -x <= -1
# -2x <= -2
A <- matrix(c(1,-1,-2),nrow=3)
b <- c(1,-1,-2)
L <- compact(A,b,neq=0,nleq=3)
expect_equal(L$neq,1)
expect_equal(L$nleq,0)

## regression test (<=0.1.5)
#  x  <= 1
# -y  <= 1 
A <- matrix(c(1,0, 0,-1), nrow=2, byrow=TRUE)
b <- c(1,1)
L <- compact(A, b, neq=0, nleq=2)
expect_equal(compact(L$A, L$b,neq=0, nleq=2), L)








