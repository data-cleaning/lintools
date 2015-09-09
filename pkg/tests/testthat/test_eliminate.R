context("variable elimination")




test_that("Gaussian elimination",{
  A <- matrix(c(
    1,1,-1
    ,0,1,-1),byrow=TRUE,nrow=2
  )
  b <- c(1,0)
  neq=2
  
  L <- eliminate(A=A,b=b,neq=neq,variable=1)
  expect_equivalent(L$A,matrix(c(0,1,-1)))
  expect_equivalent(L$b,0)
  expect_equal(L$neq,2)
  expect_equal(L$H,NULL)
  expect_equal(L$h,0)
})


test_that("Fourier-Motzkin elimination",{
  
  A <- matrix(c(
    4, -5, -3,  1,
   -1,  1, -1,  0,
    1,  1,  2,  0,
   -1,  0,  0,  0,
    0, -1,  0,  0,
    0,  0, -1,  0),byrow=TRUE,nrow=6) 
  
  b <- c(0,2,3,0,0,0)
  L <- eliminate(A=A,b=b,neq=0,variable=1)
 
  expect_equivalent(L$A,
    matrix(c(
      0, -0.25, -1.75, 0.25,
      0, -1.25, -0.75, 0.25,
      0,  2.00,  1.00, 0.00,
      0,  1.00,  2.00, 0.00,
      0, -1.00,  0.00, 0.00,
      0,  0.00, -1.00, 0.00), byrow=TRUE,nrow=6)
  )
  expect_equivalent(L$b,c(2,0,5,3,0,0))
  expect_equivalent(L$H,
    matrix(c(
      TRUE,  TRUE, FALSE, FALSE, FALSE, FALSE,
      TRUE, FALSE, FALSE,  TRUE, FALSE, FALSE,
      FALSE,  TRUE,  TRUE, FALSE, FALSE, FALSE,
      FALSE, FALSE,  TRUE,  TRUE, FALSE, FALSE,
      FALSE, FALSE, FALSE, FALSE,  TRUE, FALSE,
      FALSE, FALSE, FALSE, FALSE, FALSE,  TRUE), byrow=TRUE,nrow=6)
  )
  expect_equal(L$h,1)
})



#  
#   P <- editrules::editmatrix(c(
#     "4*x1 - 5*x2 - 3*x3 + z <= 0",
#     "-x1 + x2 -x3 <= 2",
#     "x1 + x2 + 2*x3 <= 3",
#     "-x1 <= 0",
#     "-x2 <= 0",
#     "-x3 <= 0"))
#   P1 <- editrules::eliminate(P, "x1", fancynames=TRUE)
#   matrix(c(
#     TRUE,  TRUE, FALSE, FALSE, FALSE, FALSE,
#     TRUE, FALSE, FALSE,  TRUE, FALSE, FALSE,
#     FALSE,  TRUE,  TRUE, FALSE, FALSE, FALSE,
#     FALSE, FALSE,  TRUE,  TRUE, FALSE, FALSE,
#     FALSE, FALSE, FALSE, FALSE,  TRUE, FALSE,
#     FALSE, FALSE, FALSE, FALSE, FALSE,  TRUE), byrow=TRUE,nrow=6)
#   op <- c("<=", "<=", "<=", "<=", "<=", "<=")
#   expect_true(all( Ab == getAb(P1) ))
#   expect_true(all( H  == getH(P1)  ))
#   expect_true(all( op == getOps(P1)))
#   expect_true( geth(P1) == 1 )
# })
# 
# 
# 
# e <- editrules::editmatrix(expression(x+y==1-z,y==z))
# editrules:::getH(editrules::eliminate(e,'x'))
# 
# e1 <- editrules::editmatrix(expression(x + y > z, x + z > w))
# editrules:::getH(editrules::eliminate(e1,'x'))
