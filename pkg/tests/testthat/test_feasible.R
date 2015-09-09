

context("feasibility")

test_that("is_feasible",{
  # x + y == 1 & x + y == 2
  expect_false( is_feasible(A = matrix(rep(1,4),nrow=2), b = c(1,2), neq = 2) )
  # x >= y & x <= y - 1 
  expect_false( is_feasible(A = matrix(c(-1,1,1,-1),byrow=TRUE,nrow=2), b = c(0,-1), neq = 0) )
})


