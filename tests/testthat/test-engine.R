for(seed in 1:10){
 test_that(sprintf("Heap test for seed %d",seed),{
  set.seed(seed)
  
  N<-sample(100,1)
  M<-N+10+sample(100,1) # We need M (len(B)) > N (len(A))

  A<-runif(N)
  B<-runif(M)
  B[1:N]<-pmax(A,B[1:N]) # We need B>A for 1..N

  .Call(C_heapTest,A,B,TRUE)->C

  C_exp<-c(
   sort(A,decreasing=TRUE),
   sort(B,decreasing=TRUE)
  )

  expect_equal(C,C_exp)
 })
}

test_that("Coercing works as expected",{
 expect_equal(vistla_coerce(1:3),factor(c("l1","l2","l2")))
 expect_equal(vistla_coerce(factor(1:3)),factor(c("l1","l2","l3")))
 expect_equal(vistla_coerce(pi+1:3),factor(c("l1","l1","l2")))
 expect_equal(vistla_coerce(1:2),factor(c("l1","l2")))
 expect_equal(vistla_coerce(c(pi,pi,pi)),factor(c("l1","l1","l1")))
 expect_equal(vistla_coerce(c(1L,1L,1L)),factor(c("l1","l1","l1")))
 expect_equal(vistla_coerce(c(pi,pi)),factor(c("l1","l1")))
 expect_equal(
  length(levels(vistla_coerce(1:300))),
  10L
 )
 expect_equal(
  length(levels(vistla_coerce(1:12))),
  4L
 )
 expect_equal(
  length(levels(vistla_coerce(1:12+pi))),
  4L
 )
 expect_error(vistla_coerce(c(1,NA,2)),"Non-finite values and NAs are not accepted")
 expect_error(vistla_coerce(factor(c(1,NA,2))),"NAs are not accepted")
 expect_error(vistla_coerce(c(1L,NA,2L)),"NAs are not accepted")
 expect_error(vistla_coerce(sin),"Invalid input")
})
