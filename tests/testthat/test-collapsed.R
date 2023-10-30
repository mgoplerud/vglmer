
test_that("Test block diagonal functions", {
  
  N <- 100
  P <- 5
  A <- matrix(rnorm(N * P), ncol = P)
  B <- matrix(rnorm(N * P), ncol = P)
  
  est <- block_diag_product(A = A, B = B, block_size = 2, blocks = N/2)  
  direct <- A %*% t(B)
  
  manual_2 <- t(sapply(1:(N/2), FUN=function(i){id <- 1 + (i-1) * 2 + 0:1; return(as.vector(direct[id,id]))}))
  expect_equivalent(est, manual_2)  
  
  est <- block_diag_product(A = A, B = B, block_size = 1, blocks = N)  
  
  expect_equivalent(est, diag(direct))

  BlockA <- tcrossprod(A)  
  
  two_block <- block_diag_sum(A = BlockA, d = 2, g = N/2)
  manual_block_A <- Reduce('+', lapply(1:(N/2), FUN=function(i){BlockA[1 + 2 * (i-1) + 0:1 , 1 + 2 * (i-1) + 0:1]}))
  expect_equivalent(two_block, manual_block_A)

})
