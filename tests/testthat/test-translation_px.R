context("Translation Expansion Tests")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(54321)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

print(paste0(NITER, ' for tests because ', env_test))

test_that('Translation Data Test', {

  N <- 100
  G <- 5
  x <- rnorm(N)
  x2 <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  
  for (v in c('dynamic', 'numerical', 'OSL', 'profiled')){

    est_vglmer <- vglmer(y ~ x + x2 + (1 + x | g) + (1 + x2 | g2), data = NULL,
         family = 'binomial',
         control = vglmer_control(iterations = NITER, px_method = v, debug_px = TRUE))
    
    expect_gt(min(diff(ELBO(est_vglmer, 'traj'))), -sqrt(.Machine$double.eps))
  }

})
  
test_that("Check that B_j has correct shape", {
  N <- 100
  G <- 5
  x <- rnorm(N)
  x2 <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  
  est_vglmer <- vglmer(y ~ v_s(x) + x2 + (1 + x | g) + (1 + x + x2 | g2), data = NULL,
                       family = 'linear',
                       control = vglmer_control(iterations = NITER))
  expect_gt(min(diff(ELBO(est_vglmer, 'traj'))), -sqrt(.Machine$double.eps))

  est_vglmer <- vglmer(y ~ v_s(x) + (1 + x | g) + (1 + x2 | g2), data = NULL,
                       family = 'linear',
                       control = vglmer_control(iterations = NITER))
  expect_gt(min(diff(ELBO(est_vglmer, 'traj'))), -sqrt(.Machine$double.eps))

})


