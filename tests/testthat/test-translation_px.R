context("Translation Expansion Tests")

test_that('Translation Data Test', {

  N <- 100
  G <- 5
  x <- rnorm(N)
  x2 <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  
  for (v in c('dynamic', 'numerical', 'OSL')){
    
    est_vglmer <- vglmer(y ~ x + x2 + (1 + x | g) + (1 + x2 | g2), data = NULL,
                         family = 'binomial',
                         control = vglmer_control(px_method = v, debug_px = TRUE))
    
    expect_gt(min(diff(ELBO(est_vglmer, 'traj'))), -sqrt(.Machine$double.eps))
  }

})
  

