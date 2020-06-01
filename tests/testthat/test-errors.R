context('Test vglmer *fails* in certain circumstances')


test_that('vglmer fails with incorrect initial options', {
  
  
  N <- 1000
  G <- 100
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  
  expect_error(
    vglmer(formula = y ~ x + (1 | g), data = NULL, iterations = 5000, 
   print_prog = 1000, tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
   family = 'linear', prior_variance = 'mean_exists', init = 'zero', factorization_method = v)
  )

  expect_error(
    vglmer(formula = y ~ x + (1 | g), data = NULL, iterations = 5000, 
           print_prog = 1000, tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
           family = 'logit', prior_variance = 'mean_exists', init = 'other', factorization_method = v)
  )

  expect_error(
    vglmer(formula = y ~ x + (1 | g), data = NULL, iterations = 5000, 
           print_prog = 1000, tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
           family = 'logit', prior_variance = 'random_prior', init = 'EM', factorization_method = v)
  )
  
  
  
})