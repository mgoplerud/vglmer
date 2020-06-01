context('Match glmer')

test_that('Compare against glmer', {
  
  N <- 1000
  G <- 100
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  
  est_glmer <- lme4::glmer(y ~ x + (1 | g), family = binomial)
  fmt_glmer <- format_glmer(est_glmer)
  
  for (v in c("weak", 'partial', 'strong')){
    
    example_vglmer <- vglmer(formula = y ~ x + (1 | g), data = NULL, iterations = 5000, 
                             print_prog = 1000, tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
                             family = 'logit', prior_variance = 'mean_exists', init = 'zero', factorization_method = v)
    
    fmt_vglmer <- format_vglmer(example_vglmer)
    comp_methods <- merge(fmt_glmer, fmt_vglmer, by = c('name'))
    
    cor_mean <- with(comp_methods, cor(mean.x, mean.y))
    expect_gt(cor_mean, expected = 0.99)
    

    
    
  }
  
  
})