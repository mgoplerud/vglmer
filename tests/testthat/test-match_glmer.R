context('Match glmer')

test_that('Compare against glmer', {
  
  N <- 1000
  G <- 100
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  
  est_glmer <- suppressWarnings(lme4::glmer(y ~ x + (1 | g), family = binomial))
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


test_that('Compare against glmer.nb', {
  
  skip_on_cran()
  
  N <- 1000
  G <- 50
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  alpha <- rnorm(G)
  
  y <- rnbinom(n = N, mu = exp(-1 + x + alpha[g]), size = 5)
  data <- data.frame(y = y, x = x, g = g)
  est_glmer <- suppressWarnings(glmer.nb(y ~ x + (1 | g), data = data, family = binomial))
  fmt_glmer <- format_glmer(est_glmer)
  
  for (v in c("weak", 'partial', 'strong')){
    
    example_vglmer <- vglmer(formula = y ~ x + (1 | g), data = data, iterations = 5000, 
                             print_prog = 1000, tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
                             family = 'negbin', prior_variance = 'mean_exists', init = 'zero', factorization_method = v)
    
    fmt_vglmer <- format_vglmer(example_vglmer)
    comp_methods <- merge(fmt_glmer, fmt_vglmer, by = c('name'))
    
    cor_mean <- with(comp_methods, cor(mean.x, mean.y))
    expect_gt(cor_mean, expected = 0.99)
    
  }
  
})

test_that('EM_prelim matches glm', {
  N <- 100
  p <- 10
  
  Z <- matrix(rnorm(N * p), ncol = p)  
  beta <- rnorm(p)
  
  y <- rbinom(N, 1, plogis(Z %*% beta))
  
  est_glm <- glm(y ~ Z, family = binomial)
  est_init <- EM_prelim(X = drop0(matrix(1, nrow = N)), Z = drop0(Z), s = y - 1/2, pg_b = 1, iter = 200, ridge = Inf)  
  est_init <- c(est_init$beta, est_init$alpha)
  expect_equal(as.vector(coef(est_glm)), est_init, tolerance = 1e-4)  
})

