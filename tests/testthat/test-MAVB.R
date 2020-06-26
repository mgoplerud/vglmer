context('Test MAVB')

test_that('Short test of MAVB running for CRAN', {
  
  N <- 1000
  G <- 10
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g]))
  
  est_glmer <- suppressMessages(suppressWarnings(lme4::glmer(y ~ x + (1 + x | g) + (1 | g2), family = binomial)))
  
  example_vglmer <- vglmer(formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = 'binomial',
                           control = vglmer_control(factorization_method = 'weak'))
  
  mavb_samples <- tryCatch(MAVB(object = example_vglmer, samples = 10), 
                           error = function(e){NULL})
  expect_false(is.null(mavb_samples))  
})

test_that('Longer Test of MAVB', {
  
  skip_on_cran()
  
  N <- 1000
  G <- 10
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g]))
  
  est_glmer <- suppressMessages(suppressWarnings(lme4::glmer(y ~ x + (1 + x | g) + (1 | g2), family = binomial)))
  glmer_var <- format_glmer(est_glmer)$var
  example_vglmer <- vglmer(formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = 'binomial',
                           control = vglmer_control(factorization_method = 'strong'))
  example_vglmer_weak <- vglmer(formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = 'binomial',
                           control = vglmer_control(factorization_method = 'weak'))
  
  mavb_samples <- MAVB(object = example_vglmer, samples = 2000)
  
  mavb_colvar <- colVar(mavb_samples)
  mavb_colmeans <- colMeans(mavb_samples)
  
  vi_means <- as.vector(rbind(example_vglmer$beta$mean, example_vglmer$alpha$mean))  
  vi_var <- as.vector(c(diag(vcov(example_vglmer)), example_vglmer$alpha$dia.var))  
  vi_var_weak <- as.vector(c(diag(vcov(example_vglmer_weak)), example_vglmer_weak$alpha$dia.var))  
  
  #Mean should be closer
  expect_lte(max(abs(vi_means - mavb_colmeans)), .1)
  #Variances should be larger
  expect_gte(min(mavb_colvar - vi_var), 0)
  # and closer to weak (Scheme III)
  expect_gte(median(abs(vi_var_weak - vi_var)), median(abs(vi_var_weak - mavb_colvar)))
  
})

