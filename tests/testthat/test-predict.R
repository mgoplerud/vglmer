context('Test Predict')

test_that('Prediction Matches Manual and (nearly) glmer', {
  
  N <- 1000
  G <- 10
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g]))
  
  est_glmer <- suppressWarnings(lme4::glmer(y ~ x + (1 + x | g) + (1 | g2), family = binomial))
  
  example_vglmer <- vglmer(formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, iterations = 5000, 
     print_prog = 1000, tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
     family = 'logit', prior_variance = 'mean_exists', init = 'zero', factorization_method = 'strong')
  
  glmer_predict <- predict(est_glmer)
  def_predict <- predict(example_vglmer, newdata = data.frame(y =y, x = x, g = g))
  
  expect_gt(
    cor(def_predict, glmer_predict), 0.99
  )
  
  alpha_names <- rownames(example_vglmer$alpha$mean)
  
  manual_predict <- as.vector(
    example_vglmer$alpha$mean[match(paste0('g2 @ (Intercept) @ ', g2), alpha_names)]  +
    example_vglmer$alpha$mean[match(paste0('g @ x @ ', g), alpha_names)] * x +
    example_vglmer$alpha$mean[match(paste0('g @ (Intercept) @ ', g), alpha_names)] +
    cbind(1, x) %*% example_vglmer$beta$mean  
  )
  
  expect_equal(def_predict, manual_predict)

})

test_that('Prediction Matches for New Levels in newdata', {
  
  N <- 50
  G <- 10
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  
  example_vglmer <- vglmer(formula = y ~ x + (1 | g), data = NULL, iterations = 2,
                           print_prog = 1000, tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
                           family = 'logit', prior_variance = 'mean_exists', init = 'zero', factorization_method = 'weak')
  #No old level in new level
  new_data <- data.frame(x = rep(0,5), g = 11)  
  expect_error(predict(example_vglmer, newdata = new_data)) #Error on default
  expect_equal(#zero when allow_missing_levels = TRUE
    predict(example_vglmer, newdata = new_data, allow_missing_levels = TRUE),
    rep(example_vglmer$beta$mean[1], nrow(new_data))
  )
  mixed_data <- data.frame(x = rnorm(10), g = sample(1:25, 10, replace = T))
  
  man_beta <- as.vector(cbind(1, mixed_data$x) %*% example_vglmer$beta$mean)
  man_alpha <- example_vglmer$alpha$mean[match(paste0('g @ (Intercept) @ ', mixed_data$g), rownames(example_vglmer$alpha$mean))]
  man_alpha[is.na(man_alpha)] <- 0
  expect_equivalent(man_beta + man_alpha, predict(example_vglmer, newdata = mixed_data, allow_missing_levels = T))

})

test_that('Prediction Matches for Missing in new.data', {
  
  N <- 50
  G <- 10
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  
  example_vglmer <- vglmer(formula = y ~ x + (1 | g), data = NULL, iterations = 2, 
                           print_prog = 1000, tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
                           family = 'logit', prior_variance = 'mean_exists', init = 'zero', factorization_method = 'strong')
  
  mixed_data <- data.frame(x = rnorm(20), g = rep(1:10,2))
  rownames(mixed_data) <- letters[1:20]
  
  mixed_data$x[8] <- NA
  mixed_data$g[7] <- NA
  mixed_data$x[2] <- mixed_data$g[2] <- NA
  
  man_beta <- as.vector(cbind(1, mixed_data$x) %*% example_vglmer$beta$mean)
  man_alpha <- example_vglmer$alpha$mean[match(paste0('g @ (Intercept) @ ', mixed_data$g), rownames(example_vglmer$alpha$mean))]
  expect_equivalent(man_beta + man_alpha, 
                    predict(example_vglmer, newdata = mixed_data))
  
})

test_that('Prediction Matches for Simulation', {
  
})

