context("Test Predict")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(171)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that("Prediction Matches Manual and (nearly) glmer", {
  
  if (env_test == 'local'){
    N <- 1000
  }else{
    N <- 50
  }
  G <- 3
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g]))

  est_glmer <- suppressMessages(suppressWarnings(lme4::glmer(y ~ x + (1 + x | g) + (1 | g2), family = binomial)))

  example_vglmer <- vglmer(
    formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = "binomial",
    control = vglmer_control(factorization_method = "weak")
  )
  draw_samples <- posterior_samples.vglmer(example_vglmer, samples = 10)
  draw_MAVB <- MAVB(example_vglmer, samples = 10, var_px = 1)
  expect_true(is.matrix(draw_MAVB))
  expect_true(is.matrix(draw_samples))
  
  def_predict <- predict(example_vglmer, newdata = data.frame(y = y, x = x, g = g, g2 = g2))
  
  if (env_test == 'local'){
    glmer_predict <- predict(est_glmer)
    
    expect_gt(
      cor(def_predict, glmer_predict), 0.95
    )
  }

  alpha_names <- rownames(example_vglmer$alpha$mean)
  manual_predict <- as.vector(
    example_vglmer$alpha$mean[match(paste0("g2 @ (Intercept) @ ", g2), alpha_names)] +
      example_vglmer$alpha$mean[match(paste0("g @ x @ ", g), alpha_names)] * x +
      example_vglmer$alpha$mean[match(paste0("g @ (Intercept) @ ", g), alpha_names)] +
      cbind(1, x) %*% example_vglmer$beta$mean
  )
  expect_equal(def_predict, manual_predict)
})

test_that("Prediction Matches for New Levels in newdata", {
  N <- 50
  G <- 10
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))

  example_vglmer <- vglmer(
    formula = y ~ x + (1 | g), data = NULL,
    control = vglmer_control(iterations = 2, print_prog = 10, init = "zero"),
    family = "binomial"
  )
  # No old level in new level
  new_data <- data.frame(x = rep(0, 5), g = 11)
  expect_error(predict(example_vglmer, newdata = new_data)) # Error on default
  expect_equal( # zero when allow_missing_levels = TRUE
    predict(example_vglmer, newdata = new_data, allow_missing_levels = TRUE),
    rep(example_vglmer$beta$mean[1], nrow(new_data))
  )
  mixed_data <- data.frame(x = rnorm(10), g = sample(1:25, 10, replace = T))

  man_beta <- as.vector(cbind(1, mixed_data$x) %*% example_vglmer$beta$mean)
  man_alpha <- as.vector(example_vglmer$alpha$mean)[match(paste0("g @ (Intercept) @ ", mixed_data$g), rownames(example_vglmer$alpha$mean))]
  man_alpha[is.na(man_alpha)] <- 0
  expect_equivalent(man_beta + man_alpha, predict(example_vglmer, newdata = mixed_data, allow_missing_levels = T))
})

test_that("Prediction Matches for Missing in new.data", {
  N <- 50
  G <- 10
  x <- rnorm(N + G)
  g <- c(sample(1:G, N, replace = T), 1:G)
  alpha <- rnorm(G)

  y <- rbinom(n = N + G, size = 1, prob = plogis(-1 + x + alpha[g]))

  example_vglmer <- vglmer(
    formula = y ~ x + (1 | g), data = NULL,
    control = vglmer_control(iterations = 2), family = "binomial"
  )

  mixed_data <- data.frame(x = rnorm(20), g = rep(1:10, 2))
  rownames(mixed_data) <- letters[1:20]

  mixed_data$x[8] <- NA
  mixed_data$g[7] <- NA
  mixed_data$x[2] <- mixed_data$g[2] <- NA

  man_beta <- as.vector(cbind(1, mixed_data$x) %*% example_vglmer$beta$mean)
  man_alpha <- example_vglmer$alpha$mean[match(paste0("g @ (Intercept) @ ", mixed_data$g), rownames(example_vglmer$alpha$mean))]
  expect_equivalent(
    man_beta + man_alpha,
    predict(example_vglmer, newdata = mixed_data)
  )
})

test_that("Prediction Matches for Simulation", {
  
  if (env_test == "local"){
    N <- 1000
    G <- 10
  }else{
    N <- 50
    G <- 3
  }
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g]))

  est_glmer <- suppressMessages(suppressWarnings(lme4::glmer(y ~ x + (1 + x | g) + (1 | g2), family = binomial)))

  example_vglmer <- vglmer(
    formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = "binomial",
    control = vglmer_control(factorization_method = "weak", iterations = NITER)
  )

  mixed_data <- data.frame(x = rnorm(20), g = rep(1:10, 2), g2 = sample(1:25, 20, replace = T))
  rownames(mixed_data) <- letters[1:20]

  mixed_data$x[8] <- NA
  mixed_data$g[7] <- NA
  mixed_data$x[2] <- mixed_data$g[2] <- NA

  test_data <- rbind(mixed_data, data.frame(x = x, g = g, g2 = g2)[sample(1:length(y), 100, replace = T), ])

  point_predict <- predict(example_vglmer, newdata = test_data, allow_missing_levels = T)
  terms_predict <- predict(example_vglmer, newdata = test_data, allow_missing_levels = TRUE, type = 'terms')
  
  manual_fe <- test_data$x * coef(example_vglmer)['x'] + coef(example_vglmer)[1]
  expect_equal(
    terms_predict[,'FE'], 
    ifelse(is.na(point_predict), NA, manual_fe)
  )
  
  manual_g <- rowSums(ranef(example_vglmer)$g[test_data$g,2:3] * cbind(1, test_data$x))
  manual_g[is.na(manual_g) & !is.na(test_data$g)] <- 0
  expect_equal(
    terms_predict[,'g'], 
    ifelse(is.na(point_predict), NA, manual_g)
  )

  manual_g2 <- ranef(example_vglmer)$g2[test_data$g2,2]
  manual_g2[is.na(manual_g2) & !is.na(test_data$g2)] <- 0
  expect_equal(
    terms_predict[,'g2'], 
    ifelse(is.na(point_predict), NA, manual_g2)
  )
  
  expect_equivalent(
    manual_fe + manual_g + manual_g2, point_predict
  )
  
  if (env_test == "local"){
    n_samples <- 2 * 10^4
  }else{
    n_samples <- 2
  }
  mean_predict <- predict(example_vglmer,
    newdata = test_data,
    samples = n_samples, allow_missing_levels = T
  )
  # should have "clean" rownames
  expect_equal(rownames(mean_predict), as.character(1:nrow(mean_predict)))
  
  if (env_test == "local"){
    # Should be very close
    expect_equal(mean_predict$mean, point_predict, 0.01)
  }

  if (env_test == "local"){
    matrix_predict <- predict(example_vglmer,
                              newdata = test_data,
                              samples = 2 * 10^4, allow_missing_levels = T, samples_only = TRUE
    )
    matrix_predict <- colMeans(matrix_predict)
    expect_equivalent(
      c(coef(example_vglmer), as.vector(example_vglmer$alpha$mean)),
      matrix_predict, 0.01
    )
  }
})


test_that("Prediction Matches for vglmer after MAVB", {
  
  skip_on_cran()
  
  N <- 50
  G <- 4
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g]))

  example_vglmer <- vglmer(
    formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = "binomial",
    control = vglmer_control(factorization_method = "weak")
  )


  mixed_data <- data.frame(x = rnorm(20), g = rep(1:10, 2), g2 = sample(1:25, 20, replace = T))
  rownames(mixed_data) <- letters[1:20]

  mixed_data$x[8] <- NA
  mixed_data$g[7] <- NA
  mixed_data$x[2] <- mixed_data$g[2] <- NA

  test_data <- rbind(mixed_data, data.frame(x = x, g = g, g2 = g2)[sample(1:length(y), 100, replace = T), ])

  pred.MAVB <- predict_MAVB(example_vglmer,
    newdata = test_data,
    samples = 1000, allow_missing_levels = T
  )
  base_predict <- predict(example_vglmer, newdata = test_data, allow_missing_levels = T)

  expect_equal(pred.MAVB$mean, base_predict, tolerance = 0.1)
})

test_that("Prediction with Samples", {
  
  skip_on_cran()
  
  if (env_test == "local"){
    n_samples <- 2 * 10^4
  }else{
    n_samples <- 5
  }
  N <- 50
  G <- 10
  x <- rnorm(N + G)
  g <- c(sample(1:G, N, replace = T), 1:G)
  alpha <- rnorm(G)

  y <- rbinom(n = N + G, size = 1, prob = plogis(-1 + x + alpha[g]))

  example_vglmer <- vglmer(
    formula = y ~ x + (1 | g), data = NULL,
    control = vglmer_control(iterations = 2, factorization_method = 'strong'), 
    family = "binomial"
  )

  pred_samples <- predict(example_vglmer, newdata = data.frame(x = x, g = g), samples = 10, summary = F)
  expect_equivalent(dim(pred_samples), c(10, N + G))

  draw_coef <- predict(example_vglmer,
    newdata = data.frame(x = x, g = g),
    samples = n_samples, samples_only = T
  )
  expect_equivalent(dim(draw_coef), c(n_samples, G + 2))

  if (env_test == "local"){
    expect_equivalent(
      colMeans(draw_coef), format_vglmer(example_vglmer)$mean,
      tolerance = 0.02
    )
  }
  #Confirms that it works with "weak"
  example_vglmer <- vglmer(
    formula = y ~ x + (1 | g), data = NULL,
    control = vglmer_control(iterations = 2, factorization_method = 'weak'), 
    family = "binomial"
  )
  
  pred_samples <- predict(example_vglmer, newdata = data.frame(x = x, g = g), samples = 10, summary = F)
  expect_equivalent(dim(pred_samples), c(10, N + G))
  
  draw_coef <- predict(example_vglmer,
                       newdata = data.frame(x = x, g = g),
                       samples = n_samples, samples_only = T
  )
  expect_equivalent(dim(draw_coef), c(n_samples, G + 2))
  
  if (env_test == "local"){
    expect_equivalent(
      colMeans(draw_coef), format_vglmer(example_vglmer)$mean,
      tolerance = 0.02
    )
  }
})

test_that("Prediction with factor/categorical", {
  
  N <- 50
  G <- 10
  x <- rnorm(N + G)
  g <- c(sample(1:G, N, replace = T), 1:G)
  alpha <- rnorm(G)
  
  y <- rbinom(n = N + G, size = 1, prob = plogis(-1 + x + alpha[g]))
  l <- sample(letters[1:3], length(x), replace = T)
  
  example_vglmer <- vglmer(
    formula = y ~ x + l + (1  | g), data = NULL,
    control = vglmer_control(iterations = 2, factorization_method = 'strong'), 
    family = "binomial"
  )
  
  pred1a <- predict(example_vglmer, newdata = data.frame(x = 3, g = 2, l = 'a', stringsAsFactors = TRUE))
  expect_equivalent(pred1a, sum(coef(example_vglmer)  * c(1, 3, 0, 0)) + ranef(example_vglmer)$g[2,2])
  
  pred1b <- predict(example_vglmer, newdata = data.frame(x = 3, g = 2, l = 'a', stringsAsFactors = FALSE))
  expect_equivalent(pred1b, sum(coef(example_vglmer)  * c(1, 3, 0, 0)) + ranef(example_vglmer)$g[2,2])
  
  pred2a <- predict(example_vglmer, newdata = data.frame(x = -1, g = 100, l = 'c', stringsAsFactors = TRUE), allow_missing_levels = TRUE)
  expect_equivalent(pred2a, sum(coef(example_vglmer)  * c(1, -1, 0, 1)) )
  pred2b <- predict(example_vglmer, newdata = data.frame(x = -1, g = 100, l = 'c', stringsAsFactors = FALSE), allow_missing_levels = TRUE)
  expect_equivalent(pred2b, sum(coef(example_vglmer)  * c(1, -1, 0, 1)) )
  
  expect_error(predict(example_vglmer, 
    newdata = data.frame(x = -1, g = 100, l = 'd'), 
    allow_missing_levels = TRUE), regexp = 'has new level d')
})
