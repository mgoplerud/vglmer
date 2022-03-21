context("Test vglmer robustness to certain situations")

test_that("vglmer can run with objects in environment", {
  N <- 1000
  G <- 100
  G_names <- paste(sample(letters, G, replace = T), 1:G)
  x <- rnorm(N)
  g <- sample(G_names, N, replace = T)
  alpha <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[match(g, G_names)]))

  test_nodata <- tryCatch(suppressMessages(vglmer(y ~ x + (1 | g),
    data = NULL,
    control = vglmer_control(
      init = "zero",
      iterations = 1, print_prog = 10
    ),
    family = "binomial"
  )),
  error = function(e) {
    NULL
  }
  )
  expect_false(is.null(test_nodata))

  dta <- data.frame(Y = y, X = x, G = g)
  # Inject missingness into
  dta$Y[38] <- NA
  dta$X[39] <- NA
  dta$G[190] <- NA
  dta[108, ] <- NA
  test_missing <- tryCatch(suppressMessages(vglmer(Y ~ X + (1 | G),
    data = dta,
    control = vglmer_control(
      init = "zero", return_data = T,
      iterations = 1, print_prog = 10
    ),
    family = "binomial"
  )),
  error = function(e) {
    NULL
  }
  )
  # Confirm runs
  expect_false(is.null(test_missing))
  # Confirms deletion "works"
  expect_equivalent(dta$X[-c(38, 39, 108, 190)], test_missing$data$X[, 2])
  expect_equivalent(dta$Y[-c(38, 39, 108, 190)], test_missing$data$y)
})

test_that('vglmer runs with timing and "quiet=F"', {
  N <- 25
  G <- 2
  G_names <- paste(sample(letters, G, replace = T), 1:G)
  x <- rnorm(N)
  g <- sample(G_names, N, replace = T)
  alpha <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[match(g, G_names)]))

  est_simple <- suppressMessages(vglmer(y ~ x + (1 | g),
    data = NULL,
    control = vglmer_control(do_timing = T, quiet = F),
    family = "binomial"
  ))
  expect_true(inherits(est_simple$timing, "data.frame"))
  expect_gte(min(diff(est_simple$ELBO_trajectory$ELBO)), 0)
})

test_that('vglmer parses environment correctly', {
  rm(list=ls())  
  N <- 25
  G <- 2
  G_names <- paste(sample(letters, G, replace = T), 1:G)
  
  dta <- data.frame(x = rnorm(N), g = sample(G_names, N, replace = T))
  alpha <- rnorm(G)
  
  dta$y <- rbinom(n = N, size = 1, prob = plogis(-1 + dta$x + alpha[match(dta$g, G_names)]))
  dta$size <- rpois(n = N, lambda = 2) + 1
  dta$y_b <- rbinom(n = N, size = dta$size, prob = plogis(-1 + dta$x + alpha[match(dta$g, G_names)]))
  #runs with clean environment
  est_simple <- suppressMessages(vglmer(y ~ x + (1 | g), data = dta, family = 'binomial'))
  expect_true(inherits(est_simple, 'vglmer'))
  
  est_simple <- suppressMessages(vglmer(cbind(y_b, size) ~ x + (1 | g), data = dta, family = 'binomial'))
  expect_true(inherits(est_simple, 'vglmer'))
  

})

test_that("vglmer can run with 'debug' settings", {
  N <- 1000
  G <- 100
  G_names <- paste(sample(letters, G, replace = T), 1:G)
  x <- rnorm(N)
  g <- sample(G_names, N, replace = T)
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[match(g, G_names)]))
  
  # Debug to collect parameters
  est_vglmer <- vglmer(y ~ x + (1 | g), data = data.frame(y = y, x = x, g = g),
         family = 'binomial',
         control = vglmer_control(debug_param = TRUE))  
  
  expect_true(all(c('beta', 'alpha') %in% names(est_vglmer$parameter_trajectory)))

  est_vglmer <- vglmer(y ~ x + (1 | g), 
      data = data.frame(y = y, x = x, g = g),
      family = 'binomial',
      control = vglmer_control(debug_ELBO = TRUE))
  expect_true(!is.null(est_vglmer$debug_ELBO))
  
})

test_that("vglmer can run with exactly balanced classes", {
  N <- 1000
  G <- 5
  G_names <- paste(sample(letters, G, replace = T), 1:G)
  x <- rnorm(N)
  g <- sample(G_names, N, replace = T)
  alpha <- rnorm(G)
  
  y <- c(rep(0, N/2), rep(1, N/2))
  
  # Debug to collect parameters
  est_vglmer <- vglmer(y ~ x + (1 | g), data = data.frame(y = y, x = x, g = g),
                       family = 'binomial',
                       control = vglmer_control(do_SQUAREM = TRUE,
                                                factorization_method = 'strong',
                                                parameter_expansion = 'translation'))  
  
  expect_s3_class(est_vglmer, 'vglmer')
})
