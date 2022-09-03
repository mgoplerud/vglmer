context("Test MAVB")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(161)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that("Short test of MAVB running for CRAN", {
  N <- 50
  G <- 5
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g]))

  example_vglmer <- vglmer(
    formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = "binomial",
    control = vglmer_control(factorization_method = "weak",
                             iterations = 5)
  )

  mavb_samples <- tryCatch(MAVB(object = example_vglmer, samples = 10),
    error = function(e) {
      NULL
    }
  )
  expect_false(is.null(mavb_samples))
  
  
  example_vglmer <- vglmer(
    formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = "binomial",
    control = vglmer_control(factorization_method = "strong", iterations = 5)
  )
  
  mavb_samples <- tryCatch(MAVB(object = example_vglmer, samples = 10),
                           error = function(e) {
                             NULL
                           }
  )
  expect_false(is.null(mavb_samples))
  
})