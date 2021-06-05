context("Match glmer")

test_that("Compare against lmer", {
  
  N <- 1000
  G <- 100
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  sigmasq <- abs(rcauchy(1))
  coef_scale <- rexp(1, rate = 1/2)
  y <- rnorm(n = N, mean = coef_scale * (-1 + x + alpha[g]) * sqrt(sigmasq), sd = sqrt(sigmasq))
  
  est_glmer <- suppressWarnings(lme4::lmer(y ~ x + (1 | g), REML = FALSE))
  fmt_glmer <- format_glmer(est_glmer)
  
  for (v in c("weak", "partial", "strong")) {
    example_vglmer <- vglmer(
      formula = y ~ x + (1 | g), data = NULL,
      control = vglmer_control(factorization_method = v, init = "zero"),
      family = "linear"
    )
    
    expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)
    
    fmt_vglmer <- format_vglmer(example_vglmer)
    comp_methods <- merge(fmt_glmer, fmt_vglmer, by = c("name"))
    
    cor_mean <- with(comp_methods, cor(mean.x, mean.y))
    expect_gt(cor_mean, expected = 0.99)
    
    expect_lt(with(comp_methods, mean(abs(mean.x - mean.y))), 0.1)
    
    expect_equal(sigma(example_vglmer), sigma(est_glmer), tol = 1e-1)
  }
})


test_that("Compare against glmer", {
  N <- 1000
  G <- 100
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))

  est_glmer <- suppressWarnings(lme4::glmer(y ~ x + (1 | g), family = binomial))
  fmt_glmer <- format_glmer(est_glmer)

  for (v in c("weak", "partial", "strong")) {
    example_vglmer <- vglmer(
      formula = y ~ x + (1 | g), data = NULL,
      control = vglmer_control(factorization_method = v, init = "zero"),
      family = "binomial"
    )

    expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)

    fmt_vglmer <- format_vglmer(example_vglmer)
    comp_methods <- merge(fmt_glmer, fmt_vglmer, by = c("name"))

    cor_mean <- with(comp_methods, cor(mean.x, mean.y))
    expect_gt(cor_mean, expected = 0.99)
  }
})


test_that("Compare against glmer.nb", {
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

  for (v in c("weak", "partial", "strong")) {
    example_vglmer <- vglmer(
      formula = y ~ x + (1 | g), data = data,
      family = "negbin",
      control = vglmer_control(factorization_method = v)
    )
    # Test whether it monotonically increases
    expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)

    fmt_vglmer <- format_vglmer(example_vglmer)
    comp_methods <- merge(fmt_glmer, fmt_vglmer, by = c("name"))
    cor_mean <- with(comp_methods, cor(mean.x, mean.y))
    # Test whether it is close to the truth.
    expect_gt(cor_mean, expected = 0.99)
  }
})

test_that("EM_prelim matches glm", {
  N <- 100
  p <- 5

  Z <- matrix(rnorm(N * p), ncol = p)
  beta <- runif(p, -1, 1)

  y <- rbinom(N, 1, plogis(Z %*% beta))

  est_glm <- glm(y ~ Z, family = binomial)
  est_init <- EM_prelim_logit(
    X = drop0(matrix(1, nrow = N)),
    Z = drop0(Z), s = y - 1 / 2, pg_b = rep(1, N), iter = 200, ridge = Inf
  )
  est_init <- c(est_init$beta, est_init$alpha)
  expect_equal(as.vector(coef(est_glm)), est_init, tolerance = 1e-4)
})


test_that("EM_prelim matches glm.nb", {
  quine <- MASS::quine
  N <- nrow(quine)
  quine.nb1 <- MASS::glm.nb(Days ~ Sex / (Age + Eth * Lrn), data = quine)
  X <- model.matrix(quine.nb1)
  y <- quine$Days

  est_init <- EM_prelim_nb(
    X = drop0(matrix(1, nrow = N)), Z = drop0(X[, -1]), y = y,
    est_r = quine.nb1$theta, iter = 100, ridge = Inf
  )
  est_init <- c(est_init$beta, est_init$alpha)
  expect_equal(as.vector(coef(quine.nb1)), est_init, tolerance = 1e-4)
})

test_that("Compare against glmer (binomial)", {
  N <- 1000
  G <- 100
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)

  size <- rpois(N, 1) + 1
  y <- rbinom(n = N, size = size, prob = plogis(-1 + x + alpha[g]))

  est_glmer <- suppressWarnings(lme4::glmer(cbind(y, size - y) ~ x + (1 | g), family = binomial))
  fmt_glmer <- format_glmer(est_glmer)

  example_vglmer <- vglmer(
    formula = cbind(y, size - y) ~ x + (1 | g), data = NULL,
    family = "binomial"
  )

  expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)

  fmt_vglmer <- format_vglmer(example_vglmer)
  comp_methods <- merge(fmt_glmer, fmt_vglmer, by = c("name"))

  cor_mean <- with(comp_methods, cor(mean.x, mean.y))
  expect_gt(cor_mean, expected = 0.99)
})
