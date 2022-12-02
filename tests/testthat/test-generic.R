context("Test Generic Methods")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(141)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that("Test generic methods (fixed, ranef, coef, vcov)", {
  N <- 100
  G <- 10
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g]))

  example_vglmer <- vglmer(
    formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = "binomial",
    control = vglmer_control(factorization_method = "strong")
  )
  expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)

  expect_equivalent(coef(example_vglmer), as.vector(example_vglmer$beta$mean))
  expect_equivalent(fixef(example_vglmer), coef(example_vglmer))
  expect_equivalent(vcov(example_vglmer), as.matrix(example_vglmer$beta$var))

  generic_ranef <- ranef(example_vglmer)

  expect_equivalent(example_vglmer$alpha$mean[-(1:(2 * G))], generic_ranef$g2$`(Intercept)`)
  expect_equivalent(example_vglmer$alpha$mean[(1:(2 * G))], as.vector(t(generic_ranef$g[, -1])))
  expect_equivalent(example_vglmer$alpha$dia.var[[1]], example_vglmer$alpha$var[[1]][,c(1,4)])
  expect_equivalent(as.vector(t(as.matrix(attributes(generic_ranef$g)$variance[,-1]))), 
    diag(crossprod(example_vglmer$alpha$decomp_var[,1:(2*G)])),
    tol = 1e-6)
  expect_equivalent(as.vector(t(as.matrix(attributes(generic_ranef$g2)$variance[,-1]))), 
    diag(crossprod(example_vglmer$alpha$decomp_var[,-(1:(2*G))])),
    tol = 1e-6)
  
  test_fmt <- format_vglmer(example_vglmer)
  
  expect_equivalent(test_fmt[-1:-2,][1:(2 * G),]$mean, example_vglmer$alpha$mean[(1:(2 * G))])
  expect_equivalent(test_fmt[-1:-2,][-(1:(2 * G)),]$mean, example_vglmer$alpha$mean[-(1:(2 * G))])
  expect_equivalent(test_fmt[-1:-2,][1:(2 * G),]$var, as.vector(t(as.matrix(attributes(generic_ranef$g)$variance[,-1]))))
  expect_equivalent(test_fmt[-1:-2,][-(1:(2 * G)),]$var, as.vector(t(as.matrix(attributes(generic_ranef$g2)$variance[,-1]))))
  
})

test_that("Test that print and summary run", {
  N <- 100
  G <- 10
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g]))

  example_vglmer <- vglmer(
    formula = y ~ x + (1 + x | g) + (1 | g2), data = NULL, family = "binomial",
    control = vglmer_control(factorization_method = "strong")
  )
  expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)
  save_print <- invisible(capture.output(print(example_vglmer)))
  save_summary <- invisible(capture.output(summary(example_vglmer)))
})

