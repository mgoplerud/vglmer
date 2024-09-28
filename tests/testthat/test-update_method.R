context("Test various update methods")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(191)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that("Joint vs Cyclical Update", {
  
  skip_on_cran()
  
  N <- 1000
  G <- 20
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:10, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(10)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g2]))


  for (v in c("weak", "partial", "strong")) {
    
    ex_vglmer_cyclic <- vglmer(
      formula = y ~ x + (1 | g) + (1 | g2), family = "binomial",
      data = NULL, control = vglmer_control(factorization_method = v, linpred_method = "cyclical", init = "zero")
    )

    ex_vglmer_joint <- vglmer(
      formula = y ~ x + (1 | g) + (1 | g2), family = "binomial",
      data = NULL, control = vglmer_control(factorization_method = v, linpred_method = "joint", init = "zero")
    )

    fmt_vglmer_cyclic <- format_vglmer(ex_vglmer_cyclic)
    fmt_vglmer_joint <- format_vglmer(ex_vglmer_joint)

    expect_equivalent(fmt_vglmer_cyclic, fmt_vglmer_joint, tolerance = 1e-4, scale = 1)

    if (v == "strong") {

      ex_vglmer_normal <- vglmer(
        formula = y ~ x + (1 | g) + (1 | g2), family = "binomial",
        data = NULL, 
        control = vglmer_control(factorization_method = v, 
           do_SQUAREM = FALSE,
          linpred_method = "solve_normal", init = "zero")
      )

      fmt_vglmer_normal <- format_vglmer(ex_vglmer_normal)
      expect_equivalent(fmt_vglmer_normal, fmt_vglmer_joint, tolerance = 1e-4, scale = 1)
    }
  }
})


test_that("Joint vs Cyclical Update (Nested)", {
  
  skip_on_cran()

  N <- 1000
  G <- 20
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- floor(g/3) + 1
  alpha <- rnorm(G)
  alpha2 <- rnorm(max(g2))
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g] + alpha2[g2]))
  
  
  fmla <- y ~ x + (1 | g2) + (1 | g)
  fmla_perm <- y ~ x + (1 | g) + (1 | g2)
  
  warning('This should work with *default* initialization...')
  
  ex_vglmer_cyclic_perm <- suppressMessages(vglmer(
    formula = fmla_perm, family = "binomial",
    data = NULL, control = vglmer_control(
      iterations = 100, print_prog = 500,
      init = 'EM', linpred_method = "cyclical")
  ))
  
  ex_vglmer_joint_perm <- vglmer(
    formula = fmla_perm, family = "binomial",
    control = vglmer_control(iterations = 100, print_prog = 500),
    data = NULL
  )
  
  ex_vglmer_cyclic <- suppressMessages(vglmer(
    formula = fmla, family = "binomial",
    data = NULL, control = vglmer_control(
      iterations = 100, print_prog = 500,
      init = 'EM', linpred_method = "cyclical")
  ))
  
  expect_equal(names(ranef(ex_vglmer_cyclic)), c('g2', 'g'))  
  expect_equal(names(ranef(ex_vglmer_cyclic_perm)), c('g', 'g2'))  
  
  ex_vglmer_joint <- vglmer(
    formula = fmla, family = "binomial",
    control = vglmer_control(iterations = 100,
                             print_prog = 500),
    data = NULL 
  )
  
  fmt_vglmer_cyclic_perm <- format_vglmer(ex_vglmer_cyclic_perm)
  fmt_vglmer_cyclic <- format_vglmer(ex_vglmer_cyclic)
  fmt_vglmer_joint <- format_vglmer(ex_vglmer_joint)
  
  expect_equivalent(fmt_vglmer_cyclic, fmt_vglmer_joint, tolerance = 1e-3, scale = 1)
  expect_equivalent(fmt_vglmer_cyclic, 
                    fmt_vglmer_cyclic_perm[
                      match(fmt_vglmer_cyclic$name, fmt_vglmer_cyclic_perm$name),
                    ], tolerance = 1e-2, scale = 1)
  
  expect_equivalent(ex_vglmer_joint$ELBO, ex_vglmer_joint_perm$ELBO,
                    tolerance = 1e-5, scale = 1)
  ex_vglmer_cyclic$ELBO$it <- NULL
  ex_vglmer_cyclic_perm$ELBO$it <- NULL
  expect_equivalent(ex_vglmer_cyclic$ELBO[,1], ex_vglmer_cyclic_perm$ELBO[,1],
                    tolerance = 1e-2, scale = 1)
  
})

test_that("Compare PX vs Non-PX", {
  
  skip_on_cran()
  
  N <- 1000
  G <- 20
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)

  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))


  ex_vglmer_px_t <- vglmer(
    formula = y ~ x + (1 | g), family = "binomial",
    data = NULL, control = vglmer_control(
      factorization_method = "strong", debug_px = TRUE,
      parameter_expansion = "translation", debug_ELBO = T, init = "zero"
    )
  )
  
  ex_vglmer_px <- vglmer(
    formula = y ~ x + (1 | g), family = "binomial",
    data = NULL, control = vglmer_control(
      factorization_method = "strong", debug_px = TRUE,
      parameter_expansion = "mean", debug_ELBO = T, init = "zero"
    )
  )

  ex_vglmer_no <- vglmer(
    formula = y ~ x + (1 | g), family = "binomial",
    data = NULL, control = vglmer_control(
      factorization_method = "strong",
      parameter_expansion = "none", debug_ELBO = T, init = "zero"
    )
  )

  expect_gte(min(diff(ex_vglmer_no$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gte(min(diff(ex_vglmer_px$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gte(min(diff(ex_vglmer_px_t$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  
  fmt_vglmer_px <- format_vglmer(ex_vglmer_px)
  fmt_vglmer_no <- format_vglmer(ex_vglmer_no)
  fmt_vglmer_px_t <- format_vglmer(ex_vglmer_px_t)
  
  expect_equivalent(fmt_vglmer_px, fmt_vglmer_no, tolerance = 1e-4, scale = 1)
  expect_equivalent(fmt_vglmer_px, fmt_vglmer_px_t, tolerance = 1e-4, scale = 1)
})

test_that("Compare VI r methods", {
  
  skip_on_cran()

  N <- 1000
  G <- 50
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)

  y <- rnbinom(n = N, mu = exp(-1 + x + alpha[g]), size = 5)
  data <- data.frame(y = y, x = x, g = g)

  list_output <- list()
  list_r <- list()
  
  for (v in c("VEM", "fixed")) {
    
    if (v == 'fixed'){
      v <- 1
    }
    example_vglmer <- suppressWarnings(vglmer(
      formula = y ~ x + (1 | g), data = data,
      family = "negbin",
      control = vglmer_control(parameter_expansion = 'mean',
        factorization_method = "strong",
        iterations = ifelse(v != 'VEM', 2, 1000),
        vi_r_method = v, init = "random"
      )
    ))
    # Test whether it monotonically increases
    if (v == "VEM") {
      expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
    }
    
    fmt_vglmer <- format_vglmer(example_vglmer)
    names(fmt_vglmer)[-1] <- paste0(v, "_", names(fmt_vglmer)[-1])
  }

  # Disable negative binomial tests for now
  # list_output <- Reduce(function(a, b) {
  #   merge(a, b, by = "name")
  # }, list_output)
  # expect_gte(min(as.vector(cor(list_output[, c("Laplace_mean", "delta_mean", "VEM_mean")]))), 0.95)
  # expect_gte(min(as.vector(cor(list_output[, c("Laplace_var", "delta_var", "VEM_var")]))), 0.95)
  # 
  # all_r <- sapply(list_r, FUN = function(i) {
  #   i$mu
  # })
  # # Check that mu are quite close
  # expect_lte(diff(range(all_r)), 0.02)
  # # Check that the mu standard errors are close for Laplace/delta
  # expect_lte(diff(sqrt(sapply(list_r, FUN = function(i) {
  #   i$sigma
  # }))[-1]), 0.02)
})
