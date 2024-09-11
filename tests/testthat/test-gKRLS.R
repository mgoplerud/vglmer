context("Test Generic Methods")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(151)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that("Test gKRLS runs in a few configurations", {
  
  N <- 1000
  G <- 10
  x <- rnorm(N)
  x2 <- rnorm(N)
  x3 <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- rnorm(G)
  alpha2 <- rnorm(G)
  alpha3 <- rnorm(G)
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + sin(alpha3[g] * x) + cos(x2 * x) + alpha[g] + alpha2[g]))
  
  if (requireNamespace("gKRLS", quietly = TRUE)) {
    require(gKRLS)
    
    example_vglmer <- vglmer(
      formula = y ~ x + (1 | g) + v_s(x, x2, type = 'gKRLS'), data = NULL, family = "binomial",
      control = vglmer_control(iterations = NITER, return_data = TRUE, factorization_method = "strong")
    )
    expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)
    expect_false('x2' %in% names(fixef(example_vglmer)))
    
    id_vglmer <- example_vglmer$internal_parameters$spline$attr[[1]]$knots$subsampling_id
    fit_mgcv <- gam(y ~ x + s(x, x2, bs = 'gKRLS', xt = gKRLS(sketch_method = id_vglmer)))
    # Check that the *design* aligns from gKRLS and vglmer
    expect_equivalent(
      as.matrix(example_vglmer$data$Z[,grepl(colnames(example_vglmer$data$Z), pattern='spline')]),
      model.matrix(fit_mgcv)[,-1:-2]
    )
    
    # Check with "weak" and many smooth terms but not "proper" REs
    example_vglmer <- vglmer(
      formula = y ~ v_s(x2) +
        v_s(x,type = 'gKRLS') + 
        v_s(x, x2, type = 'gKRLS',  xt = gKRLS(sketch_method = 'none', standardize = 'none')),
      data = NULL, family = "binomial",
      control = vglmer_control(iterations = NITER, factorization_method = "weak")
    )
    expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)
    expect_false('x' %in% names(fixef(example_vglmer)))
    expect_true('x2' %in% names(fixef(example_vglmer)))
    # Check that prediction works
    test_pred <- predict(example_vglmer, newdata = data.frame(
      x = rnorm(10), x2 = rnorm(10), g = '2'
    ))
    
    # Check that this works with "by"
    f_g <- factor(g)
    example_vglmer <- vglmer(
      formula = y ~ x + 
        v_s(x, x2, type = 'gKRLS',  by = f_g),
      data = NULL, family = "binomial",
      control = vglmer_control(iterations = NITER, factorization_method = "strong")
    )
    expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)
    expect_true('f_g' %in% names(ranef(example_vglmer)))
    expect_true(ncol(ranef(example_vglmer)$f_g) == 2)
    expect_false('x2' %in% names(fixef(example_vglmer)))
    
    test_pred <- predict(example_vglmer, allow_missing_levels = TRUE,
                         newdata = data.frame(x, x2, f_g = 1590)[1:15,])
    
    # Runs fine with "data" provided directly
    example_vglmer <- vglmer(
      formula = yn ~ a + v_s(a, b, type = 'gKRLS'),
      data = data.frame(yn = y, a = x, b = x2),
      family = 'linear'
    )
    
  }
  
})

test_that("Test random walk runs in a few configurations", {

  N <- 1000
  G <- 10
  x <- rnorm(N)
  x2 <- rnorm(N)
  x3 <- rnorm(N)
  g <- sample(1:G, N, replace = T)
  g2 <- sample(1:G, N, replace = T)
  alpha <- sort(rnorm(G))
  alpha2 <- rnorm(G)
  alpha3 <- rnorm(G)
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + sin(alpha3[g] * x) + cos(x2 * x) + alpha[g] + alpha2[g]))
  
  fg <- g
  example_vglmer <- expect_error(
    vglmer(
      formula = y ~ x + v_s(fg, type = 'randwalk'), data = NULL, family = "binomial",
      control = vglmer_control(iterations = NITER, return_data = TRUE, factorization_method = "strong")
    ), regexp='character or factor'
  )
  fg <- factor(g)
  example_vglmer <- expect_error(
    vglmer(
      formula = y ~ x + v_s(fg, type = 'randwalk'), data = NULL, family = "binomial",
      control = vglmer_control(iterations = NITER, return_data = TRUE, factorization_method = "strong")
    ), regexp='ordered'
  )
  fg <- factor(g, ordered = TRUE)
  example_vglmer <- vglmer(
      formula = y ~ x + v_s(fg, type = 'randwalk'), data = NULL, family = "binomial",
      control = vglmer_control(iterations = NITER, return_data = TRUE, factorization_method = "strong")
  )
  expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), 0)
  
  # Check that manual predictions work
  
  new_fg <- data.frame(fg = sample(1:5, size = 10, replace = T), x = rnorm(10))
  pred_new_fg <- predict(example_vglmer, newdata = new_fg)
  manual_pred <- sparseMatrix(i = 1:nrow(new_fg), j = match(new_fg$fg, 
      example_vglmer$internal_parameters$spline$attr[[1]]$knots$levels),
      x = 1, dims = c(nrow(new_fg), 10)
  )
  manual_pred <- as.vector(manual_pred %*% 
    example_vglmer$internal_parameters$spline$attr[[1]]$knots$transf_mat %*%
    example_vglmer$alpha$mean + 
    cbind(1, new_fg$x) %*% coef(example_vglmer)
  )
  expect_equivalent(manual_pred, pred_new_fg)
  
  fs <- sample(state.abb, N, replace = T)
  fs <- data.frame(
    y = y, state = factor(fs, ordered = TRUE)
  )
  fs$y <- fs$y + rnorm(nrow(fs)) + cos(2 * pi * ((0:50)/50))[match(fs$state, state.abb)]
  
  example_randwalk <- vglmer(
    formula = y ~ x + v_s(state, type = 'randwalk'), data = fs, family = "linear",
    control = vglmer_control(iterations = NITER, factorization_method = "weak")
  )
  expect_gte(min(diff(example_randwalk$ELBO_trajectory$ELBO)), -.Machine$double.eps)
  
  pred_randwalk <- predict(example_randwalk,
   newdata = data.frame(x = 0, state = state.abb),
   type = 'terms')
  
  # Check that prediction works with NA and that "by" works
  th <- 8
  example_vglmer <- vglmer(
    yo ~ v_s(a, type = 'randwalk', by = st),
    data = data.frame(yo = y, a = cut(x, th, ordered_result = TRUE), st = fg),
    family = 'binomial', control = vglmer_control(iterations = NITER)
  )
  expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), -.Machine$double.eps)
  # Check that REs are added for the "by" term
  expect_equivalent(
    names(ranef(example_vglmer)),
    c('st', 'spline-(a)-1-base', 'spline-(a)-1-int')
  )
  
  full_grid <- expand.grid(st = unique(fg), a = unique(cut(x, th, ordered_result = TRUE)))
  full_grid[31,]$a <- NA
  pred_grid <- predict(example_vglmer, newdata = full_grid)
  expect_true(is.na(pred_grid[31]))
  
  # Check that prediction fails
  example_vglmer <- vglmer(
    yo ~ v_s(a, type = 'randwalk') + (1 | fg),
    data = data.frame(yo = y, a = cut(x, th, ordered_result = TRUE), st = fg),
    family = 'binomial', control = vglmer_control(iterations = NITER)
  )
  expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), -.Machine$double.eps)
  expect_error(
    predict(example_vglmer, newdata = data.frame(fg = 3, a = '(1.00, 0.99]', stringsAsFactors = F)),
    regexp = 'cannot work if levels'
  )
  
})
