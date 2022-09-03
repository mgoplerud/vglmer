context("test spline functionality")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(41414)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

print(paste0(NITER, ' for tests because ', env_test))

test_that("fit with non-default options", {
  
  skip_on_cran()
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
    g = sample(state.abb[1:5], 100, replace = T),
    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  
  # Custom knots as argument
  fit_knots <- vglmer(y ~ v_s(x, knots = quantile(dat$x, c(0.25, 0.75, 0.6))),
                      data = dat, family = 'binomial',
                      control = vglmer_control(iterations = 15))
  
  fit_tpf <- vglmer(y ~ v_s(x, type = 'tpf'),
    data = dat, family = 'binomial', control = vglmer_control(iterations = NITER))
  fit_o <- vglmer(y ~ v_s(x, type = 'o'),
    data = dat, family = 'binomial',
    control = vglmer_control(iterations = NITER))
  
  
  expect_identical(fit_knots$spline$attr[[1]]$knots, quantile(dat$x, c(0.25, 0.6, 0.75)))

  expect_gt(min(diff(ELBO(fit_tpf, 'trajectory'))), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(ELBO(fit_o, 'trajectory'))), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(ELBO(fit_knots, 'trajectory'))), -sqrt(.Machine$double.eps))
  
})

test_that("fit splines with missing data", {
  
  skip_on_cran()
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
                    g = sample(state.abb[1:5], 100, replace = T),
                    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  dat$y[sample(100, 5)] <- NA
  dat$g[sample(100, 5)] <- NA  
  dat$x[sample(100, 5)] <- NA
  dat$f[sample(100, 5)] <- NA
  
  fit_spline <- vglmer(y ~ v_s(x, by = f) + (1 | g), 
               data = dat, family = 'binomial', control = vglmer_control(iterations = 5))
  expect_s3_class(fit_spline, 'vglmer')
  expect_gt(min(diff(fit_spline$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))

})

test_that("Check failures of spline fitting", {
  
  suppressWarnings(rm(f))
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
                    g = sample(state.abb[1:5], 100, replace = T),
                    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))

  nox <- dat[, !(colnames(dat) %in% 'x')]
  x <- rnorm(nrow(dat))
  nof <- dat[, !(colnames(dat) %in% 'f')]
  
  expect_error(vglmer(y ~ v_s(x, by = f) + (1 | g), 
                       data = nox, 
                      control = vglmer_control(verify_columns = TRUE),
                      family = 'binomial'))
  
  expect_error(vglmer(y ~ v_s(x, by = f) + (1 | g), 
                      data = nof, family = 'binomial'))
  
})


test_that("test spline 'by' construction", {

  x <- splines::bs(x = rnorm(10))[,]
  by_values <- sample(letters, 10, replace = T)
  u_by <- sort(unique(by_values))
  x_by <- sparseMatrix(i = 1:length(by_values), j = match(by_values, u_by), x = 1)
  
  test_m <- t(Matrix::KhatriRao(t(x_by), t(x)))
  manual_m <- do.call('cbind', lapply(u_by, FUN=function(u){
    drop0(Diagonal(x = (by_values == u)) %*% x)
  }))
  colnames(manual_m) <- NULL
  expect_equal(test_m, manual_m)  
  
})

test_that("CRAN basic spline tests", {
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
                    g = sample(state.abb[1:5], 100, replace = T),
                    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  
  
  m1 <- vglmer(y ~ x + x2 + v_s(x), 
               control = vglmer_control(iteration = 5),
               data = dat, family = 'binomial')
  expect_gt(min(diff(ELBO(m1, 'traj'))), - sqrt(.Machine$double.eps))  
  
  
  # Check runs with 2
  m2 <- vglmer(y ~ v_s(x2) + v_s(x), 
               control = vglmer_control(
                 iterations = NITER, print_prog = 20, prior_variance = 'mean_exists'),
               data = dat, family = 'binomial')
  # Check runs with "by"
  m3 <- vglmer(y ~ v_s(x2) + v_s(x, by = f), data = dat, 
               family = 'binomial',
               control = vglmer_control(iterations = NITER, print_prog = 20))
  # Check runs with RE 
  m4 <- vglmer(y ~ v_s(x, by = f) + (1 | g), 
               data = dat, family = 'binomial',
               control = vglmer_control(iterations = NITER, print_prog = 20))
  
  # Check runs with double "by"
  m5 <- vglmer(y ~ v_s(x, by = f) + v_s(x, by = g), 
               data = dat, family = 'binomial',
               control = vglmer_control(iterations = NITER, print_prog = 20))
  # Check runs with 2 "by" for single grouping
  m6 <- vglmer(y ~ v_s(x, by = f) + v_s(x2, by = f), 
     data = dat, family = 'binomial',
     control = vglmer_control(iterations = NITER, print_prog = 20,
      factorization_method = 'strong'))
  expect_equal(ncol(m6$sigma$cov[[1]]), 3)

  m7 <- vglmer(y ~ v_s(x, by = f) + v_s(x2, by = f) , 
     data = dat, family = 'binomial',
     control = vglmer_control(iterations = NITER, print_prog = 20,
      factorization_method = 'strong'))
  
  # Check ELBO increases
  expect_gt(min(diff(m1$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m2$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m3$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m4$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m5$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m6$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m7$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  
  expect_vector(predict(m1, newdata = dat[1:5,]))
  expect_vector(predict(m2, newdata = dat[1:5,]))
  expect_vector(predict(m3, newdata = dat[1:5,]))
  expect_vector(predict(m4, newdata = dat[1:5,]))
  expect_vector(predict(m5, newdata = dat[1:5,]))
  expect_vector(predict(m6, newdata = dat[1:5,]))
  expect_vector(predict(m7, newdata = dat[1:5,]))
  
})

test_that("Test order of splines", {
  
  skip_on_cran()
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
                    g = sample(state.abb[1:5], 100, replace = T),
                    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  
  m1 <- vglmer(y ~ x + x2 + v_s(x), 
     data = dat, family = 'binomial',
     control = vglmer_control(iterations = NITER))
  expect_equal(length(coef(m1)), 3)
  m1a <- vglmer(y ~ x2 + x + v_s(x), data = dat, 
    family = 'binomial', control = vglmer_control(iterations = NITER))
  
  expect_equal(ELBO(m1a), ELBO(m1))
  expect_equal(ranef(m1), ranef(m1a), tol = 1e-4, scale = 1)
  expect_equal(coef(m1), 
    coef(m1a)[match(names(coef(m1)), names(coef(m1a)))],
    tol = 1e-4, scale = 1
  )
  
})

test_that("Prediction spline test", {
  
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
                    g = sample(state.abb[1:5], 100, replace = T),
                    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  
  m1 <- vglmer(y ~ x + x2 + v_s(x), 
               data = dat, family = 'linear',
               control = vglmer_control(iterations = 3))
  # Check that prediction works for simple spline case
  expect_vector(predict(m1, newdata = dat))
  
  
})

# TO-DO tests
# check decomposition of D and then re-transformation
# test with custom knots and prediction outside of knots/etc.


