context("test spline functionality")

test_that("check mgcv", {
  
  test_one <- mgcv:::interpret.gam(y ~ v_s(x2) + v_s(x))
  test_two <- vglmer_interpret.gam0(y ~ v_s(x2) + v_s(x))
  expect_equal(test_one, test_two)  
  
})

test_that("fit with non-default options", {
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
                    g = sample(state.abb[1:5], 100, replace = T),
                    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  
  # Custom knots as argument
  fit_knots <- vglmer(y ~ v_s(x, knots = quantile(dat$x, c(0.25, 0.75, 0.6))),
                      data = dat, family = 'binomial')
  
  fit_tpf <- vglmer(y ~ v_s(x, type = 'tpf'),
                      data = dat, family = 'binomial')
  fit_o <- vglmer(y ~ v_s(x, type = 'o'),
                    data = dat, family = 'binomial')
  
  
  expect_identical(fit_knots$spline$attr[[1]]$knots, quantile(dat$x, c(0.25, 0.6, 0.75)))
  
  fit_knotsa <- vglmer(y ~ v_s(x, knots = quantile(dat$x, c(0.25, 0.6, 0.75))),
                      data = dat, family = 'binomial')
  
  expect_equal(ELBO(fit_knots), ELBO(fit_knotsa))
  expect_equal(ranef(fit_knots), ranef(fit_knotsa))
  
  expect_gt(min(diff(ELBO(fit_tpf, 'trajectory'))), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(ELBO(fit_o, 'trajectory'))), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(ELBO(fit_knots, 'trajectory'))), -sqrt(.Machine$double.eps))
  
})

test_that("fit splines with missing data", {
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
                    g = sample(state.abb[1:5], 100, replace = T),
                    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  dat$y[sample(100, 5)] <- NA
  dat$g[sample(100, 5)] <- NA  
  dat$x[sample(100, 5)] <- NA
  
  fit_spline <- vglmer(y ~ v_s(x, by = f) + (1 | g), 
               data = dat, family = 'binomial')
  expect_s3_class(fit_spline, 'vglmer')
  expect_gt(min(diff(fit_spline$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))

})

test_that("Check failures of spline fitting", {
  
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
  
  x <- splines::bs(x = rnorm(500))[,]
  by_values <- sample(letters, 500, replace = T)
  u_by <- sort(unique(by_values))
  x_by <- sparseMatrix(i = 1:length(by_values), j = match(by_values, u_by), x = 1)
  
  test_m <- t(KhatriRao(t(x_by), t(x)))
  manual_m <- do.call('cbind', lapply(u_by, FUN=function(u){
    drop0(Diagonal(x = (by_values == u)) %*% x)
  }))
  colnames(manual_m) <- NULL
  expect_equal(test_m, manual_m)  
  
})

test_that("Basic spline tests (run and predict)", {
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
                    g = sample(state.abb[1:5], 100, replace = T),
                    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  

  # Check runs with 1
  
  m1 <- vglmer(y ~ x + x2 + v_s(x), 
               data = dat, family = 'binomial')
  expect_equal(length(coef(m1)), 3)
  m1a <- vglmer(y ~ v_s(x) + x2, data = dat, 
                family = 'binomial')
  
  # Check same convergence
  
  expect_equal(ranef(m1), ranef(m1a), tol = 1e-5)
  expect_equal(coef(m1), 
    coef(m1a)[match(names(coef(m1)), names(coef(m1a)))],
    tol = 1e-5
  )
  expect_equal(ELBO(m1), ELBO(m1a))
  expect_equal(ELBO(m1, 'trajectory'), ELBO(m1a, 'trajectory'), tol = 1e-5)
  
  # Check runs with 2
  m2 <- vglmer(y ~ v_s(x2) + v_s(x), 
               control = vglmer_control(iterations = 20, 
                                        factorization_method = 'intermediate'),
               data = dat, family = 'binomial')
  # Check runs with "by"
  m3 <- vglmer(y ~ v_s(x2) + v_s(x, by = f), data = dat, 
               family = 'binomial',
               control = vglmer_control(iterations = 20))
  # Check runs with RE 
  m4 <- vglmer(y ~ v_s(x, by = f) + (1 | g), 
               data = dat, family = 'binomial',
               control = vglmer_control(iterations = 20))

  # Check runs with double "by"
  m5 <- vglmer(y ~ v_s(x, by = f) + v_s(x, by = g), 
               data = dat, family = 'binomial',
               control = vglmer_control(iterations = 20, 
                                        factorization_method = 'strong'))
  # Check runs with 2 "by" for single grouping
  m6 <- vglmer(y ~ v_s(x, by = f) + v_s(x2, by = f), 
               data = dat, family = 'binomial',
               control = vglmer_control(iterations = 20, 
                                        factorization_method = 'strong',
                                        linpred_method = 'cyclical'))
  expect_equal(ncol(m6$sigma$cov[[1]]), 3)
  
  # Check ELBO increases
  expect_gt(min(diff(m1$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m2$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m3$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m4$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m5$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m6$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  
  # Prediction tests
  # check decomposition of D and then retransformation
  # test with custom knots and prediction outside of knots/etc.
  
})


