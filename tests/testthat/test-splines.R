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

env_test <- 'local'
NITER <- 2000
print(paste0(NITER, ' for tests because ', env_test))

test_that("fit with non-default options", {
  
  skip_on_cran()
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
    g = sample(state.abb[1:5], 100, replace = T),
    f = sample(letters[1:5], 100, replace = T))
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  
  # Custom knots as argument
  custom_knots <- quantile(dat$x, c(0, 0.25, 0.75, 0.6, 1))
  fit_knots <- vglmer(y ~ v_s(x, knots = custom_knots),
                      data = dat, family = 'binomial',
                      control = vglmer_control(iterations = 15))
  
  fit_tpf <- vglmer(y ~ v_s(x, type = 'tpf'),
    data = dat, family = 'binomial', control = vglmer_control(iterations = NITER))
  fit_o <- vglmer(y ~ v_s(x, type = 'o'),
    data = dat, family = 'binomial',
    control = vglmer_control(iterations = NITER))
  
  
  expect_identical(fit_knots$internal_parameters$spline$attr[[1]]$knots, quantile(dat$x, c(0, 0.25, 0.6, 0.75, 1)))

  expect_gt(min(diff(ELBO(fit_tpf, 'trajectory'))), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(ELBO(fit_o, 'trajectory'))), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(ELBO(fit_knots, 'trajectory'))), -sqrt(.Machine$double.eps))
  
  fit_tpf <- vglmer(y ~ v_s(x, knots = 2, type = 'tpf'),
      data = dat, family = 'binomial', 
      control = vglmer_control(iterations = NITER))
  expect_length(fit_tpf$internal_parameters$spline$attr[[1]]$knots, 2)
  expect_equivalent(fit_tpf$internal_parameters$spline$attr[[1]]$knots, 
    quantile(dat$x, seq(0,1,length.out=4)[-c(1,4)]))
})

test_that("fit splines works with one knot", {
  
  
  dat <- data.frame(x = rnorm(100))
  dat$y <- rbinom(100, 1, plogis(exp(dat$x)))

  for (loop_type in c('tpf', 'o')){
    m1 <- vglmer(y ~ v_s(x, type = loop_type, knots = 1), dat = dat, 
                 control = vglmer_control(iterations = NITER),
                 family = 'linear')
    pred_m1 <- predict(m1, newdata = data.frame(x = seq(-5, 5, length.out=100)))  
    expect_equal(length(m1$internal_parameters$spline$attr[[1]]$knots), 1)
    if (loop_type == "tpf"){
      # Expect only one "bend"
      expect_equal(sum( abs(diff(diff(pred_m1))) > sqrt(.Machine$double.eps)), 2)
    }
  }
  
  expect_error(vglmer(y ~ v_s(x, knots = 0.5), data = dat, family = 'linear'), regexp = 'If an integer')
  expect_warning(est_vglmer <- vglmer(y ~ v_s(x, knots = 0.5, force_vector = T), 
    data = dat, control = vglmer_control(iterations = 2), family = 'linear'),
    regexp = 'observed data')
  expect_equivalent(0.5, est_vglmer$internal_parameters$spline$attr[[1]]$knots)
  
  est_vglmer <- expect_warning(vglmer(y ~ v_s(x, knots = 1.2, force_vector = T),
         control = vglmer_control(iterations = 2),
         data = dat, family = 'linear'), regexp = 'observed data')
  expect_equivalent(est_vglmer$internal_parameters$spline$attr[[1]]$knots, 1.2)
  expect_vector(predict(est_vglmer, newdata = data.frame(x = seq(-3,3,length.out=5))), 5)
  expect_warning(vglmer(y ~ v_s(x, knots = 1.2),
         control = vglmer_control(iterations = 2),
         data = dat, family = 'linear'), regexp = 'as.integer')
  
  est_vglmer <- vglmer(y ~ v_s(x, knots = 3), control = vglmer_control(iterations = 2), 
         data = dat, family = 'linear')
  expect_length(est_vglmer$internal_parameters$spline$attr[[1]]$knots, 3)
  est_vglmer <- vglmer(y ~ v_s(x), control = vglmer_control(iterations = 2), 
                       data = dat, family = 'linear')
  expect_length(est_vglmer$internal_parameters$spline$attr[[1]]$knots, floor(length(unique(dat$x))/4))
})

test_that("fit and predict with splines and missing data", {
  
  dat <- data.frame(x = rnorm(100), x2 = rexp(100),
    g = sample(state.abb[1:5], 100, replace = T),
    f = sample(letters[1:5], 100, replace = T), 
    stringsAsFactors = F)
  
  dat$y <- rbinom(100, 1, plogis(dat$x * runif(5)[match(dat$f, letters)]))
  dat$y[sample(100, 5)] <- NA
  dat$g[sample(100, 5)] <- NA  
  dat$x[sample(100, 5)] <- NA
  dat$f[sample(100, 5)] <- NA
  
  fit_spline <- vglmer(y ~ v_s(x, by = f) + (1 | g), 
               data = dat, family = 'binomial', control = vglmer_control(iterations = NITER))
  expect_s3_class(fit_spline, 'vglmer')
  expect_gt(min(diff(fit_spline$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  
  dat$f[6] <- 'h'
  pred_spline <- predict(fit_spline, newdata = dat, allow_missing_levels = TRUE)
  pred_spline_terms <- predict(fit_spline,
    newdata = dat, type = 'terms', allow_missing_levels = TRUE)  
  
  expect_equivalent(pred_spline, rowSums(pred_spline_terms))
  
  raw_X <- vglmer_build_spline(x = dat$x, 
    knots = fit_spline$internal_parameters$spline$attr[[1]]$knots,
    type = fit_spline$internal_parameters$spline$attr[[1]]$type,
    Boundary.knots = fit_spline$internal_parameters$spline$attr[[1]]$Boundary.knots, 
    by = NULL)[[1]]$x

  expect_equivalent(
    pred_spline_terms[, 'spline-x-1-base'],
    ifelse(is.na(pred_spline), NA, as.vector(raw_X %*% ranef(fit_spline)[['spline-x-1-base']][,2]))
  )  

  reshape_spline <- matrix(ranef(fit_spline)[['spline-x-1-int']][,2], ncol = 5)
  pred_spline_inter <- rowSums(raw_X * t(reshape_spline)[match(dat$f, letters[1:5]),])
  pred_spline_inter[!is.na(dat$f) & is.na(pred_spline_inter)] <- 0
  expect_equivalent(
    pred_spline_terms[, 'spline-x-1-int'],
    ifelse(is.na(pred_spline), NA, pred_spline_inter)
  ) 
  
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
  
  dat$y <- rbinom(100, 1, plogis(sin(sqrt(dat$x2 * 4)) + cos(dat$x * 3) * runif(5)[match(dat$f, letters)]))
  
  
  m1 <- vglmer(y ~ x + x2 + v_s(x), 
               control = vglmer_control(iteration = NITER),
               data = dat, family = 'binomial')
  expect_gt(min(diff(ELBO(m1, 'traj'))), - sqrt(.Machine$double.eps))  
  
  
  # Check runs with 2
  m2 <- vglmer(y ~ v_s(x2) + v_s(x), 
               control = vglmer_control(
                 iterations = NITER, print_prog = 20),
               data = dat, family = 'binomial')
  # Check runs with "by"
  m3 <- vglmer(y ~ v_s(x2) + v_s(x, by = f), data = dat, 
               family = 'binomial',
               control = vglmer_control(iterations = NITER, print_prog = 20, prior_variance = 'mean_exists'))
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

  m8 <- vglmer(y ~ v_s(x2, by = g) + v_s(x, by = f, by_re = FALSE), 
               control = vglmer_control(
                 iterations = NITER, print_prog = 20),
               data = dat, family = 'binomial')
  
  # Check that they are all accepted
  expect_equal(m1$ELBO$accepted_PX, 1)
  expect_equal(m2$ELBO$accepted_PX, 1)
  
  # Check ELBO increases
  expect_gt(min(diff(m1$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m2$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m3$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m4$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m5$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m6$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m7$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  expect_gt(min(diff(m8$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
  
  for (v in c(TRUE, FALSE)){
    dat_custom <- data.frame(x = mean(dat$x), x2 = mean(dat$x2), 
       g = c('CA', 'AL'), f = c('e', 'a'),
       stringsAsFactors = v)
    expect_vector(predict(m1, newdata = dat_custom))
    expect_vector(predict(m2, newdata = dat_custom))
    expect_vector(predict(m3, newdata = dat_custom))
    expect_vector(predict(m4, newdata = dat_custom))
    expect_vector(predict(m5, newdata = dat_custom))
    expect_vector(predict(m6, newdata = dat_custom))
    expect_vector(predict(m7, newdata = dat_custom))
    expect_vector(predict(m8, newdata = dat_custom))
  }
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
  
  dat$y <- rbinom(100, 1, plogis(sin(3 * dat$x) * rnorm(5)[match(dat$f, letters)]))
  
  if (env_test != 'CRAN'){
    loop_vector <- c('tpf', 'o')
  }else{
    loop_vector <- 'tpf'
  }
  
  for (loop_type in loop_vector){
    message(loop_type)
    m1 <- vglmer(y ~ x2 + v_s(x, type = loop_type), 
                 data = dat, family = 'linear',
                 control = vglmer_control(iterations = NITER))
    
    # Check that prediction works for simple spline case
    predict_dat <- expand.grid(x = seq(min(dat$x), max(dat$x), length.out=100), x2 = 0)
    predict_vglmer <- predict(m1, newdata = predict_dat)
    # Check that *fitted* aligns with *predict*  
    expect_equivalent(fitted(m1), predict(m1, newdata = dat))
    # Check that there is some slope in the function
    diff_in_diff <- range(diff(diff(predict_vglmer)))
    expect_gte(max(abs(diff_in_diff)), sqrt(.Machine$double.eps))
    
    raw_X <- vglmer_build_spline(x = predict_dat$x, knots = m1$internal_parameters$spline$attr[[1]]$knots,
                                 type = m1$internal_parameters$spline$attr[[1]]$type,
                                 Boundary.knots = m1$internal_parameters$spline$attr[[1]]$Boundary.knots, by = NULL)[[1]]$x
    term_1 <- raw_X %*% m1$alpha$mean
    term_2 <- cbind(1, 0, predict_dat$x) %*% m1$beta$mean
    direct_predict <- as.vector(term_1 + term_2)
    expect_equivalent(direct_predict, predict_vglmer)   
    terms_predict <- predict(m1, newdata = predict_dat, type = 'terms')
    expect_equivalent(term_1[,1], terms_predict[, 'spline-x-1-base'])
    expect_equivalent(term_2[,1], terms_predict[, 'FE'])
  }

  if (TRUE){
    
    for (loop_type in loop_vector){
      message(loop_type)
      m1 <- vglmer(y ~ x2 + v_s(x, by = f, type = loop_type), 
                   data = dat, family = 'linear',
                   control = vglmer_control(iterations = NITER))
      
      # Check that prediction works for simple spline case
      predict_dat <- expand.grid(x = seq(min(dat$x), max(dat$x), length.out=100), x2 = 0, f = c('a', 'd', 'e', 'c'))
      predict_vglmer <- predict(m1, newdata = predict_dat)
      # Check that *fitted* aligns with *predict*  
      expect_equivalent(fitted(m1), predict(m1, newdata = dat))
      # Check that there is some slope in the function
      diff_in_diff <- sapply(split(predict_vglmer, predict_dat$f), FUN=function(i){
        diff_in_diff <- diff(diff(i))/diff(range(i))
        return(max(abs(diff_in_diff)))        
      })
      expect_true(all(diff_in_diff > sqrt(.Machine$double.eps)))
      
      raw_X <- vglmer_build_spline(x = predict_dat$x, knots = m1$internal_parameters$spline$attr[[1]]$knots,
       type = m1$internal_parameters$spline$attr[[1]]$type,
       by = predict_dat$f,
       Boundary.knots = m1$internal_parameters$spline$attr[[1]]$Boundary.knots)[[1]]$x
      
      expand_x <- do.call('cbind', lapply(c('a', 'c', 'd', 'e'), FUN=function(i){
        drop0(Diagonal(x = predict_dat$f == i) %*% raw_X)
      }))
      all_X <- cbind(raw_X, expand_x)
      
      direct_predict_spline <- as.vector(all_X %*% m1$alpha$mean[grep(rownames(m1$alpha$mean), pattern='spline @ x @ (base|[adec])'),] + 
                cbind(1, 0, predict_dat$x) %*% m1$beta$mean)
      
      wide_alpha <- matrix(m1$alpha$mean[grepl(rownames(m1$alpha$mean), pattern='^f'),,drop=F], byrow = TRUE, ncol = 2)
      
      direct_predict <- rowSums(wide_alpha[match(predict_dat$f, c('a', 'b', 'c', 'd', 'e')),] * cbind(1, predict_dat$x))
      expect_equivalent(direct_predict + direct_predict_spline, predict_vglmer)    
    }
  }

})

test_that("Compare against mgcv", {
  
  if (env_test == 'local'){
    N <- 1000
    dat <- data.frame(x = c(-3.1, rnorm(N-1)), x2 = rexp(N),
                      g = sample(state.abb[1:5], N, replace = T),
                      f = sample(letters[1:5], N, replace = T))
    
    dat$y <- rbinom(N, 1, plogis(2 * sin(3 * dat$x) + rnorm(5)[match(dat$f, letters)]))
    
    dat$f <- factor(dat$f)
    est_gam <- mgcv::gam(y ~ s(x, bs = 'bs') + s(f, bs = 're'), method = 'REML', data = dat, family = binomial())  
    est_vglmer <- vglmer(y ~ v_s(x, type = 'o') + (1 | f), data = dat, family = 'binomial')
    
    gx <- seq(quantile(dat$x, 0.2), quantile(dat$x, 0.8), length.out=100)
    pred_estgam <- predict(est_gam, newdata = data.frame(x = gx, f = 'b'))
    pred_estvglmer <- predict(est_vglmer, newdata = data.frame(x = gx, f = 'b'))
    
    expect_gte(cor(plogis(pred_estvglmer), plogis(pred_estgam)), 0.95)
    
    est_gam <- mgcv::gam(y ~ s(x), data = dat, method = 'REML', family = gaussian())  
    est_vglmer <- vglmer(y ~ v_s(x, type = 'o'), data = dat, family = 'linear')
    
    knots_custom <- seq(-3, 3, length.out=35)
    est_vglmer_2 <- expect_warning(
      vglmer(y ~ v_s(x, knots = knots_custom, type = 'tpf'), data = dat, family = 'linear'),
      regexp = 'observed data is outside of the self-provided knots'
    )
    pred_estgam <- predict(est_gam, newdata = data.frame(x = gx, f = 'b'))
    pred_estvglmer <- predict(est_vglmer, newdata = data.frame(x = gx, f = 'b'))
    pred_estvglmer_2 <- predict(est_vglmer_2, newdata = data.frame(x = gx, f = 'b'))
    
    expect_gte(cor(pred_estvglmer, pred_estgam), 0.95)
    expect_gte(cor(pred_estgam, pred_estvglmer_2), 0.95)
    
  }
})

test_that("Predict tests for tpf; predict outside of boundary", {
  
  N <- 100
  dat <- data.frame(x = rnorm(N), x2 = rexp(N),
                    g = sample(state.abb[1:5], N, replace = T),
                    f = sample(letters[1:5], N, replace = T))
  
  dat$y <- rbinom(N, 1, plogis(2 * sin(3 * dat$x) + rnorm(5)[match(dat$f, letters)]))
  
  est_vglmer_2 <- vglmer(y ~ v_s(x, type = 'tpf'), 
                         data = dat, family = 'linear', control = vglmer_control(iterations = NITER))
  knots_custom <- est_vglmer_2$internal_parameters$spline$attr[[1]]$knots
  predict_at_knots <- predict(est_vglmer_2, newdata = data.frame(x = knots_custom))
  
  sum(est_vglmer_2$alpha$mean)
  
  diff_rhs <- diff(predict(est_vglmer_2, newdata = data.frame(x = 5:10)))
  expect_lte(max(abs(diff(diff_rhs))), sqrt(.Machine$double.eps))
  
  diff_lhs <- diff(predict(est_vglmer_2, newdata = data.frame(x = -15:-10)))
  expect_lte(max(abs(diff(diff_lhs))), sqrt(.Machine$double.eps))
  # to the "left", the extrapolated value should be the slope of the linear part
  expect_equivalent(mean(diff_lhs), coef(est_vglmer_2)[2])
  # to the right, extrapolated value should be the sum of all 
  # the RE slopes + the baseline slope
  expect_equivalent(
    coef(est_vglmer_2)[2] + sum(ranef(est_vglmer_2)[[1]][,2]),
    mean(diff_rhs)
  )
  # The slope between any two knots should be the cumulative sum
  # of the intermediate slopes from the RE
  diff_predict_at_knots <- diff(predict_at_knots)
  slope_predict_at_knots <- diff_predict_at_knots/diff(knots_custom)
  expect_equivalent(slope_predict_at_knots, coef(est_vglmer_2)[2] + cumsum(est_vglmer_2$alpha$mean[-nrow(est_vglmer_2$alpha$mean),1]))
  
  
})

# TO-DO tests
# check decomposition of D and then re-transformation when
# implemented
