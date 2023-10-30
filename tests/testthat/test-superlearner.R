context("SuperLearner tests")

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  NITER <- 2
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  NITER <- 2
  env_test <- "CRAN"
  set.seed(181)
}else{
  # If on local machine
  NITER <- 2000
  env_test <- 'local'
}

test_that("Test SuperLearner", {
  
  N <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N, size = 1, prob = .2)
  f <- sample(letters[1:3], N, replace = T)
  y <- rbinom(N, 1, plogis(x1 + runif(3)[match(f, letters)]))
  
  
  if (requireNamespace("SuperLearner", quietly = TRUE)) {
    require(SuperLearner)
    
    sl_k <- add_formula_SL('SL.glm', env = globalenv())
    sl_k <- function(...){SL.glm_f(formula = ~ x1 + x2, ...)}
    sl_m <- function(...) {
      suppressMessages(SL.vglmer(formula = ~ x1 + (1 | f), 
        control = vglmer_control(iterations = 2), ...))
    }
    sl_g <- function(...) {
      suppressMessages(SL.glmer(formula = ~ x1 + (1 | f), ...))
    }
    fit_SL <- SuperLearner::SuperLearner(
      Y = y, cvControl = list(V = 2),
      X = data.frame(x1 = x1, x2 = x2, f = f),
      SL.library = c("sl_m", "sl_g", "sl_k")
    )
    expect_s3_class(fit_SL, "SuperLearner")

    pred <- predict(fit_SL, newdata = data.frame(x1 = x1, x2 = x2, f = f))
    expect_length(pred, n = 2)
  }
})
