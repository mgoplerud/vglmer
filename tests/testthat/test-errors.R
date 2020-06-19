context('Test vglmer robustness to certain situations')

test_that('vglmer can run with objects in environment', {
  
  N <- 1000
  G <- 100
  G_names <- paste(sample(letters, G, replace = T), 1:G)
  x <- rnorm(N)
  g <- sample(G_names, N, replace = T)
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[match(g, G_names)]))
  
  test_nodata <- tryCatch(suppressMessages(vglmer(y ~ x + (1 | g), 
      data = NULL, 
      control = vglmer_control(init = 'zero', 
        iterations = 1, print_prog = 10),
      family = 'binomial')), 
    error = function(e){NULL}
  )
  expect_false(is.null(test_nodata))

  dta <- data.frame(Y = y, X = x, G = g)
  #Inject missingness into 
  dta$Y[38] <- NA
  dta$X[39] <- NA
  dta$G[190] <- NA
  dta[108,] <- NA
  test_missing <- tryCatch(suppressMessages(vglmer(Y ~ X + (1 | G), 
    data = dta, 
    control = vglmer_control(init = 'zero', return_data = T,
                             iterations = 1, print_prog = 10),
    family = 'binomial')), 
    error = function(e){NULL}
  )
  #Confirm runs
  expect_false(is.null(test_missing))
  #Confirms deletion "works"
  expect_equivalent(dta$X[-c(38, 39,108, 190)], test_missing$data$X[,2])
  expect_equivalent(dta$Y[-c(38, 39,108, 190)], test_missing$data$y)

})