context('Test various update methods')

test_that('Joint vs Cyclical Update', {
  
  N <- 1000
  G <- 20
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  

  for (v in c("weak", 'partial', 'strong')){
    
    ex_vglmer_cyclic <- vglmer(formula = y ~ x + (1 | g), family = 'binomial',
        data = NULL, control = vglmer_control(factorization_method = v, linpred_method = 'cyclical', init = 'zero')
    )

    ex_vglmer_joint <- vglmer(formula = y ~ x + (1 | g), family = 'binomial',
       data = NULL, control = vglmer_control(factorization_method = v, linpred_method = 'joint', init = 'zero')
    )
    
    fmt_vglmer_cyclic <- format_vglmer(ex_vglmer_cyclic)
    fmt_vglmer_joint <- format_vglmer(ex_vglmer_joint)
    
    expect_equivalent(fmt_vglmer_cyclic, fmt_vglmer_joint, tolerance = 1e-5)

    if (v == 'strong'){
      ex_vglmer_normal <- vglmer(formula = y ~ x + (1 | g), family = 'binomial',
        data = NULL, control = vglmer_control(factorization_method = v, linpred_method = 'solve_normal', init = 'zero')
      )
      
      fmt_vglmer_normal <- format_vglmer(ex_vglmer_normal)
      expect_equivalent(fmt_vglmer_normal, fmt_vglmer_joint, tolerance = 1e-5)
    }
    
  }
  
  
})


test_that('Compare PX vs Non-PX', {
  
  N <- 1000
  G <- 20
  x <- rnorm(N)
  g <- sample(1:G, N, replace = T) 
  alpha <- rnorm(G)
  
  y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))
  
  
  ex_vglmer_px <- vglmer(formula = y ~ x + (1 | g), family = 'binomial',
   data = NULL, control = vglmer_control(factorization_method ='strong', 
    parameter_expansion = 'mean', debug_ELBO = T, init = 'zero')
  )
  
  ex_vglmer_no <- vglmer(formula = y ~ x + (1 | g), family = 'binomial',
    data = NULL, control = vglmer_control(factorization_method = 'strong', 
      parameter_expansion =  'none', debug_ELBO = T, init = 'zero')
  )
  
  expect_gte(min(diff(ex_vglmer_no$ELBO_trajectory$ELBO)), - sqrt(.Machine$double.eps))
  expect_gte(min(diff(ex_vglmer_px$ELBO_trajectory$ELBO)),  - sqrt(.Machine$double.eps))
  
  fmt_vglmer_px <- format_vglmer(ex_vglmer_px)
  fmt_vglmer_no <- format_vglmer(ex_vglmer_no)
  
  expect_equivalent(fmt_vglmer_px, fmt_vglmer_no, tolerance = 1e-5)

})

test_that('Compare VI r methods', {
  
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
  for (v in c("VEM", "delta", "Laplace")){
    
    example_vglmer <- vglmer(formula = y ~ x + (1 | g), data = data, 
                             family = 'negbin', 
                             control = vglmer_control(factorization_method = 'strong', 
                                                      vi_r_method = v, init = 'random'))    
    #Test whether it monotonically increases   
    if (v == 'VEM'){
      expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), -sqrt(.Machine$double.eps))
    }else{
      # "Rough" for approximate VI surrogates.
      expect_gte(min(diff(example_vglmer$ELBO_trajectory$ELBO)), -0.01)
    }
    fmt_vglmer <- format_vglmer(example_vglmer)
    names(fmt_vglmer)[-1] <- paste0(v, '_', names(fmt_vglmer)[-1])
    list_output[[v]] <- fmt_vglmer
    list_r[[v]] <- example_vglmer$r
  }
  
  list_output <- Reduce(function(a,b){merge(a,b, by = 'name')}, list_output)
  expect_gte(min(as.vector(cor(list_output[,c('Laplace_mean', 'delta_mean', 'VEM_mean')]))), 0.95)
  expect_gte(min(as.vector(cor(list_output[,c('Laplace_var', 'delta_var', 'VEM_var')]))), 0.95)
  
  all_r <- sapply(list_r, FUN=function(i){i$mu})
  #Check that mu are quite close
  expect_lte(diff(range(all_r)), 0.01)
  #Check that the mu standard errors are close for Laplace/delta
  expect_lte(diff(sqrt(sapply(list_r, FUN=function(i){i$sigma}))[-1]), 0.01)
  
})

