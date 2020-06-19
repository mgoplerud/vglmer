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
