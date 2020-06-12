context('Cyclical vs Joint Comparison')

test_that('Joint vs Cyclical Update', {
  
  set.seed(02134)
  
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
