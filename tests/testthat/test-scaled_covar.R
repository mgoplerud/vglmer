scale_Fish <- Fish %>% ungroup %>% mutate(price = as.vector(scale(price)))

est_scaled <- vglmer(as.numeric(mode) ~ 1 + v_s(price) + (1 + income | choice),
   data = scale_Fish, 
   family = 'poisson', control = vglmer_control(factorization_method = 'weak', do_SQUAREM = FALSE, parameter_expansion = 'mean')
)
est_unscaled <- vglmer(as.numeric(mode) ~ 1 + v_s(price) + (1 + income | choice),
   data = Fish , 
   family = 'poisson', control = vglmer_control(factorization_method = 'weak', do_SQUAREM = FALSE, parameter_expansion = 'mean')
)


ELBO(est_scaled)
ELBO(est_unscaled)

plot(predict(est_scaled, newdata = scale_Fish), predict(est_unscaled, newdata = Fish))




est_unscaled <- vglmer(as.numeric(mode) ~ 1 + v_s(price, knots ) + (1 + income | choice),
   data = Fish , 
   family = 'poisson', 
   control = vglmer_control(factorization_method = 'weak', do_SQUAREM = FALSE, parameter_expansion = 'mean')
)

# Test that this INCREASES

m1 <- vglmer(as.numeric(mode) ~ v_s(scaled_price), 
             data = Fish , family = 'binomial',
             control = vglmer_control(return_data = F, prior_variance = 'mean_exists',
                                      factorization_method = 'strong',
                                      do_SQUAREM = T))

