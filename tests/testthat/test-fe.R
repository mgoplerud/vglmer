warning('SETUP FE TESTS')
data('Fishing', package = 'mlogit')
Fish <- dfidx(Fishing, varying = 2:9, shape = "wide", choice = "mode")
Fish$person <- Fish$idx$id1
Fish$type <- Fish$idx$id2
Fish$bin_catch <- cut(Fish$catch, 3)


fit <- vglmer(as.numeric(mode) ~ income + (1 | bin_catch) +
                v_fe(person) + v_s(catch),
              data = Fish, family = 'linear',
              control = vglmer_control(do_SQUAREM = FALSE, parameter_expansion = 'mean'))

# This should work too...
# fit <- vglmer(as.numeric(mode) ~ type + 
#                 (1 | bin_catch) + v_fe(person) + v_s(income),
#        data = Fish, family = 'linear',
#        control = vglmer_control(do_SQUAREM = FALSE))

newFish <- tidyr::crossing(
  data.frame(bin_catch = levels(Fish$bin_catch)),
  data.frame(catch = mean(Fish$catch), 
             type = unique(Fish$type)),
  data.frame(income = seq(min(Fish$income), max(Fish$income), length.out=100))
) %>% data.frame
a <- predict(fit, newdata = newFish)
b <- predict(fit, newdata = newFish %>% dplyr::mutate(income = 0))

# est_1 <- vglmer(as.numeric(mode) ~ catch + v_s(person, type = 'fe') + (1 | choice), 
#                 data = Fish, family = 'linear', 
#                 control = vglmer_control(debug_px = T, 
#                                          quiet_rho = F))
# 
# 
# est_1 <- vglmer(as.numeric(mode) ~ v_s(price) + catch, data = Fish,
#                 family = 'binomial',
#                 control = vglmer_control(quiet_rho = F, debug_px = T))
# 
# 
# 
# # FAILS to increase objective if poisson
# est_2a <- vglmer(as.numeric(mode) ~ catch + v_s(person, type = 'fe'), 
#                 data = Fish, family = 'poisson', 
#                 control = vglmer_control(debug_px = F, iterations = 10,
#                                          quiet_rho = F))
# ranef(est_2a) # FAILS
# est_2a <- vglmer(as.numeric(mode) ~ catch + v_s(person, type = 'fe'), 
#                  data = Fish, family = 'binomial', 
#                  control = vglmer_control(debug_px = F, iterations = 10,
#                                           quiet_rho = F))
# 
# est_2a <- vglmer(as.numeric(mode) ~ catch + v_s(person, type = 'fe'), 
#                  data = Fish, family = 'poisson', 
#                  control = vglmer_control(debug_px = F, iterations = 10,
#                                           quiet_rho = F))
# 
# 
# logit_vglmer <- vglmer(as.numeric(mode) ~ 0 + price + catch + 
#                          (income):choice + (1 | choice), 
#                        data = Fish %>% mutate(choice = factor(choice, levels = c('boat', 'charter', 'pier', 'beach'))), 
#                        family = 'binomial',
#                        control = vglmer_control(quiet_rho = FALSE))
# # Check that FE confirms with predictions for simple case...
# 
# 
# fit_1 <- vglmer(as.numeric(mode) ~ price + v_s(income), 
#        data = Fish,
#        family = 'binomial',
#        control = vglmer_control(quiet_rho = FALSE))
# 
# predict(fit_1, newdata = Fish)
# 
# 
# pois_vglmer <- vglmer(as.numeric(mode) ~ 1 + v_s(person, type = 'fe') + (1 | choice),
#     data = Fish %>% filter(person %in% 1:50), 
#     family = 'poisson', 
#     control = vglmer_control(do_SQUAREM = FALSE, parameter_expansion = 'mean')
# )
# 
# plot(ELBO(pois_vglmer, 'traj')[-1:-5], type = 'l')
# 
# plot(ELBO(pois_vglmer, 'traj')[800:820], type = 'l')
# 
# 
# 
# w2 <- vglmer(as.numeric(mode) ~ 1 + v_s(price) + v_s(person, type = 'fe') + (1 + income | choice),
#     data = Fish %>% filter(person %in% 1:100), 
#     family = 'poisson', 
#     control = vglmer_control(factorization_method = 'weak', debug_ELBO = TRUE,
#                              do_SQUAREM = FALSE, parameter_expansion = 'mean')
# )
# 
# plot(ELBO(w2, 'traj')[-1:-5], type = 'l')
# 
# c(1,2,3,4,6,7,8,11,12,16)
# 
# build_D <- function(M){
#     elements <- nrow(M) * (nrow(M) + 1)/2
#     D <- t(sparseMatrix(i = seq(elements), j = do.call('c', sapply(seq(nrow(M)), FUN=function(i){seq(nrow(M) + 1 - i) + (i - 1) * (1 + nrow(M))})), x = 1))
#     return(D)
# }
# 
# test_D <- 1000
# M <- matrix(1:test_D^2, nrow = test_D)
# D <- build_D(M)
# t(D) %*% kronecker(M, Diagonal(n = ncol(M))) %*% D
