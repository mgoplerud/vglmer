

# build_alpha_expansion <- function(mean_alpha, mean_beta, re_names, number_of_Rlevels_per_RE, variables_per_RE, running_per_RE){
#   adjust_alpha_mu <- c()
#   adjust_beta_mu <- c()
#
#   re_names <- unlist(re_names)
#   for (v in 1:n_REs){
#     mean_id_v <- rep(1:levels_per_RE[v], each = variables_per_RE[v])
#
#     mean_mu <- colMeans(do.call('rbind', split(mean_alpha[(1+running_per_RE[v]):running_per_RE[v+1]], mean_id_v)))
#     adjust_beta_mu <- c(adjust_beta_mu, mean_mu)
#     adjust_alpha_mu <- c(adjust_alpha_mu, rep(mean_mu, levels_per_RE[v]))
#   }
#
#   mean_alpha <- mean_alpha - adjust_alpha_mu
#
#   adjust_beta_mu <- sapply(split(adjust_beta_mu, re_names), sum)
#   adjust_beta_mu <- adjust_beta_mu[match(rownames(mean_beta), names(adjust_beta_mu))]
#   adjust_beta_mu[is.na(names(adjust_beta_mu))] <- 0
#   mean_beta <- mean_beta + adjust_beta_mu
#
#   return(list(alpha = mean_alpha, beta = mean_beta))
# }
