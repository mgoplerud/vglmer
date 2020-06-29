internal_row_means <- function(x){if (is.matrix(x)){rowMeans(x)}else{mean(x)}}

create_vi_mu_j <- function(RE_to_FE_lookup, unusual_RE, M_prime, outer_alpha_RE_positions,
                           vi_alpha_mean, vi_sigma_alpha, vi_sigma_alpha_nu, d_j){
  E_sigma_inv <- mapply(vi_sigma_alpha, vi_sigma_alpha_nu, d_j, SIMPLIFY = FALSE, FUN=function(phi, nu, d){
    inv_phi <- solve(phi)
    sigma.inv <- nu * inv_phi
    return(sigma.inv)
  })
  
  
  
  if (!unusual_RE){
    return(vi_mu_j <- t(M_prime) %*% vi_alpha_mean)
  }
  
  vi_mu_j <- mapply(E_sigma_inv, RE_to_FE_lookup, d_j, outer_alpha_RE_positions, 
    SIMPLIFY = FALSE, FUN=function(Esigmaj, RE_to_FE_j, nu, pos_j){
      na_ids <- which(is.na(RE_to_FE_j))
      non_na_ids <- setdiff(1:nu, na_ids)
      
      if (length(non_na_ids) == 0){
        return(list(adj_mu = matrix(0, nu), mu_1 = matrix(0, nu)))
      }
      
      m1 <- sapply(pos_j, FUN=function(i){vi_alpha_mean[i[non_na_ids]]})
      m1 <- internal_row_means(m1)
      
      if (length(na_ids) == 0){
        return(list(adj_mu = m1, mu_1 = m1))
      }
      
      m2 <- sapply(pos_j, FUN=function(i){vi_alpha_mean[i[na_ids]]})
      m2 <- internal_row_means(m2)
      
      S11 <- Esigmaj[non_na_ids, non_na_ids, drop = F]
      S12 <- Esigmaj[non_na_ids, na_ids, drop = F]
      
      proj_mu <- m1 + solve(S11) %*% S12 %*% m2
      mu_output <- matrix(0, nu)
      mu_output[non_na_ids,] <- proj_mu
      return(list(adj_mu = mu_output, mu_1 = m1))
    })
  
  vi_mu_j <- do.call('rbind', lapply(vi_mu_j, FUN=function(i){i$adj_mu}))
  return(vi_mu_j)
}
