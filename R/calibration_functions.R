

get_edf_uni <- function(alpha, diag_z, target){
  (sum(diag_z/(diag_z + exp(alpha))) - target)^2
}
get_edf_diag <- function(alpha, M, g){
  sum(diag(solve(Cholesky(M + kronecker(Diagonal(n = g), Diagonal(x = exp(alpha)))), M)))
}
opt_edf <- function(par, target, M, g_l){abs(get_edf_diag(alpha = par, M = M, g = g_l) - target)^2}


calibrate_init <- function(Z, weight_init, cyclical_pos, d_j, g_j, vi_sigma_alpha_nu, vi_sigmasq_a, vi_sigmasq_b){
  
  vi_sigma_alpha <- mapply(cyclical_pos, names(d_j), d_j, g_j, vi_sigma_alpha_nu, FUN=function(c_l, n_l, d_l, g_l, alpha_nu){
    
    ZtZ <- crossprod(Diagonal(x = sqrt(weight_init)) %*% Z[,c_l])
    if (d_l == 1 & isDiagonal(ZtZ)){
      diag_ZtZ <- diag(ZtZ)
      out_cl <- optim(par = 0, fn = get_edf_uni, 
                      diag_z = diag_ZtZ, target = ncol(ZtZ) * 0.90,
                      method = 'L-BFGS-B')
    }else{
      # If spline (i.e. non-diagonal), then set to relatively small to start
      if (d_l == 1){
        if (grepl(n_l, pattern='[0-9]-int$')){
          target_fraction <- 1/ncol(ZtZ)
        }else{
          target_fraction <- 1/ncol(ZtZ)
        }
      }else{target_fraction <- 0.90}
      out_cl <- optim(par = rep(0, d_l), fn = opt_edf, method = 'L-BFGS-B',
                      control = list(maxit = 10),
                      M = ZtZ, g_l = g_l, target = target_fraction * ncol(ZtZ))
    }
    out_sa <- Diagonal(x = alpha_nu * exp(-out_cl$par) * vi_sigmasq_a/vi_sigmasq_b)
    return(out_sa)
  })
  
  return(vi_sigma_alpha)
}

