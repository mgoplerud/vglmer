

update_rho <- function(XR, y, omega, prior_precision, 
                       moments_sigma_alpha,
                       prior_sigma_alpha_nu, prior_sigma_alpha_phi,
                       vi_a_a_jp, vi_a_b_jp, vi_a_nu_jp,
                       vi_a_APRIOR_jp,
                       spline_REs, vi_beta_mean,
                       p.X, d_j, stationary_rho,
                       do_huangwand,
                       px_it = NULL, init_rho = NULL,
                       method){
  
  if (do_huangwand){
    prior_weight <- vi_a_nu_jp + d_j - 1
    diag_weight <- mapply(vi_a_a_jp, vi_a_b_jp, vi_a_nu_jp, SIMPLIFY = FALSE, 
      FUN = function(tilde.a, tilde.b, nu) {
        Diagonal(x = tilde.a/tilde.b) * 2 * nu
      })
  }else{
    diag_weight <- prior_sigma_alpha_phi
    prior_weight <- prior_sigma_alpha_nu
  }
  
  ESigma <- lapply(moments_sigma_alpha[c(which(spline_REs), which(!spline_REs))], FUN=function(i){i$sigma.inv})
  Phi <- diag_weight[c(which(spline_REs), which(!spline_REs))] 
  nu <- prior_weight[c(which(spline_REs), which(!spline_REs))]
  
  y <- as.vector(y)
  sum_ysq <- sum(omega %*% y)
  tXy <- t(XR) %*% y
  tXX <- t(XR) %*% omega %*% XR
  rho_idx <- d_j[spline_REs]
  rho_idx <- rho_idx * seq_len(length(rho_idx))
  rho_idx <- c(rho_idx, rep(seq_len(sum(!spline_REs)), times = d_j[!spline_REs]^2) + sum(spline_REs))

  if (method == 'numerical'){

    
    null_rho <- c(rep(1, sum(spline_REs)), stationary_rho)
    null_rho <- c(as.vector(vi_beta_mean), null_rho)
    dim_rho <- c(rep(1, sum(spline_REs)), d_j[!spline_REs])
    
    ctrl_opt <- list(fnscale = -1)
    if (!is.null(px_it)){
      ctrl_opt$maxit <- px_it
    }
    if (is.null(init_rho)){
      init_rho <- null_rho
    }
    
    
    opt_rho <- optim(par = init_rho, fn = eval_rho, 
        gr = eval_grad_rho,
        method = 'L-BFGS-B', control = ctrl_opt,
        tXy = tXy, tXX = tXX,
        ridge = prior_precision, rho_idx = rho_idx,
        nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X
    )
    null_eval <- eval_rho(null_rho, tXy = tXy, tXX = tXX, ridge = prior_precision, rho_idx = rho_idx,
             nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X)
    if (opt_rho$value < null_eval){
      warning('Optimization failed in parameter expansion.')
      opt_rho$par <- null_rho
    }
    
    improvement <- opt_rho$value - null_eval
    opt_rho <- opt_rho$par
    names(opt_rho) <- NULL
    opt_rho <- list(rho = opt_rho, improvement = improvement)
    
  }else if (method == 'dynamic'){
    
    vec_OSL_prior <- mapply(moments_sigma_alpha[!spline_REs], 
                            diag_weight[!spline_REs], 
                            prior_weight[!spline_REs], 
                            SIMPLIFY = FALSE, FUN=function(moment_j, phi_j, nu_j){
                              as.vector(moment_j$sigma.inv %*% phi_j - nu_j * Diagonal(n = nrow(phi_j)))
                            })
    vec_OSL_prior <- do.call('c', vec_OSL_prior)
    
    if (sum(spline_REs)){
      OSL_spline_prior <- unlist(mapply(moments_sigma_alpha[spline_REs], 
                                        diag_weight[spline_REs], 
                                        prior_weight[spline_REs], 
                                        SIMPLIFY = FALSE, FUN=function(moment_j, phi_j, nu_j){
                                          as.vector(moment_j$sigma.inv %*% phi_j - nu_j * Diagonal(n = nrow(phi_j)))
                                        }))
      vec_OSL_prior <- matrix(c(rep(0, p.X), OSL_spline_prior, vec_OSL_prior))
    }else{
      vec_OSL_prior <- matrix(c(rep(0, p.X), vec_OSL_prior))
    }
    
    hw_a <- vi_a_a_jp[c(which(spline_REs), which(!spline_REs))]
    A_prior <- vi_a_APRIOR_jp[c(which(spline_REs), which(!spline_REs))]
    nu_prior <- vi_a_nu_jp[c(which(spline_REs), which(!spline_REs))]
    
    sum_d <- sum(d_j)
    
    null_rho <- c(rep(1, sum(spline_REs)), stationary_rho)
    null_rho <- c(as.vector(vi_beta_mean), null_rho)
    dim_rho <- c(rep(1, sum(spline_REs)), d_j[!spline_REs])
    
    if (is.null(init_rho)){
      init_rho <- null_rho
    }
    ctrl_opt <- list(fnscale = -1)
    if (!is.null(px_it)){
      ctrl_opt$maxit <- px_it
    }
    
    null_eval <- eval_profiled_rho(null_rho, tXy = tXy, tXX = tXX, ridge = prior_precision, rho_idx = rho_idx,
       nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X,
       sum_d = sum_d, hw_a = hw_a, A_prior = A_prior,
       nu_prior = nu_prior)
    
    OSL_rho <- vecR_fast_ridge(X = XR, 
                               omega = omega, prior_precision = prior_precision, y = y, 
                               adjust_y = as.vector(vec_OSL_prior)) 
    OSL_eval <- eval_profiled_rho(rho = OSL_rho, tXy = tXy, tXX = tXX, ridge = prior_precision, rho_idx = rho_idx,
                                  nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X, hw_a, sum_d = sum_d,
                                  A_prior = A_prior, nu_prior = nu_prior)
    OSL_improvement <- OSL_eval - null_eval
    
    if (OSL_improvement > 0){
      opt_rho <- OSL_rho
      improvement <- OSL_improvement
      # print('OSL')
      # print(c(NA, improvement))
    }else{
      # print('max')
      opt_rho <- optim(par = null_rho, fn = eval_profiled_rho, 
                       gr = eval_grad_profiled_rho,
                       method = 'L-BFGS-B', control = ctrl_opt,
                       tXy = tXy, tXX = tXX,
                       ridge = prior_precision, rho_idx = rho_idx,
                       nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X,
                       sum_d = sum_d, hw_a = hw_a, A_prior = A_prior,
                       nu_prior = nu_prior
      )
      improvement <- opt_rho$value - null_eval
      
      # compare_improvement <- c(improvement, OSL_improvement)
      # names(compare_improvement) <- c('max', 'OSL')
      # print(compare_improvement)
      # print(compare_improvement/compare_improvement['max'])
      opt_rho <- opt_rho$par
    }
    
    if (improvement < 0){
      warning('Optimization of parameter expansion failed')
      opt_rho <- null_rho
    }
    
    raw_opt_rho <- opt_rho

    if (p.X > 0){nonfe_rho <- opt_rho[-seq_len(p.X)]}else{nonfe_rho <- opt_rho}
    Rmatrix <- mapply(split(nonfe_rho, rho_idx), dim_rho, SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d, ncol = d)})
    opt_rho_hw <- mapply(Rmatrix, nu, hw_a, A_prior, ESigma, nu_prior, SIMPLIFY = FALSE, 
                         FUN=function(R_j, nu_j, hw_a_j, A_j, ESigma.inv.j, nu_prior_j){
                           inv_R_j <- solve(R_j)
                           diag_meat <- diag(t(inv_R_j) %*% ESigma.inv.j %*% inv_R_j)
                           rho_hw_j <- nu_prior_j * diag_meat + 1/A_j^2
                           return(rho_hw_j)
                         })
    names(opt_rho_hw) <- names(d_j)[c(which(spline_REs), which(!spline_REs))]
    names(opt_rho) <- NULL
    opt_rho <- list(hw = opt_rho_hw,
                    rho = opt_rho, improvement = improvement,
                    opt_par = raw_opt_rho)
  }else if (method == 'profiled'){
    
    ctrl_opt <- list(fnscale = -1)
    if (!is.null(px_it)){
      ctrl_opt$maxit <- px_it
    }
    
    hw_a <- vi_a_a_jp[c(which(spline_REs), which(!spline_REs))]
    A_prior <- vi_a_APRIOR_jp[c(which(spline_REs), which(!spline_REs))]
    nu_prior <- vi_a_nu_jp[c(which(spline_REs), which(!spline_REs))]
    
    sum_d <- sum(d_j)
    
    null_rho <- c(rep(1, sum(spline_REs)), stationary_rho)
    null_rho <- c(as.vector(vi_beta_mean), null_rho)
    dim_rho <- c(rep(1, sum(spline_REs)), d_j[!spline_REs])
    
    if (is.null(init_rho)){
      init_rho <- null_rho
    }

    opt_rho <- optim(par = null_rho, fn = eval_profiled_rho, 
                     gr = eval_grad_profiled_rho,
                     method = 'L-BFGS-B', control = ctrl_opt,
                     tXy = tXy, tXX = tXX,
                     ridge = prior_precision, rho_idx = rho_idx,
                     nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X,
                     sum_d = sum_d, hw_a = hw_a, A_prior = A_prior,
                     nu_prior = nu_prior
    )
    
    null_eval <- eval_profiled_rho(null_rho, tXy = tXy, tXX = tXX, ridge = prior_precision, rho_idx = rho_idx,
                 nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X,
                 sum_d = sum_d, hw_a = hw_a, A_prior = A_prior,
                 nu_prior = nu_prior)
    improvement <- opt_rho$value - null_eval
    if (opt_rho$value < null_eval){
      warning('Optimization of parameter expansion failed')
      opt_rho$par <- null_rho
    }
    opt_rho <- raw_opt_rho <- opt_rho$par
    if (p.X > 0){nonfe_rho <- opt_rho[-seq_len(p.X)]}else{nonfe_rho <- opt_rho}
    Rmatrix <- mapply(split(nonfe_rho, rho_idx), dim_rho, SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d, ncol = d)})
    opt_rho_hw <- mapply(Rmatrix, nu, hw_a, A_prior, ESigma, nu_prior, SIMPLIFY = FALSE, 
      FUN=function(R_j, nu_j, hw_a_j, A_j, ESigma.inv.j, nu_prior_j){
           inv_R_j <- solve(R_j)
           diag_meat <- diag(t(inv_R_j) %*% ESigma.inv.j %*% inv_R_j)
           rho_hw_j <- nu_prior_j * diag_meat + 1/A_j^2
           return(rho_hw_j)
      })
    names(opt_rho_hw) <- names(d_j)[c(which(spline_REs), which(!spline_REs))]
    names(opt_rho) <- NULL
    opt_rho <- list(hw = opt_rho_hw,
                    rho = opt_rho, improvement = improvement,
                    opt_par = raw_opt_rho)
    
  }else if (method == 'OSL'){
    
    null_rho <- c(rep(1, sum(spline_REs)), stationary_rho)
    null_rho <- c(as.vector(vi_beta_mean), null_rho)
    dim_rho <- c(rep(1, sum(spline_REs)), d_j[!spline_REs])

    vec_OSL_prior <- mapply(moments_sigma_alpha[!spline_REs], 
                            diag_weight[!spline_REs], 
                            prior_weight[!spline_REs], 
                            SIMPLIFY = FALSE, FUN=function(moment_j, phi_j, nu_j){
                              as.vector(moment_j$sigma.inv %*% phi_j - nu_j * Diagonal(n = nrow(phi_j)))
                            })
    vec_OSL_prior <- do.call('c', vec_OSL_prior)
    
    if (sum(spline_REs)){
      OSL_spline_prior <- unlist(mapply(moments_sigma_alpha[spline_REs], 
                                        diag_weight[spline_REs], 
                                        prior_weight[spline_REs], 
                                        SIMPLIFY = FALSE, FUN=function(moment_j, phi_j, nu_j){
                                          as.vector(moment_j$sigma.inv %*% phi_j - nu_j * Diagonal(n = nrow(phi_j)))
                                        }))
      vec_OSL_prior <- matrix(c(rep(0, p.X), OSL_spline_prior, vec_OSL_prior))
    }else{
      vec_OSL_prior <- matrix(c(rep(0, p.X), vec_OSL_prior))
    }
    
    OSL_rho <- vecR_fast_ridge(X = XR, 
     omega = omega, prior_precision = prior_precision, y = y, 
     adjust_y = as.vector(vec_OSL_prior)) 

    if (do_huangwand){
      sum_d <- sum(d_j)
      
      hw_a <- vi_a_a_jp[c(which(spline_REs), which(!spline_REs))]
      A_prior <- vi_a_APRIOR_jp[c(which(spline_REs), which(!spline_REs))]
      nu_prior <- vi_a_nu_jp[c(which(spline_REs), which(!spline_REs))]
      
      null_eval <- eval_profiled_rho(rho = null_rho, tXy = tXy, tXX = tXX, ridge = prior_precision, rho_idx = rho_idx,
                                     nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X, hw_a, sum_d = sum_d,
                                     A_prior = A_prior, nu_prior = nu_prior)
      OSL_eval <- eval_profiled_rho(rho = OSL_rho, tXy = tXy, tXX = tXX, ridge = prior_precision, rho_idx = rho_idx,
                                    nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X, hw_a, sum_d = sum_d,
                                    A_prior = A_prior, nu_prior = nu_prior)
    }else{
      null_eval <- eval_rho(null_rho, tXy = tXy, tXX = tXX, ridge = prior_precision, rho_idx = rho_idx,
                                     nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X)
      OSL_eval <- eval_rho(OSL_rho, tXy = tXy, tXX = tXX, ridge = prior_precision, rho_idx = rho_idx,
                                    nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X)
    }
    improvement <- OSL_eval - null_eval
    if (improvement < 0){
      OSL_rho <- null_rho
    }
    opt_rho <- OSL_rho
    opt_rho <- list(rho = opt_rho, improvement = improvement)
  }else{stop('..')}

  return(opt_rho)
}

eval_rho <- function(rho, tXy, tXX, ridge, rho_idx, nu, Phi, ESigma, dim_rho, p.X){
  ssr <- t(rho) %*% tXy - 1/2 * t(rho) %*% tXX %*% rho
  ridge <- -1/2 * as.numeric(t(rho) %*% ridge %*% rho)
  if (p.X > 0){nonfe_rho <- rho[-seq_len(p.X)]}else{nonfe_rho <- rho}
  Rmatrix <- mapply(split(nonfe_rho, rho_idx), dim_rho, SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d, ncol = d)})
  prior <- sum(mapply(Rmatrix, nu, Phi, ESigma, FUN=function(R_j, nu_j, Phi_j, ESigma.inv.j){
    inv_R_j <- solve(R_j)
    out <- - nu_j * determinant(R_j)$modulus - 1/2 * sum(Matrix::diag(inv_R_j %*% Phi_j %*% t(inv_R_j) %*% ESigma.inv.j))
    return(out)
  }))
  
  return(as.numeric( ssr + ridge + prior ) )
}

eval_grad_rho <- function(rho, tXy, tXX, ridge, rho_idx, nu, Phi, ESigma, dim_rho, p.X){
  ssr <- tXy - tXX %*% rho
  ridge <- - ridge %*% rho
  if (p.X > 0){nonfe_rho <- rho[-seq_len(p.X)]}else{nonfe_rho <- rho}
  Rmatrix <- mapply(split(nonfe_rho, rho_idx), dim_rho, SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d, ncol = d)})
  prior <- mapply(Rmatrix, nu, Phi, ESigma, SIMPLIFY = FALSE, FUN=function(R_j, nu_j, Phi_j, ESigma.inv.j){
    inv_R_j <- solve(R_j)
    inv_Phi_j <- solve(Phi_j)
    # meat <- inv_R_j %*% Phi_j %*% t(inv_R_j)
    meat <- solve(t(R_j) %*% inv_Phi_j %*% R_j)
    out <- as.vector(- nu_j * t(inv_R_j) + inv_Phi_j %*% R_j %*% meat %*% ESigma.inv.j %*% meat)
    return(out)
  })
  prior <- c(rep(0, p.X), unlist(prior))
  return(as.vector( ssr + ridge + prior ))
}

eval_rho_hw <- function(rho, tXy, tXX, ridge, rho_idx, nu, Phi, 
                        ESigma, dim_rho, p.X, dlist, sum_d, hw_a, A_prior, nu_prior){
  
  
  rho_hw <- exp(rho[seq_len(sum_d)])
  rho_hw <- split(rho_hw, dlist)
  rho <- rho[-seq_len(sum_d)]
  
  ssr <- t(rho) %*% tXy - 1/2 * t(rho) %*% tXX %*% rho
  ridge <- -1/2 * as.numeric(t(rho) %*% ridge %*% rho)
  if (p.X > 0){nonfe_rho <- rho[-seq_len(p.X)]}else{nonfe_rho <- rho}
  Rmatrix <- mapply(split(nonfe_rho, rho_idx), dim_rho, SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d, ncol = d)})
  
  prior_variance <- sum(mapply(Rmatrix, nu, hw_a, rho_hw, A_prior, ESigma, nu_prior,
     FUN=function(R_j, nu_j, hw_a_j, rho_hw_j, A_j, ESigma.inv.j, nu_prior_j){
       inv_R_j <- solve(R_j)
       # Prior from Wishart
       out <- - nu_j * determinant(R_j)$modulus - nu_j/2 * sum(log(rho_hw_j)) +
         -1/2 * sum(Matrix::diag(inv_R_j %*% Diagonal(x = 2 * nu_prior_j * hw_a_j/rho_hw_j) %*% t(inv_R_j) %*% ESigma.inv.j))
       # Prior from Inverse-Gamma and Entropy 
       out <- out + sum(-1/2 * log(rho_hw_j) - 1/A_j^2 * hw_a_j/rho_hw_j)
       return(out)
     }))
  return(as.numeric( ssr + ridge + prior_variance ) )
}



eval_profiled_rho <- function(rho, tXy, tXX, ridge, rho_idx, nu, Phi, 
                              ESigma, dim_rho, p.X, sum_d, hw_a, A_prior, nu_prior){
  
  ssr <- t(rho) %*% tXy - 1/2 * t(rho) %*% tXX %*% rho
  ridge <- -1/2 * as.numeric(t(rho) %*% ridge %*% rho)
  
  if (p.X > 0){nonfe_rho <- rho[-seq_len(p.X)]}else{nonfe_rho <- rho}
  
  Rmatrix <- mapply(split(nonfe_rho, rho_idx), dim_rho, SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d, ncol = d)})
  
  prior_variance <- sum(mapply(Rmatrix, nu, hw_a, A_prior, ESigma, nu_prior,
       FUN=function(R_j, nu_j, hw_a_j, A_j, ESigma.inv.j, nu_prior_j){
         
         inv_R_j <- solve(R_j)
         
         diag_meat <- diag(t(inv_R_j) %*% ESigma.inv.j %*% inv_R_j)
         
         rho_hw_j <- nu_prior_j * diag_meat + 1/A_j^2
         
         out <- -(nu_j + 1)/2 * sum(log(rho_hw_j)) +
           - nu_j * determinant(R_j)$modulus +
           sum(hw_a_j)
         return(out)
       }))
  return(as.numeric( ssr + ridge + prior_variance ) )
}


eval_grad_profiled_rho <- function(rho, tXy, tXX, ridge, rho_idx, nu, Phi, 
                                   ESigma, dim_rho, p.X, sum_d, hw_a, A_prior, nu_prior){
  
  ssr <- tXy - tXX %*% rho
  ridge <- - ridge %*% rho
  if (p.X > 0){nonfe_rho <- rho[-seq_len(p.X)]}else{nonfe_rho <- rho}
  
  Rmatrix <- mapply(split(nonfe_rho, rho_idx), dim_rho, SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d, ncol = d)})
  
  prior <- mapply(Rmatrix, nu, hw_a, A_prior, ESigma, nu_prior, SIMPLIFY = FALSE,
    FUN=function(R_j, nu_j, hw_a_j, A_j, ESigma.inv.j, nu_prior_j){
      inv_R_j <- solve(R_j)
      d_j <- ncol(R_j)
      diag_meat <- diag(t(inv_R_j) %*% ESigma.inv.j %*% inv_R_j)
      rho_hw_j <- diag_meat * nu_prior_j + 1/A_j^2

      term_profiled <- -(nu_j + 1)/2 * nu_prior_j * 1/rho_hw_j
      zeromat <- matrix(0, nrow = d_j, ncol = d_j)
      invRE <- t(inv_R_j) %*% ESigma.inv.j
      term_profiled <- mapply(seq_len(d_j), term_profiled, SIMPLIFY = FALSE, FUN=function(p, w){
        zeromat[,p] <- 2 * invRE[p,]
        as.vector(- t(inv_R_j) %*% zeromat %*% t(inv_R_j)) * w
      })
      term_profiled <- Reduce('+', term_profiled)
      out <- as.vector(- nu_j * t(inv_R_j) + term_profiled)
      return(out)
    })
  prior <- c(rep(0, p.X), unlist(prior))
  deriv_rho <- as.vector( ssr + ridge + prior )
  
  return(deriv_rho)
}