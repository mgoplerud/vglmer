update_rho <- function(XR, y, omega, prior_precision, 
                       moments_sigma_alpha,
                       prior_sigma_alpha_nu, prior_sigma_alpha_phi,
                       vi_a_a_jp, vi_a_b_jp, vi_a_nu_jp,
                       vi_a_APRIOR_jp, 
                       spline_REs, vi_beta_mean,
                       p.X, d_j, stationary_rho,
                       do_huangwand, offset,
                       offset_bilinear,
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
  if (offset != 0){ # Negative Binomial Offset
    tXY <- tXy + offset * matrix(colSums(omega %*% XR))
  }
  if (length(offset_bilinear) != 1){
    tXy <- tXy - t(t(offset_bilinear) %*% (omega %*% XR))
  }
  prior_precision <- prior_precision
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
      opt_rho <- tryCatch(optim(par = null_rho, fn = eval_profiled_rho, 
                       gr = eval_grad_profiled_rho,
                       method = 'L-BFGS-B', control = ctrl_opt,
                       tXy = tXy, tXX = tXX,
                       ridge = prior_precision, rho_idx = rho_idx,
                       nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X,
                       sum_d = sum_d, hw_a = hw_a, A_prior = A_prior,
                       nu_prior = nu_prior
      ), error = function(e){NULL})
      if (is.null(opt_rho)){
        message('optimization failed; trying with BFGS instead of L-BFGS-B')
        warning('optimization failed; trying with BFGS instead of L-BFGS-B')
        opt_rho <- optim(par = null_rho, fn = eval_profiled_rho, 
              gr = eval_grad_profiled_rho,
              method = 'BFGS', control = ctrl_opt,
              tXy = tXy, tXX = tXX,
              ridge = prior_precision, rho_idx = rho_idx,
              nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X,
              sum_d = sum_d, hw_a = hw_a, A_prior = A_prior,
              nu_prior = nu_prior
        )
      }
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
      improvement <- NA
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

eval_px_rotation_rho <- function(rho, SSQ_u, SSQ_v, dim_u, dim_v, nu,
                                 ESigma.inv, dim_rho, do_huangwand_mi,
                                 hw_a, A_prior, nu_prior, 
                                 iw_Phi, iw_nu, group_hier,
                                 mi_prior_type, px_sigma, simple_j,
                                 partial_fix){
  
  Rmatrix <- matrix(rho, nrow = dim_rho)
  lndet_Rmatrix <- as.numeric(determinant(Rmatrix)$modulus)
  
  inv_Rmatrix <- solve(Rmatrix)
  if (mi_prior_type %in% c('centered')){
    # u ~ N(0, R^T R) -> -1/2 tr(SSQ R^{-} R^{-T}) - 1/2 lndet(R) * 2
    prior_u <- -1/2 * sum(diag(SSQ_u %*% inv_Rmatrix %*% t(inv_Rmatrix)))
    prior_u <- prior_u + - dim_u * lndet_Rmatrix
  }else if (mi_prior_type %in% c('shared')){
    if (!px_sigma){
      prior_u <- -1/2 * sum(diag(SSQ_u %*% inv_Rmatrix %*% ESigma.inv %*% t(inv_Rmatrix)))
      prior_u <- prior_u + - dim_u * lndet_Rmatrix
    }else{
      # u ~ N(0, R2 Sigma R2)
      inv_R2 <- inv_Rmatrix %*% t(inv_Rmatrix)
      prior_u <- -1/2 * sum(diag(SSQ_u %*% inv_R2 %*% ESigma.inv %*% t(inv_R2)))
      prior_u <- prior_u + - dim_u * lndet_Rmatrix * 2
    }
  }else if (mi_prior_type %in% c('separate')){
    if (simple_j){
      # u ~ N(0, R^T Sigma_1 R) -> -1/2 tr(SSQ R^{-} Sigma_1^{-1} R^{-T})  - 1/2 lndet(R) * 2
      prior_u <- -1/2 * sum(diag(SSQ_u %*% inv_Rmatrix %*% ESigma.inv[[1]] %*% t(inv_Rmatrix)))
      prior_u <- prior_u + - dim_u * lndet_Rmatrix
    }else{
      prior_u <- sum(mapply(SSQ_u, dim_u, ESigma.inv[group_hier[[1]]], FUN=function(SSQ_l, dim_l, Einv_l){
        prior_u <- -1/2 * sum(diag(SSQ_l %*% inv_Rmatrix %*% Einv_l %*% t(inv_Rmatrix)))
        prior_u <- prior_u + - dim_l * lndet_Rmatrix
        return(prior_u)
      }))
    }
  }else{
    stop('...')
  }
  
  if (!px_sigma){
    
    prior_variance <- 0
    
    # v ~ N(0, R^- Sigma_v R^{-T}) -> -1/2 tr(SSQ R^{T} Sigma_v^{-1} R)  - 1/2 lndet(R^{-1})* 2
    if (mi_prior_type %in% c('centered', 'shared')){
      prior_v <- -1/2 * sum(diag(SSQ_v %*% t(Rmatrix) %*% ESigma.inv %*% Rmatrix))
      prior_v <- prior_v + dim_v * lndet_Rmatrix
    }else if (mi_prior_type %in% c('separate')){
      if (simple_j){
        prior_v <- -1/2 * sum(diag(SSQ_v %*% t(Rmatrix) %*% ESigma.inv[[2]] %*% Rmatrix))
        prior_v <- prior_v + dim_v * lndet_Rmatrix
      }else{
        prior_v <- sum(mapply(SSQ_v, dim_v, ESigma.inv[group_hier[[2]]], FUN=function(SSQ_l, dim_l, Einv_l){
          prior_v <- -1/2 * sum(diag(SSQ_l %*% t(Rmatrix) %*% Einv_l %*% Rmatrix))
          prior_v <- prior_v + dim_l * lndet_Rmatrix
          return(prior_v)
        }))
      }
    }else{
      stop('...')
    }
    
  }else{
    
    prior_v <- 0
    if (do_huangwand_mi){
      if (mi_prior_type %in% c('shared', 'centered')){
        
        diag_meat <- diag(t(inv_Rmatrix) %*% ESigma.inv %*% inv_Rmatrix)
        rho_hw <- nu_prior * diag_meat + 1/A_prior^2
        prior_variance <- -(nu + 1)/2 * sum(log(rho_hw)) +
          - nu * lndet_Rmatrix +
          sum(hw_a)
        
      }else if (mi_prior_type %in% c('separate')){
        
        if (simple_j){
          diag_meat <- diag(t(inv_Rmatrix) %*% ESigma.inv[[2]] %*% inv_Rmatrix)
          rho_hw <- nu_prior * diag_meat + 1/A_prior[[2]]^2
          prior_variance <- -(nu + 1)/2 * sum(log(rho_hw)) +
            - nu * lndet_Rmatrix +
            sum(hw_a[[2]])
        }else{
          prior_variance <- mapply(ESigma.inv[group_hier[[2]]], A_prior[group_hier[[2]]], 
                 hw_a[group_hier[[2]]], FUN=function(Einv_l, A_prior_l, hw_l){
            diag_meat <- diag(t(inv_Rmatrix) %*% Einv_l %*% inv_Rmatrix)
            rho_hw <- nu_prior * diag_meat + 1/A_prior_l^2
            prior_variance <- -(nu + 1)/2 * sum(log(rho_hw)) +
              - nu * lndet_Rmatrix +
              sum(hw_l)
            return(prior_variance)
          })
          prior_variance <- sum(prior_variance)
        }
        
      }else{stop("...")}
    }else{
      # p(Sigma): IW(prior_nu, R^{-T} prior_phi R^{-1})
      if (mi_prior_type %in% c('shared', 'centered')){
        prior_variance <- - iw_nu * lndet_Rmatrix - 1/2 * sum(diag(
          inv_Rmatrix %*% iw_Phi %*% t(inv_Rmatrix) %*% ESigma.inv
        ))
      }else if (mi_prior_type %in% c('separate')){
        if (simple_j){
          prior_variance <- - iw_nu * lndet_Rmatrix - 1/2 * sum(diag(
            inv_Rmatrix %*% iw_Phi[[2]] %*% t(inv_Rmatrix) %*% ESigma.inv[[2]]
          ))
        }else{
          prior_variance <- sum(mapply(ESigma.inv[group_hier[[2]]], iw_Phi[group_hier[[2]]],
             FUN=function(Einv_l, iw_Phi_l){
              prior_variance <- - iw_nu * lndet_Rmatrix - 1/2 * sum(diag(
                inv_Rmatrix %*% iw_Phi_l %*% t(inv_Rmatrix) %*% Einv_l
              ))
              return(prior_variance)       
          }))
        }
      }else{stop('invalid mi_prior_type')}
    }
  }
  return(as.numeric(prior_u + prior_v + prior_variance))
}


update_px_rotation <- function(vi_mi_SSQ,
                               vi_mi_moments,
                               vi_mi_size,
                               vi_mi_dim,
                               vi_mi_a_a_jp,
                               vi_mi_a_b_jp,
                               vi_mi_a_APRIOR_jp,
                               vi_mi_a_nu_jp,
                               Z_MI_grouping,
                               mi_prior_sigma_alpha_nu,
                               mi_prior_sigma_alpha_phi,
                               do_huangwand_mi,
                               mi_prior_type,
                               px_sigma, simple_j,
                               VEM_scalar, partial_fix
                               ){
  if (do_huangwand_mi){
    prior_weight <- vi_mi_a_nu_jp + vi_mi_dim - 1
    diag_weight <- mapply(vi_mi_a_a_jp, vi_mi_a_b_jp, vi_mi_a_nu_jp, SIMPLIFY = FALSE, 
      FUN = function(tilde.a, tilde.b, nu) {
        if (mi_prior_type %in% c('centered', 'shared')){
          Diagonal(x=tilde.a/tilde.b) * 2 * nu
        }else if (mi_prior_type %in% 'separate'){
          mapply(tilde.a, tilde.b, SIMPLIFY = FALSE, FUN=function(a_l, b_l){
            Diagonal(x = a_l/b_l) * 2 * nu
          })
        }else{
          stop('...')
        }
      })
  }else{
    diag_weight <- mi_prior_sigma_alpha_phi
    prior_weight <- mi_prior_sigma_alpha_nu
  }
  
  opt_all <- lapply(1:length(vi_mi_SSQ), FUN=function(j){
    
    Z_MI_grouping_j <- Z_MI_grouping[[j]]
    simple_j <- all(lengths(Z_MI_grouping_j) == 1)
    if (simple_j){
      dim_rho <- vi_mi_dim[j]
      dim_u <- vi_mi_size[[j]][[1]] * VEM_scalar
      dim_v <- vi_mi_size[[j]][[2]] * VEM_scalar 
      SSQ_u <- vi_mi_SSQ[[j]][[1]]
      SSQ_v <- vi_mi_SSQ[[j]][[2]]
      hw_a <- vi_mi_a_a_jp[[j]]
      A_prior <- vi_mi_a_APRIOR_jp[[j]]
      nu_prior <- vi_mi_a_nu_jp[[j]]
      nu <- prior_weight[[j]]
      Phi <- diag_weight[[j]]
      ESigma.inv <- vi_mi_moments[[j]]
      null_rho <- as.vector(Diagonal(n=ncol(SSQ_u)))
    }else{
      dim_rho <- vi_mi_dim[j]
      dim_u <- vi_mi_size[[j]][Z_MI_grouping_j[[1]]]
      dim_v <- vi_mi_size[[j]][Z_MI_grouping_j[[2]]]
      dim_u <- lapply(dim_u, FUN=function(i){i * VEM_scalar})
      dim_v <- lapply(dim_v, FUN=function(i){i * VEM_scalar})
      SSQ_u <- vi_mi_SSQ[[j]][Z_MI_grouping_j[[1]]]
      SSQ_v <- vi_mi_SSQ[[j]][Z_MI_grouping_j[[2]]]
      hw_a <- vi_mi_a_a_jp[[j]]
      A_prior <- vi_mi_a_APRIOR_jp[[j]]
      nu_prior <- vi_mi_a_nu_jp[[j]]
      nu <- prior_weight[[j]]
      Phi <- diag_weight[[j]]
      ESigma.inv <- vi_mi_moments[[j]]
      null_rho <- as.vector(Diagonal(n=dim_rho))
    }

    opt_rho <- optim(par = null_rho,
                     SSQ_u = SSQ_u, SSQ_v = SSQ_v,
                     dim_u = dim_u, dim_v = dim_v, nu = nu, 
                     ESigma.inv = ESigma.inv, dim_rho = dim_rho, 
                     hw_a = hw_a, A_prior = A_prior, nu_prior = nu_prior,
                     iw_Phi = Phi, iw_nu = nu,
                     f = eval_px_rotation_rho, do_huangwand_mi = do_huangwand_mi,
                     mi_prior_type = mi_prior_type, px_sigma = px_sigma,
                     simple_j = simple_j, group_hier = Z_MI_grouping_j,
                     partial_fix = partial_fix,
                     control = list(fnscale = -1),
                     method = 'BFGS')
    
    null_eval <- eval_px_rotation_rho(rho = null_rho,
      SSQ_u = SSQ_u, SSQ_v = SSQ_v, mi_prior_type = mi_prior_type,
      dim_u = dim_u, dim_v = dim_v, nu = nu,  do_huangwand_mi = do_huangwand_mi,
      ESigma.inv = ESigma.inv, dim_rho = dim_rho, group_hier = Z_MI_grouping_j,
      simple_j = simple_j, partial_fix = partial_fix,
      iw_Phi = Phi, iw_nu = nu, px_sigma = px_sigma,
      hw_a = hw_a, A_prior = A_prior, nu_prior = nu_prior
    )
    
    opt_Rmatrix <- matrix(opt_rho$par, nrow = dim_rho)
    
    inv_R <- solve(opt_Rmatrix)
    
    if (px_sigma){
      if (do_huangwand_mi){
        if (mi_prior_type %in% c('centered', 'shared')){
          diag_meat <- diag(t(inv_R) %*% ESigma.inv %*% inv_R)
          opt_rho_hw <- nu_prior * diag_meat + 1/A_prior^2
        }else{
          if (simple_j){
            diag_meat <- diag(t(inv_R) %*% ESigma.inv[[2]] %*% inv_R)
            opt_rho_hw <- list(NA, nu_prior * diag_meat + 1/A_prior[[2]]^2)
          }else{
            if (partial_fix){browser()}
            opt_rho_hw <- mapply(ESigma.inv[Z_MI_grouping_j[[2]]], A_prior[Z_MI_grouping_j[[2]]],
                   SIMPLIFY = FALSE,
              FUN=function(Einv_l, A_prior_l){
                diag_meat <- diag(t(inv_R) %*% Einv_l %*% inv_R)
                out <- nu_prior * diag_meat + 1/A_prior_l^2
                return(out)
            })
          }
        }
      }else{
        opt_rho_hw <- NULL
      }
    }else{
      opt_rho_hw <- NULL
    }
    
    return(list(
      R = opt_Rmatrix,
      rho_hw = opt_rho_hw,
      diff = opt_rho$value - null_eval
    ))
    
  })
  
  out <- list(
    R = lapply(opt_all, `[[`, 'R'),  
    rho_hw = lapply(opt_all, `[[`, 'rho_hw'), 
    diff = sapply(opt_all, `[[`, 'diff')
  )
  out$diff <- sum(out$diff)
  return(out)
}


fast_insert <- function(A, B, index){
  # B[index, index] <- A
  # can be *frightfully* expensive to do computationally
  A_dgT <- as(A, 'dgTMatrix')
  B_dgT <- as(B, 'dgTMatrix')
  
  zero_index <- index - 1
  # These are the positions of the non-A block of "B" that should be kept
  nonA_pos <- !( (B_dgT@i %in% zero_index) & (B_dgT@j %in% zero_index) )
  
  out <- sparseMatrix(
    # Add in triplet form the non-A block of "B" and the 
    # block of "A"
    i = c(B_dgT@i[nonA_pos], zero_index[A_dgT@i + 1]),
    j = c(B_dgT@j[nonA_pos], zero_index[A_dgT@j + 1]),
    x = c(B_dgT@x[nonA_pos], A_dgT@x),
    use.last.ij = FALSE,
    index1 = FALSE,
    repr = 'C',
    dims = dim(B)
  )  
  return(out)
  
  # Version 1: Reasonable but B[index, index] can be expensive  
  # # This is a matrix with zeros everywhere except "A"
  # aug_A <- sparseMatrix(i = index[A_dgT@i + 1],
  #                       j = index[A_dgT@j + 1],
  #                       x = A@x, dims = dim(B)
  # )
  # # This is matrix with zeros everywhere except for the old values of "B"
  # old_B <- as(B[index,index], 'dgTMatrix')
  # old_B <- sparseMatrix(i = index[old_B@i + 1],
  #                       j = index[old_B@j + 1],
  #                       x = old_B@x, dims = dim(B)
  # )
  # # This removes the old_B and adds the new A
  # out <- (B - old_B + aug_A)
  # return(out)
}
