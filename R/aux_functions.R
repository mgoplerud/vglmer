

update_rho <- function(XR, y, omega, prior_precision, 
                       vec_OSL_prior, 
                       moments_sigma_alpha,
                       diag_weight,
                       prior_weight,
                       spline_REs, vi_beta_mean,
                       p.X, d_j, stationary_rho,
                       px_it = NULL,
                       method = 'numerical'){
  
  
  ESigma <- lapply(moments_sigma_alpha[c(which(spline_REs), which(!spline_REs))], FUN=function(i){i$sigma.inv})
  Phi <- diag_weight[c(which(spline_REs), which(!spline_REs))] 
  nu <- prior_weight[c(which(spline_REs), which(!spline_REs))]
  
  y <- as.vector(y)
  
  if (method == 'numerical'){
    sum_ysq <- sum(omega %*% y)
    tXy <- t(XR) %*% y
    tXX <- t(XR) %*% omega %*% XR
    
    rho_idx <- d_j[spline_REs]
    rho_idx <- rho_idx * seq_len(length(rho_idx))
    rho_idx <- c(rho_idx, rep(seq_len(sum(!spline_REs)), times = d_j[!spline_REs]^2) + sum(spline_REs))
    
    null_rho <- c(rep(1, sum(spline_REs)), stationary_rho)
    null_rho <- c(as.vector(vi_beta_mean), null_rho)
    dim_rho <- c(rep(1, sum(spline_REs)), d_j[!spline_REs])

    ctrl_opt <- list(fnscale = -1)
    if (!is.null(px_it)){
      ctrl_opt$maxit <- px_it
    }
    
    opt_rho <- optim(par = null_rho, fn = eval_rho, 
        gr = eval_grad_rho,
        method = 'L-BFGS-B', control = ctrl_opt,
        tXy = tXy, tXX = tXX,
        ridge = prior_precision, rho_idx = rho_idx,
        nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X
    )
    null_eval <- eval_rho(null_rho, tXy = tXy, tXX = tXX, ridge = prior_precision, rho_idx = rho_idx,
             nu = nu, Phi = Phi, ESigma = ESigma, dim_rho = dim_rho, p.X = p.X)
    if (opt_rho$value < null_eval){
      message('OPTIMIZATION FAILED')
      stop()
    }
    opt_rho <- opt_rho$par
    names(opt_rho) <- NULL
  }else if (method == 'OSL'){
    opt_rho <- vecR_fast_ridge(X = XR, 
     omega = omega, prior_precision = prior_precision, y = y, 
     adjust_y = as.vector(vec_OSL_prior)) 
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


