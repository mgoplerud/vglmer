
#' Create a sparse design matrix for the fixed effects
vglmer_build_fe <- function(x, by, contrast_type, levels = NULL){
  if (is.null(levels)){
    ux <- unique(x)
  }else{
    ux <- levels
  }
  if (is.null(by)){
    x <- sparseMatrix(i = 1:length(x), j = match(x, ux), x = 1, dims = c(length(x), length(ux)))
  }else{
    x <- sparseMatrix(i = 1:length(x), j = match(x, ux), x = by, dims = c(length(x), length(ux)))
  }
  colnames(x) <- ux
  out <- list(x = x, attr = list(unique_values = ux, constrast = contrast_type(ncol(x), sparse = TRUE)))
  class(out) <- c('fe_sparse')
  return(list(out))
}

obj_poisson <- function(par, y, Textend, P, joint.XZ, diag_position, exp_diag = TRUE){
  joint_mean <- matrix(par[1:ncol(joint.XZ)])
  joint_L <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(joint.XZ), 2))
  
  lt_par <- par[-1:-ncol(joint.XZ)]
  if (exp_diag){
    lt_par[diag_position] <- exp(lt_par[diag_position])
  }
  joint_L[lower.tri(joint_L, diag = TRUE)] <- lt_par
  if (any(diag(joint_L) < 0)){stop('Invalid lower triangular matrix')}
  joint_LP <- joint_L %*% t(P)
  
  mean_lp <- as.vector(joint.XZ %*% joint_mean)
  var_lp <- rowSums((joint.XZ %*% t(joint_LP))^2)  
  poisson_weight <- exp(mean_lp + 1/2 * var_lp)
  
  ll_1 <- sum(y * mean_lp)
  ll_2 <- sum(-poisson_weight)
  
  lp_1 <- -1/2 * sum( (Textend %*% joint_mean) * joint_mean )
  lp_2 <- -1/2 * sum(diag((Textend %*% t(joint_LP)) %*% joint_LP))
  
  lndiag_L <- log(diag(joint_L)^2)
  lndet <- 1/2 * sum(lndiag_L)
  
  out <- ll_1 + ll_2 + lp_1 + lp_2 + lndet  
  return(out)
}

grad_poisson <- function(par, y, Textend, P, joint.XZ, diag_position, exp_diag = TRUE){
 
  joint_mean <- matrix(par[1:ncol(joint.XZ)])
  
  joint_L <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(joint.XZ), 2))
  lt_par <- par[-1:-ncol(joint.XZ)]
  if (exp_diag){
    lt_par[diag_position] <- exp(lt_par[diag_position])
  }
  joint_L[lower.tri(joint_L, diag = TRUE)] <- lt_par
  joint_LP <- joint_L %*% t(P)
  
  mean_lp <- as.vector(joint.XZ %*% joint_mean)
  var_lp <- rowSums((joint.XZ %*% t(joint_LP))^2)  
  poisson_weight <- exp(mean_lp + 1/2 * var_lp)
  
  grad_mean <- as.vector(t(joint.XZ) %*% y - t(joint.XZ) %*% poisson_weight +
    - Textend %*% joint_mean)
  
  grad_meat <- t(joint.XZ) %*% Diagonal(x = poisson_weight) %*% joint.XZ + Textend
  
  grad_variance <- matrixcalc::vech(
    as.matrix(
      Diagonal(x = 1/abs(diag(joint_L))) +
      - joint_L %*% t(P) %*% grad_meat %*% P
    )
  )
  if (exp_diag){
    grad_variance[diag_position] <- grad_variance[diag_position] * lt_par[diag_position]
  }
  
  return(c(grad_mean, grad_variance))
}

gradient_ascent_poisson <- function(par, diag_position, beta = 2, 
    outer_iter = 15, inner_iter = 20, y, Textend, vi_joint_LP, joint.XZ){
  
  obj_original <- obj_poisson(par = par, y = y, Textend = Textend, diag_position = diag_position,
                              P = vi_joint_LP, joint.XZ = joint.XZ)
  print('Starting')
  print(obj_original)
  out_par <- par
  out_obj <- obj_original
  store_grad <- rep(0, length(par))
  for (oi in seq_len(outer_iter)){
    cat('.')
    direction_p <- grad_poisson(par = out_par, y = y, Textend = Textend, diag_position = diag_position,
      P = vi_joint_LP, joint.XZ = joint.XZ)
    store_grad <- store_grad + direction_p^2
    print(sqrt(sum(direction_p^2)))
    prop_step <- Diagonal(x = 1/sqrt(1e-6 + store_grad)) %*% direction_p

    g <- function(alpha){obj_poisson(par = out_par + exp(alpha) * prop_step, y = y, diag_position = diag_position,
        Textend = Textend, P = vi_joint_LP, joint.XZ = joint.XZ)}
    opt_g <- optimize(f = g, interval = c(-30, 0), maximum = TRUE)
    print(opt_g)
    out_par <- out_par + exp(opt_g$maximum) * prop_step
    plot(prop_step)
    # for (iter in 0:inner_iter){
    #   
    #   prop_par <- out_par + beta^(-iter) * direction_p
    #   
    #   prop_obj <- obj_poisson(par = prop_par, y = y, diag_position = diag_position,
    #                           Textend = Textend, P = vi_joint_LP, joint.XZ = joint.XZ)  
    #   converge_inner <- prop_obj > obj_original
    #   if (converge_inner){
    #     out_par <- prop_par
    #     out_obj <- prop_obj
    #     succeed <- 1
    #     break
    #   }
    # }
    
  }
  print("output")
  print(c(iter, succeed))
  print(out_obj)
  if (out_obj < obj_original){stop('...')}
  return(out_par)
}

weight_poisson <- function(par, y, Textend, P, joint.XZ, diag_position, exp_diag = TRUE){
  joint_mean <- matrix(par[1:ncol(joint.XZ)])
  joint_L <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(joint.XZ), 2))
  
  lt_par <- par[-1:-ncol(joint.XZ)]
  if (exp_diag){
    lt_par[diag_position] <- exp(lt_par[diag_position])
  }
  joint_L[lower.tri(joint_L, diag = TRUE)] <- lt_par
  if (any(diag(joint_L) < 0)){stop('Invalid lower triangular matrix')}
  joint_LP <- joint_L %*% t(P)
  
  mean_lp <- as.vector(joint.XZ %*% joint_mean)
  var_lp <- rowSums((joint.XZ %*% t(joint_LP))^2)  
  poisson_weight <- exp(mean_lp + 1/2 * var_lp)
  return(poisson_weight)
}

poisson_NVMP <- function(y, joint.XZ, Textend, diag_position,
                         vi_pg_mean, diag_vi_pg_mean, old_param){

  chol.update.joint <- LinRegChol(X = joint.XZ,
             omega = diag_vi_pg_mean,
             prior_precision = Textend,
             y = y - vi_pg_mean, adj_y = -as.vector(Textend %*% old_param))
  Pmatrix <- sparseMatrix(i = 1:ncol(joint.XZ),
                          j = 1 + chol.update.joint$Pindex, x = 1)
  vi_joint_L_nonpermute <- drop0(solve(chol.update.joint$origL))
  chol.update.joint$mean <- old_param + chol.update.joint$mean
  vi_joint_LP <- Pmatrix
  log_det_joint_var <- -2 * sum(log(diag(chol.update.joint$origL)))

  vi_mean <- as.vector(chol.update.joint$mean)
  nvmp_par <- c(vi_mean, ks::vech(as.matrix(vi_joint_L_nonpermute)))
  nvmp_par[ncol(joint.XZ) + diag_position] <- log(nvmp_par[ncol(joint.XZ) + diag_position])
  
  return(list(par = nvmp_par, mean = vi_mean, P = vi_joint_LP, L = vi_joint_L_nonpermute))
}

damp_poisson_NVMP <- function(attempt_NVMP, y, old_param, old_weights,
                              Textend, joint.XZ, diag_position, starting_obj,
                              direct_optimize = TRUE){

  new_pois <- weight_poisson(par = attempt_NVMP$par, y = y, 
   Textend = Textend, P = attempt_NVMP$P, joint.XZ = joint.XZ,
   diag_position = diag_position)

  if (direct_optimize){
    
    f <- function(alpha_NVMP){
      weight_alpha <- new_pois * alpha_NVMP + (1 - alpha_NVMP) * old_weights
      weight_mu <- attempt_NVMP$par[1:ncol(joint.XZ)] * alpha_NVMP + (1 - alpha_NVMP) * old_param
      
      damp_chol.update.joint <- LinRegChol(X = joint.XZ, 
         omega = sparseMatrix(i = 1:nrow(joint.XZ), j = 1:nrow(joint.XZ), x = weight_alpha), 
         prior_precision = Textend,
         y = y - weight_alpha, adj_y = -as.vector(Textend %*% weight_mu))
      damp_chol.update.joint$mean <- damp_chol.update.joint$mean + weight_mu
      
      damp_Pmatrix <- sparseMatrix(i = 1:ncol(joint.XZ), j = 1 + damp_chol.update.joint$Pindex, x = 1)
      damp_Lmatrix <- drop0(solve(damp_chol.update.joint$origL))
      damp_par <- c(as.vector(damp_chol.update.joint$mean), ks::vech(as.matrix(damp_Lmatrix)))
      damp_par[ncol(joint.XZ) + diag_position] <- log(damp_par[ncol(joint.XZ) + diag_position])
      
      damp_obj <- obj_poisson(par = damp_par, 
            y = y, Textend = Textend,
            diag_position = diag_position,
            P = damp_Pmatrix, joint.XZ = joint.XZ)
      
      return(damp_obj)
    }
    
    grid_alpha <- optimize(f = f, interval = c(0,1), maximum = TRUE)$maximum
   
  } else {
    grid_alpha <- seq(0.5, 1, length.out = 11)
    
  }
  
  counter_NVMP <- 0
  while (counter_NVMP <= 10){
    
    alpha_NVMP <- grid_alpha[counter_NVMP + 1]
    
    weight_alpha <- new_pois * alpha_NVMP + (1 - alpha_NVMP) * old_weights
    weight_mu <- attempt_NVMP$par[1:ncol(joint.XZ)] * alpha_NVMP + (1 - alpha_NVMP) * old_param
    
    damp_chol.update.joint <- LinRegChol(X = joint.XZ, 
         omega = sparseMatrix(i = 1:nrow(joint.XZ), j = 1:nrow(joint.XZ), x = weight_alpha), 
         prior_precision = Textend,
         y = y - weight_alpha, adj_y = -as.vector(Textend %*% weight_mu))
    damp_chol.update.joint$mean <- damp_chol.update.joint$mean + weight_mu
    
    damp_Pmatrix <- sparseMatrix(i = 1:ncol(joint.XZ), j = 1 + damp_chol.update.joint$Pindex, x = 1)
    damp_Lmatrix <- solve(damp_chol.update.joint$origL)
    damp_par <- c(as.vector(damp_chol.update.joint$mean), ks::vech(as.matrix(damp_Lmatrix)))
    damp_par[ncol(joint.XZ) + diag_position] <- log(damp_par[ncol(joint.XZ) + diag_position])
    
    damp_obj <- obj_poisson(par = damp_par, 
                            y = y, Textend = Textend,
                            diag_position = diag_position,
                            P = damp_Pmatrix, joint.XZ = joint.XZ)
    
    if (damp_obj > starting_obj){
      break
    }else{
      counter_NVMP <- counter_NVMP + 1
    }
  }
  if (damp_obj > starting_obj){
    return(list(par = damp_par, obj = damp_obj, P = damp_Pmatrix))
  }else{
    return(NULL)
  }
}

poisson_fixedpoint <- function(par, y, Textend, P, joint.XZ, diag_position,
                               OUTER_ITER = 1, INNER_ITER = 5){
  
  init_obj <- obj_poisson(par = par, y = y, Textend = Textend, P = P, joint.XZ = joint.XZ, diag_position = diag_position)
  init_par <- par
  joint_mean <- matrix(par[1:ncol(joint.XZ)])
  joint_L <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(joint.XZ), 2))
  
  lt_par <- par[-1:-ncol(joint.XZ)]
  lt_par[diag_position] <- exp(lt_par[diag_position])
  joint_L[lower.tri(joint_L, diag = TRUE)] <- lt_par
  if (any(diag(joint_L) < 0)){stop('Invalid lower triangular matrix')}
  joint_L_decomp <- joint_L %*% t(P)
  
  mean_lp <- as.vector(joint.XZ %*% joint_mean)
  var_lp <- rowSums((joint.XZ %*% t(joint_L_decomp))^2)  
  poisson_weight <- exp(mean_lp + 1/2 * var_lp)
  
  for (oi in 1:OUTER_ITER){
    
    for (inner_mean in 1:INNER_ITER){
      # Update the mean parameters
      newton_raphson <- solve(
        Cholesky(Textend + t(joint.XZ) %*% Diagonal(x = poisson_weight) %*% joint.XZ), 
        t(joint.XZ) %*% (y - poisson_weight) - Textend %*% joint_mean)
      
      joint_mean <- par[1:ncol(joint.XZ)] + newton_raphson
      par[1:ncol(joint.XZ)] <- joint_mean
      mean_lp <- as.vector(joint.XZ %*% joint_mean)
      poisson_weight <- exp(mean_lp + 1/2 * var_lp)
    }
    
    # print('Mean')
    # print(obj_poisson(par = par, y = y, Textend = Textend, P = P, joint.XZ = joint.XZ, diag_position = diag_position))
    
    # Update the variance parameters
    meat_grad <- t(joint.XZ) %*% Diagonal(x = poisson_weight) %*% joint.XZ + Textend
    grad_vechL <- matrixcalc::vech(
      as.matrix(
        Diagonal(x = 1/abs(diag(joint_L))) +
          - joint_L %*% t(P) %*% meat_grad %*% P
      )
    )
    grad_vechL[diag_position] <- grad_vechL[diag_position] * diag(joint_L)
    current_lt <- matrixcalc::vech(as.matrix(joint_L))
    current_lt[diag_position] <- log(current_lt[diag_position])
    
    opt_lt <- optim(par = current_lt, 
      fn = function(x, ...){obj_poisson(par = c(joint_mean[,], x), ...)},
      gr = function(x, ...){grad_poisson(par = c(joint_mean[,], x), ...)[-1:-ncol(joint.XZ)]},
      method = 'BFGS', control = list(fnscale = -1, maxit = 5),
      P = P, diag_position = diag_position, joint.XZ = joint.XZ, Textend = Textend, y = y
    )
    par[-1:-ncol(joint.XZ)] <- opt_lt$par
    
    joint_L[lower.tri(joint_L, diag = T)] <- par[-1:-ncol(joint.XZ)]
    diag(joint_L) <- exp(diag(joint_L))
    # Update the weights
    joint_L_decomp <- joint_L %*% t(P)
    var_lp <- rowSums((joint.XZ %*% t(joint_L_decomp))^2)  
    poisson_weight <- exp(mean_lp + 1/2 * var_lp)
    # print('Final')
    # print(obj_poisson(par = par, y = y, Textend = Textend, P = P, joint.XZ = joint.XZ, diag_position = diag_position))
    
  }
  return(par)
}

