

fast_contr_sum <- function(N){
  N_1 <- N - 1; grid_N1 <- seq_len(N_1) 
  return(sparseMatrix(i = c(grid_N1,rep(N, N_1)), 
      j = c(grid_N1, grid_N1), x = c(rep(1,N_1), rep(-1,N_1))))
}

obj_poisson <- function(par, y, Textend, P, joint.XZ, 
  diag_position, mask_lowertri, offset_weight, exp_diag = TRUE){
  
  joint_mean <- matrix(par[1:ncol(joint.XZ)])
  joint_L <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(joint.XZ), 2))
  
  lt_par <- par[-1:-ncol(joint.XZ)]
  if (exp_diag){
    lt_par[diag_position] <- exp(lt_par[diag_position])
  }
  joint_L[mask_lowertri] <- lt_par
  if (any(diag(joint_L) < 0)){browser(); stop('Invalid lower triangular matrix')}
  joint_LP <- joint_L %*% t(P)
  
  mean_lp <- as.vector(joint.XZ %*% joint_mean)
  var_lp <- rowSums((joint.XZ %*% t(joint_LP))^2)  
  poisson_weight <- exp(mean_lp + 1/2 * var_lp + offset_weight)
  
  ll_1 <- sum(y * mean_lp)
  ll_2 <- sum(-poisson_weight)
  
  lp_1 <- -1/2 * sum( (Textend %*% joint_mean) * joint_mean )
  lp_2 <- -1/2 * sum(diag((Textend %*% t(joint_LP)) %*% joint_LP))
  
  lndiag_L <- log(diag(joint_L)^2)
  lndet <- 1/2 * sum(lndiag_L)
  
  out <- ll_1 + ll_2 + lp_1 + lp_2 + lndet  
  return(out)
}

grad_poisson <- function(par, y, Textend, P, joint.XZ,
  diag_position, mask_lowertri, offset_weight, exp_diag = TRUE){
 
  joint_mean <- matrix(par[1:ncol(joint.XZ)])
  
  joint_L <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(joint.XZ), 2))
  lt_par <- par[-1:-ncol(joint.XZ)]
  if (exp_diag){
    lt_par[diag_position] <- exp(lt_par[diag_position])
  }
  joint_L[mask_lowertri] <- lt_par
  joint_LP <- joint_L %*% t(P)
  
  mean_lp <- as.vector(joint.XZ %*% joint_mean)
  var_lp <- rowSums((joint.XZ %*% t(joint_LP))^2)  
  poisson_weight <- exp(mean_lp + 1/2 * var_lp + offset_weight)
  
  grad_mean <- as.vector(t(joint.XZ) %*% y - t(joint.XZ) %*% poisson_weight +
    - Textend %*% joint_mean)
  
  grad_meat <- t(joint.XZ) %*% Diagonal(x = poisson_weight) %*% joint.XZ + Textend
  
  grad_variance <- Diagonal(x = 1/abs(diag(joint_L))) - joint_L %*% t(P) %*% grad_meat %*% P
  grad_variance <- grad_variance[mask_lowertri]
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

weight_poisson <- function(par, y, Textend, P, joint.XZ, mask_lowertri,
                           diag_position, offset_weight, exp_diag = TRUE){
  joint_mean <- matrix(par[1:ncol(joint.XZ)])
  joint_L <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(joint.XZ), 2))
  
  lt_par <- par[-1:-ncol(joint.XZ)]
  if (exp_diag){
    lt_par[diag_position] <- exp(lt_par[diag_position])
  }
  joint_L[mask_lowertri] <- lt_par
  if (any(diag(joint_L) < 0)){stop('Invalid lower triangular matrix')}
  joint_LP <- joint_L %*% t(P)
  
  mean_lp <- as.vector(joint.XZ %*% joint_mean)
  var_lp <- rowSums((joint.XZ %*% t(joint_LP))^2)  
  poisson_weight <- exp(mean_lp + 1/2 * var_lp + offset_weight)
  return(poisson_weight)
}

poisson_NVMP <- function(y, joint.XZ, Textend, diag_position,
                         offset_weight,
                         vi_pg_mean, diag_vi_pg_mean, old_param,
                         mask_lowertri){

  chol.update.joint <- LinRegChol(X = joint.XZ,
             omega = diag_vi_pg_mean,
             prior_precision = Textend,
             y = y - vi_pg_mean, adj_y = -as.vector(Textend %*% old_param))
  Pmatrix <- sparseMatrix(i = 1:ncol(joint.XZ),
                          j = 1 + chol.update.joint$Pindex, x = 1)
  vi_joint_L_nonpermute <- drop0(solve(as(chol.update.joint$origL, 'dtCMatrix')))
  chol.update.joint$mean <- old_param + chol.update.joint$mean
  vi_joint_LP <- Pmatrix
  log_det_joint_var <- -2 * sum(log(diag(chol.update.joint$origL)))

  vi_mean <- as.vector(chol.update.joint$mean)
  nvmp_par <- c(vi_mean, vi_joint_L_nonpermute[mask_lowertri])
  nvmp_par[ncol(joint.XZ) + diag_position] <- log(nvmp_par[ncol(joint.XZ) + diag_position])
  
  return(list(par = nvmp_par, mean = vi_mean, P = vi_joint_LP, L = vi_joint_L_nonpermute))
}


damp_poisson_NVMP <- function(attempt_NVMP, y, old_param, old_weights,
                              Textend, joint.XZ, diag_position, starting_obj,
                              mask_lowertri, offset_weight,
                              direct_optimize = TRUE){

  new_pois <- weight_poisson(par = attempt_NVMP$par, y = y, 
                             offset_weight = offset_weight,
   Textend = Textend, P = attempt_NVMP$P, joint.XZ = joint.XZ,
   diag_position = diag_position, mask_lowertri = mask_lowertri)

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
      damp_par <- c(as.vector(damp_chol.update.joint$mean), damp_Lmatrix[mask_lowertri])
      damp_par[ncol(joint.XZ) + diag_position] <- log(damp_par[ncol(joint.XZ) + diag_position])
      
      damp_obj <- obj_poisson(par = damp_par, 
            y = y, Textend = Textend, offset_weight = offset_weight,
            diag_position = diag_position, mask_lowertri = mask_lowertri,
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
    if (is.na(alpha_NVMP)){break}
    
    weight_alpha <- new_pois * alpha_NVMP + (1 - alpha_NVMP) * old_weights
    weight_mu <- attempt_NVMP$par[1:ncol(joint.XZ)] * alpha_NVMP + (1 - alpha_NVMP) * old_param
    
    damp_chol.update.joint <- LinRegChol(X = joint.XZ, 
         omega = sparseMatrix(i = 1:nrow(joint.XZ), j = 1:nrow(joint.XZ), x = weight_alpha), 
         prior_precision = Textend,
         y = y - weight_alpha, adj_y = -as.vector(Textend %*% weight_mu))
    damp_chol.update.joint$mean <- damp_chol.update.joint$mean + weight_mu
    
    damp_Pmatrix <- sparseMatrix(i = 1:ncol(joint.XZ), j = 1 + damp_chol.update.joint$Pindex, x = 1)
    damp_Lmatrix <- solve(damp_chol.update.joint$origL)
    damp_par <- c(as.vector(damp_chol.update.joint$mean), damp_Lmatrix[mask_lowertri])
    damp_par[ncol(joint.XZ) + diag_position] <- log(damp_par[ncol(joint.XZ) + diag_position])
    
    damp_obj <- obj_poisson(par = damp_par, offset_weight = offset_weight,
                            y = y, Textend = Textend,
                            diag_position = diag_position, 
                            mask_lowertri = mask_lowertri,
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

update_poisson <- function(current_param, current_lowertri, 
      y, Textend, joint.XZ, diag_position, existing_P, diag_vi_pg_mean,
      vi_pg_mean, old_vi_pg_mean, old_Textend, old_param, mask_lowertri,
      starting_obj, it, force_ascent, offset_weight, allow_decrease = FALSE){

  attempt_NVMP <- poisson_NVMP(y = y, joint.XZ = joint.XZ, Textend = Textend, 
     vi_pg_mean = vi_pg_mean, diag_position = diag_position,
     mask_lowertri = mask_lowertri, offset_weight = offset_weight,
     diag_vi_pg_mean = diag_vi_pg_mean, old_param = current_param)
  
  attempt_NVMP_obj <- obj_poisson(par = attempt_NVMP$par, y = y,
    P = attempt_NVMP$P, Textend = Textend, joint.XZ = joint.XZ,
    offset_weight = offset_weight,
    diag_position = diag_position, mask_lowertri = mask_lowertri)
  
  if (allow_decrease | (abs(attempt_NVMP_obj) < 1e12 & (attempt_NVMP_obj > starting_obj) & (!force_ascent | it == 1)) ){
    vi_var_L_nonpermute <- attempt_NVMP$L 
    vi_var_P <- attempt_NVMP$P 
    vi_update_mean <- attempt_NVMP$mean
    type <- 'NVMP'
    out_obj <- attempt_NVMP_obj
  }else{
    
    if (!force_ascent & (it != 1) & (abs(attempt_NVMP_obj) < 1e12)){
      # Try dampening to address convergence issue
      damp_par <- tryCatch(damp_poisson_NVMP(
         attempt_NVMP = attempt_NVMP, y = y, old_param = old_param, 
         old_weights = old_vi_pg_mean, Textend = Textend, 
         offset_weight = offset_weight,
         joint.XZ = joint.XZ, diag_position = diag_position, 
         starting_obj = starting_obj, mask_lowertri = mask_lowertri,
         direct_optimize = TRUE), error = function(e){NULL})
    }else{
      damp_par <- NULL
    }
    
    # If damping fails or makes it worse, or force ascent, use direct
    # optimization...
    
    if (is.null(damp_par) | force_ascent){
      opt_obj <- tryCatch(optim(par = c(current_param, current_lowertri), 
        fn = obj_poisson, gr = grad_poisson,
        control = list(fnscale = -1, maxit = 30),
        offset_weight = offset_weight,
        method = 'CG', y = y, mask_lowertri = mask_lowertri,
        P = existing_P, Textend = Textend, joint.XZ = joint.XZ, 
        diag_position = diag_position), error = function(e){NULL})
      
      opt_failed <- FALSE
      if (!is.null(opt_obj)){
        if (opt_obj$value < starting_obj){
          opt_failed <- TRUE
        }
      }else{
        opt_failed <- TRUE
      }
      if (opt_failed){
        warning(paste0('Direct optimization (and NVMP and damping) failed;',
                       ' not updating q(beta,alpha) parameters at', it, 
                       '; try rescaling parameters or using calibrate_init=TRUE'))
        ascent_par <- c(current_param, current_lowertri)
        type <- "FAILED"
        browser()
        out_obj <- starting_obj
      }else{
        ascent_par <- opt_obj$par
        out_obj <- opt_obj$value
        type <- "OPTIMIZE"
      }
      vi_var_P <- existing_P
    }else{
      type <- 'DAMP'
      vi_var_P <- damp_par$P
      ascent_par <- damp_par$par
      out_obj <- damp_par$obj
    }
    
    vi_update_mean <- ascent_par[1:ncol(joint.XZ)]
    vi_var_L_nonpermute <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(joint.XZ), 2))
    ascent_lt <- ascent_par[-1:-ncol(joint.XZ)]
    ascent_lt[diag_position] <- exp(ascent_lt[diag_position])
    vi_var_L_nonpermute[mask_lowertri] <- ascent_lt
    
  }
  return(
    list(L = vi_var_L_nonpermute, P = vi_var_P, 
         mean = vi_update_mean, obj = out_obj,
         type = type)
  )
}

cyclical_update_poisson <- function(X, Z, 
    y, vi_pg_mean, diag_vi_pg_mean, cyclical_pos, 
    vi_alpha_mean, vi_beta_mean, zeromat, it, Tinv, seq_N,
    old_Tinv, old_param, old_vi_pg_mean, mask_lowertri, force_ascent,
    name_RE, number_of_RE, vi_beta_L_nonpermute, vi_alpha_L_nonpermute,
    vi_beta_LP, vi_alpha_LP, vi_alpha_decomp, vi_beta_decomp,
    allow_decrease = FALSE, adj_fe = NULL){
  
  if (it > 1){
    if (ncol(X) > 0){
      old_param_alpha <- old_param[-seq_len(ncol(X))]
    }else{
      old_param_alpha <- old_param
    }
    old_param_beta <- old_param[seq_len(ncol(X))]
  }
  log_weight <- log(vi_pg_mean)
  log_det_beta_var <- NA
  running_log_det_alpha_var <- rep(NA, number_of_RE)
  # vi_alpha_decomp <- as(vi_alpha_L_nonpermute %*% t(vi_alpha_LP), 'dgCMatrix')
  # vi_beta_decomp <- vi_beta_L_nonpermute %*% t(vi_beta_LP)

  if (number_of_RE == 0){
    loop_terms <- 0
  }else{
    loop_terms <- c(0, sample(number_of_RE))
  }
  if (ncol(X) == 0){
    loop_terms <- setdiff(loop_terms, 0)
  }
  
  for (j in loop_terms) {
    # print(j)
    if (j != 0){
      name_j <- name_RE[j]
      mask_j <- mask_lowertri[[name_j]]
      index_j <- cyclical_pos[[j]]
      Z_j <- Z[, index_j, drop = F]
      current_param <- vi_alpha_mean[index_j]
      current_Lnonperm_var <- vi_alpha_L_nonpermute[index_j, index_j]
      current_lowertri <- current_Lnonperm_var[mask_j$mask]
      current_P <- vi_alpha_LP[index_j, index_j]
      Tinv_j <- Tinv[[j]]
      if (it > 1){
        old_param_j <- old_param_alpha[index_j]
        old_Tinv_j <- old_Tinv[[j]]
      }
    }else{
      Z_j <- X
      current_param <- vi_beta_mean[,]
      mask_j <- mask_lowertri[["..fe.."]]
      current_Lnonperm_var <- vi_beta_L_nonpermute
      current_lowertri <- current_Lnonperm_var[mask_j$mask]
      Tinv_j <- zeromat
      if (it  > 1){
        old_param_j <- old_param_beta
        old_Tinv_j <- zeromat
      }
      current_P <- vi_beta_LP
    }

    current_lowertri[mask_j$diag] <- log(current_lowertri[mask_j$diag])

    current_weight_j <- as.vector(Z_j %*% current_param + 
      1/2 * rowSums( (Z_j %*% t(current_Lnonperm_var %*% t(current_P)))^2 ))
    offset_weight_j <- log_weight - current_weight_j
    
    starting_obj_j <- obj_poisson(par = c(current_param, current_lowertri),
      y = y, Textend = Tinv_j, P = current_P, offset_weight = offset_weight_j,
      joint.XZ = Z_j, diag_position = mask_j$diag, mask_lowertri = mask_j$mask)

    est_poisson <- update_poisson(current_param = current_param, 
        current_lowertri = current_lowertri, 
        y = y, Textend = Tinv_j, joint.XZ = Z_j, diag_position = mask_j$diag,
        existing_P = current_P, diag_vi_pg_mean = diag_vi_pg_mean, 
        vi_pg_mean = vi_pg_mean, 
        old_vi_pg_mean = old_vi_pg_mean, old_Textend = old_Tinv_j, 
        old_param = old_param_j, mask_lowertri = mask_j$mask, 
        starting_obj = starting_obj_j, it = it, force_ascent = force_ascent,
        allow_decrease = allow_decrease, offset_weight = offset_weight_j)
    
    update_decomp_j <- est_poisson$L %*% t(est_poisson$P)
    
    
    if (j != 0){
      vi_alpha_mean[index_j] <- est_poisson$mean
      vi_alpha_L_nonpermute[index_j, index_j] <- est_poisson$L
      vi_alpha_LP[index_j, index_j] <- est_poisson$P
      vi_alpha_decomp[index_j, index_j] <- update_decomp_j
      running_log_det_alpha_var[j] <- 2 * sum(log(diag(est_poisson$L)))
    }else{
      vi_beta_mean <- Matrix(est_poisson$mean)
      vi_beta_decomp <- update_decomp_j
      vi_beta_L_nonpermute <- est_poisson$L
      vi_beta_LP <- est_poisson$P
      log_det_beta_var <- 2 * sum(log(diag(vi_beta_L_nonpermute)))
    }
    

    diff_weight_j <- as.vector(Z_j %*% est_poisson$mean + 1/2 * rowSums( (Z_j %*% t(update_decomp_j))^2 ))
    log_weight <- offset_weight_j + diff_weight_j
    vi_pg_mean <- exp(log_weight)

    direct_update <- as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean + 1/2 * rowSums( (X %*% t(vi_beta_decomp))^2 ) +
                           1/2 * rowSums( (Z %*% t(vi_alpha_decomp))^2))
    direct_update <- exp(direct_update + adj_fe)
    vi_pg_mean_d <- as.vector(direct_update)
    
    if (!isTRUE(all.equal(as.vector(vi_pg_mean_d), as.vector(vi_pg_mean)))){
      browser('Misalignment')
    }
    
    diag_vi_pg_mean <- sparseMatrix(i = seq_N, j = seq_N, x = vi_pg_mean)
    if (est_poisson$type != 'NVMP'){
      message(est_poisson$type)
    }
    # print(c(j, est_poisson$type))
    # print(est_poisson$obj - starting_obj_j)  

  }
  
  running_log_det_alpha_var <- running_log_det_alpha_var[is.finite(running_log_det_alpha_var)]
  
  return(list(
    vi_pg_mean = vi_pg_mean,
    diag_vi_pg_mean = diag_vi_pg_mean,
    log_det_alpha_var = sum(running_log_det_alpha_var),
    log_det_beta_var = log_det_beta_var,
    vi_alpha_mean = vi_alpha_mean,
    vi_alpha_L_nonpermute = vi_alpha_L_nonpermute,
    vi_alpha_LP = vi_alpha_LP,
    vi_beta_mean = vi_beta_mean,
    vi_beta_L_nonpermute = vi_beta_L_nonpermute,
    vi_beta_LP = vi_beta_LP,
    vi_beta_decomp = vi_beta_decomp,
    vi_alpha_decomp = vi_alpha_decomp
  ))
  
}


update_poisson_FE <- function(y,
  FE_data, FE_lookup, FE_rowtens,
  vi_FE_mean, vi_FE_var, vi_FE_lndet,
  weight, it, old_FE_mean, old_weights,
  dim_fe, levels_fe, force_ascent, allow_decrease = FALSE){
  
  n_FE <- length(FE_data)
  vi_FE_raw <- as.list(rep(NA, n_FE))
  log_weight <- log(weight)
  
  for (v in 1:n_FE) {
    
    FE_rt_v <- FE_rowtens[[v]]
    FE_data_v <- FE_data[[v]]
    FE_lookup_v <- FE_lookup[[v]]
    
    init_mean_v <- rowSums(FE_data_v * (FE_lookup_v %*% vi_FE_mean[[v]]))
    init_var_v <- rowSums(FE_rt_v * (FE_lookup_v %*% vi_FE_var[[v]]))
    
    starting_obj_v <- sum(init_mean_v * y) - sum(weight) + 1/2 * vi_FE_lndet[v]
    names(starting_obj_v) <- NULL
    # NVMP update
    if (dim_fe[v] == 1){
      # NVMP_FE <- fast_1D_FE(X = FE_data_v, Z = FE_lookup_v, weights = weight, y = y - weight)
      NVMP_FE <- fast_generic_FE(
        FS_XX = FE_rt_v, X = FE_data_v, Z = FE_lookup_v,
         y = y - weight, weights = weight,
         dim_fe = dim_fe[v], levels_fe = levels_fe[v])
      
    }else{
      NVMP_FE <- fast_generic_FE(FS_XX = FE_rt_v, X = FE_data_v, Z = FE_lookup_v,
                                 y = y - weight, weights = weight,
                                 dim_fe = dim_fe[v], levels_fe = levels_fe[v])
      
    }
    NVMP_FE$mean <- NVMP_FE$mean + vi_FE_mean[[v]]
    nvmp_mean_v <- rowSums(FE_data_v * (FE_lookup_v %*% NVMP_FE$mean))
    nvmp_var_v <- rowSums(FE_rt_v * (FE_lookup_v %*% NVMP_FE$var))
    
    nvmp_log_weight <- log_weight + 
      (nvmp_mean_v - init_mean_v) + 1/2 * (nvmp_var_v - init_var_v)
    nvmp_weight <- exp(nvmp_log_weight)
    nvmp_lndet <- NVMP_FE$lndet
    obj_nvmp <- sum(nvmp_mean_v * y) - sum(nvmp_weight) + 1/2 * nvmp_lndet
    
    if ( (allow_decrease | ((abs(obj_nvmp) < 1e12) & obj_nvmp >= starting_obj_v)) & !force_ascent){
      
      log_weight <- nvmp_log_weight
      weight <- nvmp_weight
      vi_FE_mean[[v]] <- NVMP_FE$mean
      vi_FE_var[[v]] <- NVMP_FE$var
      vi_FE_lndet[v] <- nvmp_lndet
      vi_FE_raw[[v]] <- NVMP_FE$raw
      type_v <- 'NVMP'
      
    }else{
      
      if (is.finite(obj_nvmp) & it > 1 & !force_ascent){
        damp_par <- tryCatch(damp_poisson_FE_NVMP(attempt_NVMP = NVMP_FE, 
          starting_obj = starting_obj_v, y = y, 
          damp_offset_weight = log_weight - init_mean_v - 1/2 * init_var_v,
          FE_data_v = FE_data_v, FE_lookup_v = FE_lookup_v, 
          FE_rt_v = FE_rt_v, dim_fe_v = dim_fe[v], levels_fe_v = levels_fe[v],
          old_mean = vi_FE_mean[[v]], 
          old_weights = weight,
          new_weights = nvmp_weight, 
          direct_optimize = TRUE), error = function(e){NULL})
      }else{
        damp_par <- NULL
      }
      
      if (is.null(damp_par) | force_ascent){
        
        init_opt_par <- c(as.vector(vi_FE_mean[[v]]), as.vector(vi_FE_var[[v]]))
        
        opt_FE_manual <- tryCatch(optim(par = init_opt_par, 
              fn = obj_FE_poisson, gr = grad_FE_poisson,
              method = 'CG', control = list(fnscale = -1, maxit = 30),
              offset_weight = log_weight - init_mean_v - 1/2 * init_var_v,
              y = y, X = FE_data_v, P = levels_fe[v], 
              FS_XX = FE_rt_v,
              Z = FE_lookup_v), error = function(e){NULL})
        
        if (is.null(opt_FE_manual)){
          print("CG FAILED")
          diff_obj <- -Inf
        }else{
          if (opt_FE_manual$convergence != 0){
            opt_FE_manual$value <- obj_FE_poisson(par = opt_FE_manual$par, 
                           offset_weight = log_weight - init_mean_v - 1/2 * init_var_v,
                           y = y, X = FE_data_v, P = levels_fe[v], 
                           FS_XX = FE_rt_v,
                           Z = FE_lookup_v)
            
          }
          diff_obj <- opt_FE_manual$value - starting_obj_v
        }
        if (diff_obj > 0){
          print("OPTIMIZE")
          
          opt_mean_v <- matrix(opt_FE_manual$par[1:levels_fe[v]])
          opt_mean_v <- sweep(opt_mean_v, MARGIN = 2, colMeans(opt_mean_v), '-')
          ln_opt_var_v <- opt_FE_manual$par[-1:-levels_fe[v]]
          exp_opt_var_v <- exp(ln_opt_var_v)
          opt_var_v <- exp_opt_var_v - exp_opt_var_v^2/sum(exp_opt_var_v)
          opt_lndet_v <- sum(ln_opt_var_v) - log(sum(exp_opt_var_v))
          
          opt_mean_i <- rowSums(FE_data_v * (FE_lookup_v %*% opt_mean_v))
          opt_log_weight <- log_weight - init_mean_v - 1/2 * init_var_v +
            opt_mean_i + 1/2 * rowSums(FE_rt_v * (FE_lookup_v %*% opt_var_v))
          opt_weight <- exp(opt_log_weight)

          check_obj <- sum(opt_mean_i * y) - sum(opt_weight) + 1/2 * opt_lndet_v
          if (abs(opt_FE_manual$value - check_obj) > 1e-6){
            browser()
          }
          weight <- opt_weight
          log_weight <- opt_log_weight
          vi_FE_mean[[v]] <- opt_mean_v
          vi_FE_var[[v]] <- matrix(opt_var_v)
          vi_FE_lndet[v] <- opt_lndet_v
          vi_FE_raw[[v]] <- matrix(exp_opt_var_v)
          type_v <- 'OPTIMIZE'
          
        }else{
          print(diff_obj)
          log_weight <- nvmp_log_weight
          weight <- nvmp_weight
          vi_FE_mean[[v]] <- NVMP_FE$mean
          vi_FE_var[[v]] <- NVMP_FE$var
          vi_FE_lndet[v] <- nvmp_lndet
          vi_FE_raw[[v]] <- NVMP_FE$raw
          type_v <- 'FAILURE'
          print("FAILURE")
        }
      }else{
        log_weight <- damp_par$log_weight
        weight <- damp_par$weight
        vi_FE_mean[[v]] <- damp_par$par$mean
        vi_FE_var[[v]] <- damp_par$par$var
        vi_FE_lndet[v] <- damp_par$par$lndet
        vi_FE_raw[[v]] <- damp_par$raw
        type_v <- 'DAMP'
      }
    }
  }
  print(type_v)
  return(list(
    type = type_v,
    weight = weight,
    mean = vi_FE_mean,
    var = vi_FE_var,
    lndet = vi_FE_lndet,
    raw = vi_FE_raw
  ))
}


damp_poisson_FE_NVMP <- function(attempt_NVMP, starting_obj, y, 
  FE_data_v, FE_lookup_v, FE_rt_v, old_mean, old_weights,
  new_weights, dim_fe_v, levels_fe_v, damp_offset_weight,
  direct_optimize = TRUE){
  
  if (direct_optimize){
    
    f <- function(alpha_NVMP){
      
      weight_alpha <- new_weights * alpha_NVMP + (1 - alpha_NVMP) * old_weights
      weight_mu <- attempt_NVMP$mean * alpha_NVMP + (1 - alpha_NVMP) * old_mean
      
      damp_FE <- fast_generic_FE(FS_XX = FE_rt_v, X = FE_data_v, Z = FE_lookup_v,
                      y = y - weight_alpha, weights = weight_alpha,
                      dim_fe = dim_fe_v, levels_fe = levels_fe_v)
      damp_FE$mean <- damp_FE$mean + weight_mu
      damp_mean_v <- rowSums(FE_data_v * (FE_lookup_v %*% damp_FE$mean))
      damp_var_v <- rowSums(FE_rt_v * (FE_lookup_v %*% damp_FE$var))
      
      damp_log_weight <- damp_offset_weight + damp_mean_v + 1/2 * damp_var_v
      damp_weight <- exp(damp_log_weight)
      damp_lndet <- damp_FE$lndet
      
      damp_obj <- sum(damp_mean_v * y) - sum(damp_weight) + 1/2 * damp_lndet

      return(damp_obj)
    }
    
    grid_alpha <- optimize(f = f, interval = c(0,1), maximum = TRUE)$maximum
    
  } else {
    grid_alpha <- seq(0.5, 1, length.out = 11)
    
  }
  
  counter_NVMP <- 0
  while (counter_NVMP <= 10){
    
    alpha_NVMP <- grid_alpha[counter_NVMP + 1]
    if (is.na(alpha_NVMP)){break}
    
    weight_alpha <- new_weights * alpha_NVMP + (1 - alpha_NVMP) * old_weights
    weight_mu <- attempt_NVMP$mean * alpha_NVMP + (1 - alpha_NVMP) * old_mean
    
    damp_FE <- fast_generic_FE(FS_XX = FE_rt_v, X = FE_data_v, Z = FE_lookup_v,
      y = y - weight_alpha, weights = weight_alpha,
      dim_fe = dim_fe_v, levels_fe = levels_fe_v)
    damp_FE$mean <- damp_FE$mean + weight_mu
    damp_mean_v <- rowSums(FE_data_v * (FE_lookup_v %*% damp_FE$mean))
    damp_var_v <- rowSums(FE_rt_v * (FE_lookup_v %*% damp_FE$var))
    
    damp_log_weight <- damp_offset_weight + damp_mean_v + 1/2 * damp_var_v
    damp_weight <- exp(damp_log_weight)
    damp_lndet <- damp_FE$lndet
    damp_obj <- sum(damp_mean_v * y) - sum(damp_weight) + 1/2 * damp_lndet
    
    if (damp_obj >= starting_obj){
      damp_raw <- damp_FE$raw
      break
    }else{
      counter_NVMP <- counter_NVMP + 1
    }
  }
  if (damp_obj > starting_obj){
    return(list(par = damp_FE, raw = damp_raw,
                weight = damp_weight, log_weight = damp_log_weight,
                obj = damp_obj))
  }else{
    return(NULL)
  }
}

obj_FE_poisson <- function(par, y, X, FS_XX, Z, P, offset_weight){
  
  raw_beta <- par[1:P]
  beta <- matrix(raw_beta - mean(raw_beta))
  ln_diag_V <- par[-1:-P]
  raw_V <- exp(ln_diag_V)
  V <- raw_V - raw_V^2/sum(raw_V)
  
  mean_v <- rowSums(X * (Z %*% beta))
  weight <- exp(offset_weight + mean_v + 1/2 * rowSums(FS_XX * (Z %*% V)))
  
  lndet <- sum(ln_diag_V) - log(sum(raw_V))
  out <- sum( mean_v * y) - sum(weight) + 1/2 * lndet
  return(out)
}

grad_FE_poisson <- function(par, y, X, FS_XX, Z, P, offset_weight){
  
  raw_beta <- par[1:P]
  beta <- matrix(raw_beta - mean(raw_beta))
  ln_diag_V <- par[-1:-P]
  raw_V <- exp(ln_diag_V)
  V <- raw_V - raw_V^2/sum(raw_V)
  
  mean_v <- rowSums(X * (Z %*% beta))
  weight <- exp(offset_weight + mean_v + 1/2 * rowSums(FS_XX * (Z %*% V)))
  
  grad_beta <- (Diagonal(n = P) - matrix(1/P, P, P)) %*% t(Z) %*% Diagonal(x = as.vector(y - weight)) %*% X
  
  grad_V_weight <- Diagonal(x = 1 - 2 * raw_V/sum(raw_V)) + (raw_V/sum(raw_V))^2 %*% matrix(1, ncol = P)
  grad_V_weight <- grad_V_weight %*% Diagonal(x = raw_V)
  
  grad_V_weight <- -1/2 * as.vector(t(grad_V_weight) %*% t(Z) %*% Diagonal(x = weight) %*% FS_XX)
  
  grad_V_lndet <- 1/2 * (1 - raw_V/sum(raw_V))
  grad_V <- grad_V_weight + grad_V_lndet
  
  return(c(as.vector(grad_beta), as.vector(grad_V)))
}
