VB_nakajima <- function(svd_Y, c_a, c_b, sigmasq, H, Y,
                        M = NULL, L = NULL){
  
  if ( (length(c_a) != H) | (length(c_b) != H) ){
    stop('c_a and c_b must be of length "H"')
  }
  if (length(svd_Y$d) != H){
    stop('svd must only compute first H singular vectors')
  }
  
  if (L > M){
    transpose <- TRUE
    
    t_c_a <- c_a
    t_c_b <- c_b
    
    c_a <- t_c_b
    c_b <- t_c_a
    
    temp_u <- svd_Y$u
    temp_v <- svd_Y$v
    
    svd_Y$u <- temp_v
    svd_Y$v <- temp_u
    
    temp_M <- M
    temp_L <- L
    
    L <- temp_M
    M <- temp_L
    
    gc()
    
  }else{
    transpose <- FALSE
  }
  
  if (nrow(svd_Y$v) != M){stop('Mis-aligned')}
  if (nrow(svd_Y$u) != L){stop('Mis-aligned')}
  
  # Get the singular values from this SVD
  gamma <- svd_Y$d
  t1 <- (L + M)/2 * sigmasq + sigmasq^2/(2 * c_a * c_b)
  tilde_gamma <- sqrt(t1 + sqrt(t1^2 - L * M * sigmasq^2))
  # If gamma_h < tilde{gamma}_h, then zero-out this component
  zero_component <- gamma < tilde_gamma
  
  hat_gamma_VB <- sapply(1:H, FUN=function(h){
    if (zero_component[h]){
      return(0)
    }
    gamma_h <- gamma[h]
    c_ah <- c_a[h]
    c_bh <- c_b[h]
    
    etasq_h <- (1 - sigmasq * L/gamma_h^2) * (1 - sigmasq * M/gamma_h^2) * gamma_h^2
    xi_0 <- (etasq_h - sigmasq^2/(c_ah * c_bh))^2
    xi_3 <- (L - M)^2/(L * M) * gamma_h
    xi_2 <- -(xi_3 * gamma_h + (L^2 + M^2)/(L * M) * etasq_h + 2 * sigmasq^2/(c_ah * c_bh))
    xi_1 <- xi_3 * sqrt(xi_0)
    
    poly_coef <- c(xi_0, xi_1, xi_2, xi_3, 1)
    poly_root <- polyroot(poly_coef)
    real_root <- Re(poly_root)
    im_root <- Im(poly_root)
    if (!(sum(real_root > 0) == 2)){
      stop('error in nakajima; set mi_init="random"')
    }
    if (max(abs(im_root)) > 1e-4){
      stop('error in nakajima; set mi_init="random"')
      browser()
    }
    real_root <- real_root[real_root > 0]
    # Return the *second largest* real root
    sort_root <- sort(real_root, decreasing = TRUE)
    second_root <- sort_root[2]
    if (second_root < 0){
      stop('error in nakajima; set mi_init="random"')
      browser()
    }
    return(second_root)
  })
  
  if (any(gamma < hat_gamma_VB)){
    stop('...')
  }
  hat_delta_h <- (M-L) * (gamma - hat_gamma_VB) 
  hat_delta_h <- 
    (hat_delta_h + sqrt(hat_delta_h^2 + 4 * sigmasq^2 * L * M / (c_a * c_b))) *
    1/(2 * sigmasq * M * 1/c_a)
  etasq <- (1 - sigmasq * L/gamma^2) * (1 - sigmasq * M/gamma^2) * gamma^2
  hat_etasq <- ifelse(gamma > tilde_gamma, etasq, sigmasq^2 / (c_a * c_b))
  
  A_mean <- as.matrix(svd_Y$v %*% Diagonal(x=sqrt(hat_gamma_VB * hat_delta_h)))
  B_mean <- as.matrix(svd_Y$u %*% Diagonal(x=sqrt(hat_gamma_VB * 1/hat_delta_h)))
  
  t2_a <- (hat_etasq - sigmasq * (M-L))  
  A_var <- (-t2_a + sqrt(t2_a^2 + 4 * M * sigmasq * hat_etasq))/
    (2 * M * (hat_gamma_VB * 1/hat_delta_h + sigmasq * 1/c_a))
  t2_b <- (hat_etasq + sigmasq * (M-L))  
  B_var <- (-t2_b + sqrt(t2_b^2 + 4 * L * sigmasq * hat_etasq))/
    (2 * L * (hat_gamma_VB * hat_delta_h + sigmasq * 1/c_b))
  
  # if (transpose){
  #   Y <- t(Y)
  # }
  # 
  # sapply(1:H, FUN=function(h){
  #   SVB_a <- sigmasq * 1/(sum(B_mean[,h]^2) + L * B_var[h] + sigmasq/c_a[h])
  #   SVB_b <- sigmasq * 1/(sum(A_mean[,h]^2) + M * A_var[h] + sigmasq/c_b[h])
  #   mean_A <- A_var[h]/sigmasq * t(Y - B_mean[,-h] %*% t(A_mean[,-h])) %*% B_mean[,h]
  #   mean_B <- B_var[h]/sigmasq * (Y - B_mean[,-h] %*% t(A_mean[,-h])) %*% A_mean[,h]
  #   return(as.vector(mean_B))
  # })
  
  A_var <- t(array(as.vector(Diagonal(x=A_var)), c(H^2, M)))
  B_var <- t(array(as.vector(Diagonal(x=B_var)), c(H^2, L)))
  
  if (transpose){
    out <- list(
      U_mean = A_mean, V_mean = B_mean,
      U_var = A_var, V_var = B_var
    )
  }else{
    out <- list(V_mean = A_mean, U_mean = B_mean,
                V_var = A_var, U_var = B_var)
  }
  return(out)
}

# U V^T + W odot (Y - U V^T)
svd_A <- function(x, args){
  as.vector(
    (args$mean_U) %*% t(x %*% args$mean_V) +
      args$TERM_2 %*% x
  )
}
svd_Atrans <- function(x, args){
  as.vector(
    args$mean_V %*% t(x %*% args$mean_U) +
      t(x %*% args$TERM_2)
  )
}

sum_sparse <- function(A,B, f){
  sparseMatrix(i = A@i + 1,
               p = A@p,
               x = f(A@x, B@x),
               dim = dim(A),
               dimnames = dimnames(A))
}

MM_nakajima <- function(y, w, init_U, init_V,
                        Z_U, Z_V, prior_U, prior_V, D, iter, 
                        return_nonzero = FALSE,
                        tol = 1e-4){
  
  # Get the sum of weights for each (i,j) combination
  summed_weights <- t(Z_U) %*% Diagonal(x=w) %*% Z_V
  weighted_sum_Y <- (t(Z_U) %*% Diagonal(x=y * w) %*% Z_V)
  # Older but slower way
  # recip_W <- summed_weights
  # recip_W@x <- 1/recip_W@x
  # weighted_avg_Y <- weighted_sum_Y * recip_W
  weighted_avg_Y <- sum_sparse(weighted_sum_Y, summed_weights, f = `/`)
  max_weight <- max(summed_weights@x)
  mean_U <- init_U
  mean_V <- init_V
  
  safe_U <- old_U <- init_U
  safe_V <- old_V <- init_V
  
  norm_weight <- summed_weights
  norm_weight@x <- norm_weight@x/max_weight
  
  pos_nonzero <- tryCatch(as(norm_weight, 'TsparseMatrix'), error = function(e){NULL})
  if (is.null(pos_nonzero)){
    pos_nonzero <- as(norm_weight, 'dgTMatrix')
  }
  pos_nonzero_U <- pos_nonzero@i
  pos_nonzero_V <- pos_nonzero@j
  pos_nonzero_x <- pos_nonzero@x
  
  
  pos_sparse_Y <- tryCatch(as(weighted_avg_Y, 'TsparseMatrix'), error = function(e){NULL})
  if (is.null(pos_nonzero)){
    pos_sparse_Y <- as(weighted_avg_Y, 'dgTMatrix')
  }
  # Check identical positions
  c1 <- isTRUE(all.equal(pos_sparse_Y@i, pos_nonzero@i))
  c2 <- isTRUE(all.equal(pos_sparse_Y@j, pos_nonzero@j))
  if (!c1 | !c2){stop('identical positions not found...')}
  pos_sparse_Y_x <- pos_sparse_Y@x
  
  
  # sparse_WY <- norm_weight * weighted_avg_Y
  # inv_W_matrix <- 1 - norm_weight
  for (it in 1:iter){
    
    # # Direct SVD for MM-Nakajima, not advised for large problems...  
    # weight_Y <-
    #   sparse_WY +
    #   (mean_U %*% t(mean_V)) * inv_W_matrix
    # direct_svd_Y <- RSpectra::svds(A = weight_Y, k = D)
    # svd_Y <- direct_svd_Y
    
    # Use the function interface + sparse estimation
    # This does not work when Z_U or Z_V contains multiple elements
    # impute_UVt <- rowSums((Z_U %*% mean_U) * (Z_V %*% mean_V))
    # impute_UVt <-  t(Z_U) %*% Diagonal(x=impute_UVt) %*% Z_V
    # This works but is slow on big problems
    # impute_UVt <- rowSums(mean_U[pos_nonzero_U,] * mean_V[pos_nonzero_V,])
    impute_UVt <- cpp_zipped_sum(pos_U = pos_nonzero_U, pos_V = pos_nonzero_V,
                                 Ut = t(mean_U), Vt = t(mean_V))
    TERM_2 <- sparseMatrix(i = pos_nonzero_U + 1,
                           j = pos_nonzero_V + 1,
                           x = pos_nonzero_x * (pos_sparse_Y_x - impute_UVt) ,
                           dims = dim(weighted_avg_Y))
    args_svd <- list(TERM_2 = TERM_2, mean_U = mean_U, mean_V = mean_V)
    if (!requireNamespace('RSpectra', quietly = TRUE)){
      stop('RSpectra must be installed for Nakajima initialization')
    }else{
      svd_Y <- RSpectra::svds(
        A = svd_A, Atrans = svd_Atrans, k = D,
        dim = c(ncol(Z_U), ncol(Z_V)),
        args = args_svd)
    }
    
    
    # Y is L x M where the decomposition is B A^T where
    # B is L x H and A is M x H
    # So, in our notation, A is V and B is U...
    
    fit <- VB_nakajima(svd_Y = svd_Y, Y = NULL,
                       # A: Has prior c_a and M rows (corresponds to V)
                       c_a = prior_V, M = ncol(Z_V),
                       # B: Has prior c_b and L rows (corresponds to U)
                       c_b = prior_U, L = ncol(Z_U),
                       sigmasq = 1/max_weight, H = D)
    
    # \sum_i \sum_j (tilde{y}_{ij} - a_i^T b_j)^2 + prior()
    # \sum_j b_j (tilde_{y}_{ij} - b_j^T a_i)  + prior()
    # \sum_j tilde{y}_{ij} b_j - a_i^T \sum_j b_j^T
    
    # Y %*% V[,]
    grad_U <- sapply(1:D, FUN=function(h){
      svd_A(x = fit$V_mean[,h], args_svd)
    })
    grad_U <- 
      max_weight * (
        grad_U - fit$U_mean %*% crossprod(fit$V_mean) +
          - fit$U_mean %*% matrix(colSums(fit$V_var), D)
      ) +
      - fit$U_mean/prior_U
    # t(Y) %*% U[,h]
    grad_V <- sapply(1:D, FUN=function(h){
      svd_Atrans(x = fit$U_mean[,h], args_svd)
    })
    grad_V <- 
      max_weight * (
        grad_V - fit$V_mean %*% crossprod(fit$U_mean) +
          - fit$V_mean %*% matrix(colSums(fit$U_var), D)
      ) +
      - fit$V_mean/prior_V
    
    if (max(abs(grad_U)) > 1e-5 | max(abs(grad_V)) > 1e-5){
      warning('Error in NAKAJIMA ERROR init')
      # browser()
    }
    
    mean_U <- fit$U_mean  
    mean_V <- fit$V_mean  
    var_U <- fit$U_var
    var_V <- fit$V_var

    change_U <- max(abs(mean_U - old_U))
    change_V <- max(abs(mean_V - old_V))
    
    old_U <- mean_U
    old_V <- mean_V
    # print(c(change_U, change_V))
    
    if (max(abs(mean_U)) > 0){
      safe_U <- mean_U
    }
    
    if (max(abs(mean_V)) > 0){
      safe_V <- mean_V
    }
    
    if (change_U < tol & change_V < tol){
      break
    }
  }
  
  if (return_nonzero){
    if (all(mean_U == 0) | all(mean_V == 0)){
      warning('MM + Nakjaima zeroed all terms; returning last non-zeroed Nakajima solution')
    }
    mean_U <- safe_U
    mean_V <- safe_V
  }
  
  return(list(
    mean_U = mean_U, mean_V = mean_V,
    var_U = var_U, var_V = var_V
  ))
}


init_MI_from_svd <- function(data_mi, y, pg_weight, rank, prior_U, prior_V){
  
  minimum_size <- min(sapply(data_mi, ncol))
  if (rank > minimum_size){
    stop('rank > maximum size; lower to continue')
  }
  
  init_U <- array(matrix(rnorm(ncol(data_mi[[1]]) * rank, sd = sqrt(1/rank))), c(ncol(data_mi[[1]]), rank))
  init_V <- array(matrix(rnorm(ncol(data_mi[[2]]) * rank, sd = sqrt(1/rank))), c(ncol(data_mi[[2]]), rank))
  
  init_SVD <- MM_nakajima(y = y,
              w = pg_weight, 
              init_U = init_U, init_V = init_V,
              prior_U = rep(prior_U, rank),
              prior_V = rep(prior_V, rank),
              Z_U = data_mi[[1]], 
              Z_V = data_mi[[2]],
              D = rank,
              return_nonzero = TRUE,
              iter = 5)
  
  counter <- 0
  for (r in 1:rank){
    
    if (all(init_SVD$mean_U[,r] == 0)){
      init_SVD$mean_U[,r] <- rnorm(nrow(init_SVD$mean_U), sd = sqrt(1/rank))
      counter <- counter + 1
    }
    
    if (all(init_SVD$mean_V[,r] == 0)){
      init_SVD$mean_V[,r] <- rnorm(nrow(init_SVD$mean_V), sd = sqrt(1/rank))
      counter <- counter + 1
    }
    
  }
  if (counter > 0){
    warning('Nakajima initalization zeroed out; random initialization used instead')
  }
  return(init_SVD)
}
