
FS <- function(X,Y){t(KhatriRao(t(X),t(Y)))}

get_bilinear_mean <- function(Z_MI, vi_mi_mean, vi_hier_grouping, reduce = TRUE){
  
  if (length(Z_MI) != length(vi_mi_mean)){
    stop('lengths misaligned')
  }
  if (length(Z_MI) != length(vi_hier_grouping)){
    stop("lengths misaligned")
  }
  
  if (reduce){
    out <- Reduce("+", mapply(Z_MI, vi_mi_mean, vi_hier_grouping, 
      SIMPLIFY = FALSE, FUN=function(mi_data, mi_mean, mi_grouping){
        # only one group in each term
        simple_model <- all(lengths(mi_grouping) == 1)
        if (simple_model){
          out_j <- rowSums(Reduce("*", mapply(mi_data, mi_mean, SIMPLIFY = FALSE, FUN=function(i,j){
            i %*% j
          })))
        }else{
          out_j <- rowSums(Reduce("*", lapply(mi_grouping, FUN=function(g){
            # For each group, add and *then* multiply u and v
            Reduce('+', lapply(g, FUN=function(g_i){
              mi_data[[g_i]] %*% mi_mean[[g_i]]
            }))
          })))
        }
        return(out_j)
    }))
  }else{
    out <- mapply(Z_MI, vi_mi_mean, SIMPLIFY = FALSE, FUN=function(mi_data, mi_mean){
      mapply(mi_data, mi_mean, SIMPLIFY = FALSE, FUN=function(i,j){
        i %*% j
      })
    })
  }
  return(out)
}

get_bilinear_var <- function(Z_MI, vi_mi_mean, vi_mi_var, vi_hier_grouping, reduce = TRUE){
  
  if (length(Z_MI) != length(vi_mi_mean)){
    stop('lengths misaligned')
  }
  if (length(Z_MI) != length(vi_mi_var)){
    stop('lengths misaligned')
  }
  if (length(Z_MI) != length(vi_hier_grouping)){
    stop("lengths misaligned")
  }
  
  
  if (reduce){
    out <- mapply(Z_MI, vi_mi_mean, vi_mi_var, vi_hier_grouping, SIMPLIFY = FALSE, 
        FUN=function(mi_data, mi_mean, mi_var, mi_grouping){
          simple_model <- all(lengths(mi_grouping) == 1)
          if (!simple_model){
            M_mean <- mapply(mi_data, mi_mean, SIMPLIFY = FALSE, FUN=function(i,j){i %*% j})
            # Sum by group to get the mean within "u" and "v"
            M_mean <- lapply(mi_grouping, FUN=function(g){
              # For each group, add and *then* multiply u and v
              Reduce('+', lapply(g, FUN=function(g_i){
                M_mean[[g_i]]
              }))
            })
            M_mean <- lapply(M_mean, FUN=function(i){FS(i,i)})
            # For Var(u_i + u_p) given Var(u_i) and Var(u_p) are independent
            M_var <- mapply(mi_data, mi_var, SIMPLIFY = FALSE, FUN=function(i,j){i %*% j})
            M_var <- lapply(mi_grouping, FUN=function(g){
              # For each group, add together the variance
              Reduce('+', lapply(g, FUN=function(g_i){
                M_var[[g_i]]
              }))
            })
          }else{
            mi_FS_mean <- lapply(mi_mean, FUN=function(i){FS(i,i)})
            M_mean <- mapply(mi_data, mi_FS_mean, SIMPLIFY = FALSE, FUN=function(i,j){i %*% j})
            M_var <- mapply(mi_data, mi_var, SIMPLIFY = FALSE, FUN=function(i,j){i %*% j})
          }
          out_trace <- rowSums(Reduce('*', M_var))
          out_quad <- 
            rowSums( M_mean[[1]] * M_var[[2]] ) +
            rowSums( M_mean[[2]] * M_var[[1]] )
          return(out_trace + out_quad)
        
      })
    return(Reduce("+", out))
  }else{
    out <- mapply(Z_MI, vi_mi_mean, vi_mi_var, SIMPLIFY = FALSE, 
    FUN=function(mi_data, mi_mean, mi_var){
      M_var <- mapply(mi_data, mi_var, SIMPLIFY = FALSE, FUN=function(i,j){i %*% j})
      return(M_var)
    })
    return(out) 
  }
}

get_bilinear_outer <- function(vi_mi_mean, vi_mi_var, vi_mi_diag, mi_prior_type,  reduce = TRUE){

  mapply(
    vi_mi_mean, vi_mi_var, vi_mi_diag, 
    SIMPLIFY = FALSE,
    FUN=function(m_mean, m_var, m_diag, m_hier){
      if (reduce & (mi_prior_type != 'separate')){
        if (mi_prior_type %in% c('centered')){
          out_mean <- Reduce("+", lapply(m_mean[2], crossprod))
          out_var <- Reduce("+", lapply(m_var[2], colSums))
        }else if (mi_prior_type %in% c('shared')){
          out_mean <- Reduce("+", lapply(m_mean, crossprod))
          out_var <- Reduce("+", lapply(m_var, colSums))
        }
        out_var <- matrix(out_var, nrow = length(m_diag))
        return(out_mean + out_var)
      }else{
        out <- mapply(lapply(m_mean, crossprod), lapply(m_var, colSums), 
               SIMPLIFY = FALSE, FUN=function(i,j){
          i + matrix(j, nrow = length(m_diag))
        }) 
        return(out)
      }
    }
  )
}

#' @importFrom RSpectra svds
init_MI_from_svd <- function(data_mi, weight, rank, mi_prior_type, y = NULL, trials = NULL){

  f <- function(x, args){
    as.vector(t(args$Z1) %*% Diagonal(x=args$w) %*% (args$Z2 %*% x))
  }
  g <- function(x, args){
    as.vector(t(args$Z2) %*% Diagonal(x=args$w) %*% (args$Z1 %*% x))
  }

  if (init_mi_type == 'svd'){
    warning('init naive SVD: Setting d=identity; STD d by sd(vec(D))')
    if (is.null(y)){
      D <- t(data_mi[[1]]) %*% Diagonal(x=weight) %*% data_mi[[2]]
      D <- D/sd(as.vector(D@x))
    }else{
      ratio <- (y/trials > 0.5)
      ratio[is.na(ratio)] <- 1/2
      D <- t(data_mi[[1]]) %*% Diagonal(x=2 * ratio  - 1) %*% data_mi[[2]]
      D <- D/sd(as.vector(D@x))
    }
    D <- sweep(D, MARGIN = 1, STATS = rowMeans(D), FUN = '-')
    D <- sweep(D, MARGIN = 2, STATS = colMeans(D), FUN = '-')
    svd_weight <- RSpectra::svds(A = D, k = rank)
    svd_weight$d <- svd_weight$d/svd_weight$d
  }else if (init_mi_type == 'peress'){
    warning('init bad peress')
    if (is.null(y)){
      svd_peress <- ipe::ipe_start(Y = t(data_mi[[1]]) %*%
                                     Diagonal(x= as.numeric(weight < median(weight)) + 1) %*%
                                     data_mi[[2]], D = rank)
    }else{
      ratio <- (y/trials > 0.5)
      ratio[is.na(ratio)] <- 1/2
      svd_peress <- ipe::ipe_start(Y = t(data_mi[[1]]) %*%
                                     Diagonal(x= ratio) %*%
                                     data_mi[[2]], D = rank)
    }
    svd_weight <- svd_peress[c('Alpha', 'Delta')]
    names(svd_weight) <- c('u', 'v')
    svd_weight$v <- svd_weight$v[,-1,drop=T]
    svd_weight$d <- rep(1, rank)
  }else{stop('....')}
  # minimum_size <- min(sapply(data_mi, ncol))
  # if (rank >= minimum_size){
  #   warning('rank of a multiplicative interaction is above minimum size...')
  #   svd_weight <- svd(t(data_mi[[1]]) %*% Diagonal(x=weight) %*% data_mi[[2]])
  #   # Pad with zeros
  #   svd_weight$d <- c(svd_weight$d, rep(0, rank - minimum_size))
  #   svd_weight$u <- cbind(svd_weight$u, matrix(0, nrow = ncol(data_mi[[1]]), ncol = rank - minimum_size))
  #   svd_weight$v <- cbind(svd_weight$v, matrix(0, nrow = ncol(data_mi[[2]]), ncol = rank - minimum_size))
  # }else{
  #   svd_weight <- RSpectra::svds(A = f, k = rank, Atrans = g, 
  #                                dim = sapply(data_mi, ncol),
  #                                args = list(Z1 = data_mi[[1]], Z2 = data_mi[[2]], weight = weight))
  # }

  
  if (mi_prior_type %in% c('shared', 'separate')){
    weight_d <- Diagonal(x=sqrt(svd_weight$d))
    svd_weight$u <- svd_weight$u  %*% weight_d
    svd_weight$v <- svd_weight$v  %*% weight_d
  }else{
    svd_weight$u <- svd_weight$u
    svd_weight$v <- svd_weight$v %*% Diagonal(x=svd_weight$d)  
  }
  return(svd_weight)
}

bilinear_update <- function(Z_MI, Z_hier_mapping, prior_mi, d_mi, s,
                            vi_pg_mean, diag_vi_pg_mean, RFSmean, Rvar, Rmean, offset,
                            return_chol, method, UF_MI, it,
                            auxiliary_Z = NULL, auxiliary_prior = NULL){
  
  simple_j <- !inherits(Z_MI, 'list')
  
  if (!simple_j){
    
    full_Z <- do.call('cbind', Z_MI)
    
    # aug_Z <- do.call('cbind', lapply(Z_MI, FUN=function(i){
    #   FS(i, Rmean)
    # }))
    aug_Z <- FS(full_Z, Rmean)
    aug_prior <- bdiag(
      mapply(Z_MI, prior_mi, SIMPLIFY = FALSE, 
             FUN=function(i, p_i){kronecker(Diagonal(n=ncol(i)), matrix(p_i, d_mi))})
    )
    
    if (d_mi == 1){
      aug_Rvar <- t(full_Z) %*% Diagonal(x=vi_pg_mean * Rvar[,1]) %*% full_Z
    }else{
      
      tile <- lapply(1:d_mi, FUN=function(i){
        lapply(1:d_mi, FUN=function(j){
          o <- t(full_Z) %*% Diagonal(x=vi_pg_mean * Rvar[,d_mi * (i-1) + j]) %*% full_Z
          o <- as(o, 'TsparseMatrix')
          o <- cbind(o@i * d_mi + i  , o@j * d_mi + j , o@x)
          return(o)
        })
      })
      flat_tile <- do.call('rbind', lapply(tile, FUN=function(i){do.call('rbind', i)}))
      aug_Rvar <- sparseMatrix(i = flat_tile[,1], j = flat_tile[,2], x = flat_tile[,3])
      # # Correct but extremely slow method...
      # slow_aug_Rvar <- Reduce('+', lapply(1:nrow(full_Z), FUN=function(i){
      #   z_i <- full_Z[i,,drop=F]
      #   kronecker(crossprod(z_i), matrix(vi_pg_mean[i] * Rvar[i,], d_mi))
      # }))
      # Old, incorrect, code only captures diagonal terms
      # aug_Rvar <- t(full_Z) %*% diag_vi_pg_mean %*% Rvar
      # aug_Rvar <- bdiag(apply(aug_Rvar, MARGIN = 1, FUN=function(i){Matrix(i,d_mi)}))
    }
    
    if (!is.null(auxiliary_Z)){
      if (UF_MI){stop('...')}
      aug_Z <- cbind(auxiliary_Z, aug_Z)
      aug_prior <- bdiag(auxiliary_prior, aug_prior)
      aug_Rvar <- bdiag(Diagonal(x = rep(0, ncol(auxiliary_Z))), aug_Rvar)
    }
    
    aug_chol <- Cholesky(t(aug_Z) %*% Diagonal(x=vi_pg_mean) %*% aug_Z + aug_prior + aug_Rvar)
    aug_coef <- solve(
      aug_chol,
      t(aug_Z) %*% (s - vi_pg_mean * offset)
    )

    if (!is.null(auxiliary_Z)){
      auxiliary_est <- aug_coef[1:ncol(auxiliary_Z)]
      aug_coef <- aug_coef[-(1:ncol(auxiliary_Z))]    
    }else{
      auxiliary_est <- NULL
    }
    
    out_mean <- matrix(aug_coef, ncol = d_mi, byrow = T)
    split_out <- sapply(Z_MI, ncol)
    split_out <- c(0, cumsum(split_out))
    out_mean <- lapply(1:length(Z_MI), FUN=function(i){
      out_mean[seq(split_out[i] + 1, split_out[i+1]),,drop=FALSE]
    })
    
    if (!is.null(Z_hier_mapping)){
      old_mean <- out_mean
      prior_mi <- lapply(prior_mi, FUN=function(i){matrix(i, nrow = sqrt(length(i)))})
      if (do_PX_HIER & it > HIER_THRESH){
        message('hier_PX ON')
        out_mean <- hier_px(Z_hier_mapping, out_mean, d_mi, prior_mi)
      }
    }
    
    ub <- mapply(Z_MI, prior_mi, out_mean, SIMPLIFY = FALSE, FUN=function(Z_j, prior_j, om_j){
      ub_j <- invert_rowwise(
        X = as.matrix(t(Z_j) %*% diag_vi_pg_mean %*% (
          RFSmean + Rvar)
        ),
        vec_prior = prior_j,
        dim = d_mi,
        RHS = array(0, dim = c(nrow(Z_j), d_mi)),
        return_chol = return_chol
      )
      ub_j$mean <- om_j
      return(ub_j)
    })
    ub$aux <- auxiliary_est
    return(ub)
    
  }else if (method == 'rowwise'){
    
    update_bilinear <- invert_rowwise(
      X = as.matrix(
        t(Z_MI) %*% diag_vi_pg_mean %*%
          (RFSmean + Rvar)
      ),
      vec_prior = prior_mi,
      dim = as.integer(d_mi),
      RHS = as.matrix(t(Z_MI) %*%
                        Diagonal(x= s - vi_pg_mean * offset) %*% Rmean),
      return_chol = return_chol
    )

    check <- max(abs(t(apply(as.matrix(
      t(Z_MI) %*% diag_vi_pg_mean %*%
        (RFSmean + Rvar)
    ), MARGIN = 1, FUN=function(i){
      as.vector(solve(matrix(i + prior_mi, 4)))
    })) - update_bilinear$inverse))
    if (check > 1e-7){stop("...")}
    print(dim(update_bilinear$mean))    
    print('var')
    print(rowMeans(apply(update_bilinear$inverse, MARGIN = 1, FUN=function(i){diag(matrix(i, ncol = d_mi))})))
    print('mean')
    print(colMeans(update_bilinear$mean))
    print('eigen')
    print(colMeans(t(apply(as.matrix(
      t(Z_MI) %*% diag_vi_pg_mean %*%
        (RFSmean + Rvar)
    ), MARGIN = 1, FUN=function(i){
      eigen(solve(matrix(i + prior_mi, 4)))$values
    }))))
    
    if (!is.null(Z_hier_mapping)){
      
      old_mean <- update_bilinear$mean
      browser()
      stopifnot(all(sapply(prior_mi, FUN=function(i){isDiagonal(matrix(i, nrow = sqrt(length(i))))})))
      update_bilinear$mean <- hier_px(Z_hier_mapping, update_bilinear$mean, d_mi)
    }
    
    return(update_bilinear)
  }else if (method == 'sparse_direct'){

    aug_Rvar <- t(Z_MI) %*% diag_vi_pg_mean %*% Rvar
    stop('WRONG aug_Rvar')
    aug_Rvar <- bdiag(apply(aug_Rvar, MARGIN = 1, FUN=function(i){Matrix(i,d_mi)}))
    if (!is.null(auxiliary_Z)){
      orig_Z <- FS(Z_MI, Rmean)
      aug_Z <- cbind(auxiliary_Z, orig_Z)
      orig_prior <- kronecker(Diagonal(n=ncol(Z_MI)), matrix(prior_mi, d_mi))
      aug_prior <- bdiag(auxiliary_prior, orig_prior)
      orig_Rvar <- aug_Rvar
      aug_Rvar <- bdiag(Diagonal(x = rep(0, ncol(auxiliary_Z))), aug_Rvar)
    }else{
      aug_Z <- FS(Z_MI, Rmean)
      aug_prior <- kronecker(Diagonal(n=ncol(Z_MI)), matrix(prior_mi, d_mi))
    }
    aug_chol <- Cholesky(t(aug_Z) %*% Diagonal(x=vi_pg_mean) %*% aug_Z + aug_prior + aug_Rvar)
    aug_coef <- solve(
      aug_chol,
      t(aug_Z) %*% (s - vi_pg_mean * offset)
    )

    if (!is.null(auxiliary_Z)){
      auxiliary_est <- aug_coef[1:ncol(auxiliary_Z)]
      aug_coef <- aug_coef[-(1:ncol(auxiliary_Z))]    
    }else{
      auxiliary_est <- NULL
    }
    out_mean <- matrix(aug_coef, ncol = d_mi, byrow = T)
    
    if (!is.null(Z_hier_mapping)){
      
      old_mean <- out_mean
      browser()
      stopifnot(all(sapply(prior_mi, FUN=function(i){isDiagonal(matrix(i, nrow = sqrt(length(i))))})))
      
      out_mean <- hier_px(Z_hier_mapping, out_mean, d_mi)
    }
    
    ub <- invert_rowwise(
      X = as.matrix(t(Z_MI) %*% diag_vi_pg_mean %*% (
        RFSmean + Rvar)
      ),
      vec_prior = prior_mi,
      dim = d_mi,
      RHS = array(0, dim = c(nrow(Z_MI), d_mi)),
      return_chol = return_chol
    )
    
    ub$mean <- out_mean
    ub$aux <- auxiliary_est
    return(ub)
    
  }else{
    stop('invalid method for bilinear update')
  }
}

diagonalize_AB <- function(A,B, tol = 1e-8){
  eigen_A <- eigen(A, symmetric = TRUE)
  if (any(eigen_A$values < max(eigen_A$values) * tol)){
    stop('A is not invertible...')
  }else{
    eigen_A$inv_values <- 1/eigen_A$values
  }
  # Q_A L_A Q_A^T
  T_A <- Diagonal(x = sqrt(eigen_A$inv_values)) %*% t(eigen_A$vectors)
  inv_T_A <- eigen_A$vectors %*% Diagonal(x=sqrt(eigen_A$values))
  transf_A <- T_A %*% A %*% t(T_A)
  transf_B <- t(inv_T_A) %*% B %*% inv_T_A
  eigen_transf_B <- eigen(transf_B, symmetric = TRUE)
  T_out <- t(eigen_transf_B$vectors) %*% T_A 
  
  drop0(zapsmall(t(solve(T_out)) %*% B %*% solve(T_out)))
  drop0(zapsmall(T_out %*% A %*% t(T_out)))
  return(T_out)
}

obj_MI_mean <- function(Z_MI, s, omega, mi_mean, mi_var){
  lp <- get_bilinear_mean(list(Z_MI), list(mi_mean))
  term_1 <- sum(s * lp)
  term_2a <- sum(omega * rowSums( (Z_MI[[1]] %*% mi_mean[[1]]) * (Z_MI[[2]] %*% mi_mean[[2]])))
  term_2b <- sum(sapply(1:2, FUN=function(k){
    not_k <- setdiff(1:2, k)
    RFSmean <- Z_MI[[k]] %*% FS(mi_mean[[k]], mi_mean[[k]])
    Rvar <- Z_MI[[not_k]] %*% mi_var[[not_k]]
    sum(rowSums(RFSmean * omega * Rvar))
  }))
  return(term_1 - (term_2a + term_2b))
}

reformat_mean_lbfgs <- function(par, dim_j, recons_j, names_j){
  mean_j <- lapply(split(par, recons_j), FUN=function(i){matrix(i, ncol = dim_j)})
  mean_j <- mean_j[names_j]
  return(mean_j)
}

fff <- function(mean_j, data_j, s, omega, long_var_j, offset_j,
                pos_d, ridge_j, grouping_j, recons_j, dim_j, basic = FALSE, disagg = FALSE){
  out_mean <- mapply(data_j, mean_j, SIMPLIFY = FALSE, FUN=function(i,j){
    i %*% j
  })
  if (basic){
    mean_u <- out_mean[[1]]
    mean_v <- out_mean[[2]]
  }else{
    mean_u <- Reduce('+', out_mean[grouping_j[[1]]])
    mean_v <- Reduce('+', out_mean[grouping_j[[2]]])
  }
  mean_lp <- rowSums(mean_u * mean_v)
  term_1 <- sum(s * (offset_j + mean_lp))
  term_2a <- -1/2 * sum(omega * (offset_j + mean_lp)^2)
  term_2b <- -1/2 * sum(omega * rowSums(FS(mean_u, mean_u) * long_var_j[[2]]))
  term_2c <- -1/2 * sum(omega * rowSums(FS(mean_v, mean_v) * long_var_j[[1]]))
  term_2 <- term_2a + term_2b + term_2c
  term_prior <- sum(mapply(mean_j, ridge_j, FUN=function(i,j){
    -1/2 * sum(rowSums( (i %*% j) * i))
  }))
  if (disagg){
    return(c(term_1, term_2a, term_2b, term_2c, term_prior))
  }else{
    return(term_1 + term_2 + term_prior)
  }
}

f_lbfgs_mean <- function(par, data_j, s, omega, long_var_j, offset_j, pos_d, ridge_j, grouping_j, recons_j, dim_j){
  
  mean_j <- reformat_mean_lbfgs(par, dim_j, recons_j, names(data_j))
  
  out_mean <- mapply(data_j, mean_j, SIMPLIFY = FALSE, FUN=function(i,j){
      i %*% j
  })
  mean_u <- Reduce('+', out_mean[grouping_j[[1]]])
  mean_v <- Reduce('+', out_mean[grouping_j[[2]]])
  mean_lp <- rowSums(mean_u * mean_v)
  term_1 <- sum(s * (offset_j + mean_lp))
  term_2a <- -1/2 * sum(omega * (offset_j + mean_lp)^2)
  term_2b <- -1/2 * sum(omega * rowSums(FS(mean_u, mean_u) * long_var_j[[2]]))
  term_2c <- -1/2 * sum(omega * rowSums(FS(mean_v, mean_v) * long_var_j[[1]]))
  term_2 <- term_2a + term_2b + term_2c
  term_prior <- sum(mapply(mean_j, ridge_j, FUN=function(i,j){
    -1/2 * sum(rowSums( (i %*% j) * i))
  }))
  
  return(-1 * (term_1 + term_2 + term_prior))
}

gr_lbfgs_mean <- function(par, data_j, s, omega, long_var_j, offset_j, ridge_j, pos_d, grouping_j, recons_j, dim_j){
  
  mean_j <- reformat_mean_lbfgs(par, dim_j, recons_j, names(data_j))
  out_mean <- mapply(data_j, mean_j, SIMPLIFY = FALSE, FUN=function(i,j){
    i %*% j
  })
  mean_u <- Reduce('+', out_mean[grouping_j[[1]]])
  mean_v <- Reduce('+', out_mean[grouping_j[[2]]])
  mean_lp <- rowSums(mean_u * mean_v)
  adj_s <- s - omega * offset_j
  deriv_mean <- lapply(1:2, FUN=function(k){
    
    group_k <- grouping_j[[k]]
    out_k <- lapply(group_k, FUN=function(l){
      Z_l <- data_j[[l]]
      if (k == 1){
        term_1 <- t(Z_l) %*% Diagonal(x=adj_s) %*% mean_v
        term_2 <- (long_var_j[[2]] + FS(mean_v, mean_v))
        term_2 <- Reduce('+', mapply(1:ncol(mean_u), pos_d, FUN=function(i,j){
          t(Z_l) %*% Diagonal(x = omega * mean_u[,i]) %*% term_2[,j]
        }))
      }else{
        term_1 <- t(Z_l) %*% Diagonal(x=adj_s) %*% mean_u
        term_2 <- (long_var_j[[1]] + FS(mean_u, mean_u))
        term_2 <- Reduce('+', mapply(1:ncol(mean_v), pos_d, FUN=function(i,j){
          t(Z_l) %*% Diagonal(x = omega * mean_v[,i]) %*% term_2[,j]
        }))
      }
      term_prior <- mean_j[[l]] %*% ridge_j[[l]]
      return(as.vector(term_1 - term_2 - term_prior))
    })
    return(do.call('c', out_k))
  })
  deriv_mean <- do.call('c', deriv_mean)
  return(-deriv_mean)
}

flatten_mi_mean <- function(x){
  do.call('c', lapply(x, FUN=function(i){as.vector(t(i))}))
}

get_bilinear_outer_UF <- function(vi_mi_mean, vi_mi_decomp, Z_MI_grouping, outer_MI_positions){
  mapply(vi_mi_decomp, vi_mi_mean, Z_MI_grouping, outer_MI_positions,
    SIMPLIFY = FALSE, FUN=function(decomp_i, mean_i, group_i, positions_i){
      out <- mapply(decomp_i, group_i, positions_i, SIMPLIFY = FALSE, FUN=function(decomp_j, group_j, positions_j){
        out_j <- calculate_expected_outer_alpha(
           L = as(decomp_j, 'dgCMatrix'),
           alpha_mu = flatten_mi_mean(mean_i[group_j]),
           re_position_list = positions_j
          )$outer_alpha
        names(out_j) <- group_j
        return(out_j)
      })
    out <- unlist(out, recursive = FALSE)[unlist(group_i)]
  })
}
 

get_longvar_UF <- function(d_mi, Z, pos_dim, var_decomp){
  
  fmt_var <- lapply(pos_dim, FUN=function(i){var_decomp[,i,drop=F]})
  out <- array(NA, dim = c(nrow(Z), d_mi^2))
  for (i in 1:d_mi){
    for (j in 1:d_mi){
      if (i < j){
        out[,d_mi * (j-1) + i] <- rowSums( (Z %*% t(fmt_var[[i]])) * (Z %*% t(fmt_var[[j]])) ) 
      }else if (i == j){
        out[,d_mi * (j-1) + i] <- rowSums( (Z %*% t(fmt_var[[i]]))^2 ) 
      }else{
        out[,d_mi * (j-1) + i] <- out[, d_mi * (i-1) + j]
      }
    }
  }
  return(out)
}

get_bilinear_var_UF <- function(Z_MI, Z_MI_first, Z_MI_first_mapping, 
                                vi_mi_mean, vi_mi_decomp, vi_hier_grouping, 
                                dim_MI_positions, dim_MI, reduce = TRUE){
  
  if (length(Z_MI) != length(vi_mi_mean)){
    stop('lengths misaligned')
  }
  if (length(Z_MI) != length(vi_mi_decomp)){
    stop('lengths misaligned')
  }
  if (length(Z_MI) != length(vi_hier_grouping)){
    stop("lengths misaligned")
  }
  if (reduce){

    out <- mapply(Z_MI, Z_MI_first, Z_MI_first_mapping, 
        vi_mi_mean, vi_mi_decomp, vi_hier_grouping, dim_MI_positions, dim_MI,
      SIMPLIFY = FALSE,
      FUN=function(mi_data, mi_first, mi_first_mapping, mi_mean, 
                   mi_decomp, mi_grouping, mi_dim_pos, mi_dim){
        mi_FS_mean <- lapply(mi_mean, FUN=function(i){FS(i,i)})
        M_mean <- mapply(mi_data, mi_FS_mean, SIMPLIFY = FALSE, FUN=function(i,j){i %*% j})
        # Sum by group to get the mean within "u" and "v"
        M_mean <- lapply(mi_grouping, FUN=function(g){
          # For each group, add first, then multiply (see below)
          Reduce('+', lapply(g, FUN=function(g_i){
            M_mean[[g_i]]
          }))
        })
        
        M_var <- mapply(mi_first, mi_dim_pos, mi_decomp, mi_first_mapping, SIMPLIFY = FALSE,
               FUN=function(first_l, pos_l, decomp_l, mapping_l){
                 mapping_l %*% get_longvar_UF(d_mi = mi_dim, Z = first_l,
                                pos_dim = pos_l, var_decomp = decomp_l)
               })

        out_trace <- rowSums(Reduce('*', M_var))
        out_quad <- 
          rowSums( M_mean[[1]] * M_var[[2]] ) +
          rowSums( M_mean[[2]] * M_var[[1]] )
      })
    out <- Reduce("+", out)
  }else{
    stop('...')
  }
  return(out)
}

bilinear_update_UF <- function(Z_wide, Z_MI, Z_hier_mapping,
                               d_mi,
                               diag_vi_pg_mean, it,
                               s, vi_pg_mean, offset,
                               prior_mi, Rmean, Rvar){

  aug_Z <- FS(Z_wide, Rmean)  
  list_aug_prior <- (
    mapply(Z_MI, prior_mi, SIMPLIFY = FALSE,
           FUN=function(i, p_i){kronecker(Diagonal(n=ncol(i)), matrix(p_i, d_mi))})
  )
  aug_prior <- bdiag(list_aug_prior)
  aug_Rvar <- t(Z_wide) %*% diag_vi_pg_mean %*% Rvar
  
  if (d_mi == 1){
    stop('WRONG RVAR')
    aug_Rvar <- drop0(Diagonal(x=aug_Rvar[,1]))
  }else{
    stop('WRONG RVAR')
    aug_Rvar <- bdiag(apply(aug_Rvar, MARGIN = 1, FUN=function(i){Matrix(i,d_mi)}))
  }
  
  aug_chol <- Cholesky(t(aug_Z) %*% Diagonal(x=vi_pg_mean) %*% aug_Z + aug_prior + aug_Rvar)
  aug_coef <- solve(
    aug_chol,
    t(aug_Z) %*% (s - vi_pg_mean * offset)
  )
  
  out_mean <- matrix(aug_coef, ncol = d_mi, byrow = T)
  split_out <- sapply(Z_MI, ncol)
  split_out <- c(0, cumsum(split_out))
  out_mean <- lapply(1:length(Z_MI), FUN=function(i){
    out_mean[seq(split_out[i] + 1, split_out[i+1]),,drop=FALSE]
  })
  
  vi_mi_decomp <- expand(aug_chol)
  vi_mi_decomp_nonpermute <- drop0(solve(vi_mi_decomp$L))
  vi_mi_decomp_LP <- t(vi_mi_decomp$P)
  log_det_mi_var <- -2 * sum(log(diag(vi_mi_decomp$L)))
  vi_mi_decomp <- vi_mi_decomp_nonpermute %*% t(vi_mi_decomp_LP)
  vi_mi_decomp <- drop0(vi_mi_decomp)
  
  if (!is.null(Z_hier_mapping)){
    
    if (it > 50){
      old_mean <- out_mean
      browser()
      stopifnot(all(sapply(prior_mi, FUN=function(i){isDiagonal(matrix(i, nrow = sqrt(length(i))))})))
      prior_mi <- lapply(prior_mi, FUN=function(i){matrix(i, nrow = sqrt(length(i)))})
      out_mean <- hier_px(Z_hier_mapping, out_mean, d_mi, prior_mi)
    }
    
  }
  
  out <- list(
    mean = out_mean,
    lndet = log_det_mi_var,
    decomp = vi_mi_decomp
  )
  return(out)
}

hier_px <- function(Z_hier_mapping, out_mean, d_mi, prior_mi){
  
  size_mean <- sapply(out_mean, nrow)
  
  wide_mapping <- do.call('cbind', Z_hier_mapping)
  flat_X <- rbind(wide_mapping, bdiag(lapply(out_mean[-1], FUN=function(i){Diagonal(n=nrow(i))})))
  level_flat <- sapply(out_mean, nrow)[-1]
  level_flat <- rep(names(Z_hier_mapping), level_flat)
  
  px_hier <- lapply(1:d_mi, FUN=function(d){
    prior_d <- sapply(prior_mi, FUN=function(i){i[d,d]})
    W <- Diagonal(x=rep(prior_d, size_mean))
    flat_aug <- do.call('c', mapply(out_mean[-1], SIMPLIFY = FALSE,
                FUN=function(j,p){j[,d]}))
    flat_y <- c(out_mean[[1]][,d], -1.0 * flat_aug)
    flat_coef <- as.vector(solve(Cholesky(t(flat_X) %*% W %*% flat_X), t(flat_X) %*% W %*% flat_y))
    return(split(flat_coef, level_flat))
  })
  
  px_hier <- lapply(names(Z_hier_mapping), FUN=function(i){
    sapply(px_hier, `[[`, i)
  })
  
  names(px_hier) <- names(Z_hier_mapping)
  print(px_hier)
  print('Before')
  print(t(wide_mapping) %*% out_mean[[1]])
  print(out_mean[[2]])
  
  # print(lapply(px_hier, FUN=function(i){colMeans(abs(i))}))
  out_mean[[1]] <- out_mean[[1]] - as.matrix(Reduce('+', mapply(Z_hier_mapping, px_hier, SIMPLIFY = FALSE, FUN=function(i,j){
    i %*% j
  })))
  
  out_mean[-1] <- mapply(out_mean[-1], px_hier, SIMPLIFY = FALSE, FUN=function(i,j){
    i + j
  })
  print('After')
  print(t(wide_mapping) %*% out_mean[[1]])
  print(out_mean[[2]])

  return(out_mean) 
}