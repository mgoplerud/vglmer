
FS <- function(X,Y){t(KhatriRao(t(X),t(Y)))}

get_bilinear_mean <- function(Z_MI, vi_mi_mean, vi_hier_grouping, is_onehot, reduce = TRUE){
  
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

get_bilinear_var <- function(Z_MI, vi_mi_mean, vi_mi_var,
  vi_hier_grouping, is_onehot, reduce = TRUE){
  
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
    out <- mapply(Z_MI, vi_mi_mean, vi_mi_var, vi_hier_grouping, is_onehot, SIMPLIFY = FALSE, 
        FUN=function(mi_data, mi_mean, mi_var, mi_grouping, mi_onehot){
          simple_model <- all(lengths(mi_grouping) == 1) & all(mi_onehot)
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
            M_var <- mapply(mi_data, mi_var, mi_onehot, SIMPLIFY = FALSE, FUN=function(i,j,o){
              if (o){
                i %*% j
              }else{
                i^2 %*% j
              }
            })
            M_var <- lapply(mi_grouping, FUN=function(g){
              # For each group, add together the variance
              Reduce('+', lapply(g, FUN=function(g_i){
                M_var[[g_i]]
              }))
            })
          }else{
            mi_FS_mean <- lapply(mi_mean, FUN=function(i){FS(i,i)})
            M_mean <- mapply(mi_data, mi_FS_mean, SIMPLIFY = FALSE, FUN=function(i,j){i %*% j})
            M_var <- mapply(mi_data, mi_var, mi_onehot, SIMPLIFY = FALSE, FUN=function(i,j,o){
              if (o){i %*% j}else{
                if (ncol(j) > 1){browser()}
                i^2 %*% j
              }
            })
          }
          out_trace <- rowSums(Reduce('*', M_var))
          out_quad <- 
            rowSums( M_mean[[1]] * M_var[[2]] ) +
            rowSums( M_mean[[2]] * M_var[[1]] )
          return(out_trace + out_quad)
        
      })
    return(Reduce("+", out))
  }else{
    out <- mapply(Z_MI, vi_mi_mean, vi_mi_var, is_onehot, SIMPLIFY = FALSE, 
    FUN=function(mi_data, mi_mean, mi_var, mi_onehot){
      M_var <- mapply(mi_data, mi_var, mi_onehot,
        SIMPLIFY = FALSE, FUN=function(i,j,o){
          if (o){
            i %*% j
          }else{
            if (ncol(j) > 1){browser()}
            return(i^2 %*% j)
          }
        })
      return(M_var)
    })
    return(out) 
  }
}

get_bilinear_outer <- function(vi_mi_mean, vi_mi_var, vi_mi_diag, mi_prior_type,
                               is_onehot, reduce = TRUE){

  mapply(
    vi_mi_mean, vi_mi_var, vi_mi_diag, is_onehot,
    SIMPLIFY = FALSE,
    FUN=function(m_mean, m_var, m_diag, m_hier, m_onehot){
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

bilinear_update <- function(Z_MI, prior_mi, d_mi, s, is_onehot,
                            vi_pg_mean, diag_vi_pg_mean, RFSmean, Rvar, Rmean, offset,
                            return_chol, method, it,
                            auxiliary_Z = NULL, auxiliary_prior = NULL){
  
  simple_j <- !inherits(Z_MI, 'list') & all(is_onehot)

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
      aug_Rvar <- sparseMatrix(i = flat_tile[,1], j = flat_tile[,2], x = flat_tile[,3],
                               dims = rep(ncol(full_Z) * d_mi,2))
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
      aug_Z <- cbind(auxiliary_Z, aug_Z)
      aug_prior <- bdiag(auxiliary_prior, aug_prior)
      aug_Rvar <- bdiag(Diagonal(x = rep(0, ncol(auxiliary_Z))), aug_Rvar)
    }
    
    prec_aug <- t(aug_Z) %*% Diagonal(x=vi_pg_mean) %*% aug_Z + aug_prior + aug_Rvar
    aug_chol <- Cholesky(prec_aug)
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
    
    old_mean <- out_mean
    prior_mi <- lapply(prior_mi, FUN=function(i){matrix(i, nrow = sqrt(length(i)))})
    
    ub <- mapply(Z_MI, prior_mi, out_mean, is_onehot[names(Z_MI)], SIMPLIFY = FALSE, FUN=function(Z_j, prior_j, om_j, onehot_j){
      if (onehot_j){
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
      }else{
        if (d_mi > 1){
          
          Z_j <- as.matrix(Z_j)
          prec_j <- sapply(1:d_mi^2, FUN=function(d){
            colSums(Z_j * (Diagonal(x=as.vector(diag_vi_pg_mean %*% (RFSmean[,d] + Rvar[,d]))) %*% Z_j))
            # colSums((Diagonal(x=sqrt(as.vector(diag_vi_pg_mean %*% (RFSmean[,d] + Rvar[,d])))) %*% Z_j)^2)
            # diag(t(Z_j) %*%
            #        Diagonal(x=as.vector(diag_vi_pg_mean %*% (RFSmean[,d] + Rvar[,d]))) %*%
            #        Z_j)
          })
          if (ncol(Z_j) == 1){
            prec_j <- t(prec_j)
          }
          prec_j <- as.matrix(prec_j)
          ub_j <- invert_rowwise(prec_j, vec_prior = matrix(as.vector(prior_j)), 
                         RHS = matrix(rep(0, nrow(prec_j))),
                         dim = d_mi, return_chol = return_chol)
          ub_j$mean <- om_j
        }else{
          if (it == 1){warning('inefficient')}
          # prec_j <- diag(t(Z_j) %*% Diagonal(x=as.vector(diag_vi_pg_mean %*% (RFSmean + Rvar))) %*% Z_j)
          prec_j <- colSums((Diagonal(x=sqrt(as.vector(diag_vi_pg_mean %*% (RFSmean + Rvar)))) %*% Z_j)^2)
          var_j <- matrix(1/(prec_j + as.vector(prior_j)))
          ub_j <- list(inverse = var_j, det = 0.5 * log(var_j), mean = om_j)
        }
      }
      
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
    
    # print(dim(update_bilinear$mean))
    # print('var')
    # print(rowMeans(apply(update_bilinear$inverse, MARGIN = 1, FUN=function(i){diag(matrix(i, ncol = d_mi))})))
    # print('mean')
    # print(colMeans(update_bilinear$mean))
    # print('eigen')
    # print(colMeans(t(apply(as.matrix(
    #   t(Z_MI) %*% diag_vi_pg_mean %*%
    #     (RFSmean + Rvar)
    # ), MARGIN = 1, FUN=function(i){
    #   eigen(solve(matrix(i + prior_mi, 4)))$values
    # }))))
    
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

initalize_mi_mean <- function(Z_j, Z_j_onehot, j_group, j_hier, j_attr,
                              j_nesting, mi_d_j, init_method,
                              init_y = NULL, init_w = NULL){
  if (j_hier){
    if (init_method == 'random'){
      principal_group <- sapply(j_group, FUN=function(i){i[1]})
      out <- mapply(Z_j[principal_group], Z_j_onehot[principal_group], 
                    SIMPLIFY = FALSE, FUN=function(i, oh){
        if (oh){
          sd_i <- 1
        }else{
          sd_i <- sqrt(colMeans(i^2) - colMeans(i)^2)
        }
        out_i <- matrix(rnorm(mi_d_j * ncol(i), 
                     sd = 1/rep(sd_i/mi_d_j, mi_d_j)),
               ncol = mi_d_j)
        return(out_i)
      })
      names(out) <- c('u', 'v')
      out_var <- NULL
    }else{
      principal_group <- sapply(j_group, FUN=function(i){i[1]})
      out <- init_MI_from_svd(
        data_mi = Z_j[principal_group],
        rank = mi_d_j,
        y = init_y, pg_weight = init_w,
        prior_U = 1, prior_V = 1
      )
      out_var <- out[c('var_U', 'var_V')]
      out <- out[c('mean_U', 'mean_V')]
      names(out) <- c('u', 'v')
    }
    out_data <- lapply(1:2, FUN=function(k){
      if (k == 1){
        init_principal_k <- out$u
      }else{
        init_principal_k <- out$v
      }
      principal_k <- j_group[[k]][1]
      other_levels <- j_group[[k]][-1]
      levels_k <- j_attr$levels[j_group[[k]]]
      if (length(other_levels) > 0){
        
        # Estimate lambda needed to calibrate edf at around 90% of maximum
        # (light but some stabilization)
        block_prior <- sapply(Z_j[other_levels], FUN=function(i){
          iti <- crossprod(i)
          if (isDiagonal(iti)){
            diag_iti <- diag(iti)
            out_cl <- optim(par = 0, fn = get_edf_uni,
                            diag_z = diag_iti, target = ncol(iti) * 0.90,
                            method = 'L-BFGS-B')
            return(exp(out_cl$par))
          }else{
            
            f <- function(alpha, M, x, target){
              return((get_edf_diag(alpha, M, x) - target)^2)
            }
            out_cl <- optim(par = 0, 
                  fn = f, M = iti, x = ncol(iti),
                  target = ncol(iti) * 0.90,
                  method = 'L-BFGS-B')
            return(exp(out_cl$par))
          }
        })        
        n_other <- sapply(Z_j[other_levels], ncol)
        block_prior <- Diagonal(x=rep(block_prior, n_other))
        block_Z <- do.call('cbind', Z_j[other_levels])

        out_decompose <- lapply(1:ncol(init_principal_k), FUN=function(d){
          # Take the initial estimates from, e.g., Nakajima or random
          # and assign to each observed for the corresponding principal group
          initial_estimate <- as.vector(Z_j[[principal_k]] %*% init_principal_k[,d])
          # Estimate a regression of initial ~ others
          block_estimates <- 
            solve(Cholesky(crossprod(block_Z) + block_prior), 
                  t(block_Z) %*% initial_estimate)
          # Get the average residual for each group in the principal mapping
          principal_estimate <- 
            as.vector(
              solve(Cholesky(crossprod(Z_j[[principal_k]])),
                    t(Z_j[[principal_k]]) %*%
                      (initial_estimate - block_Z %*% block_estimates)
              ))
          block_estimates <- split(block_estimates, rep(1:length(n_other), n_other))
          return(c(list(principal_estimate), block_estimates))
        })
        
        out_decompose <- lapply(1:length(j_group[[k]]), FUN=function(m){
          out_iii <- do.call('cbind', lapply(out_decompose, FUN=function(iii){iii[[m]]}))
          rownames(out_iii) <- levels_k[[m]]
          colnames(out_iii) <- paste0('dim', 1:ncol(init_principal_k))
          return(out_iii)
        })
        names(out_decompose) <- j_group[[k]]
        return(out_decompose)
      }else{
        init_principal_k <- as.matrix(init_principal_k)
        rownames(init_principal_k) <- levels_k[[principal_k]]
        colnames(init_principal_k) <- paste0('dim', 1:ncol(init_principal_k))
        out_decompose <- setNames(list(init_principal_k), principal_k)
      }
      return(out_decompose)
    })
    out_data <- unlist(out_data, recursive = FALSE)
    out_data <- out_data[names(Z_j)]
    if (is.null(out_var)){
      return(list(mean = out_data))
    }else{
      names(out_var) <- principal_group
      return(list(mean = out_data, var = out_var))
    }
  }else{
    if (length(Z_j) != 2){stop('...')}
    if (control$mi_init == 'random'){
      out <- sapply(j_group, FUN=function(i){i[1]})
      out <- lapply(Z_MI[[j]][out], FUN=function(i){
        matrix(rnorm(mi_d_j[j] * ncol(i)), ncol = mi_d_j[j])
      })
      out <- list(mean = out)
    }else{
      out <- init_MI_from_svd(
        data_mi = Z_j,
        rank = mi_d_j[j], y = init_y,
        pg_weight = init_w,
        prior_U = 1, prior_V = 1)
      out <- list(mean = list(out$mean_U, out$mean_V),
                  var = list(out$var_U, out$var_V))
    }
    return(out)
  }
}

expect_mi_prior_kernel <- function(
    vi_mi_sigma_alpha, vi_mi_sigma_alpha_nu, 
    vi_mi_sigma_outer_alpha, 
    mi_d_j, mi_prior_type){
  
  
  mi_moments_sigma_alpha <- mapply(
    vi_mi_sigma_alpha, vi_mi_sigma_alpha_nu, mi_d_j,
    SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
      if (mi_prior_type %in% c('centered', 'shared')){
        inv_phi <- solve(phi)
        sigma.inv <- nu * inv_phi
        ln.det <- log(det(phi)) - sum(digamma((nu - 1:d + 1) / 2)) - d * log(2)
        return(list(sigma.inv = sigma.inv, ln.det = ln.det))
      }else{
        mapply(phi, nu, d, SIMPLIFY = FALSE, FUN=function(phi_l, nu_l, d_l){
          inv_phi_l <- solve(phi_l)
          sigma.inv_l <- nu_l * inv_phi_l
          ln.det_l <- log(det(phi_l)) - sum(digamma((nu_l - 1:d_l + 1) / 2)) - d_l * log(2)
          return(list(sigma.inv = sigma.inv_l, ln.det = ln.det_l))
        })
      }
    })

  mi_inv_sigma_alpha <- lapply(mi_moments_sigma_alpha, FUN = function(i) {
    if (mi_prior_type %in% c('centered', 'shared')){
      i$sigma.inv
    }else{
      lapply(i, `[[`, 'sigma.inv')
    }
  })
  
  if (mi_prior_type %in% c('centered', 'shared')){
    
    out <- sum(
      mapply(mi_inv_sigma_alpha, vi_mi_sigma_outer_alpha, FUN = function(a, b) {
        sum(diag(a %*% b))
    }))

  }else if (mi_prior_type %in% c('separate')){
    
    out <- sum(mapply(
      mi_inv_sigma_alpha, vi_mi_sigma_outer_alpha,
      FUN=function(invsig_l, oa_l){
        sum(mapply(invsig_l, oa_l, FUN = function(a, b) {
            sum(diag(a %*% b))
          }))
      }))
  }else{stop('set up kernel for mi_prior_type')}

  return(out)
}