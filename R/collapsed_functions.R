
FS <- function(X,Y){
  t(KhatriRao(t(X),t(Y)))
}

build_collapse_index <- function(X, Z, cyclical_pos, outer_alpha_RE_positions,
    weight, names_of_RE, k, d_j){

  mean_nonsparse_X <- colMeans(Diagonal(x = weight) %*% (X != 0))
  mean_nonsparse_Z <- colMeans(Diagonal(x = weight) %*% (Z != 0))
  # Get average of sparsity for each group "g" (relevant if d_j > 1)
  mean_nonsparse_Z <- lapply(outer_alpha_RE_positions, FUN=function(i){sapply(i, FUN=function(j){mean(mean_nonsparse_Z[j])})})
  
  if (is.character(k)){
    
    if (!all(k %in% c('FE', names(names_of_RE)))){stop('collapsed_set must be vector of named REs or number.')}
    
    if ('FE' %in% k){
      collapse_X <- 1:ncol(X)
    }else{
      collapse_X <- numeric(0)
    }
    collapse_Z <- unlist(cyclical_pos[which(names(names_of_RE) %in% k)])
    
    index_collapse_Z <- lapply(names(names_of_RE), FUN=function(i){
      if (i %in% k){
        as.numeric(names(outer_alpha_RE_positions[[i]]))
      }else{
        vector(mode = 'numeric')
      }
    })
    index_marginal_Z <- lapply(names(names_of_RE), FUN=function(i){
      if (!(i %in% k)){
        as.numeric(names(outer_alpha_RE_positions[[i]]))
      }else{
        vector(mode = 'numeric')
      }
    })
    thresh <- NULL
    
  }else if (!is.finite(k)){
    thresh <- 0
  }else if (k == 0){
    thresh <- Inf
  }else{
    if (k >= (length(mean_nonsparse_X) + sum(lengths(mean_nonsparse_Z)))){
      thresh <- 0
    }else{
      thresh <- rev(sort(c(mean_nonsparse_X, unlist(mean_nonsparse_Z))))[k]
      names(thresh) <- NULL
    }
  }

  if (!is.null(thresh)){
    collapse_Z <- mapply(outer_alpha_RE_positions, mean_nonsparse_Z, 
      SIMPLIFY = FALSE, FUN=function(i,j){
        i[which(j >= thresh)]
    })
    
    index_collapse_Z <- lapply(collapse_Z, FUN=function(i){as.numeric(names(i))})
    index_marginal_Z <- mapply(outer_alpha_RE_positions, index_collapse_Z, SIMPLIFY = FALSE, FUN=function(i,j){
      setdiff(as.numeric(names(i)), j)
    })
    collapse_Z <- unlist(collapse_Z)
    names(collapse_Z) <- NULL
    collapse_X <- which(mean_nonsparse_X >= thresh)
  }
  
  C_j <- lapply(cyclical_pos, 
    FUN=function(i){ncol(X) + base::intersect(i, collapse_Z)})
  M_j <- lapply(cyclical_pos, 
    FUN=function(i){ncol(X) + base::setdiff(i, collapse_Z)})
  
  C_j <- c(list(collapse_X), C_j)
  M_j <- c(list(base::setdiff(seq_len(ncol(X)), collapse_X)), M_j)
  
  names(C_j) <- names(M_j) <- c('Fixed Effect', names(names_of_RE))
  print(lengths(C_j)/c(1, d_j))
  
  names_of_collapsed <- list(FE = names(collapse_X), RE = names(collapse_Z))  
  return(
    list(C_j = C_j, 
         M_j = M_j, 
         index_collapse = index_collapse_Z,
         index_marginal = index_marginal_Z,
         names_of_collapsed = names_of_collapsed)
  )
}

block_diag_sum <- function(A, g, d){
  r <- 0:(d-1)
  out <- Reduce('+', lapply(0:(g-1), FUN=function(i){
    s <- 1 + d * i + r
    return(A[s,s])
  }))
  return(out)
}

prepare_total_variance <- function(
  vi_C_mean, vi_M_mean, vi_M_var_flat, aug_dj,
  vi_C_uncond,
  any_collapsed_M, any_collapsed_C,
  position_block_j,
  c_id_Mj, block_collapse){

  if (any_collapsed_M){

    ssq_M <- mapply(vi_M_mean, vi_M_var_flat, aug_dj, SIMPLIFY = FALSE,
      FUN=function(i,v,j){
        if (j == 0){
          return(NA)
        }else if (length(i) > 0){
          m <- crossprod(matrix(i,ncol=j, byrow = T))
          v <- matrix(colSums(v), ncol = j)
          return(m + v)
        }else{return(matrix(0, nrow = j, ncol = j))}
      })
    
  }else{
    ssq_M <- mapply(1:length(vi_M_mean), aug_dj, SIMPLIFY = FALSE, FUN=function(i,j){
      matrix(0, nrow = j, ncol = j)
    })
  }
  if (any_collapsed_C){
    
    ssq_C_mean <- mapply(split(vi_C_mean, position_block_j)[c_id_Mj], aug_dj, SIMPLIFY = FALSE,
       FUN=function(i,j){
         if (length(i) > 0){
           if (j != 0){
             m <- crossprod(matrix(i,ncol=j, byrow = T))
           }else{
             m <- NA
           }
           return(m)
         }else{return(matrix(0, nrow = j, ncol = j))}
       })
    ssq_C_var <- lapply(split(1:length(position_block_j), position_block_j), 
                        FUN=function(i){vi_C_uncond[i,i, drop=F]})[c_id_Mj]
    ssq_C_var <- mapply(ssq_C_var, aug_dj, FUN=function(i,j){
      if (length(i) > 0){
        if (j == 1){
          return(matrix(sum(diag(i))))
        }else if (j == 0){
          return(NA)
        }else{
          block_diag_sum(A = i, g = ncol(i)/j, d = j)
        }
      }else{return(matrix(0, nrow = j, ncol = j))}
    })
    
  }else{
    ssq_C_var <- ssq_C_mean <- mapply(1:length(vi_M_mean), aug_dj, SIMPLIFY = FALSE, FUN=function(i,j){
      matrix(0, nrow = j, ncol = j)
    })
  }
  
  if (block_collapse){
    stop('setup blockcollapse')
    block_terms <- sapply(split(diag(vi_M_var[[position_of_blocked]]) + 
                                  vi_M_mean[[position_of_blocked]]^2, position_block_j), sum)
    block_terms <- block_terms[as.character(0:number_of_RE)]
    block_terms[is.na(block_terms)] <- 0
    ssq_C <- c(NA, block_terms)
  }
  
  ssq_out <- mapply(ssq_C_var, ssq_C_mean, ssq_M, SIMPLIFY = FALSE, 
    FUN=function(i,j,k){
      if (is.null(k)){
        return(NA)
      }else if (all(is.na(k))){
        return(NA)
      }else if (ncol(k) == 0){
        return(NA) 
      }else{
        i+j+k
      }
    })
  
  return(ssq_out)
}


extract_precision <- function(object){
  
  if (!(object$family %in% c('binomial', 'linear'))){
    stop('Extracting precision not set up for this family.')
  }
  
  if (object$control$factorization_method == 'strong'){
    
    X <- object$data$X
    Z <- object$data$Z
    cyclical_pos <- object$data$cyclical_pos
    diag_vi_pg_mean <- object$data$diag_vi_pg_mean
    Tinv <- object$data$Tinv
    
    prec_Z <- bdiag(lapply(1:length(cyclical_pos), FUN=function(j){
      Z_j <- Z[, cyclical_pos[[j]], drop = F]
      return(t(Z_j) %*% diag_vi_pg_mean %*% Z_j + Tinv[[j]])
    }))
    prec_X <- t(X) %*% diag_vi_pg_mean %*% X
    vi_precision <- bdiag(prec_X, prec_Z)
    
    if (object$family == 'linear'){
      adjust_prec <- object$sigmasq$a/object$sigmasq$b
      vi_precision <- vi_precision * adjust_prec
    }
    
  }else if (object$control$factorization_method == 'weak'){
    
    X <- object$data$X
    Z <- object$data$Z
    diag_vi_pg_mean <- object$data$diag_vi_pg_mean
    
    joint.XZ <- cbind(X,Z)
    vi_precision <- t(joint.XZ) %*% diag_vi_pg_mean %*% joint.XZ + 
      bdiag(Diagonal(x = rep(0, ncol(X))), object$data$Tinv)
    
    if (object$family == 'linear'){
      adjust_prec <- object$sigmasq$a/object$sigmasq$b
      vi_precision <- vi_precision * adjust_prec
    }
    
  }else if (object$control$factorization_method == 'partially_factorized'){
    
    design_C <- object$data$C_design
    vi_M_list <- object$data$M_design
    M_j <- object$marginal$index
    diag_vi_pg_mean <- object$data$diag_vi_pg_mean
    
    Tinv_M <- object$data$Tinv
    Tinv_M <- lapply(M_j, FUN=function(i){Tinv_M[i,i,drop=F]})
    
    vi_C_var <- object$collapsed$cond_variance
    vi_C_uncond <- object$collapsed$variance
    vi_B_raw <- object$data$vi_B_raw
    vi_P <- object$data$vi_P
    aug_dj <- object$data$aug_dj

    # Get t(P) %*% bdiag(Var(alpha_j)) %*% P
    add_adjust <- Reduce('+', mapply(vi_B_raw, vi_P, aug_dj, SIMPLIFY = FALSE,
     FUN=function(b,p,d){
       b <- as.vector(t(b))
       if (length(b) > 0){
         matrix(b, ncol = ncol(p), byrow = TRUE) %*% t(p)
       }else{return(0)}
     }))
    
    if (object$family == 'linear'){
      adjust_var <- 1/(object$sigmasq$a/object$sigmasq$b)
      add_adjust <- adjust_var * add_adjust
    }
    
    if (is.null(vi_C_uncond)){
      schur_C <- matrix(nrow = 0, ncol = 0)
      vi_C_var <- matrix(nrow = 0, ncol = 0)
    }else{
      schur_C <- solve(vi_C_uncond - add_adjust)
    }
    
    if (object$family == 'linear'){
      adjust_prec <- object$sigmasq$a/object$sigmasq$b
    }else{
      adjust_prec <- 1
    }
    
    agg_P <- do.call('cbind', vi_P)
    
    vi_prec_M <- bdiag( 
      mapply(vi_M_list, vi_P, Tinv_M, SIMPLIFY = FALSE, FUN=function(M_j, P_j, Tinv_j){
        
        int_term <- t(M_j) %*% diag_vi_pg_mean %*% design_C
        
        out <- t(M_j) %*% diag_vi_pg_mean %*% M_j + 
          - int_term %*% (adjust_prec * vi_C_var) %*% t(int_term) +
          Tinv_j
        return(out)
      })
    )
    vi_prec_M <- adjust_prec * vi_prec_M
    
    vi_precision <- rbind(
      cbind(schur_C, schur_C %*% agg_P),
      cbind(t(agg_P) %*% schur_C, vi_prec_M + t(agg_P) %*% schur_C %*% agg_P)
    )
    # Permute to align with standard ordering
    vi_precision <- vi_precision[object$marginal$reverse_collapsing,object$marginal$reverse_collapsing]
    
  }else{stop('Extracting precision not set up.')}
  
  return(vi_precision)
}