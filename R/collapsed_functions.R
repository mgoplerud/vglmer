
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
    
    if (!all(k %in% c('FE', names(names_of_RE)))){stop('collapsed_size must be vector of named REs or number.')}
    
    if ('FE' %in% k){
      collapse_X <- 1:ncol(X)
    }else{
      collapse_X <- numeric(0)
    }
    collapse_Z <- unlist(cyclical_pos[which(names(names_of_RE) %in% k)])
    
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
         names_of_collapsed = names_of_collapsed)
  )
}


