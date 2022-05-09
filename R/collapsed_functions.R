
build_collapse_index <- function(X, Z, cyclical_pos, weight, names_of_RE, k){
  
  mean_nonsparse_X <- colMeans(Diagonal(x = weight) %*% (X != 0))
  mean_nonsparse_Z <- colMeans(Diagonal(x = weight) %*% (Z != 0))
  
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
    if (k >= (length(mean_nonsparse_X) + length(mean_nonsparse_Z))){
      thresh <- 0
    }else{
      thresh <- rev(sort(c(mean_nonsparse_X, mean_nonsparse_Z)))[k]
    }
  }
  
  if (!is.null(thresh)){
    collapse_Z <- which(mean_nonsparse_Z >= thresh)
    collapse_X <- which(mean_nonsparse_X >= thresh)
  }
  
  C_j <- lapply(cyclical_pos, FUN=function(i){ncol(X) + base::intersect(i, collapse_Z)})
  M_j <- lapply(cyclical_pos, FUN=function(i){ncol(X) + base::setdiff(i, collapse_Z)})
  
  C_j <- c(list(collapse_X), C_j)
  M_j <- c(list(base::setdiff(seq_len(ncol(X)), collapse_X)), M_j)
  
  names(C_j) <- names(M_j) <- c('Fixed Effect', names(names_of_RE))
  print(lengths(C_j))
  
  return(
    list(C_j = C_j, 
         M_j = M_j)
  )
}


