
prep_cholesky <- function(L){
  diag(L) <- log(diag(L))
  return(L)
}
unprep_cholesky <- function(L){
  diag(L) <- exp(diag(L))
  return(L)
}
prep_matrix <- function(M){drop0(chol(as.matrix(M)))}
unprep_matrix <- function(M){t(M) %*% M}

prep_positive <- function(x){log(x)}
unprep_positive <- function(x){exp(x)}

squarem_prep_function <- function(x, type){
  if (type == 'real'){
    x
  }else if (type == 'cholesky'){
    prep_cholesky(x)
  }else if (type == 'matrix'){
    prep_matrix(x)
  }else if (type == 'positive'){
    prep_positive(x)
  }else{stop('Invalid type')}
}

squarem_unprep_function <- function(x, type){
  if (type == 'real'){
    x
  }else if (type == 'cholesky'){
    unprep_cholesky(x)
  }else if (type == 'matrix'){
    unprep_matrix(x)
  }else if (type == 'positive'){
    unprep_positive(x)
  }else{stop('Invalid type')}
}
