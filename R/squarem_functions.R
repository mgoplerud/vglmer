
prep_lu <- function(M){
 fact_lu <- expand(Matrix::lu(M))
 return(fact_lu)
}
unprep_lu <- function(M){
  recons_M <- t(M$P) %*% M$L %*% M$U %*% M$Q
  # recons_M <- t(M$P) %*% drop0(zapsmall(M$L %*% M$U, 15)) %*% M$Q
  logdet_M <- 2 * sum(log(abs(diag(M$U))))
  return(list(M = recons_M, logdet_M = logdet_M))
}

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
  }else if (type == 'lu'){
    prep_lu(x)
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
  }else if (type == 'lu'){
    unprep_lu(x)
  }else if (type == 'cholesky'){
    unprep_cholesky(x)
  }else if (type == 'matrix'){
    unprep_matrix(x)
  }else if (type == 'positive'){
    unprep_positive(x)
  }else{stop('Invalid type')}
}
