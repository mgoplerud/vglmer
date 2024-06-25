
FS <- function(X,Y){t(KhatriRao(t(X),t(Y)))}

get_bilinear_mean <- function(Z_MI, vi_mi_mean){
  out <- Reduce("+", mapply(Z_MI, vi_mi_mean, SIMPLIFY = FALSE, FUN=function(mi_data, mi_mean){
    rowSums(Reduce("*", mapply(mi_data, mi_mean, SIMPLIFY = FALSE, FUN=function(i,j){
      i %*% j
    })))
  }))
  return(out)
}

get_bilinear_var <- function(Z_MI, vi_mi_mean, vi_mi_var){
  out <- mapply(Z_MI, vi_mi_mean, vi_mi_var, SIMPLIFY = FALSE, 
    FUN=function(mi_data, mi_mean, mi_var){
      M_mean <- mapply(mi_data, mi_mean, SIMPLIFY = FALSE, FUN=function(i,j){i %*% j})
      M_var <- mapply(mi_data, mi_var, SIMPLIFY = FALSE, FUN=function(i,j){i %*% j})
      out_trace <- rowSums(Reduce('*', M_var))
      if (length(mi_data) != 2){stop('..set up for more than two...')}
      out_quad <- 
        rowSums( FS(M_mean[[1]], M_mean[[1]]) * M_var[[2]] ) +
        rowSums( FS(M_mean[[2]], M_mean[[2]]) * M_var[[1]] )
      return(out_trace + out_quad)
  })
  return(Reduce("+", out))
}

get_bilinear_outer <- function(vi_mi_mean, vi_mi_var, vi_mi_diag){
  mapply(
    vi_mi_mean, vi_mi_var, vi_mi_diag, 
    SIMPLIFY = FALSE,
    FUN=function(m_mean, m_var, m_diag){
      out_mean <- Reduce("+", lapply(m_mean, crossprod))
      out_var <- Reduce("+", lapply(m_var, colSums))
      out_var <- matrix(out_var, nrow = length(m_diag))
      return(out_mean + out_var)
    }
  )
}

#' @importFrom RSpectra svds
init_MI_from_svd <- function(data_mi, weight, rank){
  
  # t(Z1) %*% Diag(w) %*% Z2
  f <- function(x, args){
    as.vector(t(args$Z1) %*% Diagonal(x=args$w) %*% (args$Z2 %*% x))
  }
  g <- function(x, args){
    as.vector(t(args$Z2) %*% Diagonal(x=args$w) %*% (args$Z1 %*% x))
  }
  
  minimum_size <- min(sapply(data_mi, ncol))
  if (rank >= minimum_size){
    warning('rank of a multiplicative interaction is above minimum size...')
    svd_weight <- svd(t(data_mi[[1]]) %*% Diagonal(x=weight) %*% data_mi[[2]])
    # Pad with zeros
    svd_weight$d <- c(svd_weight$d, rep(0, rank - minimum_size))
    svd_weight$u <- cbind(svd_weight$u, matrix(0, nrow = ncol(data_mi[[1]]), ncol = rank - minimum_size))
    svd_weight$v <- cbind(svd_weight$v, matrix(0, nrow = ncol(data_mi[[2]]), ncol = rank - minimum_size))
  }else{
    svd_weight <- RSpectra::svds(A = f, k = rank, Atrans = g, 
                                 dim = sapply(data_mi, ncol),
                                 args = list(Z1 = data_mi[[1]], Z2 = data_mi[[2]], weight = weight))
  }

  weight_d <- Diagonal(x=sqrt(svd_weight$d))
  svd_weight$u <- svd_weight$u  %*% weight_d
  svd_weight$v <- svd_weight$v  %*% weight_d
  return(svd_weight)
}