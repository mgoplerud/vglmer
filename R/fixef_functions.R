
# Face-Splitting Product or Row-Tensor Product
FS <- function(X,Z){t(KhatriRao(t(X), t(Z)))}

#' @export
v_fe <- function(group, interactions = ~ 1){
  # Using mgcv's syntax for "s" to make it work with "interpret.gam"
  vars <- as.list(substitute(list(group)))[-1]
  d <- length(vars)
  if (d > 1){stop('Provide a single grouping variable in each v_fe(...)')}
  
  fmt_interactions <- as.list(substitute(list(interactions)))[-1]
  
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  term[1] <- attr(terms(reformulate(term[1])), "term.labels")
  if (any(term == '.')){stop('v_fe(.) not supported')}
  
  fmt_interactions <- deparse(fmt_interactions[[1]], backtick = TRUE, width.cutoff = 500)
  fake_by <- paste(attr(terms(as.formula(fmt_interactions)), 'term.labels'), collapse = ' + ')
  if (fake_by == ""){
    fake_by <- "NA"
  }
  ret <- list(term = term, type = 'fe',
              interactions = fmt_interactions, by = fake_by, by_re = FALSE)
  class(ret) <- 'vglmer_special'
  
  return(ret)
}

#' Create a sparse design matrix for the fixed effects
vglmer_build_fe <- function(group, inter, levels = NULL){
  
  if (is.null(levels)){
    if (!is.null(levels)){
      ux <- levels(group)
      if (any(!(group %in% ux))){stop('Not all values found in factor levels.')}
    }else{
      ux <- sort(unique(group))
    }
  }else{
    ux <- levels
  }
  if (is.null(inter)){
    x <- Matrix(rep(1, length(group)))
  }else{
    x <- Matrix(as.matrix(inter))
  }
  lookup <- match(group, ux)
  lookup <- sparseMatrix(i = seq_len(nrow(x)), j = lookup, x = 1)
  
  out <- list(x = x, lookup = lookup, fe_size = ncol(x),
              attr = list(unique_values = ux, constrast = fast_contr_sum(length(ux))))
  class(out) <- c('fe_sparse')
  return(out)
}

# Sum up mean and variance of FE
calculate_FE <- function(X, Z, FS_XX, mean, var){

  calc_FE_all <- mapply(X, Z, FS_XX, mean, var, SIMPLIFY = FALSE, 
         FUN=function(X_i, Z_i, FS_XX_i, mean_i, var_i){
           est_mean <- rowSums(X_i * (Z_i %*% mean_i))
           est_var <- rowSums(FS_XX_i * (Z_i %*% var_i))  
           return(cbind(est_mean, est_var))  
  })
  
  all_mean <- rowSums(sapply(calc_FE_all, FUN=function(i){i[,1]}))
  all_var <- rowSums(sapply(calc_FE_all, FUN=function(i){i[,2]}))
  return(cbind('mean' = all_mean, 'var' = all_var))
}

# Function for updating FE with sum-to-zero constraints
fast_generic_FE <- function(FS_XX, X, Z, y, weights, dim_fe, levels_fe){
  
  lhs_meat <- as.matrix(t(Z) %*% Diagonal(x = weights) %*% FS_XX)
  inv_meat <- invert_rowwise(lhs_meat, dim_fe)
  lhs_meat <- inv_meat$inverse
  det_lhs_meat <- inv_meat$det
  rhs_meat <- t(Z) %*% Diagonal(x = y) %*% X
  raw_OLS <- matrix(NA, nrow = levels_fe, ncol = dim_fe)
  weight_matrix <- solve(matrix(colSums(lhs_meat), ncol = dim_fe))
  for (d in 1:dim_fe){
    idx_M2 <- 1 + dim_fe * (d-1) + (0:(dim_fe-1))
    select_rows <- lhs_meat[, idx_M2]
    raw_OLS[, d] <- rowSums(rhs_meat * select_rows)
  }
  adjust_restrict <- lhs_meat %*% kronecker(Diagonal(n = dim_fe), weight_matrix %*% colSums(raw_OLS))
  
  mean_FE <- raw_OLS - adjust_restrict
  
  reshape_lhs_meat <- matrix(as.vector(t(lhs_meat)), ncol = dim_fe, byrow = T)
  var_FE <- lhs_meat - t(matrix(
    colSums(matrix(as.vector(
      FS(kronecker(Diagonal(n = levels_fe), weight_matrix) %*% reshape_lhs_meat, reshape_lhs_meat)), ncol = levels_fe * dim_fe^2)), 
    ncol = levels_fe, byrow = T)
  )  
  
  det_block <- apply(lhs_meat, MARGIN = 1, FUN=function(i){
    as.numeric(determinant(matrix(i, ncol = dim_fe))$modulus)
  })
  # https://arxiv.org/pdf/2207.08038.pdf
  # lndet^*(M) = lndet*(C C^T) - lndet([FS_XZ^T weights FS_XZ]^{-1})  + lndet(sum(....))
  lndet_FE <- sum(det_block)  + as.numeric(determinant(weight_matrix)$modulus)
  # -(sum(log(diag_ZtZ)) + log(sum(1/diag_ZtZ))) [ignore]((+ log(ncol(Z)) ))
  
  return(list(mean = mean_FE, var = var_FE, lndet = lndet_FE, raw = lhs_meat))  
}

# Function for updating FE with sum-to-zero constraint, only 1-D
fast_1D_FE <- function(X, Z, weights, y){
  
  diag_ZtZ <- as.vector(t(Z) %*% Diagonal(x = weights) %*% X^2)
  inv_diag_ZtZ <- 1/diag_ZtZ
  
  OLS <- Diagonal(x = inv_diag_ZtZ) %*% (t(Z) %*% Diagonal(x = y) %*% X)
  var_FE <- inv_diag_ZtZ - inv_diag_ZtZ^2/sum(inv_diag_ZtZ)
  mean_FE <- OLS - inv_diag_ZtZ * sum(OLS)/sum(inv_diag_ZtZ)
  lndet_FE <- -(sum(log(diag_ZtZ)) + log(sum(1/diag_ZtZ)))
  
  return(list(var = matrix(var_FE), mean = Matrix(mean_FE), 
              lndet = lndet_FE, raw = inv_diag_ZtZ))
}

update_FE <- function(vi_mean, vi_var, vi_lndet, FE_data, FE_lookup, y, weights, 
                      dim_fe, levels_fe, FE_rowtens, family){

  if (family == 'poisson'){stop('set up update FE For poisson.')}
  running_FE <- calculate_FE(X = FE_data, Z = FE_lookup, FS_XX = FE_rowtens, mean = vi_mean, var = vi_var)
  
  n_FE <- length(FE_data)
  vi_var_raw <- as.list(rep(NA, n_FE))
  
  for (v in 1:n_FE){
    
    FE_rt_v <- FE_rowtens[[v]]
    FE_data_v <- FE_data[[v]]
    FE_lookup_v <- FE_lookup[[v]]
    
    init_mean_v <- rowSums(FE_data_v * (FE_lookup_v %*% vi_mean[[v]]))
    init_var_v <- rowSums(FE_rt_v * (FE_lookup_v %*% vi_var[[v]]))
    if (n_FE > 1){
      adj_FE <- running_FE[,1] - init_mean_v
    }else{
      adj_FE <- 0
    }
    
    # update_FE_v <- fast_1D_FE(X = FE_data_v, Z = FE_lookup_v, weights = weights, y = y - (running_FE[,1] - init_mean_v))
    update_FE_v <- fast_generic_FE(FS_XX = FE_rt_v, X = FE_data_v, Z = FE_lookup_v,
                    y = y - weights * adj_FE, weights = weights,
                    dim_fe = dim_fe[v], levels_fe = levels_fe[v])
    vi_mean[[v]] <- update_FE_v$mean
    vi_var[[v]] <- update_FE_v$var  
    vi_lndet[v] <- update_FE_v$lndet
    vi_var_raw[[v]] <- update_FE_v$raw
    # running_FE <- calculate_FE(X = FE_data, Z = FE_lookup, FS_XX = FE_rowtens, mean = vi_mean, var = vi_var)
    new_mean_v <- rowSums(FE_data_v * (FE_lookup_v %*% vi_mean[[v]]))
    new_var_v <- rowSums(FE_rt_v * (FE_lookup_v %*% vi_var[[v]]))

    running_FE[,1] <- running_FE[,1] + (new_mean_v - init_mean_v)
    running_FE[,2] <- running_FE[,2] + (new_var_v - init_var_v)
  }
  
  return(list(mean = vi_mean, var = vi_var, lndet = vi_lndet, raw = vi_var_raw))
  
}
