if (FALSE){
#Do speed tests and verify that exactly the same

parsed_RE_groups <- get_RE_groups(formula = formula, data = data)

mapping_new_Z <- do.call('cbind', parsed_RE_groups$design)
mapping_J <- split(1:sum(d_j^2), rep(1:length(d_j), d_j^2))
mapping_J <- lapply(mapping_J, FUN=function(i){i-1})
mapping_J <- sapply(mapping_J, min)

mapping_to_re <- parsed_RE_groups$factor
mapping_to_re <- array_branch(do.call('cbind', mapping_to_re), margin = 1)
mapping_to_re <- lapply(mapping_to_re, FUN=function(i){
  mapply(outer_alpha_RE_positions, i, SIMPLIFY = FALSE, FUN=function(a,b){a[[b]]})
})
Mmap <- t(sapply(mapping_to_re, FUN=function(i){as.integer(sapply(i, min))}))


J <- length(d_j)

start_base_Z <- cumsum(c(0,d_j))[-(J+1)]
names(start_base_Z) <- NULL

base_Rvec_ridge <- function(vi_alpha_decomp, d_j, Z, diag_vi_pg_mean, outer_alpha_RE_positions){
  J <- length(d_j)
  d_sq <- d_j^2
  list_j <- as.list(rep(NA, J))
  for (j in 1:J){
    pos_j <- outer_alpha_RE_positions[[j]]
    store_j <- array(0, dim = rep(d_sq[j],2))
    for (g in pos_j){
      store_j <- store_j + kronecker(crossprod(vi_alpha_decomp[,g]), t(Z[,g]) %*% diag_vi_pg_mean %*% Z[,g])
    }
    list_j[[j]] <- store_j
  }
  return(bdiag(list_j))
}

base_Rvec_design <- function(vi_alpha_mean, d_j){
  J <- length(d_j)
  parseZ <- array_branch(mapping_new_Z, margin = 1)  
  lapply(parseZ, FUN=function(zi){
    for (j in 1:J){
      start_base_Z[j] + 1:d_j[j]
    }
  })
  
}
vecR_ridge_general(L = vi_alpha_decomp,
                   pg_mean = vi_pg_mean,
                   mapping_to_re = mapping_to_re, 
                   mapping_z = mapping_base_z,
                   mapping_J = mapping_J,
                   d = d_j)


bR <- base_Rvec_ridge(vi_alpha_decomp = vi_alpha_decomp, 
                      outer_alpha_RE_positions = outer_alpha_RE_positions,
                      d_j = d_j, Z = Z, diag_vi_pg_mean =diag_vi_pg_mean)

sim_R <- lapply(d_j, FUN=function(i){matrix(rnorm(i^2), ncol =i)})
sim_Rj <- mapply(sim_R, g_j, SIMPLIFY = FALSE, FUN=function(i,g){kronecker(Diagonal(n = g), i)})
sim_Rall <- bdiag(sim_Rj)
vec_Rall <- as.vector(sim_Rall)
vec_R <- unlist(lapply(sim_R, as.vector), use.names = F)

direct_verify <- function(vec_Rall, vi_alpha_decomp, Z, diag_vi_pg_mean){
  t(vec_Rall) %*% kronecker(t(vi_alpha_decomp) %*% vi_alpha_decomp, t(Z) %*% diag_vi_pg_mean %*% Z) %*% (vec_Rall)
}
cpp_verify <- function(vec_R, vi_alpha_decomp, vi_pg_mean, mapping_to_re, mapping_z, mapping_J, d_j){
  R <-  vecR_ridge_general(
    L = vi_alpha_decomp,
    Z = mapping_new_Z,
    pg_mean = vi_pg_mean,   
    M = Mmap,
    mapping_J = mapping_J, start_z = start_base_Z,
    d = d_j)
  return(t(vec_R) %*% R %*% vec_R)
}


dR <- vecR_design(alpha_mu = as.vector(vi_alpha_mean), Z = mapping_new_Z, M = Mmap, mapping_J = mapping_J, d = d_j,
                  start_z = start_base_Z)

t(s) %*% (Z %*% sim_Rall %*% vi_alpha_mean)
t(s) %*% dR %*% vec_R



}
