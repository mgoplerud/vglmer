context("Translation Expansion Tests")

N <- 1000
G <- 50
x <- rnorm(N)
x2 <- rnorm(N)
g <- sample(1:G, N, replace = T)
g2 <- sample(1:G, N, replace = T)
alpha <- rnorm(G)

y <- rbinom(n = N, size = 1, prob = plogis(-1 + x + alpha[g]))

est_vglmer <- vglmer(y ~ x + x2 + (1 + x | g) + (1 + x2 | g2),
                    data = NULL, 
                    control = vglmer_control(
                      parameter_expansion = 'translation',
                      factorization_method = 'partial'),
                    family = 'binomial')


est_vglmer_kn

sum(
  sapply(purrr::array_branch(combn(1:number_of_RE, 2), margin = 2),
  FUN=function(i){
    unique(Mmap[,c(i)]) %>% nrow
  })
)
J <- number_of_RE

row_combinations <- list()
for (j in 1:J){
  cat('j')
  combo_j <- list()
  for (jprime in j:J){
    combinations <- array_branch(unique(Mmap[,c(j,jprime)]), margin = 1)
    combinations <- lapply(combinations, FUN=function(i){
      t(vi_alpha_decomp[, i[1] + 0:(d_j[j]-1), drop = F]) %*%
        vi_alpha_decomp[, i[2] + 0:(d_j[jprime]-1), drop = F]
    })
    combo_j[[jprime]] <- combinations
  }
  row_combinations[[j]] <- combo_j
}

unique_map_pairwise <- lapply(purrr::array_branch(combn(1:number_of_RE, 2), margin = 2),
  FUN=function(i){
    unique(Mmap[,c(i)]) 
  })

start_base_Z
mapping_J
re_po

# for (int k = 0; k < d_j; ++k){
#   for (int kprime = 0; kprime < d_jprime; ++kprime){
#     outer_alpha(k,kprime) = L.col(m_ij + k - 1).cwiseProduct(L.col(m_ijprime + kprime - 1)).sum();
#   }
# }

est_vglmer <- vglmer(formula = fmla,
                     data = data, 
                     control = vglmer_control(
                       parameter_expansion = 'translation',
                       factorization_method = 'partial'),
                     family = 'binomial')

#FOR SIMPLE

R_ridge <- base_Rvec_ridge(
  vi_alpha_decomp = vi_alpha_decomp,
  d_j = d_j,
  Z = Z,
  diag_vi_pg_mean = diag_vi_pg_mean,
  outer_alpha_RE_positions = outer_alpha_RE_positions
)

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


R_ridge <- vecR_ridge_general(L = vi_alpha_decomp, 
      pg_mean = diag(diag_vi_pg_mean), 
      Z = mapping_new_Z, M = Mmap, 
      mapping_J = mapping_J, d = d_j, start_z = start_base_Z,
      diag_only = FALSE  )
J <- number_of_RE


est_vglmer <- vglmer(fmla,
                     data = data, 
                     control = vglmer_control(
                       parameter_expansion = 'translation',
                       factorization_method = 'partial'),
                     family = 'binomial')


store_re_id <- store_id <- list()
for (j in 1:J){
  store_re_id_j <- store_id_j <- list()
  for (jprime in 1:j){
    print(c(j, jprime))
    umap <- unique(Mmap[, c(j, jprime)])
    store_re_id_j[[jprime]] <- purrr::array_branch(umap, margin = 1)
    id_lookup <- lapply(1:nrow(umap), FUN=function(i){
      umap_r <- umap[i,]
      id_r <- which( (Mmap[,j] %in% umap_r[1]) & (Mmap[,jprime] %in% umap_r[2]))
      return(id_r)
    })
    store_id_j[[jprime]] <- id_lookup
  }
  store_id[[j]] <- store_id_j
  store_re_id[[j]] <- store_re_id_j
}

pg_mean <- diag(diag_vi_pg_mean)
cumdsq <- cumsum(c(1, d_j^2))
all_loop <- matrix(0, nrow = sum(d_j^2), ncol = sum(d_j^2))
for (j in 1:J){
  print(j)
  for (jprime in 1:j){
      cat('|')
      length_j <- length(store_id[[j]][[jprime]])
      loop_terms <- matrix(0, nrow = d_j[j]^2, ncol = d_j[jprime]^2)
      outeralpha <- t(vi_alpha_decomp[, col_re[1] + 0:(d_j[j]-1), drop = F]) %*%
        vi_alpha_decomp[, col_re[2] + 0:(d_j[jprime]-1), drop = F]
      for (comb in 1:length_j){
        col_re <- store_re_id[[j]][[jprime]][[comb]]
        u_comb <- store_id[[j]][[jprime]][[comb]]
        
        outerz <- t(parsed_RE_groups$design[[j]][u_comb,]) %*% Diagonal(x = pg_mean[u_comb]) %*%
          parsed_RE_groups$design[[jprime]][u_comb,]
        
        loop_terms <- loop_terms + kronecker(outeralpha, outerz)  
      }
      all_loop[cumdsq[j] + 1:d_j[j]^2 - 1 , cumdsq[jprime] + 1:d_j[jprime]^2 - 1] <- as.matrix(loop_terms)
  }
  print(image(all_loop %>% drop0))
  copy_loop <- all_loop
  copy_loop[upper.tri(copy_loop, diag = F)] <- 0
  copy_loop <- drop0(copy_loop)
  
  if (max(abs(copy_loop[1:cumdsq[j], 1:cumdsq[j]] - copy_R[1:cumdsq[j],1:cumdsq[j]])) > 1e-8){
    stop('..')
  }
  
}

all.equal(copy_loop, copy_R) #(as.vector(copy_loop) - as.vector(copy_R))

all_shards <- list()
for (j in 1:J){
  print(j)
  shard_j <- list()
  for (jprime in 1:j){
    length_j <- length(store_id[[j]][[jprime]])
    counter <- 1
    list_shard <- list()
    for (comb in 1:length_j){
      col_re <- store_re_id[[j]][[jprime]][[comb]]
      u_comb <- store_id[[j]][[jprime]][[comb]]
      list_shard[[counter]] <- list(one = as.matrix(parsed_RE_groups$design[[j]][u_comb,,drop=F]),
           two = as.matrix(parsed_RE_groups$design[[jprime]][u_comb,,drop=F]))
      counter <- counter + 1
    }
    one <- do.call('rbind', lapply(list_shard, FUN=function(i){i$one}))
    two <- do.call('rbind', lapply(list_shard, FUN=function(i){i$two}))
    list_shard <- list(cbind(one, two))
    shard_j[[jprime]] <- list_shard
  }
  all_shards[[j]] <- shard_j
}

v <- lapply(all_shards, FUN=function(j){
  cat(".")
  lapply(j, FUN=function(k){
    cat("|")
    lapply(k, FUN=function(l){crossprod(l$one, l$two)})
  })
})

pg_mean <- diag(diag_vi_pg_mean)
do.call('rbind', lapply(all_shards[[4]][[1]], FUN=function(i){i$one}))

vnew <- vecR_ridge_new(L = vi_alpha_decomp, pg_mean = pg_mean,
               mapping_J = mapping_J, d = d_j, store_id = store_id,
               store_re_id = store_re_id, store_design = parsed_RE_groups$design,
               diag_only = FALSE)

diag(drop0(zapsmall(vnew, 15)))
diag(R_ridge)

plot(as.vector(R_ridge - vnew))
round(cbind(diag(R_ridge), diag(vnew)), 3)


microbenchmark::microbenchmark(
  vecR_ridge_new(L = vi_alpha_decomp, pg_mean = pg_mean,
                 mapping_J = mapping_J, d = d_j, store_id = store_id,
                 store_re_id = store_re_id, store_design = parsed_RE_groups$design,
                 diag_only = FALSE),
  vecR_ridge_general(L = vi_alpha_decomp, pg_mean = pg_mean,
                     Z = mapping_new_Z, M = Mmap, mapping_J = mapping_J,
                     d = d_j, start_z = start_base_Z, diag_only = FALSE),
  times = 2
)

saveRDS(
  list(L = vi_alpha_decomp, pg_mean = diag(diag_vi_pg_mean),
     Z = mapping_new_Z, M = Mmap, store_id = store_id,
     store_re_id = store_re_id, store_design = store_design,
     mapping_J = mapping_J, d = d_j, start_z = start_base_Z),
  'tmp_data.RDS')



vecR_ridge_new(L = vi_alpha_decomp, pg_mean = diag(diag_vi_pg_mean),
               mapping_J = mapping_J, d = d_j, store_id = store_id,
               store_re_id = store_re_id, store_design = parsed_RE_groups$design,
               diag_only = T)

data <- readRDS('tmp_data.RDS')
for (v in names(data)){ assign(v, data[[v]])}


vn <- vecR_ridge_new(L = L, pg_mean = pg_mean,
               mapping_J = mapping_J, d = d, store_id = store_id,
               store_re_id = store_re_id, store_design = store_design,
               diag_only = FALSE)

vo <- vecR_ridge_general(L = L, 
                   pg_mean = pg_mean, 
                   Z = Z, M = M, 
                   mapping_J = mapping_J, d = d, start_z = start_z,
                   diag_only = FALSE)
range(unlist(store_id))
all.equal(vn, vo)
microbenchmark::microbenchmark(
  vecR_ridge_new(L = L, pg_mean = pg_mean,
                 mapping_J = mapping_J, d = d, store_id = store_id,
                 store_re_id = store_re_id, store_design = store_design,
                 diag_only = TRUE),
  vecR_ridge_general(L = L, 
                     pg_mean = pg_mean, 
                     Z = Z, M = M, 
                     mapping_J = mapping_J, d = d, start_z = start_z,
                     diag_only = TRUE),
  vecR_ridge_new(L = L, pg_mean = pg_mean,
       mapping_J = mapping_J, d = d, store_id = store_id,
       store_re_id = store_re_id, store_design = store_design,
       diag_only = FALSE),
  vecR_ridge_general(L = L, 
         pg_mean = pg_mean, 
         Z = Z, M = M, 
         mapping_J = mapping_J, d = d, start_z = start_z,
         diag_only = FALSE), times = 3)
