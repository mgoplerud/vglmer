
safe_convert <- function(x){
  if (isDiagonal(x)){
    out <- diag(x)
    lout <- seq_len(length(out)) - 1
    out <- cbind(lout, lout, out)
    colnames(out) <- c('i', 'j', 'x')
    # out <- with(attributes(as(as(as.matrix(x), "sparseMatrix"), "dgTMatrix")), cbind(i, j, x))
  }else{
    out <- with(attributes(as(as(x, "sparseMatrix"), "dgTMatrix")), cbind(i, j, x))
  }
  return(out)
}

#' @import Matrix
#' @importFrom methods as
make_mapping_alpha <- function(sigma, px.R = FALSE) {
  
  if (!px.R) {
    lapply(sigma, FUN = function(i) {
      
      sparse_i <- safe_convert(i)
      sparse_i <- sparse_i[sparse_i[, 1] >= sparse_i[, 2], , drop = F]
      return(sparse_i)
      
    })
  } else {
    lapply(sigma, FUN = function(i) {
      
      sparse_i <- safe_convert(i)
      return(sparse_i)
      
    })
  }
}

prepare_T <- function(mapping, levels_per_RE, variables_per_RE, running_per_RE, num_REs, cyclical = FALSE, px.R = FALSE) {
  if (!cyclical) {
    RE_T <- matrix(nrow = 0, ncol = 3)
  } else {
    RE_T <- as.list(rep(NA, num_REs))
  }
  for (v in 1:num_REs) {
    mapping_v <- mapping[[v]]
    
    if (cyclical) {
      mapping_id_i <- rep(mapping_v[, 1], levels_per_RE[v]) +
        rep(seq(1, 1 + (levels_per_RE[v] - 1) * variables_per_RE[v], by = variables_per_RE[v]), each = nrow(mapping_v))
      mapping_id_j <- rep(mapping_v[, 2], levels_per_RE[v]) +
        rep(seq(1, 1 + (levels_per_RE[v] - 1) * variables_per_RE[v], by = variables_per_RE[v]), each = nrow(mapping_v))
      
      mapping_id_x <- rep(mapping_v[, 3], levels_per_RE[v])
      RE_T[[v]] <- sparseMatrix(i = mapping_id_i, j = mapping_id_j, x = mapping_id_x, symmetric = T)
    } else {
      mapping_id_i <- rep(mapping_v[, 1], levels_per_RE[v]) +
        rep(seq(1, 1 + (levels_per_RE[v] - 1) * variables_per_RE[v], by = variables_per_RE[v]), each = nrow(mapping_v))
      mapping_id_j <- rep(mapping_v[, 2], levels_per_RE[v]) +
        rep(seq(1, 1 + (levels_per_RE[v] - 1) * variables_per_RE[v], by = variables_per_RE[v]), each = nrow(mapping_v))
      
      mapping_id_x <- rep(mapping_v[, 3], levels_per_RE[v])
      mapping_id_i <- running_per_RE[v] + mapping_id_i
      mapping_id_j <- running_per_RE[v] + mapping_id_j
      RE_T <- rbind(RE_T, cbind(mapping_id_i, mapping_id_j, mapping_id_x))
    }
  }
  if (!cyclical) {
    if (px.R) {
      RE_T <- sparseMatrix(i = RE_T[, 1], j = RE_T[, 2], x = RE_T[, 3], symmetric = F)
    } else {
      RE_T <- sparseMatrix(i = RE_T[, 1], j = RE_T[, 2], x = RE_T[, 3], symmetric = T)
    }
  }
  return(RE_T)
}

multi_lgamma <- function(a, p) {
  return(lmvgamma(x = a, p = p))
  # if (length(p) != 1){stop('P must have length 1')}
  # #if (any(a + 1 < p)){stop('Undefined for a < p - 1')}
  # term.1 <- log(pi) * (p) * (p - 1) / 4
  # term.2 <- mapply(1:p, SIMPLIFY = FALSE, FUN=function(i){
  #   matrix(lgamma(a + (1 - i) / 2))
  # })
  # term.2 <- as.vector(Reduce('+', term.2))
  # return(term.1 + term.2)
}

multi_digamma <- function(a, p) {
  return(mvdigamma(x = a, p = p))
  # if (length(p) != 1){stop('P must have length 1')}
  # #if (any(a + 1 < p)){stop('Undefined for a < p - 1')}
  # term.1 <- mapply(1:p, SIMPLIFY = FALSE, FUN=function(i){
  #   digamma(a + (1 - i) / 2)
  # })
  # term.1 <- as.vector(Reduce('+', term.1))
  # return(term.1)
}


#' Simple EM algorithm for starting values.
#'
#' Use ridge penalty to prevent separation. Not be called by user!
#'
#' @param X Design matrix
#' @param Z RE design matrix
#' @param s (y_i - n_i)/2 for polya-gamma input
#' @param y Raw observed y_i
#' @param est_r Initial r value (not updated!)
#' @param pg_b n_i as vector input
#' @param iter iterations
#' @param ridge variance of ridge prior
#' @name simple_EM
#' @keywords internal
#' @importFrom stats runif
EM_prelim_logit <- function(X, Z, s, pg_b, iter, ridge = 2) {
  jointXZ <- cbind(X, Z)
  N <- nrow(X)
  
  EM_beta <- rep(0, ncol(jointXZ))
  
  if (all(jointXZ[, 1] == 1)) {
    EM_beta[1] <- qlogis(sum(s + pg_b / 2) / sum(pg_b))
  }
  if (EM_beta[1] == 0) {
    EM_beta[1] <- runif(1, -.1, .1)
  }
  EM_variance <- sparseMatrix(i = 1:ncol(jointXZ), j = 1:ncol(jointXZ), x = 1 / ridge)
  
  for (it in 1:iter) {
    EM_pg_c <- jointXZ %*% EM_beta
    EM_pg_mean <- pg_b / (2 * EM_pg_c) * tanh(EM_pg_c / 2)
    if (any(abs(EM_pg_c) < 1e-10)) {
      tiny_c <- which(abs(EM_pg_c) < 1e-10)
      EM_pg_mean[tiny_c] <- pg_b[tiny_c] / 4
    }
    EM_pg_diag_sqrt <- sparseMatrix(i = 1:N, j = 1:N, x = sqrt(EM_pg_mean))
    
    EM_beta <- solve(Matrix::Cholesky( crossprod(EM_pg_diag_sqrt %*% jointXZ) + EM_variance),
                     t(jointXZ) %*% (s) )
    
    # EM_beta <- LinRegChol(X = jointXZ, omega = EM_pg_diag, y = s, prior_precision = EM_variance)$mean
  }
  output <- list(beta = EM_beta[1:ncol(X)], alpha = EM_beta[-1:-ncol(X)])
  return(output)
}

#' @rdname simple_EM
#' @importFrom stats runif
EM_prelim_nb <- function(X, Z, y, est_r, iter, ridge = 2) {
  if (is.null(Z)) {
    jointXZ <- drop0(X)
  } else {
    jointXZ <- cbind(X, Z)
  }
  N <- nrow(jointXZ)
  
  EM_beta <- rep(0, ncol(jointXZ))
  if (all(jointXZ[, 1] == 1)) {
    EM_beta[1] <- log(mean(y))
  }
  if (EM_beta[1] == 0) {
    EM_beta[1] <- runif(1, -.1, .1)
  }
  
  EM_variance <- sparseMatrix(i = 1:ncol(jointXZ), j = 1:ncol(jointXZ), x = 1 / ridge)
  for (it in 1:iter) {
    pg_c <- as.vector(jointXZ %*% EM_beta - log(est_r))
    pg_b <- y + est_r
    pg_mean <- as.vector(pg_b / (2 * pg_c) * tanh(pg_c / 2))
    
    if (any(abs(pg_c) < 1e-10)) {
      tiny_c <- which(abs(pg_c) < 1e-10)
      pg_mean[tiny_c] <- pg_b[tiny_c] / 4
    }
    
    adj_out <- (y - est_r) / 2 + pg_mean * log(est_r)
    # EM_beta <- LinRegChol(X = jointXZ, omega = sparseMatrix(i = 1:N, j = 1:N, x = pg_mean), y = adj_out, prior_precision = EM_variance)$mean
    
    EM_beta <- solve(Matrix::Cholesky(  t(jointXZ) %*% sparseMatrix(i = 1:N, j = 1:N, x = pg_mean) %*% jointXZ + EM_variance),
                     t(jointXZ) %*% (adj_out) )
    
  }
  
  output <- list(beta = EM_beta[1:ncol(X)], alpha = EM_beta[-1:-ncol(X)])
  return(output)
}

make_log_invwishart_constant <- function(nu, Phi) {
  p <- ncol(Phi)
  output <- nu / 2 * log(det(Phi)) - (nu * p) / 2 * log(2) - multi_lgamma(a = nu / 2, p = p)
  return(output)
}

calculate_ELBO <- function(family, ELBO_type, factorization_method,
   # Fixed constants or priors
   d_j, g_j, prior_sigma_alpha_phi, prior_sigma_alpha_nu, 
   iw_prior_constant, choose_term,
   store_assignment_Z, store_design_Z, outer_alpha_RE_positions,
   # Data
   X, Z, s, y,
   # PolyaGamma Parameters
   vi_pg_b, vi_pg_mean, vi_pg_c,
   # Sigma Parameters
   vi_sigma_alpha, vi_sigma_alpha_nu, vi_sigma_outer_alpha,
   # Beta Parameters / Alpha Parameters
   vi_beta_mean, vi_beta_decomp,
   vi_alpha_mean, vi_alpha_decomp,
   log_det_beta_var , log_det_alpha_var,
   vi_alpha_var = NULL, vi_beta_var = NULL,
   cyclical_pos = NULL,
   log_det_joint_var = NULL,
   vi_joint_decomp = NULL,
   vi_P = NULL, log_det_M_var = NULL, log_det_C_var = NULL,
   vi_C_uncond = NULL, vi_C_mean = NULL, vi_M_mean = NULL,
   vi_FS_MC = NULL, lookup_marginal = NULL,
   vi_FS_MM = NULL, vi_M_list = NULL,
   design_M = NULL,  vi_M_B = NULL, vi_M_var_flat = NULL,
   design_C = NULL,
   # r Parameters
   vi_r_mu = NULL, vi_r_mean = NULL,
   vi_r_sigma = NULL,
   #linear parameters
   vi_sigmasq_a = NULL, vi_sigmasq_b = NULL,
   vi_sigmasq_prior_a = NULL, vi_sigmasq_prior_b = NULL,
   # huang_wand parameters
   do_huangwand = NULL, vi_a_a_jp = NULL, vi_a_b_jp = NULL,
   vi_a_nu_jp = NULL, vi_a_APRIOR_jp = NULL
) {
  
  ####
  ## PREPARE INTERMEDIATE QUANTITES
  ###
  
  N <- nrow(X)
  
  if (factorization_method == "partially_factorized"){
    ex_XBZA <- design_C %*% vi_C_mean 
    if (ncol(design_M) > 0){
      ex_XBZA <- ex_XBZA + design_M %*% do.call('c', vi_M_mean) 
    }
  }else{
    ex_XBZA <- (X %*% vi_beta_mean + Z %*% vi_alpha_mean)
  }
  if (family == 'negbin'){
    ex_XBZA <- ex_XBZA - vi_r_mu
  }
  # quadratic var, i.e. Var(x_i^T beta + z_i^T alpha)
  if (factorization_method == "weak") {
    
    if (is.null(vi_joint_decomp)) {
      stop("Need to provide joint decomposition for ELBO weak")
    }
    
    var_XBZA <- rowSums((cbind(X, Z) %*% t(vi_joint_decomp))^2)
    
    if (family == 'negbin'){
      var_XBZA <- var_XBZA + vi_r_sigma
    }
  } else if (factorization_method == "partially_factorized") {
    
    var_XBZA <- rowSums( (design_C %*% vi_C_uncond) * design_C)
    
    var_XBZA <- var_XBZA + rowSums(mapply(vi_FS_MM, vi_M_var_flat, lookup_marginal, names(lookup_marginal), FUN=function(xi,zi,gi, n){
      if (ncol(xi) > 0){
        rowSums( xi * (gi %*% zi))
      }else{
        return(rep(0, nrow(xi)))
      }
    }))
    # Covariance
    var_XBZA <- var_XBZA + -2 * rowSums(
      mapply(vi_FS_MC, vi_M_B, FUN=function(data_j, B_j){
        as.vector(data_j %*% B_j)})
    )
    
    # # # Variance of Marginal
    # var_XBZA <- var_XBZA + rowSums(mapply(vi_M_list, vi_M_var_flat, FUN=function(data_j, var_j){
    #   rowSums( (data_j %*% Diagonal(x = var_j)) * data_j)
    # }))
    # # # Covariance
    # var_XBZA <- var_XBZA + -2 * rowSums(
    #   mapply(vi_M_list, vi_M_B, FUN=function(data_j, B_j){
    #     rowSums((data_j %*% t(B_j)) * design_C)})
    # )
    #   
  } else {
    
    beta_quad <- rowSums((X %*% t(vi_beta_decomp))^2)
    alpha_quad <- rowSums((Z %*% t(vi_alpha_decomp))^2)
    var_XBZA <- beta_quad + alpha_quad
    if (family == 'negbin'){
      var_XBZA <- var_XBZA + vi_r_sigma 
    }
  }
  
  # Prepare vi_sigma_alpha
  moments_sigma_alpha <- mapply(vi_sigma_alpha, vi_sigma_alpha_nu, d_j, SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
    inv_phi <- solve(phi)
    
    sigma.inv <- nu * inv_phi
    
    # ln.det <- - (multi_digamma(a = nu/2, p = d) + d * log(2) + log(det(inv_phi)) )
    ln.det <- log(det(phi)) - sum(digamma((nu - 1:d + 1) / 2)) - d * log(2)
    return(list(sigma.inv = sigma.inv, ln.det = ln.det))
  })
  
  ln_det_sigma_alpha <- sapply(moments_sigma_alpha, FUN = function(i) {
    i$ln.det
  })
  inv_sigma_alpha <- lapply(moments_sigma_alpha, FUN = function(i) {
    i$sigma.inv
  })
  ## GET the terms for the expectation
  ## of the log-complete data given the variational distribution.
  if (ELBO_type == "augmented") {
    if (family == "linear") {
      e_ln_sigmasq <- log(vi_sigmasq_b) - digamma(vi_sigmasq_a)
      e_inv_sigmasq <- vi_sigmasq_a/vi_sigmasq_b
      logcomplete_1 <- sum(-1/2 * ((y - ex_XBZA)^2 + var_XBZA) * e_inv_sigmasq) +
        -1/2 * length(y) * (log(2 * pi) + e_ln_sigmasq)
      #Add log prior
      logcomplete_1 <- logcomplete_1 + 
        (-vi_sigmasq_prior_a - 1) * e_ln_sigmasq +
        -vi_sigmasq_prior_b * e_inv_sigmasq
    } else {
      # Get the terms for the p(y, w | alpha, beta, Sigma) EXCLUDING the intractable PG.
      logcomplete_1 <- -sum(vi_pg_b) * log(2) +
        as.vector(t(s) %*% ex_XBZA - 1 / 2 * t(ex_XBZA) %*% Diagonal(x = vi_pg_mean) %*% ex_XBZA) +
        -1 / 2 * sum(var_XBZA * vi_pg_mean)
    }
    # Get the terms for p(alpha | Sigma)
    
    if (family == 'linear'){
      e_ln_sigmasq <- log(vi_sigmasq_b) - digamma(vi_sigmasq_a)
      
      logcomplete_2 <- sum(-d_j * g_j / 2 * log(2 * pi) - g_j / 2 * ln_det_sigma_alpha) +
        -e_inv_sigmasq * 1 / 2 * sum(mapply(inv_sigma_alpha, vi_sigma_outer_alpha, FUN = function(a, b) {
          sum(diag(a %*% b))
        }))
      logcomplete_2 <-  logcomplete_2 +
        -1/2 * sum(d_j * g_j) * e_ln_sigmasq
    }else{
      
      logcomplete_2 <- sum(-d_j * g_j / 2 * log(2 * pi) - g_j / 2 * ln_det_sigma_alpha) +
        -1 / 2 * sum(mapply(inv_sigma_alpha, vi_sigma_outer_alpha, FUN = function(a, b) {
          sum(diag(a %*% b))
        }))
    }
    
    ## GET THE ENTROPY
    # Entropy for p(beta,alpha)
    if (factorization_method == "weak") {
      entropy_1 <- ncol(vi_joint_decomp) / 2 * log(2 * pi * exp(1)) +
        1 / 2 * log_det_joint_var
    } else if (factorization_method == "partially_factorized") {
      entropy_1 <- ncol(design_C) / 2 * log(2 * pi * exp(1)) + 1 / 2 * log_det_C_var +
        ncol(design_M) / 2 * log(2 * pi * exp(1)) + 1 / 2 * log_det_M_var
    } else {
      entropy_1 <- nrow(vi_beta_mean) / 2 * log(2 * pi * exp(1)) + 1 / 2 * log_det_beta_var +
        nrow(vi_alpha_mean) / 2 * log(2 * pi * exp(1)) + 1 / 2 * log_det_alpha_var
    }
    #ENTROPY FOR LINK SPECIFIC PARAMETERS
    if (family == 'linear'){
      entropy_2 <- vi_sigmasq_a + log(vi_sigmasq_b) + lgamma(vi_sigmasq_a) + 
        -(vi_sigmasq_a + 1) * digamma(vi_sigmasq_a)
    }else{
      # Entropy for Polya-Gamma EXCLUDING intractable term that cancels
      entropy_2 <- sum(vi_pg_b * vi_pg_c / 4 * tanh(vi_pg_c / 2) - vi_pg_b * log(cosh(vi_pg_c / 2)))
    }
    # Entropy Wisharts
    entropy_3 <- -mapply(vi_sigma_alpha_nu, vi_sigma_alpha, FUN = function(nu, Phi) {
      make_log_invwishart_constant(nu = nu, Phi = Phi)
    }) +
      (vi_sigma_alpha_nu + d_j + 1) / 2 * ln_det_sigma_alpha +
      1 / 2 * mapply(vi_sigma_alpha, inv_sigma_alpha, FUN = function(a, b) {
        sum(diag(a %*% b))
      })
    entropy_3 <- sum(entropy_3)
    
  } else if (ELBO_type == "profiled") {
    
    vi_r_var <- (exp(vi_r_sigma) - 1) * vi_r_mean^2
    
    psi <- ex_XBZA + vi_r_mu
    zVz <- var_XBZA - vi_r_sigma
    
    logcomplete_1 <- VEM.PELBO.r(ln_r = vi_r_mu, y, psi, zVz) +
      1 / 2 * VEM.PELBO.r_hessian(ln_r = vi_r_mu, y, psi, zVz) * vi_r_sigma
    
    # logcomplete_1a <- sum(lgamma(y + vi_r_hat)) +
    #   - N * lgamma(vi_r_hat) - (vi_r_hat) * N * log(2) +
    #   as.vector(t((y - exp(vi_r_hat))/2) %*% ex_XBZA)
    # logcomplete_1b <- sum(-(y + vi_r_mean) * log(cosh(1/2 * sqrt(ex_XBZA^2 + vi_r_sigma + var_XBZA))))
    # logcomplete_1c <- N/2 * vi_r_mean * vi_r_sigma
    # logcomplete_1 <- logcomplete_1a + logcomplete_1b + logcomplete_1c + choose_term
    
    logcomplete_2 <- sum(-d_j * g_j / 2 * log(2 * pi) - g_j / 2 * ln_det_sigma_alpha) +
      -1 / 2 * sum(mapply(inv_sigma_alpha, vi_sigma_outer_alpha, FUN = function(a, b) {
        sum(diag(a %*% b))
      }))
    if (factorization_method == "weak") {
      entropy_1 <- ncol(vi_joint_decomp) / 2 * log(2 * pi * exp(1)) +
        1 / 2 * log_det_joint_var
    } else {
      entropy_1 <- ncol(vi_beta_decomp) / 2 * log(2 * pi * exp(1)) + 1 / 2 * log_det_beta_var +
        ncol(vi_alpha_decomp) / 2 * log(2 * pi * exp(1)) + 1 / 2 * log_det_alpha_var
    }
    # Entropy of q(ln r)
    if (vi_r_sigma == 0) {
      entropy_2 <- 0
    } else {
      entropy_2 <- 1 / 2 * log(2 * pi * exp(1) * vi_r_sigma)
    }
    
  } else {
    stop("ELBO must be profiled or augmented")
  }
  ###############
  # Log Complete and Entropy for p(Sigma_j) or similar
  ###############
  if (do_huangwand){
    E_ln_vi_a <- mapply(vi_a_a_jp, vi_a_b_jp, FUN=function(tilde.a, tilde.b){
      sum(log(tilde.b) - digamma(tilde.a))
    })
    E_inv_v_a <- mapply(vi_a_a_jp, vi_a_b_jp, vi_a_nu_jp, SIMPLIFY = FALSE, FUN=function(tilde.a, tilde.b, nu){
      2 * nu * Diagonal(x = tilde.a/tilde.b)
    })
    logcomplete_3 <- 0 + # flat prior on beta
      sum(
        iw_prior_constant +
          - (vi_a_nu_jp + d_j - 1)/2 * (d_j * log(2 * vi_a_nu_jp) + E_ln_vi_a) +
          -(2 * d_j + vi_a_nu_jp) / 2 * ln_det_sigma_alpha +
          -1 / 2 * mapply(E_inv_v_a, inv_sigma_alpha, FUN = function(a, b) {
            sum(diag(a %*% b))
          })
      )
    # inv_sigma_alpha <<- inv_sigma_alpha
    # iw_prior_constant <<- iw_prior_constant
    # ln_det_sigma_alpha <<- ln_det_sigma_alpha
    # d_j <<- d_j; vi_a_a_jp <<- vi_a_a_jp
    # vi_a_b_jp <<- vi_a_b_jp
    # E_ln_vi_a <<- E_ln_vi_a
    # vi_a_APRIOR_jp <<- vi_a_APRIOR_jp
    logcomplete_3_a <- mapply(d_j, vi_a_a_jp, vi_a_b_jp, E_ln_vi_a, 
                              vi_a_APRIOR_jp, 
                              FUN=function(d, tilde.a, tilde.b, E_ln_vi_a.j, APRIOR.j){
                                1/2 * sum(log(1/APRIOR.j^2)) - d * lgamma(1/2) - 3/2 * E_ln_vi_a.j +
                                  sum(-1/APRIOR.j^2 * tilde.a/tilde.b)
                              })
    logcomplete_3 <- logcomplete_3 + sum(logcomplete_3_a)
  }else{
    logcomplete_3 <- 0 + # flat prior on beta
      sum(
        iw_prior_constant +
          -(prior_sigma_alpha_nu + d_j + 1) / 2 * ln_det_sigma_alpha +
          -1 / 2 * mapply(prior_sigma_alpha_phi, inv_sigma_alpha, FUN = function(a, b) {
            sum(diag(a %*% b))
          })
      )
  }
  
  entropy_3 <- -mapply(vi_sigma_alpha_nu, vi_sigma_alpha, FUN = function(nu, Phi) {
    make_log_invwishart_constant(nu = nu, Phi = Phi)
  }) +
    (vi_sigma_alpha_nu + d_j + 1) / 2 * ln_det_sigma_alpha +
    1 / 2 * mapply(vi_sigma_alpha, inv_sigma_alpha, FUN = function(a, b) {
      sum(diag(a %*% b))
    })
  entropy_3 <- sum(entropy_3)
  #########
  # Optional Entropy if using Huang and Wand (2013) prior
  #########
  if (do_huangwand){
    entropy_4 <- sum(mapply(vi_a_a_jp, vi_a_b_jp, FUN=function(tilde.a, tilde.b){
      sum(tilde.a + log(tilde.b) + lgamma(tilde.a) - (1 + tilde.a) * digamma(tilde.a))
    }))
  }else{
    entropy_4 <- 0
  }
  
  ###Combine all of the terms together
  
  logcomplete <- logcomplete_1 + logcomplete_2 + logcomplete_3 +
    choose_term
  
  entropy <- entropy_1 + entropy_2 + entropy_3 + entropy_4
  ELBO <- entropy + logcomplete
  
  return(data.frame(
    ELBO, logcomplete, entropy, logcomplete_1,
    logcomplete_2, logcomplete_3, entropy_1, entropy_2, entropy_3, entropy_4
  ))
}



update_r <- function(vi_r_mu, vi_r_sigma, y, X, Z, factorization_method,
                     vi_beta_mean, vi_alpha_mean,
                     vi_joint_decomp, vi_beta_decomp, vi_alpha_decomp,
                     vi_r_method) {
  if (vi_r_method == "fixed") {
    return(c(vi_r_mu, vi_r_sigma))
  }
  # Get intermediate quantities
  ex_XBZA <- (X %*% vi_beta_mean + Z %*% vi_alpha_mean)
  # quadratic var, i.e. Var(x_i^T beta + z_i^T alpha)
  if (factorization_method == "weak") {
    if (is.null(vi_joint_decomp)) {
      stop("Need to provide joint decomposition for ELBO weak")
    }
    var_XBZA <- rowSums((cbind(X, Z) %*% t(vi_joint_decomp))^2)
  } else {
    beta_quad <- rowSums((X %*% t(vi_beta_decomp))^2)
    alpha_quad <- rowSums((Z %*% t(vi_alpha_decomp))^2)
    var_XBZA <- beta_quad + alpha_quad
  }
  
  N <- length(y)
  
  # vi_r_mu <<- vi_r_mu
  # vi_r_sigma <<- vi_r_sigma
  # ex_XBZA <<- ex_XBZA
  # var_XBZA <<- var_XBZA
  # y <<- y
  # N <<- N
  
  if (vi_r_method == "delta") {
    opt_vi_r <- optim(
      par = c(vi_r_mu, log(vi_r_sigma)), fn = VEM.delta_method,
      y = y, psi = ex_XBZA, zVz = var_XBZA,
      control = list(fnscale = -1), method = "L-BFGS"
    )
    
    prior_vi_r <- VEM.delta_method(
      par = c(vi_r_mu, log(vi_r_sigma)), y = y,
      psi = ex_XBZA, zVz = var_XBZA
    )
    if (opt_vi_r$value < prior_vi_r) {
      warning("Optim for r decreased objective.")
      out_par <- c(vi_r_mu, vi_r_sigma)
    } else {
      out_par <- c(opt_vi_r$par[1], exp(opt_vi_r$par[2]))
    }
  } else if (vi_r_method %in% c("VEM", "Laplace")) {
    opt_vi_r <- optim(
      par = vi_r_mu, fn = VEM.PELBO.r, gr = VEM.PELBO.r_deriv,
      y = y, psi = ex_XBZA, zVz = var_XBZA,
      control = list(fnscale = -1), method = "L-BFGS"
    )
    
    if (vi_r_method == "Laplace") {
      proposed_vi_r_sigma <- -1 / VEM.PELBO.r_hessian(
        ln_r = opt_vi_r$par,
        y = y, psi = ex_XBZA, zVz = var_XBZA
      )
      # proposed_vi_r_sigma <- #as.numeric(-1/opt_vi_r$hessian)
    } else {
      proposed_vi_r_sigma <- 0
    }
    
    prior_vi_r <- VEM.PELBO.r(
      ln_r = vi_r_mu, y = y,
      psi = ex_XBZA, zVz = var_XBZA
    )
    
    if (opt_vi_r$value < prior_vi_r) {
      warning("Optim for r decreased objective.")
      out_par <- c(vi_r_mu, vi_r_sigma)
    } else {
      out_par <- c(opt_vi_r$par, proposed_vi_r_sigma)
    }
  } else {
    stop("vi_r method must be VI or VEM or fixed")
  }
  
  return(out_par)
}

VEM.PELBO.r <- function(ln_r, y, psi, zVz) {
  t1 <- -(y + exp(ln_r)) * log(2) + (y - exp(ln_r)) / 2 * (psi - ln_r) +
    -(exp(ln_r) + y) * log(cosh(1 / 2 * sqrt(zVz + (psi - ln_r)^2)))
  t2 <- lgamma(y + exp(ln_r)) - lgamma(exp(ln_r)) - lgamma(y + 1)
  return(sum(t1 + t2))
}

approx.lgamma <- function(x, mean_r, var_r) {
  input <- x + mean_r
  output <- lgamma(input) + 1 / 2 * psigamma(x = input, deriv = 1) * var_r
  return(output)
}

VEM.delta_method <- function(par, ln_r, y, psi, zVz) {
  mu <- par[1]
  sigma <- exp(par[2])
  obj <- VEM.PELBO.r(ln_r = mu, y, psi, zVz) +
    1 / 2 * VEM.PELBO.r_hessian(ln_r = mu, y, psi, zVz) * sigma +
    1 / 2 * log(sigma)
  return(obj)
}

sech <- function(x) {
  1 / cosh(x)
}

VEM.PELBO.r_deriv <- function(ln_r, y, psi, zVz) {
  N <- length(y)
  r <- exp(ln_r)
  meat <- sqrt(zVz + (psi - ln_r)^2)
  
  # -E^lnr PolyGamma[0, E^lnr] + E^lnr PolyGamma[0, E^lnr + y]
  deriv_normcon <- -N * r * psigamma(r) + r * sum(psigamma(y + r))
  # Mathematica Syntax for Derivative Ln[Cosh[1/2 * Sqrt[zVz + (psi - ln_r)^2]]]
  # -E^lnr Log[Cosh[1/2 Sqrt[(lnr - psi)^2 + v]]] - ((lnr - psi) (E^lnr + y) Tanh[
  #     1/2 Sqrt[(lnr - psi)^2 + v]])/(2 Sqrt[(lnr - psi)^2 + v])
  deriv_lncosh <- -r * log(cosh(1 / 2 * meat)) - (ln_r - psi) * (y + r) * tanh(1 / 2 * meat) / (2 * meat)
  deriv_lncosh <- sum(deriv_lncosh)
  # Mathematic Synax for
  # 1/2 (-y + E^lnr (1 + lnr - psi - Log[4]))
  deriv_prelim <- 1 / 2 * (-y + r * (1 + ln_r - psi - log(4)))
  deriv_prelim <- sum(deriv_prelim)
  return(deriv_normcon + deriv_lncosh + deriv_prelim)
}

VEM.PELBO.r_hessian <- function(ln_r, y, psi, zVz) {
  N <- length(y)
  r <- exp(ln_r)
  meat <- sqrt(zVz + (psi - ln_r)^2)
  
  # -E^lnr PolyGamma[0, E^lnr] + E^lnr PolyGamma[0, E^lnr + y] -
  #   E^(2 lnr) PolyGamma[1, E^lnr] + E^(2 lnr) PolyGamma[1, E^lnr + y]
  deriv_normcon <- -N * r * psigamma(r) + r * sum(psigamma(y + r)) +
    -N * r^2 * psigamma(r, deriv = 1) + r^2 * sum(psigamma(r + y, deriv = 1))
  # Mathematica Code
  # E^r Log[Cosh[1/2 Sqrt[(psi - r)^2 + z]]] + (
  # E^r (psi - r) Tanh[1/2 Sqrt[(psi - r)^2 + z]])/Sqrt[(psi - r)^2 +
  # z] + (-E^r - y) (((psi - r)^2 Sech[1/2 Sqrt[(psi - r)^2 + z]]^2)/(
  #   4 ((psi - r)^2 + z)) - ((psi - r)^2 Tanh[
  #     1/2 Sqrt[(psi - r)^2 + z]])/(2 ((psi - r)^2 + z)^(3/2)) +
  #     Tanh[1/2 Sqrt[(psi - r)^2 + z]]/(2 Sqrt[(psi - r)^2 + z]))
  deriv_lncosh <- -r * log(cosh(1 / 2 * meat)) + r * (psi - ln_r) * tanh(1 / 2 * meat) / meat +
    -(y + r) * ((psi - ln_r)^2 * sech(1 / 2 * meat)^2 / (4 * meat^2) - (psi - ln_r)^2 * tanh(1 / 2 * meat) / (2 * meat^3) +
                  tanh(1 / 2 * meat) / (2 * meat))
  deriv_lncosh <- sum(deriv_lncosh)
  # deriv_lncosh <- -(psi - ln_r)^2 * (r + y)/(2 * meat^2 * (1 + cosh(meat))) +
  # -r * log(cosh(1/2 * meat))  +
  # (zVz * (y + r) + r * 2 * (ln_r - psi) * meat^2) * tanh(1/2 * meat) / (2 * meat^3)
  # Mathematica
  # E^lnr - 1/2 E^lnr (-lnr + psi) - E^lnr Log[2]
  deriv_prelim <- N * r - 1 / 2 * r * sum(psi - ln_r) - N * r * log(2)
  
  return(deriv_normcon + deriv_lncosh + deriv_prelim)
}

#\sum_{j,g} E[tr(\alpha_{j,g}^T \Sigma_j^{-1} \alpha_{j,g})]
expect_alpha_prior_kernel <- function(vi_sigma_alpha, vi_sigma_alpha_nu, vi_sigma_outer_alpha, d_j){
  
  moments_sigma_alpha <- mapply(vi_sigma_alpha, vi_sigma_alpha_nu, 
                                d_j, SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
                                  inv_phi <- solve(phi)
                                  
                                  sigma.inv <- nu * inv_phi
                                  
                                  return(list(sigma.inv = sigma.inv))
                                })
  
  inv_sigma_alpha <- lapply(moments_sigma_alpha, FUN = function(i) {i$sigma.inv})
  
  out <- sum(mapply(inv_sigma_alpha, vi_sigma_outer_alpha, FUN = function(a, b) {sum(diag(a %*% b))}))
  
  return(out)
}	
