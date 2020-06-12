extract_group_memberships <- function(x, fr, drop.unused.levels){
  frloc <- factorize(x, frloc)
  if (is.null(ff <- tryCatch(eval(substitute(makeFac(fac), 
                                             list(fac = x[[3]])), frloc), error = function(e) NULL))) 
    stop("couldn't evaluate grouping factor ", deparse(x[[3]]), 
         " within model frame:", " try adding grouping factor to data ", 
         "frame explicitly if possible", call. = FALSE)
  if (all(is.na(ff))){ 
    stop("Invalid grouping factor specification, ", deparse(x[[3]]), 
         call. = FALSE)
  }
  if (drop.unused.levels){
    ff <- factor(ff, exclude = NA)
  }
  return(ff)
}

make_dgC <- function(x){
  if (!inherits(x, 'ddiMatrix')){
    x <- as(x, 'dgCMatrix')
  }else{
    x <- sparseMatrix(i = 1:nrow(x), j =1:nrow(x), x = diag(x))
  }
  return(x)
}

#' @import Matrix
#' @importFrom methods as
make_mapping_alpha <- function(sigma, px.R = FALSE){
  if (!px.R){
    lapply(sigma, FUN=function(i){
      sparse_i <- with(attributes(as(as(i, 'sparseMatrix'), 'dgTMatrix')), cbind(i,j,x))
      sparse_i <- sparse_i[sparse_i[,1] >= sparse_i[,2], , drop = F]
      return(sparse_i)
    })
  }else{
    lapply(sigma, FUN=function(i){
      sparse_i <- with(attributes(as(as(i, 'sparseMatrix'), 'dgTMatrix')), cbind(i,j,x))
      return(sparse_i)
    })
  }
}

prepare_T <- function(mapping, levels_per_RE, variables_per_RE, running_per_RE, num_REs, cyclical = FALSE, px.R = FALSE){
  if (!cyclical){
    RE_T <- matrix(nrow = 0, ncol = 3)
  }else{
    RE_T <- as.list(rep(NA, num_REs))
  }
  for (v in 1:num_REs){
    mapping_v <- mapping[[v]]
    
    if (cyclical){
      mapping_id_i <- rep(mapping_v[,1], levels_per_RE[v]) +
        rep(seq(1, 1 + (levels_per_RE[v]-1) * variables_per_RE[v], by = variables_per_RE[v]), each = nrow(mapping_v))
      mapping_id_j <- rep(mapping_v[,2], levels_per_RE[v]) +
        rep(seq(1, 1 + (levels_per_RE[v]-1) * variables_per_RE[v], by = variables_per_RE[v]), each = nrow(mapping_v))
      
      mapping_id_x <- rep(mapping_v[,3], levels_per_RE[v])
      RE_T[[v]] <- sparseMatrix(i = mapping_id_i, j = mapping_id_j, x = mapping_id_x, symmetric = T)
    }else{
      mapping_id_i <- rep(mapping_v[,1], levels_per_RE[v]) +
        rep(seq(1, 1 + (levels_per_RE[v]-1) * variables_per_RE[v], by = variables_per_RE[v]), each = nrow(mapping_v))
      mapping_id_j <- rep(mapping_v[,2], levels_per_RE[v]) +
        rep(seq(1, 1 + (levels_per_RE[v]-1) * variables_per_RE[v], by = variables_per_RE[v]), each = nrow(mapping_v))
      
      mapping_id_x <- rep(mapping_v[,3], levels_per_RE[v])
      mapping_id_i <- running_per_RE[v] + mapping_id_i
      mapping_id_j <- running_per_RE[v] + mapping_id_j
      RE_T <- rbind(RE_T, cbind(mapping_id_i, mapping_id_j, mapping_id_x))
    }
  }
  if (!cyclical){
    if (px.R){
      RE_T <- sparseMatrix(i = RE_T[,1], j = RE_T[,2], x =RE_T[,3], symmetric = F)
    }else{
      RE_T <- sparseMatrix(i = RE_T[,1], j = RE_T[,2], x =RE_T[,3], symmetric = T)
    }
  }    
  return(RE_T)
}

multi_lgamma <- function(a, p){
  
  return(lmvgamma(x = a, p = p))
  if (length(p) != 1){stop('P must have length 1')}
  #if (any(a + 1 < p)){stop('Undefined for a < p - 1')}
  term.1 <- log(pi) * (p) * (p - 1) / 4
  term.2 <- mapply(1:p, SIMPLIFY = FALSE, FUN=function(i){
    matrix(lgamma(a + (1 - i) / 2))
  })
  term.2 <- as.vector(Reduce('+', term.2))
  return(term.1 + term.2)
}

multi_digamma <- function(a, p){
  return(mvdigamma(x = a, p = p))
  if (length(p) != 1){stop('P must have length 1')}
  #if (any(a + 1 < p)){stop('Undefined for a < p - 1')}
  term.1 <- mapply(1:p, SIMPLIFY = FALSE, FUN=function(i){
    digamma(a + (1 - i) / 2)
  })
  term.1 <- as.vector(Reduce('+', term.1))
  return(term.1)
}


#' Simple EM algorithm for starting values.
#' 
#' Use ridge penalty to prevent separation. Not be called by user!
#' 
#' @param X Design matrix
#' @param Z RE design matrix
#' @param s (y_i - n_i)/2 for polya-gamma input
#' @param pg_b n_i as vector input
#' @param iter iterations
#' @param ridge variance of ridge prior
EM_prelim <- function(X, Z, s, pg_b, iter, ridge = 2){
  
  jointXZ <- cbind(X, Z)
  N <- nrow(X)
  
  EM_beta <- rep(0, ncol(jointXZ))
  
  EM_variance <- sparseMatrix(i = 1:ncol(jointXZ), j = 1:ncol(jointXZ), x = 1/ridge)
  
  for (it in 1:iter){
    
    EM_pg_c <- jointXZ %*% EM_beta
    EM_pg_mean <- pg_b/(2 * EM_pg_c) * tanh(EM_pg_c / 2)
    EM_pg_mean[which(abs(EM_pg_c) < 1e-10)] <- pg_b/4
    EM_pg_diag <- sparseMatrix(i=1:N, j=1:N, x = EM_pg_mean)
    
    EM_beta <- LinRegChol(X = jointXZ, omega = EM_pg_diag, y = s, prior_precision =  EM_variance)$mean
    
  }
  output <- list(beta = EM_beta[1:ncol(X)], alpha = EM_beta[-1:-ncol(X)])
  return(output)  
}

make_log_invwishart_constant <- function(nu, Phi){
  p <- ncol(Phi)
  output <- nu/2 * log(det(Phi)) - (nu * p)/2 * log(2) - multi_lgamma(a = nu/2, p = p)
  return(output)
}

calculate_ELBO <- function(ELBO_type, factorization_method,
                           #Fixed constants or priors
                           d_j, g_j, prior_sigma_alpha_phi, prior_sigma_alpha_nu, iw_prior_constant, choose_term,
                           #Data
                           X, Z, s, y,
                           #PolyaGamma Parameters
                           vi_pg_b, vi_pg_mean, vi_pg_c,
                           #Sigma Parameters
                           vi_sigma_alpha, vi_sigma_alpha_nu, vi_sigma_outer_alpha,
                           #Beta Parameters / Alpha Parameters
                           vi_beta_mean, vi_beta_decomp,
                           vi_alpha_mean, vi_alpha_decomp,
                           log_det_beta_var, log_det_alpha_var, 
                           log_det_joint_var = NULL,
                           vi_joint_decomp = NULL,
                           #r Parameters
                           vi_r_mu = NULL, vi_r_mean = NULL, vi_r_sigma = NULL
){
  ####
  ##PREPARE INTERMEDIATE QUANTITES
  ###
  N <- nrow(X)
  #linear predictor: E[XB + ZA - log(r)]
  ex_XBZA <- (X %*% vi_beta_mean + Z %*% vi_alpha_mean) - vi_r_mu
  #quadratic var, i.e. Var(x_i^T beta + z_i^T alpha)
  if (factorization_method == 'weak'){
    if (is.null(vi_joint_decomp) | is.null(log_det_joint_var)){stop('Need to provide joint decomposition for ELBO weak')}
    var_XBZA <- rowSums( (cbind(X,Z) %*% t(vi_joint_decomp))^2 )
  }else{
    beta_quad <- rowSums( (X %*% t(vi_beta_decomp))^2 )
    alpha_quad <- rowSums( (Z %*% t(vi_alpha_decomp))^2 )
    var_XBZA <- beta_quad + alpha_quad
  }
  #Prepare vi_sigma_alpha
  moments_sigma_alpha <- mapply(vi_sigma_alpha, vi_sigma_alpha_nu, d_j, SIMPLIFY = FALSE, FUN=function(phi, nu, d){
    inv_phi <- solve(phi)
    
    sigma.inv <- nu * inv_phi
    
    #ln.det <- - (multi_digamma(a = nu/2, p = d) + d * log(2) + log(det(inv_phi)) )
    ln.det <- log(det(phi)) - sum(digamma((nu - 1:d + 1)/2)) - d * log(2)
    return(list(sigma.inv = sigma.inv, ln.det = ln.det))
  })
  
  ln_det_sigma_alpha <- sapply(moments_sigma_alpha, FUN=function(i){i$ln.det})
  inv_sigma_alpha <- lapply(moments_sigma_alpha, FUN=function(i){i$sigma.inv})
  ##GET the terms for the expectation
  ##of the log-complete data given the variational distribution.
  if (ELBO_type == 'augmented'){
    #Get the terms for the p(y, w | alpha, beta, Sigma) EXCLUDING the intractable PG.
    logcomplete_1 <- -sum(vi_pg_b) * log(2) +
      as.vector(t(s) %*% ex_XBZA - 1/2 * t(ex_XBZA) %*% Diagonal(x = vi_pg_mean) %*% ex_XBZA) +
      -1/2 * sum(var_XBZA * vi_pg_mean)
    #Get the terms for p(alpha | Sigma)
    
    logcomplete_2 <- sum(-d_j * g_j / 2 * log(2 * pi) - g_j/2 * ln_det_sigma_alpha) +
      -1/2 * sum(mapply(inv_sigma_alpha, vi_sigma_outer_alpha, FUN=function(a,b){sum(diag(a %*% b))}))
    
    #Get the terms for the final priors!
    logcomplete_3 <- 0 + #flat prior on beta
      sum(
        iw_prior_constant +
          -(prior_sigma_alpha_nu + d_j + 1)/2 * ln_det_sigma_alpha +
          -1/2 * mapply(prior_sigma_alpha_phi, inv_sigma_alpha, FUN=function(a,b){sum(diag(a %*% b))})
      )
    ##GET THE ENTROPY
    #Entropy for p(beta,alpha)
    if (factorization_method == 'weak'){
      entropy_1 <-  ncol(vi_joint_decomp)/2 * log(2 * pi * exp(1)) + 
        1/2 * log_det_joint_var 
    }else{
      entropy_1 <- ncol(vi_beta_decomp)/2 * log(2 * pi * exp(1)) + 1/2 * log_det_beta_var +
        ncol(vi_alpha_decomp)/2 * log(2 * pi * exp(1)) + 1/2 * log_det_alpha_var
    }
    #Entropy for Polya-Gamma EXCLUDING intractable term that cancels
    entropy_2 <- sum(vi_pg_b * vi_pg_c / 4 * tanh(vi_pg_c/2) - vi_pg_b * log(cosh(vi_pg_c/2)))
    #Entropy Wisharts
    entropy_3 <- -mapply(vi_sigma_alpha_nu, vi_sigma_alpha, FUN=function(nu,Phi){make_log_invwishart_constant(nu = nu, Phi = Phi)}) +
      (vi_sigma_alpha_nu + d_j + 1)/2 * ln_det_sigma_alpha +
      1/2 * mapply(vi_sigma_alpha, inv_sigma_alpha, FUN=function(a,b){sum(diag(a %*% b))})
    entropy_3 <- sum(entropy_3)
    
    logcomplete <- logcomplete_1 + logcomplete_2 + logcomplete_3
    logcomplete <- logcomplete + choose_term
    
  }else if (ELBO_type == 'profiled'){
    vi_r_var <- ( exp(vi_r_sigma) - 1 ) * vi_r_mean^2 
    
    logcomplete_1a <- sum(approx.lgamma(y, mean_r = vi_r_mean, var_r = vi_r_var)) +
      - N * approx.lgamma(x = 0, mean_r = vi_r_mean, var_r = vi_r_var) - (vi_r_mean) * N * log(2) +
      as.vector(t((y - vi_r_mean)/2) %*% ex_XBZA)
    logcomplete_1b <- sum(-(y + vi_r_mean) * log(cosh(1/2 * sqrt(ex_XBZA^2 + vi_r_sigma + var_XBZA))))
    logcomplete_1c <- N/2 * vi_r_mean * vi_r_sigma
    logcomplete_1 <- logcomplete_1a + logcomplete_1b + logcomplete_1c + choose_term
    
    logcomplete_2 <- sum(-d_j * g_j / 2 * log(2 * pi) - g_j/2 * ln_det_sigma_alpha) +
      -1/2 * sum(mapply(inv_sigma_alpha, vi_sigma_outer_alpha, FUN=function(a,b){sum(diag(a %*% b))}))
    logcomplete_3 <- 0 + #flat prior on beta
      sum(
        iw_prior_constant +
          -(prior_sigma_alpha_nu + d_j + 1)/2 * ln_det_sigma_alpha +
          -1/2 * mapply(prior_sigma_alpha_phi, inv_sigma_alpha, FUN=function(a,b){sum(diag(a %*% b))})
      )
    
    if (factorization_method == 'weak'){
      entropy_1 <-  ncol(vi_joint_decomp)/2 * log(2 * pi * exp(1)) + 
        1/2 * log_det_joint_var 
    }else{
      entropy_1 <- ncol(vi_beta_decomp)/2 * log(2 * pi * exp(1)) + 1/2 * log_det_beta_var +
        ncol(vi_alpha_decomp)/2 * log(2 * pi * exp(1)) + 1/2 * log_det_alpha_var
    }
    entropy_2 <- 0
    #Entropy Wisharts
    entropy_3 <- -mapply(vi_sigma_alpha_nu, vi_sigma_alpha, FUN=function(nu,Phi){make_log_invwishart_constant(nu = nu, Phi = Phi)}) +
      (vi_sigma_alpha_nu + d_j + 1)/2 * ln_det_sigma_alpha +
      1/2 * mapply(vi_sigma_alpha, inv_sigma_alpha, FUN=function(a,b){sum(diag(a %*% b))})
    entropy_3 <- sum(entropy_3)
    
    logcomplete <- logcomplete_1 + logcomplete_2 + logcomplete_3
  }else{stop('ELBO must be profiled or augmented')}
  
  entropy <- entropy_1 + entropy_2 + entropy_3
  ELBO <- entropy + logcomplete
  
  return(data.frame(ELBO, logcomplete, entropy, logcomplete_1, 
                    logcomplete_2, logcomplete_3, entropy_1, entropy_2, entropy_3))
}



update_r <- function(vi_r_mu, vi_r_sigma, y, X, Z, factorization_method,
                     vi_beta_mean, vi_alpha_mean, 
                     vi_joint_decomp, vi_beta_decomp, vi_alpha_decomp,
                     vi_r_method){
  
  if(vi_r_method == 'fixed'){
    return(c(vi_r_mu, vi_r_sigma))
  }
  #Get intermediate quantities
  ex_XBZA <- (X %*% vi_beta_mean + Z %*% vi_alpha_mean)
  #quadratic var, i.e. Var(x_i^T beta + z_i^T alpha)
  if (factorization_method == 'weak'){
    if (is.null(vi_joint_decomp)){stop('Need to provide joint decomposition for ELBO weak')}
    var_XBZA <- rowSums( (cbind(X,Z) %*% t(vi_joint_decomp))^2 )
  }else{
    beta_quad <- rowSums( (X %*% t(vi_beta_decomp))^2 )
    alpha_quad <- rowSums( (Z %*% t(vi_alpha_decomp))^2 )
    var_XBZA <- beta_quad + alpha_quad
  }
  
  N <- length(y)
  
  vi_r_mu <<- vi_r_mu
  vi_r_sigma <<- vi_r_sigma
  ex_XBZA <<- ex_XBZA
  var_XBZA <<- var_XBZA
  y <<- y
  N <<- N
  
  if (vi_r_method == 'VI'){
    
    
    opt_vi_r <- optim(par = vi_r_mu, fn = VEM.PELBO.r, y = y, psi = ex_XBZA, zVz = var_XBZA,
                      control = list(fnscale = -1), method = 'L-BFGS', hessian = T)
    
    prior_vi_r <- VEM.PELBO.r(ln_r = vi_r_mu, y = y,
                              psi = ex_XBZA, zVz = var_XBZA)
    
    
    if (opt_vi_r$value < prior_vi_r){
      warning('Optim for r decreased objective.')
      out_par <- c(vi_r_mu, vi_r_sigma)
    }else{
      vi_r_sigma <- as.numeric(-1/opt_vi_r$hessian)
      out_par <- c(opt_vi_r$par, vi_r_sigma)
    }
    
    # opt_vi_r <- optim(par = c(vi_r_mu, log(vi_r_sigma)), fn = PELBO.r,
    #                   y = y, psi = ex_XBZA, zVz = var_XBZA, N = N,
    #                   control = list(fnscale = - 1), method = 'L-BFGS')
    # prior_vi_r <- PELBO.r(par = c(vi_r_mu, log(vi_r_sigma)), y = y,
    #                       psi = ex_XBZA, zVz = var_XBZA, N = N)
    # if (opt_vi_r$value < prior_vi_r){
    #   warning('Optim for r decreased objective.')
    #   out_par <- c(vi_r_mu, vi_r_sigma)
    # }else{
    #   out_par <- c(opt_vi_r$par[1], exp(opt_vi_r$par[2]))
    # }
  }else if (vi_r_method == 'VEM'){
    opt_vi_r <- optim(par = vi_r_mu, fn = VEM.PELBO.r,
                      y = y, psi = ex_XBZA, zVz = var_XBZA, 
                      control = list(fnscale = - 1), method = 'L-BFGS')
    prior_vi_r <- VEM.PELBO.r(ln_r = vi_r_mu, y = y, 
                              psi = ex_XBZA, zVz = var_XBZA)
    
    if (opt_vi_r$value < prior_vi_r){
      warning('Optim for r decreased objective.')
      out_par <- c(vi_r_mu, 0)
    }else{
      out_par <- c(opt_vi_r$par, 0)
    }
  }else{stop('vi_r method must be VI or VEM or fixed')}
  
  return(out_par)
}

VEM.PELBO.r <- function(ln_r, y, psi, zVz){
  t1 <- -(y+exp(ln_r)) * log(2) + (y-exp(ln_r))/2 * (psi  - ln_r) +
    -(exp(ln_r) + y) * log(cosh(1/2 * sqrt(zVz + (psi - ln_r)^2)))
  t2 <- lgamma(y + exp(ln_r)) - lgamma(exp(ln_r)) - lgamma(y+1)
  return(sum(t1 + t2))
}

approx.lgamma <- function(x, mean_r, var_r){
  input <- x + mean_r
  output <- lgamma(input) + 1/2 * psigamma(x = input, deriv = 1) * var_r
  return(output)
}

PELBO.r <- function(par, y, psi, zVz, N){
  mu_r <- par[1]
  log_sigma_r <- par[2]
  
  sigma_r <- exp(log_sigma_r)
  mean_r <- exp(mu_r + sigma_r/2)
  var_r <- ( exp(sigma_r) - 1 ) * mean_r^2
  cov_r <- mean_r * sigma_r
  entropy_r <- mu_r + 1/2 * log_sigma_r
  
  #Gamma Approximation Term
  t1 <- sum( approx.lgamma(x = y, mean_r = mean_r, var_r = var_r) ) +
    - N * approx.lgamma(x = 0, mean_r = mean_r, var = var_r) +
    - N * log(2) * mean_r
  #E[s^T (psi - ln r)] ....
  t2 <- sum( (y - mean_r)/2 * (psi - mu_r) ) + N/2 * cov_r
  
  t3 <- sum(
    -(y + mean_r) * log(cosh(1/2 * sqrt(zVz + (psi - mu_r)^2 + sigma_r)))
  )
  
  return(t1 + t2 + t3 + entropy_r)
}
