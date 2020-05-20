
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

#' Format Existing Objects
#' 
#' Takes a regression output from glmer, stan, or vglmer and extracts the fixed
#' and random effects.
#'
#' @name format_obj
#' 
#' @param object Object from glmer, vglmer, or stan to the respective "format_"
#'   function
#' @importFrom dplyr bind_rows mutate
#' @importFrom lme4 fixef ranef
#' @importFrom stats vcov
#' @importFrom rlang .data
#' @export
format_glmer <- function(object){
  output <- bind_rows(lapply(ranef(object), FUN=function(i){
    obj <- data.frame(var = as.vector(apply(attributes(i)$postVar, MARGIN = 3, FUN=function(i){diag(i)})),
                      mean = as.vector(t(as.matrix(i))),
                      name = paste0(rep(colnames(i), nrow(i)), ' @ ', rep(rownames(i), each = ncol(i))), stringsAsFactors = F)
    return(obj)
  }), .id = '.re')
  output$name <- with(output, paste0(.re, ' @ ', name))
  output_fe <- data.frame(mean = fixef(object), var = diag(vcov(object)))
  output_fe$name <- rownames(output_fe)
  
  output <- bind_rows(output, output_fe)
  output <- output[, (names(output) != '.re')]
  
  rownames(output) <- NULL
  
  return(output)
}

#' Format vglmer
#' @rdname format_obj
#' @export
format_vglmer <- function(object){
  beta.output <- data.frame(name = rownames(object$beta$mean), mean = as.vector(object$beta$mean), var = diag(object$beta$var), stringsAsFactors = F)
  alpha.output <- data.frame(name = rownames(object$alpha$mean), mean = as.vector(object$alpha$mean), var = as.vector(object$alpha$dia.var), stringsAsFactors = F)
  output <- bind_rows(beta.output, alpha.output)
  return(output)
}

#' Format Stan
#' @param useSigma Return variance component parameters from STAN? Default "FALSE".
#' @rdname format_obj
#' @importFrom stringr str_split
#' @importFrom stats var
#' @export
format_stan <- function(object, useSigma = FALSE){
  post_stan <- as.matrix(object)
  if (!useSigma){
    post_stan <- post_stan[,!grepl(colnames(post_stan), pattern='^Sigma')]
  }
  parse_stan_names <- str_split(colnames(post_stan), pattern='^b\\[| |\\]')

  fmt_stan_names <- sapply(parse_stan_names, FUN=function(i){
    if (length(i) == 1){
      return(i)
    }else{
      i_one <- unlist(str_split(i[3], pattern=':'))
      return(paste(i_one[1], i[2], i_one[2], sep=' @ '))
    }
  })
  colnames(post_stan) <- fmt_stan_names
  output <- data.frame(var = apply(post_stan, MARGIN = 2, var), mean = colMeans(post_stan))
  output$name <- rownames(output)
  rownames(output) <- NULL
  return(output)
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

loop_outer_alpha <- function(vi_alpha_mean, vi_alpha_decomp, outer_alpha_RE_positions){
  #Must be such that t(vi_alpha_decomp) %*% vi_alpha_decomp = VAR
  store_oa <- as.list(rep(NA, length(outer_alpha_RE_positions)))
  store_alpha_outer <- as.list(rep(NA, length(outer_alpha_RE_positions)))
  counter_j <- 1
  for (j in outer_alpha_RE_positions){
    #cat('.')
    summed_oa <- summed_alpha_outer <- array(0, dim = rep(length(j[[1]]), 2))
    for (g in j){
      summed_oa <- summed_oa + crossprod(vi_alpha_decomp[,g])      
      summed_alpha_outer <- summed_alpha_outer + tcrossprod(vi_alpha_mean[g])
    }
    store_oa[[counter_j]] <- summed_oa
    store_alpha_outer[[counter_j]] <- summed_alpha_outer
    counter_j <- counter_j + 1
  }
  return(list(outer_alpha = store_oa, alpha_mu_outer = store_alpha_outer))
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
  EM_beta <- rep(1/ncol(jointXZ), ncol(jointXZ))
  EM_variance <- sparseMatrix(i = 1:ncol(jointXZ), j = 1:ncol(jointXZ), x = 1/ridge)
  
  for (it in 1:iter){
    
    EM_pg_c <- jointXZ %*% EM_beta
    EM_pg_mean <- pg_b/(2 * EM_pg_c) * tanh(EM_pg_c / 2)
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

calculate_ELBO <- function(factorization_method,
   #Fixed constants or priors
   d_j, g_j, prior_sigma_alpha_phi, prior_sigma_alpha_nu, iw_prior_constant, choose_term,
   #Data
   X, Z, s, 
   #PolyaGamma Parameters
   vi_pg_b, vi_pg_mean, vi_pg_c,
   #Sigma Parameters
   vi_sigma_alpha, vi_sigma_alpha_nu, vi_sigma_outer_alpha,
   #Beta Parameters / Alpha Parameters
   vi_beta_mean, vi_beta_decomp,
   vi_alpha_mean, vi_alpha_decomp,
   log_det_beta_var, log_det_alpha_var, 
   log_det_joint_var = NULL,
   vi_joint_decomp = NULL){
  ####
  ##PREPARE INTERMEDIATE QUANTITES
  ###
    #linear predictor: E[XB + ZA]
    ex_XBZA <- (X %*% vi_beta_mean + Z %*% vi_alpha_mean)
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
  logcomplete <- logcomplete_1 + logcomplete_2 + logcomplete_3
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
  
  logcomplete <- logcomplete + choose_term
  
  entropy <- entropy_1 + entropy_2 + entropy_3
  ELBO <- entropy + logcomplete
  return(data.frame(ELBO, logcomplete, entropy, logcomplete_1, logcomplete_2, logcomplete_3, entropy_1, entropy_2, entropy_3))
}

#' Variational Inference for Non-Linear Hierarchical Models
#' 
#' Estimate a hierarchical model (logistic) using mean-field variational
#' inference. Accepts identical syntax to glmer. Options are described below.
#' 
#' @param formula Standard glmer-style formula for random effects.
#' @param data data.frame containing the outcome and variables.
#' @param iterations Number of iterations for the model.
#' @param factorization_method The factorization method to use. Described in
#'   detail in the dissertation. strong, partial, and weak correspond to Schemes
#'   I, II, and III respectively. "weak" should have best performance but is
#'   slowest.
#' @param prior_variance Options are jeffreys, mcmcglmm, mvD, mean_exists,
#'   limit, and uniform. They are defined as an Inverse Wishart with the
#'   following parameters where d is the dimensionality of the random effect. At
#'   present, all models use "mean_exists" for a well-defined proper prior. 
#'   \itemize{
#'   \item jeffreys: IW(0, 0)
#'   \item mcmcglmm: IW(0, I)
#'   \item mvD: IW(-d, I)
#'   \item mean_exists: IW(d + 1, I)
#'   \item limit: IW(d - 1, 0)
#'   \item uniform: IW(-[d+1], 0)
#'   }
#'   The model may fail to converge   
#' @param tolerance_elbo Change in ELBO to stop algorithm.
#' @param tolerance_parameters Change in value of any parameter to stop algorithm.
#' 
#' @param parameter_expansion At moment, accepts 'mean' or 'none'. 'mean' is
#'   costless and should always be used!
#' @param prevent_degeneracy Ignored for the moment.
#' @param force_whole Require whole numbers. Set to FALSE to allow "quasi-binomial".
#' 
#' @param init Initialization method can be EM, zero, or random at moment.
#' @param random_seed Set seed for initialization.
#' 
#' @param debug_param Debug parameter convergence.
#' @param debug_ELBO Debug ELBO trajectory.
#' 
#' @param linpred_method Method for updating means of beta and alpha. "joint" is best.
#' @param print_prog Print after print_prog iterations to show progress.
#' @param quiet Don't print noisy intermediate output.
#' @param return_data Return the design (X,Z) for debugging afterwards.
#' 
#' @import CholWishart purrr stringr
#' @importFrom dplyr select group_by group_by_at summarize
#' @importFrom lme4 mkReTrms findbars subbars
#' @importFrom stats model.response model.matrix model.frame rnorm rWishart qlogis
#' @importFrom graphics plot
#' @importFrom rlang .data
#' @useDynLib vglmer
#' @export
vglmer_logit <- function(formula, data, iterations, prior_variance, factorization_method = 'weak', 
     tolerance_elbo = 1e-4, tolerance_parameters = 1e-4,
     prevent_degeneracy = FALSE, force_whole = TRUE,
     parameter_expansion = 'mean', random_seed = 1,
     outcome = 'logit',
     debug_param = FALSE, return_data = FALSE, linpred_method = 'joint',
     debug_ELBO = FALSE, print_prog = NULL, quiet = T, init = 'EM'){
  
  if (!(factorization_method %in% c('weak', 'strong', 'partial'))){
    stop('factorization_method must be in weak, strong, or partial.')
  }
  if (is.null(print_prog)){
    print_prog <- max(c(1, floor(iterations / 20)))
  }
  #Extract outcome y
  y <- model.response(model.frame(nobars(formula), data = data))
  names(y) <- NULL
  if (!(class(y) %in% c('numeric', 'integer'))){
    if (min(y) < 0){
      stop('Negative Numbers not Permitted in outcome')
    }
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if (any(is.wholenumber(y) == FALSE)){
      if (force_whole){
        stop('If force_whole = TRUE, must provide whole numbers')
      }else{
        warning('Non-integer numbers in y')
      }
    }
    #Total trials (Success + Failure)
    trials <- rowSums(y)
    rownames(trials) <- NULL
    #Successes
    y <- y[,1]
    rownames(y) <- NULL
    
  }else{
    if (!all(y %in% c(0,1)) & outcome == 'logit'){
      stop('Only {0,1} outcomes permitted for numeric y.')
    }
    trials <- rep(1, length(y))
  }
  N <- length(y)
  
  #Extract X (FE design matrix)
  X <- model.matrix(nobars(formula), data = data)
  
  #Extract the Z (Random Effect) design matrix.
  mk_Z <- mkReTrms(findbars(formula), model.frame(subbars(formula), data = data), reorder.terms = FALSE, reorder.vars = FALSE)
  Z <- t(mk_Z$Zt)

  p.X <- ncol(X)
  p.Z <- ncol(Z)
  
  ####
  #Process the REs to get various useful terms.
  ####
    #RE names and names of variables included for each.
    names_of_RE <- mk_Z$cnms
  
    number_of_RE <- length(mk_Z$Gp) - 1
    
    if (number_of_RE < 1){stop('Need to provide at least one random effect...')}
    #The position that demarcates each random effect.
    #That is, breaks_for_RE[2] means at that position + 1 does RE2 start.
    breaks_for_RE <- c(0, cumsum(diff(mk_Z$Gp)))
    #Dimensionality of \alpha_{j,g}, i.e. 1 if random intercept
    #2 if random intercept + random slope
    d_j <- lengths(names_of_RE)
    #Number of GROUPs for each random effect.
    g_j <- diff(mk_Z$Gp)/d_j
    
    #Empty vector to build the formatted names for each random effect.
    fmt_names_Z <- c()
    init_Z_names <- colnames(Z)
    for (v in 1:number_of_RE){
      name_of_effects_v <- names_of_RE[[v]]
      
      mod_name <- rep(name_of_effects_v, g_j[v])
      
      levels_of_re <- init_Z_names[(1+breaks_for_RE[v]):breaks_for_RE[v+1]]
      
      fmt_names_Z <- c(fmt_names_Z, paste0(names(names_of_RE)[v], ' @ ', mod_name, ' @ ', levels_of_re))
      
    }
    colnames(Z) <- fmt_names_Z
    cyclical_pos <- lapply(1:number_of_RE, FUN=function(i){seq(breaks_for_RE[i] + 1, breaks_for_RE[i+1])})
    
    M.names <- cbind(unlist(mapply(names_of_RE, g_j, FUN=function(i,j){rep(i,j)})))
    M <- cbind(match(M.names[,1], colnames(X)), rep(1/g_j, d_j * g_j))
    M <- sparseMatrix(i = 1:ncol(Z), j = M[,1], x = M[,2], dims = c(ncol(Z), ncol(X)))
    
    M_prime.names <- paste0(rep(names(names_of_RE), g_j * d_j), ' @ ', M.names)
    M_prime <- cbind(match(M_prime.names, unique(M_prime.names)), rep(1/g_j, d_j * g_j))
    M_prime <- sparseMatrix(i = 1:ncol(Z), j = M_prime[,1], x = M_prime[,2])
    colnames(M_prime) <- unique(M_prime.names)
  
    M_prime_one <- M_prime
    M_prime_one@x <- rep(1, length(M_prime_one@x))
    
    stopifnot(identical(paste0(rep(names(names_of_RE), d_j), ' @ ', unlist(names_of_RE)), colnames(M_prime)))
    
    M_mu_to_beta <- sparseMatrix(i = 1:sum(d_j), j = match(unlist(names_of_RE), colnames(X)), x = 1, dims = c(sum(d_j), p.X))
    colnames(M_mu_to_beta) <- colnames(X)
    rownames(M_mu_to_beta) <- colnames(M_prime)
    #List of Lists
    #Outer list: one for RE
    #Inner List: One for each GROUP with its row positions.
    outer_alpha_RE_positions <- mapply(d_j, g_j, breaks_for_RE[-length(breaks_for_RE)], SIMPLIFY = FALSE, FUN=function(a,b,m){
      split(m + seq(1, a * b), rep(1:b, each = a))
    })
    
    if (anyDuplicated(unlist(outer_alpha_RE_positions)) != 0 | max(unlist(outer_alpha_RE_positions)) != ncol(Z)){
      stop('Issue with greating OA positions')
    }
  ####
  #Prepare Initial Values
  ###

    if (outcome == 'logit'){
      s <- y - trials/2
      vi_pg_b <- trials
      vi_r <- 1
      choose_term <- sum(lchoose(n = round(trials), k = round(y)))
    }else{
      #Initialize
      vi_r <- 1/MASS::glm.nb(y ~ 1)$theta
      s <- (y - vi_r)/2
      vi_pg_b <- y + vi_r
      choose_term <- 0
    }

    #Initalize variational parameters.
    #Note that we keep a sparse matrix or lowertri such that
    #t(vi_beta_decomp) %*% vi_beta_decomp = VARIANCE

    vi_beta_decomp <- Diagonal(x = rep(0, ncol(X)))
    vi_alpha_decomp <- Diagonal(x = rep(0, ncol(Z)))
    
    vi_sigma_alpha_nu <-  g_j
    
    if (prior_variance == 'jeffreys'){
      prior_sigma_alpha_nu <- rep(0, number_of_RE)
      prior_sigma_alpha_phi <- lapply(d_j, FUN=function(i){diag(x = 0, nrow = i, ncol = i)})
    }else if (prior_variance == 'mcmcglmm'){
      prior_sigma_alpha_nu <- rep(0, number_of_RE)
      prior_sigma_alpha_phi <- lapply(d_j, FUN=function(i){diag(x = 0, nrow = i, ncol = i)})
    }else if (prior_variance == 'mvD'){
      prior_sigma_alpha_nu <- -d_j
      prior_sigma_alpha_phi <- lapply(d_j, FUN=function(i){diag(x = 0, nrow = i, ncol = i)})
    }else if (prior_variance == 'mean_exists'){
      prior_sigma_alpha_nu <- d_j + 1 #Ensures the mean exists...
      prior_sigma_alpha_phi <- lapply(d_j, FUN=function(i){diag(x = 1, nrow = i, ncol = i)})
    }else if (prior_variance == 'limit'){
      prior_sigma_alpha_nu <- d_j - 1
      prior_sigma_alpha_phi <- lapply(d_j, FUN=function(i){diag(x = 0, nrow = i, ncol = i)})
    }else if (prior_variance == 'uniform'){
      prior_sigma_alpha_nu <- -(d_j + 1)
      prior_sigma_alpha_phi <- lapply(d_j, FUN=function(i){diag(x = 0, nrow = i, ncol = i)})
    }else{
      stop('Options for prior variance are jeffreys and mean_exists')
    }

    #normalizingly constant for wishart to make ELBO have right value to compare models.
    iw_prior_constant <- mapply(prior_sigma_alpha_nu, prior_sigma_alpha_phi, 
        FUN=function(nu,Phi){
          if (nu <= (ncol(Phi) - 1)){
            return(0)
          }else{
            return(make_log_invwishart_constant(nu,Phi))
          }
        })
    
    
    vi_sigma_alpha_nu <- vi_sigma_alpha_nu + prior_sigma_alpha_nu
        
    if (init == 'EM'){
      if (outcome == 'negbin'){stop('EM init not yet enabled for NB.')}
      EM_init <- EM_prelim(X = X, Z = Z, s = s, pg_b = vi_pg_b, iter = 15, ridge = 4)

      vi_beta_mean <- matrix(EM_init$beta)
      vi_alpha_mean <- matrix(EM_init$alpha)
      
      vi_sigma_alpha <- calculate_expected_outer_alpha(alpha_mu = matrix(EM_init$alpha), 
        L = sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(Z), 2)), 
        re_position_list = outer_alpha_RE_positions)
      vi_sigma_alpha <- mapply(vi_sigma_alpha$outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN=function(i,j){i+j})
    }else if (init == 'random'){
      set.seed(random_seed)
      vi_beta_mean <- rnorm(ncol(X))
      vi_alpha_mean <- rep(0, ncol(Z))
      
      vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN=function(d,g){
        rWishart(n = 1, df = g, Sigma = diag(d))[,,1]
      })
    }else if (init == 'zero'){
      
      vi_beta_mean <- rep(0, ncol(X))
      
      if (outcome == 'logit'){
        vi_beta_mean[1] <- qlogis(sum(y)/sum(trials))
      }else if (outcome == 'negbin'){
        vi_beta_mean[1] <- log(mean(y))
      }else{stop('Set up init')}
      
      vi_alpha_mean <- rep(0, ncol(Z))
      
      vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN=function(d,g){
        diag(x = 1, ncol = d, nrow = d)
      })
      
    }else{stop('Provide init = EM or random')}
    
    zero_mat <- sparseMatrix(i=1,j=1,x=0,dims = c(ncol(X), ncol(X)))
    zero_mat <- drop0(zero_mat)
    
    if (factorization_method == 'weak'){
      vi_joint_decomp <- bdiag(vi_beta_decomp, vi_alpha_decomp)
      joint.XZ <- cbind(X, Z)
    }else{
      vi_joint_decomp = NULL
      log_det_joint_var <- NULL
    }
    
    if (!quiet){warning('Check vi sigma alpha nu')}
    
    #Create mapping for this to allow sparse implementations.
    mapping_sigma_alpha <- make_mapping_alpha(vi_sigma_alpha)

    running_log_det_alpha_var <- rep(NA, number_of_RE)
    
    lagged_alpha_mean <- rep(-Inf, ncol(Z))
    lagged_beta_mean <- rep(-Inf, ncol(X))
    lagged_sigma_alpha <- vi_sigma_alpha
    if (factorization_method == 'weak'){
      lagged_joint_decomp <- vi_joint_decomp
    }else{
      lagged_alpha_decomp <- vi_alpha_decomp
      lagged_beta_decomp <- vi_beta_decomp
    }
    lagged_vi_r <- -Inf
    
    lagged_ELBO <- -Inf
    accepted_times <- NA
    
    if (parameter_expansion == 'translation'){
      stop('parameter_expansion_translation not allowed yet.')
      accepted_times <- 0
      zeromat_beta <- drop0(Diagonal(x = rep(0, ncol(X))))
      
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
      
      start_base_Z <- cumsum(c(0,d_j))[-(number_of_RE+1)]
      names(start_base_Z) <- NULL
      
      rm(parsed_RE_groups, mapping_to_re)
    }
    store_ELBO <- data.frame()
    
    if (debug_param){
      store_beta <- array(NA, dim = c(iterations, ncol(X)))
    }
    ##Begin VI algorithm:
    if (!quiet){message('Begin Regression')}
    for (it in 1:iterations){
      if (it %% print_prog == 0){cat('.')}
      ###
      ##Polya-Gamma Updates
      ###
      #Get the x_i^T Var(beta) x_i terms.
      log_vi_r <- log(vi_r)
      if (factorization_method == 'weak'){
        joint_quad <- rowSums( (joint.XZ %*% t(vi_joint_decomp))^2 )
        vi_pg_c <- sqrt(as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean - log_vi_r)^2 + joint_quad)
      }else{
        beta_quad <- rowSums( (X %*% t(vi_beta_decomp))^2 )
        alpha_quad <- rowSums( (Z %*% t(vi_alpha_decomp))^2 )
        vi_pg_c <- sqrt(as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean - log_vi_r)^2 + beta_quad + alpha_quad)
      }
      vi_pg_mean <- vi_pg_b/(2 * vi_pg_c) * tanh(vi_pg_c / 2)

      diag_vi_pg_mean <- sparseMatrix(i = 1:N, j = 1:N, x = vi_pg_mean)
      
      if (debug_ELBO & it != 1){
        debug_ELBO.1 <- calculate_ELBO(
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = vi_sigma_outer_alpha,
          vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
          vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
          log_det_joint_var = log_det_joint_var
        )
      }

      #Process Sigma_j for manipulation      
      #if Sigma_{j} is InverseWishart(a,Phi)
      #Then E[Sigma^{-1}_j] = a * Phi^{-1}
      if (factorization_method == 'strong'){
        cyclical_T <- TRUE
      }else{
        cyclical_T <- FALSE
      }
      inv_mapping_alpha <- mapply(vi_sigma_alpha_nu, lapply(vi_sigma_alpha, solve), SIMPLIFY = FALSE, FUN=function(a,b){a * b})
      inv_mapping_alpha <- make_mapping_alpha(inv_mapping_alpha)
      Tinv <- prepare_T(mapping = inv_mapping_alpha, levels_per_RE = g_j, num_REs = number_of_RE,
                        variables_per_RE = d_j, running_per_RE = breaks_for_RE, cyclical = cyclical_T)
      if (!cyclical_T){
        Tinv <- as(Tinv, 'dgCMatrix')
      }else{
        Tinv <- lapply(Tinv, FUN=function(i){as(i, 'dgCMatrix')})  
      }
      
      if (factorization_method == 'weak'){
        ##Update <beta, alpha> jointly
        chol.update.joint <- LinRegChol(X = joint.XZ, omega = diag_vi_pg_mean, 
                                        prior_precision =  bdiag(zero_mat, Tinv), 
                                        y = s + vi_pg_mean * log(vi_r))
        Pmatrix <- sparseMatrix(i = 1:ncol(joint.XZ), j = 1 + chol.update.joint$Pindex, x = 1)
        
        
        vi_joint_decomp <- drop0(solve(chol.update.joint$origL) %*% t(Pmatrix))
        
        vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))
        
        vi_alpha_decomp <- vi_joint_decomp[,-1:-p.X, drop = F]
        vi_beta_decomp <- vi_joint_decomp[,1:p.X, drop = F]
        
        log_det_joint_var <- - 2 * sum(log(diag(chol.update.joint$origL)))
        
      }else if (factorization_method == 'partial'){
        
        if (linpred_method == 'cyclical'){
          # Do not run except as backup
          # ###Non optimized
          # precision_beta <- t(X) %*% diag_vi_pg_mean %*% X
          # nonopt_beta <- solve(precision_beta, t(X) %*% (s - diag_vi_pg_mean %*% Z %*% vi_alpha_mean))
          # precision_alpha <- t(Z) %*% diag_vi_pg_mean %*% Z + Tinv
          # nonopt_alpha <- solve(precision_alpha, t(Z) %*% (s - diag_vi_pg_mean %*% X %*% nonopt_beta))
          
          chol.update.beta <- LinRegChol(X = as(X, 'sparseMatrix'), omega = diag_vi_pg_mean, prior_precision =  zero_mat, 
                                         y = as.vector(s - diag_vi_pg_mean %*% Z %*% vi_alpha_mean))
          Pmatrix <- sparseMatrix(i = 1:p.X, j = 1 + chol.update.beta$Pindex, x = 1)
          
          #P origL oriL^T P^T = PRECISION
          #t(decompVar) %*%  decompVar = VARIANCE = (origL^{-1} t(P))^T (origL^{-1} t(P))
          
          vi_beta_decomp <- solve(chol.update.beta$origL) %*% t(Pmatrix)
          vi_beta_mean <- chol.update.beta$mean
          log_det_beta_var <- - 2 * sum(log(diag(chol.update.beta$origL)))
          
          chol.update.alpha <- LinRegChol(X = Z, omega = diag_vi_pg_mean, prior_precision =  Tinv, 
                                          y = as.vector(s - diag_vi_pg_mean %*% X %*% vi_beta_mean))
          Pmatrix <- sparseMatrix(i = 1:p.Z, j = 1 + chol.update.alpha$Pindex, x = 1)
          
          vi_alpha_decomp <- solve(chol.update.alpha$origL) %*% t(Pmatrix)
          vi_alpha_decomp <- drop0(vi_alpha_decomp)
          vi_alpha_mean <- chol.update.alpha$mean
          log_det_alpha_var <- - 2 * sum(log(diag(chol.update.alpha$origL)))
          
          vi_beta_mean <- Matrix(vi_beta_mean, dimnames = list(colnames(X), NULL))
          vi_alpha_mean <- Matrix(vi_alpha_mean, dimnames = list(fmt_names_Z, NULL))
          
        }else if (linpred_method == 'joint'){
          joint.XZ <- cbind(X,Z)  
          chol.update.joint <- LinRegChol(X = joint.XZ, omega = diag_vi_pg_mean, prior_precision =  bdiag(zero_mat, Tinv),  
                                          y= s + vi_pg_mean * log(vi_r))
          vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
          vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))
          
          vi_beta_decomp <- solve(t(chol(as.matrix(t(X) %*% diag_vi_pg_mean %*% X))))
          log_det_beta_var <- 2 * sum(log(diag(vi_beta_decomp)))
          
          chol.update.alpha <- LinRegChol(X = Z, omega = diag_vi_pg_mean, prior_precision =  Tinv, 
                                          y =  s + vi_pg_mean * log(vi_r))
          Pmatrix <- sparseMatrix(i = 1:p.Z, j = 1 + chol.update.alpha$Pindex, x = 1)
          
          vi_alpha_decomp <- solve(chol.update.alpha$origL) %*% t(Pmatrix)
          vi_alpha_decomp <- drop0(vi_alpha_decomp)
          log_det_alpha_var <- - 2 * sum(log(diag(chol.update.alpha$origL)))
          
        }else{stop('Invalid linpred method for partial scheme')}
        
      }else if (factorization_method == 'strong'){
        
        running_log_det_alpha_var <- rep(NA, number_of_RE)
        vi_alpha_decomp <- sparseMatrix(i = 1, j =1, x = 0, dims = rep(p.Z, 2))
        
        # 
        # vi_alpha_mean <<- vi_alpha_mean
        # vi_beta_mean <<- vi_beta_mean
        # vi_alpha_decomp <<- vi_alpha_decomp
        # X <<- X
        # Z <<- Z
        # fmt_names_Z <<- fmt_names_Z
        # Tinv <<- Tinv
        # cyclical_pos <<- cyclical_pos
        # outcome_s <<- s
        # diag_vi_pg_mean <<- diag_vi_pg_mean
        # running_log_det_alpha_var <<- running_log_det_alpha_var
        # number_of_RE <<- number_of_RE
        # zero_mat <<- zero_mat
        # p.X <<- p.X
        
        
        if (linpred_method == 'joint'){
          
          joint.XZ <- cbind(X,Z)  
          chol.update.joint <- LinRegChol(X = joint.XZ, omega = diag_vi_pg_mean, prior_precision =  bdiag(zero_mat, bdiag(Tinv)),  y = s + vi_pg_mean * log(vi_r))
          vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
          vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))
          
          vi_beta_decomp <- solve(t(chol(as.matrix(t(X) %*% diag_vi_pg_mean %*% X))))
          log_det_beta_var <- 2 * sum(log(diag(vi_beta_decomp)))
            #-log(det(t(X) %*% diag_vi_pg_mean %*% X))
          
          running_log_det_alpha_var <- rep(NA, number_of_RE)
          
          for (j in 1:number_of_RE){
            
            index_j <- cyclical_pos[[j]]
            Z_j <- Z[,index_j, drop = F]
            prec_j <- t(Z_j) %*% diag_vi_pg_mean %*% Z_j + Tinv[[j]]
            
            chol_var_j <- solve(t(chol(prec_j)))
            running_log_det_alpha_var[j] <- 2 * sum(log(diag(chol_var_j)))
            
            vi_alpha_decomp[index_j, index_j] <- as(chol_var_j, 'dgTMatrix')
            
          }
          
          log_det_alpha_var <- sum(running_log_det_alpha_var)
          
        }else if (linpred_method == 'solve_normal'){
          
          bind_rhs_j <- list()
          bind_lhs_j <- list()
          
          for (j in 1:number_of_RE){
            
            index_j <- cyclical_pos[[j]]
            Z_j <- Z[,index_j, drop = F]
            Z_negj <- Z[,-index_j, drop= F]
            prec_j <- t(Z_j) %*% diag_vi_pg_mean %*% Z_j + Tinv[[j]]
            
            chol_prec_j <- t(chol(prec_j))
            chol_var_j <- solve(chol_prec_j)
            
            mod_j <- solve(prec_j)
            
            term_j <- mod_j %*% t(Z_j) %*% diag_vi_pg_mean %*% Z
            term_j[,index_j, drop = F] <- Diagonal(n = ncol(Z_j))
            term_j <- cbind(term_j, mod_j %*% t(Z_j) %*% diag_vi_pg_mean %*% X)
            
            bind_lhs_j[[j]] <- term_j
            bind_rhs_j[[j]] <- mod_j %*% t(Z_j) %*% s
            
            running_log_det_alpha_var[j] <- 2 * sum(log(diag(chol_var_j)))
            vi_alpha_decomp[index_j, index_j] <- as(chol_var_j, 'dgTMatrix')
            
          }
          
          log_det_alpha_var <- sum(running_log_det_alpha_var)
          
          bind_lhs_j <- drop0(do.call('rbind', bind_lhs_j))
          bind_rhs_j <- do.call('rbind', bind_rhs_j)
          
          vi_beta_decomp <- solve(t(chol(as.matrix(t(X) %*% diag_vi_pg_mean %*% X))))
          vi_beta_var <- solve(t(X) %*% diag_vi_pg_mean %*% X)
          log_det_beta_var <- 2 * sum(log(diag(vi_beta_decomp)))
          
          # vi_beta_mean <- vi_beta_var %*% t(X) %*% (s - diag_vi_pg_mean %*% Z %*% vi_alpha_mean)
          # vi_alpha_mean <- solve(bind_lhs_j[,1:ncol(Z)], bind_rhs_j)
          # 
          # vi_alpha_mean <- Matrix(vi_alpha_mean)
          # vi_beta_mean <- Matrix(vi_beta_mean)
          
          bind_lhs_j <- drop0(rbind(bind_lhs_j, cbind(vi_beta_var %*% t(X) %*% diag_vi_pg_mean %*% Z, Diagonal(n = ncol(X)))))
          bind_rhs_j <- rbind(bind_rhs_j, vi_beta_var %*% t(X) %*% s)
          # 
          bind_solution <- solve(bind_lhs_j) %*% bind_rhs_j
          #print(cbind(bind_solution, rbind(vi_alpha_mean, vi_beta_mean)))
          # 
          vi_beta_mean <- Matrix(bind_solution[-1:-ncol(Z)], dimnames = list(colnames(X), NULL))
          vi_alpha_mean <- Matrix(bind_solution[1:ncol(Z)], dimnames = list(fmt_names_Z, NULL))
          
          
        }else if (linpred_method == 'cyclical'){
          
          for (j in 1:number_of_RE){
            index_j <- cyclical_pos[[j]]
            Z_j <- Z[,index_j, drop = F]
            Z_negj <- Z[,-index_j, drop= F]
            
            chol.j <- LinRegChol(X = Z_j, omega = diag_vi_pg_mean, prior_precision =  Tinv[[j]], 
                                 y = as.vector(s + vi_pg_mean * log(vi_r) - diag_vi_pg_mean %*% (X %*% vi_beta_mean + Z_negj %*% vi_alpha_mean[-index_j])))
            vi_alpha_mean[index_j] <- chol.j$mean 
            
            Pmatrix <- sparseMatrix(i = 1:ncol(Z_j), j = 1 + chol.j$Pindex, x = 1)
            
            running_log_det_alpha_var[j] <- - 2 * sum(log(diag(chol.j$origL)))
            vi_alpha_decomp[index_j, index_j] <- solve(chol.j$origL) %*% t(Pmatrix)
          }
          
          #vi_alpha_decomp <- bdiag(vi_alpha_decomp)
          log_det_alpha_var <- sum(running_log_det_alpha_var)
          
          chol.update.beta <- LinRegChol(X = as(X, 'sparseMatrix'), omega = diag_vi_pg_mean, prior_precision =  zero_mat, 
                                         y = as.vector(s + vi_pg_mean * log(vi_r) - diag_vi_pg_mean %*% Z %*% vi_alpha_mean))
          Pmatrix <- sparseMatrix(i = 1:p.X, j = 1 + chol.update.beta$Pindex, x = 1)
          
          vi_beta_decomp <- solve(chol.update.beta$origL) %*% t(Pmatrix)
          vi_beta_mean <- chol.update.beta$mean
          log_det_beta_var <- - 2 * sum(log(diag(chol.update.beta$origL)))
          
          vi_beta_mean <- Matrix(vi_beta_mean, dimnames = list(colnames(X), NULL))
          vi_alpha_mean <- Matrix(vi_alpha_mean, dimnames = list(fmt_names_Z, NULL))
        }else{stop('Invalid linpred method')}
        


      }else{stop('Invalid factorization method.')}

      if (debug_ELBO & it != 1){
        variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
        vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha

        debug_ELBO.2 <- calculate_ELBO(
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = vi_sigma_outer_alpha,
          vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
          vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
          log_det_joint_var = log_det_joint_var
        )
      }
      ###
      # Update \Sigma_j
      ##
      
      variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
      vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
      vi_sigma_alpha <- mapply(vi_sigma_outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN=function(i,j){i+j})
      # 
      # # #Slow non-optimized function:
      # old_expectations_alpha_outer <- loop_outer_alpha(vi_alpha_mean, vi_alpha_decomp, outer_alpha_RE_positions)
      # old_vi_sigma_outer_alpha <- mapply(old_expectations_alpha_outer[[1]], old_expectations_alpha_outer[[2]], SIMPLIFY = FALSE, FUN=function(a,b){a+b})
      # old_vi_sigma_alpha <- mapply(old_vi_sigma_outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN=function(i,j){i+j})
      # 
      # vi_sigma_alpha <- old_vi_sigma_alpha
      # vi_sigma_outer_alpha <- old_vi_sigma_outer_alpha
      # 
      # comp <- mapply(vi_sigma_outer_alpha, old_vi_sigma_outer_alpha, FUN=function(i,j){all.equal(as.vector(i), as.vector(j))})
      # if (any(!sapply(comp, isTRUE))){print(comp); stop()}
      
      #Update the auxilary parameters
      if (outcome == 'negbin'){
      
        vi_r <- update_r(vi_r = vi_r, y = y, X = X, Z = Z, factorization_method = factorization_method,
                 vi_beta_mean = vi_beta_mean, vi_beta_decomp = vi_beta_decomp,
                 vi_alpha_mean = vi_alpha_mean, vi_alpha_decomp = vi_alpha_decomp,
                 vi_joint_decomp = vi_joint_decomp)
        
        s <- (y - vi_r)/2
        vi_pg_b <- y + vi_r
        
      }
      
      ###PARAMETER EXPANSIONS!
      if (debug_ELBO){
        debug_ELBO.3 <- calculate_ELBO(
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = vi_sigma_outer_alpha,
          vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
          vi_joint_decomp = vi_joint_decomp,
          log_det_joint_var = log_det_joint_var, choose_term = choose_term
        )
      }
      
      if (parameter_expansion == 'mean'){
        #Do a simple mean adjusted expansion.
        
        #Get the mean of each random effect.
        vi_mu_j <- t(M_prime) %*% vi_alpha_mean

        #Remove the "excess mean" mu_j from each random effect \alpha_{j,g}
        #and add the summd mass back to the betas.
        vi_alpha_mean <- vi_alpha_mean - M_prime_one %*% vi_mu_j
        vi_beta_mean <- vi_beta_mean + t(M_mu_to_beta) %*% vi_mu_j
        
        variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
        vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
        
        accept.PX <- TRUE
      }else if (parameter_expansion == 'translation'){
        
        prior.ELBO <- calculate_ELBO(
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = vi_sigma_outer_alpha,
          vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
          vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
          log_det_joint_var = log_det_joint_var
        )
        
        if (factorization_method == 'weak'){
          stop('no Translation PX for weak yet...')
        }
        R_ridge <- vecR_ridge_general(
            L = vi_alpha_decomp,
            Z = mapping_new_Z,
            pg_mean = vi_pg_mean,   
            M = Mmap,
            mapping_J = mapping_J, start_z = start_base_Z,
            d = d_j)
        
        
        R_design <- vecR_design(alpha_mu = as.vector(vi_alpha_mean), Z = mapping_new_Z, M = Mmap, mapping_J = mapping_J, d = d_j,
          start_z = start_base_Z)
        
        XR <- cbind(X, R_design)
        R_ridge <- bdiag(zeromat_beta, R_ridge)
        
        update_expansion_XR <- solve(t(XR) %*% diag_vi_pg_mean %*% XR + R_ridge, t(XR) %*% s)
        
        update_expansion_bX <- update_expansion_XR[1:p.X,,drop=F]
        
        update_expansion_R <- mapply(split(update_expansion_XR[-1:-p.X], rep(1:number_of_RE, d_j^2)), d_j, SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d)})

        prop_vi_sigma_alpha <- mapply(vi_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE, 
                                      FUN=function(Phi, R){R %*% Phi %*% t(R)})
        
        if (prevent_degeneracy){
          update_expansion_R <- mapply(prop_vi_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE, FUN=function(prop, r){
            #If determinant is negative or very small, do not allow!
            det_prop <- det(prop)
            if (det_prop < 1e-10 | det(r) < 0){
              r <- diag(ncol(r))
            }
            return(r)
          })
          
          prop_vi_sigma_alpha <- mapply(vi_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE, 
                                        FUN=function(Phi, R){R %*% Phi %*% t(R)})
        }
        
        
        mapping_for_R_block <- make_mapping_alpha(update_expansion_R, px.R = TRUE)
        update_expansion_Rblock <- prepare_T(mapping = mapping_for_R_block, levels_per_RE = g_j, num_REs = number_of_RE,
                          variables_per_RE = d_j, running_per_RE = breaks_for_RE, cyclical = FALSE, px.R = TRUE)
        
        #all.equal(update_expansion_Rblock, bdiag(mapply(update_expansion_R, g_j, FUN=function(i,g){bdiag(lapply(1:g, FUN=function(k){i}))})))
        
        update_expansion_mu <- t(M_prime) %*% vi_alpha_mean
        
        update_expansion_R_logdet <- log(sapply(update_expansion_R, det))

        
        prop_vi_beta_mean <- t(M_mu_to_beta) %*% bdiag(update_expansion_R) %*% update_expansion_mu +
          update_expansion_bX 
        prop_vi_alpha_mean <- update_expansion_Rblock %*% (vi_alpha_mean - M_prime_one %*% update_expansion_mu)
        
        
        #L^T L = Variance
        #R Var R^T --->
        # L %*% R^T
        prop_vi_alpha_decomp <- vi_alpha_decomp %*% t(update_expansion_Rblock)
        prop_log_det_alpha_var <- log_det_alpha_var + 2 * sum(update_expansion_R_logdet * g_j)
        
        prop_variance_by_alpha_jg <- calculate_expected_outer_alpha(L = prop_vi_alpha_decomp, alpha_mu = as.vector(prop_vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
        prop_vi_sigma_outer_alpha <- prop_variance_by_alpha_jg$outer_alpha
        
        prop.ELBO <- calculate_ELBO(
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = prop_vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = prop_vi_sigma_outer_alpha, 
          vi_beta_mean = prop_vi_beta_mean, vi_alpha_mean = prop_vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = prop_log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = prop_vi_alpha_decomp, 
          choose_term = choose_term
        )
        
        if (prop.ELBO$ELBO > prior.ELBO$ELBO){
          #Accept the PX-VB adjustment (OSL).
          vi_beta_mean <- prop_vi_beta_mean
          vi_alpha_mean <- prop_vi_alpha_mean
          vi_sigma_alpha <- prop_vi_sigma_alpha
          vi_alpha_decomp <- prop_vi_alpha_decomp
          log_det_alpha_var <- prop_log_det_alpha_var
          variance_by_alpha_jg <- prop_variance_by_alpha_jg
          vi_sigma_outer_alpha <- prop_vi_sigma_outer_alpha
          accept.PX <- TRUE
        }else{
          accept.PX <- FALSE
        }
        accepted_times <- accept.PX + accepted_times
        
        rm(prop_vi_beta_mean, prop_vi_alpha_mean, prop_vi_sigma_alpha, prop_vi_alpha_decomp,
           prop_log_det_alpha_var, prop_variance_by_alpha_jg, prop_vi_sigma_outer_alpha)
        
        
        rownames(vi_alpha_mean) <- fmt_names_Z
        
      }else if (parameter_expansion == 'none'){
        accept.PX <- TRUE
      }else{stop('Invalid PX method.')}
      
      #Adjust the terms in the ELBO calculation that are different.

      if (accept.PX){
        final.ELBO <- calculate_ELBO(
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = vi_sigma_outer_alpha,
          vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
          vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
          log_det_joint_var = log_det_joint_var
        )
      }else{
        if (accept.PX){
          final.ELBO <- prop.ELBO
        }else{
          final.ELBO <- prior.ELBO
        }
      }
      
      if (debug_ELBO & it != 1){
        debug_ELBO.1$step <- 1
        debug_ELBO.2$step <- 2
        debug_ELBO.3$step <- 3
        final.ELBO$step <- 4
        update_ELBO <- bind_rows(debug_ELBO.1, debug_ELBO.2, debug_ELBO.3, final.ELBO)
        update_ELBO$it <- it
        store_ELBO <- bind_rows(store_ELBO, update_ELBO)
      }else{
        final.ELBO$it <- it
        store_ELBO <- bind_rows(store_ELBO, final.ELBO)
      }
      
      ##Change diagnostics
      change_elbo <- final.ELBO - lagged_ELBO
      
      change_alpha_mean <- max(abs(vi_alpha_mean - lagged_alpha_mean))
      change_beta_mean <- max(abs(vi_beta_mean - lagged_beta_mean))
      change_sigma_mean <- mapply(vi_sigma_alpha, lagged_sigma_alpha, FUN=function(i,j){max(abs(i-j))})
      
      if (factorization_method == 'weak'){
        change_joint_var <- 0 #change_joint_var <- max(abs(vi_joint_decomp - lagged_joint_decomp))
        change_alpha_var <- change_beta_var <- 0
      }else{
        change_joint_var <- 0
        change_alpha_var <- max(abs(vi_alpha_decomp - lagged_alpha_decomp))
        change_beta_var <- max(abs(vi_beta_decomp - lagged_beta_decomp))
      }
      
      change_vi_r <- vi_r - lagged_vi_r
      
      if (debug_param){
        store_beta[it,] <- as.vector(vi_beta_mean)
      }
      change_all <- data.frame(change_alpha_mean, change_beta_mean, t(change_sigma_mean), change_alpha_var, change_beta_var, change_joint_var, change_vi_r)
      if ( (max(change_all) < tolerance_parameters) | (change_elbo$ELBO > 0 & change_elbo$ELBO < tolerance_elbo) ){
        if (!quiet){
          message(paste0('Converged after ', it, ' iterations with ELBO change of ', round(change_elbo[1], 1 + abs(floor(log(tolerance_elbo)/log(10))))))
          message(paste0('The largest change in any variational parameter was ', round(max(change_all), 1 + abs(floor(log(tolerance_parameters)/log(10))))))
        }
        break
      }
      if (!quiet & (it %% print_prog == 0)){
        plot(vi_alpha_mean)
        message(paste0('ELBO Change: ', round(change_elbo$ELBO, 10)))
        message('Other Parameter Changes')
        print(data.frame(change_all, accept.PX))
      }
      
      lagged_alpha_mean <- vi_alpha_mean
      lagged_beta_mean <- vi_beta_mean
      lagged_alpha_decomp <- vi_alpha_decomp
      lagged_beta_decomp <- vi_beta_decomp
      lagged_sigma_alpha <- vi_sigma_alpha
      lagged_vi_r <- vi_r
      lagged_ELBO <- final.ELBO
    }
    if (it == iterations){
      message(paste0('Ended without Convergence after', it, ' iterations : ELBO change of ', round(change_elbo[1], abs(floor(log(tolerance_elbo)/log(10))))))
    }
    
    if (debug_ELBO){
      d.ELBO <- with(store_ELBO, ELBO - dplyr::lag(ELBO))
      sum.ELBO <- store_ELBO
      sum.ELBO$diff <- d.ELBO
      sum.ELBO <- sum.ELBO %>% group_by_at(.vars = 'step') %>% summarize(negative = mean(.data$diff < 0, na.rm=T))  
    }else{
      d.ELBO <- NULL
    }
    if (parameter_expansion == 'translation'){
      final.ELBO$accepted_PX <- accepted_times / it
    }
    output <- list(beta = list(mean = vi_beta_mean),
                   ELBO = final.ELBO, debug_ELBO = d.ELBO,
                   ELBO_trajectory = store_ELBO,
                   parameter.change = change_all,
                   sigma = list(cov = vi_sigma_alpha, df = vi_sigma_alpha_nu),
                   alpha = list(mean = vi_alpha_mean)
    )
    
    if (debug_param){
      output$parameter_trajectory <- list(beta = store_beta)
    }
    if (factorization_method == 'weak'){
      output$joint <- vi_joint_decomp
    }
    if (return_data){
      output$data <- list(X = X, Z = Z, y = y, trials = trials)
    }
    output$formula <- formula
    output$factorization_method <- factorization_method
    output$alpha$dia.var <- unlist(lapply(variance_by_alpha_jg$variance_jg, FUN=function(i){as.vector(sapply(i, diag))}))
    output$beta$var <- t(vi_beta_decomp) %*% vi_beta_decomp
    output$beta$decomp_var <- vi_beta_decomp
    
    if (outcome == 'negbin'){
      output$r <- list(mean = vi_r, variance = NA)
    }

    output$internal.parameters <- list(names_of_RE = names_of_RE, d_j = d_j, g_j = g_j)
    
    output$alpha$var <- variance_by_alpha_jg$variance_jg
    output$alpha$decomp_var <- vi_alpha_decomp
    
    output$MAVB_parameters <- list(
      M_mu_to_beta = M_mu_to_beta,
      M_prime = M_prime,
      M_prime_one = M_prime_one,
      outer_alpha_RE_positions = outer_alpha_RE_positions,
      d_j = d_j, g_j = g_j
    )
    class(output) <- 'vglmer'
    return(output)
}

#' Predict after vglmer
#' 
#' Get linear predictor for new observations after using vglmer.
#' 
#' @param model Object from vglmer.
#' @param newdata Data to get predictions on.
#' @param samples How many samples to draw? 0, default, gets the expected value.
#'   Two methods: 
#'   \itemize{
#'   \item Number - draw that many samples.
#'   \item Matrix - use previous samples for prediction. Currently used for MAVB
#'   but will be depreciated in updated version.
#'   }
#' @param samples_only Return only samples *not* linear predictor.
#' @param summary Return summary of linear predictor, not full posterior.
#' @param allow_missing_levels Allow prediction for random effects not in model.
#'   As is standard, give an estimate of "0" for that effect.
#' @export
vglmer_predict <- function(model, newdata, 
                           samples = 0, samples_only = FALSE, 
                           summary = TRUE, allow_missing_levels = FALSE){
  
  fmla <- model$formula
  #Extract X (FE design matrix)
  X <- model.matrix(nobars(fmla), data = newdata)
  
  orig_X_names <- rownames(model$beta$mean)
  if (!identical(colnames(X), orig_X_names)){
    print(all.equal(colnames(X), orig_X_names))
    stop('Misaligned Fixed Effects')
  }
  
  #Extract the Z (Random Effect) design matrix.
  mk_Z <- mkReTrms(findbars(fmla), model.frame(subbars(fmla), data = newdata), reorder.terms = FALSE, reorder.vars = FALSE)
  Z <- t(mk_Z$Zt)
  
  #RE names and names of variables included for each.
  names_of_RE <- mk_Z$cnms
  
  number_of_RE <- length(mk_Z$Gp) - 1
  #The position that demarcates each random effect.
  #That is, breaks_for_RE[2] means at that position + 1 does RE2 start.
  breaks_for_RE <- c(0, cumsum(diff(mk_Z$Gp)))
  #Dimensionality of \alpha_{j,g}, i.e. 1 if random intercept
  #2 if random intercept + random slope
  d_j <- lengths(names_of_RE)
  #Number of GROUPs for each random effect.
  g_j <- diff(mk_Z$Gp)/d_j
  
  #Empty vector to build the formatted names for each random effect.
  fmt_names_Z <- c()
  init_Z_names <- colnames(Z)
  for (v in 1:number_of_RE){
    name_of_effects_v <- names_of_RE[[v]]
    
    mod_name <- rep(name_of_effects_v, g_j[v])
    
    levels_of_re <- init_Z_names[(1+breaks_for_RE[v]):breaks_for_RE[v+1]]
    
    fmt_names_Z <- c(fmt_names_Z, paste0(names(names_of_RE)[v], ' @ ', mod_name, ' @ ', levels_of_re))
    
  }
  colnames(Z) <- fmt_names_Z
  #####
  ###Confirm Alignment of the Z
  #####
  orig_Z_names <- rownames(model$alpha$mean)
  
  not_in_original_Z <- setdiff(fmt_names_Z, orig_Z_names)
  not_in_new_Z <- setdiff(orig_Z_names, fmt_names_Z)
  
  if (length(not_in_original_Z) > 0){
    if (!allow_missing_levels){stop('New levels not allowed unless allow_missing_levels = TRUE')}
  }
  
  in_both <- intersect(fmt_names_Z, orig_Z_names)
  recons_Z <- drop0(sparseMatrix(i = 1, j = 1, x = 0, dims = c(nrow(Z), length(orig_Z_names))))
  colnames(recons_Z) <- orig_Z_names
  rownames(recons_Z) <- 1:nrow(recons_Z)
  
  recons_Z[,match(in_both, orig_Z_names)] <- Z[,match(in_both, fmt_names_Z)] 
  
  #Check that the entirely missing columns match those not in the original
  checksum_align <- setdiff(not_in_new_Z, sort(names(which(colSums(recons_Z != 0) == 0))))
  if (length(checksum_align) > 0){
    stop('Alignment Error')
  }
  
  Z <- recons_Z
  rm(recons_Z)
  ####
  
  XZ <- cbind(X, Z)
  factorization_method <- model$factorization_method
  if (is.matrix(samples)){
    if (ncol(samples) != ncol(XZ)){
      stop('Samples must be {m, ncol(Z) + ncol(X)}')
    }
    samples <- t(samples)
    only.lp <- FALSE
    print('Using Provided Samples')
  }else{
    if (samples == 0){
      only.lp <- TRUE
    }else{
      only.lp <- FALSE
    }
    if (factorization_method %in% c('strong', 'partial')){
      
      vi_alpha_mean <- model$alpha$mean
      vi_alpha_decomp <- model$alpha$decomp_var
      
      p.Z <- nrow(vi_alpha_mean)
      
      vi_beta_mean <- model$beta$mean
      vi_beta_decomp <- model$beta$decomp_var
      
      p.X <- nrow(vi_beta_mean)
      
      if (!only.lp){
        sim_init_alpha <- matrix(rnorm(samples * p.Z), ncol = samples)
        sim_init_alpha <- t(vi_alpha_decomp) %*% sim_init_alpha
        sim_init_alpha <- sim_init_alpha + kronecker(vi_alpha_mean, t(matrix(1, samples)))
        
        sim_init_beta <- matrix(rnorm(samples * p.X), ncol = samples)
        sim_init_beta <- t(vi_beta_decomp) %*% sim_init_beta
        sim_init_beta <- sim_init_beta + kronecker(vi_beta_mean, t(matrix(1,samples)))
      }else{
        sim_init_alpha<- vi_alpha_mean
        sim_init_beta <- vi_beta_mean
      }
      
    }else if (factorization_method == 'weak'){
      
      vi_alpha_mean <- model$alpha$mean
      p.Z <- nrow(vi_alpha_mean)
      
      vi_beta_mean <- model$beta$mean
      p.X <- nrow(vi_beta_mean)
      
      if (!only.lp){
        vi_joint_decomp <- model$joint
        sim_init_joint <- matrix(rnorm(samples * (p.X + p.Z)), ncol = samples)
        sim_init_joint <- t(vi_joint_decomp) %*% sim_init_joint
        
        sim_init_beta <- sim_init_joint[1:p.X, , drop = F]
        sim_init_alpha <- sim_init_joint[-1:-p.X, ,drop = F]
        
        rm(sim_init_joint)
        
        sim_init_alpha <- sim_init_alpha + kronecker(vi_alpha_mean, t(matrix(1, samples)))
        sim_init_beta <- sim_init_beta + kronecker(vi_beta_mean, t(matrix(1,samples)))
      }else{
        sim_init_alpha <- vi_alpha_mean
        sim_init_beta <- vi_beta_mean
      }
      
    }else{stop('')}
    
    samples <- rbind(sim_init_beta, sim_init_alpha)
    rm(sim_init_beta, sim_init_alpha)
  }
  
  if (samples_only){
    return(samples)
  }
  
  lp <- XZ %*% samples
  if (summary){
    if (!only.lp){
      lp <- t(apply(lp, MARGIN = 1, FUN=function(i){c(mean(i), var(i))}))
      lp <- data.frame(mean = lp[,1], var = lp[,2])  
    }else{
      lp <- as.vector(t(apply(lp, MARGIN = 1, FUN=function(i){mean(i)})))
    }
    return(lp)
  }else{
    return(lp)
  }
}



#' Perform MAVB after fitting vglmer
#' 
#' Given a model from vglmer, perform MAVB to improve approximation quality. See
#' dissertation for details. This function uses a naive loop so is slower than
#' necessary. This will be updated shortly.
#' 
#' @param model Model fit using vglmer
#' @param samples Samples to draw from MAVB distribution.
#' @param var_px Default (Inf); variance of working prior. Higher is more
#'   diffuse and thus likely better.
#' 
#' @importFrom purrr array_branch
#' @importFrom mvtnorm rmvnorm
#' @export
MAVB <- function(model, samples, var_px = Inf){
  
  # if (!inherits(model, 'vglmer')){
  #   stop('Must provide model from vglmer')
  # }
  
  M_prime <- model$MAVB_parameters$M_prime
  M_prime_one <- model$MAVB_parameters$M_prime_one
  M_mu_to_beta <- model$MAVB_parameters$M_mu_to_beta
  d_j <- model$MAVB_parameters$d_j
  g_j <- model$MAVB_parameters$g_j
  outer_alpha_RE_positions <- model$MAVB_parameters$outer_alpha_RE_positions
  factorization_method <- model$factorization_method
  
  if (model$factorization_method == 'weak'){
    decomp_joint <- model$joint
    joint_mean <- rbind(model$beta$mean, model$alpha$mean)
    p.XZ <- ncol(decomp_joint)
    p.X <- nrow(model$beta$mean)
  }else{
    decomp_varA <- model$alpha$decomp_var
    decomp_varB <- model$beta$decomp_var
    p.XZ <- ncol(decomp_varA) + ncol(decomp_varB)
    p.X <- nrow(model$beta$mean)
  }
  
  MAVB_sims <- matrix(NA, nrow = samples, ncol = p.XZ)
  regen_store <- MAVB_diff <- matrix(NA, nrow = samples, ncol = sum(d_j))
  
  n_MAVB <- samples

  for (it in 1:n_MAVB){
    if (it %% 1000 == 0){cat('.')}
    #Sim from VARIATIONAL approximation to posterior.
    if (factorization_method == 'weak'){
      sim_joint <- joint_mean + t(decomp_joint) %*% rnorm(p.XZ)
      sim_beta <- sim_joint[1:p.X,,drop=F]
      sim_a <- sim_joint[-1:-p.X,,drop=F]
    }else{
      sim_a <- model$alpha$mean + t(decomp_varA) %*% rnorm(ncol(decomp_varA))
      sim_beta <- model$beta$mean + t(decomp_varB) %*% rnorm(ncol(decomp_varB))
    }
    sim_sigma <- mapply(model$sigma$df, model$sigma$cov, SIMPLIFY = FALSE, FUN=function(i,j){rInvWishart(n = 1, df = i, Sigma = j)[,,1]})
    #Working Prior
    if (var_px == Inf){
      sim_px <- rep(0, sum(d_j))
    }else{
      sim_px <- rnorm(sum(d_j), 0, sd = sqrt(var_px))
    }
    #Transform t^{-1}_a(z) = w
    sim_atilde <- sim_a + M_prime_one %*% sim_px
    sim_btilde <- sim_beta - t(M_mu_to_beta) %*% sim_px
    #Draw sim_px AGAIN from its full conditional
    if (var_px == Inf){
      #Use the LIMITING transition.
      var_redux <- as.matrix(bdiag(mapply(sim_sigma, g_j, SIMPLIFY = FALSE, FUN=function(S, g){S/g})))
      #Get the MEAN not the SUM
      mean_redux <- t(M_prime) %*% sim_atilde
    }else{
      var_redux <- solve(diag(x = 1/var_px, ncol = sum(d_j), nrow = sum(d_j)) + as.matrix(bdiag(mapply(sim_sigma, g_j, SIMPLIFY = FALSE, FUN=function(S, g){solve(S) * g}))))
      #Use the SUM
      mean_redux <- var_redux %*% solve(bdiag(sim_sigma)) %*% t(M_prime_one) %*% sim_atilde
    }
    regen_px <- rmvnorm(1, mean_redux, var_redux)
    regen_px <- t(regen_px)
    regen_store[it,] <- regen_px  
    MAVB_diff[it,] <- regen_px - sim_px
    final_atilde <- sim_atilde - M_prime_one %*% regen_px
    final_btilde <- sim_btilde + t(M_mu_to_beta) %*% regen_px
    
    MAVB_sims[it,] <- c(as.vector(final_btilde), as.vector(final_atilde))
  }
  return(MAVB_sims)
}

#' Linear Predictor Following MAVB
#' 
#' Combine prediction and MAVB into a single function.
#' 
#' @inheritParams MAVB
#' @inheritParams vglmer_predict
#' @export
MAVB_linpred <- function(model, newdata, samples, var_px = Inf, summary = TRUE){
  pxSamples <- MAVB(model = model, samples = samples, var_px = var_px)
  lp <- vglmer_predict(model, newdata = newdata, samples = pxSamples, summary = summary)
  return(lp)
}


#' @import lme4
get_RE_groups <- function(formula, data){
  stop('Figure out workaround for lme4:::')
  # bars <- findbars(formula)
  # names(bars) <- lme4:::barnames(bars)
  # fr <- model.frame(subbars(formula), data = data)
  # blist <- lapply(bars, simple_blist, fr, drop.unused.levels = F, reorder.vars = FALSE)
  # blist <- lapply(blist, FUN=function(i){i[c('ff', 'mm')]})
  # 
  # ff <- lapply(blist, FUN=function(i){i$ff})
  # ff <- lapply(ff, FUN=function(i){match(i, levels(i))})
  # mm <- lapply(blist, FUN=function(i){i$mm})
  # return(list(factor = ff, design = mm))
}

#' @import lme4
simple_blist <- function(x, frloc, drop.unused.levels = TRUE, reorder.vars = FALSE){
  stop('figure out workaround for lme4:::')
  # frloc <- factorize(x, frloc)
  # if (is.null(ff <- tryCatch(eval(substitute(lme4:::makeFac(fac), 
  #                                            list(fac = x[[3]])), frloc), error = function(e) NULL))) 
  #   stop("couldn't evaluate grouping factor ", deparse(x[[3]]), 
  #        " within model frame:", " try adding grouping factor to data ", 
  #        "frame explicitly if possible", call. = FALSE)
  # if (all(is.na(ff))) 
  #   stop("Invalid grouping factor specification, ", deparse(x[[3]]), 
  #        call. = FALSE)
  # if (drop.unused.levels) 
  #   ff <- factor(ff, exclude = NA)
  # nl <- length(levels(ff))
  # mm <- model.matrix(eval(substitute(~foo, list(foo = x[[2]]))), 
  #                    frloc)
  # if (reorder.vars) {
  #   mm <- mm[colSort(colnames(mm)), ]
  # }
  # list(ff = ff, nl = nl, mm = mm, cnms = colnames(mm))
}


#' Variance of Rows or Columns of Matrices
#' 
#' Base R implementation for variance. Analogue of rowMeans.
#' @name var_mat
#' @param matrix Matrix of numeric inputs.
rowVar <- function(matrix){apply(matrix, MARGIN = 1, var)}

#' @importFrom stats var
#' @rdname var_mat
colVar <- function(matrix){apply(matrix, MARGIN = 2, var)}


#' @param data Data to get predictions on.
#' @rdname hmc_samples
custom_HMC_linpred <- function(HMC, data){
  if (!requireNamespace('rstanarm', quietly = TRUE)){
    stop('rstanarm must be installed to analyze HMC objects.')
  }
  hmc_samples <- as.matrix(HMC)
  hmc_samples <- hmc_samples[, !grepl(colnames(hmc_samples), pattern='^Sigma')]
  
  parse_stan_names <- str_split(colnames(hmc_samples), pattern='^b\\[| |\\]')
  
  fmt_stan_names <- sapply(parse_stan_names, FUN=function(i){
    if (length(i) == 1){
      return(i)
    }else{
      i_one <- unlist(str_split(i[3], pattern=':'))
      return(paste(i_one[1], i[2], i_one[2], sep=' @ '))
    }
  })
  colnames(hmc_samples) <- fmt_stan_names
  
  hmc.XZ <- rstanarm::posterior_linpred(HMC, data = data, XZ = TRUE)
  
  hmc.linpred <- hmc.XZ %*% t(hmc_samples)
  
  return(hmc.linpred)
}

#' Get HMC samples but ordered to match vector of names provided
#' @export
#' @name hmc_samples
#' @param HMC Object from rstanarm
#' @param ordering vector of names to order the posterior samples.
#' @importFrom mvtnorm rmvnorm
custom_HMC_samples <- function(HMC, ordering){
  if (!requireNamespace('rstanarm', quietly = TRUE)){
    stop('rstanarm must be installed to analyze HMC objects.')
  }
  
  hmc_samples <- as.matrix(HMC)
  hmc_samples <- hmc_samples[, !grepl(colnames(hmc_samples), pattern='^Sigma')]
  
  parse_stan_names <- str_split(colnames(hmc_samples), pattern='^b\\[| |\\]')
  
  fmt_stan_names <- sapply(parse_stan_names, FUN=function(i){
    if (length(i) == 1){
      return(i)
    }else{
      i_one <- unlist(str_split(i[3], pattern=':'))
      return(paste(i_one[1], i[2], i_one[2], sep=' @ '))
    }
  })
  colnames(hmc_samples) <- fmt_stan_names
  hmc_samples <- hmc_samples[, match(ordering, colnames(hmc_samples))]
  return(return(hmc_samples))
}

#' Get samples from GLMER
#' 
#' Order samples from glmer to match names from vglmer.
#' 
#' @param glmer object fitted using glmer
#' @param samples number of samples to draw
#' @param ordering order of output
#' 
#' @export
#' @importFrom stats rnorm
custom_glmer_samples <- function(glmer, samples, ordering){
  
  fmt_glmer <- format_glmer(glmer)
  
  glmer_samples <- mapply(fmt_glmer$mean, fmt_glmer$var, FUN=function(m,v){rnorm(samples, mean = m, sd = sqrt(v))})
  colnames(glmer_samples) <- fmt_glmer$name
  
  glmer_samples <- glmer_samples[, match(ordering, colnames(glmer_samples))]
  return(glmer_samples)  
}


update_r <- function(vi_r, y, X, Z, factorization_method,
    vi_beta_mean, vi_alpha_mean, 
    vi_joint_decomp, vi_beta_decomp, vi_alpha_decomp){
  
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
  
  # opt.r <- optimize(f = function(lr,y,psi,zVz){
  #   Q.adjust.r(r = exp(lr), y = y, psi = psi, zVz = zVz)},
  #   lower = log(sqrt(.Machine$double.eps)), upper = 10,
  #   y = y, psi = ex_XBZA, zVz = var_XBZA,
  #   maximum = T)
  
  opt_ln_r <- optim(par = log(vi_r), fn = Q.adjust.r,
                    y = y, psi = ex_XBZA, zVz = var_XBZA,
                    control = list(fnscale = - 1), method = 'L-BFGS-B')
  if (opt_ln_r$value < Q.adjust.r(log(vi_r), y = y, psi = ex_XBZA, zVz = var_XBZA)){
    warning('Optim for r decreased objective.')
  }else{
    vi_r <- exp(opt_ln_r$par)
  }
  # vi_r <- exp(opt.r$maximum)
  
  return(vi_r)
}

# Q.adjust.r <- function(r, y, psi, zVz){
#   t1 <- -(y+r) * log(2) + (y-r)/2 * (psi  - log(r)) - (r + y) * log(cosh(1/2 * sqrt(zVz + (psi - log(r))^2)))
#   t2 <- lgamma(y + r) - lgamma(r) - lgamma(y+1)
#   return(sum(t1 + t2))
# }

Q.adjust.r <- function(ln_r, y, psi, zVz){
  t1 <- -(y+exp(ln_r)) * log(2) + (y-exp(ln_r))/2 * (psi  - ln_r) +
    -(exp(ln_r) + y) * log(cosh(1/2 * sqrt(zVz + (psi - ln_r)^2)))
  t2 <- lgamma(y + exp(ln_r)) - lgamma(exp(ln_r)) - lgamma(y+1)
  return(sum(t1 + t2))
}









