#' Variational Inference for Non-Linear Hierarchical Models
#' 
#' Estimate a hierarchical model using mean-field variational
#' inference. Accepts standard syntax to glmer: y ~ X + (1 + Z | g). Options are
#' described below.
#' 
#' @param formula Standard glmer-style formula for random effects.
#' @param data data.frame containing the outcome and variables.
#' @param family Options are "binomial" or "negbin". If "binomial", outcome must
#'   be either {0,1} (binary) or cbind(success, failure) as per standard glm(er)
#'   syntax. Non-integer values are permitted for binomial if "force_whole" is
#'   set to FALSE in vglmer_control.
#' @param control Control additional arguments. Must be made using
#'   vglmer_control(); see for documentation for additional details.
#'
#' @examples
#' 
#' sim_data <- data.frame(x = rnorm(100), 
#' y = rbinom(100, 1, 0.5), 
#' g =sample(letters, 100, replace = TRUE))
#' 
#' # Run with defaults
#' est_vglmer <- vglmer(y ~ x + (x|g), data = sim_data, family = "binomial")
#' 
#' # Simple prediction 
#' predict(est_vglmer, newdata = sim_data)
#' 
#' # Summarize results
#' summary(est_vglmer)
#' 
#' # Extract parameters
#' coef(est_vglmer); vcov(est_vglmer)
#' 
#' # Comparability with lme4, 
#' # although ranef is formatted differently.
#' ranef(est_vglmer); fixef(est_vglmer)
#' 
#' #' # Run with stronger (i.e. less good) approximation
#' \dontrun{
#' vglmer(y ~ x + (x|g), data = sim_data, 
#' control = vglmer_control(factorization_method = "strong"),
#' family = "binomial")
#' }
#'
#' @importFrom dplyr select group_by group_by_at summarize n lag
#' @importFrom lme4 mkReTrms findbars subbars
#' @importFrom stats model.response model.matrix model.frame rnorm rWishart
#'   qlogis optim
#' @importFrom rlang .data
#' @importFrom graphics plot
#' @importFrom checkmate assert assert_formula assert_choice
#'   check_data_frame
#' @useDynLib vglmer
#' @export
vglmer <- function(formula, data, family, control = vglmer_control()){
  
  #Verify integrity of parameter arguments
  checkdf <- check_data_frame(data, null.ok = TRUE)
  if (checkdf != TRUE){
    warning(paste0('data is not a data.frame? Behavior may be unexpected: ', checkdf))
  }
  assert_formula(formula)
  # Delete the missing data 
  # (i.e. sub out the random effects, do model.frame)
  #
  nobs_init <- nrow(data)
  data <- model.frame(subbars(formula), data)
  nobs_complete <- nrow(data)
  missing_obs <- nobs_init - nobs_complete
  if (length(missing_obs) == 0){
    missing_obs <- '??'
  }
  assert_choice(family, c('negbin', 'binomial'))
  if (!inherits(control, 'vglmer_control')){
    stop('control must be object from vglmer_control().')
  }
  
  do_timing <- control$do_timing
  factorization_method <- control$factorization_method
  print_prog <- control$print_prog
  iterations <- control$iterations
  quiet <- control$quiet
  parameter_expansion <- control$parameter_expansion
  tolerance_elbo <- control$tolerance_elbo
  tolerance_parameters <- control$tolerance_parameters
  prevent_degeneracy <- control$prevent_degeneracy
  debug_param <- control$debug_param
  linpred_method <- control$linpred_method
  vi_r_method <- control$vi_r_method
  vi_r_val <- control$vi_r_val
  debug_ELBO <- control$debug_ELBO
  verbose_time <- control$verbose_time
  
  if (do_timing){
    if (!requireNamespace('tictoc', quietly = TRUE)){
      stop('tictoc must be installe to do timing')
    }
    tic <- tictoc::tic
    toc <- tictoc::toc
    tic.clear <- tictoc::tic.clear
    tic.clearlog <- tictoc::tic.clearlog
    
    tic.clear(); tic.clearlog();
    tic('Prepare Model')
  }
  if (!(factorization_method %in% c('weak', 'strong', 'partial'))){
    stop('factorization_method must be in weak, strong, or partial.')
  }
  if (is.null(print_prog)){
    print_prog <- max(c(1, floor(iterations / 20)))
  }
  if (!(family %in% c('binomial', 'negbin'))){
    stop('family must be "binomial" or "negbin".')
  }
  #Extract outcome y
  y <- model.response(model.frame(nobars(formula), data = data))
  names(y) <- NULL
  
  if (family == 'binomial'){
    if (is.matrix(y)){
    # if (!(class(y) %in% c('numeric', 'integer'))){
      if (min(y) < 0){
        stop('Negative Numbers not Permitted in outcome')
      }
      is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
      if (any(is.wholenumber(y) == FALSE)){
        if (control$force_whole){
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
      if (!all(y %in% c(0,1)) & family == 'binomial'){
        stop('Only {0,1} outcomes permitted for numeric y.')
      }
      trials <- rep(1, length(y))
    }
  }else{
    
    if (!(class(y) %in% c('numeric', 'integer'))){
      stop("Must provide vector of numbers with negbin.")
    }
    
    if (min(y) < 0){
      stop('Negative numbers not Permitted in outcome')
    }
    
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if (any(is.wholenumber(y) == FALSE)){
      if (control$force_whole){
        stop('If force_whole = TRUE, must provide whole numbers')
      }else{
        warning('Non-integer numbers in y')
      }
    }
    
  }
  
  if (family == 'binomial'){
    ELBO_type <- 'augmented'
  }else if (family == 'negbin'){
    ELBO_type <- 'profiled'
  }else{stop('Check ELBO_type')}
  
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

    if (family == 'binomial'){
      s <- y - trials/2
      vi_pg_b <- trials
      vi_r_mu <- 0
      vi_r_sigma <- 0
      choose_term <- sum(lchoose(n = round(trials), k = round(y)))
    }else{
      #Initialize
      if (vi_r_method == 'fixed'){
        vi_r_mu <- vi_r_val
        vi_r_mean <- exp(vi_r_mu)
        vi_r_sigma <- 0
      }else if (vi_r_method == 'VEM'){
        if (!requireNamespace('MASS', quietly = TRUE)){
          stop('Install MASS for negbin')          
        }
        vi_r_mean <- MASS::glm.nb(y ~ 1)$theta
        vi_r_mu <- log(vi_r_mean)
        vi_r_sigma <- 0
      }else{
        init_r <- optim(par = 0, fn = VEM.PELBO.r, method = 'L-BFGS', hessian = T,
                        control = list(fnscale = -1), y = y, psi = rep(log(mean(y)), length(y)), zVz = 0)
        vi_r_mu <- init_r$par
        vi_r_sigma <- as.numeric(-1/init_r$hessian)
        
        vi_r_mean <- exp(vi_r_mu + vi_r_sigma/2)
      }
      s <- (y - vi_r_mean)/2
      vi_pg_b <- y + vi_r_mean
      
      choose_term <- NA #-sum(lgamma(y + 1)) - sum(y) * log(2)
    }

    #Initalize variational parameters.
    #Note that we keep a sparse matrix or lowertri such that
    #t(vi_beta_decomp) %*% vi_beta_decomp = VARIANCE

    vi_beta_decomp <- Diagonal(x = rep(0, ncol(X)))
    vi_alpha_decomp <- Diagonal(x = rep(0, ncol(Z)))
    
    vi_sigma_alpha_nu <-  g_j
    
    prior_variance <- control$prior_variance
    
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
        
    if (control$init == 'EM'){
      if (family == 'negbin'){
        EM_init <- EM_prelim_nb(X = X, Z = Z, y = y, est_r = exp(vi_r_mu), iter = 15, ridge = 4)
      }else{
        EM_init <- EM_prelim_logit(X = X, Z = Z, s = s, pg_b = vi_pg_b, iter = 15, ridge = 4)
      }
      
      vi_beta_mean <- matrix(EM_init$beta)
      vi_alpha_mean <- matrix(EM_init$alpha)
      
      vi_sigma_alpha <- calculate_expected_outer_alpha(alpha_mu = matrix(EM_init$alpha), 
        L = sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(Z), 2)), 
        re_position_list = outer_alpha_RE_positions)
      vi_sigma_alpha <- mapply(vi_sigma_alpha$outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN=function(i,j){i+j})
    }else if (control$init == 'random'){
      set.seed(control$random_seed)
      vi_beta_mean <- rnorm(ncol(X))
      vi_alpha_mean <- rep(0, ncol(Z))
      
      vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN=function(d,g){
        rWishart(n = 1, df = g, Sigma = diag(d))[,,1]
      })
    }else if (control$init == 'zero'){
      
      vi_beta_mean <- rep(0, ncol(X))
      
      if (family == 'binomial'){
        vi_beta_mean[1] <- qlogis(sum(y)/sum(trials))
        #Deal with edge case of 0.5
        if (vi_beta_mean[1] == 0){
          vi_beta_mean[1] <- runif(1, -.1, .1)
        }
        
      }else if (family == 'negbin'){
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
      log_det_beta_var <- log_det_alpha_var <- NULL
    }else{
      vi_joint_decomp = NULL
      log_det_joint_var <- NULL
    }
    
    # if (!quiet){warning('Check vi sigma alpha nu')}
    
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
    lagged_vi_r_mu <- -Inf
    
    lagged_ELBO <- -Inf
    accepted_times <- NA
    
    if (parameter_expansion == 'translation'){
      stop('parameter_expansion_translation not allowed yet.')
      # accepted_times <- 0
      # zeromat_beta <- drop0(Diagonal(x = rep(0, ncol(X))))
      # 
      # parsed_RE_groups <- get_RE_groups(formula = formula, data = data)
      # 
      # mapping_new_Z <- do.call('cbind', parsed_RE_groups$design)
      # mapping_J <- split(1:sum(d_j^2), rep(1:length(d_j), d_j^2))
      # mapping_J <- lapply(mapping_J, FUN=function(i){i-1})
      # mapping_J <- sapply(mapping_J, min)
      # 
      # mapping_to_re <- parsed_RE_groups$factor
      # mapping_to_re <- array_branch(do.call('cbind', mapping_to_re), margin = 1)
      # mapping_to_re <- lapply(mapping_to_re, FUN=function(i){
      #   mapply(outer_alpha_RE_positions, i, SIMPLIFY = FALSE, FUN=function(a,b){a[[b]]})
      # })
      # Mmap <- t(sapply(mapping_to_re, FUN=function(i){as.integer(sapply(i, min))}))
      # 
      # start_base_Z <- cumsum(c(0,d_j))[-(number_of_RE+1)]
      # names(start_base_Z) <- NULL
      # 
      # rm(parsed_RE_groups, mapping_to_re)
    }
    store_ELBO <- data.frame()
    
    if (debug_param){
      store_beta <- array(NA, dim = c(iterations, ncol(X)))
    }
    if (do_timing){
      toc(quiet = verbose_time, log = TRUE)
      tic.clear()
    }
    ##Begin VI algorithm:
    if (!quiet){message('Begin Regression')}
    for (it in 1:iterations){
      if (it %% print_prog == 0){cat('.')}
      ###
      ##Polya-Gamma Updates
      ###
      #Get the x_i^T Var(beta) x_i terms.
      if (do_timing){
        tic('Update PG')
      }
      
      if (factorization_method == 'weak'){
        # joint_quad <- rowSums( (joint.XZ %*% t(vi_joint_decomp))^2 )
        # vi_joint_decomp <<- vi_joint_decomp
        # joint.XZ <<- joint.XZ
        joint_quad <- cpp_zVz(Z = joint.XZ, V = as(vi_joint_decomp, 'dgCMatrix')) + vi_r_sigma
        vi_pg_c <- sqrt(as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean - vi_r_mu)^2 + joint_quad)
      }else{
        # vi_beta_decomp <<- vi_beta_decomp
        # vi_alpha_decomp <<- vi_alpha_decomp
        # X <<- X
        # Z <<- Z
        # beta_quad <- cpp_zVz(Z = drop0(X), V = make_dgC(vi_beta_decomp))
        # alpha_quad <- cpp_zVz(Z = Z, V = make_dgC(vi_alpha_decomp))
        beta_quad <- rowSums( (X %*% t(vi_beta_decomp))^2 )
        alpha_quad <- rowSums( (Z %*% t(vi_alpha_decomp))^2 )
        vi_pg_c <- sqrt(as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean - vi_r_mu)^2 + beta_quad + alpha_quad + vi_r_sigma)
      }
      vi_pg_mean <- vi_pg_b/(2 * vi_pg_c) * tanh(vi_pg_c / 2)
      diag_vi_pg_mean <- sparseMatrix(i = 1:N, j = 1:N, x = vi_pg_mean)
      if (debug_ELBO & it != 1){
        debug_ELBO.1 <- calculate_ELBO(ELBO_type = ELBO_type,
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s, y = y,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = vi_sigma_outer_alpha,
          vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
          vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
          log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma
        )
        
      }
      
      if (do_timing){
        toc(quiet = verbose_time, log = TRUE)
        tic('Prepare Sigma')
      }
      
      #Process Sigma_j for manipulation      
      #if Sigma_{j} is InverseWishart(a,Phi)
      #Then E[Sigma^{-1}_j] = a * Phi^{-1}
      if (factorization_method == 'strong'){
        cyclical_T <- TRUE
      }else{
        cyclical_T <- FALSE
      }
      inv_mapping_alpha <- mapply(vi_sigma_alpha_nu, lapply(vi_sigma_alpha, solve), 
          SIMPLIFY = FALSE, FUN=function(a,b){a * b})
      inv_mapping_alpha <- make_mapping_alpha(inv_mapping_alpha)
      
      Tinv <- prepare_T(mapping = inv_mapping_alpha, levels_per_RE = g_j, num_REs = number_of_RE,
        variables_per_RE = d_j, running_per_RE = breaks_for_RE, cyclical = cyclical_T)
      if (!cyclical_T){
        Tinv <- as(Tinv, 'dgCMatrix')
      }else{
        Tinv <- lapply(Tinv, FUN=function(i){as(i, 'dgCMatrix')})  
      }
      if (do_timing){
        toc(quiet = verbose_time, log = T)
        tic('Update Beta')
      }
      if (factorization_method == 'weak'){
        ##Update <beta, alpha> jointly
        chol.update.joint <- LinRegChol(X = joint.XZ, omega = diag_vi_pg_mean, 
                                        prior_precision =  bdiag(zero_mat, Tinv), 
                                        y = s + vi_pg_mean * vi_r_mu)
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
                                          y= s + vi_pg_mean * vi_r_mu)
          vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
          vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))
          
          vi_beta_decomp <- solve(t(chol(as.matrix(t(X) %*% diag_vi_pg_mean %*% X))))
          log_det_beta_var <- 2 * sum(log(diag(vi_beta_decomp)))
          
          chol.update.alpha <- LinRegChol(X = Z, omega = diag_vi_pg_mean, prior_precision =  Tinv, 
                                          y =  s + vi_pg_mean * vi_r_mu)
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
        # family_s <<- s
        # diag_vi_pg_mean <<- diag_vi_pg_mean
        # running_log_det_alpha_var <<- running_log_det_alpha_var
        # number_of_RE <<- number_of_RE
        # zero_mat <<- zero_mat
        # p.X <<- p.X
        
        
        if (linpred_method == 'joint'){
          
          joint.XZ <- cbind(X,Z)  
          chol.update.joint <- LinRegChol(X = joint.XZ, omega = diag_vi_pg_mean, prior_precision =  bdiag(zero_mat, bdiag(Tinv)),  y = s + vi_pg_mean * vi_r_mu)
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
                                 y = as.vector(s + vi_pg_mean * vi_r_mu - diag_vi_pg_mean %*% (X %*% vi_beta_mean + Z_negj %*% vi_alpha_mean[-index_j])))
            vi_alpha_mean[index_j] <- chol.j$mean 
            
            Pmatrix <- sparseMatrix(i = 1:ncol(Z_j), j = 1 + chol.j$Pindex, x = 1)
            
            running_log_det_alpha_var[j] <- - 2 * sum(log(diag(chol.j$origL)))
            vi_alpha_decomp[index_j, index_j] <- solve(chol.j$origL) %*% t(Pmatrix)
          }
          
          #vi_alpha_decomp <- bdiag(vi_alpha_decomp)
          log_det_alpha_var <- sum(running_log_det_alpha_var)
          
          chol.update.beta <- LinRegChol(X = as(X, 'sparseMatrix'), omega = diag_vi_pg_mean, prior_precision =  zero_mat, 
                                         y = as.vector(s + vi_pg_mean * vi_r_mu - diag_vi_pg_mean %*% Z %*% vi_alpha_mean))
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

        debug_ELBO.2 <- calculate_ELBO(ELBO_type = ELBO_type,
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s, y = y,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = vi_sigma_outer_alpha,
          vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
          vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
          log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, 
          vi_r_sigma = vi_r_sigma
        )
      }
      if (do_timing){
        toc(quiet = verbose_time, log = T)
        tic('Update Sigma')
      }
      ###
      # Update \Sigma_j
      ##
      
      variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
      vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
      vi_sigma_alpha <- mapply(vi_sigma_outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN=function(i,j){i+j})
  
      if (do_timing){
        toc(quiet = verbose_time, log = T)
        tic('Update Aux')
      }
      
      #Update the auxilary parameters
      if (family == 'negbin'){
        
        vi_r_param <- update_r(vi_r_mu = vi_r_mu, vi_r_sigma = vi_r_sigma,
         y = y, X = X, Z = Z, factorization_method = factorization_method,
         vi_beta_mean = vi_beta_mean, vi_beta_decomp = vi_beta_decomp,
         vi_alpha_mean = vi_alpha_mean, vi_alpha_decomp = vi_alpha_decomp,
         vi_joint_decomp = vi_joint_decomp, vi_r_method = vi_r_method,
         vi_pg_mean = vi_pg_mean)
        
        vi_r_mu <- vi_r_param[1]
        vi_r_sigma <- vi_r_param[2]
        vi_r_mean <- exp(vi_r_mu + vi_r_sigma/2)

        s <- (y - vi_r_mean)/2
        vi_pg_b <- y + vi_r_mean
        
      }
      
      if (do_timing){
        toc(quiet = verbose_time, log = T)
      }
      ###PARAMETER EXPANSIONS!
      if (debug_ELBO){
        debug_ELBO.3 <- calculate_ELBO(ELBO_type = ELBO_type,
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s, y = y,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = vi_sigma_outer_alpha,
          vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
          vi_joint_decomp = vi_joint_decomp,
          log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma, choose_term = choose_term
        )
      }
      
      if (parameter_expansion == 'mean'){
        if (do_timing){
          tic('Update PX')
        }
        #Do a simple mean adjusted expansion.
        #Get the mean of each random effect.
        vi_mu_j <- t(M_prime) %*% vi_alpha_mean

        #Remove the "excess mean" mu_j from each random effect \alpha_{j,g}
        #and add the summd mass back to the betas.
        vi_alpha_mean <- vi_alpha_mean - M_prime_one %*% vi_mu_j
        vi_beta_mean <- vi_beta_mean + t(M_mu_to_beta) %*% vi_mu_j
        
        variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, 
                                                               alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
        vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
        accept.PX <- TRUE
        if (do_timing){
          toc(quiet = verbose_time, log = TRUE)
        }
      }else if (parameter_expansion == 'translation'){
        stop('translation not yet implemented')
        # prior.ELBO <- calculate_ELBO(ELBO_type = ELBO_type,
        #   factorization_method = factorization_method,
        #   d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
        #   prior_sigma_alpha_nu = prior_sigma_alpha_nu,
        #   iw_prior_constant = iw_prior_constant,
        #   X = X, Z = Z, s = s, y = y,
        #   vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
        #   vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
        #   vi_sigma_outer_alpha = vi_sigma_outer_alpha,
        #   vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
        #   log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
        #   vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
        #   vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
        #   log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma
        # )
        # 
        # if (factorization_method == 'weak'){
        #   stop('no Translation PX for weak yet...')
        # }
        # R_ridge <- vecR_ridge_general(
        #     L = vi_alpha_decomp,
        #     Z = mapping_new_Z,
        #     pg_mean = vi_pg_mean,   
        #     M = Mmap,
        #     mapping_J = mapping_J, start_z = start_base_Z,
        #     d = d_j)
        # 
        # 
        # R_design <- vecR_design(alpha_mu = as.vector(vi_alpha_mean), Z = mapping_new_Z, M = Mmap, mapping_J = mapping_J, d = d_j,
        #   start_z = start_base_Z)
        # 
        # XR <- cbind(X, R_design)
        # R_ridge <- bdiag(zeromat_beta, R_ridge)
        # 
        # update_expansion_XR <- solve(t(XR) %*% diag_vi_pg_mean %*% XR + R_ridge, t(XR) %*% s)
        # 
        # update_expansion_bX <- update_expansion_XR[1:p.X,,drop=F]
        # 
        # update_expansion_R <- mapply(split(update_expansion_XR[-1:-p.X], rep(1:number_of_RE, d_j^2)), d_j, SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d)})
        # 
        # prop_vi_sigma_alpha <- mapply(vi_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE, 
        #                               FUN=function(Phi, R){R %*% Phi %*% t(R)})
        # 
        # if (prevent_degeneracy){
        #   update_expansion_R <- mapply(prop_vi_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE, FUN=function(prop, r){
        #     #If determinant is negative or very small, do not allow!
        #     det_prop <- det(prop)
        #     if (det_prop < 1e-10 | det(r) < 0){
        #       r <- diag(ncol(r))
        #     }
        #     return(r)
        #   })
        #   
        #   prop_vi_sigma_alpha <- mapply(vi_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE, 
        #                                 FUN=function(Phi, R){R %*% Phi %*% t(R)})
        # }
        # 
        # 
        # mapping_for_R_block <- make_mapping_alpha(update_expansion_R, px.R = TRUE)
        # update_expansion_Rblock <- prepare_T(mapping = mapping_for_R_block, levels_per_RE = g_j, num_REs = number_of_RE,
        #                   variables_per_RE = d_j, running_per_RE = breaks_for_RE, cyclical = FALSE, px.R = TRUE)
        # 
        # #all.equal(update_expansion_Rblock, bdiag(mapply(update_expansion_R, g_j, FUN=function(i,g){bdiag(lapply(1:g, FUN=function(k){i}))})))
        # 
        # update_expansion_mu <- t(M_prime) %*% vi_alpha_mean
        # 
        # update_expansion_R_logdet <- log(sapply(update_expansion_R, det))
        # 
        # 
        # prop_vi_beta_mean <- t(M_mu_to_beta) %*% bdiag(update_expansion_R) %*% update_expansion_mu +
        #   update_expansion_bX 
        # prop_vi_alpha_mean <- update_expansion_Rblock %*% (vi_alpha_mean - M_prime_one %*% update_expansion_mu)
        # 
        # 
        # #L^T L = Variance
        # #R Var R^T --->
        # # L %*% R^T
        # prop_vi_alpha_decomp <- vi_alpha_decomp %*% t(update_expansion_Rblock)
        # prop_log_det_alpha_var <- log_det_alpha_var + 2 * sum(update_expansion_R_logdet * g_j)
        # 
        # prop_variance_by_alpha_jg <- calculate_expected_outer_alpha(L = prop_vi_alpha_decomp, alpha_mu = as.vector(prop_vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
        # prop_vi_sigma_outer_alpha <- prop_variance_by_alpha_jg$outer_alpha
        # 
        # prop.ELBO <- calculate_ELBO(ELBO_type = ELBO_type,
        #   factorization_method = factorization_method,
        #   d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
        #   prior_sigma_alpha_nu = prior_sigma_alpha_nu,
        #   iw_prior_constant = iw_prior_constant,
        #   X = X, Z = Z, s = s, y = y,
        #   vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
        #   vi_sigma_alpha = prop_vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
        #   vi_sigma_outer_alpha = prop_vi_sigma_outer_alpha, 
        #   vi_beta_mean = prop_vi_beta_mean, vi_alpha_mean = prop_vi_alpha_mean,
        #   log_det_beta_var = log_det_beta_var, log_det_alpha_var = prop_log_det_alpha_var,
        #   vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = prop_vi_alpha_decomp, 
        #   choose_term = choose_term
        # )
        # 
        # if (prop.ELBO$ELBO > prior.ELBO$ELBO){
        #   #Accept the PX-VB adjustment (OSL).
        #   vi_beta_mean <- prop_vi_beta_mean
        #   vi_alpha_mean <- prop_vi_alpha_mean
        #   vi_sigma_alpha <- prop_vi_sigma_alpha
        #   vi_alpha_decomp <- prop_vi_alpha_decomp
        #   log_det_alpha_var <- prop_log_det_alpha_var
        #   variance_by_alpha_jg <- prop_variance_by_alpha_jg
        #   vi_sigma_outer_alpha <- prop_vi_sigma_outer_alpha
        #   accept.PX <- TRUE
        # }else{
        #   accept.PX <- FALSE
        # }
        # accepted_times <- accept.PX + accepted_times
        # 
        # rm(prop_vi_beta_mean, prop_vi_alpha_mean, prop_vi_sigma_alpha, prop_vi_alpha_decomp,
        #    prop_log_det_alpha_var, prop_variance_by_alpha_jg, prop_vi_sigma_outer_alpha)
        # 
        # 
        # rownames(vi_alpha_mean) <- fmt_names_Z
        
      }else if (parameter_expansion == 'none'){
        accept.PX <- TRUE
      }else{stop('Invalid PX method.')}
      
      #Adjust the terms in the ELBO calculation that are different.

      # ELBO_type <<- ELBO_type;
      #                factorization_method <<- factorization_method;
      #                d_j <<- d_j; g_j <<- g_j; prior_sigma_alpha_phi <<- prior_sigma_alpha_phi; 
      #                prior_sigma_alpha_nu <<- prior_sigma_alpha_nu;
      #                iw_prior_constant <<- iw_prior_constant;
      #                X <<- X; Z <<- Z; s <<- s; y <<- y;
      #                vi_pg_b <<- vi_pg_b; vi_pg_mean <<- vi_pg_mean; vi_pg_c <<- vi_pg_c;
      #                vi_sigma_alpha <<- vi_sigma_alpha; vi_sigma_alpha_nu <<- vi_sigma_alpha_nu; 
      #                vi_sigma_outer_alpha <<- vi_sigma_outer_alpha;
      #                vi_beta_mean <<- vi_beta_mean; vi_alpha_mean <<- vi_alpha_mean;
      #                log_det_beta_var <<- log_det_beta_var; log_det_alpha_var <<- log_det_alpha_var;
      #                vi_beta_decomp <<- vi_beta_decomp; vi_alpha_decomp <<- vi_alpha_decomp; 
      #                vi_joint_decomp <<- vi_joint_decomp; choose_term <<- choose_term;
      #                log_det_joint_var <<- log_det_joint_var; vi_r_mu <<- vi_r_mu
      
      if (accept.PX){
        final.ELBO <- calculate_ELBO(ELBO_type = ELBO_type,
          factorization_method = factorization_method,
          d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi, 
          prior_sigma_alpha_nu = prior_sigma_alpha_nu,
          iw_prior_constant = iw_prior_constant,
          X = X, Z = Z, s = s, y = y,
          vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
          vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu, 
          vi_sigma_outer_alpha = vi_sigma_outer_alpha,
          vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
          log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
          vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp, 
          vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
          log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma
        )
      }else{
        stop('PX should not fail in current version')
        # Should be no case
        # if (accept.PX){
        #   final.ELBO <- prop.ELBO
        # }else{
        #   final.ELBO <- prior.ELBO
        # }
      }
      
      if (do_timing){
        toc(quiet = verbose_time, log = T)
        tic('Final Cleanup')
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
      
      change_vi_r_mu <- vi_r_mu - lagged_vi_r_mu
      
      if (do_timing){
        toc(quiet = verbose_time, log =T)
      }
      if (debug_param){
        store_beta[it,] <- as.vector(vi_beta_mean)
      }
      change_all <- data.frame(change_alpha_mean, change_beta_mean, t(change_sigma_mean), change_alpha_var, change_beta_var, change_joint_var, change_vi_r_mu)
      if ( (max(change_all) < tolerance_parameters) | (change_elbo$ELBO > 0 & change_elbo$ELBO < tolerance_elbo) ){
        if (!quiet){
          message(paste0('Converged after ', it, ' iterations with ELBO change of ', round(change_elbo[1], 1 + abs(floor(log(tolerance_elbo)/log(10))))))
          message(paste0('The largest change in any variational parameter was ', round(max(change_all), 1 + abs(floor(log(tolerance_parameters)/log(10))))))
        }
        break
      }
      if (!quiet & (it %% print_prog == 0)){
        message(paste0('ELBO Change: ', round(change_elbo$ELBO, 10)))
        message(paste0('Other Parameter Changes: ', max(change_all)))
      }
      
      lagged_alpha_mean <- vi_alpha_mean
      lagged_beta_mean <- vi_beta_mean
      lagged_alpha_decomp <- vi_alpha_decomp
      lagged_beta_decomp <- vi_beta_decomp
      lagged_sigma_alpha <- vi_sigma_alpha
      lagged_vi_r_mu <- vi_r_mu
      lagged_ELBO <- final.ELBO
    }
    if (it == iterations){
      message(paste0('Ended without Convergence after', it, ' iterations : ELBO change of ', round(change_elbo[1], abs(floor(log(tolerance_elbo)/log(10))))))
    }
    
    if (debug_ELBO){
      d.ELBO <- with(store_ELBO, ELBO - dplyr::lag(ELBO))
      sum.ELBO <- store_ELBO
      sum.ELBO$diff <- d.ELBO
      sum.ELBO
      sum.ELBO <- summarize(group_by_at(sum.ELBO, .vars = 'step'), 
                negative = mean(.data$diff < 0, na.rm=T)) 
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
    output$family <- family
    output$control <- control
    
    if (do_timing){
      tic_log <- tictoc::tic.log(format = FALSE)
      tic_log <- data.frame(stage = sapply(tic_log, FUN=function(i){i$msg}), time = sapply(tic_log, FUN=function(i){i$toc - i$tic}), stringsAsFactors = F)
      
      tic.clear()
      tic.clearlog()
      
      tic_summary <- summarize(.data = group_by_at(.vars = 'stage', tic_log), 
                .groups = 'keep',
                n = dplyr::n(), mean = mean(.data$time), 
                min = min(.data$time), max = max(.data$time),
                total = sum(.data$time))
    }else{
      tic_summary <- NULL
    }
    if (debug_param){
      output$parameter_trajectory <- list(beta = store_beta)
    }
    if (factorization_method == 'weak'){
      output$joint <- vi_joint_decomp
    }
    if (control$return_data){
      output$data <- list(X = X, Z = Z, y = y, trials = trials)
    }
    output$formula <- formula
    output$alpha$dia.var <- unlist(lapply(variance_by_alpha_jg$variance_jg, FUN=function(i){as.vector(sapply(i, diag))}))
    output$beta$var <- t(vi_beta_decomp) %*% vi_beta_decomp
    output$beta$decomp_var <- vi_beta_decomp
    
    if (family == 'negbin'){
      output$r <- list(mu = vi_r_mu, sigma = vi_r_sigma, method = vi_r_method)
    }

    output$internal_parameters <- list(it_used = it, it_max = iterations, 
                                       missing_obs = missing_obs, N = nrow(X),
                                       names_of_RE = names_of_RE, d_j = d_j, g_j = g_j)

    output$alpha$var <- variance_by_alpha_jg$variance_jg
    output$alpha$decomp_var <- vi_alpha_decomp
    output$timing <- tic_summary
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

#' Control for vglmer
#' 
#' Provides a set of control arguments to vglmer
#' 
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
#' @param vi_r_method Type of estimate for "r"; at moment, accepts VEM (Variational EM) or fixed
#' @param vi_r_val For fixed "r", which value?
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
#' @param verbose_time Print time for each step (debugging only)
#' @param do_timing Estimate timing with tictoc
#' 
#' @importFrom checkmate assert check_double check_logical check_choice check_int check_integerish
#' @export
vglmer_control <- function(iterations = 1000,
  prior_variance = 'mean_exists', factorization_method = 'weak', 
  tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
  prevent_degeneracy = FALSE, force_whole = TRUE, verbose_time = TRUE,
  parameter_expansion = 'mean', random_seed = 1, do_timing = FALSE,
  debug_param = FALSE, return_data = FALSE, linpred_method = 'joint',
  vi_r_method = 'VEM', vi_r_val = NA,
  debug_ELBO = FALSE, print_prog = NULL, quiet = T, init = 'EM'
){
  #use checkmate package to verify arguments
  assert(
    check_integerish(iterations, lower = 1),
    check_double(c(tolerance_elbo, tolerance_parameters), len = 2, lower = 0), 
    check_logical(c(prevent_degeneracy, force_whole, verbose_time, do_timing, 
        debug_param, return_data, debug_ELBO, quiet), len = 8),
    check_choice(factorization_method, c("weak", "strong", "partial")),
    check_choice(prior_variance, c('mean_exists', 'jeffreys', 'mcmcglmm', 'mvD', 'limit', 'uniform')),
    check_choice(linpred_method, c('joint', 'cyclical', 'solve_normal')),
    check_choice(vi_r_method, c('VEM', 'fixed', 'Laplace', 'delta', 'temp')),
    check_double(vi_r_val, all.missing = TRUE),
    check_int(print_prog, null.ok = TRUE),
    check_choice(init, c('EM', 'random', 'zero')),
    check_double(random_seed),
    combine = 'and')

  if (vi_r_method == 'fixed' & is.na(vi_r_val)){
    stop('vi_r_val must not be NA if vi_r_method = "fixed"')
  }
  
  output <- namedList(iterations, prior_variance, factorization_method, 
                 tolerance_elbo, tolerance_parameters,
                 prevent_degeneracy, force_whole, verbose_time,
                 parameter_expansion, random_seed, do_timing, debug_param, return_data,
                 linpred_method, vi_r_method, vi_r_val, debug_ELBO, print_prog, quiet, init)
  
  class(output) <- c('vglmer_control')
  return(output)
}

# Simple function to create named list 
# https://stackoverflow.com/questions/16951080/can-lists-be-created-that-name-themselves-based-on-input-object-names
# Identical to version used in lme4:::namedList, see also loo::nlist
#' @importFrom stats setNames
namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}
