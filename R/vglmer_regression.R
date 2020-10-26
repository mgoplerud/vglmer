#' Variational Inference for Non-Linear Hierarchical Models
#'
#' Estimate a hierarchical model using mean-field variational
#' inference. Accepts standard syntax to glmer: \code{y ~ X + (1 + Z | g)}.
#' Options are described below.
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
#' sim_data <- data.frame(
#'   x = rnorm(100),
#'   y = rbinom(100, 1, 0.5),
#'   g = sample(letters, 100, replace = TRUE)
#' )
#'
#' # Run with defaults
#' est_vglmer <- vglmer(y ~ x + (x | g), data = sim_data, family = "binomial")
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
#' vglmer(y ~ x + (x | g),
#'   data = sim_data,
#'   control = vglmer_control(factorization_method = "strong"),
#'   family = "binomial"
#' )
#' }
#' @return Returns an object of class vglmer: See the available methods (e.g.
#'   \code{coef}) using \code{methods(class="vglmer")}. A few of the internal
#'   outputs are described below.
#' \describe{
#'
#' \item{alpha}{Contains the posterior distribution of each random effect.
#' \code{mean} contains the posterior means; \code{dia.var} contains the
#' variance of each random effect. \code{var} contains the variance-covariance
#' matrix of each random effect (j,g). \code{decomp_var} contains a matrix L
#' such that L^T L the full variance of the entire set of random effects.}
#'
#' \item{sigma}{Contains the posterior distribution of each random effect error
#' variance.}
#'
#' \item{ELBO}{Contains the ELBO at the termination of the algorithm.}
#'
#' \item{ELBO_trajectory}{A data.frame tracking the ELBO per iteration.}
#'
#' }
#' @importFrom dplyr select group_by group_by_at summarize n lag
#' @importFrom lme4 mkReTrms findbars subbars
#' @importFrom stats model.response model.matrix model.frame rnorm rWishart
#'   qlogis optim residuals lm
#' @importFrom rlang .data
#' @importFrom graphics plot
#' @importFrom checkmate assert assert_formula assert_choice
#'   check_data_frame
#' @useDynLib vglmer
#' @export
vglmer <- function(formula, data, family, control = vglmer_control()) {

  # Verify integrity of parameter arguments
  checkdf <- check_data_frame(data, null.ok = TRUE)
  if (checkdf != TRUE) {
    warning(paste0("data is not a data.frame? Behavior may be unexpected: ", checkdf))
  }
  assert_formula(formula)
  # Delete the missing data
  # (i.e. sub out the random effects, do model.frame)
  #
  nobs_init <- nrow(data)
  
  if (control$parameter_expansion %in% c('translation', 'diagonal')){
    parsed_RE_groups <- get_RE_groups(formula = formula, data = data)
  }
  data <- model.frame(subbars(formula), data)
  nobs_complete <- nrow(data)
  missing_obs <- nobs_init - nobs_complete
  if (length(missing_obs) == 0) {
    missing_obs <- "??"
  }

  #Extract the Outcome
  y <- model.response(data)
  if (is.matrix(y)){
    N <- nrow(y)
    rownames(y) <- NULL
  }else{
    N <- length(y)
    y <- as.vector(y)
    names(y) <- NULL
  }
  
  
  assert_choice(family, c("negbin", "binomial", "linear"))
  if (!inherits(control, "vglmer_control")) {
    stop("control must be object from vglmer_control().")
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

  if (do_timing) {
    if (!requireNamespace("tictoc", quietly = TRUE)) {
      stop("tictoc must be installe to do timing")
    }
    tic <- tictoc::tic
    toc <- tictoc::toc
    tic.clear <- tictoc::tic.clear
    tic.clearlog <- tictoc::tic.clearlog

    tic.clear()
    tic.clearlog()
    tic("Prepare Model")
  }
  if (!(factorization_method %in% c("weak", "strong", "partial"))) {
    stop("factorization_method must be in weak, strong, or partial.")
  }
  if (is.null(print_prog)) {
    print_prog <- max(c(1, floor(iterations / 20)))
  }
  if (!(family %in% c("binomial", "negbin", "linear"))) {
    stop('family must be "linear", "binomial", "negbin".')
  }

  if (family == "binomial") {
    if (is.matrix(y)) {
      # if (!(class(y) %in% c('numeric', 'integer'))){
      if (min(y) < 0) {
        stop("Negative Numbers not Permitted in outcome")
      }
      is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
      if (any(is.wholenumber(y) == FALSE)) {
        if (control$force_whole) {
          stop("If force_whole = TRUE, must provide whole numbers")
        } else {
          warning("Non-integer numbers in y")
        }
      }
      # Total trials (Success + Failure)
      trials <- rowSums(y)
      rownames(trials) <- NULL
      # Successes
      y <- y[, 1]
      rownames(y) <- NULL
    } else {
      if (!all(y %in% c(0, 1)) & family == "binomial") {
        stop("Only {0,1} outcomes permitted for numeric y.")
      }
      trials <- rep(1, length(y))
    }
  } else if (family == 'negbin') {
    if (!(class(y) %in% c("numeric", "integer"))) {
      stop("Must provide vector of numbers with negbin.")
    }

    if (min(y) < 0) {
      stop("Negative numbers not Permitted in outcome")
    }

    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
    if (any(is.wholenumber(y) == FALSE)) {
      if (control$force_whole) {
        stop("If force_whole = TRUE, must provide whole numbers")
      } else {
        warning("Non-integer numbers in y")
      }
    }
  } else if (family == 'linear') {
    y <- as.numeric(y)
    #Do nothing if linear
  } else {
    stop('family is invalid.')
  }

  if (family %in% c("binomial", "linear")) {
    ELBO_type <- "augmented"
  } else if (family == "negbin") {
    ELBO_type <- "profiled"
  } else {
    stop("Check ELBO_type")
  }

  # Extract X (FE design matrix)
  X <- model.matrix(nobars(formula), data = data)

  # Extract the Z (Random Effect) design matrix.
  mk_Z <- mkReTrms(findbars(formula), data, reorder.terms = FALSE, reorder.vars = FALSE)
  Z <- t(mk_Z$Zt)

  p.X <- ncol(X)
  p.Z <- ncol(Z)

  ####
  # Process the REs to get various useful terms.
  ####
  # RE names and names of variables included for each.
  names_of_RE <- mk_Z$cnms

  number_of_RE <- length(mk_Z$Gp) - 1

  if (number_of_RE < 1) {
    stop("Need to provide at least one random effect...")
  }
  # The position that demarcates each random effect.
  # That is, breaks_for_RE[2] means at that position + 1 does RE2 start.
  breaks_for_RE <- c(0, cumsum(diff(mk_Z$Gp)))
  # Dimensionality of \alpha_{j,g}, i.e. 1 if random intercept
  # 2 if random intercept + random slope
  d_j <- lengths(names_of_RE)
  # Number of GROUPs for each random effect.
  g_j <- diff(mk_Z$Gp) / d_j

  # Empty vector to build the formatted names for each random effect.
  fmt_names_Z <- c()
  init_Z_names <- colnames(Z)
  for (v in 1:number_of_RE) {
    name_of_effects_v <- names_of_RE[[v]]

    mod_name <- rep(name_of_effects_v, g_j[v])

    levels_of_re <- init_Z_names[(1 + breaks_for_RE[v]):breaks_for_RE[v + 1]]

    fmt_names_Z <- c(fmt_names_Z, paste0(names(names_of_RE)[v], " @ ", mod_name, " @ ", levels_of_re))
  }
  colnames(Z) <- fmt_names_Z
  cyclical_pos <- lapply(1:number_of_RE, FUN = function(i) {
    seq(breaks_for_RE[i] + 1, breaks_for_RE[i + 1])
  })

  M.names <- cbind(unlist(mapply(names_of_RE, g_j, FUN = function(i, j) {
    rep(i, j)
  })))
  M <- cbind(match(M.names[, 1], colnames(X)), rep(1 / g_j, d_j * g_j))
  M <- sparseMatrix(i = 1:ncol(Z), j = M[, 1], x = M[, 2], dims = c(ncol(Z), ncol(X)))

  M_prime.names <- paste0(rep(names(names_of_RE), g_j * d_j), " @ ", M.names)
  M_prime <- cbind(match(M_prime.names, unique(M_prime.names)), rep(1 / g_j, d_j * g_j))
  M_prime <- sparseMatrix(i = 1:ncol(Z), j = M_prime[, 1], x = M_prime[, 2])
  colnames(M_prime) <- unique(M_prime.names)

  M_prime_one <- M_prime
  M_prime_one@x <- rep(1, length(M_prime_one@x))

  stopifnot(identical(paste0(rep(names(names_of_RE), d_j), " @ ", unlist(names_of_RE)), colnames(M_prime)))

  M_mu_to_beta <- sparseMatrix(i = 1:sum(d_j), j = match(unlist(names_of_RE), colnames(X)), x = 1, dims = c(sum(d_j), p.X))
  colnames(M_mu_to_beta) <- colnames(X)
  rownames(M_mu_to_beta) <- colnames(M_prime)
  # List of Lists
  # Outer list: one for RE
  # Inner List: One for each GROUP with its row positions.
  outer_alpha_RE_positions <- mapply(d_j, g_j, breaks_for_RE[-length(breaks_for_RE)], SIMPLIFY = FALSE, FUN = function(a, b, m) {
    split(m + seq(1, a * b), rep(1:b, each = a))
  })

  if (anyDuplicated(unlist(outer_alpha_RE_positions)) != 0 | max(unlist(outer_alpha_RE_positions)) != ncol(Z)) {
    stop("Issue with greating OA positions")
  }
  ####
  # Prepare Initial Values
  ###

  vi_sigmasq_prior_a <- 0
  vi_sigmasq_a <- (length(y) + sum(d_j * g_j))/2 + vi_sigmasq_prior_a
  
  vi_sigmasq_prior_b <- 1
  vi_sigmasq_b <- 0
  
  if (family == "linear") {
    s <- y
    vi_pg_b <- 1
    vi_r_mu <- 0
    vi_r_sigma <- 0
    vi_r_mean <- 0
    
    vi_sigmasq_b <- sum(residuals(lm(y ~ 1))^2)/2 + vi_sigmasq_prior_b
    

    choose_term <- -length(y)/2 * log(2 * pi)
  } else if (family == "binomial") {
    s <- y - trials / 2
    vi_pg_b <- trials
    vi_r_mu <- 0
    vi_r_mean <- 0
    vi_r_sigma <- 0
    choose_term <- sum(lchoose(n = round(trials), k = round(y)))
  } else {
    # Initialize
    if (vi_r_method == "fixed") {
      vi_r_mu <- vi_r_val
      vi_r_mean <- exp(vi_r_mu)
      vi_r_sigma <- 0
    } else if (vi_r_method == "VEM") {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Install MASS for negbin")
      }
      vi_r_mean <- MASS::glm.nb(y ~ 1)$theta
      vi_r_mu <- log(vi_r_mean)
      vi_r_sigma <- 0
    } else if (vi_r_method %in% c("Laplace", "delta")) {
      init_r <- optim(
        par = 0, fn = VEM.PELBO.r, method = "L-BFGS", hessian = T,
        control = list(fnscale = -1), y = y, psi = rep(log(mean(y)), length(y)), zVz = 0
      )
      vi_r_mu <- init_r$par
      vi_r_sigma <- as.numeric(-1 / init_r$hessian)

      vi_r_mean <- exp(vi_r_mu + vi_r_sigma / 2)
    } else {
      stop("vi_r_method must be VEM, fixed, delta, or Laplace.")
    }
    s <- (y - vi_r_mean) / 2
    vi_pg_b <- y + vi_r_mean

    choose_term <- NA #-sum(lgamma(y + 1)) - sum(y) * log(2)
  }

  # Initalize variational parameters.
  # Note that we keep a sparse matrix or lowertri such that
  # t(vi_beta_decomp) %*% vi_beta_decomp = VARIANCE

  vi_beta_decomp <- Diagonal(x = rep(0, ncol(X)))
  vi_alpha_decomp <- Diagonal(x = rep(0, ncol(Z)))

  vi_sigma_alpha_nu <- g_j

  prior_variance <- control$prior_variance
  do_huangwand <- FALSE
  vi_a_APRIOR_jp <- vi_a_nu_jp <- vi_a_a_jp <- vi_a_b_jp <- NULL
  prior_sigma_alpha_nu <- prior_sigma_alpha_phi <- NULL
  
  if (prior_variance == 'hw') {
    do_huangwand <- TRUE
    vi_a_nu_jp <- rep(3, length(d_j))
    names(vi_a_nu_jp) <- names(names_of_RE)
    vi_a_APRIOR_jp <- lapply(d_j, FUN=function(i){rep(2.5, i)})
    vi_a_a_jp <- mapply(d_j, vi_a_nu_jp, SIMPLIFY = FALSE, 
                        FUN=function(i,nu){1/2 * (nu + rep(i, i))})
    vi_a_b_jp <- lapply(vi_a_APRIOR_jp, FUN=function(i){1/i^2})
  } else if (prior_variance == 'kn') {
    
    if (family != 'binomial'){stop('Set up kn for non-binomial.')}
    fe_only <- EM_prelim_logit(X = drop0(X), Z = NULL, s = s, pg_b =  vi_pg_b, iter = 15)
    weight_p <- plogis(as.vector(X %*% fe_only$beta))
    weight_p <- weight_p * (1-weight_p) * trials
        
    weight_kn <- Diagonal(x = sqrt(weight_p)) %*% Z
    weight_kn <- lapply(outer_alpha_RE_positions, FUN=function(i){
      D_matrix <- Reduce('+', lapply(i, FUN=function(j){
        crossprod(weight_kn[,j])
      }))
      D_matrix <- solve(1/length(i) * D_matrix)
    })
    prior_sigma_alpha_nu <- d_j
    prior_sigma_alpha_phi <- lapply(weight_kn, FUN=function(i){i * nrow(i)})

  } else if (prior_variance == "jeffreys") {
    prior_sigma_alpha_nu <- rep(0, number_of_RE)
    prior_sigma_alpha_phi <- lapply(d_j, FUN = function(i) {
      diag(x = 0, nrow = i, ncol = i)
    })
  } else if (prior_variance == "mcmcglmm") {
    prior_sigma_alpha_nu <- rep(0, number_of_RE)
    prior_sigma_alpha_phi <- lapply(d_j, FUN = function(i) {
      diag(x = 0, nrow = i, ncol = i)
    })
  } else if (prior_variance == "mvD") {
    prior_sigma_alpha_nu <- -d_j
    prior_sigma_alpha_phi <- lapply(d_j, FUN = function(i) {
      diag(x = 0, nrow = i, ncol = i)
    })
  } else if (prior_variance == "mean_exists") {
    prior_sigma_alpha_nu <- d_j + 1 # Ensures the mean exists...
    prior_sigma_alpha_phi <- lapply(d_j, FUN = function(i) {
      diag(x = 1, nrow = i, ncol = i)
    })
  } else if (prior_variance == "limit") {
    prior_sigma_alpha_nu <- d_j - 1
    prior_sigma_alpha_phi <- lapply(d_j, FUN = function(i) {
      diag(x = 0, nrow = i, ncol = i)
    })
  } else if (prior_variance == "uniform") {
    prior_sigma_alpha_nu <- -(d_j + 1)
    prior_sigma_alpha_phi <- lapply(d_j, FUN = function(i) {
      diag(x = 0, nrow = i, ncol = i)
    })
  } else {
    stop("Options for prior variance are jeffreys and mean_exists")
  }

  if (do_huangwand){
    iw_prior_constant <- mapply(vi_a_nu_jp, d_j,
      FUN = function(nu, d) {
        nu <- nu + d - 1
        return(- (nu * d) / 2 * log(2) - multi_lgamma(a = nu / 2, p = d))
      }
    )
    vi_sigma_alpha_nu <- vi_sigma_alpha_nu + vi_a_nu_jp + d_j - 1
  }else{
    # normalizingly constant for wishart to make ELBO have right value to compare models.
    iw_prior_constant <- mapply(prior_sigma_alpha_nu, prior_sigma_alpha_phi,
                                FUN = function(nu, Phi) {
                                  if (nu <= (ncol(Phi) - 1)) {
                                    return(0)
                                  } else {
                                    return(make_log_invwishart_constant(nu, Phi))
                                  }
                                }
    )
    
    vi_sigma_alpha_nu <- vi_sigma_alpha_nu + prior_sigma_alpha_nu
  }

  if (control$init == "EM") {
    if (family == "linear"){
      jointXZ <- cbind(X,Z)
      EM_init <- LinRegChol(X = jointXZ, 
        omega = sparseMatrix(i = 1:nrow(jointXZ), j = 1:nrow(jointXZ), x = 1), 
        y = y, prior_precision = sparseMatrix(i = 1:ncol(jointXZ), j = 1:ncol(jointXZ), x = 1/4))$mean
      EM_init <- list('beta' = EM_init[1:ncol(X)], 'alpha' = EM_init[-1:-ncol(X)])
      rm(jointXZ)
    } else if (family == "negbin") {
      EM_init <- EM_prelim_nb(X = X, Z = Z, y = y, est_r = exp(vi_r_mu), iter = 15, ridge = 4)
    } else {
      EM_init <- EM_prelim_logit(X = X, Z = Z, s = s, pg_b = vi_pg_b, iter = 15, ridge = 4)
    }

    vi_beta_mean <- matrix(EM_init$beta)
    vi_alpha_mean <- matrix(EM_init$alpha)

    vi_sigma_alpha <- calculate_expected_outer_alpha(
      alpha_mu = matrix(EM_init$alpha),
      L = sparseMatrix(i = 1, j = 1, x = 0, dims = rep(ncol(Z), 2)),
      re_position_list = outer_alpha_RE_positions
    )
    if (!do_huangwand){
      vi_sigma_alpha <- mapply(vi_sigma_alpha$outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN = function(i, j) {
        i + j
      })
    }else{
      #Update Inverse-Wishart
      vi_a_b_jp <<- vi_a_b_jp
      vi_a_nu_jp <<- vi_a_nu_jp
      vi_a_a_jp <<- vi_a_a_jp
      vi_sigma_alpha <<- vi_sigma_alpha
      vi_a_APRIOR_jp <<- vi_a_APRIOR_jp
      vi_sigma_alpha <- mapply(vi_sigma_alpha$outer_alpha, vi_a_a_jp, 
       vi_a_b_jp, vi_a_nu_jp, SIMPLIFY = FALSE, FUN = function(i, tilde.a, tilde.b, nu) {
         i + Diagonal(n = nrow(i))
       })
      #Update a_{j,p}
      diag_Einv_sigma <- mapply(vi_sigma_alpha, 
          vi_sigma_alpha_nu, d_j, SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
            inv_phi <- solve(phi)
            sigma.inv <- nu * inv_phi
            return(diag(sigma.inv))
          })
      # vi_a_b_jp <- lapply(vi_sigma_alpha, FUN=function(i){diag(solve(i))})
      vi_a_b_jp <- mapply(vi_a_nu_jp, vi_a_APRIOR_jp, diag_Einv_sigma,
        SIMPLIFY = FALSE,
        FUN=function(nu, APRIOR, diag_j){
          1/APRIOR^2 + nu * diag_j
        })
    }
  } else if (control$init == "random") {
    set.seed(control$random_seed)
    vi_beta_mean <- rnorm(ncol(X))
    vi_alpha_mean <- rep(0, ncol(Z))

    if (do_huangwand){
      vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN = function(d, g) {
        rWishart(n = 1, df = g, Sigma = diag(d))[, , 1]
      })
    }else{
      vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN = function(d, g) {
        rWishart(n = 1, df = g, Sigma = diag(d))[, , 1]
      })
    }
  } else if (control$init == "zero") {
    vi_beta_mean <- rep(0, ncol(X))

    if (family == "binomial") {
      vi_beta_mean[1] <- qlogis(sum(y) / sum(trials))
    } else if (family == "negbin") {
      vi_beta_mean[1] <- log(mean(y))
    } else if (family == 'linear'){
      vi_beta_mean[1] <- mean(y)
    } else {
      stop('Set up init')
    }

    vi_alpha_mean <- rep(0, ncol(Z))

    vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN = function(d, g) {
      diag(x = 1, ncol = d, nrow = d)
    })
    # if (do_huangwand){stop('Setup init for zero')}
  } else {
    stop("Provide init = EM or random")
  }

  zero_mat <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(X), ncol(X)))
  zero_mat <- drop0(zero_mat)

  if (factorization_method == "weak") {
    vi_joint_decomp <- bdiag(vi_beta_decomp, vi_alpha_decomp)
    joint.XZ <- cbind(X, Z)
    log_det_beta_var <- log_det_alpha_var <- NULL
  } else {
    vi_joint_decomp <- NULL
    log_det_joint_var <- NULL
  }

  # if (!quiet){warning('Check vi sigma alpha nu')}

  # Create mapping for this to allow sparse implementations.
  mapping_sigma_alpha <- make_mapping_alpha(vi_sigma_alpha)

  running_log_det_alpha_var <- rep(NA, number_of_RE)

  lagged_alpha_mean <- rep(-Inf, ncol(Z))
  lagged_beta_mean <- rep(-Inf, ncol(X))
  lagged_sigma_alpha <- vi_sigma_alpha
  if (factorization_method == "weak") {
    lagged_joint_decomp <- vi_joint_decomp
  } else {
    lagged_alpha_decomp <- vi_alpha_decomp
    lagged_beta_decomp <- vi_beta_decomp
  }
  lagged_vi_r_mu <- -Inf
  lagged_vi_sigmasq_a <- lagged_vi_sigmasq_b <- -Inf
  lagged_ELBO <- -Inf
  accepted_times <- NA

  skip_translate <- FALSE
  if (parameter_expansion %in%  c("translation", "diagonal")) {
    accepted_times <- 0
    attempted_expansion <- 0
    
    stationary_rho <- do.call('c', lapply(d_j, FUN=function(i){as.vector(diag(x = i))}))
    diag_rho <- which(stationary_rho == 1)
    zeromat_beta <- drop0(Diagonal(x = rep(0, ncol(X))))

    # parsed_RE_groups <- get_RE_groups(formula = formula, data = data)
    parsed_RE_groups <<- parsed_RE_groups
    
    mapping_new_Z <- do.call('cbind', parsed_RE_groups$design)
    
    mapping_J <- split(1:sum(d_j^2), rep(1:length(d_j), d_j^2))
    mapping_J <- lapply(mapping_J, FUN=function(i){i-1})
    mapping_J <- sapply(mapping_J, min)

    mapping_to_re <- parsed_RE_groups$factor
    mapping_to_re <- purrr::array_branch(do.call('cbind', mapping_to_re), margin = 1)
    mapping_to_re <- lapply(mapping_to_re, FUN=function(i){
      mapply(outer_alpha_RE_positions, i, SIMPLIFY = FALSE, FUN=function(a,b){a[[b]]})
    })
    Mmap <- do.call('rbind', lapply(mapping_to_re, FUN=function(i){as.integer(sapply(i, min))}))

    start_base_Z <- cumsum(c(0,d_j))[-(number_of_RE+1)]
    names(start_base_Z) <- NULL

    rm(parsed_RE_groups, mapping_to_re)
    
    store_re_id <- store_id <- list()
    for (j in 1:number_of_RE){
      store_re_id_j <- store_id_j <- list()
      for (jprime in 1:j){
        # print(c(j, jprime))
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
    store_design <- parsed_RE_groups$design
    
    gc()
  }
  store_parameter_traj <- store_vi <- store_ELBO <- data.frame()

  if (debug_param) {
    store_beta <- array(NA, dim = c(iterations, ncol(X)))
  }
  if (do_timing) {
    toc(quiet = verbose_time, log = TRUE)
    tic.clear()
  }
  ## Begin VI algorithm:
  if (!quiet) {
    message("Begin Regression")
  }
  do_SQUAREM <- TRUE
  if (do_SQUAREM){
    squarem_list <- list()
    squarem_counter <- 1
  }
  for (it in 1:iterations) {
    if (it %% print_prog == 0) {
      cat(".")
    }
    ###
    ## Polya-Gamma Updates
    ###
    # Get the x_i^T Var(beta) x_i terms.
    if (do_timing) {
      tic("Update PG")
    }

    if (factorization_method == "weak") {
      # joint_quad <- rowSums( (joint.XZ %*% t(vi_joint_decomp))^2 )
      # vi_joint_decomp <<- vi_joint_decomp
      # joint.XZ <<- joint.XZ
      joint_quad <- cpp_zVz(Z = joint.XZ, V = as(vi_joint_decomp, "dgCMatrix")) 
      if (family == 'negbin'){
        joint_quad <- joint_quad + vi_r_sigma
      }
      vi_pg_c <- sqrt(as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean - vi_r_mu)^2 + joint_quad)
    } else {
      # vi_beta_decomp <<- vi_beta_decomp
      # vi_alpha_decomp <<- vi_alpha_decomp
      # X <<- X
      # Z <<- Z
      # beta_quad <- cpp_zVz(Z = drop0(X), V = make_dgC(vi_beta_decomp))
      # alpha_quad <- cpp_zVz(Z = Z, V = make_dgC(vi_alpha_decomp))
      beta_quad <- rowSums((X %*% t(vi_beta_decomp))^2)
      alpha_quad <- rowSums((Z %*% t(vi_alpha_decomp))^2)
      joint_var <- beta_quad + alpha_quad
      if (family == 'negbin'){
        joint_var <- joint_var + vi_r_sigma
      }
      vi_pg_c <- sqrt(as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean - vi_r_mu)^2 + joint_var)
    }
    if (family == 'linear'){
      vi_pg_mean <- rep(1, nrow(X))
      diag_vi_pg_mean <- sparseMatrix(i = 1:N, j = 1:N, x = vi_pg_mean)
    }else{
      vi_pg_mean <- vi_pg_b / (2 * vi_pg_c) * tanh(vi_pg_c / 2)
      diag_vi_pg_mean <- sparseMatrix(i = 1:N, j = 1:N, x = vi_pg_mean)
    }
    if (debug_ELBO & it != 1) {
      debug_ELBO.1 <- calculate_ELBO(family = family,
        ELBO_type = ELBO_type,
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
        log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma,
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
      )
    }

    if (do_timing) {
      toc(quiet = verbose_time, log = TRUE)
      tic("Prepare Sigma")
    }

    # Process Sigma_j for manipulation
    # if Sigma_{j} is InverseWishart(a,Phi)
    # Then E[Sigma^{-1}_j] = a * Phi^{-1}
    if (factorization_method == "strong") {
      cyclical_T <- TRUE
    } else {
      cyclical_T <- FALSE
    }
    inv_mapping_alpha <- mapply(vi_sigma_alpha_nu, lapply(vi_sigma_alpha, solve),
      SIMPLIFY = FALSE, FUN = function(a, b) {
        a * b
      }
    )
    inv_mapping_alpha <- make_mapping_alpha(inv_mapping_alpha)

    Tinv <- prepare_T(
      mapping = inv_mapping_alpha, levels_per_RE = g_j, num_REs = number_of_RE,
      variables_per_RE = d_j, running_per_RE = breaks_for_RE, cyclical = cyclical_T
    )
    if (!cyclical_T) {
      Tinv <- as(Tinv, "dgCMatrix")
    } else {
      Tinv <- lapply(Tinv, FUN = function(i) {
        as(i, "dgCMatrix")
      })
    }
    if (do_timing) {
      toc(quiet = verbose_time, log = T)
      tic("Update Beta")
    }
    if (factorization_method == "weak") {
      ## Update <beta, alpha> jointly
      chol.update.joint <- LinRegChol(
        X = joint.XZ, omega = diag_vi_pg_mean,
        prior_precision = bdiag(zero_mat, Tinv),
        y = s + vi_pg_mean * vi_r_mu
      )
      Pmatrix <- sparseMatrix(i = 1:ncol(joint.XZ), j = 1 + chol.update.joint$Pindex, x = 1)


      vi_joint_decomp <- drop0(solve(chol.update.joint$origL) %*% t(Pmatrix))

      vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
      vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))

      vi_alpha_decomp <- vi_joint_decomp[, -1:-p.X, drop = F]
      vi_beta_decomp <- vi_joint_decomp[, 1:p.X, drop = F]

      log_det_joint_var <- -2 * sum(log(diag(chol.update.joint$origL)))
    } else if (factorization_method == "partial") {
      if (linpred_method == "cyclical") {
        # Do not run except as backup
        # ###Non optimized
        # precision_beta <- t(X) %*% diag_vi_pg_mean %*% X
        # nonopt_beta <- solve(precision_beta, t(X) %*% (s - diag_vi_pg_mean %*% Z %*% vi_alpha_mean))
        # precision_alpha <- t(Z) %*% diag_vi_pg_mean %*% Z + Tinv
        # nonopt_alpha <- solve(precision_alpha, t(Z) %*% (s - diag_vi_pg_mean %*% X %*% nonopt_beta))

        chol.update.beta <- LinRegChol(
          X = as(X, "sparseMatrix"), omega = diag_vi_pg_mean, prior_precision = zero_mat,
          y = as.vector(s - diag_vi_pg_mean %*% Z %*% vi_alpha_mean)
        )
        Pmatrix <- sparseMatrix(i = 1:p.X, j = 1 + chol.update.beta$Pindex, x = 1)

        # P origL oriL^T P^T = PRECISION
        # t(decompVar) %*%  decompVar = VARIANCE = (origL^{-1} t(P))^T (origL^{-1} t(P))

        vi_beta_decomp <- solve(chol.update.beta$origL) %*% t(Pmatrix)
        vi_beta_mean <- chol.update.beta$mean
        log_det_beta_var <- -2 * sum(log(diag(chol.update.beta$origL)))

        chol.update.alpha <- LinRegChol(
          X = Z, omega = diag_vi_pg_mean, prior_precision = Tinv,
          y = as.vector(s - diag_vi_pg_mean %*% X %*% vi_beta_mean)
        )
        Pmatrix <- sparseMatrix(i = 1:p.Z, j = 1 + chol.update.alpha$Pindex, x = 1)

        vi_alpha_decomp <- solve(chol.update.alpha$origL) %*% t(Pmatrix)
        vi_alpha_decomp <- drop0(vi_alpha_decomp)
        vi_alpha_mean <- chol.update.alpha$mean
        log_det_alpha_var <- -2 * sum(log(diag(chol.update.alpha$origL)))

        vi_beta_mean <- Matrix(vi_beta_mean, dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(vi_alpha_mean, dimnames = list(fmt_names_Z, NULL))
      } else if (linpred_method == "joint") {
        joint.XZ <- cbind(X, Z)
        chol.update.joint <- LinRegChol(
          X = joint.XZ, omega = diag_vi_pg_mean, prior_precision = bdiag(zero_mat, Tinv),
          y = s + vi_pg_mean * vi_r_mu
        )
        vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))

        vi_beta_decomp <- solve(t(chol(as.matrix(t(X) %*% diag_vi_pg_mean %*% X))))
        log_det_beta_var <- 2 * sum(log(diag(vi_beta_decomp)))

        chol.update.alpha <- LinRegChol(
          X = Z, omega = diag_vi_pg_mean, prior_precision = Tinv,
          y = s + vi_pg_mean * vi_r_mu
        )
        Pmatrix <- sparseMatrix(i = 1:p.Z, j = 1 + chol.update.alpha$Pindex, x = 1)

        vi_alpha_decomp <- solve(chol.update.alpha$origL) %*% t(Pmatrix)
        vi_alpha_decomp <- drop0(vi_alpha_decomp)
        log_det_alpha_var <- -2 * sum(log(diag(chol.update.alpha$origL)))
      } else {
        stop("Invalid linpred method for partial scheme")
      }
    } else if (factorization_method == "strong") {
      running_log_det_alpha_var <- rep(NA, number_of_RE)
      vi_alpha_decomp <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(p.Z, 2))

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


      if (linpred_method == "joint") {
        joint.XZ <- cbind(X, Z)
        chol.update.joint <- LinRegChol(X = joint.XZ, omega = diag_vi_pg_mean, prior_precision = bdiag(zero_mat, bdiag(Tinv)), y = s + vi_pg_mean * vi_r_mu)
        vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))

        vi_beta_decomp <- solve(t(chol(as.matrix(t(X) %*% diag_vi_pg_mean %*% X))))
        log_det_beta_var <- 2 * sum(log(diag(vi_beta_decomp)))
        #-log(det(t(X) %*% diag_vi_pg_mean %*% X))

        running_log_det_alpha_var <- rep(NA, number_of_RE)

        for (j in 1:number_of_RE) {
          index_j <- cyclical_pos[[j]]
          Z_j <- Z[, index_j, drop = F]
          prec_j <- t(Z_j) %*% diag_vi_pg_mean %*% Z_j + Tinv[[j]]

          chol_var_j <- solve(t(chol(prec_j)))
          running_log_det_alpha_var[j] <- 2 * sum(log(diag(chol_var_j)))

          vi_alpha_decomp[index_j, index_j] <- as(chol_var_j, "dgTMatrix")
        }

        log_det_alpha_var <- sum(running_log_det_alpha_var)
      } else if (linpred_method == "solve_normal") {
        bind_rhs_j <- list()
        bind_lhs_j <- list()

        for (j in 1:number_of_RE) {
          index_j <- cyclical_pos[[j]]
          Z_j <- Z[, index_j, drop = F]
          Z_negj <- Z[, -index_j, drop = F]
          prec_j <- t(Z_j) %*% diag_vi_pg_mean %*% Z_j + Tinv[[j]]

          chol_prec_j <- t(chol(prec_j))
          chol_var_j <- solve(chol_prec_j)

          mod_j <- solve(prec_j)

          term_j <- mod_j %*% t(Z_j) %*% diag_vi_pg_mean %*% Z
          term_j[, index_j, drop = F] <- Diagonal(n = ncol(Z_j))
          term_j <- cbind(term_j, mod_j %*% t(Z_j) %*% diag_vi_pg_mean %*% X)

          bind_lhs_j[[j]] <- term_j
          bind_rhs_j[[j]] <- mod_j %*% t(Z_j) %*% s

          running_log_det_alpha_var[j] <- 2 * sum(log(diag(chol_var_j)))
          vi_alpha_decomp[index_j, index_j] <- as(chol_var_j, "dgTMatrix")
        }

        log_det_alpha_var <- sum(running_log_det_alpha_var)

        bind_lhs_j <- drop0(do.call("rbind", bind_lhs_j))
        bind_rhs_j <- do.call("rbind", bind_rhs_j)

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
        # print(cbind(bind_solution, rbind(vi_alpha_mean, vi_beta_mean)))
        #
        vi_beta_mean <- Matrix(bind_solution[-1:-ncol(Z)], dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(bind_solution[1:ncol(Z)], dimnames = list(fmt_names_Z, NULL))
      } else if (linpred_method == "cyclical") {
        for (j in 1:number_of_RE) {
          index_j <- cyclical_pos[[j]]
          Z_j <- Z[, index_j, drop = F]
          Z_negj <- Z[, -index_j, drop = F]

          chol.j <- LinRegChol(
            X = Z_j, omega = diag_vi_pg_mean, prior_precision = Tinv[[j]],
            y = as.vector(s + vi_pg_mean * vi_r_mu - diag_vi_pg_mean %*% (X %*% vi_beta_mean + Z_negj %*% vi_alpha_mean[-index_j]))
          )
          vi_alpha_mean[index_j] <- chol.j$mean

          Pmatrix <- sparseMatrix(i = 1:ncol(Z_j), j = 1 + chol.j$Pindex, x = 1)

          running_log_det_alpha_var[j] <- -2 * sum(log(diag(chol.j$origL)))
          vi_alpha_decomp[index_j, index_j] <- solve(chol.j$origL) %*% t(Pmatrix)
        }

        # vi_alpha_decomp <- bdiag(vi_alpha_decomp)
        log_det_alpha_var <- sum(running_log_det_alpha_var)

        chol.update.beta <- LinRegChol(
          X = as(X, "sparseMatrix"), omega = diag_vi_pg_mean, prior_precision = zero_mat,
          y = as.vector(s + vi_pg_mean * vi_r_mu - diag_vi_pg_mean %*% Z %*% vi_alpha_mean)
        )
        Pmatrix <- sparseMatrix(i = 1:p.X, j = 1 + chol.update.beta$Pindex, x = 1)

        vi_beta_decomp <- solve(chol.update.beta$origL) %*% t(Pmatrix)
        vi_beta_mean <- chol.update.beta$mean
        log_det_beta_var <- -2 * sum(log(diag(chol.update.beta$origL)))

        vi_beta_mean <- Matrix(vi_beta_mean, dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(vi_alpha_mean, dimnames = list(fmt_names_Z, NULL))
      } else {
        stop("Invalid linpred method")
      }
    } else {
      stop("Invalid factorization method.")
    }
    
    if (family == 'linear'){
      adjust_var <- 1/sqrt(vi_sigmasq_a/vi_sigmasq_b)
      vi_beta_decomp <- vi_beta_decomp * adjust_var
      vi_alpha_decomp <- vi_alpha_decomp * adjust_var
      vi_joint_decomp <- vi_joint_decomp * adjust_var
      
      e_ln_sigmasq <- log(vi_sigmasq_b) - digamma(vi_sigmasq_a)
      
      log_det_joint_var <- log_det_joint_var + ncol(vi_joint_decomp) * e_ln_sigmasq
      log_det_beta_var <- log_det_beta_var + ncol(vi_beta_decomp) * e_ln_sigmasq
      log_det_alpha_var <- log_det_alpha_var + ncol(vi_alpha_decomp) * e_ln_sigmasq
    }

    if (debug_ELBO & it != 1) {
      variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
      vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
      
      debug_ELBO.2 <- calculate_ELBO(family = family,
        ELBO_type = ELBO_type,
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
        vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
        vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
        log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean,
        vi_r_sigma = vi_r_sigma,
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
      )
    }
    if (do_timing) {
      toc(quiet = verbose_time, log = T)
      tic("Update Sigma")
    }
    ###
    # Update \Sigma_j
    ##

    if (!do_huangwand){#Update standard Inverse-Wishart
      variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
      vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
      vi_sigma_alpha <- mapply(vi_sigma_outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN = function(i, j) {
        i + j
      })
    }else{
        #Update Inverse-Wishart
        # vi_sigma_outer_alpha <<- vi_sigma_outer_alpha
        # vi_a_a_jp <<- vi_a_a_jp
        # vi_a_b_jp <<- vi_a_b_jp
        variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
        vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
        
        for (inner_it in 1:10){
          
          vi_sigma_alpha <- mapply(vi_sigma_outer_alpha, vi_a_a_jp, 
           vi_a_b_jp, vi_a_nu_jp, SIMPLIFY = FALSE, 
           FUN = function(i, tilde.a, tilde.b, nu) {
             i + Diagonal(x = tilde.a/tilde.b) * 2 * nu
           })
          
          # vi_sigma_alpha <<- vi_sigma_alpha
          #Update a_{j,p}
          diag_Einv_sigma <- mapply(vi_sigma_alpha, 
                                    vi_sigma_alpha_nu, d_j, SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
                                      inv_phi <- solve(phi)
                                      sigma.inv <- nu * inv_phi
                                      return(diag(sigma.inv))
                                    })
          vi_a_b_jp <- mapply(vi_a_nu_jp, vi_a_APRIOR_jp, diag_Einv_sigma,
                              SIMPLIFY = FALSE,
                              FUN=function(nu, APRIOR, diag_j){
                                1/APRIOR^2 + nu * diag_j
                              })
          
        }
        # d_j <<- d_j
        # vi_alpha_decomp <<- vi_alpha_decomp
        # Tinv <<- Tinv
        # vi_alpha_mean <<- vi_alpha_mean
    }
    
    if (do_timing) {
      toc(quiet = verbose_time, log = T)
      tic("Update Aux")
    }
    # print(unlist(lapply(vi_sigma_alpha, as.vector)))
    # Update the auxilary parameters
    if (family == "negbin") {
      vi_r_param <- update_r(
        vi_r_mu = vi_r_mu, vi_r_sigma = vi_r_sigma,
        y = y, X = X, Z = Z, factorization_method = factorization_method,
        vi_beta_mean = vi_beta_mean, vi_beta_decomp = vi_beta_decomp,
        vi_alpha_mean = vi_alpha_mean, vi_alpha_decomp = vi_alpha_decomp,
        vi_joint_decomp = vi_joint_decomp, vi_r_method = vi_r_method
      )

      vi_r_mu <- vi_r_param[1]
      vi_r_sigma <- vi_r_param[2]
      vi_r_mean <- exp(vi_r_mu + vi_r_sigma / 2)
      
      s <- (y - vi_r_mean) / 2
      vi_pg_b <- y + vi_r_mean
    } else if (family == 'linear') {
      
      if (factorization_method == 'weak'){
        joint_quad <- cpp_zVz(Z = joint.XZ, V = as(vi_joint_decomp, "dgCMatrix"))
        vi_lp <- (s - as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean))^2 + joint_quad
      } else{
        beta_quad <- rowSums((X %*% t(vi_beta_decomp))^2)
        alpha_quad <- rowSums((Z %*% t(vi_alpha_decomp))^2)
        vi_lp <- (s - as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean))^2 + beta_quad + alpha_quad
      }
      
      vi_kernel <- expect_alpha_prior_kernel(vi_sigma_alpha = vi_sigma_alpha, 
          vi_sigma_alpha_nu = vi_sigma_alpha_nu, d_j = d_j,
          vi_sigma_outer_alpha = vi_sigma_outer_alpha)
      
      vi_sigmasq_b_trad <- (sum(vi_lp) + vi_kernel)/2 + vi_sigmasq_prior_b

      vi_sigmasq_b <- vi_sigmasq_b_trad
    }

    if (do_timing) {
      toc(quiet = verbose_time, log = T)
    }
    ### PARAMETER EXPANSIONS!
    if (debug_ELBO) {
      debug_ELBO.3 <- calculate_ELBO(family = family,
        ELBO_type = ELBO_type,
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
        log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma, choose_term = choose_term,
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
      )
    }
    if (do_timing) {
      tic("Update PX")
    }
    
    if (parameter_expansion %in% c("mean", "translation", "diagonal")) {
      
      # Do a simple mean adjusted expansion.
      # Get the mean of each random effect.
      vi_mu_j <- t(M_prime) %*% vi_alpha_mean
      
      # Remove the "excess mean" mu_j from each random effect \alpha_{j,g}
      # and add the summd mass back to the betas.
      vi_alpha_mean <- vi_alpha_mean - M_prime_one %*% vi_mu_j
      vi_beta_mean <- vi_beta_mean + t(M_mu_to_beta) %*% vi_mu_j

      variance_by_alpha_jg <- calculate_expected_outer_alpha(
        L = vi_alpha_decomp,
        alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions
      )
      vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
      if (parameter_expansion == "mean"){accept.PX <- TRUE}
    }  
    if (parameter_expansion %in% c("translation", "diagonal") & skip_translate == FALSE) {
      
      attempted_expansion <- attempted_expansion + 1
      
      prior.ELBO <- calculate_ELBO(family = family, ELBO_type = ELBO_type,
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
        log_det_joint_var = log_det_joint_var, 
        vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma,
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
      )

      Mmap <<- Mmap
      start_base_Z <<- start_base_Z
      mapping_new_Z <<- mapping_new_Z
      mapping_J <<- mapping_J
      zeromat_beta <<- zeromat_beta
      diag_vi_pg_mean <<- diag_vi_pg_mean
      p.X <<- ncol(X)
      number_of_RE <<- length(d_j)
      breaks_for_RE <<- breaks_for_RE
      
      cat('r')
      
      
      vi_alpha_decomp <<- vi_alpha_decomp
      diag_vi_pg_mean <<- diag_vi_pg_mean
      Z <<- Z
      store_design <<- store_design
      store_re_id <<- store_re_id
      store_id <<- store_id
      
      raw_R <- R_ridge <- vecR_ridge_new(L = vi_alpha_decomp, pg_mean = diag(diag_vi_pg_mean),
        mapping_J = mapping_J, d = d_j,
        store_id = store_id, store_re_id = store_re_id,
        store_design = store_design, 
        diag_only = (factorization_method == 'strong'))

      if (factorization_method == 'weak'){
        stop('no Translation PX for weak yet...')
      }
      cat('r')

      vi_alpha_mean <<- vi_alpha_mean
      mapping_new_Z <<- mapping_new_Z
      mapping_J <<- mapping_J
      Mmap <<- Mmap
      mapping_J <<- mapping_J
      start_base_Z <<- start_base_Z
      d_j <<- d_j
      R_design <- vecR_design(alpha_mu = as.vector(vi_alpha_mean), Z = mapping_new_Z, 
        M = Mmap, mapping_J = mapping_J, d = d_j,
        start_z = start_base_Z)

      vi_sigma_alpha <<- vi_sigma_alpha
      vi_sigma_alpha_nu <<- vi_sigma_alpha_nu
      
      moments_sigma_alpha <- mapply(vi_sigma_alpha, vi_sigma_alpha_nu, d_j, SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
        inv_phi <- solve(phi)
        
        sigma.inv <- nu * inv_phi
        
        ln.det <- log(det(phi)) - sum(digamma((nu - 1:d + 1) / 2)) - d * log(2)
        return(list(sigma.inv = sigma.inv, ln.det = ln.det))
      })
      
      diag_weight <- prior_sigma_alpha_phi
      prior_weight <- prior_sigma_alpha_nu
      
      if (do_huangwand){
        diag_weight <- mapply(vi_a_a_jp, vi_a_b_jp, vi_a_nu_jp, 
         SIMPLIFY = FALSE, 
         FUN = function(tilde.a, tilde.b, nu) {
           Diagonal(x = tilde.a/tilde.b) * 2 * nu
         })
        prior_weight <- vi_a_nu_jp + d_j - 1
      }
      
      moments_sigma_alpha <<- moments_sigma_alpha 
      diag_weight <<- diag_weight
      prior_weight <<- prior_weight
      
      vec_OSL_prior <- mapply(moments_sigma_alpha, diag_weight, prior_weight, 
        SIMPLIFY = FALSE, FUN=function(moment_j, phi_j, nu_j){
        as.vector(moment_j$sigma.inv %*% phi_j - nu_j * Diagonal(n = nrow(phi_j)))
      })
      vec_OSL_prior <- do.call('c', vec_OSL_prior)
      vec_OSL_prior <- matrix(c(rep(0, ncol(X)), vec_OSL_prior))
      
      R_design <<- R_design
      R_ridge <<- R_ridge
      X <<- X
      s <<- s
      

      #If a DIAGONAL expansion, then only update the diagonal elements
      if (parameter_expansion == "diagonal"){
        XR <- cbind(X, R_design[, diag_rho])
        R_ridge <- bdiag(zeromat_beta, R_ridge[diag_rho, diag_rho])
        if (do_huangwand){
          vec_OSL_prior <- do.call('c', mapply(vi_a_APRIOR_jp, vi_a_a_jp, vi_a_b_jp,
                                               SIMPLIFY = FALSE,
            FUN=function(i,a,b){1-2/i^2 * a/b}))
          vec_OSL_prior <- c(rep(0, p.X), vec_OSL_prior)
        }else{
          vec_OSL_prior <- vec_OSL_prior[c(1:p.X, p.X + diag_rho),,drop=F]
        }
        update_expansion_XR <- vecR_fast_ridge(X = drop0(XR), 
         omega = diag_vi_pg_mean, prior_precision = R_ridge, y = as.vector(s), 
         adjust_y = as.vector(vec_OSL_prior))
        update_expansion_bX <- Matrix(update_expansion_XR[1:p.X])
        update_expansion_R <- mapply(split(update_expansion_XR[-1:-p.X], 
          rep(1:number_of_RE, d_j)), d_j, SIMPLIFY = FALSE, 
          FUN=function(i,d){
            dg <- diag(x = d)
            diag(dg) <- i
            return(dg)
          })
         update_diag_R <- split(update_expansion_XR[-1:-p.X], rep(1:number_of_RE, d_j))
         rownames(update_expansion_bX) <- colnames(X)
      }else{#Do the FULL update
        XR <- cbind(X, R_design)
        R_ridge <- bdiag(zeromat_beta, R_ridge)
        update_expansion_XR <- vecR_fast_ridge(X = drop0(XR), 
         omega = diag_vi_pg_mean, prior_precision = R_ridge, y = as.vector(s), 
         adjust_y = as.vector(vec_OSL_prior))
        
        update_expansion_bX <- Matrix(update_expansion_XR[1:p.X])
        
        update_expansion_R <- mapply(split(update_expansion_XR[-1:-p.X], 
          rep(1:number_of_RE, d_j^2)), d_j, 
          SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d)})
        
      }
      
      est_rho <<- update_expansion_XR[-1:-p.X]
      print(round(est_rho, 4))
      if (parameter_expansion == 'diagonal'){
        if (max(abs(est_rho - 1)) < 1e-4){
          print('No further improvements')
          skip_translate <- TRUE
        }
      }else{
        if (max(abs(est_rho - stationary_rho)) < 1e-4){
          print('No further improvements')
          skip_translate <- TRUE
        }
      }
      

      prop_vi_sigma_alpha <- mapply(vi_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE,
                                    FUN=function(Phi, R){R %*% Phi %*% t(R)})
      # cat('r')
      mapping_for_R_block <- make_mapping_alpha(update_expansion_R, px.R = TRUE)
      update_expansion_Rblock <- prepare_T(mapping = mapping_for_R_block, levels_per_RE = g_j, num_REs = number_of_RE,
                        variables_per_RE = d_j, running_per_RE = breaks_for_RE, cyclical = FALSE, px.R = TRUE)
      # cat('r')
      stopifnot(all.equal(update_expansion_Rblock, bdiag(mapply(update_expansion_R, g_j, FUN=function(i,g){bdiag(lapply(1:g, FUN=function(k){i}))}))))

      update_expansion_R_logdet <- log(sapply(update_expansion_R, det))
      
      prop_vi_beta_mean <- update_expansion_bX
      prop_vi_alpha_mean <- update_expansion_Rblock %*% vi_alpha_mean
      cat('r')
      
      if (factorization_method != 'weak'){
        prop_log_det_joint_var <- prop_vi_joint_decomp <- NULL
        prop_vi_alpha_decomp <- vi_alpha_decomp %*% t(update_expansion_Rblock)
        prop_log_det_alpha_var <- log_det_alpha_var + 2 * sum(update_expansion_R_logdet * g_j)
        prop_log_det_beta_var <- log_det_beta_var
        prop_vi_beta_decomp <- vi_beta_decomp
        
        prop_variance_by_alpha_jg <- calculate_expected_outer_alpha(L = prop_vi_alpha_decomp, alpha_mu = as.vector(prop_vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
        prop_vi_sigma_outer_alpha <- prop_variance_by_alpha_jg$outer_alpha
      }else{
        stop('...')
      }
      
      if (do_huangwand){
        if (parameter_expansion == "diagonal"){
          prop_vi_a_b_jp <- mapply(vi_a_b_jp, update_diag_R, SIMPLIFY = FALSE,
                              FUN=function(i,j){i / j^2})
        }else{
          update_expansion_R <<- update_expansion_R
          
          prop_moments <- mapply(moments_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE,
             FUN=function(Phi, R){
               inv_R <- solve(R)
               return(diag(t(inv_R) %*% Phi$sigma.inv %*% inv_R))
              })
          prop_vi_a_b_jp <- mapply(vi_a_nu_jp, vi_a_APRIOR_jp, prop_moments,
             SIMPLIFY = FALSE,
             FUN=function(nu, APRIOR, diag_j){
               1/APRIOR^2 + nu * diag_j
             })
          
          # prop_diag_Einv_sigma <- mapply(prop_vi_sigma_alpha, 
          #   vi_sigma_alpha_nu, d_j, SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
          #     inv_phi <- solve(phi)
          #     sigma.inv <- nu * inv_phi
          #     return(diag(sigma.inv))
          #   })
          # prop_vi_a_b_jp <- mapply(vi_a_nu_jp, vi_a_APRIOR_jp, prop_diag_Einv_sigma,
          #   SIMPLIFY = FALSE,
          #   FUN=function(nu, APRIOR, diag_j){
          #     1/APRIOR^2 + nu * diag_j
          #   })
          
        }
      }else{
        prop_vi_a_b_jp <- NULL
      }
      # #L^T L = Variance
      # #R Var R^T --->
      # # L %*% R^T

      prop.ELBO <- calculate_ELBO(family = family,
        ELBO_type = ELBO_type,
        factorization_method = factorization_method,
        d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi,
        prior_sigma_alpha_nu = prior_sigma_alpha_nu,
        iw_prior_constant = iw_prior_constant,
        X = X, Z = Z, s = s, y = y,
        vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
        vi_sigma_alpha = prop_vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu,
        vi_sigma_outer_alpha = prop_vi_sigma_outer_alpha,
        vi_beta_mean = prop_vi_beta_mean, vi_alpha_mean = prop_vi_alpha_mean,
        
        log_det_beta_var = prop_log_det_beta_var, 
        log_det_alpha_var = prop_log_det_alpha_var,
        log_det_joint_var = prop_vi_log_det_joint_var,
        
        vi_beta_decomp = prop_vi_beta_decomp, 
        vi_alpha_decomp = prop_vi_alpha_decomp,
        vi_joint_decomp = prop_vi_joint_decomp,
        
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, 
        vi_a_b_jp = prop_vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
        choose_term
      )
      cat('d')
      if (prop.ELBO$ELBO > prior.ELBO$ELBO){
        #Accept the PX-VB adjustment (OSL).
        vi_beta_mean <- prop_vi_beta_mean
        vi_alpha_mean <- prop_vi_alpha_mean
        vi_sigma_alpha <- prop_vi_sigma_alpha
        if (factorization_method == 'weak'){
          stop('Setup reassignment weak')
        }else{
          vi_alpha_decomp <- prop_vi_alpha_decomp
          log_det_alpha_var <- prop_log_det_alpha_var
        }
        variance_by_alpha_jg <- prop_variance_by_alpha_jg
        vi_sigma_outer_alpha <- prop_vi_sigma_outer_alpha
        
        if (do_huangwand){
          vi_a_b_jp <- prop_vi_a_b_jp
        }
        accept.PX <- TRUE
      }else{
        accept.PX <- FALSE
      }
      print(accept.PX)
      accepted_times <- accept.PX + accepted_times

      rm(prop_vi_beta_mean, prop_vi_alpha_mean, prop_vi_sigma_alpha, prop_vi_alpha_decomp,
         prop_log_det_alpha_var, prop_variance_by_alpha_jg, prop_vi_sigma_outer_alpha)


      rownames(vi_alpha_mean) <- fmt_names_Z
    }
    
    if (parameter_expansion == "none") {
      accept.PX <- TRUE
    }

    if (do_timing) {
      toc(quiet = verbose_time, log = TRUE)
    }
    
    # Adjust the terms in the ELBO calculation that are different.

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

    if (accept.PX) {
      final.ELBO <- calculate_ELBO(family = family,
        ELBO_type = ELBO_type,
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
        log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma,
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
      )
    } else {
      # Should be no case
      if (accept.PX){
        final.ELBO <- prop.ELBO
      }else{
        final.ELBO <- prior.ELBO
      }
    }
    
    if (do_SQUAREM){
      squarem_list[[squarem_counter]] <- namedList(vi_sigma_alpha_nu, 
           vi_sigma_alpha, vi_alpha_mean, vi_beta_mean,
           vi_pg_c,
           vi_a_a_jp, vi_a_b_jp)
      if (squarem_counter %% 3 == 0){
        final.ELBO <<- final.ELBO
        squarem_list <<- squarem_list
        ELBOargs <<- list(family = family,
           ELBO_type = ELBO_type,
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
           log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma,
           do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
           vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
        )
        
        squarem_par <- c('vi_a_b_jp', 'vi_sigma_alpha', 'vi_pg_c',
                         'vi_alpha_mean', 'vi_beta_mean')
        # squarem_par <- c('vi_alpha_mean', 'vi_beta_mean')
        
        prep_SQUAREM <- mapply(squarem_par, SIMPLIFY = FALSE, FUN=function(i){
          if (class(squarem_list[[1]][[i]]) == 'list'){
            r <- mapply(squarem_list[[2]][[i]], squarem_list[[1]][[i]], FUN=function(i,j){i - j})
            d2 <- mapply(squarem_list[[3]][[i]], squarem_list[[2]][[i]], FUN=function(i,j){i - j})
            v <- mapply(d2, r, FUN=function(i,j){i - j})
            norm_sq_r <- sum(unlist(lapply(r, as.vector))^2)
            norm_sq_v <- sum(unlist(lapply(v, as.vector))^2)
          }else{
            r <- squarem_list[[2]][[i]] - squarem_list[[1]][[i]]
            d2 <- squarem_list[[3]][[i]] - squarem_list[[2]][[i]]
            v <- d2 - r
            norm_sq_r <- sum(r^2)
            norm_sq_v <- sum(v^2)
          }
          return(list(first = squarem_list[[1]][[i]], 
                      second = squarem_list[[2]][[i]], r = r, v = v, norm_sq_r = norm_sq_r, norm_sq_v = norm_sq_v))
        })
        
        alpha <- -sqrt(sum(sapply(prep_SQUAREM, FUN=function(i){i$norm_sq_r}))) /
          sqrt(sum(sapply(prep_SQUAREM, FUN=function(i){i$norm_sq_v})))

        alpha <- alpha / 2
        
        if (alpha > -1){
          alpha <- -1.01
        }

        prop_squarem <- lapply(prep_SQUAREM, FUN=function(i){
          if (class(i$first) == 'list'){
            prop_squarem <- mapply(i$first, i$second, 
                                   i$r, i$v, SIMPLIFY = FALSE, FUN=function(i_1, s_1, r_1, v_1){
              i_1 - 2 * alpha * r_1 + alpha^2 * v_1
              # s_1 - alpha * (v_1 + r_1)
            })
            names(prop_squarem) <- names(i$first)
          }else{
            # prop_squarem <- i$second - alpha * (i$r + i$v)
            prop_squarem <- i$first - 2 * alpha * i$r + alpha^2 * i$v
          }
          return(prop_squarem)
        })
        
        names(prop_squarem) <- squarem_par

        prop_ELBOargs <- ELBOargs
        
        if ('vi_a_b_jp' %in% squarem_par){
          prop_squarem$vi_a_b_jp <- mapply(prop_squarem$vi_a_b_jp, squarem_list[[3]]$vi_a_b_jp, FUN=function(i,j){
            i[i < 0] <- j[i < 0]
            return(i)
          })
        }
        if ('vi_pg_c' %in% squarem_par){
          prop_squarem$vi_pg_c[prop_squarem$vi_pg_c < 0] <- squarem_list[[3]]$vi_pg_c[prop_squarem$vi_pc_c < 0]
        }
        if ('vi_sigma_alpha' %in% squarem_par){
          prop_squarem$vi_sigma_alpha <- mapply(prop_squarem$vi_sigma_alpha, squarem_list[[3]]$vi_sigma_alpha, FUN=function(i,j){
            if (det(i) < 0){
              i <- j 
            }
            return(i)
          })
        }
        
        for (v in names(prop_squarem)){
          prop_ELBOargs[[v]] <- prop_squarem[[v]]
        }
        if ('vi_alpha_mean' %in% squarem_par){
          prop_variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, 
              alpha_mu = as.vector(prop_squarem$vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
          prop_ELBOargs[['vi_sigma_outer_alpha']] <- prop_variance_by_alpha_jg$outer_alpha
          squarem_par <- c(squarem_par, 'vi_sigma_outer_alpha')
        }
        if ('vi_pg_c' %in% squarem_par){
          if (family != 'binomial'){stop('check squarem for non-binomial case')}
          prop_vi_pg_mean <- vi_pg_b / (2 * prop_ELBOargs$vi_pg_c) * tanh(prop_ELBOargs$vi_pg_c / 2)
          prop_ELBOargs[['vi_pg_mean']] <- prop_vi_pg_mean
          squarem_par <- c(squarem_par, 'vi_pg_c')
        }
        
        
        elbo_init <- do.call("calculate_ELBO", ELBOargs)$ELBO
        elbo_squarem <- do.call("calculate_ELBO", prop_ELBOargs)$ELBO
        print(c(elbo_squarem, elbo_init))
        
        if (elbo_squarem > elbo_init){
          cat('SUCCESS')
          for (v in squarem_par){
            assign(v, prop_ELBOargs[[v]])
          }
        }else{
          cat('FAIL')
        }
        squarem_list <- list()
        squarem_counter <- 1
      }else{
        
        squarem_counter <- squarem_counter + 1
        
      }
      
    }
    if (do_timing) {
      toc(quiet = verbose_time, log = T)
      tic("Final Cleanup")
    }

    if (debug_ELBO & it != 1) {
      debug_ELBO.1$step <- 1
      debug_ELBO.2$step <- 2
      debug_ELBO.3$step <- 3
      final.ELBO$step <- 4
      update_ELBO <- bind_rows(debug_ELBO.1, debug_ELBO.2, debug_ELBO.3, final.ELBO)
      update_ELBO$it <- it
      store_ELBO <- bind_rows(store_ELBO, update_ELBO)
    } else {
      final.ELBO$it <- it
      store_ELBO <- bind_rows(store_ELBO, final.ELBO)
    }

    ## Change diagnostics
    change_elbo <- final.ELBO - lagged_ELBO

    change_alpha_mean <- max(abs(vi_alpha_mean - lagged_alpha_mean))
    change_beta_mean <- max(abs(vi_beta_mean - lagged_beta_mean))
    
    unlist_vi <- c(unlist(lapply(vi_sigma_alpha, as.vector)), unlist(vi_a_b_jp))
    if (debug_ELBO){
      store_vi <- rbind(store_vi, data.frame(t(as.vector(unlist_vi))) %>% mutate(it = it))
      print(unlist_vi)
    }
    
    change_sigma_mean <- mapply(vi_sigma_alpha, lagged_sigma_alpha, FUN = function(i, j) {
      max(abs(i - j))
    })

    if (factorization_method == "weak") {
      change_joint_var <- 0 # change_joint_var <- max(abs(vi_joint_decomp - lagged_joint_decomp))
      change_alpha_var <- change_beta_var <- 0
    } else {
      change_joint_var <- 0
      change_alpha_var <- max(abs(vi_alpha_decomp - lagged_alpha_decomp))
      change_beta_var <- max(abs(vi_beta_decomp - lagged_beta_decomp))
    }

    change_vi_r_mu <- vi_r_mu - lagged_vi_r_mu

    if (do_timing) {
      toc(quiet = verbose_time, log = T)
    }
    if (debug_param) {
      store_beta[it, ] <- as.vector(vi_beta_mean)
    }
    change_all <- data.frame(change_alpha_mean, change_beta_mean, 
        t(change_sigma_mean), change_alpha_var, change_beta_var, change_joint_var, change_vi_r_mu)
    if ((max(change_all) < tolerance_parameters) | (change_elbo$ELBO > 0 & change_elbo$ELBO < tolerance_elbo)) {
      if (!quiet) {
        message(paste0("Converged after ", it, " iterations with ELBO change of ", round(change_elbo[1], 1 + abs(floor(log(tolerance_elbo) / log(10))))))
        message(paste0("The largest change in any variational parameter was ", round(max(change_all), 1 + abs(floor(log(tolerance_parameters) / log(10))))))
      }
      break
    }
    if (debug_ELBO){
      change_all$it <- it
      store_parameter_traj <- rbind(store_parameter_traj, change_all)
    }
    if (!quiet & (it %% print_prog == 0)) {
      message(paste0("ELBO Change: ", round(change_elbo$ELBO, 10)))
      message(paste0("Other Parameter Changes: ", max(change_all)))
    }

    lagged_alpha_mean <- vi_alpha_mean
    lagged_beta_mean <- vi_beta_mean
    lagged_alpha_decomp <- vi_alpha_decomp
    lagged_beta_decomp <- vi_beta_decomp
    lagged_sigma_alpha <- vi_sigma_alpha
    lagged_vi_r_mu <- vi_r_mu
    lagged_ELBO <- final.ELBO
  }
  if (it == iterations) {
    message(paste0("Ended without Convergence after", it, " iterations : ELBO change of ", round(change_elbo[1], abs(floor(log(tolerance_elbo) / log(10))))))
  }

  if (debug_ELBO) {
    d.ELBO <- with(store_ELBO, ELBO - dplyr::lag(ELBO))
    sum.ELBO <- store_ELBO
    sum.ELBO$diff <- d.ELBO
    sum.ELBO <- summarize(group_by_at(sum.ELBO, .vars = "step"),
      negative = mean(.data$diff < 0, na.rm = T)
    )
  } else {
    store_parameter_traj <- NULL
    d.ELBO <- NULL
  }
  if (parameter_expansion == "translation") {
    final.ELBO$accepted_PX <- accepted_times / attempted_expansion
  }
  output <- list(
    beta = list(mean = vi_beta_mean),
    ELBO = final.ELBO, debug_ELBO = d.ELBO,
    ELBO_trajectory = store_ELBO,
    parameter.change = change_all,
    parameter.vi = store_vi,
    parameter.path = store_parameter_traj,
    sigma = list(cov = vi_sigma_alpha, df = vi_sigma_alpha_nu),
    alpha = list(mean = vi_alpha_mean)
  )
  if (family == 'linear'){
    output$sigmasq <- list(a = vi_sigmasq_a, b = vi_sigmasq_b)
  }else if (family == 'negbin'){
    
  }
  output$family <- family
  output$control <- control

  if (do_timing) {
    tic_log <- tictoc::tic.log(format = FALSE)
    tic_log <- data.frame(stage = sapply(tic_log, FUN = function(i) {
      i$msg
    }), time = sapply(tic_log, FUN = function(i) {
      i$toc - i$tic
    }), stringsAsFactors = F)

    tic.clear()
    tic.clearlog()

    tic_summary <- summarize(
      .data = group_by_at(.vars = "stage", tic_log),
      .groups = "keep",
      n = dplyr::n(), mean = mean(.data$time),
      min = min(.data$time), max = max(.data$time),
      total = sum(.data$time)
    )
  } else {
    tic_summary <- NULL
  }
  if (debug_param) {
    output$parameter_trajectory <- list(beta = store_beta)
  }
  if (factorization_method == "weak") {
    output$joint <- vi_joint_decomp
  }
  if (control$return_data) {
    output$data <- list(X = X, Z = Z, y = y, trials = trials)
  }
  output$formula <- formula
  output$alpha$dia.var <- unlist(lapply(variance_by_alpha_jg$variance_jg, FUN = function(i) {
    as.vector(sapply(i, diag))
  }))
  output$beta$var <- t(vi_beta_decomp) %*% vi_beta_decomp
  output$beta$decomp_var <- vi_beta_decomp

  if (family == "negbin") {
    output$r <- list(mu = vi_r_mu, sigma = vi_r_sigma, method = vi_r_method)
  }
  if (do_huangwand){
    output$hw <- list(a = vi_a_a_jp, b = vi_a_b_jp)
  }

  output$internal_parameters <- list(
    it_used = it, it_max = iterations,
    missing_obs = missing_obs, N = nrow(X),
    names_of_RE = names_of_RE, d_j = d_j, g_j = g_j
  )

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
  class(output) <- "vglmer"
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
#'   \item kn: See Kass and Natarajan (2006)
#'   \item jeffreys: IW(0, 0)
#'   \item mcmcglmm: IW(0, I)
#'   \item mvD: IW(-d, I)
#'   \item mean_exists: IW(d + 1, I)
#'   \item limit: IW(d - 1, 0)
#'   \item uniform: IW(-[d+1], 0)
#'   }
#'   The model may fail to converge if an improper prior is used.
#' @param tolerance_elbo Change in ELBO to stop algorithm.
#' @param tolerance_parameters Change in value of any parameter to stop algorithm.
#'
#' @param parameter_expansion At moment, accepts 'mean' or 'none'. 'mean' is
#'   costless and should always be used!
#' @param prevent_degeneracy Ignored for the moment.
#' @param force_whole Require whole numbers. Set to FALSE to allow "quasi-binomial".
#'
#' @param vi_r_method Type of estimate for "r"; at moment, "fixed" (provide r),
#'   "VEM" (treat r as point estimate; default);
#'   "Laplace" (estimate using Laplace approximation described in the Addendum on GitHub); or "delta"
#'   (experimential).
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
                           prior_variance = "mean_exists", factorization_method = "weak",
                           tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
                           prevent_degeneracy = FALSE, force_whole = TRUE, verbose_time = TRUE,
                           parameter_expansion = "mean", random_seed = 1, do_timing = FALSE,
                           debug_param = FALSE, return_data = FALSE, linpred_method = "joint",
                           vi_r_method = "VEM", vi_r_val = NA,
                           debug_ELBO = FALSE, print_prog = NULL, quiet = T, init = "EM") {
  # use checkmate package to verify arguments
  assert(
    check_integerish(iterations, lower = 1),
    check_double(c(tolerance_elbo, tolerance_parameters), len = 2, lower = 0),
    check_logical(c(
      prevent_degeneracy, force_whole, verbose_time, do_timing,
      debug_param, return_data, debug_ELBO, quiet
    ), len = 8),
    check_choice(factorization_method, c("weak", "strong", "partial")),
    check_choice(prior_variance, c("kn", "hw", "mean_exists", "jeffreys", "mcmcglmm", "mvD", "limit", "uniform")),
    check_choice(linpred_method, c("joint", "cyclical", "solve_normal")),
    check_choice(vi_r_method, c("VEM", "fixed", "Laplace", "delta")),
    check_double(vi_r_val, all.missing = TRUE),
    check_int(print_prog, null.ok = TRUE),
    check_choice(init, c("EM", "random", "zero")),
    check_double(random_seed),
    combine = "and"
  )

  if (vi_r_method == "fixed" & is.na(vi_r_val)) {
    stop('vi_r_val must not be NA if vi_r_method = "fixed"')
  }

  output <- namedList(
    iterations, prior_variance, factorization_method,
    tolerance_elbo, tolerance_parameters,
    prevent_degeneracy, force_whole, verbose_time,
    parameter_expansion, random_seed, do_timing, debug_param, return_data,
    linpred_method, vi_r_method, vi_r_val, debug_ELBO, print_prog, quiet, init
  )

  class(output) <- c("vglmer_control")
  return(output)
}

# Simple function to create named list
# https://stackoverflow.com/questions/16951080/can-lists-be-created-that-name-themselves-based-on-input-object-names
# Identical to version used in lme4:::namedList, see also loo::nlist
#' @importFrom stats setNames
namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm == "")) nm[nonames] <- snm[nonames]
  setNames(L, nm)
}
