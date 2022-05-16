#' Variational Inference for Non-Linear Hierarchical Models
#'
#' Estimate a hierarchical model using mean-field variational
#' inference. Accepts standard syntax for lme4: \code{y ~ X + (1 + Z | g)}.
#' Options are described below.
#'
#' @param formula lme4-style formula that can include multiple random effects.
#'   Splines can be specified via `s_v(x)`, see the examples for details.
#' @param data data.frame containing the outcome and variables.
#' @param family Options are "binomial" or "negbin". If "binomial", outcome must
#'   be either {0,1} (binary) or cbind(success, failure). Non-integer values are
#'   permitted for binomial if "force_whole" is set to FALSE in vglmer_control.
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
#' @importFrom lme4 mkReTrms findbars subbars
#' @importFrom stats model.response model.matrix model.frame rnorm rWishart
#'   qlogis optim residuals lm plogis
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
  
  parse_formula <- vglmer_interpret.gam0(subbars(formula), extra.special = 'v_s')

  if (any(!sapply(parse_formula$smooth.spec, inherits, what = 'vglmer_spline'))){
    stop('gam specials are not permitted; use v_s(...) and see documentation.')
  }
  
  if (control$verify_columns){
    if (!all(parse_formula$pred.names %in% colnames(data))){
      missing_columns <- setdiff(parse_formula$pred.names, colnames(data))
      stop(
        paste0('The following columns are missing from "data". Can override with vglmer_control (not usually desirable): ', 
               paste(missing_columns, collapse =', '))
      )
    }
  }
  
  data <- model.frame(parse_formula$fake.formula, data)
  
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
  fe_fmla <- vglmer_interpret.gam0(nobars(formula), extra.special = 'v_s')
  
  if (length(fe_fmla$smooth.spec) > 0){
    # Add the linear spline terms to the main effect.
    
    fe_fmla <- update.formula(fe_fmla$pf,
        paste0('. ~ . + 1 + ',
        paste0(sapply(fe_fmla$smooth.spec, FUN=function(i){i$term}),
        collapse = " + "))
    )
    
    
  }else{
    fe_fmla <- fe_fmla$pf
  }
  
  X <- model.matrix(fe_fmla, data = data)
  
  # Extract the Z (Random Effect) design matrix.
  re_fmla <- findbars(formula)
  
  # If using splines by group, add random effects to
  # the main level.
  if (!all(sapply(parse_formula$smooth.spec, 
      FUN=function(i){i$by}) %in% c('NA'))){
    
    by_splines <- parse_formula$smooth.spec[
      which(sapply(parse_formula$smooth.spec, FUN=function(i){i$by != "NA"}))
    ]
    
    character_re <- lapply(re_fmla, FUN=function(i){strsplit(deparse(i), split = ' \\| ')[[1]]})
    character_re_group <- sapply(character_re, FUN=function(i){i[2]})
    
    if (any(duplicated(character_re_group))){
      stop('Some random groups are duplicated. Reformulate initial formula.')
    }
    
    for (v in sapply(character_re, FUN=function(i){i[2]})){
      if (!(is.factor(data[[v]]) | is.character(data[[v]]))){
        data[[v]] <- as.character(data[[v]])
      }
    } 
    
    for (b in by_splines){
      
      b_term <- b$term
      b_by <- b$by
      
      if (!(is.factor(data[[b_by]]) | is.character(data[[b_by]]))){
       stop('For now, all v_s spline "by" factors must be characters or factors.') 
      }
      
      # If "by" grouping already used, then add to the RE
      if (b_by %in% character_re_group){
        
        position_b_by <- which(b_by == character_re_group)
        existing_re_b_by <- character_re[[position_b_by]][1]
        new_re_b_by <- paste0(unique(c('1', strsplit(existing_re_b_by, split=' \\+ ')[[1]], b_term)), collapse = ' + ')
        character_re[[position_b_by]][1] <- new_re_b_by
      }else{
        # If not, then add a new RE group with a 
        # random intercept and random slope.
        character_re <- c(character_re, list(c(paste0('1 + ', b_term), b_by)))
        character_re_group <- sapply(character_re, FUN=function(i){i[2]})
      }
    }
    
    character_re_fmla <- paste(sapply(character_re, FUN=function(i){paste0('(', i[1], ' | ', i[2], ' )')}), collapse = " + ")
    
    old_re <- re_fmla
    re_fmla <- lapply(character_re, FUN=function(i){str2lang(paste0(i[1], ' | ', i[2]))})
    
  }
  
  if (!is.null(re_fmla)){
    mk_Z <- mkReTrms(re_fmla, data, reorder.terms = FALSE, reorder.vars = FALSE)
    Z <- t(mk_Z$Zt)
    
    p.X <- ncol(X)
    p.Z <- ncol(Z)
    
    ####
    # Process the REs to get various useful terms.
    ####
    # RE names and names of variables included for each.
    names_of_RE <- mk_Z$cnms
    
    number_of_RE <- length(mk_Z$Gp) - 1
    
    if ( (number_of_RE < 1) & (length(parse_formula$smooth.spec) == 0) ) {
      stop("Need to provide at least one random effect or spline...")
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
  }else{
    
    Z <- drop0(Matrix(nrow = nrow(X), ncol = 0))
    p.X <- ncol(X)
    p.Z <- 0
    names_of_RE <- c()
    number_of_RE <- 0
    breaks_for_RE <- c(0)
    d_j <- c()
    g_j <- c()
    fmt_names_Z <- c()
    cyclical_pos <- list()
    
    if ( (length(parse_formula$smooth.spec) == 0) ) {
      stop("Need to provide at least one random effect or spline...")
    }

  }
  
  M.names <- cbind(unlist(mapply(names_of_RE, g_j, FUN = function(i, j) {
    rep(i, j)
  })))
  M <- cbind(match(M.names[, 1], colnames(X)), rep(1 / g_j, d_j * g_j))
  

  if (nrow(M) > 0){
    M <- sparseMatrix(i = 1:p.Z, j = M[, 1], x = M[, 2], dims = c(p.Z, ncol(X)))
  }else{
    M <- drop0(matrix(0, nrow = 0, ncol = ncol(X)))
  }
  
  if (!is.null(names_of_RE)){
    any_Mprime <- TRUE
    M_prime.names <- paste0(rep(names(names_of_RE), g_j * d_j), " @ ", M.names)
    M_prime <- cbind(match(M_prime.names, unique(M_prime.names)), rep(1 / g_j, d_j * g_j))
    M_prime <- sparseMatrix(i = 1:p.Z, j = M_prime[, 1], x = M_prime[, 2])
    colnames(M_prime) <- unique(M_prime.names)
    
    M_prime_one <- M_prime
    M_prime_one@x <- rep(1, length(M_prime_one@x))
    
    stopifnot(identical(paste0(rep(names(names_of_RE), d_j), " @ ", unlist(names_of_RE)), colnames(M_prime)))
    
    M_mu_to_beta <- sparseMatrix(i = 1:sum(d_j), j = match(unlist(names_of_RE), colnames(X)), x = 1, dims = c(sum(d_j), p.X))
    
    
  }else{
    any_Mprime <- FALSE
    M_prime_one <- M_prime <- drop0(matrix(0, nrow = 0, ncol = 0))
    M_mu_to_beta <- drop0(matrix(0, nrow = 0, ncol = p.X))
  }
  
  colnames(M_mu_to_beta) <- colnames(X)
  rownames(M_mu_to_beta) <- colnames(M_prime)
  
  
  # Extract the Specials
  

  if (length(parse_formula$smooth.spec) > 0){
    
    base_specials <- length(parse_formula$smooth.spec)
    
    # Number of splines + one for each "by"...
    
    n.specials <- base_specials +
      sum(sapply(parse_formula$smooth.spec, FUN=function(i){i$by}) != "NA")
    
    Z.spline.attr <- as.list(rep(NA, base_specials))

    Z.spline <- as.list(rep(NA, n.specials))
    Z.spline.size <- rep(NA, n.specials)
    special_counter <- 1
    
    
    for (i in 1:base_specials){
      
      special_i <- parse_formula$smooth.spec[[i]]

      all_splines_i <- vglmer_build_spline(x = data[[special_i$term]], 
          by = data[[special_i$by]],
          knots = special_i$knots, type = special_i$type,
          outer_okay = special_i$outer_okay)
      
      Z.spline.attr[[i]] <- c(all_splines_i[[1]]$attr, 
        list(type = special_i$type, by = special_i$by))
      
      spline_counter <- 1
      for (spline_i in all_splines_i){
        
        stopifnot(spline_counter %in% 1:2)
        
        colnames(spline_i$x) <- paste0('spline @ ', special_i$term, ' @ ', colnames(spline_i$x))
        
        if (spline_counter > 1){
          spline_name <- paste0('spline-',special_i$term,'-', i, '-int')
        }else{
          spline_name <- paste0('spline-', special_i$term, '-', i, '-base')
        }
        
        Z.spline[[special_counter]] <- spline_i$x
        Z.spline.size[special_counter] <- ncol(spline_i$x)
        
        names_of_RE[[spline_name]] <- spline_name
        number_of_RE <- number_of_RE + 1
        d_j <- setNames(c(d_j, 1), c(names(d_j), spline_name))
        g_j <- setNames(c(g_j, ncol(spline_i$x)), c(names(g_j), spline_name))
        breaks_for_RE <- c(breaks_for_RE, max(breaks_for_RE) + ncol(spline_i$x))
        fmt_names_Z <- c(fmt_names_Z, colnames(spline_i$x))
        p.Z <- p.Z + ncol(spline_i$x)
        spline_counter <- spline_counter + 1
        special_counter <- special_counter + 1
      }
    }
    
    cyclical_pos <- lapply(1:number_of_RE, FUN = function(i) {
      seq(breaks_for_RE[i] + 1, breaks_for_RE[i + 1])
    })
    
    Z.spline <- drop0(do.call('cbind', Z.spline))
    
    Z <- drop0(cbind(Z, Z.spline))
    
    if (ncol(M_prime) == 0){
      M_prime <- rbind(M_prime, 
                       drop0(matrix(0, nrow = ncol(Z.spline), ncol = 0)))
      M_prime_one <- rbind(M_prime_one, 
                           drop0(matrix(0, nrow = ncol(Z.spline), ncol = 0)))
    }else{
      M_prime <- rbind(M_prime, 
         drop0(sparseMatrix(i = 1, j = 1, x = 0, 
         dims = c(ncol(Z.spline), ncol(M_prime)))))
      M_prime_one <- rbind(M_prime_one, 
         drop0(sparseMatrix(i = 1, j = 1, x = 0, 
         dims = c(ncol(Z.spline), ncol(M_prime_one)))))
    }
    
  }else{
    n.specials <- 0
    Z.spline.attr <- NULL
    Z.spline <- NULL
    Z.spline.size <- NULL
  }

  if (factorization_method == 'collapsed' | control$parameter_expansion %in% c('translation', 'diagonal')){
    parsed_RE_groups <- get_RE_groups(formula = re_fmla, data = data)
    store_design_Z <- parsed_RE_groups$design
    store_design_Z <- lapply(store_design_Z, as.matrix)
    store_assignment_Z <- parsed_RE_groups$factor

    rm(parsed_RE_groups); gc()
  }
  

  
  # List of Lists
  # Outer list: one for RE
  # Inner List: One for each GROUP with its row positions.
  outer_alpha_RE_positions <- mapply(d_j, g_j, breaks_for_RE[-length(breaks_for_RE)], 
    SIMPLIFY = FALSE, FUN = function(a, b, m) {
      split(m + seq(1, a * b), rep(1:b, each = a))
  })

  if (anyDuplicated(unlist(outer_alpha_RE_positions)) != 0 | max(unlist(outer_alpha_RE_positions)) != p.Z) {
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
    
  } else if (family == 'negbin') {
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

    choose_term <- -sum(lgamma(y + 1)) - sum(y) * log(2)
  }else{
    stop('family must be linear, binomial, or negative binomial.')
  }

  # Initalize variational parameters.
  # Note that we keep a sparse matrix or lowertri such that
  # t(vi_beta_decomp) %*% vi_beta_decomp = VARIANCE

  vi_beta_decomp <- Diagonal(x = rep(0, ncol(X)))
  vi_alpha_decomp <- Diagonal(x = rep(0, p.Z))

  vi_sigma_alpha_nu <- g_j

  prior_variance <- control$prior_variance
  do_huangwand <- FALSE
  vi_a_APRIOR_jp <- vi_a_nu_jp <- vi_a_a_jp <- vi_a_b_jp <- NULL
  prior_sigma_alpha_nu <- prior_sigma_alpha_phi <- NULL
  
  if (prior_variance == 'hw') {
    
    if (family == 'linear'){stop('Set up Huang-Wand for non-binomial.')}
    
    do_huangwand <- TRUE
    vi_a_nu_jp <- rep(2, length(d_j))
    names(vi_a_nu_jp) <- names(names_of_RE)
    vi_a_APRIOR_jp <- lapply(d_j, FUN=function(i){rep(5, i)})
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

  if ( control$init %in% c("EM", "EM_FE") ){
    if (family == "linear"){
      jointXZ <- cbind(X,Z)
      if (control$init == 'EM_FE'){
        EM_init <- LinRegChol(X = X,
           omega = sparseMatrix(i = 1:nrow(X), j = 1:nrow(X), x = 1),
           y = y, prior_precision = sparseMatrix(i = 1:p.X, j = 1:p.X, x = 1e-5))$mean
        # stop('Setup EM init for linear')
        # solve(Matrix::Cholesky(  t(joint.XZ) %*% sparseMatrix(i = 1:N, j = 1:N, x = pg_mean) %*% joint.XZ + EM_variance),
        #       t(joint.XZ) %*% (adj_out) )
        EM_init <- list('beta' = EM_init, 'alpha' = rep(0, p.Z))
      }else{
        stop('Setup EM init')
        
        EM_init <- LinRegChol(X = jointXZ, 
                              omega = sparseMatrix(i = 1:nrow(jointXZ), j = 1:nrow(jointXZ), x = 1), 
                              y = y, prior_precision = sparseMatrix(i = 1:ncol(jointXZ), j = 1:ncol(jointXZ), x = 1/4))$mean
        EM_init <- list('beta' = EM_init[1:p.X], 'alpha' = EM_init[-1:-p.X])
      }
      rm(jointXZ)
    } else if (family == "negbin") {
      if (control$init == 'EM_FE'){
        EM_init <- EM_prelim_nb(X = X, Z = drop0(matrix(0, nrow = nrow(X), ncol = 0)), y = y, est_r = exp(vi_r_mu), iter = 15, ridge = 10^5)
        EM_init <- list('beta' = EM_init$beta, 'alpha' = rep(0, p.Z))
      }else{
        EM_init <- EM_prelim_nb(X = X, Z = Z, y = y, est_r = exp(vi_r_mu), iter = 15, ridge = 4)
      }
    } else {
      if (control$init == 'EM_FE'){
        EM_init <- EM_prelim_logit(X = X, Z = drop0(matrix(0, nrow = nrow(X), ncol = 0)), s = s, pg_b = vi_pg_b, iter = 15, ridge = 10^5)
        EM_init <- list('beta' = EM_init$beta, 'alpha' = rep(0, p.Z))
      }else{
        EM_init <- EM_prelim_logit(X = X, Z = Z, s = s, pg_b = vi_pg_b, iter = 15, ridge = 4)
      }
    }

    vi_beta_mean <- matrix(EM_init$beta)
    vi_alpha_mean <- matrix(EM_init$alpha)

    vi_sigma_alpha <- calculate_expected_outer_alpha(
      alpha_mu = vi_alpha_mean,
      alpha_decomp_var = sparseMatrix(i = 1, j = 1, x = 1e-4, dims = rep(p.Z, 2)),
      factorization_method = 'init',
      re_position_list = outer_alpha_RE_positions, do_adjustment = FALSE,
      tP = matrix(0, 1, 1), L_beta = matrix(0, 1, 1)
    )
    if (!do_huangwand){
      vi_sigma_alpha <- mapply(vi_sigma_alpha$outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN = function(i, j) {
        i + j
      })
    }else{
      #Update Inverse-Wishart
      # vi_a_b_jp <<- vi_a_b_jp
      # vi_a_nu_jp <<- vi_a_nu_jp
      # vi_a_a_jp <<- vi_a_a_jp
      # vi_sigma_alpha <<- vi_sigma_alpha
      # vi_a_APRIOR_jp <<- vi_a_APRIOR_jp
      
      vi_sigma_alpha <- mapply(vi_sigma_alpha$outer_alpha, vi_a_a_jp, 
       vi_a_b_jp, vi_a_nu_jp, SIMPLIFY = FALSE, FUN = function(i, tilde.a, tilde.b, nu) {
         i + sparseMatrix(i = seq_len(nrow(i)), j = seq_len(nrow(i)), x = 1)
       })
      
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
  } else if (control$init == "random") {
    set.seed(control$random_seed)
    vi_beta_mean <- rnorm(p.X)
    vi_alpha_mean <- rep(0, p.Z)

    vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN = function(d, g) {
      
      out <- rWishart(n = 1, df = ifelse(g >= d, g, d), Sigma = diag(d))[ , , 1]
      
      if (d == 1){
        out <- matrix(out)
      }
      
      return(out)
      
    })

  } else if (control$init == "zero") {
    vi_beta_mean <- rep(0, p.X)

    if (family == "binomial") {
      vi_beta_mean[1] <- qlogis(sum(y) / sum(trials))
    } else if (family == "negbin") {
      vi_beta_mean[1] <- log(mean(y))
    } else if (family == 'linear'){
      vi_beta_mean[1] <- mean(y)
    } else {
      stop('Set up init')
    }

    vi_alpha_mean <- rep(0, p.Z)

    vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN = function(d, g) {
      diag(x = 1, ncol = d, nrow = d)
    })
    # if (do_huangwand){stop('Setup init for zero')}
  } else {
    stop("Provide init = EM or random")
  }

  zero_mat <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(p.X, p.X))
  zero_mat <- drop0(zero_mat)

  if (factorization_method %in% c("weak")) {
    vi_joint_decomp <- bdiag(vi_beta_decomp, vi_alpha_decomp)
    joint.XZ <- cbind(X, Z)
    log_det_beta_var <- log_det_alpha_var <- NULL
  } else {
    vi_joint_decomp <- NULL
    log_det_joint_var <- NULL
  }

  if (family == 'linear'){
    exclude_collapsed <- seq_len(nrow(Z))
  }else{
    exclude_collapsed <- which(vi_pg_b != 0)
  }
  
  # if (!quiet){warning('Check vi sigma alpha nu')}

  # Create mapping for this to allow sparse implementations.

  mapping_sigma_alpha <- make_mapping_alpha(vi_sigma_alpha)

  running_log_det_alpha_var <- rep(NA, number_of_RE)

  if (factorization_method != "collapsed"){
    lagged_alpha_mean <- rep(-Inf, p.Z)
    lagged_beta_mean <- rep(-Inf, p.X)
    lagged_sigma_alpha <- vi_sigma_alpha
    if (factorization_method %in% c("weak")) {
      lagged_joint_decomp <- vi_joint_decomp
    } else {
      lagged_alpha_decomp <- vi_alpha_decomp
      lagged_beta_decomp <- vi_beta_decomp
    }
  }
  lagged_vi_r_mu <- -Inf
  lagged_vi_sigmasq_a <- lagged_vi_sigmasq_b <- -Inf
  lagged_ELBO <- -Inf
  accepted_times <- NA

  skip_translate <- FALSE
  if (factorization_method == "collapsed" | (parameter_expansion %in%  c("translation", "diagonal") & any_Mprime)) {
    
    if (do_timing){
      tic('Build PX R Terms')
    }
    
    accepted_times <- 0
    attempted_expansion <- 0
    
    spline_REs <- grepl(names(d_j), pattern='^spline-')
    
    nonspline_positions <- sort(unlist(outer_alpha_RE_positions[!spline_REs]))
    
    size_splines <- sum((d_j * g_j)[spline_REs])
    
    stationary_rho <- do.call('c', lapply(d_j[!spline_REs], FUN=function(i){as.vector(diag(x = i))}))
    diag_rho <- which(stationary_rho == 1)
    
    zeromat_beta <- drop0(Diagonal(x = rep(0, p.X)))
    
    mapping_new_Z <- do.call('cbind', store_design_Z)
    
    mapping_J <- split(1:sum(d_j[!spline_REs]^2), rep(1:length(d_j[!spline_REs]), d_j[!spline_REs]^2))
    mapping_J <- lapply(mapping_J, FUN=function(i){i-1})
    mapping_J <- sapply(mapping_J, min)

    mapping_to_re <- purrr::array_branch(do.call('cbind', store_assignment_Z), margin = 1)
    
    mapping_to_re <- lapply(mapping_to_re, FUN=function(i){
      mapply(outer_alpha_RE_positions[!spline_REs], i, SIMPLIFY = FALSE, 
          FUN=function(a,b){a[[b]]})
    })
    Mmap <- do.call('rbind', lapply(mapping_to_re, FUN=function(i){as.integer(sapply(i, min))}))

    start_base_Z <- cumsum(c(0,d_j[!spline_REs]))[-(number_of_RE - sum(spline_REs) +1)]
    names(start_base_Z) <- NULL

    store_re_id <- store_id <- list()
    id_range <- 1:nrow(Mmap)
    for (j in 1:(number_of_RE - sum(spline_REs))){
      store_re_id_j <- store_id_j <- list()
      for (jprime in 1:j){
        # print(c(j, jprime))
        umap <- unique(Mmap[, c(j, jprime)])
        store_re_id_j[[jprime]] <- unlist(apply(umap, MARGIN = 1, list), recursive = F)
        # store_re_id_j[[jprime]] <- purrr::array_branch(umap, margin = 1)
        
        id_lookup <- split(id_range, paste(Mmap[,j], Mmap[,jprime]))
        id_lookup <- id_lookup[paste(umap[,1], umap[,2])]
        names(id_lookup) <- NULL
        
        # id_lookup <- lapply(1:nrow(umap), FUN=function(i){
        #   umap_r <- umap[i,]
        #   id_r <- which( (Mmap[,j] %in% umap_r[1]) & (Mmap[,jprime] %in% umap_r[2]))
        #   return(id_r)
        # })
        store_id_j[[jprime]] <- id_lookup
      }
      store_id[[j]] <- store_id_j
      store_re_id[[j]] <- store_re_id_j
    }
    
    gc()
    
    if (do_timing){
      toc(quiet = verbose_time, log = T)
    }
  }
  store_parameter_traj <- store_vi <- store_ELBO <- data.frame()

  if (do_timing) {
    toc(quiet = verbose_time, log = TRUE)
    tic.clear()
  }
  ## Begin VI algorithm:
  if (!quiet) {
    message("Begin Regression")
  }
  
  do_SQUAREM <- control$do_SQUAREM
  
  vi_collapsed_P <- matrix(0, nrow = 1, ncol = 0)
  prop_vi_collapsed_P <- vi_collapsed_P
  adjust_ceoa <- FALSE
  vi_alpha_var <- vi_beta_var <- NULL
  
  if (factorization_method == "collapsed"){
    
    warning('Turning off SQUAREM for "collapsed')
    do_SQUAREM <- FALSE
    
    if (family == 'binomial'){
      init_c <- as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean)
      init_mean <- vi_pg_b / (2 * init_c) * tanh(init_c / 2)
      fill_zero <- which(abs(init_c) < 1e-6)
      if (length(fill_zero) > 0){
        init_mean[fill_zero] <- vi_pg_b[fill_zero] / 4
      }
    }else{
      init_c <- NULL
      init_mean <- rep(1, nrow(Z)) 
    }

    list2env(build_collapse_index(X = X, Z = Z, weight = vi_pg_b,
              cyclical_pos = cyclical_pos, 
              names_of_RE = names_of_RE, k = control$collapse_size),
             envir = base::environment())
    block_collapse <- control$block_collapse
    position_block_j <- rep(0:number_of_RE, lengths(C_j))
    
    if (block_collapse){
      orig_Cj <- C_j
      block_j <- unlist(C_j)
      if (length(block_j) == 0){
        block_collapse <- FALSE
      }else{
        names(block_j) <- NULL
        M_j <- c(list('...blocked' = block_j), M_j)
        C_j <- lapply(c('...blocked', names(C_j)), FUN=function(i){numeric(0)})
        names(C_j) <- names(M_j)
        
        position_of_blocked <- which(names(M_j) == '...blocked')
      }
    }

    joint_XZ <- cbind(X,Z)
    stopifnot(all(joint_XZ[,1] == 1))
    
    if (1 %in% C_j[[1]]){
      if (C_j[[1]][1] != 1){stop('Intercept in wrong position.')}
      intercept_in_Cj <- TRUE
    }else{
      if (M_j[[1]][1] != 1){stop('Intercept in wrong position.')}
      intercept_in_Cj <- FALSE
    }
    
    design_C <- joint_XZ[, unlist(C_j), drop = F]
    design_M <- joint_XZ[, unlist(M_j), drop = F]
    
    any_collapsed_C <- ncol(design_C) > 0
    any_collapsed_M <- ncol(design_M) > 0
    stopifnot(all(d_j == 1))
    stopifnot(all(spline_REs == FALSE))
    
    vi_C_mean <- c(vi_beta_mean, vi_alpha_mean)[unlist(C_j)]
    vi_M_mean <- c(vi_beta_mean, vi_alpha_mean)[unlist(M_j)]
    
    id_Mj <- rep(seq_len(length(M_j)) - 1, lengths(M_j))
    c_id_Mj <- as.character(seq_len(length(M_j)) - 1)
    keep_id_Mj <- rep(TRUE, length(M_j))
    names(keep_id_Mj) <- names(M_j)
    keep_id_Mj[names(keep_id_Mj) %in% c('...blocked', 'Fixed Effect')] <- FALSE
    vi_M_mean <- split(as.vector(vi_M_mean), id_Mj)
    vi_M_mean <- vi_M_mean[c_id_Mj]
    names(vi_M_mean) <- NULL
    
    rm(vi_alpha_mean, vi_beta_mean)
    vi_alpha_mean <- vi_beta_mean <- NULL
    
    joint_V_init <- bdiag(crossprod(vi_beta_decomp), crossprod(vi_alpha_decomp))
    vi_C_var <- joint_V_init[unlist(C_j), unlist(C_j), drop = F]
    vi_M_var <- lapply(M_j, FUN=function(i){as(joint_V_init[i,i,drop=F], 'dgCMatrix')})
    vi_M_list <- lapply(M_j, FUN=function(i){joint_XZ[,i, drop = F]})
    rm(joint_V_init, joint_XZ)
    
    vi_P <- lapply(vi_M_list, FUN=function(i){
      vi_C_var %*% 
        t(design_C) %*% Diagonal(x = init_mean) %*% i
    })
    
    vi_C_uncond <- vi_C_var + 
      Reduce('+', mapply(vi_P, vi_M_var, SIMPLIFY = FALSE, 
        FUN=function(pi_m, Mi){pi_m %*% Mi %*% t(pi_m)}))
    vi_C_uncond <- as(vi_C_uncond, 'dgCMatrix')
    
    vi_var_M_list <- lapply(vi_M_var, FUN=function(i){
      if (ncol(i) == 0){return(NULL)}
      o <- diag(i)
      o <- sparseMatrix(i = seq_len(length(o)), j = seq_len(length(o)), x = o)
      return(o)
    })


    lagged_M_var <- vi_M_var
    lagged_C_var <- vi_C_var
    
    lagged_C_mean <- rep(-Inf, length(unlist(C_j)))
    lagged_M_mean <- rep(-Inf, length(unlist(M_j)))
    lagged_sigma_alpha <- vi_sigma_alpha
    
    rm(init_c, init_mean)
    
  }
  
  if (debug_param) {
    if (factorization_method == "collapsed"){
      store_beta <- array(NA, dim = c(iterations, sum(lengths(C_j))))
      store_alpha <- array(NA, dim = c(iterations, sum(lengths(M_j))))
      store_sigma <- array(NA, dim = c(iterations, sum(d_j^2)))
    }else{
      store_beta <- array(NA, dim = c(iterations, p.X))
      store_alpha <- array(NA, dim = c(iterations, p.Z))
      store_sigma <- array(NA, dim = c(iterations, sum(d_j^2)))
    }
    if (do_huangwand){
      store_hw <- array(NA, dim = c(iterations, length(unlist(vi_a_b_jp))))
    }
  }
  
  if (family %in% c('negbin', 'linear')){
    if (do_SQUAREM){warning('Turning off SQUAREM for negbin/linear temporarily.')}
    do_SQUAREM <- FALSE
  }
  if (family == 'negbin' & !(control$vi_r_method %in% c('VEM', 'fixed'))){
    if (do_SQUAREM){warning('Turning off SQUAREM if "negbin" and not VEM/fixed.')}
    do_SQUAREM <- FALSE
  }

  if (do_SQUAREM){
    squarem_success <- c(0, 0)
    squarem_list <- list()
    squarem_counter <- 1
  }else{
    squarem_success <- NA
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
    
    if (family != 'linear'){
      if (factorization_method == "weak") {
        # joint_quad <- rowSums( (joint.XZ %*% t(vi_joint_decomp))^2 )
        # vi_joint_decomp <<- vi_joint_decomp
        # joint.XZ <<- joint.XZ
        joint_quad <- cpp_zVz(Z = joint.XZ, V = as(vi_joint_decomp, "dgCMatrix")) 
        if (family == 'negbin'){
          joint_quad <- joint_quad + vi_r_sigma
        }
        vi_pg_c <- sqrt(as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean - vi_r_mu)^2 + joint_quad)
      } else if (factorization_method == "collapsed"){
      
        # joint_quad <- rowSums( (design_C %*% vi_C_uncond) * design_C)
        # for (j in seq_len(length(M_j))){
        #   var_Mj <- vi_M_var[[j]]
        #   data_Mj <- vi_M_list[[j]]
        #   jq1 <- rowSums((data_Mj %*% var_Mj) * data_Mj)
        #   jq2 <- rowSums( (data_Mj %*% t(vi_P[[j]] %*% var_Mj)) * design_C)
        #   joint_quad <- joint_quad + jq1 - 2 * jq2
        # }
        
        joint_quad <- cpp_var_lp(design_C = design_C,
          vi_C_uncond = vi_C_uncond,
          vi_M_var = vi_M_var,
          vi_M_list = vi_M_list,
          vi_P = vi_P,
          sparse_input = class(vi_M_var[[1]]) == 'dgCMatrix'
        )
        
        if (family == 'negbin'){
          joint_quad <- joint_quad + vi_r_sigma
        }
        
        vi_pg_c <- design_C %*% vi_C_mean
        if (any_collapsed_M){
          vi_pg_c <- vi_pg_c + design_M %*% do.call('c', vi_M_mean)
        }
        vi_pg_c <- sqrt(as.vector(vi_pg_c - vi_r_mu)^2 + joint_quad)
        
      } else {
        
        beta_quad <- rowSums((X %*% t(vi_beta_decomp))^2)
        alpha_quad <- rowSums((Z %*% t(vi_alpha_decomp))^2)
        joint_var <- beta_quad + alpha_quad
        
        if (family == 'negbin'){
          joint_var <- joint_var + vi_r_sigma
        }
        vi_pg_c <- sqrt(as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean - vi_r_mu)^2 + joint_var)
      }
    }
    
    if (family == 'linear'){
      vi_pg_mean <- rep(1, nrow(X))
      diag_vi_pg_mean <- sparseMatrix(i = 1:N, j = 1:N, x = vi_pg_mean)
    }else{
      vi_pg_mean <- vi_pg_b / (2 * vi_pg_c) * tanh(vi_pg_c / 2)
      fill_zero <- which(abs(vi_pg_c) < 1e-6)
      if (length(fill_zero) > 0){
        vi_pg_mean[fill_zero] <- vi_pg_b[fill_zero] / 4
      }
      diag_vi_pg_mean <- sparseMatrix(i = 1:N, j = 1:N, x = vi_pg_mean)
    }
    
    if (debug_ELBO & it != 1) {
      debug_ELBO.1 <- calculate_ELBO(family = family,
        ELBO_type = ELBO_type, 
        factorization_method = factorization_method,
        d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi,
        prior_sigma_alpha_nu = prior_sigma_alpha_nu,
        iw_prior_constant = iw_prior_constant,
        store_assignment_Z = store_assignment_Z, 
        store_design_Z = store_design_Z, 
        X = X, Z = Z, s = s, y = y, outer_alpha_RE_positions = outer_alpha_RE_positions,
        vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
        vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu,
        vi_sigma_outer_alpha = vi_sigma_outer_alpha,
        vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
        log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
        vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp,
        vi_alpha_var = vi_alpha_var, vi_beta_var = vi_beta_var, cyclical_pos = cyclical_pos,
        vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
        vi_P = vi_P, vi_C_uncond = vi_C_uncond, vi_M_var = vi_M_var, vi_M_list = vi_M_list,
        vi_M_mean = vi_M_mean, vi_C_mean = vi_C_mean, 
        design_M = design_M,
        design_C = design_C, log_det_M_var = log_det_M_var, log_det_C_var = log_det_C_var,
        vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
        vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
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
    if (factorization_method %in% c("strong")) {
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
    
    if (factorization_method == "collapsed"){
      Tinv <- bdiag(Diagonal(x = rep(0, p.X)), Tinv)
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

      vi_joint_L_nonpermute <- drop0(solve(chol.update.joint$origL))
      vi_joint_LP <- Pmatrix
      vi_joint_decomp <- vi_joint_L_nonpermute %*% t(vi_joint_LP)

      vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
      vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))

      vi_alpha_decomp <- vi_joint_decomp[, -1:-p.X, drop = F]
      vi_beta_decomp <- vi_joint_decomp[, 1:p.X, drop = F]

      log_det_joint_var <- -2 * sum(log(diag(chol.update.joint$origL)))
      if (do_SQUAREM){
        vi_joint_L_nonpermute <- vi_joint_decomp
        vi_joint_LP <- Diagonal(n = ncol(vi_joint_decomp))
      }
    } else if (factorization_method == "collapsed") {
      
      running_log_det_M_var <- rep(0, length(M_j))
      
      if (do_timing){
        tic('beta C')
      }
      if (any_collapsed_C){
        
        
        Tinv_C <- Tinv[unlist(C_j), unlist(C_j), drop = F]
        
        # chol.update.C <- Cholesky(t(design_C) %*% diag_vi_pg_mean %*% design_C + Tinv_C)
        # C_hat <- as.vector(solve(chol.update.C, t(design_C) %*% s))
        # vi_P <- lapply(vi_M_list, FUN=function(i){
        #   solve(chol.update.C, t(design_C) %*% diag_vi_pg_mean %*% i)
        # })
        # vi_C_var <- with(expand(chol.update.C), crossprod(solve(L) %*% P))
        # # vi_C_var_alt <- solve(t(design_C) %*% diag_vi_pg_mean %*% design_C + Tinv_C)
        # log_det_C_var <- -2 * as.numeric(determinant(chol.update.C)$modulus)
        
        update_C <- cpp_update_c_var(diag_vi_pg_mean = diag_vi_pg_mean,
                         design_C = design_C, Tinv_C = Tinv_C, s = s,
                         vi_M_list = vi_M_list)
        
        list2env(update_C, envir = environment())
        rm(update_C)
        
      }else{
        Tinv_C <- vi_C_var <- sparseMatrix(i = numeric(0), j = numeric(0), x = numeric(0), dims = c(0,0))
        vi_P <- lapply(vi_M_list, FUN=function(i){
          sparseMatrix(i = numeric(0), j = numeric(0), x = numeric(0), dims = c(0, ncol(i)))
        })
        C_hat <- numeric(0)
        log_det_C_var <- 0
      }
      
      if (do_timing){
        toc(quiet = verbose_time, log = TRUE)
      }
      if (linpred_method == 'cyclical'){
        
        pg_lp_C_hat <- diag_vi_pg_mean %*% design_C %*% C_hat

        running_Cadjust <- rowSums(do.call('cbind',
           mapply(vi_P, vi_M_mean, SIMPLIFY = FALSE,
                  FUN=function(i,j){if (ncol(i) > 0){i %*% j}})))

        running_lp <- as.vector(
          design_M %*% do.call('c', vi_M_mean) - design_C %*% running_Cadjust
        )
        
        running_log_det_M_var <- rep(0, length(M_j))
        vi_M_var <- vector(mode = "list", length(M_j))
        
        for (j in seq_len(length(M_j))) {
          index_Mj <- M_j[[j]]
          if (length(index_Mj) == 0){
            vi_M_var[[j]] <- sparseMatrix(i = numeric(0), j = numeric(0), x = numeric(0), dims = c(0,0))
            next
          }
          data_M_j <- vi_M_list[[j]]
          
          Ainv <- solve(t(data_M_j) %*% diag_vi_pg_mean %*% data_M_j + bdiag(Tinv)[index_Mj, index_Mj])
          
          if (any_collapsed_C){
            bread_j <- t(data_M_j) %*% diag_vi_pg_mean %*% design_C
            meat_woodbury <- t(design_C) %*% diag_vi_pg_mean %*% design_C -
              t(bread_j) %*% Ainv %*% bread_j + Tinv_C
            var_M_j <- Ainv %*% bread_j %*% solve(meat_woodbury) %*% t(bread_j) %*% Ainv
            running_log_det_M_var[j] <- -as.numeric(determinant(meat_woodbury)$modulus +
                                                      determinant(vi_C_var)$modulus - determinant(Ainv)$modulus)
          }else{
            var_M_j <- Ainv 
            running_log_det_M_var[j] <- as.numeric(determinant(Ainv)$modulus)
          }
          
          
          # data_projM_j <- data_M_j - design_C %*% vi_P[[j]]
          # var_M_j <- solve(t(data_projM_j) %*% diag_vi_pg_mean %*% data_projM_j +
          #     bdiag(Tinv)[index_Mj, index_Mj] + t(vi_P[[j]]) %*% Tinv_C %*% vi_P[[j]])
          # 
          # running_log_det_M_var[j] <- as.numeric(determinant(var_M_j)$modulus)
          
          adj_s <- s - pg_lp_C_hat -  diag_vi_pg_mean %*% (running_lp - (data_M_j %*% vi_M_mean[[j]] - design_C %*% vi_P[[j]] %*% vi_M_mean[[j]]))
          adj_C <- C_hat - (running_Cadjust - vi_P[[j]] %*% vi_M_mean[[j]])

          # other_termRHS1 <- Reduce('+', mapply(vi_M_list[-j], vi_M_mean[-j], SIMPLIFY = FALSE, FUN=function(i,j){if (ncol(i) > 0){i %*% j}else{0}}))
          # other_termRHS1 <- other_termRHS1 - design_C %*% Reduce('+', mapply(vi_P[-j], vi_M_mean[-j], SIMPLIFY = FALSE, FUN=function(i,j){if (ncol(i) > 0){i %*% j}else{0}}))
          # 
          # adj_salt <- s - diag_vi_pg_mean %*% design_C %*% C_hat - diag_vi_pg_mean %*% other_termRHS1
          # adj_Calt <- C_hat - Reduce('+', mapply(vi_P[-j], vi_M_mean[-j], FUN=function(i,j){if (ncol(i) > 0){i %*% j}else{0}}))
          # 
          # adj_s <- adj_salt
          # adj_C <- adj_Calt
          
          RHS_1 <- t(data_M_j) %*% adj_s - t(design_C %*% vi_P[[j]]) %*% adj_s
          RHS_2 <- t(vi_P[[j]]) %*% Tinv_C %*% adj_C
          
          update_Mj <- var_M_j %*% (RHS_1 + RHS_2)
          
          running_lp <- running_lp + data_M_j %*% (update_Mj - vi_M_mean[[j]]) +
            - design_C %*% vi_P[[j]] %*% (update_Mj - vi_M_mean[[j]])
          running_Cadjust <- running_Cadjust + vi_P[[j]] %*% (update_Mj - vi_M_mean[[j]])
          
          vi_M_var[[j]] <- as(var_M_j, 'dgCMatrix')
          vi_M_mean[[j]] <- as.vector(update_Mj)
          
        }
        
        vi_M_mean <- do.call('c', vi_M_mean)
        
      }else if (linpred_method == 'joint'){
        
        all_P <- do.call('cbind', vi_P)
        Tinv_M <- Tinv[unlist(M_j), unlist(M_j), drop = F]
        
        if (do_timing){
          
          tic('beta mean')
        }
        
        if (TRUE){
          
          if (any_collapsed_M){
            cg_alpha <- cg_custom(Z = design_M[exclude_collapsed,], 
                                  P = as.matrix(all_P), X = drop0(design_C[exclude_collapsed,]),
                                  it_max = 25, 
                                  tol = sqrt(.Machine$double.eps), ridge_X = Tinv_C,
                                  omega = vi_pg_mean[exclude_collapsed], ridge_Z = Tinv_M, 
                                  s = as.vector(s - diag_vi_pg_mean %*% design_C %*% C_hat)[exclude_collapsed], 
                                  old_alpha = do.call('c', vi_M_mean),
                                  offset_ridge_X = as.vector(t(all_P) %*% Tinv_C %*% C_hat))
            
            vi_M_mean <- matrix(cg_alpha$alpha)
          }
        }else{
          Mproj <- design_M - design_C %*% do.call('cbind', vi_P)
          
          direct_M <- solve(
            t(all_P) %*% Tinv_C %*% all_P + t(Mproj) %*% diag_vi_pg_mean %*% Mproj + Tinv_M,
            t(Mproj) %*% (s - diag_vi_pg_mean %*% design_C %*% C_hat) +
              + t(all_P) %*% Tinv_C %*% C_hat
          )
          vi_M_mean <- direct_M
        }
        
        if (do_timing){
          toc(quiet = verbose_time, log = TRUE)
          tic('beta var')
        }
        if (any_collapsed_M){
          
          # # Base R Loop Version
          # for (j in seq_len(length(M_j))) {
          # 
          #   index_Mj <- M_j[[j]]
          # 
          #   if (length(index_Mj) == 0){
          #     vi_M_var[[j]] <- drop0(matrix(data = 1, nrow = 0, ncol = 0))
          #     running_log_det_M_var[j] <- 0
          #     next
          #   }
          #   # Z_j^T Omega Z_j + R - Z_j^T Omega X (X^T Omega X)^{-1} X^T Omega Z_j
          #   # A = Z_j^T Omega Z_j + R
          #   # U = - Z_j^T Omega X
          #   # C = X^T Omega X^{-1}
          #   # V = X^T Omega Z_j
          #   # A^{-1} - A^{-1} U (C^{-1} + V A^{-1} U) V A^{-1}
          #   data_M_j <- vi_M_list[[j]]
          # 
          #   Ainv <- solve(t(data_M_j) %*% diag_vi_pg_mean %*% data_M_j + bdiag(Tinv)[index_Mj, index_Mj])
          #   if (any_collapsed_C){
          #     bread_j <- t(data_M_j) %*% diag_vi_pg_mean %*% design_C
          #     meat_woodbury <- t(design_C) %*% diag_vi_pg_mean %*% design_C -
          #       t(bread_j) %*% Ainv %*% bread_j + Tinv_C
          #     var_M_j <- Ainv + Ainv %*% bread_j %*% solve(meat_woodbury) %*% t(bread_j) %*% Ainv
          #     running_log_det_M_var[j] <- -as.numeric(determinant(meat_woodbury)$modulus +
          #                                               log_det_C_var - determinant(Ainv)$modulus)
          #   }else{
          #     var_M_j <- Ainv
          #     running_log_det_M_var[j] <- as.numeric(determinant(Ainv)$modulus)
          #   }
          #   # data_projM_j <- data_M_j - design_C %*% vi_P[[j]]
          #   # var_M_j <- solve(t(data_projM_j) %*% diag_vi_pg_mean %*% data_projM_j +
          #   #                    bdiag(Tinv)[index_Mj, index_Mj] + t(vi_P[[j]]) %*% Tinv_C %*% vi_P[[j]])
          #   # running_log_det_M_var[j] <- as.numeric(determinant(var_M_j)$modulus)
          # 
          #   vi_M_var[[j]] <- as(var_M_j, 'dgCMatrix')
          # 
          # }
          # 

          Tinv_M <- bdiag(Tinv)
          Tinv_M <- lapply(M_j, FUN=function(i){Tinv_M[i,i, drop = F]})
          
          update_vi_M <- cpp_update_m_var(diag_vi_pg_mean = diag_vi_pg_mean,
            design_C = design_C, Tinv_C = Tinv_C, list_Tinv_M = Tinv_M,
            vi_M_list = vi_M_list, any_collapsed_C = any_collapsed_C,
            lndet_C = log_det_C_var)

          running_log_det_M_var <- update_vi_M$running_log_det_M_var
          vi_M_var <- update_vi_M$vi_M_var
          rm(update_vi_M)

        }
        
        if (do_timing){
          toc(quiet = verbose_time, log = TRUE)
        }
        
      }else{stop('linpred_method must be "joint" or "cyclical" if collapsed.')}

      # Turn vi_M_mean into "list"
      if (any_collapsed_M){
        vi_M_mean <- split(as.vector(vi_M_mean), id_Mj)
        vi_M_mean <- vi_M_mean[c_id_Mj]
        names(vi_M_mean) <- NULL
      }
      # Get E[alpha_C] = E[E[alpha_C | alpha_M]]
      if (any_collapsed_C){
        vi_C_mean <- C_hat - Reduce('+', mapply(vi_P, vi_M_mean, FUN=function(i,j){if (!is.null(j)){i %*% j}else{0}}))
      }
      
      vi_C_uncond <- vi_C_var + 
        Reduce('+', mapply(vi_P, vi_M_var, SIMPLIFY = FALSE, 
                           FUN=function(pi_m, Mi){pi_m %*% Mi %*% t(pi_m)}))
      vi_C_uncond <- as(vi_C_uncond, 'dgCMatrix')
      log_det_M_var <- sum(running_log_det_M_var)
      
      log_det_joint_var <- NA
      
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

        vi_beta_L_nonpermute <- drop0(solve(chol.update.beta$origL))
        vi_beta_LP <- Pmatrix
        vi_beta_decomp <- vi_beta_L_nonpermute %*% t(vi_beta_LP)
        vi_beta_mean <- chol.update.beta$mean
        log_det_beta_var <- -2 * sum(log(diag(chol.update.beta$origL)))

        chol.update.alpha <- LinRegChol(
          X = Z, omega = diag_vi_pg_mean, prior_precision = Tinv,
          y = as.vector(s - diag_vi_pg_mean %*% X %*% vi_beta_mean)
        )
        Pmatrix <- sparseMatrix(i = 1:p.Z, j = 1 + chol.update.alpha$Pindex, x = 1)

        vi_alpha_L_nonpermute <- drop0(solve(chol.update.alpha$origL))
        vi_alpha_LP <- Pmatrix
        vi_alpha_decomp <- vi_alpha_L_nonpermute  %*% t(vi_alpha_LP)
        vi_alpha_decomp <- drop0(vi_alpha_decomp)
        vi_alpha_decomp <- as(vi_alpha_decomp, 'dgCMatrix')
        vi_alpha_mean <- chol.update.alpha$mean
        log_det_alpha_var <- -2 * sum(log(diag(chol.update.alpha$origL)))

        vi_beta_mean <- Matrix(vi_beta_mean, dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(vi_alpha_mean, dimnames = list(fmt_names_Z, NULL))
      } else if (linpred_method == "joint") {
        joint.XZ <- cbind(X, Z)

        
        chol.update.joint <- solve(Matrix::Cholesky(  
          crossprod(Diagonal(x = sqrt(vi_pg_mean)) %*% joint.XZ) + 
            bdiag(zero_mat, bdiag(Tinv)) ),
               t(joint.XZ) %*% (s + vi_pg_mean * vi_r_mu) )
        vi_beta_mean <- Matrix(chol.update.joint[1:p.X,], dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(chol.update.joint[-1:-p.X,], dimnames = list(fmt_names_Z, NULL))
        
        # chol.update.joint <- LinRegChol(
        #   X = joint.XZ, omega = diag_vi_pg_mean, prior_precision = bdiag(zero_mat, Tinv),
        #   y = s + vi_pg_mean * vi_r_mu,
        #   save_chol = FALSE
        # )
        # vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
        # vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))
        
        vi_beta_decomp <- solve(t(chol(as.matrix(t(X) %*% diag_vi_pg_mean %*% X))))
        
        vi_beta_L_nonpermute <- vi_beta_decomp
        vi_beta_LP <- Diagonal(n = ncol(vi_beta_decomp))
        
        log_det_beta_var <- 2 * sum(log(diag(vi_beta_decomp)))

        chol.update.alpha <- LinRegChol(
          X = Z, omega = diag_vi_pg_mean, prior_precision = Tinv,
          y = s + vi_pg_mean * vi_r_mu
        )
        Pmatrix <- sparseMatrix(i = 1:p.Z, j = 1 + chol.update.alpha$Pindex, x = 1)

        vi_alpha_L_nonpermute <- drop0(solve(chol.update.alpha$origL))
        vi_alpha_LP <- Pmatrix
        
        vi_alpha_decomp <- vi_alpha_L_nonpermute %*% t(vi_alpha_LP)
        vi_alpha_decomp <- drop0(vi_alpha_decomp)
        log_det_alpha_var <- -2 * sum(log(diag(chol.update.alpha$origL)))
        
        if (do_SQUAREM){
          vi_alpha_L_nonpermute <- vi_alpha_decomp
          vi_alpha_LP <- Diagonal(n = ncol(vi_alpha_decomp))
        }
        
      } else {
        stop("Invalid linpred method for partial scheme")
      }
    } else if (factorization_method == "strong") {
      running_log_det_alpha_var <- rep(NA, number_of_RE)
      vi_alpha_decomp <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(p.Z, 2))
      
      if (linpred_method == "joint") {
        if (it == 1){
          joint.XZ <- cbind(X, Z)
        }
        if (do_timing) {
          tic("ux_mean")
        }
        
        sqrt_pg_weights <- Diagonal(x = sqrt(vi_pg_mean))
        chol.update.joint <- solve(Matrix::Cholesky(  
          crossprod(sqrt_pg_weights %*% joint.XZ) + 
            bdiag(zero_mat, bdiag(Tinv)) ),
          t(joint.XZ) %*% (s + vi_pg_mean * vi_r_mu) )
        
        vi_beta_mean <- Matrix(chol.update.joint[1:p.X,], dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(chol.update.joint[-1:-p.X,], dimnames = list(fmt_names_Z, NULL))

        # chol.update.joint <- LinRegChol(X = joint.XZ,
        #   omega = diag_vi_pg_mean,
        #   prior_precision = bdiag(zero_mat, bdiag(Tinv)),
        #   y = s + vi_pg_mean * vi_r_mu,
        #   save_chol = FALSE)
        # vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
        # vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))

        if (do_timing) {
          toc(quiet = verbose_time, log = T)
          tic("ux_var")
        }
        
        vi_beta_decomp <- solve(t(chol(as.matrix(t(X) %*% diag_vi_pg_mean %*% X))))
        
        vi_beta_L_nonpermute <- vi_beta_decomp
        vi_beta_LP <- Diagonal(n = nrow(vi_beta_decomp))
        log_det_beta_var <- 2 * sum(log(diag(vi_beta_decomp)))
        #-log(det(t(X) %*% diag_vi_pg_mean %*% X))

        running_log_det_alpha_var <- rep(NA, number_of_RE)

        # cyclical_pos <<- cyclical_pos
        # number_of_RE <<- number_of_RE
        # Z <<- Z
        # X <<- X
        # diag_vi_pg_mean <<- diag_vi_pg_mean
        # Tinv <<- Tinv
        # running_log_det_alpha_var <<- running_log_det_alpha_var
        # vi_alpha_decomp <<- vi_alpha_decomp
        # stop()
        
        for (j in 1:number_of_RE) {
          index_j <- cyclical_pos[[j]]
          Z_j <- Z[, index_j, drop = F]
          prec_j <- crossprod(sqrt_pg_weights %*% Z_j) + Tinv[[j]]

          chol_var_j <- solve(t(chol(prec_j)))
          running_log_det_alpha_var[j] <- 2 * sum(log(diag(chol_var_j)))

          vi_alpha_decomp[index_j, index_j] <- as(chol_var_j, "dgTMatrix")
        }
        vi_alpha_L_nonpermute <- vi_alpha_decomp
        vi_alpha_LP <- Diagonal(n = nrow(vi_alpha_L_nonpermute))
        
        log_det_alpha_var <- sum(running_log_det_alpha_var)
        
        if (do_timing){
          toc(quiet = verbose_time, log = T)
        }
      } else if (linpred_method == "solve_normal") {
        bind_rhs_j <- list()
        bind_lhs_j <- list()

        for (j in 1:number_of_RE) {
          index_j <- cyclical_pos[[j]]
          Z_j <- Z[, index_j, drop = F]
          Z_negj <- Z[, -index_j, drop = F]
          prec_j <- crossprod(sqrt_pg_weights %*% Z_j) + Tinv[[j]]

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
        # vi_alpha_mean <- solve(bind_lhs_j[,1:p.Z], bind_rhs_j)
        #
        # vi_alpha_mean <- Matrix(vi_alpha_mean)
        # vi_beta_mean <- Matrix(vi_beta_mean)

        bind_lhs_j <- drop0(rbind(bind_lhs_j, cbind(vi_beta_var %*% t(X) %*% diag_vi_pg_mean %*% Z, Diagonal(n = p.X))))
        bind_rhs_j <- rbind(bind_rhs_j, vi_beta_var %*% t(X) %*% s)
        #
        bind_solution <- solve(bind_lhs_j) %*% bind_rhs_j
        # print(cbind(bind_solution, rbind(vi_alpha_mean, vi_beta_mean)))
        #
        vi_beta_mean <- Matrix(bind_solution[-1:-p.Z], dimnames = list(colnames(X), NULL))
        vi_alpha_mean <- Matrix(bind_solution[1:p.Z], dimnames = list(fmt_names_Z, NULL))
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

        vi_alpha_L_nonpermute <- vi_alpha_decomp
        vi_alpha_LP <- Diagonal(n = ncol(vi_alpha_L_nonpermute))
        # vi_alpha_decomp <- bdiag(vi_alpha_decomp)
        log_det_alpha_var <- sum(running_log_det_alpha_var)

        chol.update.beta <- LinRegChol(
          X = as(X, "sparseMatrix"), omega = diag_vi_pg_mean, prior_precision = zero_mat,
          y = as.vector(s + vi_pg_mean * vi_r_mu - diag_vi_pg_mean %*% Z %*% vi_alpha_mean)
        )
        Pmatrix <- sparseMatrix(i = 1:p.X, j = 1 + chol.update.beta$Pindex, x = 1)

        vi_beta_L_nonpermute <- drop0(solve(chol.update.beta$origL))
        vi_beta_LP <- Pmatrix
        
        vi_beta_decomp <- vi_beta_L_nonpermute %*% t(Pmatrix)
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
      
      vi_alpha_var <- vi_alpha_var * adjust_var^2
      vi_beta_var <- vi_beta_var * adjust_var^2
      if (factorization_method == "collapsed"){
        vi_alpha_var_list <- lapply(vi_alpha_var_list, FUN=function(i){adjust_var^2 * i})
      }
      
      e_ln_sigmasq <- log(vi_sigmasq_b) - digamma(vi_sigmasq_a)

      log_det_joint_var <- log_det_joint_var + (p.X + p.Z) * e_ln_sigmasq
      log_det_beta_var <- log_det_beta_var + p.X * e_ln_sigmasq
      log_det_alpha_var <- log_det_alpha_var + p.Z * e_ln_sigmasq
    }

    if (debug_ELBO & it != 1) {
      if (factorization_method == "collapsed"){
        
        ssq_M <- sapply(vi_M_var, FUN=function(i){sum(diag(i))}) + 
          sapply(vi_M_mean, FUN=function(i){sum(i^2)})
        ssq_C <- sapply(
          split((as.vector(vi_C_mean)^2 + diag(vi_C_uncond)), position_block_j)[c_id_Mj], 
          sum)
        
        if (block_collapse){
          block_terms <- sapply(split(diag(vi_M_var[[position_of_blocked]]) + vi_M_mean[[position_of_blocked]]^2, position_block_j), sum)
          block_terms <- block_terms[as.character(0:number_of_RE)]
          block_terms[is.na(block_terms)] <- 0
          ssq_C <- c(NA, block_terms)
        }
        
        variance_by_alpha_jg <- list()
        variance_by_alpha_jg$outer_alpha <- as.list(ssq_M + ssq_C)[keep_id_Mj]
        names(variance_by_alpha_jg$outer_alpha) <- names(names_of_RE)
        vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
      }else{
        variance_by_alpha_jg <- calculate_expected_outer_alpha(
          factorization_method = factorization_method,
          alpha_decomp_var = vi_alpha_decomp, 
          alpha_mu = as.vector(vi_alpha_mean), alpha_var = vi_alpha_var,
          re_position_list = outer_alpha_RE_positions,
          tP = t(vi_collapsed_P), L_beta = as.matrix(vi_beta_decomp), do_adjustment = adjust_ceoa)
        vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
      }
      
      debug_ELBO.2 <- calculate_ELBO(family = family,
        ELBO_type = ELBO_type,
        factorization_method = factorization_method,
        d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi,
        prior_sigma_alpha_nu = prior_sigma_alpha_nu,
        iw_prior_constant = iw_prior_constant,
        store_assignment_Z = store_assignment_Z, 
        store_design_Z = store_design_Z, 
        X = X, Z = Z, s = s, y = y, outer_alpha_RE_positions = outer_alpha_RE_positions,
        vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
        vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu,
        vi_sigma_outer_alpha = vi_sigma_outer_alpha,
        vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
        vi_alpha_var = vi_alpha_var, vi_beta_var = vi_beta_var, cyclical_pos = cyclical_pos,
        log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
        vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp,
        vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
        vi_P = vi_P, vi_C_uncond = vi_C_uncond, vi_M_var = vi_M_var, vi_M_list = vi_M_list,
        vi_M_mean = vi_M_mean, vi_C_mean = vi_C_mean, 
        design_M = design_M,
        design_C = design_C, log_det_M_var = log_det_M_var, log_det_C_var = log_det_C_var,
        vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
        vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
        log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean,
        vi_r_sigma = vi_r_sigma,
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
      )
      
      if (it > 1){
        if (debug_ELBO.2$ELBO < debug_ELBO.1$ELBO){browser()}
      }
      
    }
    if (do_timing) {
      toc(quiet = verbose_time, log = T)
      tic("Update Sigma")
    }
    
    ###
    # Update \Sigma_j
    ##

    if (factorization_method == "collapsed"){
      
      stopifnot(all(spline_REs == FALSE))
      stopifnot(all(d_j == 1))
      
      ssq_M <- sapply(vi_M_var, FUN=function(i){sum(diag(i))}) + 
        sapply(vi_M_mean, FUN=function(i){sum(i^2)})
      
      raw_Csq <- as.vector(vi_C_mean)^2 + diag(vi_C_uncond)
      
      ssq_C <- split(raw_Csq, position_block_j)[c_id_Mj]
      ssq_C <- sapply(ssq_C, sum)
      
      if (block_collapse){
        block_terms <- sapply(split(diag(vi_M_var[[position_of_blocked]]) + 
            vi_M_mean[[position_of_blocked]]^2, position_block_j), sum)
        block_terms <- block_terms[as.character(0:number_of_RE)]
        block_terms[is.na(block_terms)] <- 0
        ssq_C <- c(NA, block_terms)
      }
      
      variance_by_alpha_jg <- list()
      variance_by_alpha_jg$outer_alpha <- as.list(ssq_M + ssq_C)[keep_id_Mj]
      names(variance_by_alpha_jg$outer_alpha) <- names(names_of_RE)
      vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
    }else{
      
      variance_by_alpha_jg <- calculate_expected_outer_alpha(
        factorization_method = factorization_method,
        alpha_decomp_var = vi_alpha_decomp, 
        alpha_mu = as.vector(vi_alpha_mean), alpha_var = vi_alpha_var,
        re_position_list = outer_alpha_RE_positions,
        L_beta = as.matrix(vi_beta_decomp), tP = t(vi_collapsed_P), 
        do_adjustment = adjust_ceoa)
      
    }
    
    
    if (!do_huangwand){#Update standard Inverse-Wishart

      vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
      vi_sigma_alpha <- mapply(vi_sigma_outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN = function(i, j) {
        i + j
      })
      
    }else{
        #Update Inverse-Wishart
        # vi_sigma_outer_alpha <<- vi_sigma_outer_alpha
        # vi_a_a_jp <<- vi_a_a_jp
        # vi_a_b_jp <<- vi_a_b_jp
        
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
      } else if (factorization_method == "collapsed"){
        
        # vi_collapsed_P_list <- lapply(cyclical_pos, FUN=function(i){vi_collapsed_P[,i]})
        # joint_quad_old <- apply(X, MARGIN = 1, FUN=function(i){as.numeric(t(i) %*% vi_beta_var %*% i)})
        # vi_collapsed_P_list <- lapply(cyclical_pos, FUN=function(i){vi_collapsed_P[,i]})
        # for (j in 1:number_of_RE){
        #   var_aj <- vi_alpha_var_list[[j]]
        #   M_j <- vi_Z_list[[j]] - X %*% vi_collapsed_P_list[[j]]
        #   joint_quad_old <- joint_quad_old + apply(M_j, MARGIN = 1, FUN=function(i){as.numeric(t(i) %*% var_aj %*% i)})
        # }
          
        joint_quad <- cpp_quad_collapsed(V = vi_alpha_var, 
                                         re_position_list = outer_alpha_RE_positions,
                                         Z_list_raw = store_design_Z,
                                         individual_assignments = store_assignment_Z,
                                         vi_beta_var = as.matrix(vi_beta_var), 
                                         P = vi_collapsed_P, X = X
        )
        joint_quad <- as.vector(joint_quad)
        
        vi_lp <- (s - as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean))^2 + joint_quad

      } else {
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
        store_assignment_Z = store_assignment_Z, 
        store_design_Z = store_design_Z, 
        X = X, Z = Z, s = s, y = y, outer_alpha_RE_positions = outer_alpha_RE_positions,
        vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
        vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu,
        vi_sigma_outer_alpha = vi_sigma_outer_alpha,
        vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
        vi_alpha_var = vi_alpha_var, vi_beta_var = vi_beta_var, cyclical_pos = cyclical_pos,
        log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
        vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp,
        vi_joint_decomp = vi_joint_decomp,
        vi_P = vi_P, vi_C_uncond = vi_C_uncond, vi_M_var = vi_M_var, vi_M_list = vi_M_list,
        vi_M_mean = vi_M_mean, vi_C_mean = vi_C_mean, 
        design_M = design_M, design_C = design_C, 
        log_det_M_var = log_det_M_var, log_det_C_var = log_det_C_var,
        log_det_joint_var = log_det_joint_var, 
        vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma, 
        choose_term = choose_term,
        vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
        vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
      )
      if (it > 1){
        if (debug_ELBO.3$ELBO < debug_ELBO.2$ELBO){browser()}
      }
    }

    if (do_timing) {
      tic("Update PX")
    }
    
    if (parameter_expansion == "none" | !any_Mprime) {
      
      accept.PX <- TRUE
      
    } else {
      
      if (parameter_expansion == "mean"){accept.PX <- TRUE}
      
      # Do a simple mean adjusted expansion.
      # Get the mean of each random effect.
      
      # Remove the "excess mean" mu_j from each random effect \alpha_{j,g}
      # and add the summd mass back to the betas.
      if (factorization_method == "collapsed"){

        stopifnot(all(spline_REs == FALSE))
        stopifnot(all(d_j == 1))
        
        if (!block_collapse){
          sum_M <- sapply(vi_M_mean, sum)
          if (!any_collapsed_C){
            sum_all <- sum_M
          }else{
            sum_C <- sapply(split(vi_C_mean, position_block_j)[c_id_Mj], sum)
            sum_all <- sum_M + sum_C
          }
        }else{
          sum_all <- sapply(split(vi_M_mean[[position_of_blocked]], position_block_j), sum)
          sum_all <- sum_all[as.character(0:number_of_RE)]
          sum_all[is.na(sum_all)] <- 0
          sum_all <- sum_all + sapply(vi_M_mean, sum)[-position_of_blocked]
        }
        sum_all <- sum_all[-1]
        mean_RE <- sum_all/g_j
        
        vi_M_mean[keep_id_Mj] <- mapply(vi_M_mean[keep_id_Mj], mean_RE, 
            SIMPLIFY = FALSE, 
            FUN=function(i,j){i-j})
        
        if (!block_collapse){
          if (max(lengths(C_j)) != 0){
            vi_C_mean <- vi_C_mean - rep(c(0, mean_RE), lengths(C_j))
          }
          if (intercept_in_Cj){
            vi_C_mean[1] <- vi_C_mean[1] + sum(mean_RE)
          }else{
            vi_M_mean[[1]][1] <- vi_M_mean[[1]][1] + sum(mean_RE)
          }
        }else{
          
          vi_M_mean[[position_of_blocked]] <- vi_M_mean[[position_of_blocked]] + 
            -c(0, mean_RE)[position_block_j + 1]
          
          vi_M_mean[[position_of_blocked]][1] <- vi_M_mean[[position_of_blocked]][1] +
            sum(mean_RE)
        }
        
        ssq_M <- sapply(vi_M_var, FUN=function(i){sum(diag(i))}) + sapply(vi_M_mean, FUN=function(i){sum(i^2)})
        ssq_C <- sapply(split((vi_C_mean^2 + diag(vi_C_uncond)), position_block_j)[c_id_Mj], sum)
        
        if (block_collapse){
          block_terms <- sapply(split(diag(vi_M_var[[position_of_blocked]]) + vi_M_mean[[position_of_blocked]]^2, position_block_j), sum)
          block_terms <- block_terms[as.character(0:number_of_RE)]
          block_terms[is.na(block_terms)] <- 0
          ssq_C <- c(NA, block_terms)
        }
        variance_by_alpha_jg <- list()
        variance_by_alpha_jg$outer_alpha <- as.list(ssq_M + ssq_C)[keep_id_Mj]
        
      }else{
        vi_mu_j <- t(M_prime) %*% vi_alpha_mean
        vi_alpha_mean <- vi_alpha_mean - M_prime_one %*% vi_mu_j
        vi_beta_mean <- vi_beta_mean + t(M_mu_to_beta) %*% vi_mu_j
        
        variance_by_alpha_jg <- calculate_expected_outer_alpha(
          factorization_method = factorization_method,
          alpha_decomp_var = vi_alpha_decomp, alpha_var = vi_alpha_var,
          alpha_mu = as.vector(vi_alpha_mean), 
          re_position_list = outer_alpha_RE_positions,
          L_beta = as.matrix(vi_beta_decomp), tP = t(vi_collapsed_P), 
          do_adjustment = adjust_ceoa
        )
        vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
        
      }
    }  

    quiet_rho <- control$quiet_rho
    
    if (parameter_expansion %in% c("translation", "diagonal") & skip_translate == FALSE & any_Mprime) {
      
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
        vi_alpha_var = vi_alpha_var, vi_beta_var = vi_beta_var, cyclical_pos = cyclical_pos,
        log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
        vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp,
        vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
        vi_P = vi_P, vi_C_uncond = vi_C_uncond, vi_M_var = vi_M_var, vi_M_list = vi_M_list,
        vi_M_mean = vi_M_mean, vi_C_mean = vi_C_mean, 
        design_M = design_M,
        design_C = design_C, log_det_M_var = log_det_M_var, log_det_C_var = log_det_C_var,
        vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
        vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
        log_det_joint_var = log_det_joint_var, 
        vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma,
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
      )

      if (!quiet_rho){cat('r')}

      raw_R <- R_ridge <- vecR_ridge_new(L = vi_alpha_decomp, pg_mean = diag(diag_vi_pg_mean),
        mapping_J = mapping_J, d = d_j[!spline_REs],
        store_id = store_id, store_re_id = store_re_id,
        store_design = store_design_Z, 
        diag_only = (factorization_method == 'strong'))

      if (factorization_method == 'weak'){
        stop('no Translation PX for weak yet...')
      }
      
      if (!quiet_rho){cat('r')}

      R_design <- vecR_design(alpha_mu = as.vector(vi_alpha_mean), Z = mapping_new_Z, 
        M = Mmap, mapping_J = mapping_J, d = d_j[!spline_REs],
        start_z = start_base_Z)

      moments_sigma_alpha <- mapply(vi_sigma_alpha, 
                                    vi_sigma_alpha_nu, 
                                    d_j, SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
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
      
      diag_weight <- diag_weight[!spline_REs]
      prior_weight <- prior_weight[!spline_REs]

      vec_OSL_prior <- mapply(moments_sigma_alpha[!spline_REs], diag_weight, prior_weight, 
        SIMPLIFY = FALSE, FUN=function(moment_j, phi_j, nu_j){
        as.vector(moment_j$sigma.inv %*% phi_j - nu_j * Diagonal(n = nrow(phi_j)))
      })
      vec_OSL_prior <- do.call('c', vec_OSL_prior)
      vec_OSL_prior <- matrix(c(rep(0, p.X), rep(0, size_splines), vec_OSL_prior))
      
      if (sum(spline_REs)){
        if (factorization_method == 'strong'){
          R_spline_ridge <- bdiag(Tinv[spline_REs])
        }else{
          R_spline_ridge <- Tinv[-nonspline_positions, -nonspline_positions]
        }
      }else{
        R_spline_ridge <- drop0(matrix(0, nrow = 0, ncol = 0))
      }
      #If a DIAGONAL expansion, then only update the diagonal elements
      if (parameter_expansion == "diagonal"){
        
        XR <- cbind(X, Z[,-nonspline_positions], R_design[, diag_rho])
        R_ridge <- bdiag(zeromat_beta, R_spline_ridge, R_ridge[diag_rho, diag_rho])
        
        if (do_huangwand){
          vec_OSL_prior <- do.call('c', mapply(vi_a_APRIOR_jp[!spline_REs], 
                                               vi_a_a_jp[!spline_REs], 
                                               vi_a_b_jp[!spline_REs],
                                               SIMPLIFY = FALSE,
            FUN=function(i,a,b){1-2/i^2 * a/b}))
          vec_OSL_prior <- c(rep(0, p.X + size_splines), vec_OSL_prior)
        }else{
          vec_OSL_prior <- vec_OSL_prior[c(seq_len(p.X + size_splines), p.X + size_splines + diag_rho),,drop=F]
        }
        if (length(vec_OSL_prior) != ncol(XR)){stop('MISALIGNED DIMENSIONS')}
        
        update_expansion_XR <- vecR_fast_ridge(X = drop0(XR), 
         omega = diag_vi_pg_mean, prior_precision = R_ridge, y = as.vector(s), 
         adjust_y = as.vector(vec_OSL_prior))
        
        update_expansion_bX <- Matrix(update_expansion_XR[1:p.X])
        update_expansion_splines <- Matrix(update_expansion_XR[-(1:p.X)][seq_len(size_splines)])
        
        update_expansion_R <- mapply(split(update_expansion_XR[-seq_len(p.X + size_splines)], 
          rep(1:(number_of_RE - sum(spline_REs)), d_j[!spline_REs])), d_j[!spline_REs], SIMPLIFY = FALSE, 
          FUN=function(i,d){
            dg <- diag(x = d)
            diag(dg) <- i
            return(dg)
          })
         update_diag_R <- split(update_expansion_XR[-seq_len(p.X + size_splines)], 
                                rep(1:(number_of_RE - sum(spline_REs)), d_j[!spline_REs]))
         rownames(update_expansion_bX) <- colnames(X)
      }else{#Do the FULL update
        
        
        browser()
        
        XR <- drop0(cbind(drop0(X), Z[,-nonspline_positions], drop0(R_design)))
        R_ridge <- bdiag(zeromat_beta, R_spline_ridge, R_ridge)
        
        update_expansion_XR <- vecR_fast_ridge(X = XR, 
         omega = diag_vi_pg_mean, prior_precision = R_ridge, y = as.vector(s), 
         adjust_y = as.vector(vec_OSL_prior))
        
        update_expansion_bX <- Matrix(update_expansion_XR[1:p.X])
        update_expansion_splines <- Matrix(update_expansion_XR[-(1:p.X)][seq_len(size_splines)])
        update_expansion_R <- mapply(split(update_expansion_XR[-1:-(p.X + size_splines)], 
          rep(1:(number_of_RE - sum(spline_REs)), d_j[!spline_REs]^2)), d_j[!spline_REs], 
          SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d)})
        
      }
      
      gc()

      est_rho <- update_expansion_XR[-seq_len(p.X + size_splines)]
      
      if (!quiet_rho){print(round(est_rho, 6))}
      if (parameter_expansion == 'diagonal'){
        if (max(abs(est_rho - 1)) < 1e-6){
          if (!quiet_rho){print('No further improvements')}
          skip_translate <- TRUE
        }
      }else{
        if (max(abs(est_rho - stationary_rho)) < 1e-6){
          if (!quiet_rho){print('No further improvements')}
          skip_translate <- TRUE
        }
      }
      
      if (sum(spline_REs) > 0){

        if (parameter_expansion == 'diagonal'){
          old_update_diag_R <- update_diag_R
          update_diag_R <- lapply(d_j, FUN=function(i){rep(1, i)})
          update_diag_R[!spline_REs] <- old_update_diag_R
        }
        
        old_update_expansion_R <- update_expansion_R
        update_expansion_R <- lapply(d_j, FUN=function(i){diag(i)})
        update_expansion_R[!spline_REs] <- old_update_expansion_R
        rm(old_update_expansion_R)
      }

      prop_vi_sigma_alpha <- mapply(vi_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE,
                                    FUN=function(Phi, R){R %*% Phi %*% t(R)})
      # cat('r')
      mapping_for_R_block <- make_mapping_alpha(update_expansion_R, px.R = TRUE)
      update_expansion_Rblock <- prepare_T(mapping = mapping_for_R_block, levels_per_RE = g_j, num_REs = number_of_RE,
                        variables_per_RE = d_j, running_per_RE = breaks_for_RE, cyclical = FALSE, px.R = TRUE)
      # cat('r')
      stopifnot(all.equal(update_expansion_Rblock, bdiag(mapply(update_expansion_R, g_j, FUN=function(i,g){bdiag(lapply(1:g, FUN=function(k){i}))}))))

      update_expansion_R_logdet <- sapply(update_expansion_R, FUN=function(i){determinant(i)$modulus})
      
      prop_vi_beta_mean <- update_expansion_bX
      prop_vi_alpha_mean <- update_expansion_Rblock %*% vi_alpha_mean
      
      if (!quiet_rho){cat('r')}
      
      if (factorization_method != 'weak'){
        prop_log_det_joint_var <- prop_vi_joint_decomp <- NULL
        prop_vi_alpha_decomp <- vi_alpha_decomp %*% t(update_expansion_Rblock)
        prop_log_det_alpha_var <- log_det_alpha_var + 2 * sum(update_expansion_R_logdet * g_j)
        prop_log_det_beta_var <- log_det_beta_var
        prop_vi_beta_decomp <- vi_beta_decomp

        prop_variance_by_alpha_jg <- calculate_expected_outer_alpha(
          factorization_method = factorization_method,
          alpha_decomp_var = prop_vi_alpha_decomp, alpha_var = prop_vi_alpha_var,
          alpha_mu = as.vector(prop_vi_alpha_mean),
          re_position_list = outer_alpha_RE_positions,
          L_beta = prop_vi_beta_decomp, tP = t(prop_vi_collapsed_P), 
          do_adjustment = adjust_ceoa)
        prop_vi_sigma_outer_alpha <- prop_variance_by_alpha_jg$outer_alpha
      }else{
        stop('...')
      }
      
      if (do_huangwand){
        if (parameter_expansion == "diagonal"){
          prop_vi_a_b_jp <- mapply(vi_a_b_jp, update_diag_R, SIMPLIFY = FALSE,
                              FUN=function(i,j){i / j^2})
        }else{
          # update_expansion_R <<- update_expansion_R
          
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
        vi_alpha_var = prop_vi_alpha_var, vi_beta_var = prop_vi_beta_var, cyclical_pos = cyclical_pos,
        log_det_beta_var = prop_log_det_beta_var, 
        log_det_alpha_var = prop_log_det_alpha_var,
        log_det_joint_var = prop_log_det_joint_var,
        
        vi_beta_decomp = prop_vi_beta_decomp, 
        vi_alpha_decomp = prop_vi_alpha_decomp,
        vi_joint_decomp = prop_vi_joint_decomp,
        
        vi_P = vi_P, vi_C_uncond = vi_C_uncond, vi_M_var = vi_M_var, vi_M_list = vi_M_list,
        vi_M_mean = vi_M_mean, vi_C_mean = vi_C_mean, 
        design_M = design_M,
        design_C = design_C, log_det_M_var = log_det_M_var, log_det_C_var = log_det_C_var,
        
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, 
        vi_a_b_jp = prop_vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
        choose_term
      )
      if (!quiet_rho){cat('d')}
      if (prop.ELBO$ELBO > prior.ELBO$ELBO){
        #Accept the PX-VB adjustment (OSL).
        vi_beta_mean <- prop_vi_beta_mean
        vi_alpha_mean <- prop_vi_alpha_mean
        vi_sigma_alpha <- prop_vi_sigma_alpha
        if (factorization_method == 'weak'){
          stop('Setup reassignment weak')
          if (do_SQUAREM){stop('...')}
        }else{
          vi_alpha_decomp <- prop_vi_alpha_decomp
          log_det_alpha_var <- prop_log_det_alpha_var
          if (do_SQUAREM){
            vi_alpha_L_nonpermute <- vi_alpha_decomp
            vi_alpha_LP <- Diagonal(n = ncol(vi_alpha_decomp))
          }
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
      if (!quiet_rho){print(accept.PX)}
      accepted_times <- accept.PX + accepted_times

      rm(prop_vi_beta_mean, prop_vi_alpha_mean, prop_vi_sigma_alpha, prop_vi_alpha_decomp,
         prop_log_det_alpha_var, prop_variance_by_alpha_jg, prop_vi_sigma_outer_alpha)


      rownames(vi_alpha_mean) <- fmt_names_Z
    }
    
    # Adjust the terms in the ELBO calculation that are different.

    if (accept.PX) {

      final.ELBO <- calculate_ELBO(family = family,
        ELBO_type = ELBO_type,
        factorization_method = factorization_method,
        store_assignment_Z = store_assignment_Z,
        store_design_Z = store_design_Z,
        outer_alpha_RE_positions = outer_alpha_RE_positions,
        d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi,
        prior_sigma_alpha_nu = prior_sigma_alpha_nu,
        iw_prior_constant = iw_prior_constant,
        X = X, Z = Z, s = s, y = y,
        vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
        vi_sigma_alpha = vi_sigma_alpha, vi_sigma_alpha_nu = vi_sigma_alpha_nu,
        vi_sigma_outer_alpha = vi_sigma_outer_alpha,
        vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
        vi_alpha_var = vi_alpha_var, vi_beta_var = vi_beta_var, cyclical_pos = cyclical_pos,
        log_det_beta_var = log_det_beta_var, log_det_alpha_var = log_det_alpha_var,
        vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp,
        vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
        vi_P = vi_P, vi_C_uncond = vi_C_uncond, vi_M_var = vi_M_var, vi_M_list = vi_M_list,
        vi_M_mean = vi_M_mean, vi_C_mean = vi_C_mean, 
        design_M = design_M,
        design_C = design_C, log_det_M_var = log_det_M_var, log_det_C_var = log_det_C_var,
        vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b,
        vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
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
    
    if (do_timing) {
      toc(quiet = verbose_time, log = TRUE)
      tic("Update Squarem")
    }
    
    if (do_SQUAREM){
      
      if (factorization_method %in% c('weak', 'collapsed')){
        vi_alpha_L_nonpermute <- vi_beta_L_nonpermute <- NULL
        vi_alpha_LP <- vi_beta_LP <- NULL
      }else{
        vi_joint_L_nonpermute <- vi_joint_LP <- NULL
      }
      
      squarem_list[[squarem_counter]] <- namedList(vi_sigma_alpha_nu, 
           vi_sigma_alpha, vi_alpha_mean, vi_beta_mean,
           vi_pg_c, vi_alpha_L_nonpermute, vi_alpha_LP,
           vi_beta_L_nonpermute, vi_beta_LP,
           vi_alpha_L_nonpermute,
           vi_joint_L_nonpermute, vi_joint_LP,
           vi_a_a_jp, vi_a_b_jp,
           vi_r_mu, vi_r_sigma, vi_r_mean)
      if (squarem_counter %% 3 == 0){
        
        # final.ELBO <<- final.ELBO
        # squarem_list <<- squarem_list
        # outer_alpha_RE_positions <<- outer_alpha_RE_positions
        # 
        
        ELBOargs <- list(family = family,
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
           log_det_joint_var = log_det_joint_var, 
           vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma,
           do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
           vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
        )
        
        if (factorization_method %in% c('weak', 'collapsed')){
          squarem_par <- c('vi_a_b_jp', 'vi_sigma_alpha', 'vi_pg_c',
                           'vi_alpha_mean', 'vi_beta_mean', 'vi_joint_L_nonpermute')
          squarem_type <- c('positive', 'matrix', 'positive',
                            'real', 'real', 'cholesky')
          squarem_structure <- c('list', 'list', 'vector', 'vector', 'vector',
                                 'vector')
        }else{
          squarem_par <- c('vi_a_b_jp', 'vi_sigma_alpha', 'vi_pg_c',
                           'vi_alpha_mean', 'vi_beta_mean', 'vi_beta_L_nonpermute',
                           'vi_alpha_L_nonpermute')
          squarem_type <- c('positive', 'matrix', 'positive',
                            'real', 'real', 'cholesky', 'cholesky')
          squarem_structure <- c('list', 'list', 'vector', 'vector', 'vector',
                                 'vector', 'vector')
          
        }
        
        if (family == 'negbin'){
          
          stop('Setup SQUAREM For negbin')
          if (vi_r_method == 'VEM'){
            
            squarem_par <- c(squarem_par, 'vi_r_mu')
            squarem_type <- c(squarem_type, 'real')
            squarem_structure <- c(squarem_structure, 'vector')
            
          } else if (vi_r_method %in% c('Laplace', 'delta')) {
            
            stop('Set up Laplace/delta for SQUAREM')
            
            squarem_par <- c(squarem_par, 'vi_r_mu', 'vi_r_sigma')
            squarem_type <- c(squarem_type, 'real', 'positive')
            squarem_structure <- c(squarem_structure, 'vector', 'vector')
            
          } else if (vi_r_method == 'fixed') {
            
          }
        }
        if (family %in% 'linear'){
          stop('Set up SQUAREM for linear model.')
        }
        
        if (!do_huangwand){
          remove_hw_b <- FALSE
          squarem_par <- squarem_par[-1]
          squarem_type <- squarem_type[-1]
          squarem_structure <- squarem_structure[-1]
        }else{
          remove_hw_b <- FALSE
          if (remove_hw_b){
            squarem_type <- squarem_type[!grepl(squarem_par, pattern='vi_a_b_jp')]
            squarem_structure <- squarem_structure[!grepl(squarem_par, pattern='vi_a_b_jp')]
            squarem_par <- squarem_par[!grepl(squarem_par, pattern='vi_a_b_jp')]
          }
        }
        
        remove_c <- FALSE
        if (remove_c){
          squarem_type <- squarem_type[!grepl(squarem_par, pattern='vi_pg_c')]
          squarem_structure <- squarem_structure[!grepl(squarem_par, pattern='vi_pg_c')]
          squarem_par <- squarem_par[!grepl(squarem_par, pattern='vi_pg_c')]
        }
        
        check_tri <- sapply(squarem_par[squarem_type == 'cholesky'], FUN=function(nm_i){
          si <- sapply(squarem_list, FUN=function(i){isTriangular(i[[nm_i]])})
          return(all(si))
        })
        squarem_type[squarem_type == 'cholesky'][check_tri == FALSE] <- 'lu'
        
        if ('vi_pg_c' %in% squarem_par){
          # Address possibility of "zero" for vi_pg_c
          squarem_list <- lapply(squarem_list, FUN=function(i){
            i$vi_pg_c <- ifelse(abs(i$vi_pg_c) < 1e-6, 1e-6, i$vi_pg_c)
            return(i)
          })
        }
        
        squarem_list <- lapply(squarem_list, FUN=function(i){
          
          i[squarem_par] <- mapply(squarem_par, squarem_type, 
            squarem_structure, SIMPLIFY = FALSE, FUN=function(s_par, s_type, s_str){
              if (s_str == 'vector'){
                out <- squarem_prep_function(i[[s_par]], s_type) 
              }else{
                out <- lapply(i[[s_par]], FUN=function(j){squarem_prep_function(j, s_type)})
              }
              return(out)
          })
          return(i)
        })
        
        prep_SQUAREM <- mapply(squarem_par, squarem_structure, squarem_type,
          SIMPLIFY = FALSE, FUN=function(s_par, s_str, s_type){
            
            if (s_type == 'lu'){
              
              r <- list(
                'L' = squarem_list[[2]][[s_par]]$L - squarem_list[[1]][[s_par]]$L,
                'U' = squarem_list[[2]][[s_par]]$U - squarem_list[[1]][[s_par]]$U
              )
              d2 <- list(
                'L' = squarem_list[[3]][[s_par]]$L - squarem_list[[2]][[s_par]]$L,
                'U' = squarem_list[[3]][[s_par]]$U - squarem_list[[2]][[s_par]]$U
              )
              v <- list("L" = d2$L - r$L, 'U' = d2$U - r$U)
              norm_sq_r <- sum(sapply(r, FUN=function(i){sum(i@x^2)}))
              norm_sq_v <- sum(sapply(v, FUN=function(i){sum(i@x^2)}))
              max_d <- max(sapply(d2, FUN=function(i){max(abs(i@x))}))
              P <- squarem_list[[3]][[s_par]]$P
              Q <- squarem_list[[3]][[s_par]]$Q
            }else if (s_str == 'list'){
              r <- mapply(squarem_list[[2]][[s_par]], squarem_list[[1]][[s_par]], SIMPLIFY = FALSE, FUN=function(i,j){i - j})
              d2 <- mapply(squarem_list[[3]][[s_par]], squarem_list[[2]][[s_par]], SIMPLIFY = FALSE, FUN=function(i,j){i - j})
              v <- mapply(d2, r, SIMPLIFY = FALSE, FUN=function(i,j){i - j})
              norm_sq_r <- sum(unlist(lapply(r, as.vector))^2)
              norm_sq_v <- sum(unlist(lapply(v, as.vector))^2)
              max_d <- max(abs(sapply(d2, FUN=function(j){max(abs(j))})))
              P <- NULL
              Q <- NULL
            }else{
              r <- squarem_list[[2]][[s_par]] - squarem_list[[1]][[s_par]]
              d2 <- squarem_list[[3]][[s_par]] - squarem_list[[2]][[s_par]]
              v <- d2 - r
              norm_sq_r <- sum(r^2)
              norm_sq_v <- sum(v^2)
              max_d = max(abs(d2))
              P <- NULL
              Q <- NULL
            }
          return(list(first = squarem_list[[1]][[s_par]], 
                      second = squarem_list[[2]][[s_par]], 
                      max_d = max_d, P = P, Q = Q,
                      r = r, v = v, norm_sq_r = norm_sq_r, norm_sq_v = norm_sq_v))
        })
        
        ind_alpha <- FALSE
        
        if (ind_alpha){
          alpha <- -sqrt((sapply(prep_SQUAREM, FUN=function(i){i$norm_sq_r}))) /
            sqrt((sapply(prep_SQUAREM, FUN=function(i){i$norm_sq_v})))
          if (any(alpha > -1)){
            alpha[which(alpha > -1)] <- -1.01
          }
          if (any(alpha < -10)){
            alpha[which(alpha < -10)] <- -10
          }
          max_d <- sapply(prep_SQUAREM, FUN=function(i){i$max_d})
          if (any(max_d < tolerance_parameters)){
            alpha[which(max_d < tolerance_parameters)] <- -1.01
          }
        }else{
          
          alpha <- -sqrt(sum(sapply(prep_SQUAREM, FUN=function(i){i$norm_sq_r}))) /
            sqrt(sum(sapply(prep_SQUAREM, FUN=function(i){i$norm_sq_v})))
          
          if (alpha > -1){
            alpha <- -1.01
          }
          if (alpha < -10){
            alpha <- -10
          }
          
          alpha <- rep(alpha, length(prep_SQUAREM))
          names(alpha) <- names(prep_SQUAREM)
        }
        if (!quiet_rho){print(alpha)}
        
        orig_squarempar <- squarem_par
        orig_alpha <- alpha
        for (attempt_SQUAREM in 1:10){
          
          if (!quiet_rho){print(mean(alpha))}
          
          squarem_par <- orig_squarempar
          if (attempt_SQUAREM > 1){
            alpha <- (alpha - 1)/2
          }
          
          prop_squarem <- mapply(prep_SQUAREM, squarem_structure, squarem_type, alpha, SIMPLIFY = FALSE,
             FUN=function(i, s_str, s_type, s_alpha){
               if (s_type == 'lu'){
                 prop_squarem <- lapply(c('L', 'U'), FUN=function(k){
                   i$first[[k]] - 2 * s_alpha * i$r[[k]] + s_alpha^2 * i$v[[k]] 
                 })
                 names(prop_squarem) <- c('L', 'U')
                 prop_squarem$P <- i$P
                 prop_squarem$Q <- i$Q
                 if (!quiet_rho){if (!isTRUE(all.equal(i$second$P, i$P))){print('MISALIGNED at P')}}
                 if (!quiet_rho){if (!isTRUE(all.equal(i$second$Q, i$Q))){print('MISALIGNED at Q')}}
               }else if (s_str == 'list'){
                 prop_squarem <- mapply(i$first, i$second, 
                    i$r, i$v, SIMPLIFY = FALSE, FUN=function(i_1, s_1, r_1, v_1){
                      out <- i_1 - 2 * s_alpha * r_1 + s_alpha^2 * v_1
                      return(out)
                    })
                 names(prop_squarem) <- names(i$first)
               }else{
                 prop_squarem <- i$first - 2 * s_alpha * i$r + s_alpha^2 * i$v
               }
               return(prop_squarem)
             })
          
          names(prop_squarem) <- squarem_par
          
          prop_ELBOargs <- ELBOargs
          
          prop_squarem <- mapply(prop_squarem, squarem_type, 
           squarem_structure, SIMPLIFY = FALSE, FUN=function(i, s_type, s_str){
             if (s_str == 'vector'){
               out <- squarem_unprep_function(i, s_type) 
             }else{
               out <- lapply(i, FUN=function(j){squarem_unprep_function(j, s_type)})
             }
             return(out)
           })
          
          if (factorization_method == 'weak'){
            
            if (squarem_type[squarem_par == 'vi_joint_L_nonpermute'] == 'lu'){
              prop_squarem$vi_joint_decomp <- prop_squarem$vi_joint_L_nonpermute$M              
              prop_ELBOargs$log_det_joint_var <- prop_squarem$vi_joint_L_nonpermute$logdet_M
            }else{
              prop_squarem$vi_joint_decomp <- prop_squarem$vi_joint_L_nonpermute %*% t(squarem_list[[1]]$vi_joint_LP)
              prop_ELBOargs$log_det_joint_var <- 2 * sum(log(diag(prop_squarem$vi_joint_L_nonpermute)))
            }
            prop_squarem$vi_alpha_decomp <- prop_squarem$vi_joint_decomp[, -1:-p.X, drop = F]
            prop_squarem$vi_beta_decomp <- prop_squarem$vi_joint_decomp[, 1:p.X, drop = F]
            
            squarem_par <- c(squarem_par, 'log_det_joint_var')
            squarem_par <- c(squarem_par, 'vi_joint_decomp')
            
          }else if (factorization_method == 'collapsed'){
            stop('Setup squarem for collapsed')
            if (squarem_type[squarem_par == 'vi_joint_L_nonpermute'] == 'lu'){
              prop_squarem$vi_joint_decomp <- prop_squarem$vi_joint_L_nonpermute$M              
              prop_ELBOargs$log_det_joint_var <- prop_squarem$vi_joint_L_nonpermute$logdet_M
            }else{
              prop_squarem$vi_joint_decomp <- prop_squarem$vi_joint_L_nonpermute %*% t(squarem_list[[1]]$vi_joint_LP)
              prop_ELBOargs$log_det_joint_var <- 2 * sum(log(diag(prop_squarem$vi_joint_L_nonpermute)))
            }
            prop_squarem$vi_alpha_decomp <- prop_squarem$vi_joint_decomp[, -1:-p.X, drop = F]
            prop_squarem$vi_beta_decomp <- prop_squarem$vi_joint_decomp[, 1:p.X, drop = F]
            
            squarem_par <- c(squarem_par, 'log_det_joint_var')
            squarem_par <- c(squarem_par, 'vi_joint_decomp')
            
            
          }else{

            if (squarem_type[squarem_par == 'vi_beta_L_nonpermute'] == 'lu'){
              prop_ELBOargs$log_det_beta_var <- prop_squarem$vi_beta_L_nonpermute$logdet_M
              prop_squarem$vi_beta_decomp <- prop_squarem$vi_beta_L_nonpermute$M              
            }else{
              prop_ELBOargs$log_det_beta_var <- 2 * sum(log(diag(prop_squarem$vi_beta_L_nonpermute)))
              prop_squarem$vi_beta_decomp <- prop_squarem$vi_beta_L_nonpermute %*% t(squarem_list[[1]]$vi_beta_LP)
            }

            if (squarem_type[squarem_par == 'vi_alpha_L_nonpermute'] == 'lu'){
              prop_ELBOargs$log_det_alpha_var <- prop_squarem$vi_alpha_L_nonpermute$logdet_M
              prop_squarem$vi_alpha_decomp <- prop_squarem$vi_alpha_L_nonpermute$M              
            }else{
              prop_ELBOargs$log_det_alpha_var <- 2 * sum(log(diag(prop_squarem$vi_alpha_L_nonpermute)))
              prop_squarem$vi_alpha_decomp <- prop_squarem$vi_alpha_L_nonpermute %*% t(squarem_list[[1]]$vi_alpha_LP)
            }
            
            squarem_par <- c(squarem_par, 'log_det_alpha_var', 'log_det_beta_var')
            squarem_par <- c(squarem_par, 'vi_alpha_decomp', 'vi_beta_decomp')
            
          }
          
          if (family == 'negbin'){
            
            if (vi_r_method == 'VEM'){
              
              prop_ELBOargs$vi_r_mean <- exp(prop_squarem$vi_r_mu)

            } else if (vi_r_method %in% c('Laplace', 'delta')){
              
              prop_ELBOargs$vi_r_mean <- exp(prop_squarem$vi_r_mu + prop_squarem$vi_r_sigma/2)
            
            } 
            
            if (factorization_method != 'weak'){
              prop_joint_var <- rowSums((X %*% t(prop_squarem$vi_beta_decomp))^2) + 
                rowSums((Z %*% t(prop_squarem$vi_alpha_decomp))^2)
            }else{
              prop_joint_var <-  cpp_zVz(Z = joint.XZ, V = as(prop_squarem$vi_joint_decomp, "dgCMatrix")) 
            }
            
            if (vi_r_method %in% c('Laplace', 'delta')){
              prop_joint_var <- prop_joint_var + vi_r_sigma
            }

            prop_ELBOargs$vi_pg_c <- sqrt(as.vector(X %*% prop_squarem$vi_beta_mean + Z %*% prop_squarem$vi_alpha_mean - prop_squarem$vi_r_mu)^2 + prop_joint_var)
            prop_ELBOargs$vi_pg_b <- y + prop_ELBOargs$vi_r_mean
            
          }
          
          for (v in names(prop_squarem)){
            prop_ELBOargs[[v]] <- prop_squarem[[v]]
          }
          
          prop_ELBOargs$vi_alpha_L_nonpermute <- NULL
          prop_ELBOargs$vi_beta_L_nonpermute <- NULL
          prop_ELBOargs$vi_joint_L_nonpermute <- NULL
          
          prop_squarem$vi_collapsed_P <- matrix(0, 1, 1)
          prop_squarem$vi_beta_decomp <- as.matrix(prop_squarem$vi_beta_decomp)  

          if ('vi_alpha_mean' %in% names(prop_squarem)){
            prop_variance_by_alpha_jg <- calculate_expected_outer_alpha(
              factorization_method = factorization_method,
              alpha_decomp_var = (prop_squarem$vi_alpha_decomp), 
              alpha_var = prop_squarem$vi_alpha_var,
              alpha_mu = as.vector(prop_squarem$vi_alpha_mean), 
              re_position_list = outer_alpha_RE_positions,
              L_beta = prop_squarem$vi_beta_decomp, tP = t(prop_squarem$vi_collapsed_P), 
              do_adjustment = adjust_ceoa)
            prop_ELBOargs[['vi_sigma_outer_alpha']] <- prop_variance_by_alpha_jg$outer_alpha
            squarem_par <- c(squarem_par, 'vi_sigma_outer_alpha')
          }
          
          if (remove_hw_b){
            prop_diag_Einv_sigma <- mapply(prop_ELBOargs$vi_sigma_alpha, 
                                           vi_sigma_alpha_nu, d_j, SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
                                             inv_phi <- solve(phi)
                                             sigma.inv <- nu * inv_phi
                                             return(diag(sigma.inv))
                                           })
            prop_ELBOargs$vi_a_b_jp <- mapply(vi_a_nu_jp, vi_a_APRIOR_jp, prop_diag_Einv_sigma,
                                              SIMPLIFY = FALSE,
                                              FUN=function(nu, APRIOR, diag_j){
                                                1/APRIOR^2 + nu * diag_j
                                              })
            squarem_par <- c(squarem_par, 'vi_a_b_jp')
          }
          if ('vi_pg_c' %in% squarem_par){
            
            if (family != 'binomial'){stop('check squarem for linear case')}
            
            prop_vi_pg_mean <- prop_ELBOargs$vi_pg_b / (2 * prop_ELBOargs$vi_pg_c) * tanh(prop_ELBOargs$vi_pg_c / 2)
            fill_zero <- which(abs(prop_ELBOargs$vi_pg_c) < 1e-6)
            if (length(fill_zero) > 0){
              prop_vi_pg_mean[fill_zero] <- prop_ELBOargs$vi_pg_b[fill_zero]/4
            }
            prop_ELBOargs[['vi_pg_mean']] <- prop_vi_pg_mean
            squarem_par <- c(squarem_par, 'vi_pg_mean')
            
          }else{
            
            if (family != 'binomial'){stop('check squarem for non-binomial case')}
            
            if (factorization_method %in% c("weak")) {
              joint_quad <- cpp_zVz(Z = joint.XZ, V = as(prop_ELBOargs$vi_joint_decomp, "dgCMatrix")) 
              if (family == 'negbin'){
                joint_quad <- joint_quad + prop_ELBOargs$vi_r_sigma
              }
              prop_ELBOargs$vi_pg_c <- sqrt(as.vector(X %*% prop_ELBOargs$vi_beta_mean + Z %*% prop_ELBOargs$vi_alpha_mean - prop_ELBOargs$vi_r_mu)^2 + joint_quad)
            } else if (grepl(factorization_method, pattern='collapsed')) {
              stop('Setup Squarem for "collapsed"')
            } else {
              beta_quad <- rowSums((X %*% t(prop_ELBOargs$vi_beta_decomp))^2)
              alpha_quad <- rowSums((Z %*% t(prop_ELBOargs$vi_alpha_decomp))^2)
              joint_var <- beta_quad + alpha_quad
              if (family == 'negbin'){
                joint_var <- joint_var + prop_ELBOargs$vi_r_sigma
              }
              prop_ELBOargs$vi_pg_c <- sqrt(as.vector(X %*% prop_ELBOargs$vi_beta_mean + Z %*% prop_ELBOargs$vi_alpha_mean - prop_ELBOargs$vi_r_mu)^2 + joint_var)
            }
            
            prop_vi_pg_mean <- prop_ELBOargs$vi_pg_b / (2 * prop_ELBOargs$vi_pg_c) * tanh(prop_ELBOargs$vi_pg_c / 2)
            fill_zero <- which(abs(prop_ELBOargs$vi_pg_c) < 1e-6)
            if (length(fill_zero) > 0){
              prop_vi_pg_mean[fill_zero] <- prop_ELBOargs$vi_pg_b[fill_zero]/4
            }
            
            prop_ELBOargs[['vi_pg_mean']] <- prop_vi_pg_mean
            squarem_par <- c(squarem_par, 'vi_pg_c', 'vi_pg_mean')
            
          }
          
          elbo_init <- do.call("calculate_ELBO", ELBOargs)
          elbo_squarem <- do.call("calculate_ELBO", prop_ELBOargs)
          if (!quiet_rho){print(c(elbo_squarem$ELBO, elbo_init$ELBO))}
          if (elbo_squarem$ELBO >= elbo_init$ELBO){break}
            
        }
        
        # squarem_par <<- squarem_par
        # prop_ELBOargs <<- prop_ELBOargs

        if (elbo_squarem$ELBO >= elbo_init$ELBO){
          if (!quiet_rho){cat('SUCCESS')}
          squarem_success <- squarem_success + 1
          squarem.ELBO <- elbo_squarem
          final.ELBO <- elbo_squarem
          
          for (v in squarem_par){
            assign(v, prop_ELBOargs[[v]])
          }
          
          test_ELBO <- calculate_ELBO(family = family,
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
               vi_alpha_var = vi_alpha_var, vi_beta_var = vi_beta_var, cyclical_pos = cyclical_pos,
               vi_beta_decomp = vi_beta_decomp, vi_alpha_decomp = vi_alpha_decomp,
               vi_joint_decomp = vi_joint_decomp, choose_term = choose_term,
               vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
               vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
               log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma,
               do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
               vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
          )
          if (test_ELBO$ELBO != elbo_squarem$ELBO){stop('....')}
        }else{
          if (!quiet_rho){cat('FAIL')}
          squarem_success[1] <- squarem_success[1] + 1
          final.ELBO <- squarem.ELBO <- final.ELBO
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
      if (do_SQUAREM & (it %% 3 == 0)){
        squarem.ELBO$step <- 4
        final.ELBO$step <- 5
        update_ELBO <- rbind(debug_ELBO.1, debug_ELBO.2, debug_ELBO.3, squarem.ELBO, final.ELBO)
      }else{
        final.ELBO$step <- 4
        update_ELBO <- rbind(debug_ELBO.1, debug_ELBO.2, debug_ELBO.3, final.ELBO)
      }
      update_ELBO$it <- it
      store_ELBO <- rbind(store_ELBO, update_ELBO)
    } else {
      final.ELBO$it <- it
      store_ELBO <- rbind(store_ELBO, final.ELBO)
    }
    

    # if (!quiet_rho){
    #   if (factorization_method == 'weak'){
    #     print('NonsparseA')
    #     print(length(vi_joint_decomp@x))
    #   }else{
    #     print('NonsparseA')
    #     print(length(vi_alpha_decomp@x))
    #   }
    # }

    ## Change diagnostics
    
    change_elbo <- final.ELBO$ELBO - lagged_ELBO

    if (factorization_method == "collapsed"){
      if (any_collapsed_M){
        change_alpha_mean <- max(abs(do.call('c', vi_M_mean) - lagged_M_mean))
      }else{
        change_alpha_mean <- 0
      }
      if (any_collapsed_C){
        change_beta_mean <- max(abs(vi_C_mean - lagged_C_mean))
      }else{
        change_beta_mean <- 0
      }
    }else{
      change_alpha_mean <- max(abs(vi_alpha_mean - lagged_alpha_mean))
      change_beta_mean <- max(abs(vi_beta_mean - lagged_beta_mean))
    }
    
    unlist_vi <- c(unlist(lapply(vi_sigma_alpha, as.vector)), unlist(vi_a_b_jp))
    if (debug_ELBO){
      svi <- data.frame(t(as.vector(unlist_vi)))
      svi$it <- it
      store_vi <- rbind(store_vi,  svi)
    }
    
    change_sigma_mean <- mapply(vi_sigma_alpha, lagged_sigma_alpha, FUN = function(i, j) {
      max(abs(i - j))
    })

    if (factorization_method == "weak") {
      change_joint_var <- 0 # change_joint_var <- max(abs(vi_joint_decomp - lagged_joint_decomp))
      change_alpha_var <- change_beta_var <- 0
    } else if (factorization_method == "collapsed") {
      
      if (any_collapsed_M){
        change_alpha_var <- max(mapply(vi_M_var, lagged_M_var, FUN=function(i,j){if(ncol(i) > 0){max(abs(i-j))}else{0}}))
      }else{
        change_alpha_var <- 0
      }
      if (any_collapsed_C){
        change_beta_var <- max(abs(vi_C_var - lagged_C_var))
      }else{change_beta_var <- 0}
      change_joint_var <- 0
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
      if (factorization_method == "collapsed"){
        store_beta[it, ] <- as.vector(vi_C_mean)
        store_alpha[it, ] <- unlist(vi_M_mean)
        store_sigma[it,] <- do.call('c', lapply(vi_sigma_alpha, as.vector))
      }else{
        store_beta[it, ] <- as.vector(vi_beta_mean)
        store_alpha[it, ] <- as.vector(vi_alpha_mean)
        store_sigma[it,] <- do.call('c', lapply(vi_sigma_alpha, as.vector))
      }
      if (do_huangwand){
        store_hw[it,] <- do.call('c', vi_a_b_jp)
      }
    }
    
    change_all <- data.frame(change_alpha_mean, change_beta_mean, 
        t(change_sigma_mean), change_alpha_var, change_beta_var, change_joint_var, change_vi_r_mu)

    if ((max(change_all) < tolerance_parameters) | (change_elbo > 0 & change_elbo < tolerance_elbo)) {
      if (!quiet) {
        message(paste0("Converged after ", it, " iterations with ELBO change of ", round(change_elbo, 1 + abs(floor(log(tolerance_elbo) / log(10))))))
        message(paste0("The largest change in any variational parameter was ", round(max(change_all), 1 + abs(floor(log(tolerance_parameters) / log(10))))))
      }
      break
    }
    if (debug_ELBO){
      change_all$it <- it
      store_parameter_traj <- rbind(store_parameter_traj, change_all)
    }
    if (!quiet & (it %% print_prog == 0)) {
      message(paste0("ELBO Change: ", round(change_elbo, 10)))
      message(paste0("Other Parameter Changes: ", max(change_all)))
    }

    if (factorization_method != "collapsed"){
      lagged_alpha_mean <- vi_alpha_mean
      lagged_beta_mean <- vi_beta_mean
      lagged_alpha_decomp <- vi_alpha_decomp
      lagged_beta_decomp <- vi_beta_decomp
      
      lagged_alpha_var <- vi_alpha_var
      lagged_beta_var <- vi_beta_var
      
      lagged_sigma_alpha <- vi_sigma_alpha
      lagged_vi_r_mu <- vi_r_mu
      
      lagged_ELBO <- final.ELBO$ELBO
    }else{

      lagged_C_mean <- vi_C_mean
      lagged_M_mean <- do.call('c', vi_M_mean)
      lagged_C_var <- vi_C_var
      lagged_M_var <- vi_M_var
      
      lagged_sigma_alpha <- vi_sigma_alpha
      lagged_vi_r_mu <- vi_r_mu
      
      lagged_ELBO <- final.ELBO$ELBO
      
    }
  }
  if (it == iterations) {
    message(paste0("Ended without Convergence after", it, " iterations : ELBO change of ", round(change_elbo[1], abs(floor(log(tolerance_elbo) / log(10))))))
  }

  if (parameter_expansion == "translation") {
    final.ELBO$accepted_PX <- accepted_times / attempted_expansion
  }
  
  output <- list(
    beta = list(mean = vi_beta_mean),
    ELBO = final.ELBO, 
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

    tic_summary <- lapply(split(tic_log$time, tic_log$stage),
           FUN=function(i){
             data.frame(n = length(i), mean = mean(i), min = min(i), max = max(i),
                        total = sum(i))
           }
    )
    tic_summary <- do.call('rbind', tic_summary)
    tic_summary$variable <- rownames(tic_summary)
    rownames(tic_summary) <- NULL
  } else {
    tic_summary <- NULL
  }
  if (debug_param) {
    
    store_beta <- store_beta[1:it,,drop=F]
    store_alpha <- store_alpha[1:it,,drop=F]
    store_sigma <- store_sigma[1:it,,drop=F]
    if (do_huangwand){store_hw <- store_hw[1:it,,drop=F]}else{store_hw <- NULL}
    output$parameter_trajectory <- list(beta = store_beta,
                                        alpha = store_alpha,
                                        sigma = store_sigma,
                                        hw = store_hw)
  }
  if (factorization_method == "weak") {
    output$joint <- vi_joint_decomp
  }else if (factorization_method == "collapsed") {

    vi_C_var <- solve(drop0(
      t(design_C) %*% diag_vi_pg_mean %*% design_C + Tinv_C
    ))
    C_hat <- vi_C_var %*% t(design_C) %*% s
    output$collapsed$cond_variance <- vi_C_var
    output$collapsed$C.hat <- C_hat
    output$collapsed$index <- C_j
    vi_P <- lapply(vi_M_list, FUN=function(i){
      vi_C_var %*% 
        t(design_C) %*% Diagonal(x = vi_pg_mean) %*% i
    })
    vi_P <- do.call('cbind', vi_P)

    if (any_collapsed_C){
      vi_C_uncond <- vi_C_var + vi_P %*% bdiag(vi_M_var) %*% t(vi_P)
      output$collapsed$variance <- vi_C_uncond
      output$collapsed$mean <- vi_C_mean
      output$collapsed$names_of_collapsed <- names_of_collapsed
    }else{
      output$collapsed <- list(variance = NULL, mean = NULL, decomp = NULL)
    }
    output$marginal <- list(variance = vi_M_var, mean = vi_M_mean, index = M_j)
    
    if (block_collapse){

      vi_C_mean <- vi_M_mean[[position_of_blocked]]
      vi_M_mean <- vi_M_mean[-position_of_blocked]
      
      vi_C_var <- vi_M_var[[position_of_blocked]]
      vi_M_var <- vi_M_var[-position_of_blocked]
      
      C_j <- orig_Cj
      M_j <- M_j[names(M_j) != '...blocked']
      
    }
    
    fmt_mean <- split(as.vector(vi_C_mean), rep(seq_len(length(C_j)), lengths(C_j)))
    names(fmt_mean) <- names(C_j[lengths(C_j) != 0])
    fmt_mean <- fmt_mean[names(C_j)]
    names(fmt_mean) <- names(C_j)
    fmt_mean <- mapply(fmt_mean, vi_M_mean, C_j, M_j, SIMPLIFY = FALSE,
        FUN=function(ci,mi,index_c, index_m){
          c(ci, mi)[order(c(index_c, index_m))]
    })

    var_M <- lapply(vi_M_var, diag)
    if (block_collapse){
      var_C <- split(diag(vi_C_var), position_block_j)
    }else{
      var_C <- split(diag(vi_C_uncond), position_block_j)
    }
    names(var_C) <- names(C_j[lengths(C_j) != 0])
    var_C <- var_C[names(C_j)]
    names(var_C) <- names(C_j)
    
    fmt_var <- mapply(var_M, var_C, C_j, M_j, SIMPLIFY = FALSE,
       FUN=function(ci,mi,index_c, index_m){
         c(ci, mi)[order(c(index_c, index_m))]
       })
    
    output$beta$mean <- matrix(fmt_mean[[1]])
    output$beta$var <- Diagonal(x = fmt_var[[1]])
    output$alpha$mean <- matrix(do.call('c', fmt_mean[-1]))
    output$alpha$dia.var <- unlist(fmt_var[-1])
    rownames(output$beta$mean) <- colnames(X)
    rownames(output$alpha$mean) <- colnames(Z)
    rownames(output$alpha$dia.var) <- NULL
  }
  
  if (control$return_data) {
    output$data <- list(X = X, Z = Z, y = y, trials = trials)
  }
  
  output$spline <- list(attr = Z.spline.attr, size = Z.spline.size)

  output$formula <- list(formula = formula, 
     re = re_fmla, fe = fe_fmla,
     interpret_gam = parse_formula)
  
  
  if (factorization_method != "collapsed"){
    output$alpha$dia.var <- unlist(lapply(variance_by_alpha_jg$variance_jg, FUN = function(i) {
      as.vector(sapply(i, diag))
    }))
    
    output$beta$var <- t(vi_beta_decomp) %*% vi_beta_decomp
    output$beta$decomp_var <- vi_beta_decomp
    
    rownames(output$beta$mean) <- rownames(output$beta$var) <- colnames(X)
    colnames(output$beta$var) <- colnames(X)
    
    output$alpha$var <- variance_by_alpha_jg$variance_jg
    output$alpha$decomp_var <- vi_alpha_decomp
    
  }
  if (family == "negbin") {
    output$r <- list(mu = vi_r_mu, sigma = vi_r_sigma, method = vi_r_method)
  }
  if (do_huangwand){
    output$hw <- list(a = vi_a_a_jp, b = vi_a_b_jp)
  }

  output$internal_parameters <- list(
    vi_collapsed_P = vi_collapsed_P,
    it_used = it, it_max = iterations,
    missing_obs = missing_obs, N = nrow(X),
    acceleration = list(accept.PX = accept.PX, 
      squarem_success = squarem_success),
    names_of_RE = names_of_RE, d_j = d_j, g_j = g_j
  )

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
#' @param do_SQUAREM Accelerate method using SQUAREM
#' 
#' @importFrom checkmate assert check_double check_logical check_choice check_int check_integerish
#' @export
vglmer_control <- function(iterations = 1000, 
   collapse_size = "FE", block_collapse = FALSE,
   prior_variance = "mean_exists", factorization_method = "weak",
   tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
   prevent_degeneracy = FALSE, force_whole = TRUE, verbose_time = TRUE,
   parameter_expansion = "mean", random_seed = 1, do_timing = FALSE,
   debug_param = FALSE, return_data = FALSE, linpred_method = "joint",
   vi_r_method = "VEM", vi_r_val = NA, do_SQUAREM = TRUE, verify_columns = FALSE,
   debug_ELBO = FALSE, print_prog = NULL, quiet = T, quiet_rho = TRUE,
   init = "EM_FE") {
  
  # use checkmate package to verify arguments
  assert(
    check_integerish(iterations, lower = 1),
    check_double(c(tolerance_elbo, tolerance_parameters), len = 2, lower = 0),
    check_logical(c(
      prevent_degeneracy, force_whole, verbose_time, do_timing,
      debug_param, return_data, debug_ELBO, quiet
    ), len = 8),
    check_choice(factorization_method, c("weak", "strong", "collapsed", "partial")),
    check_choice(prior_variance, c("kn", "hw", "mean_exists", "jeffreys", "mcmcglmm", "mvD", "limit", "uniform")),
    check_choice(linpred_method, c("joint", "cyclical", "solve_normal")),
    check_choice(vi_r_method, c("VEM", "fixed", "Laplace", "delta")),
    check_double(vi_r_val, all.missing = TRUE),
    check_int(print_prog, null.ok = TRUE),
    check_choice(init, c("EM", "random", "zero", "EM_FE")),
    check_double(random_seed),
    combine = "and"
  )

  if (vi_r_method == "fixed" & is.na(vi_r_val)) {
    stop('vi_r_val must not be NA if vi_r_method = "fixed"')
  }

  output <- namedList(
    iterations, prior_variance, factorization_method,
    tolerance_elbo, tolerance_parameters, quiet_rho,verify_columns,
    prevent_degeneracy, force_whole, verbose_time, do_SQUAREM,
    parameter_expansion, random_seed, do_timing, debug_param, return_data,
    linpred_method, vi_r_method, vi_r_val, debug_ELBO, print_prog, quiet, init,
    collapse_size, block_collapse
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

calculate_expected_outer_alpha <- function(alpha_decomp_var, alpha_var, alpha_mu, re_position_list, 
      tP, L_beta, do_adjustment, factorization_method){
  
  if (factorization_method != 'collapsed'){
    out <- decomp_calculate_expected_outer_alpha(L = alpha_decomp_var, alpha_mu = alpha_mu,
      re_position_list = re_position_list,
      tP = tP, L_beta = L_beta, do_adjustment = do_adjustment)
  }else{
    out <- direct_calculate_expected_outer_alpha(V = alpha_var, alpha_mu = alpha_mu,
                                                 re_position_list = re_position_list)
  }

  return(out)
}
