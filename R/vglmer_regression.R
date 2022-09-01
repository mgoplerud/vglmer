#' Variational Inference for Non-Linear Hierarchical Models
#'
#' Estimate hierarchicals model using mean-field variational inference. Accepts
#' standard syntax used for \code{lme4}, e.g., \code{y ~ X + (1 + Z | g)}. Options are
#' described below. Goplerud (2022a; 2022b) provides details on the variational
#' algorithm.
#'
#' @param formula \code{lme4} formula for random effects. Options involving
#'   \code{||} have not been tested. Typically, \code{(1 + Z | G)} indicates a
#'   random effect for each level of variable \code{"G"} with a differing slope
#'   for the effect of variable \code{"Z"} and an intercept (\code{1}). See,
#'   e.g., Gelman and Hill (2006) for a discussion of these models. Splines can
#'   be estimated as described in the "Details" section.
#' @param data data.frame containing the outcome and variables.
#' @param family Options are "binomial", "negbin", or "linear". If "binomial",
#'   outcome must be either {0,1} (binary) or cbind(success, failure) as per
#'   standard glm(er) syntax. Non-integer values are permitted for binomial if
#'   \code{force_whole} is set to FALSE in vglmer_control.
#' @param control Adjust internal options for estimation. Must use an object
#'   created by \link{vglmer_control}.
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
#' # Run with weaker (i.e. better) approximation
#' \donttest{
#' vglmer(y ~ x + (x | g),
#'   data = sim_data,
#'   control = vglmer_control(factorization_method = "weak"),
#'   family = "binomial"
#' )
#' }
#' 
#' @details XX
#' 
#' \bold{Estimation Syntax:} The syntax 
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
#'   qlogis optim residuals lm plogis setNames
#' @importFrom graphics plot
#' @importFrom Rcpp sourceCpp
#' 
#' @references
#' Goplerud, Max. 2022. "Fast and Accurate Estimation of Non-Nested Binomial
#' Hierarchical Models Using Variational Inference." Bayesian Analysis. 17(2):
#' 623-650.
#' @useDynLib vglmer
#' @export
vglmer <- function(formula, data, family, control = vglmer_control()) {

  # Verify integrity of parameter arguments
  family <- match.arg(family, choices = c("negbin", "binomial", "linear"))
  
  checkdf <- inherits(data, 'data.frame')
  if (is.null(data)){
    checkdf <- TRUE
  }
  if (checkdf != TRUE) {
    warning(paste0("data is not a data.frame? Behavior may be unexpected: ", checkdf))
  }
  if (!inherits(formula, 'formula')){
    stop('"formula" must be a formula.')
  }
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
  
  data <- model.frame(parse_formula$fake.formula, data,
                      drop.unused.levels = TRUE)
  
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
  if (!(factorization_method %in% c("weak", "strong", "partial", "collapsed"))) {
    stop("factorization_method must be in weak, strong, partial or collapsed.")
  }
  if (is.null(print_prog)) {
    print_prog <- max(c(1, floor(iterations / 20)))
  }
  if (!(family %in% c("binomial", "negbin", "linear"))) {
    stop('family must be "linear", "binomial", "negbin".')
  }
  
  if (family %in% c('binomial', 'linear')){
    if (control$prior_variance == 'hw' & control$prior_variance %in% c('diagonal', 'translation')){
      message('hw and negative binomial or linear not yet implemented.')
      control$parameter_expansion <- 'mean'
    }
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
    fe_update <- sapply(fe_fmla$smooth.spec, FUN=function(i){
      if (i$by != "NA" & i$by_re == FALSE){
        fe_i <- paste0(i$term, ' * ', i$by)
      }else{
        fe_i <- i$term
      }
    })
    fe_update <- paste0(fe_update, collapse = ' + ')
    
    fe_fmla <- update.formula(fe_fmla$pf,
        paste0('. ~ . + 1 + ', fe_update)
    )
    
  }else{
    fe_fmla <- fe_fmla$pf
  }
  
  # Create the FE design
  X <- model.matrix(fe_fmla, data = data)
  
  # Extract the Z (Random Effect) design matrix.
  re_fmla <- findbars(formula)
  
  # If using splines by group, add random effects to
  # the main level.
  if (!all(sapply(parse_formula$smooth.spec, 
      FUN=function(i){i$by}) %in% c('NA'))){
    
    by_splines <- parse_formula$smooth.spec[
      which(sapply(parse_formula$smooth.spec, FUN=function(i){(i$by != "NA" & i$by_re == TRUE)}))
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
  id_M <- seq_len(ncol(Z))
  # Remove FE
  if (any(is.na(M[,1]))){
    which_is_na_M <- which(is.na(M[,1]))
    id_M <- id_M[-which_is_na_M]
    M <- M[-which_is_na_M, , drop = F]
  }
  
  if (nrow(M) > 0){
    M <- sparseMatrix(i = id_M, j = M[, 1], x = M[, 2], dims = c(ncol(Z), ncol(X)))
  }else{
    M <- drop0(matrix(0, nrow = 0, ncol = ncol(X)))
  }
  
  if (!is.null(names_of_RE)){
    any_Mprime <- TRUE
    M_prime.names <- paste0(rep(names(names_of_RE), g_j * d_j), " @ ", M.names)
    M_prime <- cbind(match(M_prime.names, unique(M_prime.names)), rep(1 / g_j, d_j * g_j))
    M_prime <- sparseMatrix(i = seq_len(ncol(Z))[id_M],
                            j = M_prime[id_M, 1], 
                            x = M_prime[id_M, 2],
                            dims = c(ncol(Z), max(M_prime[,1])))
    colnames(M_prime) <- unique(M_prime.names)
    
    M_prime_one <- M_prime
    M_prime_one@x <- rep(1, length(M_prime_one@x))
    
    stopifnot(identical(paste0(rep(names(names_of_RE), d_j), " @ ", unlist(names_of_RE)), colnames(M_prime)))
    
    mu_to_beta_names <- match(unlist(names_of_RE), colnames(X))
    
    id_mu_to_beta <- seq_len(sum(d_j))
    which_is_na_mu_to_beta <- which(is.na(mu_to_beta_names))
    if (length(which_is_na_mu_to_beta) > 0){
      mu_to_beta_names <- mu_to_beta_names[-which_is_na_mu_to_beta]
      id_mu_to_beta <- id_mu_to_beta[-which_is_na_mu_to_beta]
    }
    
    M_mu_to_beta <- sparseMatrix(
      i = id_mu_to_beta, j = mu_to_beta_names, 
      x = 1, dims = c(sum(d_j), p.X))
    
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
          outer_okay = special_i$outer_okay, by_re = special_i$by_re)
      
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

  debug_px <- control$debug_px
  if (control$parameter_expansion %in% c('translation', 'diagonal')){
    px_method <- control$px_method
    px_it <- control$px_numerical_it
    opt_prior_rho <- NULL
    parsed_RE_groups <- get_RE_groups(formula = re_fmla, data = data)
  }
  
  # List of Lists
  # Outer list: one for RE
  # Inner List: One for each GROUP with its row positions.
  outer_alpha_RE_positions <- mapply(d_j, g_j, breaks_for_RE[-length(breaks_for_RE)], 
    SIMPLIFY = FALSE, FUN = function(a, b, m) {
      split(m + seq(1, a * b), rep(1:b, each = a))
  })

  if (anyDuplicated(unlist(outer_alpha_RE_positions)) != 0 | max(unlist(outer_alpha_RE_positions)) != ncol(Z)) {
    stop("Issue with greating OA positions")
  }
  ####
  # Prepare Initial Values
  ###

  vi_sigmasq_prior_a <- 0
  vi_sigmasq_prior_b <- 0
  
  vi_sigmasq_a <- vi_sigmasq_b <- 1
  
  if (family == "linear") {
    
    vi_sigmasq_a <- (length(y) + sum(d_j * g_j))/2 + vi_sigmasq_prior_a
    vi_sigmasq_b <- 0
    
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
  vi_alpha_decomp <- Diagonal(x = rep(0, ncol(Z)))

  vi_sigma_alpha_nu <- g_j

  prior_variance <- control$prior_variance
  do_huangwand <- FALSE
  vi_a_APRIOR_jp <- vi_a_nu_jp <- vi_a_a_jp <- vi_a_b_jp <- NULL
  prior_sigma_alpha_nu <- prior_sigma_alpha_phi <- NULL
  
  if (prior_variance == 'hw') {
    
    do_huangwand <- TRUE
    INNER_IT <- control$hw_INNER
    vi_a_nu_jp <- rep(2, length(d_j))
    names(vi_a_nu_jp) <- names(names_of_RE)
    vi_a_APRIOR_jp <- lapply(d_j, FUN=function(i){rep(5, i)})
    vi_a_a_jp <- mapply(d_j, vi_a_nu_jp, SIMPLIFY = FALSE, 
                        FUN=function(i,nu){1/2 * (nu + rep(i, i))})
    vi_a_b_jp <- lapply(vi_a_APRIOR_jp, FUN=function(i){1/i^2})
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
  } else if (prior_variance == "gamma") {
    prior_sigma_alpha_nu <- rep(0.001 * 2, length(d_j))
    prior_sigma_alpha_phi <- lapply(d_j, FUN = function(i) {
      diag(x = 0.001 * 2, nrow = i, ncol = i)
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
        EM_init <- LinRegChol(X = drop0(X),
           omega = sparseMatrix(i = 1:nrow(X), j = 1:nrow(X), x = 1),
           y = y, prior_precision = sparseMatrix(i = 1:ncol(X), j = 1:ncol(X), x = 1e-5))$mean
        # stop('Setup EM init for linear')
        # solve(Matrix::Cholesky(  t(joint.XZ) %*% sparseMatrix(i = 1:N, j = 1:N, x = pg_mean) %*% joint.XZ + EM_variance),
        #       t(joint.XZ) %*% (adj_out) )
        EM_init <- list('beta' = EM_init, 'alpha' = rep(0, ncol(Z)))
      }else{
        stop('Setup EM init')
        
        EM_init <- LinRegChol(X = jointXZ, 
                              omega = sparseMatrix(i = 1:nrow(jointXZ), j = 1:nrow(jointXZ), x = 1), 
                              y = y, prior_precision = sparseMatrix(i = 1:ncol(jointXZ), j = 1:ncol(jointXZ), x = 1/4))$mean
        EM_init <- list('beta' = EM_init[1:ncol(X)], 'alpha' = EM_init[-1:-ncol(X)])
      }
      rm(jointXZ)
    } else if (family == "negbin") {
      if (control$init == 'EM_FE'){
        EM_init <- EM_prelim_nb(X = X, Z = drop0(matrix(0, nrow = nrow(X), ncol = 0)), y = y, est_r = exp(vi_r_mu), iter = 15, ridge = 10^5)
        EM_init <- list('beta' = EM_init$beta, 'alpha' = rep(0, ncol(Z)))
      }else{
        EM_init <- EM_prelim_nb(X = X, Z = Z, y = y, est_r = exp(vi_r_mu), iter = 15, ridge = 4)
      }
    } else {
      if (control$init == 'EM_FE'){
        EM_init <- EM_prelim_logit(X = X, Z = drop0(matrix(0, nrow = nrow(X), ncol = 0)), s = s, pg_b = vi_pg_b, iter = 15, ridge = 10^5)
        EM_init <- list('beta' = EM_init$beta, 'alpha' = rep(0, ncol(Z)))
      }else{
        EM_init <- EM_prelim_logit(X = X, Z = Z, s = s, pg_b = vi_pg_b, iter = 15, ridge = 4)
      }
    }

    vi_beta_mean <- matrix(EM_init$beta)
    vi_alpha_mean <- matrix(EM_init$alpha)

    vi_sigma_alpha <- calculate_expected_outer_alpha(
      alpha_mu = vi_alpha_mean,
      L = sparseMatrix(i = 1, j = 1, x = 1e-4, dims = rep(ncol(Z), 2)),
      re_position_list = outer_alpha_RE_positions
    )
    
    if (!do_huangwand){
      vi_sigma_alpha <- mapply(vi_sigma_alpha$outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN = function(i, j) {
        i + j
      })
    }else{
      
      #Update Inverse-Wishart

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
    vi_beta_mean <- rnorm(ncol(X))
    vi_alpha_mean <- rep(0, ncol(Z))

    vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN = function(d, g) {
      
      out <- rWishart(n = 1, df = ifelse(g >= d, g, d), Sigma = diag(d))[ , , 1]
      
      if (d == 1){
        out <- matrix(out)
      }
      
      return(out)
      
    })

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

  if (factorization_method %in% c("weak", "collapsed")) {
    vi_joint_decomp <- bdiag(vi_beta_decomp, vi_alpha_decomp)
    joint.XZ <- cbind(X, Z)
    log_det_beta_var <- log_det_alpha_var <- NULL
  } else {
    vi_joint_decomp <- NULL
    log_det_joint_var <- NULL
  }

  # Create mapping for this to allow sparse implementations.

  mapping_sigma_alpha <- make_mapping_alpha(vi_sigma_alpha)

  running_log_det_alpha_var <- rep(NA, number_of_RE)

  lagged_alpha_mean <- rep(-Inf, ncol(Z))
  lagged_beta_mean <- rep(-Inf, ncol(X))
  lagged_sigma_alpha <- vi_sigma_alpha
  if (factorization_method %in% c("weak", "collapsed")) {
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
  
  accepted_times <- 0
  attempted_expansion <- 0
  
  if (parameter_expansion %in%  c("translation", "diagonal") & any_Mprime) {
    
    if (do_timing){
      tic('Build PX R Terms')
    }
    
    
    spline_REs <- grepl(names(d_j), pattern='^spline-')
    
    nonspline_positions <- sort(unlist(outer_alpha_RE_positions[!spline_REs]))
    
    size_splines <- sum((d_j * g_j)[spline_REs])
    
    stationary_rho <- do.call('c', lapply(d_j[!spline_REs], FUN=function(i){as.vector(diag(x = i))}))
    est_rho <- stationary_rho
    diag_rho <- which(stationary_rho == 1)
    
    zeromat_beta <- drop0(Diagonal(x = rep(0, ncol(X))))

    # parsed_RE_groups <- get_RE_groups(formula = formula, data = data)
    # parsed_RE_groups <- parsed_RE_groups
    
    mapping_new_Z <- do.call('cbind', parsed_RE_groups$design)
    
    mapping_J <- split(1:sum(d_j[!spline_REs]^2), rep(1:length(d_j[!spline_REs]), d_j[!spline_REs]^2))
    mapping_J <- lapply(mapping_J, FUN=function(i){i-1})
    mapping_J <- sapply(mapping_J, min)

    mapping_to_re <- parsed_RE_groups$factor
    mapping_to_re <- unlist(apply(do.call('cbind', mapping_to_re), MARGIN = 1, list), recursive = F)
    # mapping_to_re <- purrr::array_branch(do.call('cbind', mapping_to_re), margin = 1)
    
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
    store_design <- parsed_RE_groups$design
    
    rm(parsed_RE_groups, mapping_to_re)
    
    gc()
    if (do_timing){
      toc(quiet = verbose_time, log = T)
    }
  }
  store_parameter_traj <- store_vi <- store_ELBO <- data.frame()

  if (debug_param) {
    store_beta <- array(NA, dim = c(iterations, ncol(X)))
    store_alpha <- array(NA, dim = c(iterations, ncol(Z)))
    store_sigma <- array(NA, dim = c(iterations, sum(d_j^2)))
    if (do_huangwand){
      store_hw <- array(NA, dim = c(iterations, sum(d_j)))
    }
  }
  if (do_timing) {
    toc(quiet = verbose_time, log = TRUE)
    tic.clear()
  }
  ## Begin VI algorithm:
  if (!quiet) {
    message("Begin Regression")
  }
  do_SQUAREM <- control$do_SQUAREM
  if (factorization_method == 'collapsed'){
    warning('Turning off SQUAREM for "collapsed')
    do_SQUAREM <- FALSE
  }
  if (family %in% c('negbin')){
    if (do_SQUAREM){warning('Turning off SQUAREM for negbin temporarily.')}
    do_SQUAREM <- FALSE
  }
  if (family == 'negbin' & !(control$vi_r_method %in% c('VEM', 'fixed'))){
    if (do_SQUAREM){warning('Turning off SQUAREM if "negbin" and not VEM/fixed.')}
    do_SQUAREM <- FALSE
  }

  if (do_SQUAREM){
    namedList <- utils::getFromNamespace('namedList', 'lme4')
    squarem_success <- c(0, 0)
    squarem_list <- list()
    squarem_counter <- 1
  }else{
    squarem_success <- NA
  }
  if (debug_px){
    debug_PX_ELBO <- rep(NA, iterations)
  }else{
    debug_PX_ELBO <- NULL
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

    if (factorization_method %in% c("weak", "collapsed")) {
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
      fill_zero <- which(abs(vi_pg_c) < 1e-6)
      if (length(fill_zero) > 0){
        vi_pg_mean[fill_zero] <- vi_pg_b[fill_zero] / 4
      }
      diag_vi_pg_mean <- sparseMatrix(i = 1:N, j = 1:N, x = vi_pg_mean)
    }
    sqrt_pg_weights <- Diagonal(x = sqrt(vi_pg_mean))
    
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

    if (factorization_method == "collapsed"){
      cyclical_T <- TRUE
    }
    
    Tinv <- prepare_T(
      mapping = inv_mapping_alpha, levels_per_RE = g_j, num_REs = number_of_RE,
      variables_per_RE = d_j, running_per_RE = breaks_for_RE, cyclical = cyclical_T
    )

    if (!cyclical_T & factorization_method != "collapsed") {
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
      
      if (family != 'binomial'){stop('"collapsed" not set up.')}
      
      beta_var <- solve(t(X) %*% diag_vi_pg_mean %*% X)
      beta_hat <- beta_var %*% t(X) %*% s
      
      P <- beta_var %*% t(X) %*% diag_vi_pg_mean %*% Z
      M <- Z - X %*% P
      
      vi_alpha_mean <- solve(t(M) %*% diag_vi_pg_mean %*% M + bdiag(Tinv),
         t(M) %*% (s - diag_vi_pg_mean %*% X %*% beta_hat)
      )
      vi_beta_mean <- beta_hat - P %*% vi_alpha_mean
      
      sqrt_pg_weights <- Diagonal(x = sqrt(vi_pg_mean))
      
      for (j in 1:number_of_RE) {
        index_j <- cyclical_pos[[j]]
        M_j <- as(M[, index_j, drop = F], 'dgCMatrix')
        prec_j <- crossprod(sqrt_pg_weights %*% M_j) + Tinv[[j]]
        
        chol_var_j <- solve(t(chol(prec_j)))
        running_log_det_alpha_var[j] <- 2 * sum(log(diag(chol_var_j)))
        
        vi_alpha_decomp[index_j, index_j] <- as(chol_var_j, "dgTMatrix")
      }
      vi_alpha_L_nonpermute <- vi_alpha_decomp
      vi_alpha_LP <- Diagonal(n = nrow(vi_alpha_L_nonpermute))
      vi_alpha_decomp <- vi_alpha_L_nonpermute  %*% t(vi_alpha_LP)
      vi_alpha_decomp <- drop0(vi_alpha_decomp)
      vi_alpha_decomp <- as(vi_alpha_decomp, 'dgCMatrix')
      
      log_det_alpha_var <- sum(running_log_det_alpha_var)
      
      var_ALPHA <- t(vi_alpha_decomp) %*% vi_alpha_decomp
      vi_joint_all <- bdiag(beta_var, var_ALPHA)
      
      vi_joint_all[seq_len(nrow(beta_var)), seq_len(nrow(beta_var))] <- 
        P %*% var_ALPHA %*% t(P) + vi_joint_all[seq_len(nrow(beta_var)), seq_len(nrow(beta_var))]
      vi_joint_all[seq_len(nrow(beta_var)),-seq_len(nrow(beta_var)), drop = F] <- - P %*% var_ALPHA
      vi_joint_all[-seq_len(nrow(beta_var)),seq_len(nrow(beta_var)),drop=F] <- - t(P %*% var_ALPHA)
      
      vi_joint_decomp <- chol(vi_joint_all)

      vi_beta_decomp <- vi_joint_decomp[,1:p.X,drop=F]
      # vi_beta_decomp <- chol(beta_var)
      # vi_beta_L_nonpermute <- vi_beta_decomp
      # vi_beta_LP <- Diagonal(n = nrow(vi_beta_mean))
      # vi_joint_LP <- Diagonal(n = nrow(vi_joint_decomp))
      # vi_joint_L_nonpermute <- vi_joint_decomp
      
      log_det_joint_var <- NA
      log_det_beta_var <- as.numeric(determinant(beta_var)$modulus)
      
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
        i * vi_sigmasq_a/vi_sigmasq_b + j
      })
      
    }else{
        #Update Inverse-Wishart
        variance_by_alpha_jg <- calculate_expected_outer_alpha(L = vi_alpha_decomp, alpha_mu = as.vector(vi_alpha_mean), re_position_list = outer_alpha_RE_positions)
        vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
        
        for (inner_it in 1:INNER_IT){
          
          vi_sigma_alpha <- mapply(vi_sigma_outer_alpha, vi_a_a_jp, 
           vi_a_b_jp, vi_a_nu_jp, SIMPLIFY = FALSE, 
           FUN = function(i, tilde.a, tilde.b, nu) {
             i * vi_sigmasq_a/vi_sigmasq_b + Diagonal(x = tilde.a/tilde.b) * 2 * nu
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
        # d_j <<- d_j
        # vi_alpha_decomp <<- vi_alpha_decomp
        # Tinv <<- Tinv
        # vi_alpha_mean <<- vi_alpha_mean
    }

    if (do_timing) {
      toc(quiet = verbose_time, log = T)
      tic("Update Aux")
    }

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
        log_det_joint_var = log_det_joint_var, 
        vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma, 
        choose_term = choose_term,
        vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
        vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
        do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
      )
    }

    if (parameter_expansion == "none" | !any_Mprime) {
      
      accept.PX <- TRUE
      
    } else {
      
      if (do_timing) {
        tic("Update PX")
      }
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

    quiet_rho <- control$quiet_rho
    
    if (parameter_expansion %in% c("translation", "diagonal") & skip_translate == FALSE & any_Mprime) {
      
      attempted_expansion <- attempted_expansion + 1
      
      if (debug_px){

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
           vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
           vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
           log_det_joint_var = log_det_joint_var, 
           vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma,
           do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
           vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
        )
      }

      if (!quiet_rho){cat('r')}
      
      if (do_timing){
        tic('px_r')
      }
      
      raw_R <- R_ridge <- vecR_ridge_new(L = vi_alpha_decomp[,nonspline_positions], pg_mean = diag(diag_vi_pg_mean),
        mapping_J = mapping_J, d = d_j[!spline_REs],
        store_id = store_id, store_re_id = store_re_id,
        store_design = store_design, 
        diag_only = (factorization_method == 'strong'))

      if (factorization_method == 'weak'){
        stop('no Translation PX for weak yet...')
      }
      
      if (!quiet_rho){cat('r')}

      R_design <- vecR_design(alpha_mu = as.vector(vi_alpha_mean), Z = mapping_new_Z, 
        M = Mmap, mapping_J = mapping_J, d = d_j[!spline_REs],
        start_z = start_base_Z)

      if (sum(spline_REs)){
        R_spline_design <- sapply(cyclical_pos[spline_REs], FUN=function(i){
          as.vector(Z[,i] %*% vi_alpha_mean[i,])
        })
        
        R_spline_ridge <- sapply(cyclical_pos[spline_REs], FUN=function(s){vi_alpha_decomp[,s, drop = F]})
        R_spline_ridge <- Diagonal(x =mapply(R_spline_ridge, cyclical_pos[spline_REs], FUN=function(V, pos){
          sum(vi_pg_mean * cpp_zVz(Z = drop0(Z[,pos,drop=F]), V = as(V, 'dgCMatrix')))
        }))
        # Manually convert "ddiMatrix" to "dgCMatrix" so doesn't fail on
        # old versions of "Matrix" package.
        if (inherits(R_spline_ridge, 'ddiMatrix')){
          R_spline_ridge <- diag(R_spline_ridge)
          R_spline_ridge <- sparseMatrix(
            i = seq_len(length(R_spline_ridge)),
            j = seq_len(length(R_spline_ridge)),
            x = R_spline_ridge)
        }else{
          R_spline_ridge <- as(R_spline_ridge, 'dgCMatrix')
        }
      }else{
        R_spline_ridge <- drop0(matrix(0, nrow = 0, ncol = 0))
        R_spline_design <- matrix(nrow = nrow(X), ncol = 0)
      }
      
      
      if (do_timing){
        toc(quiet = verbose_time, log = TRUE)
        tic('px_fit')
      }
      #If a DIAGONAL expansion, then only update the diagonal elements
      if (parameter_expansion == "diagonal"){
        stop('parameter_expansion "diagonal" turned off.')
        # XR <- cbind(X, R_spline_design, R_design[, diag_rho])
        # R_ridge <- bdiag(zeromat_beta, R_spline_ridge, R_ridge[diag_rho, diag_rho])
        # 
        # if (do_huangwand){
        #   vec_OSL_prior <- do.call('c', mapply(vi_a_APRIOR_jp[!spline_REs], 
        #                                        vi_a_a_jp[!spline_REs], 
        #                                        vi_a_b_jp[!spline_REs],
        #                                        SIMPLIFY = FALSE,
        #     FUN=function(i,a,b){1-2/i^2 * a/b}))
        #   vec_OSL_prior <- c(rep(0, p.X), OSL_spline_prior, vec_OSL_prior)
        # }else{
        #   vec_OSL_prior <- vec_OSL_prior[c(seq_len(p.X + sum(spline_REs)), p.X + sum(spline_REs) + diag_rho),,drop=F]
        # }
        # if (length(vec_OSL_prior) != ncol(XR)){stop('MISALIGNED DIMENSIONS')}
        # 
        # update_expansion_XR <- vecR_fast_ridge(X = drop0(XR), 
        #  omega = diag_vi_pg_mean, prior_precision = R_ridge, y = as.vector(s), 
        #  adjust_y = as.vector(vec_OSL_prior))
        # 
        # update_expansion_bX <- Matrix(update_expansion_XR[1:p.X])
        # update_expansion_splines <- Matrix(update_expansion_XR[-(1:p.X)][seq_len(size_splines)])
        # 
        # update_expansion_R <- mapply(split(update_expansion_XR[-seq_len(p.X + size_splines)], 
        #   rep(1:(number_of_RE - sum(spline_REs)), d_j[!spline_REs])), d_j[!spline_REs], SIMPLIFY = FALSE, 
        #   FUN=function(i,d){
        #     dg <- diag(x = d)
        #     diag(dg) <- i
        #     return(dg)
        #   })
        #  update_diag_R <- split(update_expansion_XR[-seq_len(p.X + size_splines)], 
        #                         rep(1:(number_of_RE - sum(spline_REs)), d_j[!spline_REs]))
        #  rownames(update_expansion_bX) <- colnames(X)
      }else{
        
        XR <- drop0(cbind(drop0(X), drop0(R_spline_design), drop0(R_design)))
        R_ridge <- bdiag(zeromat_beta, R_spline_ridge, R_ridge)
        
        moments_sigma_alpha <- mapply(vi_sigma_alpha, vi_sigma_alpha_nu, d_j, 
            SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
              inv_phi <- solve(phi)
              sigma.inv <- nu * inv_phi
              ln.det <- log(det(phi)) - sum(digamma((nu - 1:d + 1) / 2)) - d * log(2)
              return(list(sigma.inv = sigma.inv, ln.det = ln.det))
            })
        
        if (family == 'linear'){# Rescale for linear
          XR <- XR * sqrt(vi_sigmasq_a/vi_sigmasq_b)
          adj_s <- s * sqrt(vi_sigmasq_a/vi_sigmasq_b)
          R_ridge <- R_ridge * vi_sigmasq_a/vi_sigmasq_b
          offset <- 0
        }else if (family == 'negbin'){
          adj_s <- s
          offset <- vi_r_mu
          stop('translation not set up for negative binomial.')
        }else if (family == 'binomial'){
          adj_s <- s
          offset <- 0
        }else{stop("family not set up for translation expansion.")}
        
        update_expansion_XR <- update_rho(
          XR = XR, y = adj_s, omega = diag_vi_pg_mean, 
          prior_precision = R_ridge, vi_beta_mean = vi_beta_mean,
          moments_sigma_alpha = moments_sigma_alpha,
          prior_sigma_alpha_nu = prior_sigma_alpha_nu, prior_sigma_alpha_phi = prior_sigma_alpha_phi,
          vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp, vi_a_nu_jp = vi_a_nu_jp,
          vi_a_APRIOR_jp = vi_a_APRIOR_jp, 
          stationary_rho = stationary_rho,
          spline_REs = spline_REs, d_j = d_j,
          do_huangwand = do_huangwand, offset = offset,
          p.X = p.X, method = px_method, px_it = px_it,
          init_rho = opt_prior_rho
        )
        
        if (px_method %in% c('numerical_hw', 'profiled', 'dynamic')){
          px_improve <- update_expansion_XR$improvement
          opt_prior_rho <- update_expansion_XR$opt_par
          update_expansion_hw <- update_expansion_XR$hw
          update_expansion_XR <- update_expansion_XR$rho
        }else if (px_method %in% c('numerical', 'OSL')){
          px_improve <- update_expansion_XR$improvement
          opt_prior_rho <- update_expansion_XR <- update_expansion_XR$rho
        }
        opt_prior_rho <- NULL
        
        update_expansion_bX <- Matrix(update_expansion_XR[1:p.X])
        update_expansion_splines <- as.list(update_expansion_XR[-(1:p.X)][seq_len(sum(spline_REs))])
        update_expansion_R <- mapply(split(update_expansion_XR[-1:-(p.X + sum(spline_REs))], 
          rep(1:(number_of_RE - sum(spline_REs)), d_j[!spline_REs]^2)), d_j[!spline_REs], 
          SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d)})
        
      }
      
      if (do_timing){
        toc(quiet = verbose_time, log = TRUE)
        tic('px_propose')
      }
      
      est_rho_all <- update_expansion_XR[-(1:p.X)]
      if (sum(spline_REs)){
        est_rho_spline <- est_rho_all[seq_len(sum(spline_REs))]
        est_rho <- est_rho_all[-seq_len(sum(spline_REs))]
      }else{
        est_rho <- est_rho_all
        est_rho_spline <- 1
      }
      
      if (px_method %in% c('numerical_hw', 'profiled', 'dynamic')){
        
        check_rho_hw <- unlist(vi_a_b_jp[c(which(spline_REs), which(!spline_REs))])
        check_rho_hw <- check_rho_hw - unlist(update_expansion_hw)
        names(check_rho_hw) <- NULL
        
      }else{
        check_rho_hw <- 0
      }
      if (!quiet_rho){
        print(round(c(est_rho_spline, est_rho, check_rho_hw), 5))
      }
      if (parameter_expansion == 'diagonal'){
        if (!is.na(px_improve) & (max(abs(est_rho - 1)) < 1e-6) & (max(abs(est_rho_spline - 1)) < 1e-6) ){
          if (!quiet_rho){print('No further improvements')}
          skip_translate <- TRUE
        }
      }else{
        if (!is.na(px_improve) & (max(abs(est_rho - stationary_rho)) < 1e-6) & (max(abs(est_rho_spline - 1)) < 1e-6) ){
          if (!quiet_rho){print('No further improvements')}
          skip_translate <- TRUE
        }
        if (!is.na(px_improve)){
          if (abs(px_improve) < 1e-7){
            if (!quiet_rho){print('No further improvements (ELBO)')}
            skip_translate <- TRUE
          }
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
        update_expansion_R[spline_REs] <- lapply(update_expansion_splines, FUN=function(i){matrix(i)})
        
        rm(old_update_expansion_R)
      }
      
      prop_vi_sigma_alpha <- mapply(vi_sigma_alpha, update_expansion_R, SIMPLIFY = FALSE,
        FUN=function(Phi, R){R %*% Phi %*% t(R)})
      
      # cat('r')
      # Are any of the estimated "R_j" have a negative determinant?
      sign_detRj <- sign(sapply(update_expansion_R, det))
      any_neg_det <- any(sign_detRj < 0)
      
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
        
        if (!any_neg_det){
          prop_vi_alpha_decomp <- vi_alpha_decomp %*% t(update_expansion_Rblock)
        }else{
          warning(paste0('Manually corrected R_j with negative determinant at iteration ', it))
          if (all(d_j == 1)){
            if (!isDiagonal(update_expansion_Rblock)){
              stop('Correction failed as R_j is not diagonal. Try requiring optimization of PX objective.')
            }
            diag(update_expansion_Rblock) <- abs(diag(update_expansion_Rblock))
            prop_vi_alpha_decomp <- vi_alpha_decomp %*% t(update_expansion_Rblock)
          }else{
            prop_vi_alpha_decomp <- update_expansion_Rblock %*% t(vi_alpha_decomp) %*% 
              vi_alpha_decomp %*% t(update_expansion_Rblock)
            prop_vi_alpha_decomp <- Matrix::Cholesky(prop_vi_alpha_decomp)
            prop_vi_alpha_decomp <- with(expand(prop_vi_alpha_decomp), t(L) %*% P)
          }
        }
        
        prop_log_det_alpha_var <- log_det_alpha_var + 2 * sum(update_expansion_R_logdet * g_j)
        prop_log_det_beta_var <- log_det_beta_var
        prop_vi_beta_decomp <- vi_beta_decomp
        
        prop_variance_by_alpha_jg <- calculate_expected_outer_alpha(
            L = prop_vi_alpha_decomp, 
            alpha_mu = as.vector(prop_vi_alpha_mean), 
            re_position_list = outer_alpha_RE_positions)
        prop_vi_sigma_outer_alpha <- prop_variance_by_alpha_jg$outer_alpha
      }else{
        stop('...')
        # Be sure to set up linear case here too..
      }
      
      if (do_huangwand){
        if (parameter_expansion == "diagonal"){
          prop_vi_a_b_jp <- mapply(vi_a_b_jp, update_diag_R, SIMPLIFY = FALSE,
                              FUN=function(i,j){i / j^2})
          if (px_method != 'OSL'){stop('Double check diagonal expansion')}
        }else{

          if (px_method %in% c('OSL')){
            
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
            
          }else if (px_method == 'numerical'){
            prop_vi_a_b_jp <- vi_a_b_jp
          }else{
            prop_vi_a_b_jp <- update_expansion_hw[names(vi_a_b_jp)]
          }
          

        }
      }else{
        prop_vi_a_b_jp <- NULL
      }
      
      # #L^T L = Variance
      # #R Var R^T --->
      # # L %*% R^T

      if (debug_px){
        prop.ELBO <- calculate_ELBO(family = family,
            ELBO_type = ELBO_type,
            factorization_method = factorization_method,
            d_j = d_j, g_j = g_j, prior_sigma_alpha_phi = prior_sigma_alpha_phi,
            prior_sigma_alpha_nu = prior_sigma_alpha_nu,
            iw_prior_constant = iw_prior_constant,
            X = X, Z = Z, s = s, y = y,
            vi_pg_b = vi_pg_b, vi_pg_mean = vi_pg_mean, vi_pg_c = vi_pg_c,
            vi_sigma_alpha_nu = vi_sigma_alpha_nu,
            
            vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
            vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
            
            vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma, vi_r_mu = vi_r_mu,
            
            vi_sigma_alpha = prop_vi_sigma_alpha, 
            vi_a_b_jp = prop_vi_a_b_jp,
            vi_sigma_outer_alpha = prop_vi_sigma_outer_alpha,
            vi_beta_mean = prop_vi_beta_mean, vi_alpha_mean = prop_vi_alpha_mean,
            
            log_det_beta_var = prop_log_det_beta_var, 
            log_det_alpha_var = prop_log_det_alpha_var,
            log_det_joint_var = prop_log_det_joint_var,
            
            vi_beta_decomp = prop_vi_beta_decomp, 
            vi_alpha_decomp = prop_vi_alpha_decomp,
            vi_joint_decomp = prop_vi_joint_decomp,
            
            do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, 
            vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
            choose_term
        )
      }
      if (!quiet_rho){cat('d')}
      
      # If debugging, check whether the change in ELBO
      # from the profiled objective agrees with the 
      # change from the actual ELBO.
      if (debug_px){
        ELBO_diff <- prop.ELBO$ELBO - prior.ELBO$ELBO
        if (!is.na(px_improve)){
          if (abs(ELBO_diff - px_improve) > sqrt(.Machine$double.eps)){
            stop('PX does not agree with debug.')
            # browser()
            # stop()
          }
        }else{
          if (!isTRUE(all.equal(ELBO_diff, 0))){stop('PX altered parameters when NA.')}
        }
        debug_PX_ELBO[it] <- ELBO_diff
      }
      
      if (is.na(px_improve)){
        accept.PX <- FALSE
      }else if (px_improve > 0){
        accept.PX <- TRUE
      }else{
        accept.PX <- FALSE
      }

      if (accept.PX){
        
        # Accept the PX-VB adjustment
        
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
      }
      
      if (!quiet_rho){
        print(accept.PX)
        if (debug_px){
          out_px <- c(prop.ELBO$ELBO, prior.ELBO$ELBO)
          names(out_px) <- c('PX', 'prior')
          print(out_px)
        }
      }
      
      if (isFALSE(accept.PX) & (px_method %in% c('numerical', 'profiled', 'numerical_hw'))){stop("PX SHOULD NOT FAIL")}
      
      accepted_times <- accept.PX + accepted_times

      if (do_timing){
        toc(quiet = verbose_time, log = TRUE)
      }
      rm(prop_vi_beta_mean, prop_vi_alpha_mean, prop_vi_sigma_alpha, prop_vi_alpha_decomp,
         prop_log_det_alpha_var, prop_variance_by_alpha_jg, prop_vi_sigma_outer_alpha)


      rownames(vi_alpha_mean) <- fmt_names_Z
    }
    
    # Adjust the terms in the ELBO calculation that are different.

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
      vi_sigmasq_a = vi_sigmasq_a, vi_sigmasq_b = vi_sigmasq_b, 
      vi_sigmasq_prior_a = vi_sigmasq_prior_a, vi_sigmasq_prior_b = vi_sigmasq_prior_b,
      log_det_joint_var = log_det_joint_var, vi_r_mu = vi_r_mu, vi_r_mean = vi_r_mean, vi_r_sigma = vi_r_sigma,
      do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
      vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp
    )

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

          if ('vi_alpha_mean' %in% names(prop_squarem)){
            prop_variance_by_alpha_jg <- calculate_expected_outer_alpha(
              L = (prop_squarem$vi_alpha_decomp), 
              alpha_mu = as.vector(prop_squarem$vi_alpha_mean), 
              re_position_list = outer_alpha_RE_positions)
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
            
            if (factorization_method %in% c("weak", "collapsed")) {
              joint_quad <- cpp_zVz(Z = joint.XZ, V = as(prop_ELBOargs$vi_joint_decomp, "dgCMatrix")) 
              if (family == 'negbin'){
                joint_quad <- joint_quad + prop_ELBOargs$vi_r_sigma
              }
              prop_ELBOargs$vi_pg_c <- sqrt(as.vector(X %*% prop_ELBOargs$vi_beta_mean + Z %*% prop_ELBOargs$vi_alpha_mean - prop_ELBOargs$vi_r_mu)^2 + joint_quad)
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
      store_ELBO$step <- NA
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

    change_alpha_mean <- max(abs(vi_alpha_mean - lagged_alpha_mean))
    change_beta_mean <- max(abs(vi_beta_mean - lagged_beta_mean))
    
    unlist_vi <- c(unlist(lapply(vi_sigma_alpha, as.vector)), unlist(vi_a_b_jp))
    if (debug_ELBO){
      svi <- data.frame(t(as.vector(unlist_vi)))
      svi$it <- it
      store_vi <- rbind(store_vi,  svi)
    }
    
    change_sigma_mean <- mapply(vi_sigma_alpha, lagged_sigma_alpha, FUN = function(i, j) {
      max(abs(i - j))
    })

    if (factorization_method %in% c("weak", "collapsed")) {
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
      store_alpha[it, ] <- as.vector(vi_alpha_mean)
      if (do_huangwand){
        store_hw[it,] <- unlist(vi_a_b_jp)
        colnames(store_hw) <- names(unlist(vi_a_b_jp))
      }
      store_sigma[it,] <- unlist(lapply(vi_sigma_alpha, as.vector))
      colnames(store_sigma) <- names(unlist(lapply(vi_sigma_alpha, as.vector)))
      
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

    lagged_alpha_mean <- vi_alpha_mean
    lagged_beta_mean <- vi_beta_mean
    lagged_alpha_decomp <- vi_alpha_decomp
    lagged_beta_decomp <- vi_beta_decomp
    lagged_sigma_alpha <- vi_sigma_alpha
    lagged_vi_r_mu <- vi_r_mu
    lagged_ELBO <- final.ELBO$ELBO
  }
  if (it == iterations) {
    message(paste0("Ended without Convergence after", it, " iterations : ELBO change of ", round(change_elbo[1], abs(floor(log(tolerance_elbo) / log(10))))))
  }

  if (parameter_expansion %in% c("translation", "diagonal")) {
    final.ELBO$accepted_PX <- accepted_times / attempted_expansion
  }

  rownames(vi_beta_mean) <- colnames(X)
  
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
    if (do_huangwand){
      store_hw <- store_hw[1:it,,drop=F]
    }else{store_hw <- NULL}
    store_sigma <- store_sigma[1:it,,drop=F]
    output$parameter_trajectory <- list(beta = store_beta,
                                        alpha = store_alpha,
                                        sigma = store_sigma,
                                        hw = store_hw)
  }
  if (factorization_method %in% c("weak", "collapsed")) {
    output$joint <- vi_joint_decomp
  }
  if (control$return_data) {
    output$data <- list(X = X, Z = Z, y = y, trials = trials)
  }
  
  output$spline <- list(attr = Z.spline.attr, size = Z.spline.size)

  output$formula <- list(formula = formula, 
     re = re_fmla, fe = fe_fmla,
     interpret_gam = parse_formula)
  
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
    acceleration = list(accept.PX = accept.PX, 
      squarem_success = squarem_success, debug_PX_ELBO = debug_PX_ELBO),
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
#' @param factorization_method The factorization method to use. Default of
#'   \code{strong}. Described in detail in Goplerud (2022a). \code{strong},
#'   \code{partial}, and \code{weak} correspond to Schemes I, II, and III
#'   respectively. "weak" should
#'   have best performance but is slowest.
#' @param prior_variance Options are \code{hw}, \code{jeffreys},
#'   \code{mcmcglmm}, \code{mvD}, \code{mean_exists}, limit, and uniform. The
#'   default (\code{hw}) is the Huang-Wand (2013) prior whose hyper-parameters
#'   are nu = 2 and A = 5.
#'   
#'   Otherwise, the prior is an Inverse Wishart with the
#'   following parameters where d is the dimensionality of the random effect.
#'   \itemize{
#'   \item hw: Huang and Wand (2013)
#'   \item mean_exists: IW(d + 1, I)
#'   \item jeffreys: IW(0, 0)
#'   \item mcmcglmm: IW(0, I)
#'   \item mvD: IW(-d, I)
#'   \item limit: IW(d - 1, 0)
#'   \item uniform: IW(-[d+1], 0)
#'   }
#'   The model may fail to be estimable if an improper prior is used. In that
#'   case, use either \code{hw} or \code{mean_exists}.
#' @param tolerance_elbo Change in ELBO to stop algorithm.
#' @param tolerance_parameters Change in value of any parameter to stop
#'   algorithm.
#' @param parameter_expansion Default of \code{translation}  (see Goplerud 2022b).
#'   Accepts either \code{translation}, \code{mean}, or \code{none}. \code{mean}
#'   should be employed if \code{translation} is not enabled or is too
#'   computationally expensive.
#' @param px_method For \code{translation} expansion, how to update? "dynamic"
#'   tries OSL and then backup numerical improvement via L-BFGS-B.
#' @param px_numerical_it How many steps of L-BFGS-B are used to improve?
#' @param hw_INNER For HW prior, how many "loops" between optimizing the
#'   Inverse-Wishart and Inverse-Gamma parameters are done at each iteration?
#' @param prevent_degeneracy Ignored for the moment.
#' @param force_whole Require whole numbers. Set to FALSE to allow "quasi-binomial".
#' @param vi_r_method Type of estimate for "r"; at moment, "fixed" (provide r),
#'   "VEM" (treat r as point estimate; default);
#'   "Laplace" (estimate using Laplace approximation described in the Addendum on GitHub); or "delta"
#'   (experimential).
#' @param vi_r_val For fixed "r", which value?
#'
#' @param init Initialization method can be one of four options: "EM_FE" sets
#'   the random effects to zero, estimates the fixed effects and initializes the
#'   model. "EM" initializes the model with a ridge regression with a guess as
#'   to the random effect variance. "zero" initializes the variational means at
#'   zero. "random" initializes randomly.
#'
#' @param debug_param Debug parameter convergence.
#' @param debug_ELBO Debug ELBO trajectory.
#' @param quiet_rho Debug parameter expansion by printing updates
#' @param debug_px Debug parameter expansion by verifying ELBO
#'
#' @param linpred_method Method for updating means of beta and alpha. "joint" is best.
#' @param print_prog Print after print_prog iterations to show progress.
#' @param quiet Don't print noisy intermediate output.
#' @param return_data Return the design (X,Z) for debugging afterwards.
#'
#' @param verbose_time Print time for each step (debugging only)
#' @param do_timing Estimate timing with tictoc
#' @param do_SQUAREM Accelerate method using SQUAREM
#' @param verify_columns Verify that all columns are drawn from the data.frame itself.
#' 
#' @references 
#' Goplerud, Max. 2022a. "Fast and Accurate Estimation of Non-Nested Binomial
#' Hierarchical Models Using Variational Inference." Bayesian Analysis. 17(2):
#' 623-650.
#'
#' Goplerud, Max. 2022b. "Re-Evaluating Machine Learning for MRP Given the
#' Comparable Performance of (Deep) Hierarchical Models." Working Paper.
#' @export
vglmer_control <- function(iterations = 1000,
   prior_variance = "hw",
   factorization_method = c("strong", "partial", "weak"),
   tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
   prevent_degeneracy = FALSE, force_whole = TRUE, verbose_time = TRUE,
   parameter_expansion = "translation", do_timing = FALSE,
   debug_param = FALSE, return_data = FALSE, linpred_method = "joint",
   vi_r_method = "VEM", vi_r_val = NA, do_SQUAREM = TRUE, verify_columns = FALSE,
   debug_ELBO = FALSE, print_prog = NULL, quiet = T, quiet_rho = TRUE,
   debug_px = FALSE, px_method = 'dynamic', px_numerical_it = 10,
   hw_INNER = 10,
   init = "EM_FE") {
  
  factorization_method <- match.arg(factorization_method)
  prior_variance <- match.arg(prior_variance, 
    choices = c("hw", "mean_exists", "jeffreys", "mcmcglmm", "mvD", "limit", "uniform", "gamma"))
  linpred_method <- match.arg(linpred_method, choices = c("joint", "cyclical", "solve_normal"))    
  parameter_expansion <- match.arg(parameter_expansion, choices = c("translation", "mean", "none"))
  vi_r_method <- match.arg(vi_r_method, choices = c("VEM", "fixed", "Laplace", "delta"))
  init <- match.arg(init, choices = c("EM_FE", "EM", "random", "zero"))
  if (!is.null(print_prog)){
    if (print_prog < 0){stop('print_prog must be non-negative integer or NULL.')}
  }
  
  if (iterations < 0){stop('iterations must be positive integer')}
  if (tolerance_elbo < 0 | tolerance_parameters < 0){
    stop('tolerance for ELBO and parameters must be non-negative.')
  }
  
  if (factorization_method != "strong" & parameter_expansion != "mean"){
    message('Setting parameter_expansion to mean for non-strong factorization')
    parameter_expansion <- 'mean'
  }
  if (prior_variance != 'hw' & px_method != 'OSL' & parameter_expansion %in% c('diagonal', 'translation')){
    px_method <- 'OSL'
    message('Setting px_method to "OSL" if translation & non-HW prior.')
  }
  if (vi_r_method == "fixed" & is.na(vi_r_val)) {
    stop('vi_r_val must not be NA if vi_r_method = "fixed"')
  }

  output <- mget(ls())
  
  class(output) <- c("vglmer_control")
  return(output)
}
