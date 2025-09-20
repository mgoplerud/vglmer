#' Variational Inference for Hierarchical Generalized Linear Models
#'
#' This function estimates hierarchical models using mean-field variational
#' inference. \code{vglmer} accepts standard syntax used for \code{lme4}, e.g.,
#' \code{y ~ x + (x | g)}. Options are described below. Goplerud (2022a; 2022b)
#' provides details on the variational algorithms.
#'
#' @param formula \code{lme4} style-formula for random effects. Typically,
#'   \code{(1 + z | g)} indicates a random effect for each level of variable
#'   \code{"g"} with a differing slope for the effect of variable \code{"z"} and
#'   an intercept (\code{1}); see "Details" for further discussion and how to
#'   incorporate splines.
#' @param data \code{data.frame} containing the outcome and predictors.
#' @param family Options are "binomial", "linear", or "negbin" (experimental).
#'   If "binomial", outcome must be either binary (\eqn{\{0,1\}}) or
#'   \code{cbind(success, failure)} as per standard \code{glm(er)} syntax.
#'   Non-integer values are permitted for binomial if \code{force_whole} is set
#'   to \code{FALSE} in \code{vglmer_control}.
#' @param control Adjust internal options for estimation. Must use an object
#'   created by \link{vglmer_control}.
#'
#' @examples
#' 
#' set.seed(234)
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
#' \donttest{
#' # Run with weaker (i.e. better) approximation
#' vglmer(y ~ x + (x | g),
#'   data = sim_data,
#'   control = vglmer_control(factorization_method = "weak"),
#'   family = "binomial")
#' }
#' 
#' \donttest{
#' # Use a spline on x with a linear outcome
#' vglmer(y ~ v_s(x),
#'   data = sim_data,
#'   family = "linear")
#' }
#' 
#' @details
#' 
#' \bold{Estimation Syntax:} The \code{formula} argument takes syntax designed
#' to be a similar as possible to \code{lme4}. That is, one can specify models
#' using \code{y ~ x + (1 | g)} where \code{(1 | g)} indicates a random intercept. While
#' not tested extensively, terms of \code{(1 | g / f)} should work as expected. Terms
#' of \code{(1 + x || g)} may work, although will raise a warning about duplicated
#' names of random effects. \code{(1 + x || g)} terms may not work with spline
#' estimation. To get around this, one can might copy the column \code{g} to
#' \code{g_copy} and then write \code{(1 | g) + (0 + x | g_copy)}.
#' 
#' \bold{Splines:} Splines can be added using the term \code{v_s(x)} for a
#' spline on the variable \code{x}. These are transformed into hierarchical
#' terms in a standard fashion (e.g. Ruppert et al. 2003) and then estimated
#' using the variational algorithms. At the present, only truncated linear
#' functions (\code{type = "tpf"}; the default) and O'Sullivan splines (Wand and
#' Ormerod 2008) are included. The options are described in more detail at
#' \link{v_s}.
#'
#' It is possible to have the spline vary across some categorical predictor by
#' specifying the \code{"by"} argument such as \code{v_s(x, by = g)}. In effect,
#' this adds additional hierarchical terms for the group-level deviations from
#' the "global" spline. \emph{Note:} In contrast to the typical presentation of
#' these splines interacted with categorical variables (e.g., Ruppert et al.
#' 2003), the default use of \code{"by"} includes the lower order interactions
#' that are regularized, i.e. \code{(1 + x | g)}, versus their unregularized
#' version (e.g., \code{x * g}); this can be changed using the \code{by_re}
#' argument described in \link{v_s}. Further, all group-level deviations from
#' the global spline share the same smoothing parameter (same prior
#' distribution).
#' 
#' \bold{Default Settings:} By default, the model is estimated using the
#' "strong" (i.e. fully factorized) variational assumption. Setting
#' \code{vglmer_control(factorization_method = "weak")} will improve the quality
#' of the variance approximation but may take considerably more time to
#' estimate. See Goplerud (2022a) for discussion. 
#' 
#' By default, the prior on each random effect variance (\eqn{\Sigma_j}) uses a Huang-Wand prior (Huang
#' and Wand 2013) with hyper-parameters \eqn{\nu_j = 2} and \eqn{A_{j,k} = 5}.
#' This is designed to be proper but weakly informative. Other options are
#' discussed in \link{vglmer_control} under the \code{prior_variance} argument.
#' 
#' By default, estimation is accelerated using SQUAREM (Varadhan and Roland
#' 2008) and (one-step-late) parameter expansion for variational Bayes. Under
#' the default \code{"strong"} factorization, a "translation" expansion is used;
#' under other factorizations a "mean" expansion is used. These can be adjusted
#' using \link{vglmer_control}. See Goplerud (2022b) for more discussion of
#' these methods.
#' 
#' @return This returns an object of class \code{vglmer}. The available methods
#'   (e.g. \code{coef}) can be found using \code{methods(class="vglmer")}.
#' \describe{
#' \item{beta}{Contains the estimated distribution of the fixed effects
#' (\eqn{\beta}). It is multivariate normal. \code{mean} contains the means;
#' \code{var} contains the variance matrix; \code{decomp_var} contains a matrix
#' \eqn{L} such that \eqn{L^T L} equals the full variance matrix.}
#' \item{alpha}{Contains the estimated distribution of the random effects
#' (\eqn{\alpha}). They are all multivariate normal. \code{mean} contains the
#' means; \code{dia.var} contains the variance of each random effect. \code{var}
#' contains the variance matrix of each random effect (j,g). \code{decomp_var}
#' contains a matrix \eqn{L} such that \eqn{L^T L} equals the full variance of
#' the entire set of random effects.}
#' \item{joint}{If \code{factorization_method="weak"}, this is a list with one
#' element (\code{decomp_var}) that contains a matrix \eqn{L} such that \eqn{L^T
#' L} equals the full variance matrix between the fixed and random effects
#' \eqn{q(\beta,\alpha)}. The marginal variances are included in \code{beta} and
#' \code{alpha}. If the factorization method is not \code{"weak"}, this is
#' \code{NULL}.}
#' \item{sigma}{Contains the estimated distribution of each random
#' effect covariance \eqn{\Sigma_j}; all distributions are Inverse-Wishart.
#' \code{cov} contains a list of the estimated scale matrices. \code{df}
#' contains a list of the degrees of freedom.}
#' \item{hw}{If a Huang-Wand prior is used (see Huang and Wand 2013 or Goplerud
#' 2022b for more details), then the estimated distribution. Otherwise, it is
#' \code{NULL}. All distributions are Inverse-Gamma. \code{a} contains a list of
#' the scale parameters. \code{b} contains a list of the shape parameters.}
#' \item{sigmasq}{If \code{family="linear"}, this contains a list of the
#' estimated parameters for \eqn{\sigma^2}; its distribution is Inverse-Gamma.
#' \code{a} contains the scale parameter; \code{b} contains the shape
#' parameter.}
#' \item{ln_r}{If \code{family="negbin"}, this contains the variational
#' parameters for the log dispersion parameter \eqn{\ln(r)}. \code{mu} contains
#' the mean; \code{sigma} contains the variance.}
#' \item{family}{Family of outcome.}
#' \item{ELBO}{Contains the ELBO at the termination of the algorithm.}
#' \item{ELBO_trajectory}{\code{data.frame} tracking the ELBO per iteration.}
#' \item{control}{Contains the control parameters from \code{vglmer_control}
#' used in estimation.}
#' \item{internal_parameters}{Variety of internal parameters used in
#' post-estimation functions.}
#' \item{formula}{Contains the formula used for estimation; contains the
#' original formula, fixed effects, and random effects parts separately for
#' post-estimation functions. See \code{formula.vglmer} for more details.}
#' }
#' @importFrom lme4 mkReTrms findbars subbars
#' @importFrom stats model.response model.matrix model.frame rnorm rWishart
#'   qlogis optim residuals lm plogis setNames .getXlevels
#' @importFrom graphics plot
#' @importFrom Rcpp sourceCpp
#' @importFrom mgcv interpret.gam
#' @references
#' Goplerud, Max. 2022a. "Fast and Accurate Estimation of Non-Nested Binomial
#' Hierarchical Models Using Variational Inference." \emph{Bayesian Analysis}. 17(2):
#' 623-650.
#' 
#' Goplerud, Max. 2022b. "Re-Evaluating Machine Learning for MRP Given the
#' Comparable Performance of (Deep) Hierarchical Models." Working paper.
#'
#' Huang, Alan, and Matthew P. Wand. 2013. "Simple Marginally Noninformative
#' Prior Distributions for Covariance Matrices." \emph{Bayesian Analysis}.
#' 8(2):439-452.
#' 
#' Ruppert, David, Matt P. Wand, and Raymond J. Carroll. 2003.
#' \emph{Semiparametric Regression}. Cambridge University Press.
#' 
#' Varadhan, Ravi, and Christophe Roland. 2008. "Simple and Globally Convergent
#' Methods for Accelerating the Convergence of any EM Algorithm." \emph{Scandinavian
#' Journal of Statistics}. 35(2): 335-353.
#' 
#' Wand, Matt P. and Ormerod, John T. 2008. "On Semiparametric Regression with
#' O'Sullivan Penalized Splines". \emph{Australian & New Zealand Journal of Statistics}.
#' 50(2): 179-198.
#' 
#' @useDynLib vglmer
#' @export
vglmer <- function(formula, data, family, control = vglmer_control()) {

  # Verify integrity of parameter arguments
  family <- match.arg(family, choices = c("negbin", "binomial", "linear"))
  if (family == "negbin" & !(control$parameter_expansion %in% c('none', 'mean'))){
    message('Setting parameter_expansion to mean for negative binomial estimation')
    control$parameter_expansion <- 'mean'
  }
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

  # Interpret gam using mgcv::interpret.gam
  parse_formula <- NULL
  # parse_formula <- tryCatch(
  #   interpret.gam(subbars(formula), extra.special = c('v_s', 'v_mi')), error = function(e){NULL})
  if (is.null(parse_formula)){
    # If this fails, usually when there is custom argument in environment, use this instead
    parse_formula <- fallback_interpret.gam0(subbars(formula), extra.special = c('v_s', 'v_mi', 'v_hier_mi'))
  }
  
  if (any(!sapply(parse_formula$smooth.spec, inherits, what = c('vglmer_spline', 'vglmer_multiplicative')))){
    stop('gam specials are not permitted; use v_s(...) or v_mi(...) and see documentation.')
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
  
  tt <- terms(data)

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
  debug_param <- control$debug_param
  linpred_method <- control$linpred_method
  vi_r_method <- control$vi_r_method
  freeze_mi_var <- control$freeze_mi_var
  mi_prior_type <- control$mi_prior_type
  if (any(control$mi_parameter_expansion %in% 'none')){
    do_PX_MI <- FALSE
  }else{
    do_PX_MI <- TRUE
    px_mi_type <- control$mi_parameter_expansion
  }
    
  if (is.numeric(vi_r_method)){
    if (length(vi_r_method) > 1){stop('If "vi_r_method" is numeric, it must be a single number.')}
    vi_r_val <- as.numeric(vi_r_method)
    vi_r_method <- "fixed"
  }else{
    vi_r_val <- NA
  }
  debug_ELBO <- control$debug_ELBO
  # Flip given that "tictoc" accepts "quiet=quiet_time"
  quiet_time <- !control$verbose_time

  if (do_timing) {
    if (!requireNamespace("tictoc", quietly = TRUE)) {
      stop("tictoc must be installed to do timing")
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
    stop("factorization_method must be 'weak', 'strong', or 'partial'.")
  }
  if (is.null(print_prog)) {
    print_prog <- max(c(1, floor(iterations / 20)))
  }
  if (!(family %in% c("binomial", "negbin", "linear"))) {
    stop('family must be one of "linear", "binomial", "negbin".')
  }
  
  vi_alpha_L_nonpermute <- vi_alpha_LP <- NULL
  vi_beta_L_nonpermute <- vi_beta_LP <- NULL
  vi_alpha_L_nonpermute <- variance_by_alpha_jg <- NULL
  vi_joint_L_nonpermute <- vi_joint_LP <- NULL
  
  
  if (family == "binomial") {
    if (is.matrix(y)) {
      # if (!(class(y) %in% c('numeric', 'integer'))){
      if (min(y) < 0) {
        stop("Negative numbers not permitted in outcome")
      }
      is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
      if (any(is.wholenumber(y) == FALSE)) {
        if (control$force_whole) {
          stop("If force_whole = TRUE, must provide whole numbers as outcome")
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
    
    if (is.matrix(y)) {
      stop('"linear" family requires a vector outcome.')
    }
    
    if (!(class(y) %in% c("numeric", "integer"))) {
      stop("Must provide vector of numbers with negbin.")
    }

    if (min(y) < 0) {
      stop("Negative numbers not permitted in outcome")
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
    
    if (is.matrix(y)) {
      stop('"linear" family requires a vector outcome.')
    }
    if (!(class(y) %in% c("numeric", "integer"))) {
      stop("Must provide vector of numbers with linear.")
    }
    
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
  fe_fmla <- NULL
  # fe_fmla <- tryCatch(
  #   interpret.gam(nobars(formula), extra.special = c('v_s', 'v_mi')), error = function(e){NULL})
  if (is.null(fe_fmla)){
    # If this fails, usually when there is custom argument in environment, use this instead
    fe_fmla <- fallback_interpret.gam0(nobars(formula), extra.special = c('v_s', 'v_mi', 'v_hier_mi'))
  }

  if (length(fe_fmla$smooth.spec) > 0){
    
    # Add the linear spline terms to the main effect.
    fe_update <- sapply(fe_fmla$smooth.spec, FUN=function(i){
      if (i$mi){return(NULL)}
      if (i$by != "NA" & i$by_re == FALSE){
        fe_i <- paste0(i$term, ' * ', i$by)
      }else{
        fe_i <- i$term
      }
    })
    if (!all(sapply(fe_update, is.null))){
      fe_update <- paste0(fe_update, collapse = ' + ')
      
      fe_fmla <- update.formula(fe_fmla$pf,
                                paste0('. ~ . + 1 + ', fe_update)
      )
    }else{
      fe_fmla <- fe_fmla$pf
    }
    
  }else{
    fe_fmla <- fe_fmla$pf
  }
  
  # Create the FE design
  X <- model.matrix(fe_fmla, data = data)
  fe_terms <- terms(fe_fmla)
  fe_Xlevels <- .getXlevels(fe_terms, data)
  fe_contrasts <- attr(X, 'contrasts')
  
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
      stop('Some grouping factors for random effects are duplicated. Reformulate initial formula.')
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
  
  if (!is.null(re_fmla) & (length(re_fmla) > 0)){
    mk_Z <- mkReTrms(re_fmla, data, reorder.terms = FALSE, reorder.vars = FALSE)
    Z <- t(mk_Z$Zt)
    
    p.X <- ncol(X)
    p.Z <- ncol(Z)
    
    ####
    # Process the REs to get various useful terms.
    ####
    # RE names and names of variables included for each.
    names_of_RE <- mk_Z$cnms
    
    if (anyDuplicated(names(names_of_RE)) > 0){
      warning('Some random effects names are duplicated. Re-naming for stability by adding "-[0-9]" at end.')
      nre <- names(names_of_RE)
      unre <- unique(nre)
      for (u in unre){
        nre_u <- which(nre == u)
        if (length(nre_u) > 1){
          nre[nre_u] <- paste0(nre[nre_u], '-', seq_len(length(nre_u)))
        }
      }
      names(names_of_RE) <- nre
      if (anyDuplicated(names(names_of_RE)) > 0){
        stop('Renaming duplicates failed. Please rename random effects to proceed.')
      }
    }

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
  
  M.names <- cbind(unlist(mapply(names_of_RE, g_j, SIMPLIFY = FALSE, FUN = function(i, j) {
    rep(i, j)
  })))
  
  if (!is.null(M.names)){
    U_names <- unique(cbind(rep(names(names_of_RE), g_j * d_j), M.names))
    B_j <- lapply(split(U_names[,2], U_names[,1]), FUN=function(j){
      B_jj <- which(j %in% colnames(X))
      if (length(B_jj) == 0){
        return(matrix(nrow = length(j), ncol = 0))
      }
      sparseMatrix(i = B_jj, j = 1:length(B_jj), 
        x = 1, dims = c(length(j), length(B_jj)))   
    })
    U_names_Bj <- unique(U_names[,2])
    M <- cbind(match(M.names[, 1], U_names_Bj), rep(1 / g_j, d_j * g_j))
    M <- sparseMatrix(i = 1:nrow(M), j = M[, 1], x = M[, 2], dims = c(ncol(Z), length(U_names_Bj)))
  }else{
    B_j <- list()
    M <- drop0(matrix(0, nrow = 0, ncol = ncol(X)))
  }

  if (!is.null(names_of_RE)){
    any_Mprime <- TRUE
    M_prime.names <- paste0(rep(names(names_of_RE), g_j * d_j), " @ ", M.names)
    M_prime <- cbind(match(M_prime.names, unique(M_prime.names)), rep(1 / g_j, d_j * g_j))
    M_prime <- sparseMatrix(i = seq_len(ncol(Z)),
                            j = M_prime[, 1], 
                            x = M_prime[, 2],
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
  re_type <- setNames(rep('RE', length(names_of_RE)), names(names_of_RE))
  
  if (length(parse_formula$smooth.spec) > 0){
    
    any_Mprime <- TRUE
    
    base_specials <- length(parse_formula$smooth.spec)
    
    # Number of specials (# multiplicative, # splines + one for each "by")
    
    n.specials <- base_specials +
      sum(sapply(parse_formula$smooth.spec, 
        FUN=function(i){!i$mi & (i$by != "NA")}))
    
    Z.special.attr <- as.list(rep(NA, base_specials))

    Z.special.type <- rep(NA, n.specials)
    Z.special <- as.list(rep(NA, n.specials))
    Z.special.size <- rep(NA, n.specials)
    special_counter <- 0

    for (i in 1:base_specials){
      
      special_i <- parse_formula$smooth.spec[[i]]

      if (special_i$mi){
        
        all_special_i <- vglmer_build_mi(
          data = data,
          object = special_i
        )
        
        Z.special.attr[[i]] <- c(all_special_i[[1]], 
            list(type = 'mi', prior_var = special_i$prior_var,
            by = special_i$by))
        Z.special.attr[[i]]$attr$xt <- special_i$xt
        # Find REs with intercept that correspond to factors in MI
        re_mapping_mi <- lapply(match(special_i$term, names(names_of_RE)), FUN=function(i){
          if (!is.na(i)){
            if ('(Intercept)' %in% names_of_RE[[i]]){
              stem_i <- paste0('^', names(names_of_RE)[i], ' @ \\(Intercept\\) @')
              pos_i <- grep(colnames(Z), pattern=stem_i)
              val_i <- colnames(Z)[pos_i]
              return(list(index = i, position = pos_i, value = val_i))
            }else{
              return(NULL)
            }
          }else{
            return(NULL)
          }
        })
        names(re_mapping_mi) <- special_i$term
        Z.special.attr[[i]]$attr$mapping_RE <- re_mapping_mi
      }else{

        all_special_i <- vglmer_build_spline(
          object = special_i, data = data
        )
        
        Z.special.attr[[i]] <- c(all_special_i[[1]]$attr, 
                                 list(type = special_i$type, 
                                      names = special_i$term,
                                      by = special_i$by))
      }
      
      inner_counter <- 1
      for (si in all_special_i){
        
        special_counter <- special_counter + 1
        stopifnot(inner_counter %in% 1:2)
        
        if (special_i$mi){

          if (si$attr$hier & !(mi_prior_type %in% c('separate', 'partial_fixed', 'fixed'))){
            stop('v_mi_hier requires "separate", "fixed", or "partial_fixed" prior type.')
          }
          
          Z.special[[special_counter]] <- si$x
          if (mi_prior_type %in% c('centered', 'shared', 'separate')){
            if (!is.null(special_i$xt)){
              stop('xt must not be provided for "centered", "shared", or "separate"')
            }
          }
          if (mi_prior_type %in% c('centered')){
            Z.special.size[special_counter] <- sapply(si$attr$storage, nrow)[2]
          }else if (mi_prior_type %in% 'shared'){
            Z.special.size[special_counter] <- sum(
              sapply(si$attr$storage, nrow)
            )
          }else if (mi_prior_type %in% c('separate')){
            Z.special.size[special_counter] <- NA
          }else if (mi_prior_type %in% c('fixed', 'partial_fixed')){
            Z.special.size[special_counter] <- NA
            if (special_i$hier){
              # Is the first fixed, but all others omitteid
              check_1 <- all(sapply(special_i$hier_fmla, FUN=function(i){
                all(isTRUE(i[1] %in% names(special_i$xt)))
              }))            
              check_2 <- isTRUE(length(setdiff(names(special_i$xt), 
                sapply(special_i$hier_fmla, FUN=function(i){i[1]})))==0)
            }else{
              check_1 <- all(names(special_i$xt) %in% names(special_i$term))
              check_2 <- all(names(special_i$term) %in% names(special_i$xt))
            }
            if (!(check_1 & check_2)){
              if (special_i$hier){
                stop('If mi_prior_type == "partial_fixed", then variances must be provided for only the first term of q(u) and q(v) using the xt argument.')
              }else{
                stop('If mi_prior_type == "fixed", then variances must be provided for each term using the xt argument.')
              }
            }
          }else{stop('invalid mi prior type')}

          if (si$attr$hier){
            spline_name <- sapply(special_i$hier_fmla, FUN=function(i){i[1]})
            spline_name <- paste0('mi-', paste(spline_name, collapse=','))
          }else{
            spline_name <- paste0('mi-', paste(special_i$term, collapse=','))
          }
          
          names_of_RE[[spline_name]] <- spline_name
          d_j <- setNames(c(d_j, special_i$rank), c(names(d_j), spline_name))
          g_j <- setNames(c(g_j, Z.special.size[special_counter]), c(names(g_j), spline_name))
          Z.special.type[special_counter] <- 'mi'
          re_type <- setNames(c(re_type, 'mi'), c(names(re_type), spline_name))
          
        }else{
          
          colnames(si$x) <- paste0('spline @ ', special_i$term, ' @ ', colnames(si$x))
          
          if (inner_counter > 1){
            spline_name <- paste0('spline-', special_i$term,'-', i, '-int')
          }else{
            spline_name <- paste0('spline-', special_i$term, '-', i, '-base')
          }
          
          Z.special[[special_counter]] <- si$x
          Z.special.size[special_counter] <- ncol(si$x)
          Z.special.type[special_counter] <- 'spline'
          
          names_of_RE[[spline_name]] <- spline_name
          number_of_RE <- number_of_RE + 1
          d_j <- setNames(c(d_j, 1), c(names(d_j), spline_name))
          g_j <- setNames(c(g_j, ncol(si$x)), c(names(g_j), spline_name))
          re_type <- setNames(c(re_type, 'spline'), c(names(re_type), spline_name))
          breaks_for_RE <- c(breaks_for_RE, max(breaks_for_RE) + ncol(si$x))
          fmt_names_Z <- c(fmt_names_Z, colnames(si$x))
          p.Z <- p.Z + ncol(si$x)
        }
        
      }
    }
    
    # Set up options for MI prior variance
    if (control$mi_prior_variance == 'hw'){
      do_huangwand_mi <- TRUE
      mi_prior_variance <- 'huangwand'
    }else{
      do_huangwand_mi <- FALSE
      mi_prior_variance <- 'mean_exists'
    }
    partial_fix <- FALSE
    
    if (mi_prior_type == 'fixed'){
      mi_prior_type <- 'separate'
      mi_prior_variance <- 'fixed'
      do_huangwand_mi <- FALSE
      freeze_mi_var <- TRUE
    }else if (mi_prior_type == 'partial_fixed'){
      mi_prior_type <- 'separate'
      partial_fix <- TRUE
    }
    
    if (number_of_RE > 0){
      cyclical_pos <- lapply(1:number_of_RE, FUN = function(i) {
        seq(breaks_for_RE[i] + 1, breaks_for_RE[i + 1])
      })
    }
    
    Z_MI <- Z.special[which(Z.special.type == 'mi')]
    Z_MI_attr <- lapply(Z.special.attr[which(Z.special.type == 'mi')], `[[`, 'attr')
    if (any(Z.special.type %in% 'spline')){
      Z.special <- drop0(do.call('cbind', Z.special[which(Z.special.type == 'spline')]))
      Z <- drop0(cbind(Z, Z.special))
      any_spline <- TRUE
    }else{
      any_spline <- FALSE
    }
    if (ncol(Z) != p.Z){stop('...')}
    if (max(breaks_for_RE) != p.Z){stop('...')}
    if (any_spline){

      if (ncol(M_prime) == 0){
        M_prime <- rbind(M_prime, 
           drop0(matrix(0, nrow = ncol(Z.special), ncol = 0)))
        M_prime_one <- rbind(M_prime_one, 
           drop0(matrix(0, nrow = ncol(Z.special), ncol = 0)))
      }else{
        M_prime <- rbind(M_prime, 
           drop0(sparseMatrix(i = 1, j = 1, x = 0, 
                            dims = c(ncol(Z.special), ncol(M_prime)))))
        M_prime_one <- rbind(M_prime_one, 
           drop0(sparseMatrix(i = 1, j = 1, x = 0, 
                              dims = c(ncol(Z.special), ncol(M_prime_one)))))
      }

      M_prime <- cbind(M_prime, drop0(sparseMatrix(i = 1, j = 1, x = 0,
                                                   dims = c(nrow(M_prime), special_counter - sum(re_type == 'mi'))
      )))
      M_prime_one <- cbind(M_prime_one, drop0(sparseMatrix(i = 1, j = 1, x = 0,
                                                           dims = c(nrow(M_prime_one), special_counter  - sum(re_type == 'mi'))
      )))
     
      M_mu_to_beta <- rbind(M_mu_to_beta, 
                            drop0(sparseMatrix(i = 1, j = 1, x = 0, 
                                               dims = c(special_counter - sum(re_type == 'mi'), p.X)))
      )
      
    }
  
    extra_Bj <- lapply(setdiff(names(names_of_RE[re_type != 'mi']), names(B_j)), FUN=function(i){Diagonal(n = d_j[i])})
    names(extra_Bj) <- setdiff(names(names_of_RE[re_type != 'mi']), names(B_j))
    B_j <- c(B_j, extra_Bj)[names(names_of_RE[re_type != 'mi'])]
    
  }else{
    n.specials <- 0
    Z.special.attr <- NULL
    Z.special <- NULL
    Z.special.size <- NULL
    Z_MI <- NULL
  }
  
  if (length(B_j) > 0){
    B_j <- bdiag(B_j)
  }else{B_j <- NULL}
  
  debug_px <- control$debug_px
  if (control$parameter_expansion %in% c('translation', 'diagonal')){
    px_method <- control$px_method
    px_it <- control$px_numerical_it
    opt_prior_rho <- NULL
    parsed_RE_groups <- get_RE_groups(formula = re_fmla, data = data)
  }
  
  # Adjust certain things if multiplicative interactions are used
  any_mi <- any(re_type == 'mi')
  any_RE <- any(re_type != 'mi')
  
  # Update MI from first iteration; do not use VEM
  MI_VEM_THRESH <- 0
  MI_STAB <- 0
  # Turn off unusual initalization options but don't fully remove them yet...
  INIT_MCMC <- FALSE
  ADD_INT <- FALSE
  freeze_RE <- FALSE
  do_huangwand_mi <- NA
  if (any_mi){
    
    # Other estimation options to hide but not remove
    MAX_K <- 1
    bi_method <- 'rowwise'
    
    easy_message <- function(x){
      message(paste0(capture.output(x), collapse = "\n"))
    }
    mi_d_j <- d_j[re_type == 'mi']
    mi_g_j <- g_j[re_type == 'mi']
    mi_names_of_RE <- names_of_RE[re_type == 'mi']
    names(Z_MI) <- mi_names_of_RE
    names(Z_MI_attr) <- mi_names_of_RE
    if (mi_prior_type %in% c('separate')){
      mi_g_j <- lapply(Z_MI, FUN=function(i){sapply(i, ncol)})
      names(mi_g_j) <- mi_names_of_RE
    }
    
    d_j <- d_j[re_type != 'mi']
    g_j <- g_j[re_type != 'mi']
    names_of_RE <- names_of_RE[re_type != 'mi']
  }
  if (any_RE){
    # List of Lists
    # Outer list: one for RE
    # Inner List: One for each GROUP with its row positions.
    outer_alpha_RE_positions <- mapply(d_j, g_j, breaks_for_RE[-length(breaks_for_RE)], 
                                       SIMPLIFY = FALSE, FUN = function(a, b, m) {
                                         split(m + seq(1, a * b), rep(1:b, each = a))
                                       })
    
    if (anyDuplicated(unlist(outer_alpha_RE_positions)) != 0 | max(unlist(outer_alpha_RE_positions)) != ncol(Z)) {
      stop("Issue with creating OA positions")
    }
  }else{
    outer_alpha_RE_positions <- NULL
    vi_sigma_outer_alpha <- NULL
  }
  
  ####
  # Prepare Initial Values
  ###

  vi_sigmasq_prior_a <- 0
  vi_sigmasq_prior_b <- 0
  
  vi_sigmasq_a <- vi_sigmasq_b <- 1
  
  if (family == "linear") {

    vi_sigmasq_a <- (nrow(X) + sum(d_j * g_j))/2 + vi_sigmasq_prior_a
    vi_sigmasq_b <- sum(residuals(lm(y ~ 1))^2)/2 + vi_sigmasq_prior_b
    
    s <- y
    vi_pg_b <- 1
    vi_pg_c <- NULL
    vi_r_mu <- 0
    vi_r_sigma <- 0
    vi_r_mean <- 0
    

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
        stop("Install MASS to use negbin")
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
      stop("vi_r_method must be 'VEM' or 'fixed'.")
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
    INNER_IT <- control$hw_inner
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
    stop("Invalid option for prior variance provided.")
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
        if (ncol(X) == 0){
          EM_init <- list(beta = double(), alpha = rep(0, ncol(Z)))
        }else{
          EM_init <- EM_prelim_logit(X = X, Z = Z, s = s, pg_b = vi_pg_b, iter = 15, ridge = 4)
        }
      }
    }

    vi_beta_mean <- matrix(EM_init$beta)

    if (any_RE){
      vi_alpha_mean <- matrix(EM_init$alpha)
      vi_sigma_alpha <- calculate_expected_outer_alpha(
        alpha_mu = vi_alpha_mean,
        L = sparseMatrix(i = 1, j = 1, x = 1e-4, dims = rep(ncol(Z), 2)),
        re_position_list = outer_alpha_RE_positions
      )
      
      # Update Inverse-Wishart
      if (!do_huangwand){
        vi_sigma_alpha <- mapply(vi_sigma_alpha$outer_alpha, prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN = function(i, j) {
          i + j
        })
      }else{
        
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
      
    }else{
      vi_sigma_alpha <- NULL
      ln_det_sigma_alpha <- NULL
      log_det_alpha_var <- 0
      vi_alpha_mean <- matrix(nrow =0,ncol=1)
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

    if (ncol(X) > 0){
      if (family == "binomial") {
        vi_beta_mean[1] <- qlogis(sum(y) / sum(trials))
      } else if (family == "negbin") {
        vi_beta_mean[1] <- log(mean(y))
      } else if (family == 'linear'){
        vi_beta_mean[1] <- mean(y)
      } else {
        stop('Set up init')
      }
    }

    vi_alpha_mean <- rep(0, ncol(Z))

    vi_sigma_alpha <- mapply(d_j, g_j, SIMPLIFY = FALSE, FUN = function(d, g) {
      diag(x = 1, ncol = d, nrow = d)
    })
    # if (do_huangwand){stop('Setup init for zero')}
  } else {
    stop("Invalid initialization method")
  }

  if (ncol(X) == 0){
    zero_mat <- matrix(ncol = 0, nrow = 0)
  }else{
    zero_mat <- sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(X), ncol(X)))
    zero_mat <- drop0(zero_mat)
  }
  
  if (factorization_method %in% c("weak", "collapsed")) {
    vi_joint_decomp <- bdiag(vi_beta_decomp, vi_alpha_decomp)
    joint.XZ <- cbind(X, Z)
    log_det_beta_var <- log_det_alpha_var <- NULL
  } else {
    vi_joint_decomp <- NULL
    log_det_joint_var <- NULL
  }

  
  Z_MI_grouping <- NULL
  if (any_mi){
    
    N_MI <- length(Z_MI)
    names_of_MI <- names(Z_MI)
    
    if (family == 'linear'){
      init_y <- (s - as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean))
      init_w <- rep(1, length(init_y))
    }else if (family == 'binomial'){
      init_lp <- as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean)
      init_w <- vi_pg_b / (2 * init_lp) * tanh(init_lp / 2)
      fill_zero <- which(abs(init_w) < 1e-6)
      if (length(fill_zero) > 0){
        init_w[fill_zero] <- vi_pg_b[fill_zero] / 4
      }
      
      init_y <- s/init_w - init_lp
      rm(init_lp); gc()
    }else{stop('set up weight init for MI')}
    
    Z_MI_grouping <- lapply(Z_MI_attr, `[[`, 'hier_grouping')
    Z_MI_hier <- sapply(Z_MI_attr, `[[`, 'hier')
    Z_MI_hier_size <- sapply(Z_MI, length)
    Z_MI_nested_hier <- lapply(Z_MI_attr, `[[`, 'nested_hier')
    
    init_vi_mi <- mapply(1:N_MI, Z_MI_hier, Z_MI_grouping, Z_MI_nested_hier, SIMPLIFY = FALSE, FUN=function(j, j_hier, j_group, j_nesting){
      if (j_hier){
        if (control$mi_init == 'random'){
          principal_group <- sapply(j_group, FUN=function(i){i[1]})
          out <- lapply(Z_MI[[j]][principal_group], FUN=function(i){
            matrix(rnorm(mi_d_j[j] * ncol(i)), ncol = mi_d_j[j])
          })
          names(out) <- c('u', 'v')
          out_var <- NULL
        }else{
          principal_group <- sapply(j_group, FUN=function(i){i[1]})
          out <- init_MI_from_svd(
            data_mi = Z_MI[[j]][principal_group],
            rank = mi_d_j[j],
            y = init_y, pg_weight = init_w,
            prior_U = 1, prior_V = 1
          )
          out_var <- out[c('var_U', 'var_V')]
          out <- out[c('mean_U', 'mean_V')]
          names(out) <- c('u', 'v')
        }
        out_data <- lapply(1:2, FUN=function(k){
          if (k == 1){
            init_principal_k <- out$u
          }else{
            init_principal_k <- out$v
          }
          principal_k <- j_group[[k]][1]
          other_levels <- j_group[[k]][-1]
          levels_k <- Z_MI_attr[[j]]$levels[j_group[[k]]]
          if (length(other_levels) > 0){
            nest_k <- j_nesting[[k]]
            init_nest_k <- init_principal_k[match(nest_k[,1], levels_k[[principal_k]]),,drop=FALSE]
            dat_nest_k <- data.frame(nest_k, stringsAsFactors = FALSE)
            fmla_k <- paste(paste0('(1 | ', colnames(dat_nest_k)[-1], ')'), collapse = " + ")
            fmla_k <- as.formula(paste('ideal ~ ', fmla_k))
            fit_nest_k <- lapply(1:ncol(init_nest_k), FUN=function(d){
              dat_nest_k$ideal <- as.vector(init_nest_k[,d])
              
              if ('const' %in% names(dat_nest_k)){
                fit_nest_k <- lmer(update.formula(fmla_k, '. ~ . - (1 | const)'), data = dat_nest_k)
                resid_k <- dat_nest_k$ideal - fitted(fit_nest_k) 
                const_k <- fixef(fit_nest_k)
                fit_nest_k <- ranef(fit_nest_k)
                fit_nest_k <- lapply(fit_nest_k, FUN=function(i){
                  data.frame(id = rownames(i), dim = i[,1], stringsAsFactors = F)
                })
                fit_nest_k <- c(
                  list("const" = data.frame(id = "1", dim = const_k, stringsAsFactors = FALSE)),
                  fit_nest_k
                )
              }else{
                fit_nest_k <- lmer(update.formula(fmla_k, '. ~ 0 + .'), data = dat_nest_k)
                resid_k <- dat_nest_k$ideal - fitted(fit_nest_k) 
                fit_nest_k <- ranef(fit_nest_k)
                fit_nest_k <- lapply(fit_nest_k, FUN=function(i){
                  data.frame(id = rownames(i), dim = i[,1], stringsAsFactors = F)
                })
              }
              # fit_nest_k <- vglmer(fmla_k, data = dat_nest_k, 
              #                      control = vglmer_control(iterations = 15),
              #                      family = 'linear')
              # Add in the residuals as the initial estimates
              # for the principal grouping
              fit_nest_k <- c(setNames(list(
                data.frame(
                  id = dat_nest_k[,1],
                  dim = resid_k
                )), names(dat_nest_k)[1]), 
                fit_nest_k)
              
              fit_nest_k <- lapply(fit_nest_k, FUN=function(i){
                if (control$mi_init == 'random'){
                  i[,2] <- rnorm(nrow(i))
                }
                names(i)[2] <- paste0('dim_', d)
                return(i)
              })
              fit_nest_k <- fit_nest_k[setdiff(names(dat_nest_k), 'ideal')]
              return(fit_nest_k)
            })
            fit_nest_k <- lapply(j_group[[k]], FUN=function(i){
              Reduce(merge, lapply(fit_nest_k, `[[`, i))
            })
            fit_nest_k <- mapply(fit_nest_k, levels_k, SIMPLIFY = FALSE, FUN=function(i,levels_i){
              out <- as.matrix(i[match(levels_i, i[,1]),-1])
              rownames(out) <- levels_i
              return(out)
            })
            names(fit_nest_k) <- names(levels_k)
          }else{
            init_principal_k <- as.matrix(init_principal_k)
            rownames(init_principal_k) <- levels_k[[principal_k]]
            colnames(init_principal_k) <- paste0('dim', 1:ncol(init_principal_k))
            fit_nest_k <- setNames(list(init_principal_k), principal_k)
          }
          return(fit_nest_k)
        })
        out_data <- unlist(out_data, recursive = FALSE)
        out_data <- out_data[names(Z_MI[[j]])]
        if (is.null(out_var)){
          return(list(mean = out_data))
        }else{
          names(out_var) <- principal_group
          return(list(mean = out_data, var = out_var))
        }
      }else{
        if (length(Z_MI[[j]]) != 2){stop('...')}
        if (control$mi_init == 'random'){
          out <- sapply(j_group, FUN=function(i){i[1]})
          out <- lapply(Z_MI[[j]][out], FUN=function(i){
            matrix(rnorm(mi_d_j[j] * ncol(i)), ncol = mi_d_j[j])
          })
          out <- list(mean = out)
        }else{
          out <- init_MI_from_svd(
            data_mi = Z_MI[[j]],
            rank = mi_d_j[j], y = init_y,
            pg_weight = init_w,
            prior_U = 1, prior_V = 1)
          out <- list(mean = list(out$mean_U, out$mean_V),
                      var = list(out$var_U, out$var_V))
        }
        return(out)
      }
    })

    vi_mi_mean <- lapply(init_vi_mi, `[[`, "mean")
    save_init <<- vi_mi_mean
    
    Z_MI_hier_mapping <- mapply(Z_MI_grouping, lapply(Z_MI_attr, `[[`, 'levels'), 
          lapply(Z_MI_attr, `[[`, "id"), SIMPLIFY = FALSE, FUN=function(group_j, levels_j, id_j){
      mapply(group_j, SIMPLIFY = FALSE, FUN=function(group_l){
        if (length(group_l) == 1){
          return(NULL)
        }else{
          warning('clean up principal mapping...')
          p_level <- levels_j[[group_l[1]]]
          non_principal <- group_l[-1]
          pos_princip <- match(1:length(p_level), id_j[,group_l[1]])
          mapping_non_principal <- lapply(non_principal, FUN=function(v){
            sparseMatrix(i = 1:length(p_level),
                         j = id_j[,v][pos_princip],
                         x = 1,
                         dims = c(length(p_level), length(levels_j[[v]]))
            )
          })
          names(mapping_non_principal) <- non_principal
          return(mapping_non_principal)
        }
      })
    })
    
    vi_mi_var <- mapply(Z_MI_attr, init_vi_mi, SIMPLIFY = FALSE, FUN=function(i,v_i){
      o <- lapply(i$storage, FUN=function(j){
        drop0(matrix(data = 0 * as.vector(diag(ncol(j))), 
                     byrow = T,
                     nrow = nrow(j), ncol = ncol(j)^2))
      })
      for (v in names(v_i$var)){
        o[[v]] <- v_i$var[[v]]
      }
      return(o)
    })
    
    # Get number of parameters
    vi_mi_sigma_alpha_nu <- lapply(
      Z_MI_attr, FUN=function(i){
        if (mi_prior_type %in% c('centered')){
          nrow(i$storage[[2]])
        }else if (mi_prior_type %in% c('shared')){
          sum(sapply(i$storage, nrow))
        }else if (mi_prior_type %in% c('separate')){
          sapply(i$storage, nrow)
        }
      })
    
    
    if (do_huangwand_mi){

      vi_mi_a_nu_jp <- rep(2, N_MI)
      names(vi_mi_a_nu_jp) <- names_of_MI
      mi_d_j <- sapply(Z_MI_attr, FUN=function(i){ncol(i$storage[[1]])})
      vi_mi_diag <- lapply(mi_d_j, FUN=function(j){
        cumsum(c(1, rep(j + 1, j-1)))
      })

      vi_mi_diag <- lapply(mi_d_j, FUN=function(j){
        cumsum(c(1, rep(j + 1, j-1)))
      })
      
      if (mi_prior_type %in% c('shared', 'centered')){
        vi_mi_a_APRIOR_jp <- lapply(mi_d_j, FUN=function(i){rep(5, i)})
        vi_mi_a_a_jp <- mapply(mi_d_j, vi_mi_a_nu_jp, SIMPLIFY = FALSE, 
                               FUN=function(i,nu){1/2 * (nu + rep(i, i))})
        vi_mi_a_b_jp <- lapply(vi_mi_a_a_jp, FUN=function(i){1/i^2})
        
      }else{
        vi_mi_a_APRIOR_jp <- mapply(mi_d_j, Z_MI_hier_size, lapply(Z_MI, names), SIMPLIFY = FALSE,
            FUN=function(i, size_i, names_i){
              setNames(lapply(1:size_i, FUN=function(k){rep(5, i)}), names_i)
            })
        vi_mi_a_a_jp <- mapply(mi_d_j, Z_MI_hier_size, vi_mi_a_nu_jp, lapply(Z_MI, names), SIMPLIFY = FALSE, 
         FUN=function(i,size_i, nu, names_i){
           setNames(lapply(1:size_i,FUN=function(k){1/2 * (nu + rep(i, i))}), names_i)})
        vi_mi_a_b_jp <- lapply(vi_mi_a_a_jp, FUN=function(i){lapply(i, FUN=function(k){1/k^2})})

      }
      
      vi_mi_sigma_alpha <- get_bilinear_outer(vi_mi_mean, vi_mi_var, 
                                              vi_mi_diag, mi_prior_type)
      vi_mi_sigma_alpha <- mapply(vi_mi_sigma_alpha, vi_mi_a_a_jp, 
        vi_mi_a_b_jp, vi_mi_a_nu_jp, SIMPLIFY = FALSE, FUN = function(i, tilde.a, tilde.b, nu) {
          if (mi_prior_type %in% c('shared', 'centered')){
            i + sparseMatrix(i = seq_len(nrow(i)), j = seq_len(nrow(i)), x = 1)
          }else{
            mapply(i, tilde.a, tilde.b, nu, FUN=function(i_l, tilde.a_l, tilde.b_l, nu_l){
              i_l + sparseMatrix(i = seq_len(nrow(i_l)), j = seq_len(nrow(i_l)), x = 1)
            })
          }
        })
      #Update a_{j,p}
      mi_diag_Einv_sigma <- mapply(vi_mi_sigma_alpha, 
       vi_mi_sigma_alpha_nu, mi_d_j, SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
         if (mi_prior_type %in% c('centered', 'shared')){
           inv_phi <- solve(phi)
           sigma.inv <- nu * inv_phi
           return(diag(sigma.inv))
         }else{
           mapply(phi, nu, d, SIMPLIFY = FALSE, FUN=function(phi_l, nu_l, d_l){
             inv_phi_l <- solve(phi_l)
             sigma.inv <- nu_l * inv_phi_l
             return(diag(sigma.inv))
             
           })
         }
       })
      
      vi_mi_a_b_jp <- mapply(vi_mi_a_nu_jp, vi_mi_a_APRIOR_jp, mi_diag_Einv_sigma,
         SIMPLIFY = FALSE,
         FUN=function(nu, APRIOR, diag_j){
           if (mi_prior_type %in% c('centered', 'shared')){
             1/APRIOR^2 + nu * diag_j
           }else{
             mapply(APRIOR, diag_j, SIMPLIFY = FALSE, FUN=function(APRIOR_l, diag_j_l){
               1/APRIOR_l^2 + nu * diag_j_l
             })
           }
         })
      
      mi_iw_prior_constant <- mapply(vi_mi_a_nu_jp, mi_d_j,
         FUN = function(nu, d) {
           nu <- nu + d - 1
           return(- (nu * d) / 2 * log(2) - multi_lgamma(a = nu / 2, p = d))
         }
      )
      if (mi_prior_type %in% 'separate'){
        mi_iw_prior_constant <- mi_iw_prior_constant * Z_MI_hier_size
      }
      
      vi_mi_sigma_alpha_nu <- mapply(vi_mi_sigma_alpha_nu, vi_mi_a_nu_jp, mi_d_j, SIMPLIFY = FALSE,
        FUN=function(l,m,d){l + m + d - 1})
      
      if (partial_fix & mi_prior_type %in% 'separate' | mi_prior_variance == 'fixed'){
        vi_mi_sigma_alpha <- mapply(vi_mi_sigma_alpha, Z_MI_grouping, 
                vi_mi_sigma_alpha_nu, mi_d_j, lapply(Z_MI_attr, `[[`, 'xt'),
          SIMPLIFY = FALSE, FUN=function(sigma_j, group_j, nu_j, d_j, xt_j){
            sigma_j[[group_j[[1]][1]]] <- Diagonal(n=d_j) * xt_j[[group_j[[1]][1]]] * nu_j[group_j[[1]][1]]
            sigma_j[[group_j[[2]][1]]] <- Diagonal(n=d_j) * xt_j[[group_j[[2]][1]]] * nu_j[group_j[[2]][1]]
            return(sigma_j)
          })
      }
      
      mi_prior_sigma_alpha_nu <- NULL
      mi_prior_sigma_alpha_phi <- NULL
    }else{
     
      if (mi_prior_variance == 'mean_exists'){
        mi_prior_sigma_alpha_nu <- mi_d_j + 1
        mi_prior_sigma_alpha_phi <- lapply(mi_d_j, FUN = function(i) {
          diag(x = 1, nrow = i, ncol = i)
        })
      }else if (mi_prior_variance == 'fixed'){
        mi_prior_sigma_alpha_nu <- mi_d_j + 1
        mi_prior_sigma_alpha_phi <- lapply(mi_d_j, FUN = function(i) {
          diag(x = 1, nrow = i, ncol = i)
        })
        freeze_mi_var <- TRUE        
      }else{
        stop('invalid value for mi_prior_variance')
      }
      
      if (mi_prior_type %in% c('separate')){
        mi_prior_sigma_alpha_phi <- mapply(mi_prior_sigma_alpha_phi, Z_MI_hier_size, lapply(Z_MI, names),
          SIMPLIFY = FALSE, FUN=function(i, size_i, names_i){
            out <- setNames(lapply(1:size_i, FUN=function(j){i}), names_i)
            return(out)
        })
      }
      
      vi_mi_diag <- lapply(mi_d_j, FUN=function(j){
        cumsum(c(1, rep(j + 1, j-1)))
      })
      vi_mi_sigma_alpha <- get_bilinear_outer(vi_mi_mean, vi_mi_var, vi_mi_diag, mi_prior_type)
      
      mi_iw_prior_constant <- mapply(mi_prior_sigma_alpha_nu, mi_prior_sigma_alpha_phi, SIMPLIFY = FALSE,
         FUN = function(nu, Phi) {
           if (mi_prior_type %in% c('shared', 'centered')){
             if (nu <= (ncol(Phi) - 1)) {
               return(0)
             } else {
               return(make_log_invwishart_constant(nu, Phi))
             }
          }else if (mi_prior_type %in% c('separate')){
               out <- mapply(nu, Phi, FUN=function(nu_l, Phi_l){
                 if (nu_l <= (ncol(Phi_l) - 1)) {
                   return(0)
                 } else {
                   return(make_log_invwishart_constant(nu_l, Phi_l))
                 }
               })
               return(out)
          }else{stop('invalid mi_prior_type')}
      })
      if (mi_prior_type %in% c('centered', 'shared')){
        vi_mi_sigma_alpha_nu <- vi_mi_sigma_alpha_nu + mi_prior_sigma_alpha_nu
      }else{
        vi_mi_sigma_alpha_nu <- mapply(vi_mi_sigma_alpha_nu, mi_prior_sigma_alpha_nu, SIMPLIFY = FALSE, FUN=function(i,j){
          i + j
        })
      }
      
      if (mi_prior_variance == 'fixed'){
        vi_mi_sigma_alpha <- mapply(vi_mi_sigma_alpha_nu, 
                                    vi_mi_sigma_alpha,
                                    mi_prior_sigma_alpha_phi, Z_MI_grouping,
                                    lapply(Z_MI, names), lapply(Z_MI_attr, `[[`, "xt"),
                                    SIMPLIFY = FALSE, FUN = function(nu, i, j, group_j, names_j, xt_j) {
          if (mi_prior_type %in% c('shared', 'centered')){
            j * 1/nu
          }else if (mi_prior_type %in% c('separate')){
            if (all(lengths(group_j) == 1)){
              return(mapply(nu, i, j, group_j, SIMPLIFY = FALSE, 
                  FUN=function(nu_l, a_l, b_l, g_l){xt_j[[g_l]] * nu_l * b_l}))
            }else{
              out <- as.list(rep(NA, length(names_j)))
              out <- setNames(out, names_j)
              
              
              out_group <- mapply(group_j, SIMPLIFY = FALSE, FUN=function(k){
                if (!all(k %in% names(xt_j))){
                  stop(paste0('xt must include all terms.'))
                }
                return(
                  mapply(nu[k], i[k], j[k], xt_j[k], SIMPLIFY = FALSE, FUN=function(nu_l, a_l, b_l, w_l){w_l * nu_l * b_l})
                )
              })
              out[group_j[[1]]] <- out_group[[1]]
              out[group_j[[2]]] <- out_group[[2]]
              return(out)
            }
          }else{stop('invalid mi_prior_type')}
        })
      }else{
        vi_mi_sigma_alpha <- mapply(vi_mi_sigma_alpha, mi_prior_sigma_alpha_phi, SIMPLIFY = FALSE, FUN = function(i, j) {
          if (mi_prior_type %in% c('shared', 'centered')){
            i + j
          }else if (mi_prior_type %in% c('separate')){
            return(mapply(i, j, SIMPLIFY = FALSE, FUN=function(a,b){a+b}))
          }else{stop('invalid mi_prior_type')}
        })
      }
      
      vi_mi_a_b_jp <- vi_mi_a_a_jp <- vi_mi_a_APRIOR_jp <- vi_mi_a_nu_jp <- NULL 
    }

    vi_mi_sigma_outer_alpha <- get_bilinear_outer(
      vi_mi_mean, vi_mi_var, vi_mi_diag, mi_prior_type)
    vi_mi_lndet <- lapply(Z_MI, FUN=function(i){
      setNames(rep(0, length(i)), names(i))
    })

    lagged_vi_mi_mean <- lapply(Z_MI_attr, FUN=function(i){
      lapply(i$storage, FUN=function(j){
        matrix(data = -Inf, nrow = nrow(j), ncol = ncol(j))
      })
    })
    lagged_vi_mi_var <- lapply(Z_MI_attr, 
                               FUN=function(i){
                                 lapply(i$storage, FUN=function(j){
                                   drop0(matrix(data = -Inf, nrow = nrow(j), ncol = ncol(j)^2))
                                 })
                               })
    lagged_vi_mi_sigma_alpha <- mapply(Z_MI_attr, Z_MI_hier_size, SIMPLIFY = FALSE, FUN=function(i,s_i){
      array(-Inf, rep(ncol(i$storage[[1]]),s_i))
    })
    lagged_vi_mi_lndet <- lapply(vi_mi_lndet, FUN=function(i){
      setNames(rep(-Inf, length(i)), names(i))
    })
    
  }else{
    vi_mi_sigma_alpha <- NULL
    vi_mi_sigma_alpha_nu <- NULL
    vi_mi_sigma_outer_alpha <- NULL
    vi_mi_mean <- NULL
    vi_mi_var <- NULL
    vi_mi_lndet <- NULL
    mi_d_j <- NULL
    mi_g_j <- NULL
    vi_mi_a_a_jp <- vi_mi_a_b_jp <- vi_mi_a_APRIOR_jp <- vi_mi_a_nu_jp <- NULL
    mi_prior_sigma_alpha_nu <- mi_prior_sigma_alpha_phi <- NULL
    mi_iw_prior_constant <- NULL
    mi_prior_type <- NULL
    Z_MI <- NULL
    vi_mi_diag <- NULL
  }
  
  # Create mapping for this to allow sparse implementations.
  if (any_RE){
    mapping_sigma_alpha <- make_mapping_alpha(vi_sigma_alpha)
  }

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
  
  spline_REs <- grepl(names(d_j), pattern='^spline-')
  zeromat_beta <- drop0(Diagonal(x = rep(0, ncol(X))))
  stationary_rho <- do.call('c', lapply(d_j[!spline_REs], FUN=function(i){as.vector(diag(x = i))}))
  
  if (parameter_expansion %in%  c("translation", "diagonal") & any_Mprime & any(!spline_REs)) {
    
    if (do_timing){
      tic('Build PX R Terms')
    }
    
    
    
    nonspline_positions <- sort(unlist(outer_alpha_RE_positions[!spline_REs]))
    
    size_splines <- sum((d_j * g_j)[spline_REs])
    
    est_rho <- stationary_rho
    diag_rho <- which(stationary_rho == 1)
    

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

    id_range <- 1:nrow(Mmap)
    store_re_id <- store_id <- list()
    for (j in 1:(number_of_RE - sum(spline_REs))){
      store_re_id_j <- store_id_j <- list()
      for (jprime in j){
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
    store_assignment_Z <- parsed_RE_groups$factor
    store_levels_Z <- parsed_RE_groups$nl
    
    lookup_Z <- mapply(store_assignment_Z, store_levels_Z, FUN=function(i,nl){
      sparseMatrix(i = 1:length(i), j = i, x = 1, dims = c(length(i), nl))
    })

    rm(parsed_RE_groups, mapping_to_re)
    
    gc()
    if (do_timing){
      toc(quiet = quiet_time, log = T)
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
    if (any_mi){
      n_param_mean <- sapply(Z_MI_attr, FUN=function(i){
        sum(sapply(i$storage, FUN=function(j){prod(dim(j))}))
      })
      n_name_mean <- unlist(mapply(Z_MI_attr, mi_d_j, SIMPLIFY = FALSE, FUN=function(i,dim_i){
        do.call('c', mapply(i$levels, names(i$levels), SIMPLIFY = FALSE, FUN=function(j, n){
          paste(n, as.vector(outer(j, 
                          1:dim_i, FUN=function(x,y){paste(x,y,sep='@')})), sep='@')
        }))
      }))
      n_name_var <- unlist(mapply(names_of_MI, mi_d_j, Z_MI_attr, SIMPLIFY = FALSE, FUN=function(a,b,i){
        prefix <- paste(a, '', 1:b^2, sep = '@')
        do.call('c', mapply(i$levels, names(i$levels), SIMPLIFY = FALSE, FUN=function(j,n){
          paste(n, as.vector(outer(j, 
                          prefix, FUN=function(x,y){paste(x,y, sep='@')})), sep = '@')
        }))
      }))
      
      store_mi_mean <- array(NA, dim = c(iterations, sum(n_param_mean)))
      colnames(store_mi_mean) <- n_name_mean
      store_mi_var <- array(NA, dim = c(iterations, sum(n_param_mean * mi_d_j)))
      colnames(store_mi_var) <- n_name_var
      store_mi_sigma <- array(NA, dim = c(iterations, sum(lengths(lapply(Z_MI_attr, `[[`, 'levels')) * mi_d_j^2)))
      if (do_huangwand_mi){
        store_mi_hw <- array(NA, dim = c(iterations, 
         sum(sapply(vi_mi_a_b_jp, FUN=function(i){sum(lengths(i))}))))
      }
      rm(n_param_mean, n_name_mean)
    }
  }
  if (do_timing) {
    toc(quiet = quiet_time, log = TRUE)
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
  # Create terms needed for MI SQUAREM
  if (any_mi){
    type_invert <- do_SQUAREM & (mi_d_j > 1)
    keep_mi_decomp <- do_SQUAREM
    if (keep_mi_decomp){
      vi_mi_marg_decomp <- lapply(vi_mi_var, FUN=function(i){
        lapply(i, FUN=function(j){
          return(j * NA)
        })
      })
    }else{
      vi_mi_marg_decomp <- NULL
    }
    

    f_add <- function(rho, 
                    mean_k, vc_k, 
                    alpha_mean, mean_l, alpha_ESigma.inv,
                    pg, long_mean_k, long_var_l, simple){
    rho_intercept <- rho[1]
    rho <- rho[-1]
    if (simple){
      demean_k <- sweep(mean_k, MARGIN = 2, FUN = '-', STATS = rho)
      t1 <- -1/2 * sum(rowSums( (demean_k %*% vc_k) * demean_k))
      t2 <- alpha_mean - rho_intercept + mean_l %*% rho
      t2 <- -1/2 * sum(t2^2) * as.numeric(alpha_ESigma.inv)
      
      demean_long_k <- sweep(long_mean_k, MARGIN = 2, FUN = '-', STATS = rho)
      t3 <- -1/2 * sum(pg * rowSums(FS(demean_long_k, demean_long_k) * long_var_l))
      return(t1 + t2 + t3)
    }else{
      rho <- matrix(rho, nrow = ncol(vc_k[[1]]))
      rho <- lapply(1:ncol(rho), FUN=function(i){rho[,i]})
      demean_k <- mapply(mean_k, rho, SIMPLIFY = FALSE, FUN=function(i,j){
        sweep(i, MARGIN = 2, FUN = '-', STATS = j)
      })
      t1 <- -1/2 * sum(mapply(demean_k, vc_k, FUN=function(i,j){sum(rowSums( (i %*% j) * i))}))
      t2 <- alpha_mean - rho_intercept + mean_l %*% Reduce('+', rho)
      t2 <- -1/2 * sum(t2^2) * as.numeric(alpha_ESigma.inv)
      demean_long_k <- mapply(long_mean_k, rho, SIMPLIFY = FALSE, FUN=function(i,j){
        sweep(i,MARGIN=2,FUN='-', STATS = j)
      })
      t3 <- -1/2 * sum(sapply(demean_long_k, FUN=function(i){
        sum(pg * rowSums(FS(i,i) * long_var_l))
      }))
      return(t1 + t2 + t3)
    }
    # out <- as.vector(
    #   # Prior on x_i: -1/2 E[(x_i - mu)^2] E[1/sigma^2_x]
    #   -1/2 * sum( (mean_k - rho)^2) * vc_k +
    #   # Prior on alpha_j: -1/2 E[(alpha_j - mu_INT + E[beta_j] mu)^2] E[1/sigma^2_alpha]
    #   -1/2 * sum( (alpha_mean - rho_intercept +  mean_l * rho)^2) * alpha_ESigma.inv
    #   # Prior on likelihood: -1/2 E[omega_i] Var(beta_j) * (E[x_i] - mu)^2
    #   -1/2 * sum( pg * long_var_l * (long_mean_k - rho)^2)
    # )
    # return(out)
  }
  
    outer_px_size <- lapply(Z_MI_attr, FUN=function(i){sapply(i$storage, nrow)})
    outer_px_add_map <- lapply(Z_MI_attr, FUN=function(i){i$mapping_RE}) 
    outer_px_add_map <- mapply(outer_px_add_map, Z_MI_attr, 
     SIMPLIFY = FALSE,FUN=function(i, attr_i){
       mapply(i, attr_i$levels, names(attr_i$storage), SIMPLIFY = FALSE, 
              FUN=function(j, j_levels, j_term){
                if (!is.null(j)){
                  mi_value <- paste0(j_term, ' @ (Intercept) @ ', j_levels)
                  if (!all(mi_value == j$value)){
                    warning('Alignment Issue')
                    browser()
                  }else{
                    return(list(position = j$position, index = j$index))
                  }
                }
              })
     })

  }
  
  if (debug_px){
    debug_PX_ELBO <- rep(NA, iterations)
  }else{
    debug_PX_ELBO <- NULL
  }
  
  if (freeze_RE & any_mi){
    list2env(init_model, base::environment())
  }
  
  if (INIT_MCMC & any_mi){
    vi_beta_mean[,] <- mean(MCMC_estimates$alpha[,3])
    vi_alpha_mean[,] <- MCMC_estimates$alpha[,3] - mean(MCMC_estimates$alpha[,3])
    diag(vi_alpha_decomp) <- MCMC_estimates$alpha[,2]
    vi_mi_mean[[1]]$bioname[,1] <- MCMC_estimates$bioname[,3]
    vi_mi_mean[[1]]$state[,1] <- MCMC_estimates$state[,3]
    vi_mi_mean[[1]]$party[,1] <- MCMC_estimates$party[,3]
    vi_mi_mean[[1]]$rollnumber[,1] <- MCMC_estimates$beta[,3]
    vi_mi_var[[1]]$bioname[,1] <- MCMC_estimates$bioname[,2]
    vi_mi_var[[1]]$state[,1] <- MCMC_estimates$state[,2]
    vi_mi_var[[1]]$party[,1] <- MCMC_estimates$party[,2]
    vi_mi_var[[1]]$rollnumber[,1] <- MCMC_estimates$beta[,2]
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

    if (family %in% 'linear'){# Ignore Polya-Gamma or Similar Updates
      
      vi_pg_mean <- rep(1, nrow(X))
      diag_vi_pg_mean <- sparseMatrix(i = 1:N, j = 1:N, x = vi_pg_mean)
      
    }else{# Estimate Polya-Gamma or Similar Updates

      if (any_mi){
        bilinear_mean <- get_bilinear_mean(Z_MI, vi_mi_mean, Z_MI_grouping)
        bilinear_var <- get_bilinear_var(Z_MI, vi_mi_mean, vi_mi_var, Z_MI_grouping)
      }else{
        bilinear_mean <- 0
        bilinear_var <- 0
      }
      if (factorization_method %in% c("weak", "collapsed")) {
        # joint_var <- rowSums( (joint.XZ %*% t(vi_joint_decomp))^2 )
        # vi_joint_decomp <<- vi_joint_decomp
        # joint.XZ <<- joint.XZ
        joint_var <- cpp_zVz(Z = joint.XZ, V = as(vi_joint_decomp, "generalMatrix")) 
        if (family == 'negbin'){
          joint_var <- joint_var + vi_r_sigma
        }
        joint_var <- joint_var + bilinear_var
      } else {
        
        # beta_quad <- rowSums((X %*% t(vi_beta_decomp))^2)
        beta_quad <- cpp_dense_zVz(X, as.matrix(vi_beta_decomp))
        alpha_quad <- rowSums((Z %*% t(vi_alpha_decomp))^2)
        joint_var <- beta_quad + alpha_quad
        if (family == 'negbin'){
          joint_var <- joint_var + vi_r_sigma
        }
        joint_var <- joint_var + bilinear_var
      }
      
      vi_pg_c <- sqrt(as.vector(bilinear_mean + X %*% vi_beta_mean + Z %*% vi_alpha_mean - vi_r_mu)^2 + joint_var)
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
        do_huangwand = do_huangwand, do_huangwand_mi = do_huangwand_mi,
        vi_a_a_jp = vi_a_a_jp, vi_a_b_jp = vi_a_b_jp,
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
        # Multiplicative Interaction
        any_RE = any_RE, mi_prior_type = mi_prior_type,
        any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
        vi_mi_diag = vi_mi_diag, 
        vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
        vi_mi_sigma_alpha = vi_mi_sigma_alpha,
        vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
        vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
        vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
        vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
        mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
        mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
        mi_iw_prior_constant = mi_iw_prior_constant,
        mi_d_j = mi_d_j, mi_g_j = mi_g_j
      )
      if (debug_ELBO.1$ELBO < final.ELBO$ELBO - sqrt(.Machine$double.eps)){
        browser()
      }
    }

    loop_RE <- c()
    loop_RE <- c(loop_RE, 're')
    update_mi <- it > MI_STAB
    if (any_mi & update_mi){
      loop_RE <- c(loop_RE, 'mi')
    }
    loop_RE <- sort(loop_RE, decreasing = T)
    print(loop_RE)
    
    vi_alpha_L_nonpermute <- NULL
    vi_beta_L_nonpermute <- NULL
    vi_alpha_LP <- NULL
    vi_beta_LP <- NULL
    
    for (lll in loop_RE){
      
      if (lll == 're'){
        
        if (do_timing) {
          toc(quiet = quiet_time, log = TRUE)
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
        
        if (any_RE){
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
            Tinv <- as(Tinv, "generalMatrix")
          } else {
            Tinv <- lapply(Tinv, FUN = function(i) {
              as(i, "generalMatrix")
            })
          }
          
        }else{
          Tinv <- matrix(nrow=0,ncol=0)
        }
        
        if (do_timing) {
          toc(quiet = quiet_time, log = T)
          tic("Update Beta")
        }
        

        if (any_mi){
          offset_bilinear <- get_bilinear_mean(Z_MI, vi_mi_mean, Z_MI_grouping)
        }else{
          offset_bilinear <- 0
        }
        
        if (!(freeze_RE & any_mi)){
          
          if (factorization_method == "weak") {
            ## Update <beta, alpha> jointly
            chol.update.joint <- LinRegChol(
              X = joint.XZ, omega = diag_vi_pg_mean,
              prior_precision = bdiag(zero_mat, Tinv),
              y = s + vi_pg_mean * (vi_r_mu - offset_bilinear)
            )
            Pmatrix <- sparseMatrix(i = 1:ncol(joint.XZ), j = 1 + chol.update.joint$Pindex, x = 1)
            
            vi_joint_L_nonpermute <- drop0(solve(chol.update.joint$origL))
            vi_joint_LP <- Pmatrix
            vi_joint_decomp <- vi_joint_L_nonpermute %*% t(vi_joint_LP)
            
            if (ncol(X) > 0){
              vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
              vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))
            }else{
              vi_alpha_mean <- Matrix(chol.update.joint$mean, dimnames = list(fmt_names_Z, NULL))
            }
            
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
                                   t(M) %*% (s - diag_vi_pg_mean %*% (X %*% beta_hat + offset_bilinear))
            )
            vi_beta_mean <- beta_hat - P %*% vi_alpha_mean
            
            sqrt_pg_weights <- Diagonal(x = sqrt(vi_pg_mean))
            
            for (j in 1:number_of_RE) {
              index_j <- cyclical_pos[[j]]
              M_j <- as(M[, index_j, drop = F], 'generalMatrix')
              prec_j <- crossprod(sqrt_pg_weights %*% M_j) + Tinv[[j]]
              
              chol_var_j <- solve(t(chol(prec_j)))
              running_log_det_alpha_var[j] <- 2 * sum(log(diag(chol_var_j)))
              
              vi_alpha_decomp <- fast_insert(
                  A = drop0(chol_var_j), 
                  B = vi_alpha_decomp, 
                  index = index_j)
              # vi_alpha_decomp[index_j, index_j] <- drop0(chol_var_j)
              # as(
              #   as(chol_var_j, "generalMatrix"), "TsparseMatrix"
              # )
            }
            vi_alpha_L_nonpermute <- vi_alpha_decomp
            vi_alpha_LP <- Diagonal(n = nrow(vi_alpha_L_nonpermute))
            vi_alpha_decomp <- vi_alpha_L_nonpermute  %*% t(vi_alpha_LP)
            vi_alpha_decomp <- drop0(vi_alpha_decomp)
            vi_alpha_decomp <- as(vi_alpha_decomp, 'generalMatrix')
            
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
              
              browser() # Use LinRegChol_fe
              chol.update.beta <- LinRegChol(
                X = as(X, "sparseMatrix"), omega = diag_vi_pg_mean, prior_precision = zero_mat,
                y = as.vector(s - diag_vi_pg_mean %*% (offset_bilinear + Z %*% vi_alpha_mean))
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
                y = as.vector(s - diag_vi_pg_mean %*% (X %*% vi_beta_mean + offset_bilinear))
              )
              Pmatrix <- sparseMatrix(i = 1:p.Z, j = 1 + chol.update.alpha$Pindex, x = 1)
              
              vi_alpha_L_nonpermute <- drop0(solve(chol.update.alpha$origL))
              vi_alpha_LP <- Pmatrix
              vi_alpha_decomp <- vi_alpha_L_nonpermute  %*% t(vi_alpha_LP)
              vi_alpha_decomp <- drop0(vi_alpha_decomp)
              vi_alpha_decomp <- as(vi_alpha_decomp, 'generalMatrix')
              vi_alpha_mean <- chol.update.alpha$mean
              log_det_alpha_var <- -2 * sum(log(diag(chol.update.alpha$origL)))
              
              vi_beta_mean <- Matrix(vi_beta_mean, dimnames = list(colnames(X), NULL))
              vi_alpha_mean <- Matrix(vi_alpha_mean, dimnames = list(fmt_names_Z, NULL))
            } else if (linpred_method == "joint") {
              
              joint.XZ <- cbind(X, Z)
              
              chol.update.joint <- solve(Matrix::Cholesky(  
                crossprod(Diagonal(x = sqrt(vi_pg_mean)) %*% joint.XZ) + 
                  bdiag(zero_mat, bdiag(Tinv)) ),
                t(joint.XZ) %*% (s + vi_pg_mean * (vi_r_mu - offset_bilinear)) )
              
              if (ncol(X) > 0){
                vi_beta_mean <- Matrix(chol.update.joint[1:p.X,], dimnames = list(colnames(X), NULL))
                vi_alpha_mean <- Matrix(chol.update.joint[-1:-p.X,], dimnames = list(fmt_names_Z, NULL))
              }else{
                vi_alpha_mean <- Matrix(chol.update.joint, dimnames = list(fmt_names_Z, NULL))
              }
              
              # chol.update.joint <- LinRegChol(
              #   X = joint.XZ, omega = diag_vi_pg_mean, prior_precision = bdiag(zero_mat, Tinv),
              #   y = s + vi_pg_mean * (vi_r_mu - offset_bilinear),
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
                y = s + vi_pg_mean * (vi_r_mu - offset_bilinear)
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
            if (any_RE){
              vi_alpha_decomp <- sparseMatrix(i = 1, j = 1, x = 0, dims = rep(p.Z, 2))
            }else{
              vi_alpha_decomp <- matrix(nrow=0,ncol=0)
            }
            
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
                t(joint.XZ) %*% (s + vi_pg_mean * (vi_r_mu - offset_bilinear)) )
              
              if (ncol(X) > 0){
                vi_beta_mean <- Matrix(chol.update.joint[1:p.X,], dimnames = list(colnames(X), NULL))
                vi_alpha_mean <- Matrix(chol.update.joint[-1:-p.X,], dimnames = list(fmt_names_Z, NULL))
              }else{
                vi_alpha_mean <- Matrix(chol.update.joint, dimnames = list(fmt_names_Z, NULL))
              }
              
              # chol.update.joint <- LinRegChol(X = joint.XZ,
              #   omega = diag_vi_pg_mean,
              #   prior_precision = bdiag(zero_mat, bdiag(Tinv)),
              #   y = s + vi_pg_mean * (vi_r_mu - offset_bilinear),
              #   save_chol = FALSE)
              # vi_beta_mean <- Matrix(chol.update.joint$mean[1:p.X], dimnames = list(colnames(X), NULL))
              # vi_alpha_mean <- Matrix(chol.update.joint$mean[-1:-p.X], dimnames = list(fmt_names_Z, NULL))
              
              if (do_timing) {
                toc(quiet = quiet_time, log = T)
                tic("ux_var")
              }
              
              vi_beta_decomp <- solve(t(chol(as.matrix(t(X) %*% diag_vi_pg_mean %*% X))))
              
              vi_beta_L_nonpermute <- vi_beta_decomp
              vi_beta_LP <- Diagonal(n = nrow(vi_beta_decomp))
              log_det_beta_var <- 2 * sum(log(diag(vi_beta_decomp)))
              #-log(det(t(X) %*% diag_vi_pg_mean %*% X))
              
              running_log_det_alpha_var <- rep(NA, number_of_RE)
              
              if (any_RE){
                
                for (j in 1:number_of_RE) {
                  index_j <- cyclical_pos[[j]]
                  Z_j <- Z[, index_j, drop = F]
                  prec_j <- crossprod(sqrt_pg_weights %*% Z_j) + Tinv[[j]]
                  
                  chol_var_j <- solve(t(chol(prec_j)))
                  running_log_det_alpha_var[j] <- 2 * sum(log(diag(chol_var_j)))
                  
                  vi_alpha_decomp <- fast_insert(A = drop0(chol_var_j), B = vi_alpha_decomp, 
                                                 index = index_j)
                  # vi_alpha_decomp[index_j, index_j] <- drop0(chol_var_j)

                }
                
              }
              
              vi_alpha_L_nonpermute <- vi_alpha_decomp
              vi_alpha_LP <- Diagonal(n = nrow(vi_alpha_L_nonpermute))
              
              log_det_alpha_var <- sum(running_log_det_alpha_var)
              
              if (do_timing){
                toc(quiet = quiet_time, log = T)
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
                
                if (any_mi){stop('setup for mi')}
                bind_lhs_j[[j]] <- term_j
                bind_rhs_j[[j]] <- mod_j %*% t(Z_j) %*% (s)
                
                running_log_det_alpha_var[j] <- 2 * sum(log(diag(chol_var_j)))
                vi_alpha_decomp <- fast_insert(
                  A = drop0(chol_var_j),
                  B = vi_alpha_decomp, 
                  index = index_j)
                # vi_alpha_decomp[index_j, index_j] <- drop0(chol_var_j)
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
              
              if (any_mi){stop('setup for mi')}
              
              bind_lhs_j <- drop0(rbind(bind_lhs_j, cbind(vi_beta_var %*% t(X) %*% diag_vi_pg_mean %*% Z, Diagonal(n = ncol(X)))))
              bind_rhs_j <- rbind(bind_rhs_j, vi_beta_var %*% t(X) %*% s)
              #
              bind_solution <- solve(bind_lhs_j) %*% bind_rhs_j
              # print(cbind(bind_solution, rbind(vi_alpha_mean, vi_beta_mean)))
              #
              vi_beta_mean <- Matrix(bind_solution[-1:-ncol(Z)], dimnames = list(colnames(X), NULL))
              vi_alpha_mean <- Matrix(bind_solution[1:ncol(Z)], dimnames = list(fmt_names_Z, NULL))
            } else if (linpred_method == "cyclical") {
              
              
              if (any_RE){
                
                for (j in 1:number_of_RE) {
                  index_j <- cyclical_pos[[j]]
                  Z_j <- Z[, index_j, drop = F]
                  Z_negj <- Z[, -index_j, drop = F]
                  
                  chol.j <- LinRegChol(
                    X = Z_j, omega = diag_vi_pg_mean, prior_precision = Tinv[[j]],
                    y = as.vector(s + vi_pg_mean * vi_r_mu - diag_vi_pg_mean %*% (offset_bilinear + X %*% vi_beta_mean + Z_negj %*% vi_alpha_mean[-index_j]))
                  )
                  vi_alpha_mean[index_j] <- chol.j$mean
                  
                  Pmatrix <- sparseMatrix(
                    i = 1:ncol(Z_j), 
                    j = 1 + chol.j$Pindex,
                    x = 1)
                  
                  running_log_det_alpha_var[j] <- -2 * sum(log(diag(chol.j$origL)))
                  vi_alpha_decomp <- fast_insert(A = solve(chol.j$origL) %*% t(Pmatrix), B = vi_alpha_decomp, index = index_j)
                  # vi_alpha_decomp[index_j, index_j] <- solve(chol.j$origL) %*% t(Pmatrix)
                }
                
                vi_alpha_L_nonpermute <- vi_alpha_decomp
                vi_alpha_LP <- Diagonal(n = ncol(vi_alpha_L_nonpermute))
                # vi_alpha_decomp <- bdiag(vi_alpha_decomp)
                log_det_alpha_var <- sum(running_log_det_alpha_var)
                
              }
              
              chol.update.beta <- LinRegChol_fe(
                X = X, omega = diag_vi_pg_mean, 
                y = as.vector((s + vi_pg_mean * (vi_r_mu - offset_bilinear) - diag_vi_pg_mean %*% Z %*% vi_alpha_mean))
              )
              Pmatrix <- sparseMatrix(i = 1 + chol.update.beta$Pindex, j = 1:p.X, x = 1)
              
              # old <- LinRegChol(
              #   X = as(X, "sparseMatrix"), omega = diag_vi_pg_mean, prior_precision = zero_mat,
              #   y = as.vector(s + vi_pg_mean * (vi_r_mu - offset_bilinear) - diag_vi_pg_mean %*% Z %*% vi_alpha_mean)
              # )
              # Pmatrix <- sparseMatrix(i = 1:p.X, j = 1 + chol.update.beta$Pindex, x = 1)
              
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
          
        }else{
          # Do not update q(beta, alpha) if freezing RE      
        }
        
        ###
        # Update \Sigma_j
        ###
        
        if (any_RE){
          
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
            
          }
          
        }
        
      }else if (lll == 'mi'){

        # Update parameters for bilinear terms
        if (any_mi){
          
          # Get initial linear predictor
          running_billinear_lp <- 
            as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean) +
            get_bilinear_mean(Z_MI, vi_mi_mean, Z_MI_grouping)
          
          mi_inv_mapping_alpha <- mapply(
            vi_mi_sigma_alpha_nu, vi_mi_sigma_alpha,
            SIMPLIFY = FALSE, FUN = function(a, b) {
              if (mi_prior_type %in% c('centered', 'shared')){
                a * solve(b)
              }else{
                mapply(a, b, SIMPLIFY = FALSE, FUN=function(a_l, b_l){
                  a_l * solve(b_l)
                })
              }
            }
          )
          
          if (bi_method == 'joint_direct'){
            if (any_RE){
              Tinv <- mapply(vi_sigma_alpha_nu, lapply(vi_sigma_alpha, solve),
                             SIMPLIFY = FALSE, FUN = function(a, b) {
                               a * b
                             }
              )
              Tinv <- make_mapping_alpha(Tinv)
              Tinv <- prepare_T(
                mapping = Tinv, levels_per_RE = g_j, num_REs = number_of_RE,
                variables_per_RE = d_j, running_per_RE = breaks_for_RE, cyclical = F
              )
              Tinv <- bdiag(Diagonal(x=rep(0,ncol(X))), Tinv)
            }else{
              Tinv <- Diagonal(x=rep(0, ncol(X)))
            }
            
            if (it %in% c(1, MI_STAB, MI_STAB + 1)){
              XZ <- cbind(X,Z)
            }
          }

          if (exists('inner_vi_lndet')){
            old_lndet_all <- vi_mi_lndet
            temp_vi_lndet_all <- vi_mi_lndet
          }else{
            temp_vi_lndet_all <- lapply(Z_MI, FUN=function(i){
              out <- rep(0, length(i))
              names(out) <- names(i)
              return(out)
            })
          }
          vi_mi_lndet <- lapply(Z_MI, FUN=function(i){
            setNames(rep(NA, length(i)), names(i))
          })
          
          for (j in 1:N_MI){
            
            if (exists('inner_vi_lndet')){
              old_lndet <- old_lndet_all[[j]]
              temp_vi_lndet <- temp_vi_lndet_all[[j]]
            }else{
              temp_vi_lndet <- rep(0, length(Z_MI[[j]]))
              names(temp_vi_lndet) <- names(Z_MI[[j]])
            }
            
            Z_MI_grouping_j <- Z_MI_grouping[[j]]
            if (length(Z_MI_grouping_j) != 2){stop('...')}
            
            offset_j <- running_billinear_lp - 
              get_bilinear_mean(Z_MI[j], vi_mi_mean[j], list(Z_MI_grouping_j))
            
            inner_vi_lndet <- setNames(rep(NA, length(Z_MI[[j]])), names(Z_MI[[j]]))
            
            if (mi_prior_type %in% c('centered', 'shared')){
              mi_j_vec_prior <- matrix(as.vector(mi_inv_mapping_alpha[[j]]))
            }else{
              mi_j_vec_prior <- lapply(mi_inv_mapping_alpha[[j]], FUN=function(l){
                matrix(as.vector(l))
              })
            }
            type_invert_j <- type_invert[j]
            if (it == 1 & MAX_K > 1){
              warning("max K > 1")
            }
            if (MAX_K > 0){
              init_mi <- vi_mi_mean[[j]]
            }
            position_j <- outer_px_add_map[[j]]
            simple_j <- all(lengths(Z_MI_grouping_j) == 1)        
            
            if (bi_method == 'joint'){
              
              long_var_j <- lapply(1:2, FUN=function(k){
                Reduce('+', mapply(Z_MI[[j]][Z_MI_grouping_j[[k]]], vi_mi_var[[j]][Z_MI_grouping_j[[k]]], FUN=function(i,j){
                  i %*% j
                }))
              })
              
              ridge_j <- mi_inv_mapping_alpha[[j]]
              dim_j <- mi_d_j[j]
              grouping_j <- Z_MI_grouping_j
              data_j <- Z_MI[[j]]
              recons_j <- lengths(Z_MI_attr[[j]]$levels)
              recons_j <- rep(names(recons_j), recons_j * dim_j)
              omega <- vi_pg_mean
              pos_d <- lapply(1:dim_j, FUN=function(d){dim_j * (d-1) + seq(dim_j)})
              init_par <- unlist(lapply(grouping_j, FUN=function(i){
                unlist(lapply(i, FUN=function(o){as.vector(vi_mi_mean[[j]][[o]])}))
              }))
              
              fit_lbfgs <- lbfgs::lbfgs(
                max_iterations = 10,
                call_eval = f_lbfgs_mean, 
                call_grad = gr_lbfgs_mean,
                recons_j = recons_j,
                grouping_j = grouping_j, pos_d = pos_d,
                s = s, omega = omega, offset_j = offset_j,
                long_var_j = long_var_j, ridge_j = ridge_j,
                data_j = data_j, dim_j = dim_j,
                vars = init_par
              )

              vi_mi_mean[[j]] <- reformat_mean_lbfgs(fit_lbfgs$par, dim_j, recons_j, names(data_j))
              
            }
      
            if (it <= MI_VEM_THRESH){
              VEM_ONLY <- TRUE
            }else{
              VEM_ONLY <- FALSE
            }
            if (VEM_ONLY){
              message(paste0('VEM_ONLY: ', VEM_ONLY))
              VEM_scalar <- 0
            }else{
              VEM_scalar <- 1
            }
            
            g <- function(){calculate_ELBO(family = family,
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
                           vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
                           # Multiplicative Interaction
                           do_huangwand_mi = do_huangwand_mi,
                           any_RE = any_RE, mi_prior_type = mi_prior_type,
                           any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
                           vi_mi_diag = vi_mi_diag, 
                           vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
                           vi_mi_sigma_alpha = vi_mi_sigma_alpha,
                           vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
                           vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = temp_vi_lndet,
                           vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
                           vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
                           mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
                           mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
                           mi_iw_prior_constant = mi_iw_prior_constant,
                           mi_d_j = mi_d_j, mi_g_j = mi_g_j
            )}

            hhh <- function(xxx){
              xxx_oa <- get_bilinear_outer(
                xxx, vi_mi_var, vi_mi_diag, mi_prior_type)
              return(calculate_ELBO(family = family,
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
                                           vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
                                           # Multiplicative Interaction
                                           do_huangwand_mi = do_huangwand_mi,
                                           any_RE = any_RE, mi_prior_type = mi_prior_type,
                                           any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
                                           vi_mi_diag = vi_mi_diag, 
                                           vi_mi_sigma_outer_alpha = xxx_oa,
                                           vi_mi_sigma_alpha = vi_mi_sigma_alpha,
                                           vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
                                           vi_mi_mean = xxx, vi_mi_var = vi_mi_var, vi_mi_lndet = temp_vi_lndet,
                                           vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
                                           vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
                                           mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
                                           mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
                                           mi_iw_prior_constant = mi_iw_prior_constant,
                                           mi_d_j = mi_d_j, mi_g_j = mi_g_j
            ))}
            
            for (k in 1:MAX_K){
              ELBO_start <- g()
              # Freezing v, update u
              if (simple_j){
                RFSmean <- Z_MI[[j]][[2]] %*% FS(vi_mi_mean[[j]][[2]], vi_mi_mean[[j]][[2]])
                Rmean <- Z_MI[[j]][[2]] %*% vi_mi_mean[[j]][[2]]
                Rvar <- Z_MI[[j]][[2]] %*% vi_mi_var[[j]][[2]]
              }else{
                
                Rmean <- Reduce('+', mapply(
                  Z_MI[[j]][Z_MI_grouping_j[[2]]],
                  vi_mi_mean[[j]][Z_MI_grouping_j[[2]]],
                  SIMPLIFY = FALSE,
                  FUN=function(i,j){
                    i %*% j
                  })) 
                
                # Can this be improved for speed?
                RFSmean <- FS(Rmean, Rmean)
                
                Rvar <- Reduce('+', mapply(
                  Z_MI[[j]][Z_MI_grouping_j[[2]]],
                  vi_mi_var[[j]][Z_MI_grouping_j[[2]]],
                  SIMPLIFY = FALSE,
                  FUN=function(i,j){
                    i %*% j
                  }))
              }
              if (mi_prior_type %in% c('centered')){
                prior_u <- matrix(as.vector(Diagonal(n = ncol(Rmean))))
              }else if (mi_prior_type %in% c('shared')){
                prior_u <- mi_j_vec_prior
              }else if (mi_prior_type %in% c('separate')){
                if (simple_j){
                  prior_u <- mi_j_vec_prior[[1]]
                }else{
                  prior_u <- mi_j_vec_prior[Z_MI_grouping_j[[1]]]
                }
              }
              
              if (!simple_j){
                if (bi_method %in% c('joint', 'rowwise')){
                  
                  update_u <- bilinear_update(
                    Z_MI = Z_MI[[j]][Z_MI_grouping_j[[1]]],
                    Z_hier_mapping = Z_MI_hier_mapping[[j]][[1]],
                    prior_mi = prior_u, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                    diag_vi_pg_mean = diag_vi_pg_mean, 
                    RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                    offset = offset_j, it = it,
                    return_chol = type_invert_j & (k == MAX_K),
                    method = 'sparse_direct'
                  )
                  
                  # long_var_j <- lapply(1:2, FUN=function(k){
                  #   Reduce('+', mapply(Z_MI[[j]][Z_MI_grouping_j[[k]]], vi_mi_var[[j]][Z_MI_grouping_j[[k]]], FUN=function(i,j){
                  #     i %*% j
                  #   }))
                  # })
                  # ridge_j <- mi_inv_mapping_alpha[[j]]
                  # dim_j <- mi_d_j[j]
                  # grouping_j <- Z_MI_grouping_j
                  # data_j <- Z_MI[[j]]
                  # recons_j <- lengths(Z_MI_attr[[j]]$levels)
                  # recons_j <- rep(names(recons_j), recons_j * dim_j)
                  # omega <- vi_pg_mean
                  # pos_d <- lapply(1:dim_j, FUN=function(d){dim_j * (d-1) + seq(dim_j)})
                  # 
                  # val_1 <- fff(mean_j = vi_mi_mean[[j]],
                  #   recons_j = recons_j,
                  #   grouping_j = grouping_j, pos_d = pos_d,
                  #   s = s, omega = omega, offset_j = offset_j,
                  #   long_var_j = long_var_j, ridge_j = ridge_j,
                  #   data_j = data_j, dim_j = dim_j
                  # )
                  # 
                  # fit_update_BASE <- bilinear_update(
                  #   Z_MI = Z_MI[[j]][Z_MI_grouping_j[[1]]],
                  #   Z_hier_mapping = Z_MI_hier_mapping[[j]][[1]],
                  #   prior_mi = prior_u, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                  #   diag_vi_pg_mean = diag_vi_pg_mean,
                  #   RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean,
                  #   offset = offset_j, it = 5,
                  #   return_chol = type_invert_j & (k == MAX_K),
                  #   method = 'sparse_direct'
                  # )
                  # 
                  # update_BASE <- vi_mi_mean[[j]]
                  # update_BASE[names(fit_update_BASE)] <- lapply(fit_update_BASE, `[[`, 'mean')
                  # val_BASE <- fff(mean_j = update_BASE,
                  #     recons_j = recons_j,
                  #     grouping_j = grouping_j, pos_d = pos_d,
                  #     s = s, omega = omega, offset_j = offset_j,
                  #     long_var_j = long_var_j, ridge_j = ridge_j,
                  #     data_j = data_j, dim_j = dim_j
                  # )
                  # 
                  # fit_update_PX <- bilinear_update(
                  #   Z_MI = Z_MI[[j]][Z_MI_grouping_j[[1]]],
                  #   Z_hier_mapping = Z_MI_hier_mapping[[j]][[1]],
                  #   prior_mi = prior_u, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                  #   diag_vi_pg_mean = diag_vi_pg_mean,
                  #   RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean,
                  #   offset = offset_j, it = Inf,
                  #   return_chol = type_invert_j & (k == MAX_K),
                  #   method = 'sparse_direct'
                  # )
                  # update_PX <- vi_mi_mean[[j]]
                  # update_PX[names(fit_update_PX)] <- lapply(fit_update_PX, `[[`, 'mean')
                  # val_PX <- fff(mean_j = update_PX,
                  #              recons_j = recons_j,
                  #              grouping_j = grouping_j, pos_d = pos_d,
                  #              s = s, omega = omega, offset_j = offset_j,
                  #              long_var_j = long_var_j, ridge_j = ridge_j,
                  #              data_j = data_j, dim_j = dim_j
                  # )
                  # 
                  # print('------------u---------------')
                  # print(c('prior' = val_1, 'standard' = val_BASE, 'px' = val_PX))
                  # if (val_PX < val_BASE & abs(val_PX - val_BASE) > 1e-6){
                  #   browser()
                  # }
                  # if (val_1 - val_BASE > 0){
                  #   browser()
                  #   print('UNUSUAL: NO HIER WORSE THAN PRIOR...')
                  #   warning('UNUSUAL: NO HIER WORSE THAN PRIOR...')
                  # }
                  # message(round(val_PX - val_BASE, 2))
                }else if (bi_method == 'joint_direct'){
                  update_u <- bilinear_update(
                    Z_MI = Z_MI[[j]][Z_MI_grouping_j[[1]]],
                    Z_hier_mapping = Z_MI_hier_mapping[[j]][[1]],
                    prior_mi = prior_u, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                    diag_vi_pg_mean = diag_vi_pg_mean,  it = it,
                    RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                    offset = 0, auxiliary_Z = XZ, auxiliary_prior = Tinv,
                    return_chol = type_invert_j & (k == MAX_K),
                    method = 'sparse_direct'
                  )
                  vi_beta_mean <- Matrix(update_u$aux[1:ncol(X)])
                  rownames(vi_beta_mean) <- colnames(X)
                  vi_alpha_mean <- Matrix(update_u$aux[-(1:ncol(X))])
                  rownames(vi_alpha_mean) <- colnames(Z)
                  update_u$aux <- NULL
                }else{
                  stop('invalid bi_method')
                }
              }else if (bi_method %in% c('joint', 'rowwise')){
                
                update_u <- bilinear_update(
                  Z_MI = Z_MI[[j]][[1]], Z_hier_mapping = Z_MI_hier_mapping[[j]][[1]],
                  prior_mi = prior_u, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                  diag_vi_pg_mean = diag_vi_pg_mean, 
                  RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                  offset = offset_j, it = it,
                  return_chol = type_invert_j & (k == MAX_K),
                  method = 'rowwise'
                )
              }else if (bi_method == 'sparse_direct'){
                update_u <- bilinear_update(
                  Z_MI = Z_MI[[j]][[1]],
                  prior_mi = prior_u, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                  diag_vi_pg_mean = diag_vi_pg_mean, 
                  RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                  offset = offset_j, it = it,
                  return_chol = type_invert_j & (k == MAX_K),
                  method = 'sparse_direct'
                )
              }else if (bi_method == 'joint_direct'){
                update_u <- bilinear_update(
                  Z_MI = Z_MI[[j]][[1]],
                  prior_mi = prior_u, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                  diag_vi_pg_mean = diag_vi_pg_mean, it = it,
                  RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                  offset = 0, auxiliary_Z = XZ, auxiliary_prior = Tinv,
                  return_chol = type_invert_j & (k == MAX_K),
                  method = 'sparse_direct'
                )
                vi_beta_mean <- Matrix(update_u$aux[1:ncol(X)])
                rownames(vi_beta_mean) <- colnames(X)
                vi_alpha_mean <- Matrix(update_u$aux[-(1:ncol(X))])
                rownames(vi_alpha_mean) <- colnames(Z)
              }else{
                stop('invalid bi_method')
              }
              
              if (simple_j){
                vi_mi_mean[[j]][[1]] <- update_u$mean
                vi_mi_var[[j]][[1]] <- update_u$inverse * VEM_scalar
                inner_vi_lndet[1] <- 2 * sum(update_u$det) * VEM_scalar
                temp_vi_lndet[1] <- 2 * sum(update_u$det) * VEM_scalar
                
                if (keep_mi_decomp & (k == MAX_K)){
                  if (type_invert_j){
                    vi_mi_marg_decomp[[j]][[1]] <- update_u$Linv
                  }else{
                    vi_mi_marg_decomp[[j]][[1]] <- sqrt(update_u$inverse)
                  }
                }
                
              }else{
                vi_mi_mean[[j]][Z_MI_grouping_j[[1]]] <- lapply(update_u, `[[`, 'mean')
                vi_mi_var[[j]][Z_MI_grouping_j[[1]]] <- lapply(update_u, FUN=function(i){i$inverse * VEM_scalar})
                inner_vi_lndet[Z_MI_grouping_j[[1]]] <- 2 * sapply(update_u, FUN=function(i){sum(i$det) * VEM_scalar})
                temp_vi_lndet[Z_MI_grouping_j[[1]]] <- 2 * sapply(update_u, FUN=function(i){sum(i$det) * VEM_scalar})
                
                if (keep_mi_decomp & (k == MAX_K)){
                  if (type_invert_j){
                    vi_mi_marg_decomp[[j]][Z_MI_grouping_j[[1]]] <- lapply(update_u, `[[`, "Linv")
                  }else{
                    vi_mi_marg_decomp[[j]][Z_MI_grouping_j[[1]]] <- lapply(update_u, FUN=function(i){sqrt(i$inverse)})
                  }
                }
                
              }
              
              vi_mi_sigma_outer_alpha <- get_bilinear_outer(
                vi_mi_mean, vi_mi_var, vi_mi_diag, mi_prior_type)
              ELBO_u <- g()
              
              rm(Rmean, RFSmean, Rvar, update_u)
              # Freezing u, update v
              if (simple_j){
                RFSmean <- Z_MI[[j]][[1]] %*% FS(vi_mi_mean[[j]][[1]], vi_mi_mean[[j]][[1]])
                Rmean <- Z_MI[[j]][[1]] %*% vi_mi_mean[[j]][[1]]
                Rvar <- Z_MI[[j]][[1]] %*% vi_mi_var[[j]][[1]]
              }else{
               
                Rmean <- Reduce('+', mapply(
                  Z_MI[[j]][Z_MI_grouping_j[[1]]],
                  vi_mi_mean[[j]][Z_MI_grouping_j[[1]]],
                  SIMPLIFY = FALSE,
                  FUN=function(i,j){
                    i %*% j
                  })) 
                
                # old_RFSmean <- Reduce('+', mapply(
                #   Z_MI[[j]][Z_MI_grouping_j[[1]]],
                #   vi_mi_mean[[j]][Z_MI_grouping_j[[1]]],
                #   SIMPLIFY = FALSE,
                #   FUN=function(i,j){
                #     i %*% FS(j,j)
                #   })) 

                # Can this be improved for speed?                
                RFSmean <- FS(Rmean, Rmean)

                Rvar <- Reduce('+', mapply(
                  Z_MI[[j]][Z_MI_grouping_j[[1]]],
                  vi_mi_var[[j]][Z_MI_grouping_j[[1]]],
                  SIMPLIFY = FALSE,
                  FUN=function(i,j){
                    i %*% j
                  }))
                
              }
              if (mi_prior_type %in% c('centered')){
                prior_v <- mi_j_vec_prior
              }else if (mi_prior_type %in% c('shared')){
                prior_v <- mi_j_vec_prior
              }else if (mi_prior_type %in% c('separate')){
                if (simple_j){
                  prior_v <- mi_j_vec_prior[[2]]
                }else{
                  prior_v <- mi_j_vec_prior[Z_MI_grouping_j[[2]]]
                }
              }
              
              if (!simple_j){
                
                if (bi_method %in% c('joint', 'rowwise')){
                  
                  update_v <- bilinear_update(
                    Z_MI = Z_MI[[j]][Z_MI_grouping_j[[2]]],
                    Z_hier_mapping = Z_MI_hier_mapping[[j]][[2]],
                    prior_mi = prior_v, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                    diag_vi_pg_mean = diag_vi_pg_mean, 
                    RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                    offset = offset_j, it = it,
                    return_chol = type_invert_j & (k == MAX_K),
                    method = 'sparse_direct'
                  )
                  
                  # long_var_j <- lapply(1:2, FUN=function(k){
                  #   Reduce('+', mapply(Z_MI[[j]][Z_MI_grouping_j[[k]]], vi_mi_var[[j]][Z_MI_grouping_j[[k]]], FUN=function(i,j){
                  #     i %*% j
                  #   }))
                  # })
                  # ridge_j <- mi_inv_mapping_alpha[[j]]
                  # dim_j <- mi_d_j[j]
                  # grouping_j <- Z_MI_grouping_j
                  # data_j <- Z_MI[[j]]
                  # recons_j <- lengths(Z_MI_attr[[j]]$levels)
                  # recons_j <- rep(names(recons_j), recons_j * dim_j)
                  # omega <- vi_pg_mean
                  # pos_d <- lapply(1:dim_j, FUN=function(d){dim_j * (d-1) + seq(dim_j)})
                  # 
                  # val_1 <- fff(mean_j = vi_mi_mean[[j]],
                  #              recons_j = recons_j,
                  #              grouping_j = grouping_j, pos_d = pos_d,
                  #              s = s, omega = omega, offset_j = offset_j,
                  #              long_var_j = long_var_j, ridge_j = ridge_j,
                  #              data_j = data_j, dim_j = dim_j
                  # )
                  # update_mean <- vi_mi_mean[[j]]
                  # update_mean[names(update_v)] <- lapply(update_v, `[[`, 'mean')
                  # val_2 <- fff(mean_j = update_mean,
                  #              recons_j = recons_j,
                  #              grouping_j = grouping_j, pos_d = pos_d,
                  #              s = s, omega = omega, offset_j = offset_j,
                  #              long_var_j = long_var_j, ridge_j = ridge_j,
                  #              data_j = data_j, dim_j = dim_j
                  # )
                  # 
                  # update_2 <- bilinear_update(
                  #   Z_MI = Z_MI[[j]][Z_MI_grouping_j[[2]]],
                  #   Z_hier_mapping = Z_MI_hier_mapping[[j]][[2]],
                  #   prior_mi = prior_v, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                  #   diag_vi_pg_mean = diag_vi_pg_mean,
                  #   RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean,
                  #   offset = offset_j, it = Inf,
                  #   return_chol = type_invert_j & (k == MAX_K),
                  #   method = 'sparse_direct'
                  # )
                  # update_mean_2 <- vi_mi_mean[[j]]
                  # update_mean_2[names(update_2)] <- lapply(update_2, `[[`, 'mean')
                  # val_3 <- fff(mean_j = update_mean_2,
                  #              recons_j = recons_j,
                  #              grouping_j = grouping_j, pos_d = pos_d,
                  #              s = s, omega = omega, offset_j = offset_j,
                  #              long_var_j = long_var_j, ridge_j = ridge_j,
                  #              data_j = data_j, dim_j = dim_j
                  # )
                  # print('------------v---------------')
                  # print(c('base' = val_1, 'standard' = val_2, 'px' = val_3))
                  # 
                  # if (val_2 > val_3 & abs(val_2 - val_3) > 1e-6){
                  #   browser()
                  # }
                  # if (val_2 - val_1 < 0){
                  #   browser()
                  # }else{
                  #   message(round(val_2 - val_1, 2))
                  # }
                  
                }else if (bi_method == 'joint_direct'){
                  update_v <- bilinear_update(
                    Z_MI = Z_MI[[j]][Z_MI_grouping_j[[2]]],
                    Z_hier_mapping = Z_MI_hier_mapping[[j]][[2]],
                    prior_mi = prior_v, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                    diag_vi_pg_mean = diag_vi_pg_mean,  it = it,
                    RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                    offset = 0, auxiliary_Z = XZ, auxiliary_prior = Tinv,
                    return_chol = type_invert_j & (k == MAX_K),
                    method = 'sparse_direct'
                  )
                  vi_beta_mean <- Matrix(update_v$aux[1:ncol(X)])
                  rownames(vi_beta_mean) <- colnames(X)
                  vi_alpha_mean <- Matrix(update_v$aux[-(1:ncol(X))])
                  rownames(vi_alpha_mean) <- colnames(Z)
                  update_v$aux <- NULL
                }else{stop('...')}
                
              }else if (bi_method %in% c('joint', 'rowwise')){
                
                update_v <- bilinear_update(
                  Z_MI = Z_MI[[j]][[2]], Z_hier_mapping = Z_MI_hier_mapping[[j]][[2]],
                  prior_mi = prior_v, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                  diag_vi_pg_mean = diag_vi_pg_mean, 
                  RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                  offset = offset_j, it = it,
                  return_chol = type_invert_j & (k == MAX_K),
                  method = 'rowwise'
                )
                
                # long_var_j <- lapply(1:2, FUN=function(k){
                #   Z_MI[[j]][[k]] %*% vi_mi_var[[j]][[k]]
                # })
                # ridge_j <- mi_inv_mapping_alpha[[j]]
                # dim_j <- mi_d_j[j]
                # grouping_j <- Z_MI_grouping_j
                # data_j <- Z_MI[[j]]
                # recons_j <- lengths(Z_MI_attr[[j]]$levels)
                # recons_j <- rep(names(recons_j), recons_j * dim_j)
                # omega <- vi_pg_mean
                # pos_d <- lapply(1:dim_j, FUN=function(d){dim_j * (d-1) + seq(dim_j)})
                # 
                # if (mi_prior_type == 'shared'){
                #   ridge_j <- list(ridge_j, ridge_j)
                # }
                # 
                # orig_mean <- vi_mi_mean[[j]]
                # val_1 <- fff(mean_j = vi_mi_mean[[j]],
                #              recons_j = recons_j,
                #              grouping_j = grouping_j, pos_d = pos_d,
                #              s = s, omega = omega, offset_j = offset_j,
                #              long_var_j = long_var_j, ridge_j = ridge_j,
                #              data_j = data_j, dim_j = dim_j
                # )
                # update_mean <- vi_mi_mean[[j]]
                # update_mean[[2]] <- update_v$mean
                # val_2 <- fff(mean_j = update_mean,
                #              recons_j = recons_j,
                #              grouping_j = grouping_j, pos_d = pos_d,
                #              s = s, omega = omega, offset_j = offset_j,
                #              long_var_j = long_var_j, ridge_j = ridge_j,
                #              data_j = data_j, dim_j = dim_j
                # )
                # 
                # print('------------v---------------')
                # print(c('base' = val_1, 'standard' = val_2))
                # 
                # if (val_2 - val_1 < 0){
                #   browser()
                # }else{
                #   message(round(val_2 - val_1, 2))
                # }
                
              }else if (bi_method == 'sparse_direct'){
                update_v <- bilinear_update(
                  Z_MI = Z_MI[[j]][[2]],
                  prior_mi = prior_v, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                  diag_vi_pg_mean = diag_vi_pg_mean, 
                  RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                  offset = offset_j, it = it,
                  return_chol = type_invert_j & (k == MAX_K),
                  method = 'sparse_direct'
                )
              }else{
                stop('...')
                update_v <- bilinear_update(
                  Z_MI = Z_MI[[j]][[2]],
                  prior_mi = prior_v, d_mi = mi_d_j[j], s = s, vi_pg_mean = vi_pg_mean,
                  diag_vi_pg_mean = diag_vi_pg_mean, it = it,
                  RFSmean = RFSmean, Rvar = Rvar, Rmean = Rmean, 
                  offset = 0, auxiliary_Z = XZ, auxiliary_prior = Tinv, 
                  return_chol = type_invert_j & (k == MAX_K),
                  method = 'sparse_direct'
                )
                vi_beta_mean <- Matrix(update_v$aux[1:ncol(X)])
                rownames(vi_beta_mean) <- colnames(X)
                vi_alpha_mean <- Matrix(update_v$aux[-(1:ncol(X))])
                rownames(vi_alpha_mean) <- colnames(Z)
              }
              
              if (simple_j){
                
                vi_mi_mean[[j]][[2]] <- update_v$mean
                vi_mi_var[[j]][[2]] <- update_v$inverse * VEM_scalar
                inner_vi_lndet[2] <- 2 * sum(update_v$det) * VEM_scalar
                temp_vi_lndet[2] <- 2 * sum(update_v$det) * VEM_scalar
                if (keep_mi_decomp & (k == MAX_K)){
                  if (type_invert_j){
                    vi_mi_marg_decomp[[j]][[2]] <- update_v$Linv
                  }else{
                    vi_mi_marg_decomp[[j]][[2]] <- sqrt(update_v$inverse)
                  }
                }
                
              }else{
                
                vi_mi_mean[[j]][Z_MI_grouping_j[[2]]] <- lapply(update_v, `[[`, 'mean')
                vi_mi_var[[j]][Z_MI_grouping_j[[2]]] <- lapply(update_v, FUN=function(i){i$inverse * VEM_scalar})
                inner_vi_lndet[Z_MI_grouping_j[[2]]] <- 2 * sapply(update_v, FUN=function(i){sum(i$det) * VEM_scalar})
                temp_vi_lndet[Z_MI_grouping_j[[2]]] <- 2 * sapply(update_v, FUN=function(i){sum(i$det) * VEM_scalar})
                
                if (keep_mi_decomp & (k == MAX_K)){
                  if (type_invert_j){
                    vi_mi_marg_decomp[[j]][Z_MI_grouping_j[[2]]] <- lapply(update_v, `[[`, "Linv")
                  }else{
                    vi_mi_marg_decomp[[j]][Z_MI_grouping_j[[2]]] <- lapply(update_v, FUN=function(i){sqrt(i$inverse)})
                  }
                }
              }
              
              vi_mi_sigma_outer_alpha <- get_bilinear_outer(
                vi_mi_mean, vi_mi_var, vi_mi_diag, mi_prior_type)
              
              ELBO_v <- g()
              debug_MI <- rbind(ELBO_start, ELBO_u, ELBO_v)[,1,drop=F]
              print(
                debug_MI
              )
              
              if (it %in% c(1, MI_STAB, MI_STAB + 1)){
                
              }else{
                if (any(diff(debug_MI[,1]) < -sqrt(.Machine$double.eps))){browser()}
              }
              if (MAX_K > 0){

                print(mapply(vi_mi_mean[[j]], init_mi, FUN=function(i,j){
                  max(abs(i-j))
                  # sqrt(sum(abs(i-j)^2))
                }))
                init_mi <- vi_mi_mean[[j]]
              }
              
            }
           
          }
          
          vi_mi_lndet[[j]] <- inner_vi_lndet
          # Update running linear predictor
          running_billinear_lp <- offset_j +
            get_bilinear_mean(Z_MI[j], vi_mi_mean[j], Z_MI_grouping[j])

          if (simple_j & (bi_method == 'sparse_direct')){
            variance_by_alpha_jg <- calculate_expected_outer_alpha(
              L = vi_alpha_decomp, 
              alpha_mu = as.vector(vi_alpha_mean),
              re_position_list = outer_alpha_RE_positions)
            
            vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
          }
          
          vi_mi_sigma_outer_alpha <- get_bilinear_outer(
            vi_mi_mean, vi_mi_var, vi_mi_diag, mi_prior_type)

        }
        # Update the Sigma_j for bilinear terms
        if (any_mi & update_mi){
          
          if (freeze_mi_var){
            # Do Nothing
          }else{
            
            vi_mi_sigma_outer_alpha <- get_bilinear_outer(
              vi_mi_mean, vi_mi_var, vi_mi_diag, mi_prior_type)
            
            if (!do_huangwand_mi){

              vi_mi_sigma_alpha <- mapply(
                vi_mi_sigma_outer_alpha, mi_prior_sigma_alpha_phi,
                SIMPLIFY = FALSE, FUN = function(i, j) {
                  if (mi_prior_type %in% c('centered', 'shared')){
                    i * vi_sigmasq_a/vi_sigmasq_b + j
                  }else{
                    mapply(i, j, SIMPLIFY = FALSE, FUN=function(a,b){
                      a * vi_sigmasq_a/vi_sigmasq_b + b
                    })
                  }
                })
              
            }else{
              

              if (it < (MI_STAB + 5)){
                INNER_IT_MI <- 1
              }else{
                INNER_IT_MI <- INNER_IT
              }
              
              for (inner_it in 1:INNER_IT_MI){
                
                vi_mi_sigma_alpha <- mapply(
                  vi_mi_sigma_outer_alpha, vi_mi_a_a_jp,
                  vi_mi_a_b_jp, vi_mi_a_nu_jp, Z_MI_grouping, 1:N_MI, SIMPLIFY = FALSE,
                  FUN = function(i, tilde.a, tilde.b, nu, group.j, pos.j) {
                    if (mi_prior_type %in% c('centered', 'shared')){
                      i * vi_sigmasq_a/vi_sigmasq_b +
                        Diagonal(x = tilde.a/tilde.b) * 2 * nu
                    }else{
                      prop_out <- mapply(i, tilde.a, tilde.b, SIMPLIFY = FALSE,
                             FUN=function(i_l, tilde.a_l, tilde.b_l){
                               i_l * vi_sigmasq_a/vi_sigmasq_b +
                                 Diagonal(x=tilde.a_l/tilde.b_l) * 2 * nu
                             })
                      if (partial_fix){
                        for (k in 1:2){
                          k_name <- group.j[[k]][1]
                          prop_out[[k_name]] <- vi_mi_sigma_alpha[[pos.j]][[k_name]]
                        }
                      }
                      return(prop_out)
                    }
                  })
                #Update a_{j,p}
                mi_diag_Einv_sigma <- mapply(
                  vi_mi_sigma_alpha, vi_mi_sigma_alpha_nu, mi_d_j,
                  SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
                    if (mi_prior_type %in% c('centered', 'shared')){
                      inv_phi <- solve(phi)
                      sigma.inv <- nu * inv_phi
                      return(diag(sigma.inv))
                    }else{
                      mapply(phi, nu, SIMPLIFY = FALSE, FUN=function(phi_l, nu_l){
                        inv_phi_l <- solve(phi_l)
                        sigma.inv <- nu_l * inv_phi_l
                        return(diag(sigma.inv))
                        
                      })
                    }
                  })
                
                vi_mi_a_b_jp <- mapply(
                  vi_mi_a_nu_jp, vi_mi_a_APRIOR_jp, mi_diag_Einv_sigma,
                  SIMPLIFY = FALSE,
                  FUN=function(nu, APRIOR, diag_j){
                    if (mi_prior_type %in% c('centered', 'shared')){
                      1/APRIOR^2 + nu * diag_j
                    }else{
                      mapply(APRIOR, diag_j, SIMPLIFY = FALSE, FUN=function(APRIOR_l, diag_j_l){
                        1/APRIOR_l^2 + nu * diag_j_l
                      })
                    }
                  })
                
              }
              
            }
            
          }
        }
      }else{stop('lll must be mi or re')}
      
      if (debug_ELBO & it != 1) {
        
        if (any_RE){
          variance_by_alpha_jg <- calculate_expected_outer_alpha(
            L = vi_alpha_decomp, 
            alpha_mu = as.vector(vi_alpha_mean),
            re_position_list = outer_alpha_RE_positions)
          
          vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
        }
        
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
                                       vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
                                       # Multiplicative Interaction
                                       do_huangwand_mi = do_huangwand_mi,
                                       any_RE = any_RE, mi_prior_type = mi_prior_type,
                                       any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
                                       vi_mi_diag = vi_mi_diag, 
                                       vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
                                       vi_mi_sigma_alpha = vi_mi_sigma_alpha,
                                       vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
                                       vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
                                       vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
                                       vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
                                       mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
                                       mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
                                       mi_iw_prior_constant = mi_iw_prior_constant,
                                       mi_d_j = mi_d_j, mi_g_j = mi_g_j
        )
        
      }
      
    }

    
    if (family == 'linear'){
      
      if (any_mi){browser()}
      
      adjust_var <- 1/sqrt(vi_sigmasq_a/vi_sigmasq_b)
      
      vi_beta_decomp <- vi_beta_decomp * adjust_var
      vi_alpha_decomp <- vi_alpha_decomp * adjust_var
      vi_joint_decomp <- vi_joint_decomp * adjust_var
      
      if (factorization_method == 'weak'){
        vi_joint_L_nonpermute <- vi_joint_L_nonpermute * adjust_var
      }else{
        vi_beta_L_nonpermute <- vi_beta_L_nonpermute * adjust_var
        vi_alpha_L_nonpermute <- vi_alpha_L_nonpermute * adjust_var
      }

      ln_sigmasq <- log(vi_sigmasq_b) - log(vi_sigmasq_a)
      log_det_joint_var <- log_det_joint_var + ncol(vi_joint_decomp) * ln_sigmasq
      log_det_beta_var <- log_det_beta_var + ncol(vi_beta_decomp) * ln_sigmasq
      log_det_alpha_var <- log_det_alpha_var + ncol(vi_alpha_decomp) * ln_sigmasq
    }

    if (do_timing) {
      toc(quiet = quiet_time, log = T)
      tic("Update Sigma")
    }
    


    if (debug_ELBO & it > 1) {
      
      temp_ELBO <- calculate_ELBO(family = family,
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
          vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
          # Multiplicative Interaction
          do_huangwand_mi = do_huangwand_mi,
          any_RE = any_RE, mi_prior_type = mi_prior_type,
          any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
          vi_mi_diag = vi_mi_diag, 
          vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
          vi_mi_sigma_alpha = vi_mi_sigma_alpha,
          vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
          vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
          vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
          vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
          mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
          mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
          mi_iw_prior_constant = mi_iw_prior_constant,
          mi_d_j = mi_d_j, mi_g_j = mi_g_j
        )
      
      if (any_mi & it > 1){
        if (temp_ELBO$ELBO - debug_ELBO.2$ELBO < -sqrt(.Machine$double.eps) ){
          browser()
        }
      }
      
    }
    
    if (do_timing) {
      toc(quiet = quiet_time, log = T)
      tic("Update Aux")
    }
    
    
    # Update the auxilary parameters
    
    if (family == "negbin") {
      
      if (any_mi){browser()}

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
      
      if (any_mi){
        bilinear_mean <- get_bilinear_mean(Z_MI, vi_mi_mean, Z_MI_grouping)
        bilinear_var <- get_bilinear_var(Z_MI, vi_mi_mean, vi_mi_var, Z_MI_grouping)
      }else{
        bilinear_mean <- 0
        bilinear_var <- 0
      }
      
      if (factorization_method == 'weak'){
        joint_quad <- bilinear_var + cpp_zVz(Z = joint.XZ, V = as(vi_joint_decomp, "generalMatrix"))
        vi_lp <- (s - as.vector(bilinear_mean + X %*% vi_beta_mean + Z %*% vi_alpha_mean))^2 + joint_quad
      } else{
        # beta_quad <- rowSums((X %*% t(vi_beta_decomp))^2)
        beta_quad <- cpp_dense_zVz(X, as.matrix(vi_beta_decomp))
        alpha_quad <- rowSums((Z %*% t(vi_alpha_decomp))^2)
        vi_lp <- (s - as.vector(bilinear_mean + X %*% vi_beta_mean + Z %*% vi_alpha_mean))^2 + 
          beta_quad + alpha_quad + bilinear_var
      }
      
      if (any_mi){browser()}
      
      vi_kernel <- expect_alpha_prior_kernel(vi_sigma_alpha = vi_sigma_alpha, 
          vi_sigma_alpha_nu = vi_sigma_alpha_nu, d_j = d_j,
          vi_sigma_outer_alpha = vi_sigma_outer_alpha)
      vi_sigmasq_b <- (sum(vi_lp) + vi_kernel)/2 + vi_sigmasq_prior_b
      
    }

    if (do_timing) {
      toc(quiet = quiet_time, log = T)
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
        vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
        # Multiplicative Interaction
        do_huangwand_mi = do_huangwand_mi,
        any_RE = any_RE, mi_prior_type = mi_prior_type,
        any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
        vi_mi_diag = vi_mi_diag, 
        vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
        vi_mi_sigma_alpha = vi_mi_sigma_alpha,
        vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
        vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
        vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
        vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
        mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
        mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
        mi_iw_prior_constant = mi_iw_prior_constant,
        mi_d_j = mi_d_j, mi_g_j = mi_g_j
      )
      
      if (any_mi & it > 1){
        if (debug_ELBO.3$ELBO - temp_ELBO$ELBO < -sqrt(.Machine$double.eps)){
          browser()
        }
      }
      
    }

    accept.PX <- NA
    
    if (parameter_expansion == "none" | !any_Mprime) {
      
      accept.PX <- TRUE
      
    } else {
      
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
           vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
           # Multiplicative Interaction
           do_huangwand_mi = do_huangwand_mi,
           any_RE = any_RE, mi_prior_type = mi_prior_type,
           any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
           vi_mi_diag = vi_mi_diag, 
           vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
           vi_mi_sigma_alpha = vi_mi_sigma_alpha,
           vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
           vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
           vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
           vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
           mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
           mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
           mi_iw_prior_constant = mi_iw_prior_constant,
           mi_d_j = mi_d_j, mi_g_j = mi_g_j
        )
      }
      
      if (do_timing) {
        tic("Update PX")
      }

      if (any_mi & do_PX_MI & update_mi){
        
        change_px_mi_rotation <- 0
        change_px_mi_add <- 0
        
        for (px_mi_order in px_mi_type){

          if (px_mi_order == 'rotate'){
            px_mi_moments <- mapply(
              vi_mi_sigma_alpha, vi_mi_sigma_alpha_nu, mi_d_j,
              SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
                if (mi_prior_type %in% c('centered', 'shared')){
                  inv_phi <- solve(phi)
                  sigma.inv <- nu * inv_phi
                  return(sigma.inv)
                }else if (mi_prior_type %in% c('separate')){
                  mapply(nu, phi, SIMPLIFY = FALSE, FUN=function(nu_l, phi_l){
                    inv_phi_l <- solve(phi_l)
                    sigma.inv_l <- nu_l * inv_phi_l
                    return(sigma.inv_l)
                  })
                }else{stop('...')}
              })
            
            # Rotation Expansion for PX
            # Get E[u_g u_g^T] and E[v_g v_g^T]
            outer_px <- get_bilinear_outer(
              vi_mi_mean = vi_mi_mean,
              vi_mi_var = vi_mi_var,
              vi_mi_diag = vi_mi_diag,
              mi_prior_type = 'blank',
              reduce = FALSE)

            px_sigma <- !freeze_mi_var
            if (partial_fix){
              px_sigma <- FALSE
            }
            
            px_mi_rotation <- update_px_rotation(
              vi_mi_SSQ = outer_px,
              vi_mi_moments = px_mi_moments,
              vi_mi_size = outer_px_size,
              vi_mi_dim = mi_d_j,
              vi_mi_a_a_jp = vi_mi_a_a_jp,
              vi_mi_a_b_jp = vi_mi_a_b_jp,
              vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
              vi_mi_a_nu_jp = vi_mi_a_nu_jp,
              Z_MI_grouping = Z_MI_grouping,
              mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
              mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu,
              do_huangwand_mi = do_huangwand_mi,
              mi_prior_type = mi_prior_type,
              px_sigma = px_sigma, VEM_scalar = VEM_scalar,
              partial_fix = partial_fix
            )
            
            rm(outer_px, px_mi_moments)
            
            easy_message(px_mi_rotation$R)
            # easy_message(lapply(px_mi_rotation, `[[`, 'rho'))
            if (px_sigma){
              vi_mi_sigma_alpha <- mapply(
                vi_mi_sigma_alpha, px_mi_rotation$R, Z_MI_grouping, SIMPLIFY = FALSE,
                FUN=function(Phi, R, group_hier){
                  if (mi_prior_type %in% c('centered', 'shared')){
                    return(R %*% Phi %*% t(R))
                  }else if (mi_prior_type %in% 'separate'){
                    
                    simple_j <- all(lengths(group_hier) == 1)
                    if (simple_j){
                      new_Phi <- Phi
                      new_Phi[[2]] <- R %*% Phi[[2]] %*% t(R)
                    }else{
                      new_Phi <- Phi
                      if (partial_fix){
                        for (v in group_hier[[2]][-1]){
                          new_Phi[[v]] <- R %*% new_Phi[[v]] %*% t(R)
                        }
                      }else{
                        new_Phi[group_hier[[2]]] <- lapply(
                          Phi[group_hier[[2]]], FUN=function(Phi_l){
                            R %*% Phi_l %*% t(R)
                          })
                      }
                    }
                    return(new_Phi)
                  }else{stop('...')}
                })
              
              if (do_huangwand_mi){
                if (mi_prior_type %in% c('centered', 'shared')){
                  vi_mi_a_b_jp <- px_mi_rotation$rho_hw
                  names(vi_mi_a_b_jp) <- names(vi_mi_a_a_jp)
                }else if (mi_prior_type %in% c('separate')){
                  vi_mi_a_b_jp <- mapply(vi_mi_a_b_jp, 
                    px_mi_rotation$rho_hw, Z_MI_grouping, 
                    SIMPLIFY = FALSE, FUN=function(a,b, group_hier){
                      simple_j <- all(lengths(group_hier) == 1)
                      if (is.null(names(b)) & simple_j){
                        a[[2]] <- b[[2]]
                      }else{
                        a[group_hier[[2]]] <- b[group_hier[[2]]]
                      }
                      return(a)
                  })
                }else{
                  stop('...')
                }
              }
              
            }
            
            if (debug_px){
              change_px_mi_rotation <- change_px_mi_rotation + sum(px_mi_rotation$diff)
            }
            
            for (j in 1:N_MI){
              rho_j <- t(px_mi_rotation$R[[j]])
              Z_MI_grouping_j <- Z_MI_grouping[[j]]
              simple_j <- all(lengths(Z_MI_grouping_j) == 1)
              if (mi_d_j[j] == 1){
                rho_j <- as.vector(rho_j)
                if (simple_j){
                  vi_mi_mean[[j]][[1]] <- vi_mi_mean[[j]][[1]] * 1/rho_j
                  vi_mi_var[[j]][[1]] <- vi_mi_var[[j]][[1]] * 1/rho_j^2 * VEM_scalar
                  vi_mi_lndet[[j]][[1]] <- vi_mi_lndet[[j]][[1]] + -outer_px_size[[j]][1] * 2 * log(abs(rho_j)) * VEM_scalar
                  
                  vi_mi_mean[[j]][[2]] <- vi_mi_mean[[j]][[2]] * rho_j
                  vi_mi_var[[j]][[2]] <- vi_mi_var[[j]][[2]] * rho_j^2 * VEM_scalar
                  vi_mi_lndet[[j]][[2]] <- vi_mi_lndet[[j]][[2]] + outer_px_size[[j]][2] * 2 * log(abs(rho_j)) * VEM_scalar
                  
                  if (keep_mi_decomp){
                    vi_mi_marg_decomp[[j]][[1]] <- vi_mi_marg_decomp[[j]][[1]] * 1/rho_j * VEM_scalar
                    vi_mi_marg_decomp[[j]][[2]] <- vi_mi_marg_decomp[[j]][[2]] * rho_j * VEM_scalar
                  }
                }else{
                  for (u_name in Z_MI_grouping_j[[1]]){
                    vi_mi_mean[[j]][[u_name]] <- vi_mi_mean[[j]][[u_name]] * 1/rho_j
                    vi_mi_var[[j]][[u_name]] <- vi_mi_var[[j]][[u_name]] * 1/rho_j^2 * VEM_scalar
                    vi_mi_lndet[[j]][[u_name]] <- vi_mi_lndet[[j]][[u_name]] + -outer_px_size[[j]][u_name] * 2 * log(abs(rho_j)) * VEM_scalar
                    if (keep_mi_decomp){
                      vi_mi_marg_decomp[[j]][[u_name]] <- vi_mi_marg_decomp[[j]][[u_name]] * 1/rho_j * VEM_scalar
                    }
                  }
                  for (v_name in Z_MI_grouping_j[[2]]){
                    vi_mi_mean[[j]][[v_name]] <- vi_mi_mean[[j]][[v_name]] * rho_j
                    vi_mi_var[[j]][[v_name]] <- vi_mi_var[[j]][[v_name]] * rho_j^2 * VEM_scalar
                    vi_mi_lndet[[j]][[v_name]] <- vi_mi_lndet[[j]][[v_name]] + outer_px_size[[j]][v_name] * 2 * log(abs(rho_j)) * VEM_scalar
                    if (keep_mi_decomp){
                      vi_mi_marg_decomp[[j]][[v_name]] <- vi_mi_marg_decomp[[j]][[v_name]] * rho_j * VEM_scalar
                    }
                    
                  }
                }
                
              }else{
                
                I_j <- Diagonal(n = ncol(rho_j))
                inv_rho_j <- solve(rho_j)
                det_rho_j <- determinant(rho_j)
                # if (sign(det_rho_j$sign) != 1){
                #   stop('Negative Determinant')
                # }
                lndet_rho_j <- as.numeric(det_rho_j$modulus)
                
                if (simple_j){
                  vi_mi_mean[[j]][[1]] <- vi_mi_mean[[j]][[1]] %*% t(inv_rho_j)
                  vi_mi_var[[j]][[1]] <- vi_mi_var[[j]][[1]] %*%
                    t(kronecker(inv_rho_j, inv_rho_j)) * VEM_scalar
                  vi_mi_lndet[[j]][[1]] <- vi_mi_lndet[[j]][[1]] + -outer_px_size[[j]][1] * 2 * lndet_rho_j * VEM_scalar
                  
                  vi_mi_mean[[j]][[2]] <- vi_mi_mean[[j]][[2]] %*% t(t(rho_j))
                  vi_mi_var[[j]][[2]] <- vi_mi_var[[j]][[2]] %*%
                    t(kronecker(t(rho_j), t(rho_j))) * VEM_scalar
                  vi_mi_lndet[[j]][[2]] <- vi_mi_lndet[[j]][[2]] + outer_px_size[[j]][2] * 2 * lndet_rho_j * VEM_scalar
                  
                  if (keep_mi_decomp){
                    vi_mi_marg_decomp[[j]][[1]] <- as.matrix(vi_mi_marg_decomp[[j]][[1]] %*%
                                                          t(kronecker(inv_rho_j, I_j))) * VEM_scalar
                    vi_mi_marg_decomp[[j]][[2]] <- as.matrix(vi_mi_marg_decomp[[j]][[2]] %*%
                                                          t(kronecker(t(rho_j), I_j))) * VEM_scalar
                  }
                }else{
                  for (u_name in Z_MI_grouping_j[[1]]){
                    vi_mi_mean[[j]][[u_name]] <- vi_mi_mean[[j]][[u_name]] %*% t(inv_rho_j)
                    vi_mi_var[[j]][[u_name]] <- vi_mi_var[[j]][[u_name]] %*%
                      t(kronecker(inv_rho_j, inv_rho_j)) * VEM_scalar
                    vi_mi_lndet[[j]][[u_name]] <- vi_mi_lndet[[j]][[u_name]] + -outer_px_size[[j]][u_name] * 2 * lndet_rho_j * VEM_scalar
                    if (keep_mi_decomp){
                      vi_mi_marg_decomp[[j]][[u_name]] <- as.matrix(vi_mi_marg_decomp[[j]][[u_name]] %*%
                                                            t(kronecker(inv_rho_j, I_j))) * VEM_scalar
                    }
                  }
                  for (v_name in Z_MI_grouping_j[[2]]){
                    vi_mi_mean[[j]][[v_name]] <- vi_mi_mean[[j]][[v_name]] %*% t(t(rho_j))
                    vi_mi_var[[j]][[v_name]] <- vi_mi_var[[j]][[v_name]] %*%
                      t(kronecker(t(rho_j), t(rho_j))) * VEM_scalar
                    vi_mi_lndet[[j]][[v_name]] <- vi_mi_lndet[[j]][[v_name]] + outer_px_size[[j]][v_name] * 2 * lndet_rho_j * VEM_scalar
                    if (keep_mi_decomp){
                      vi_mi_marg_decomp[[j]][[v_name]] <- as.matrix(vi_mi_marg_decomp[[j]][[v_name]] %*%
                                                            t(kronecker(t(rho_j), I_j))) * VEM_scalar
                    }
                    
                  }
                }
              }
              rm(simple_j)
            }
          }
          
          # ### Additive Expansion for PX
          if (px_mi_order == 'add'){
            
            px_mi_moments <- mapply(
              vi_mi_sigma_alpha, vi_mi_sigma_alpha_nu, mi_d_j,
              SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
                if (mi_prior_type %in% c('centered', 'shared')){
                  inv_phi <- solve(phi)
                  sigma.inv <- nu * inv_phi
                  return(sigma.inv)
                }else{
                  mapply(nu, phi, SIMPLIFY = FALSE, FUN=function(nu_l, phi_l){
                    inv_phi_l <- solve(phi_l)
                    sigma.inv_l <- nu_l * inv_phi_l
                    return(sigma.inv_l)
                  })
                }
              })
            
            flat_mean <- get_bilinear_mean(Z_MI, vi_mi_mean, Z_MI_grouping, reduce = FALSE)
            flat_var <- get_bilinear_var(Z_MI, vi_mi_mean, vi_mi_var, Z_MI_grouping, reduce = FALSE)
            
            moments_RE <- mapply(vi_sigma_alpha, vi_sigma_alpha_nu, d_j,
             SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
               inv_phi <- solve(phi)
               sigma.inv <- nu * inv_phi
               return(sigma.inv)
             })
            
            px_mi_add <- mapply(outer_px_add_map, outer_px_size,
            mi_d_j, px_mi_moments, Z_MI, vi_mi_mean, flat_mean, flat_var, Z_MI_grouping, SIMPLIFY = FALSE,
            FUN=function(pos_j, g_j, d_j, vc_j, Z_j, mean_j, flat_mean_j, flat_var_j, grouping_j){
              
              if (mi_prior_type %in% c('centered')){
                vc_j <- list(Diagonal(n=d_j), vc_j)
              }else if (mi_prior_type %in% c('shared')){
                vc_j <- lapply(1:2, FUN=function(i){vc_j})
              }else{
                vc_j <- vc_j
              }
              simple_j <- all(lengths(grouping_j) == 1)
              out <- lapply(1:2, FUN=function(k){
                if (simple_j){

                  # Terms from "l"
                  pos_l <- pos_j[[-k]]
                  if (is.null(pos_l)){
                    return(NULL)
                  }
                  mean_l <- mean_j[[-k]]
                  alpha_mean <- vi_alpha_mean[pos_l$position]
                  alpha_ESigma.inv <- moments_RE[[pos_l$index]]
                  if (ncol(alpha_ESigma.inv) != 1){
                    stop('additive PX not set up yet..')
                  }
                  # Terms from "k"
                  g_k <- g_j[k]
                  vc_k <- vc_j[[k]]
                  mean_k <- mean_j[[k]]
                  
                  flat_mean_k <- flat_mean_j[[k]]
                  flat_var_l <- flat_var_j[[-k]]
                  
                  
                  aug_l <- cbind(-1,mean_l)
                  denom <- g_k * bdiag(0, vc_k) +
                    as.numeric(alpha_ESigma.inv) * crossprod(aug_l) +
                    bdiag(0, matrix(
                      colSums(Diagonal(x=vi_pg_mean) %*% flat_var_l),
                      ncol = d_j))
                  numer <- rbind(0, vc_k %*% colSums(mean_k)) +
                    -as.numeric(alpha_ESigma.inv) * as.vector(t(alpha_mean) %*% aug_l)
                  numer <- numer + 
                    c(0, rowSums(matrix(
                      colSums(flat_var_l * 
                                kronecker(Diagonal(x=vi_pg_mean) %*% flat_mean_k,
                                          matrix(1, ncol = d_j))),
                      nrow = d_j
                    )))
                  rho_add <- as.vector(solve(denom, numer))
                  name_receive_l <- NULL
                }else{
                  group_l <- grouping_j[[-k]]
                  group_k <- grouping_j[[k]]
                  # Terms from "l"
                  pos_l <- pos_j[group_l]
                  if (all(sapply(pos_l, is.null))){
                    return(NULL)
                  }
                  
                  # Get the receiving "l" terms that are not null
                  reciving_l <- pos_l[!sapply(pos_l, is.null)]
                  if (length(reciving_l) == 1){
                    name_receive_l <- names(reciving_l)
                    reciving_l <- reciving_l[[1]]
                  }else{
                    # If multiple, take the principal one
                    name_receive_l <- group_l[1]
                    reciving_l <- reciving_l[[group_l[1]]]
                    if (is.null(reciving_l)){
                      browser()
                    }
                  }
                  mean_l <- mean_j[[name_receive_l]]
                  alpha_mean <- vi_alpha_mean[reciving_l$position]
                  alpha_ESigma.inv <- moments_RE[[reciving_l$index]]
                  if (ncol(alpha_ESigma.inv) != 1){
                    stop('Not set up for random slopes on alpha_j')
                  }
                  # Terms from "k"
                  g_k <- g_j[grouping_j[[k]]]
                  vc_k <- vc_j[grouping_j[[k]]]
                  mean_k <- mean_j[grouping_j[[k]]]
                  
                  flat_mean_k <- flat_mean_j[grouping_j[[k]]]
                  full_flat_mean_k <- flat_mean_k
                  flat_mean_k <- Reduce('+', flat_mean_k)
                  flat_var_l <- flat_var_j[grouping_j[[-k]]]
                  flat_var_l <- Reduce('+', flat_var_l)
                  
                  aug_l <- cbind(-1,mean_l)
                  
                  K <- kronecker(t(matrix(1, length(group_k))), Diagonal(d_j))
                  K_aug <- bdiag(1, K)
                  
                  denom_1 <- bdiag(0,
                        bdiag(mapply(g_k, vc_k, SIMPLIFY = FALSE, FUN=function(i,j){i * j}))
                  )
                  denom_2 <- as.numeric(alpha_ESigma.inv) * (t(K_aug) %*% crossprod(aug_l) %*% K_aug)
                  denom_3 <- matrix(
                    colSums(Diagonal(x=vi_pg_mean) %*% flat_var_l),
                    ncol = d_j)
                  denom_3 <- bdiag(0, t(K) %*% denom_3 %*% K)
                  
                  denom <- denom_1 + denom_2 + denom_3
                  # denom <- g_k * bdiag(0, vc_k) +
                  #   as.numeric(alpha_ESigma.inv) * crossprod(aug_l) +
                  #   bdiag(0, matrix(
                  #     colSums(Diagonal(x=vi_pg_mean) %*% flat_var_l),
                  #     ncol = d_j))
                  
                  numer_1 <- rbind(0, do.call('rbind', mapply(vc_k, mean_k, SIMPLIFY = FALSE, FUN=function(i,j){
                    i %*% colSums(j)
                  })))
                  numer_2 <- -as.numeric(alpha_ESigma.inv) * (t(K_aug) %*% t(t(alpha_mean) %*% aug_l))
                  numer_3 <- c(0, 
                               as.vector(t(K) %*% matrix(
                                 colSums(flat_var_l * 
                                           kronecker(Diagonal(x=vi_pg_mean) %*% flat_mean_k,
                                                     matrix(1, ncol = d_j))),
                                 nrow = d_j
                               )))
                  numer <- numer_1 + numer_2 + numer_3
                  # numer <- rbind(0, vc_k %*% colSums(mean_k)) +
                  #   -as.numeric(alpha_ESigma.inv) * as.vector(t(alpha_mean) %*% aug_l)
                  # numer <- numer + 
                  #   c(0, rowSums(matrix(
                  #     colSums(flat_var_l * 
                  #               kronecker(Diagonal(x=vi_pg_mean) %*% flat_mean_k,
                  #                         matrix(1, ncol = d_j))),
                  #     nrow = d_j
                  #   )))
                  rho_add <- as.vector(solve(denom, numer))
                  attr(rho_add, 'receive') <- name_receive_l
                }
                warning('OPTIM FOR RHO ADD')
                opt_rho_add <- optim(par = rep(0, length(rho_add)),
                    f = f_add, mean_k = mean_k, vc_k = vc_k,
                    alpha_mean = alpha_mean, mean_l = mean_l,
                    alpha_ESigma.inv = alpha_ESigma.inv,
                    pg = vi_pg_mean, long_mean_k = full_flat_mean_k,
                    long_var_l = flat_var_l, simple = simple_j,
                    method = 'L-BFGS',
                    control = list(fnscale = -1, maxit = 15))
                rho_add[] <- opt_rho_add$par
                
                f_diff <- 
                  f_add(
                    rho = rho_add, mean_k = mean_k, vc_k = vc_k,
                    alpha_mean = alpha_mean, mean_l = mean_l,
                    alpha_ESigma.inv = alpha_ESigma.inv, 
                    pg = vi_pg_mean, long_mean_k = full_flat_mean_k,
                    long_var_l = flat_var_l, simple = simple_j) -
                  f_add(
                    rho = rep(0, length(rho_add)), mean_k = mean_k, vc_k = vc_k,
                    alpha_mean = alpha_mean, mean_l = mean_l,
                    alpha_ESigma.inv = alpha_ESigma.inv, 
                    pg = vi_pg_mean, long_mean_k = full_flat_mean_k,
                    long_var_l = flat_var_l, simple = simple_j)  
                
                if (f_diff < -sqrt(.Machine$double.eps)){browser()}
                
                # optimize(f = f_add, mean_k = mean_k, vc_k = vc_k,
                #       alpha_mean = alpha_mean, mean_l = mean_l, 
                #       alpha_ESigma.inv = alpha_ESigma.inv, 
                #       pg = vi_pg_mean, long_mean_k = flat_mean_k,
                #       long_var_l = flat_var_l,
                #       interval = c(-5,5), maximum = TRUE)
                
                if (debug_px){
                  return(list(diff = f_diff, rho = rho_add))
                }else{
                  return(rho_add)
                }
              })
              return(out)
            })
            message(px_mi_add)
            
            if (debug_px){
              change_px_mi_add <- change_px_mi_add + sum(sapply(px_mi_add, FUN=function(i){
                sum(sapply(i, FUN=function(j){
                  if (is.null(j)){
                    return(0)
                  }else{
                    return(j$diff)
                  }
                }))
              }))
              
              px_mi_add <- lapply(px_mi_add, FUN=function(i){
                lapply(i, FUN=function(j){
                  if (!is.null(j)){
                    return(j$rho)
                  }
                })
              })
            }
            
            for (j in 1:N_MI){
              Z_MI_grouping_j <- Z_MI_grouping[[j]]
              simple_j <- all(lengths(Z_MI_grouping_j) == 1)
              add_j <- px_mi_add[[j]]
              receive_j <- lapply(add_j, FUN=function(i){attr(i, 'receive')})
              if (all(sapply(px_mi_add[[j]], is.null))){
                next
              }
              add_j_int <- lapply(add_j, FUN=function(i){i[1]})
              add_j <- lapply(add_j, FUN=function(i){i[-1]})
              prior_vi_mean_j <- vi_mi_mean[[j]]
              prior_vi_var_j <- vi_mi_var[[j]]
              for (k in 1:2){
                if (!is.null(add_j[[k]])){
                  if (simple_j){
                    # Subtract off the adjustment from the mean of the MI
                    vi_mi_mean[[j]][[k]] <- sweep(
                      prior_vi_mean_j[[k]], MARGIN = 2, FUN='-',
                      STATS = add_j[[k]])
                    # *Add* the adjustment to the RE of "l"
                    pos_l <- outer_px_add_map[[j]][[-k]]$position
                    vi_alpha_mean[pos_l] <- vi_alpha_mean[pos_l] +
                      as.vector(prior_vi_mean_j[[-k]] %*% add_j[[k]]) -
                      add_j_int[[k]]
                    vi_beta_mean['(Intercept)',] <- vi_beta_mean['(Intercept)',] +
                      add_j_int[[k]]
                  }else{
                    add_j_rho <- matrix(add_j[[k]], nrow = mi_d_j[j])
                    add_j_rho <- lapply(1:ncol(add_j_rho), FUN=function(i){
                      add_j_rho[,i]
                    })
                    # Subtract off the adjustment from the mean of the MI
                    vi_mi_mean[[j]][Z_MI_grouping_j[[k]]] <- mapply(
                      prior_vi_mean_j[Z_MI_grouping_j[[k]]], add_j_rho, SIMPLIFY = FALSE, FUN=function(i,j){
                      sweep(i, MARGIN = 2, FUN = '-', STATS = j)
                    })
                    # *Add* the adjustment to the RE of "l"
                    pos_l <- outer_px_add_map[[j]][[receive_j[[k]]]]$position
                    vi_alpha_mean[pos_l] <- vi_alpha_mean[pos_l] +
                      as.vector(prior_vi_mean_j[[receive_j[[k]]]] %*% Reduce('+', add_j_rho)) -
                      add_j_int[[k]]
                    vi_beta_mean['(Intercept)',] <- vi_beta_mean['(Intercept)',] +
                      add_j_int[[k]]
                  }
                }
              }
              if (all(sapply(add_j, is.null) == FALSE)){
                browser()
                
                vi_beta_mean['(Intercept)',] <- vi_beta_mean['(Intercept)',] -
                  as.numeric(t(add_j[[1]]) %*% add_j[[2]])
              }
              rm(simple_j)
            }
            if (nrow(vi_alpha_mean) > 0){
              variance_by_alpha_jg <- calculate_expected_outer_alpha(
                L = vi_alpha_decomp,
                alpha_mu = as.vector(vi_alpha_mean),
                re_position_list = outer_alpha_RE_positions
              )
              vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
            }
            
          }
          
        }

        vi_mi_sigma_outer_alpha <- get_bilinear_outer(
          vi_mi_mean, vi_mi_var, vi_mi_diag, mi_prior_type)
        
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
            
            vi_sigma_alpha = vi_sigma_alpha, 
            vi_a_b_jp = vi_a_b_jp,
            vi_sigma_outer_alpha = vi_sigma_outer_alpha,
            vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
            
            log_det_beta_var = log_det_beta_var, 
            log_det_alpha_var = log_det_alpha_var,
            log_det_joint_var = log_det_joint_var,
            
            vi_beta_decomp = vi_beta_decomp, 
            vi_alpha_decomp = vi_alpha_decomp,
            vi_joint_decomp = vi_joint_decomp,
            
            do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, 
            vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
            choose_term,
            # Multiplicative Interaction
            do_huangwand_mi = do_huangwand_mi,
            any_RE = any_RE, mi_prior_type = mi_prior_type,
            any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
            vi_mi_diag = vi_mi_diag, 
            vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
            vi_mi_sigma_alpha = vi_mi_sigma_alpha,
            vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
            vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
            vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
            vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
            mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
            mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
            mi_iw_prior_constant = mi_iw_prior_constant,
            mi_d_j = mi_d_j, mi_g_j = mi_g_j
          )

          diff_ELBO <- prop.ELBO$ELBO - prior.ELBO$ELBO
          verify_diff <- 
            abs(diff_ELBO - change_px_mi_add - change_px_mi_rotation) > 
            1e-6
          print('MI-PX')
          print(round(change_px_mi_add + change_px_mi_rotation, 4))
          
          if (diff_ELBO < -1e-6){
            browser()
          }else if (verify_diff){
            browser()
          }
          prior.ELBO <- prop.ELBO
        }
      }
      
      if (any_RE){

        # Do a simple mean adjusted expansion.
        # Get the mean of each random effect.
        
        vi_mu_j <- t(M_prime) %*% vi_alpha_mean
        
        meat_Bj <- bdiag(mapply(vi_sigma_alpha, vi_sigma_alpha_nu, d_j, 
                                SIMPLIFY = FALSE, FUN = function(phi, nu, d) {
                                  inv_phi <- solve(phi)
                                  sigma.inv <- nu * inv_phi
                                  return(sigma.inv)
                                }))
        
        proj_vi_mu_j <- B_j %*% solve(t(B_j) %*% meat_Bj %*% B_j) %*% t(B_j) %*% meat_Bj %*% vi_mu_j
        
        # Remove the "excess mean" mu_j from each random effect \alpha_{j,g}
        # and add the summd mass back to the betas.
        vi_alpha_mean <- vi_alpha_mean - M_prime_one %*% proj_vi_mu_j
        vi_beta_mean <- vi_beta_mean + t(M_mu_to_beta) %*% proj_vi_mu_j
        
        variance_by_alpha_jg <- calculate_expected_outer_alpha(
          L = vi_alpha_decomp,
          alpha_mu = as.vector(vi_alpha_mean), 
          re_position_list = outer_alpha_RE_positions
        )
        
        vi_sigma_outer_alpha <- variance_by_alpha_jg$outer_alpha
        
      }
      if (parameter_expansion == "mean"){accept.PX <- TRUE}
    }  

    quiet_rho <- control$quiet_rho

    if (any_RE & parameter_expansion %in% c("translation", "diagonal") & skip_translate == FALSE & any_Mprime) {
      
      attempted_expansion <- attempted_expansion + 1
      
      if (debug_px){
        prior.ELBO <- calculate_ELBO(family = family,
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
                                    
                                    vi_sigma_alpha = vi_sigma_alpha, 
                                    vi_a_b_jp = vi_a_b_jp,
                                    vi_sigma_outer_alpha = vi_sigma_outer_alpha,
                                    vi_beta_mean = vi_beta_mean, vi_alpha_mean = vi_alpha_mean,
                                    
                                    log_det_beta_var = log_det_beta_var, 
                                    log_det_alpha_var = log_det_alpha_var,
                                    log_det_joint_var = log_det_joint_var,
                                    
                                    vi_beta_decomp = vi_beta_decomp, 
                                    vi_alpha_decomp = vi_alpha_decomp,
                                    vi_joint_decomp = vi_joint_decomp,
                                    
                                    do_huangwand = do_huangwand, vi_a_a_jp = vi_a_a_jp, 
                                    vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
                                    choose_term,
                                    # Multiplicative Interaction
                                    do_huangwand_mi = do_huangwand_mi,
                                    any_RE = any_RE, mi_prior_type = mi_prior_type,
                                    any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
                                    vi_mi_diag = vi_mi_diag, 
                                    vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
                                    vi_mi_sigma_alpha = vi_mi_sigma_alpha,
                                    vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
                                    vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
                                    vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
                                    vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
                                    mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
                                    mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
                                    mi_iw_prior_constant = mi_iw_prior_constant,
                                    mi_d_j = mi_d_j, mi_g_j = mi_g_j
                                    
        )
      }
      
      if (!quiet_rho){cat('r')}
      
      if (do_timing){
        tic('px_r')
      }
      
      if (any(!spline_REs)){
        raw_R <- R_ridge <- vecR_ridge_new(L = vi_alpha_decomp[,nonspline_positions], pg_mean = diag(diag_vi_pg_mean),
                                           mapping_J = mapping_J, d = d_j[!spline_REs],
                                           store_id = store_id, store_re_id = store_re_id,
                                           store_design = store_design, 
                                           diag_only = (factorization_method == 'strong'))
      }else{
        raw_R <- R_ridge <- matrix(0, ncol = 0, nrow = 0)
      }

      if (factorization_method == 'weak'){
        stop('no Translation PX for weak yet...')
      }
      
      if (!quiet_rho){cat('r')}

      if (any(!spline_REs)){
        R_design <- vecR_design(alpha_mu = as.vector(vi_alpha_mean), Z = mapping_new_Z, 
                                M = Mmap, mapping_J = mapping_J, d = d_j[!spline_REs],
                                start_z = start_base_Z)
      }else{
        R_design <- matrix(0, nrow = N, ncol = 0)
      }

      if (sum(spline_REs)){
        R_spline_design <- sapply(cyclical_pos[spline_REs], FUN=function(i){
          as.vector(Z[,i,drop=F] %*% vi_alpha_mean[i,])
        })
        
        R_spline_ridge <- sapply(cyclical_pos[spline_REs], FUN=function(s){vi_alpha_decomp[,s, drop = F]})
        R_spline_ridge <- Diagonal(x =mapply(R_spline_ridge, cyclical_pos[spline_REs], FUN=function(V, pos){
          sum(vi_pg_mean * cpp_zVz(Z = drop0(Z[,pos,drop=F]), V = as(V, 'generalMatrix')))
        }))
        # Manually convert "ddiMatrix" to "generalMatrix" so doesn't fail on
        # old versions of "Matrix" package.
        if (inherits(R_spline_ridge, 'ddiMatrix')){
          R_spline_ridge <- diag(R_spline_ridge)
          R_spline_ridge <- sparseMatrix(
            i = seq_len(length(R_spline_ridge)),
            j = seq_len(length(R_spline_ridge)),
            x = R_spline_ridge)
        }else{
          R_spline_ridge <- as(R_spline_ridge, 'generalMatrix')
        }
      }else{
        R_spline_ridge <- drop0(matrix(0, nrow = 0, ncol = 0))
        R_spline_design <- matrix(nrow = nrow(X), ncol = 0)
      }
      
      
      if (do_timing){
        toc(quiet = quiet_time, log = TRUE)
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
        
        if (any_mi){
          offset_bilinear <- get_bilinear_mean(Z_MI, vi_mi_mean, Z_MI_grouping)
        }else{
          offset_bilinear <- 0
        }

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
          offset_bilinear = offset_bilinear,
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
        
        if (any(!spline_REs)){
          update_expansion_R <- mapply(split(update_expansion_XR[-1:-(p.X + sum(spline_REs))], 
                                             rep(1:(number_of_RE - sum(spline_REs)), d_j[!spline_REs]^2)), d_j[!spline_REs], 
                                       SIMPLIFY = FALSE, FUN=function(i,d){matrix(i, nrow = d)})
        }
        
      }
      
      if (do_timing){
        toc(quiet = quiet_time, log = TRUE)
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
        if (length(est_rho) > 0){
          diff_rho <- max(abs(est_rho - stationary_rho))
        }else{diff_rho <- 0}
        if (!is.na(px_improve) & (diff_rho < 1e-6) & (max(abs(est_rho_spline - 1)) < 1e-6) ){
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
        
        if (any(!spline_REs)){
          old_update_expansion_R <- update_expansion_R
          update_expansion_R <- lapply(d_j, FUN=function(i){diag(i)})
          update_expansion_R[!spline_REs] <- old_update_expansion_R
          rm(old_update_expansion_R)
        }else{
          update_expansion_R <- as.list(rep(NA, length(spline_REs)))
        }
        update_expansion_R[spline_REs] <- lapply(update_expansion_splines, FUN=function(i){matrix(i)})
        
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

      check_Rblock <- bdiag(mapply(update_expansion_R, g_j, FUN=function(i,g){bdiag(lapply(1:g, FUN=function(k){i}))}))
      if (max(abs(check_Rblock - update_expansion_Rblock)) != 0){
        warning('Error in creating parameter expansion; check that ELBO increases monotonically.')
      }
      
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
            choose_term,
            # Multiplicative Interaction
            do_huangwand_mi = do_huangwand_mi,
            any_RE = any_RE, mi_prior_type = mi_prior_type,
            any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
            vi_mi_diag = vi_mi_diag, 
            vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
            vi_mi_sigma_alpha = vi_mi_sigma_alpha,
            vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
            vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
            vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
            vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
            mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
            mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
            mi_iw_prior_constant = mi_iw_prior_constant,
            mi_d_j = mi_d_j, mi_g_j = mi_g_j
        )
      }
      if (!quiet_rho){cat('d')}

      # If debugging, check whether the change in ELBO
      # from the profiled objective agrees with the 
      # change from the actual ELBO.
      if (debug_px){
        ELBO_diff <- prop.ELBO$ELBO - prior.ELBO$ELBO
        if (!is.na(px_improve)){
          if (abs(ELBO_diff - px_improve) > 1e-6){
            browser()
            warning('PX does not agree with debug.')
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
        toc(quiet = quiet_time, log = TRUE)
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
      vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
      # Multiplicative Interaction
      do_huangwand_mi = do_huangwand_mi,
      any_RE = any_RE, mi_prior_type = mi_prior_type,
      any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
      vi_mi_diag = vi_mi_diag, 
      vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
      vi_mi_sigma_alpha = vi_mi_sigma_alpha,
      vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
      vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
      vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
      vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
      mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
      mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
      mi_iw_prior_constant = mi_iw_prior_constant,
      mi_d_j = mi_d_j, mi_g_j = mi_g_j
    )

    if (do_timing) {
      toc(quiet = quiet_time, log = TRUE)
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
      
      if (family == 'linear'){
        squarem_list[[squarem_counter]]$vi_sigmasq_a <- vi_sigmasq_a
        squarem_list[[squarem_counter]]$vi_sigmasq_b <- vi_sigmasq_b
      }
      if (any_mi & update_mi){
        if (it == 1){warning('Slow vi_mi_marg_decomp for SQUAREM')}
        squarem_list[[squarem_counter]]$vi_mi_mean <- vi_mi_mean
        squarem_list[[squarem_counter]]$vi_mi_marg_decomp <- vi_mi_marg_decomp
        vi_mi_marg_decomp <- lapply(vi_mi_marg_decomp, FUN=function(i){
          lapply(i, FUN=function(j){
            o <- t(apply(j, MARGIN = 1, FUN=function(k){
              m <- matrix(k, sqrt(length(k)))
              if (all(m == 0)){
               return(as.vector(rep(0, nrow(m)^2))) 
              }else{
                as.vector(chol(crossprod(m)))
              }
            }))
            if (ncol(j) == 1){
              o <- t(o)
            }
            return(o)
          })
        })
        squarem_list[[squarem_counter]]$vi_mi_marg_decomp <- mapply(vi_mi_marg_decomp, vi_mi_diag, SIMPLIFY = FALSE, FUN=function(i,diag_i){
          lapply(i, FUN=function(j){
            j[,diag_i] <- log(j[,diag_i])
            if (any(is.na(j))){browser()}
            return(j * VEM_scalar)
          })
        })
        squarem_list[[squarem_counter]]$vi_mi_a_b_jp <- vi_mi_a_b_jp
        squarem_list[[squarem_counter]]$vi_mi_sigma_alpha <- vi_mi_sigma_alpha
      }
      
      if (squarem_counter %% 3 == 0){
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
           vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
           any_RE = any_RE, mi_prior_type = mi_prior_type,
           any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
           do_huangwand_mi = do_huangwand_mi,
           vi_mi_diag = vi_mi_diag, 
           vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
           vi_mi_sigma_alpha = vi_mi_sigma_alpha,
           vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
           vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
           vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
           vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
           mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
           mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
           mi_iw_prior_constant = mi_iw_prior_constant,
           mi_d_j = mi_d_j, mi_g_j = mi_g_j
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
        
        if (!any_RE){
          re_remove <- match(c('vi_a_b_jp', 'vi_sigma_alpha', 'vi_alpha_mean', 'vi_alpha_L_nonpermute'), squarem_par)
          squarem_par <- squarem_par[-re_remove]
          squarem_type <- squarem_type[-re_remove]
          squarem_structure <- squarem_structure[-re_remove]
          
        }
        
        if (any_mi & update_mi){
          squarem_par <- c(squarem_par, 'vi_mi_mean', 'vi_mi_marg_decomp',
                           'vi_mi_a_b_jp', 'vi_mi_sigma_alpha')
          squarem_type <- c(squarem_type, 'real', 'chol_vec', 'positive', 'matrix')
          if (mi_prior_type %in% c('centered', 'shared')){
            squarem_structure <- c(squarem_structure, 'nested_list', 
                                   'nested_list', 'list', 'list')
          }else if (mi_prior_type %in% c('separate')){
            squarem_structure <- c(squarem_structure, 'nested_list', 
                                   'nested_list', 'nested_list', 'nested_list')
            
          }else{stop('invalid mi_prior_type')}
          if (!do_huangwand_mi){
            squarem_type <- squarem_type[!grepl(squarem_par, pattern='vi_mi_a_b_jp')]
            squarem_structure <- squarem_structure[!grepl(squarem_par, pattern='vi_mi_a_b_jp')]
            squarem_par <- squarem_par[!grepl(squarem_par, pattern='vi_mi_a_b_jp')]
          }
          if (VEM_ONLY){
            squarem_type <- squarem_type[!grepl(squarem_par, pattern='vi_mi_marg_decomp')]
            squarem_structure <- squarem_structure[!grepl(squarem_par, pattern='vi_mi_marg_decomp')]
            squarem_par <- squarem_par[!grepl(squarem_par, pattern='vi_mi_marg_decomp')]
          }
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
          squarem_par <- c(squarem_par, 'vi_sigmasq_b')
          squarem_type <- c(squarem_type, 'positive')
          squarem_structure <- c(squarem_structure, 'vector')
        }
        
        remove_hw_b <- FALSE
        if (!do_huangwand){
          squarem_type <- squarem_type[!grepl(squarem_par, pattern='vi_a_b_jp')]
          squarem_structure <- squarem_structure[!grepl(squarem_par, pattern='vi_a_b_jp')]
          squarem_par <- squarem_par[!grepl(squarem_par, pattern='vi_a_b_jp')]
        }else{
          if (remove_hw_b){
            squarem_type <- squarem_type[!grepl(squarem_par, pattern='vi_a_b_jp')]
            squarem_structure <- squarem_structure[!grepl(squarem_par, pattern='vi_a_b_jp')]
            squarem_par <- squarem_par[!grepl(squarem_par, pattern='vi_a_b_jp')]
          }
        }
        
        remove_c <- FALSE
        if (remove_c | family %in% 'linear'){
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
              }else if (s_str == 'nested_list'){
                out <- lapply(i[[s_par]], FUN=function(j){
                  lapply(j, FUN=function(k){squarem_prep_function(k,s_type)})
                })
              }else if (s_str == 'list'){
                out <- lapply(i[[s_par]], FUN=function(j){squarem_prep_function(j, s_type)})
              }else{stop('unrecognized structure')}
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
            }else if (s_str == 'nested_list'){
              r <- mapply(squarem_list[[2]][[s_par]], 
                          squarem_list[[1]][[s_par]], SIMPLIFY = FALSE,
                          FUN=function(i,j){mapply(i,j, SIMPLIFY = FALSE,
                                                   FUN=function(a,b){
                            a-b
                          })})
              d2 <- mapply(squarem_list[[3]][[s_par]], 
                          squarem_list[[2]][[s_par]], SIMPLIFY = FALSE,
                          FUN=function(i,j){mapply(i,j, SIMPLIFY = FALSE,
                                                   FUN=function(a,b){
                                                     a-b
                                                   })})
              v <- mapply(d2, r, SIMPLIFY = FALSE, FUN=function(i,j){
                mapply(i,j, SIMPLIFY = FALSE, FUN=function(a,b){a-b})
              })
              norm_sq_r <- sum(sapply(r, FUN=function(i){sapply(i, FUN=function(j){sum(j^2)})}))
              norm_sq_v <- sum(sapply(v, FUN=function(i){sapply(i, FUN=function(j){sum(j^2)})}))
              max_d <- max(sapply(v, FUN=function(i){sapply(i, FUN=function(j){max(abs(j))})}))

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
               }else if (s_str == 'nested_list'){
                 
                 prop_squarem <- mapply(i$first, i$second, 
                    i$r, i$v, SIMPLIFY = FALSE, FUN=function(i_1, s_1, r_1, v_1){
                      out <- mapply(i_1, s_1, r_1, v_1, SIMPLIFY = FALSE,
                        FUN=function(i__,s__,r__,v__){
                          return(i__ - 2 * s_alpha * r__ + s_alpha^2 * v__)
                      })
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
             }else if (s_str == 'nested_list'){
               out <- lapply(i, FUN=function(j){
                 lapply(j, FUN=function(k){squarem_unprep_function(k,s_type)})
               })
             }else if (s_str == 'list'){
               out <- lapply(i, FUN=function(j){squarem_unprep_function(j, s_type)})
             }else{stop('unrecognized structure')}
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
              prop_ELBOargs$log_det_beta_var <- 2 * 
                sum(log(diag(prop_squarem$vi_beta_L_nonpermute)))
              prop_squarem$vi_beta_decomp <- 
                prop_squarem$vi_beta_L_nonpermute %*% t(squarem_list[[1]]$vi_beta_LP)
            }

            if (!any_RE){
              # pass
            }else if (squarem_type[squarem_par == 'vi_alpha_L_nonpermute'] == 'lu'){
              prop_ELBOargs$log_det_alpha_var <- prop_squarem$vi_alpha_L_nonpermute$logdet_M
              prop_squarem$vi_alpha_decomp <- prop_squarem$vi_alpha_L_nonpermute$M              
            }else{
              prop_ELBOargs$log_det_alpha_var <- 2 * 
                sum(log(diag(prop_squarem$vi_alpha_L_nonpermute)))
              prop_squarem$vi_alpha_decomp <- 
                prop_squarem$vi_alpha_L_nonpermute %*% t(squarem_list[[1]]$vi_alpha_LP)
            }
            
            squarem_par <- c(squarem_par, 'log_det_alpha_var', 'log_det_beta_var')
            squarem_par <- c(squarem_par, 'vi_alpha_decomp', 'vi_beta_decomp')
            if (!any_RE){
              squarem_par <- setdiff(squarem_par,
                c('log_det_alpha_var', 'vi_alpha_decomp'))
            }
          }
          
          if (any_mi & update_mi & ('vi_mi_marg_decomp' %in% squarem_par)){
            recons_var <- mapply(prop_squarem$vi_mi_marg_decomp, mi_d_j, vi_mi_diag, SIMPLIFY = FALSE, 
              FUN=function(i, dim_i, diag_i){
                
              if (dim_i == 1){
                i <- lapply(i, FUN=function(j){
                  exp(j)^2 * VEM_scalar
                })
                lndet_i <- sapply(i, FUN=function(l){sum(log(l))})
                if (VEM_scalar == 0){
                  lndet_i <- 0
                }
                # i[[1]] <- exp(i[[1]])
                # i[[2]] <- exp(i[[2]])
                # i[[1]] <- i[[1]]^2
                # i[[2]] <- i[[2]]^2
                # lndet_i <- sum(sapply(i, FUN=function(l){sum(log(l))}))
              }else{
                i <- lapply(i, FUN=function(j){
                  j[,diag_i] <- exp(j[,diag_i])
                  j <- decomp_to_var_rowwise(j, dim_i, get_lndet = TRUE)
                  return(j)
                })
                lndet_i <- sapply(i, `[[`, 'lndet') * VEM_scalar
                i <- lapply(i, FUN=function(i){i$var * VEM_scalar})
                # i[[1]][,diag_i] <- exp(i[[1]][,diag_i])
                # i[[2]][,diag_i] <- exp(i[[2]][,diag_i])
                # i[[1]] <- decomp_to_var_rowwise(i[[1]], dim_i, get_lndet = TRUE)
                # i[[2]] <- decomp_to_var_rowwise(i[[2]], dim_i, get_lndet = TRUE)
                # lndet_i <- i[[1]]$lndet + i[[2]]$lndet
                # i <- list(i[[1]]$var, i[[2]]$var)
              }
              return(list(var = i, lndet = lndet_i))
            })
            
            prop_squarem$vi_mi_lndet <- lapply(recons_var, `[[`, 'lndet')
            prop_squarem$vi_mi_var <- lapply(recons_var, `[[`, "var")
            prop_squarem$vi_mi_marg_decomp <- NULL
            squarem_par <- c(squarem_par, 'vi_mi_var', 'vi_mi_lndet')
          }else{
            prop_squarem$vi_mi_lndet <- vi_mi_lndet
            prop_squarem$vi_mi_var <- vi_mi_var
            prop_squarem$vi_mi_marg_decomp <- NULL
            squarem_par <- c(squarem_par, 'vi_mi_var', 'vi_mi_lndet')
          }
          if (family == 'negbin'){
            
            if (vi_r_method == 'VEM'){
              
              prop_ELBOargs$vi_r_mean <- exp(prop_squarem$vi_r_mu)

            } else if (vi_r_method %in% c('Laplace', 'delta')){
              
              prop_ELBOargs$vi_r_mean <- exp(prop_squarem$vi_r_mu + prop_squarem$vi_r_sigma/2)
            
            } 
            
            if (factorization_method != 'weak'){
              prop_joint_var <- cpp_dense_zVz(X, as.matrix(prop_squarem$vi_beta_decomp)) +
                rowSums((Z %*% t(prop_squarem$vi_alpha_decomp))^2)
              # prop_joint_var <- rowSums((X %*% t(prop_squarem$vi_beta_decomp))^2) + 
              #   rowSums((Z %*% t(prop_squarem$vi_alpha_decomp))^2)
            }else{
              prop_joint_var <-  cpp_zVz(Z = joint.XZ, 
                V = as(prop_squarem$vi_joint_decomp, "generalMatrix")) 
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
          
          if (any_RE){
            prop_variance_by_alpha_jg <- calculate_expected_outer_alpha(
              L = (prop_squarem$vi_alpha_decomp), 
              alpha_mu = as.vector(prop_squarem$vi_alpha_mean), 
              re_position_list = outer_alpha_RE_positions)
            prop_ELBOargs[['vi_sigma_outer_alpha']] <- prop_variance_by_alpha_jg$outer_alpha
            squarem_par <- c(squarem_par, 'vi_sigma_outer_alpha')
          }
          
          if (any_mi & update_mi){
            prop_ELBOargs[['vi_mi_sigma_outer_alpha']] <-
              get_bilinear_outer(vi_mi_mean = prop_squarem$vi_mi_mean,
                                 vi_mi_var = prop_squarem$vi_mi_var,
                                 vi_mi_diag = vi_mi_diag,
                                 mi_prior_type = mi_prior_type)
            squarem_par <- c(squarem_par, 'vi_mi_sigma_outer_alpha')
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
            
            if (family %in% 'binomial'){
              prop_vi_pg_mean <- prop_ELBOargs$vi_pg_b / (2 * prop_ELBOargs$vi_pg_c) * tanh(prop_ELBOargs$vi_pg_c / 2)
              fill_zero <- which(abs(prop_ELBOargs$vi_pg_c) < 1e-6)
              if (length(fill_zero) > 0){
                prop_vi_pg_mean[fill_zero] <- prop_ELBOargs$vi_pg_b[fill_zero]/4
              }
              prop_ELBOargs[['vi_pg_mean']] <- prop_vi_pg_mean
              squarem_par <- c(squarem_par, 'vi_pg_mean')
            }else{stop('Set up SQUAREM for other family')}

          }else if (!(family %in% 'linear')){
            
            if (any_mi){stop('....')}
            if (family != 'binomial'){stop('check squarem for non-binomial case')}
            
            if (factorization_method %in% c("weak", "collapsed")) {
              joint_quad <- cpp_zVz(Z = joint.XZ, 
                  V = as(prop_ELBOargs$vi_joint_decomp, "generalMatrix")) 
              if (family == 'negbin'){
                joint_quad <- joint_quad + prop_ELBOargs$vi_r_sigma
              }
              prop_ELBOargs$vi_pg_c <- sqrt(as.vector(X %*% prop_ELBOargs$vi_beta_mean + Z %*% prop_ELBOargs$vi_alpha_mean - prop_ELBOargs$vi_r_mu)^2 + joint_quad)
            } else {
              # beta_quad <- rowSums((X %*% t(prop_ELBOargs$vi_beta_decomp))^2)
              beta_quad <- cpp_dense_zVz(X, as.matrix(prop_ELBOargs$vi_beta_decomp))
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
        
        if (elbo_squarem$ELBO >= elbo_init$ELBO){
          if (!quiet_rho){cat('SUCCESS')}
          squarem_success <- squarem_success + 1
          squarem.ELBO <- elbo_squarem
          final.ELBO <- elbo_squarem
          
          if (freeze_RE & any_mi){
            squarem_par <- squarem_par[!grepl(squarem_par, pattern='nonperm')]
          }
          
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
               vi_a_nu_jp = vi_a_nu_jp, vi_a_APRIOR_jp = vi_a_APRIOR_jp,
               # Multiplicative Interaction
               do_huangwand_mi = do_huangwand_mi,
               any_RE = any_RE, mi_prior_type = mi_prior_type,
               any_mi = any_mi, Z_MI = Z_MI, Z_MI_grouping = Z_MI_grouping,
               vi_mi_diag = vi_mi_diag, 
               vi_mi_sigma_outer_alpha = vi_mi_sigma_outer_alpha,
               vi_mi_sigma_alpha = vi_mi_sigma_alpha,
               vi_mi_sigma_alpha_nu = vi_mi_sigma_alpha_nu,
               vi_mi_mean = vi_mi_mean, vi_mi_var = vi_mi_var, vi_mi_lndet = vi_mi_lndet,
               vi_mi_a_a_jp = vi_mi_a_a_jp,  vi_mi_a_APRIOR_jp = vi_mi_a_APRIOR_jp,
               vi_mi_a_b_jp = vi_mi_a_b_jp, vi_mi_a_nu_jp = vi_mi_a_nu_jp,
               mi_prior_sigma_alpha_nu = mi_prior_sigma_alpha_nu, 
               mi_prior_sigma_alpha_phi = mi_prior_sigma_alpha_phi,
               mi_iw_prior_constant = mi_iw_prior_constant,
               mi_d_j = mi_d_j, mi_g_j = mi_g_j
          )
          if (test_ELBO$ELBO != elbo_squarem$ELBO){stop('SQUAREM misalignment')}
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
      toc(quiet = quiet_time, log = T)
      tic("Final Cleanup")
    }

    if (debug_ELBO & it != 1) {
      if (any_mi){
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
        
      }else{
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
      }
      update_ELBO$it <- it
      store_ELBO <- rbind(store_ELBO, update_ELBO)
    } else {
      final.ELBO$step <- NA
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

    change_beta_mean <- max(abs(vi_beta_mean - lagged_beta_mean))
    
    if (any_RE){
      change_alpha_mean <- max(abs(vi_alpha_mean - lagged_alpha_mean))
      unlist_vi <- c(unlist(lapply(vi_sigma_alpha, as.vector)), unlist(vi_a_b_jp))
    }else{
      unlist_vi <- double(0)
      change_alpha_mean <- 0
    }
    
    if (debug_ELBO){
      svi <- data.frame(t(as.vector(unlist_vi)))
      svi$it <- it
      store_vi <- rbind(store_vi,  svi)
    }
    
    if (any_RE){
      change_sigma_mean <- mapply(vi_sigma_alpha, lagged_sigma_alpha, FUN = function(i, j) {
        max(abs(i - j))
      })
    }else{
      change_sigma_mean <- 0
    }

    if (factorization_method %in% c("weak", "collapsed")) {
      change_joint_var <- 0 # change_joint_var <- max(abs(vi_joint_decomp - lagged_joint_decomp))
      change_alpha_var <- change_beta_var <- 0
    } else {
      change_joint_var <- 0
      if (any_RE){
        change_alpha_var <- max(abs(vi_alpha_decomp - lagged_alpha_decomp))
      }else{
        change_alpha_var <- 0
      }
      change_beta_var <- max(abs(vi_beta_decomp - lagged_beta_decomp))
    }

    if (any_mi){
      change_mi_vi_mean <- sum(mapply(vi_mi_mean, lagged_vi_mi_mean, FUN=function(i,j){
        sum(mapply(i,j, FUN=function(i,j){max(abs(i-j))}))
      }))
      change_mi_vi_var <- sum(mapply(vi_mi_var, lagged_vi_mi_var, FUN=function(i,j){
        sum(mapply(i,j, FUN=function(i,j){max(abs(i-j))}))
      }))
      if (mi_prior_type %in% c('separate')){
        change_mi_sigma_mean <- sum(mapply(vi_mi_sigma_alpha, lagged_vi_mi_sigma_alpha, FUN = function(i, j) {
          max(mapply(i,j, FUN=function(i_l, j_l){max(abs(i_l - j_l))}))
        }))
      }else{
        change_mi_sigma_mean <- sum(mapply(vi_mi_sigma_alpha, lagged_vi_mi_sigma_alpha, FUN = function(i, j) {
          max(abs(i - j))
        }))
      }
    }else{
      change_mi_vi_mean <- 0
      change_mi_vi_var <- 0
      change_mi_sigma_mean <- 0
    }
    change_vi_r_mu <- vi_r_mu - lagged_vi_r_mu

    if (do_timing) {
      toc(quiet = quiet_time, log = T)
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
      if (any_mi){
        store_mi_mean[it,] <- unlist(vi_mi_mean)
        store_mi_var[it,] <- unlist(lapply(vi_mi_var, FUN=function(i){lapply(i, as.matrix)}))
        store_mi_sigma[it,] <- unlist(lapply(vi_mi_sigma_alpha, FUN=function(i){
          lapply(i, as.matrix)
        }))
        if (do_huangwand_mi){
          store_mi_hw[it,] <- unlist(vi_mi_a_b_jp)
          colnames(store_mi_hw) <- names(unlist(vi_mi_a_b_jp))
        }
      }
    }
    change_all <- data.frame(change_alpha_mean, change_beta_mean, 
        t(change_sigma_mean), change_alpha_var, 
        change_beta_var, change_joint_var,
        change_mi_vi_mean, change_mi_vi_var,
        change_mi_sigma_mean,
        change_vi_r_mu)
    
    print(it)
    print(change_all)
    
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
    
    if (any_mi){
      lagged_vi_mi_mean <- vi_mi_mean
      lagged_vi_mi_var <- vi_mi_var
      lagged_vi_mi_sigma_alpha <- vi_mi_sigma_alpha
      lagged_vi_mi_lndet <- vi_mi_lndet
    }
    
    lagged_vi_r_mu <- vi_r_mu
    lagged_ELBO <- final.ELBO$ELBO
  }
  if (it == iterations) {
    message(paste0("Ended without Convergence after ", it, " iterations : ELBO change of ", round(change_elbo[1], abs(floor(log(tolerance_elbo) / log(10))))))
  }

  if (parameter_expansion %in% c("translation", "diagonal")) {
    final.ELBO$accepted_PX <- accepted_times / attempted_expansion
  }

  rownames(vi_beta_mean) <- colnames(X)
  
  output <- list(
    beta = list(mean = vi_beta_mean),
    ELBO = final.ELBO, 
    ELBO_trajectory = store_ELBO,
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
    if (any_mi){
      store_mi_mean <- store_mi_mean[1:it,,drop=F]
      store_mi_var <- store_mi_var[1:it,,drop=F]
      store_mi_sigma <- store_mi_sigma[1:it,,drop=FALSE]
      if (do_huangwand_mi){
        store_mi_hw <- store_mi_hw[1:it,,drop=F]
      }else{
        store_mi_hw <- NULL
      }
    }else{
      store_mi_mean <- store_mi_var <- store_mi_hw <- store_mi_sigma <- NULL
    }
    output$parameter_trajectory <- list(beta = store_beta,
                                        alpha = store_alpha,
                                        sigma = store_sigma,
                                        hw = store_hw,
                                        mi_mean = store_mi_mean,
                                        mi_var = store_mi_var,
                                        mi_sigma = store_mi_sigma,
                                        mi_hw = store_mi_hw)
  }
  if (factorization_method %in% c("weak", "collapsed")) {
    output$joint <- list(decomp_var = vi_joint_decomp)
  }
  if (control$return_data) {
    output$data <- list(X = X, Z = Z, y = y, trials = trials)
  }
  
  output$formula <- list(formula = formula, 
     re = re_fmla, fe = fe_fmla,
     interpret_gam = parse_formula,
     tt = tt, fe_Xlevels = fe_Xlevels,
     fe_contrasts = fe_contrasts, fe_terms = fe_terms)
  
  if (any_RE){
    output$alpha$dia.var <- unlist(lapply(variance_by_alpha_jg$variance_jg, FUN = function(i) {
      as.vector(sapply(i, diag))
    }))
  }
  
  output$beta$var <- t(vi_beta_decomp) %*% vi_beta_decomp
  output$beta$decomp_var <- vi_beta_decomp

  if (any_mi){
    
    vi_mi_mean <- mapply(vi_mi_mean, Z_MI_attr, SIMPLIFY = FALSE, FUN=function(i,j){
      mapply(i, j$levels, SIMPLIFY = FALSE, FUN=function(i_dat, i_level){
        rownames(i_dat) <- i_level
        return(i_dat)
      })
    })
    
    vi_mi_var <- mapply(vi_mi_var, Z_MI_attr, SIMPLIFY = FALSE, FUN=function(i,j){
      mapply(i, j$levels, SIMPLIFY = FALSE, FUN=function(i_dat, i_level){
        rownames(i_dat) <- i_level
        return(i_dat)
      })
    })
    

    output$mi <- list(
      mean = vi_mi_mean,
      var = vi_mi_var,
      sigma = list(cov = vi_mi_sigma_alpha, df = vi_mi_sigma_alpha_nu)
    )
    if (do_huangwand_mi){
      output$mi$hw <- list(a = vi_mi_a_a_jp, b = vi_mi_a_b_jp) 
    }
    
  }
  if (family == "negbin") {
    output$ln_r <- list(mu = vi_r_mu, sigma = vi_r_sigma, method = vi_r_method)
  }
  if (do_huangwand){
    output$hw <- list(a = vi_a_a_jp, b = vi_a_b_jp)
  }

  output_lp <- 0
  if (any_RE){
    output_lp <- as.vector(X %*% vi_beta_mean + Z %*% vi_alpha_mean - vi_r_mu)
  }
  if (any_mi){
    output_lp <- output_lp + get_bilinear_mean(Z_MI, vi_mi_mean, Z_MI_grouping)
  }
  
  output$internal_parameters <- list(
    it_used = it, it_max = iterations,
    cyclical_pos = cyclical_pos,
    lp = output_lp,
    parameter.change = change_all,
    parameter.vi = store_vi,
    parameter.path = store_parameter_traj,
    spline = list(attr = Z.special.attr, size = Z.special.size),
    missing_obs = missing_obs, N = nrow(X),
    acceleration = list(accept.PX = accept.PX, 
      squarem_success = squarem_success, debug_PX_ELBO = debug_PX_ELBO),
    names_of_RE = names_of_RE, d_j = d_j, g_j = g_j
  )

  if (any_RE){
    MAVB_parameters <- list(
      M_mu_to_beta = M_mu_to_beta,
      M_prime = M_prime,
      M_prime_one = M_prime_one,
      B_j = B_j,
      outer_alpha_RE_positions = outer_alpha_RE_positions,
      d_j = d_j, g_j = g_j
    )
    output$internal_parameters$MAVB_parameters <- MAVB_parameters
    output$alpha$var <- variance_by_alpha_jg$variance_jg
    output$alpha$decomp_var <- vi_alpha_decomp
  }
  output$timing <- tic_summary
  class(output) <- "vglmer"
  
  output$init_model <- lme4:::namedList(
    vi_sigma_alpha, vi_sigma_outer_alpha,
    vi_beta_mean, vi_alpha_mean, log_det_beta_var,
    log_det_alpha_var, vi_beta_decomp, vi_alpha_decomp,
    vi_joint_decomp, log_det_joint_var, vi_a_b_jp,
    vi_alpha_L_nonpermute, vi_alpha_LP,
    vi_beta_L_nonpermute, vi_beta_LP,
    vi_alpha_L_nonpermute, variance_by_alpha_jg,
    vi_joint_L_nonpermute, vi_joint_LP
  )
  return(output)
}

#' Control for vglmer estimation
#'
#' This function controls various estimation options for \code{vglmer}.
#'
#' @param iterations Default of 1000; this sets the maximum number of iterations
#'   used in estimation.
#' @param factorization_method Factorization assumption for the variational
#'   approximation. Default of \code{"strong"}, i.e. a fully factorized model.
#'   Described in detail in Goplerud (2022a). \code{"strong"}, \code{"partial"},
#'   and \code{"weak"} correspond to Schemes I, II, and III respectively in that
#'   paper.
#' @param prior_variance Prior distribution on the random effect variance
#'   \eqn{\Sigma_j}. Options are \code{hw}, \code{jeffreys}, \code{mean_exists},
#'   \code{uniform}, and \code{gamma}. The default (\code{hw}) is the Huang-Wand
#'   (2013) prior whose hyper-parameters are \eqn{\nu_j} = 2 and \eqn{A_{j,k}} =
#'   5. Otherwise, the prior is an Inverse Wishart with the following parameters
#'   where \eqn{d_j} is the dimensionality of the random effect \eqn{j}.
#'   \itemize{
#'   \item mean_exists: \eqn{IW(d_j + 1, I)}
#'   \item jeffreys: \eqn{IW(0, 0)}
#'   \item uniform: \eqn{IW(-[d_j+1], 0)}
#'   \item limit: \eqn{IW(d_j - 1, 0)}
#'   }
#'   Estimation may fail if an improper prior (\code{jeffreys}, \code{uniform},
#'   \code{limit}) is used.
#' @param tolerance_elbo Default (\code{1e-8}) sets a convergence threshold if
#'   the change in the ELBO is below the tolerance.
#' @param tolerance_parameters Default (\code{1e-5}) sets a convergence
#'   threshold that is achieved if no parameter changes by more than the
#'   tolerance from the prior estimated value.
#' @param parameter_expansion Default of \code{"translation"}  (see Goplerud
#'   2022b). Valid options are \code{"translation"}, \code{"mean"}, or
#'   \code{"none"}. \code{"mean"} should be employed if \code{"translation"} is
#'   not enabled or is too computationally expensive. For negative binomial
#'   estimation or any estimation where \code{factorization_method != "strong"},
#'   only \code{"mean"} and \code{"none"} are available.
#' @param px_method When code \code{parameter_expansion="translation"}, default
#'   (\code{"dynamic"}) tries a one-step late update and, if this fails, a
#'   numerical improvement by L-BFGS-B. For an Inverse-Wishart prior on
#'   \eqn{\Sigma_j}, this is set to \code{"osl"} that only attempts a
#'   one-step-late update.
#' @param px_numerical_it Default of 10; if L-BFGS_B is needed for a parameter
#'   expansion, this sets the number of steps used.
#' @param hw_inner If \code{prior_variance="hw"}, this sets the number of
#'   repeated iterations between estimating \eqn{\Sigma_j} and \eqn{a_{j,k}}
#'   variational distributions at each iteration. A larger number approximates
#'   jointly updating both parameters. Default (10) typically performs well.
#' @param force_whole Default (\code{TRUE}) requires integers for observed
#'   outcome for binomial or count models. \code{FALSE} allows for fractional
#'   responses.
#' @param vi_r_method Default (\code{"VEM"}) uses a variational EM algorithm for
#'   updating \eqn{r} if \code{family="negbin"}. This assumes a point mass
#'   distribution on \eqn{r}. A number can be provided to fix \eqn{r}. These are
#'   the only available options.
#' @param init Default (\code{"EM_FE"}) initializes the mean variational
#'   parameters for \eqn{q(\beta, \alpha)} by setting the random effects to zero
#'   and estimating the fixed effects using a short-running EM algorithm.
#'   \code{"EM"} initializes the model with a ridge regression with a guess as
#'   to the random effect variance. \code{"random"} initializes the means
#'   randomly. \code{"zero"} initializes them at zero.
#' @param debug_param Default (\code{FALSE}) does not store parameters before
#'   the final iteration. Set to \code{TRUE} to debug convergence issues.
#' @param debug_ELBO Default (\code{FALSE}) does not store the ELBO after each
#'   parameter update. Set to \code{TRUE} to debug convergence issues.
#' @param quiet_rho Default (\code{FALSE}) does not print information about
#'   parameter expansions. Set to \code{TRUE} to debug convergence issues.
#' @param debug_px Default (\code{FALSE}) does not store information about
#'   whether parameter expansion worked. Set to \code{TRUE} to convergence
#'   issues.
#' @param linpred_method Default (\code{"joint"}) updates the mean parameters
#'   for the fixed and random effects simultaneously. This can improve the speed
#'   of estimation but may be costly for large datasets; use \code{"cyclical"}
#'   to update each parameter block separately.
#' @param print_prog Default (\code{NULL}) prints a \code{"."} to indicate once
#'   5\% of the total iterations have elapsed. Set to a positive integer
#'   \code{int} to print a \code{"."} every \code{int} iterations.
#' @param quiet Default (\code{FALSE}) does not print intermediate output about
#'   convergence. Set to \code{TRUE} to debug.
#' @param return_data Default (\code{FALSE}) does not return the original
#'   design. Set to \code{TRUE} to debug convergence issues.
#' @param verbose_time Default (\code{FALSE}) does not print the time elapsed
#'   for each parameter update. Set to \code{TRUE}, in conjunction with
#'   \code{do_timing=TRUE}, to see the time taken for each parameter update.
#' @param do_timing Default (\code{FALSE}) does not estimate timing of each
#'   variational update; \code{TRUE} requires the package \code{tictoc}.
#' @param do_SQUAREM Default (\code{TRUE}) accelerates estimation using SQUAREM
#'   (Varadhan and Roland 2008).
#' @param verify_columns Default (\code{FALSE}) \bold{does not} verify that all
#'   columns are drawn from the data.frame itself versus the environment. Set to
#'   \code{TRUE} to debug potential issues.
#' @param mi_init Initialization for MI: "random" or "nakajima"
#' @param mi_parameter_expansion Parameter-Expansion for MI: "none" or "rotate"
#' @param mi_prior_type Prior type: "centered", "separated", "shared", "fixed",
#'   or "partial_fixed"
#' @param mi_prior_variance "hw" or "mean_exists"
#' @return This function returns a named list with class \code{vglmer_control}.
#'   It is passed to \code{vglmer} in the argument \code{control}. This argument
#'   only accepts objects created using \code{vglmer_control}.
#' 
#' @references 
#' Goplerud, Max. 2022a. "Fast and Accurate Estimation of Non-Nested Binomial
#' Hierarchical Models Using Variational Inference." \emph{Bayesian Analysis}.
#' 17(2): 623-650.
#'
#' Goplerud, Max. 2022b. "Re-Evaluating Machine Learning for MRP Given the
#' Comparable Performance of (Deep) Hierarchical Models." Working Paper.
#'
#' Huang, Alan, and Matthew P. Wand. 2013. "Simple Marginally Noninformative
#' Prior Distributions for Covariance Matrices." \emph{Bayesian Analysis}.
#' 8(2):439-452.
#'
#' Varadhan, Ravi, and Christophe Roland. 2008. "Simple and Globally Convergent
#' Methods for Accelerating the Convergence of any EM Algorithm."
#' \emph{Scandinavian Journal of Statistics}. 35(2): 335-353.
#' @export
vglmer_control <- function(iterations = 1000,
   prior_variance = "hw",
   factorization_method = c("strong", "partial", "weak"),
   parameter_expansion = "translation", do_SQUAREM = TRUE, 
   tolerance_elbo = 1e-8, tolerance_parameters = 1e-5,
   force_whole = TRUE, print_prog = NULL,
   do_timing = FALSE, verbose_time = FALSE,
   return_data = FALSE, linpred_method = "joint",
   vi_r_method = "VEM", verify_columns = FALSE,
   debug_param = FALSE, debug_ELBO = FALSE, debug_px = FALSE, 
   quiet = TRUE, quiet_rho = TRUE,
   freeze_mi_var = FALSE, 
   mi_init = c('nakajima', 'random'),
   mi_parameter_expansion = 'rotate',
   mi_prior_type = c('centered', 'separate', 'shared', 'fixed', 'partial_fixed'),
   mi_prior_variance = c('hw', 'mean_exists'),
   px_method = 'dynamic', px_numerical_it = 10,
   hw_inner = 10,
   init = "EM_FE") {
  
  mi_parameter_expansion <- match.arg(mi_parameter_expansion, choices = c('none', 'add', 'rotate'))
  if (any(mi_parameter_expansion == 'none')){
    if (length(mi_parameter_expansion) > 1){
      stop('mi_parameter_expansion must be "none" or some combination of "add" and "rotate".')
    }
    mi_parameter_expansion <- 'none'
  }
  mi_init <- match.arg(mi_init)
  mi_prior_variance <- match.arg(mi_prior_variance)
  mi_prior_type <- match.arg(mi_prior_type)
  factorization_method <- match.arg(factorization_method)
  prior_variance <- match.arg(prior_variance, 
    choices = c("hw", "mean_exists", "jeffreys", "limit", "uniform"))
  linpred_method <- match.arg(linpred_method, choices = c("joint", "cyclical", "solve_normal"))    
  parameter_expansion <- match.arg(parameter_expansion, choices = c("translation", "mean", "none"))
  # vi_r_method <- match.arg(vi_r_method, choices = c("VEM", "fixed", "Laplace", "delta"))
  init <- match.arg(init, choices = c("EM_FE", "EM", "random", "zero"))
  if (!is.null(print_prog)){
    if (print_prog < 0){stop('print_prog must be non-negative integer or NULL.')}
  }
  
  if (iterations < 0){stop('iterations must be positive integer')}
  if (tolerance_elbo < 0 | tolerance_parameters < 0){
    stop('tolerance for ELBO and parameters must be non-negative.')
  }
  
  if (factorization_method != "strong" & !(parameter_expansion %in% c("mean", "none"))){
    message('Setting parameter_expansion to mean for non-strong factorization')
    parameter_expansion <- 'mean'
  }
  if (prior_variance != 'hw' & px_method != 'OSL' & parameter_expansion %in% c('diagonal', 'translation')){
    px_method <- 'OSL'
    message('Setting px_method to "OSL" if translation & non-HW prior.')
  }

  output <- mget(ls())
  
  class(output) <- c("vglmer_control")
  return(output)
}
