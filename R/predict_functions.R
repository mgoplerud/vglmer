#' Predict after vglmer
#'
#' @description These functions calculate the estimated linear predictor using
#'   the variational distributions. \code{predict.vglmer} draws predictions
#'   using the estimated variational distributions; \code{predict_MAVB} does so
#'   using the MAVB procedure described in Goplerud (2022a).
#' @name vglmer_predict
#' @param object Model fit using \code{vglmer}.
#' @param newdata Dataset to use for predictions. It cannot be missing.
#' @param samples Number of samples to draw. Using \code{0} (default) gives the
#'   expectation of the linear predictor. A positive integer draws
#'   \code{samples} samples from the variational distributions and calculates
#'   the linear predictor.
#' @param type Default (\code{"link"}) returns the linear predictor;
#'   \code{"terms"} returns the predicted value for each random effect (or
#'   spline) separately as well as one that collects all fixed effects. At the
#'   moment, other options are not enabled.
#' @param samples_only Default (\code{FALSE}) returns the samples from the
#'   variational distributions, \bold{not} the prediction. Each row is a sample and
#'   each column is a parameter.
#' @param summary Default (\code{TRUE}) returns the mean and variance of the
#'   samples for each observation. \code{FALSE} returns a matrix of the sampled
#'   linear predictor for each observation. Each row is a sample and each column
#'   is an observation.
#' @param allow_missing_levels Default (\code{FALSE}) does not allow prediction
#'   for levels not observed in the original data. \code{TRUE} allows for
#'   prediction on unseen levels; the value of \code{0} (with no uncertainty) is
#'   used for the corresponding random effect.
#' @param ... Not used; included to maintain compatibility with existing
#'   methods.
#'
#' @examples
#' 
#' set.seed(123)
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
#' # Return 10 posterior draws of the linear predictor for each observation.
#' predict_MAVB(est_vglmer, newdata = sim_data, summary = FALSE, samples = 10)
#' # Predict with a new level; note this would fail if 
#' # allow_missing_levels = FALSE (the default)
#' predict(est_vglmer,
#'   newdata = data.frame(g = "AB", x = 0),
#'   allow_missing_levels = TRUE
#' )
#' @return This function returns an estimate of the linear predictor. The
#'   default returns the expected mean, i.e. \eqn{E_{q(\alpha,\beta)}[x_i^T
#'   \beta + z_i^T\alpha]}. If \code{samples > 0}, these functions return a
#'   summary of the prediction for each observation, i.e. the estimated mean and
#'   variance. If \code{summary = FALSE}, the sampled values of the linear
#'   predictor are returned as a matrix. \code{predict_MAVB} performs MAVB as
#'   described in Goplerud (2022a) before returning the linear predictor.
#'   
#'   If \code{allow_missing_levels = TRUE}, then observations with a new
#'   (unseen) level for the random effect are given a value of zero for that
#'   term of the prediction.
#' @importFrom stats delete.response terms na.pass
#' @export
predict.vglmer <- function(object, newdata, type = 'link',
                           samples = 0, samples_only = FALSE,
                           summary = TRUE, allow_missing_levels = FALSE, ...) {
  if (length(list(...)) > 0) {
    stop("... not used for predict.vglmer")
  }
  if (!(type %in% c('link', 'terms'))){
    stop('vglmer only uses "terms" and "link" for "type" in predict.')
  }
  newdata <- as.data.frame(newdata)
  rownames(newdata) <- as.character(1:nrow(newdata))

  parse_formula <- object$formula$interpret_gam
  if (!all(parse_formula$pred.names %in% colnames(newdata))){
    missing_columns <- setdiff(parse_formula$pred.names, colnames(newdata))
    stop(
      paste0('The following columns are missing from "newdata": ', 
        paste(missing_columns, collapse =', '))
    )
  }
  fmla <- formula(object, form = 'original')
  
  newdata_FE <- model.frame(delete.response(object$formula$fe_terms), 
      data = newdata, xlev = object$formula$fe_Xlevels, na.action = na.pass)
  X <- model.matrix(
    delete.response(object$formula$fe_terms), newdata_FE, 
      contrasts.arg = object$formula$fe_contrasts)
  
  orig_X_names <- rownames(object$beta$mean)
  if (!identical(colnames(X), orig_X_names)) {
    print(all.equal(colnames(X), orig_X_names))
    stop("Misaligned Fixed Effects")
  }

  mk_Z <- model.frame(delete.response(terms(object$formula$interpret_gam$fake.formula)), 
                      data = newdata, drop.unused.levels = TRUE)
  rownames_Z <- rownames(mk_Z)
  

  if (!is.null(object$formula$re) & (length(object$formula$re) > 0) ){
    
    # Extract the Z (Random Effect) design matrix.
    mk_Z <- mkReTrms(formula(object, form = 're'), mk_Z, reorder.terms = FALSE, reorder.vars = FALSE)
    Z <- t(mk_Z$Zt)
    
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
    
  }
  
  # Extract the Specials
  if (length(parse_formula$smooth.spec) > 0){
    base_specials <- length(parse_formula$smooth.spec)
    # Number of splines + one for each "by"...
    n.specials <- base_specials +
      sum(sapply(parse_formula$smooth.spec, FUN=function(i){i$by}) != "NA")
    
    
    Z.spline <- as.list(rep(NA, n.specials))
    Z.spline.size <- rep(NA, n.specials)
    Z.spline.attr <- object$internal_parameters$spline$attr
    
    special_counter <- 1
    store_spline_type <- rep(NA, n.specials)
    for (i in 1:base_specials){
      
      special_i <- parse_formula$smooth.spec[[i]]
      
      all_splines_i <- vglmer_build_spline(x = newdata[[special_i$term]], 
         knots = Z.spline.attr[[i]]$knots, 
         Boundary.knots =  Z.spline.attr[[i]]$Boundary.knots,
         by = newdata[[Z.spline.attr[[i]]$by]], outer_okay = TRUE,
         type = Z.spline.attr[[i]]$type, override_warn = TRUE,
         force_vector = TRUE)
      
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
        
        store_spline_type[special_counter] <- spline_counter
        
        spline_counter <- spline_counter + 1
        special_counter <- special_counter + 1
        
      }
      
    }

    Z.spline <- drop0(do.call('cbind', Z.spline))
    rownames(Z.spline) <- rownames(newdata)
    
    if (ncol(Z) > 0){
      Z.spline <- Z.spline[match(rownames(Z), rownames(Z.spline)),, drop = FALSE]
      Z <- drop0(cbind(Z, Z.spline))
    }else{
      Z <- Z.spline
    }
    
    if (!isTRUE(all.equal(names_of_RE, object$internal_parameters$names_of_RE))){
      stop('Names of REs do not match estimation data. This may occur when REs have to be re-named.')
    }
    
    if (!isTRUE(identical(object$internal_parameters$spline$size[store_spline_type %in% 1], 
                          Z.spline.size[store_spline_type  %in% 1]))){
      stop('Misalignment of splines in prediction.')
    }
    if (!isTRUE(identical(names_of_RE, object$internal_parameters$names_of_RE))){
      stop('Misalignment of spline names in prediction.')
    }
    
  }else{
    n.specials <- 0
    Z.spline.attr <- NULL
    Z.spline <- NULL
    Z.spline.size <- NULL
  }
  
  #####
  ### Confirm Alignment of the Z
  #####
  orig_Z_names <- rownames(object$alpha$mean)

  not_in_original_Z <- setdiff(fmt_names_Z, orig_Z_names)
  not_in_new_Z <- setdiff(orig_Z_names, fmt_names_Z)
  
  if (length(not_in_original_Z) > 0) {
    if (!allow_missing_levels) {
      stop("New levels not allowed unless allow_missing_levels = TRUE")
    }
  }

  # Select overlapping columns
  in_both <- intersect(fmt_names_Z, orig_Z_names)
  # Find the ones that are missing
  missing_cols <- setdiff(orig_Z_names, in_both)

  recons_Z <- Z[, match(in_both, fmt_names_Z), drop = F]
  if (length(missing_cols) > 0){
    # Create a matrix of zeros to pad the missing columns
    pad_zero <- sparseMatrix(i = 1, j = 1, x = 0, 
                             dims = c(nrow(Z), length(missing_cols)))
    colnames(pad_zero) <- missing_cols
    # Combine and then reorder to be lined-up correctly
    recons_Z <- cbind(recons_Z, pad_zero)
  }
  recons_Z <- recons_Z[, match(orig_Z_names, colnames(recons_Z)), drop = F]
    
  # Old method for prediction
  # in_both <- intersect(fmt_names_Z, orig_Z_names)
  # recons_Z <- drop0(sparseMatrix(i = 1, j = 1, x = 0, dims = c(nrow(Z), length(orig_Z_names))))
  # colnames(recons_Z) <- orig_Z_names
  # rownames(recons_Z) <- rownames_Z
  # recons_Z[, match(in_both, orig_Z_names)] <- Z[, match(in_both, fmt_names_Z)]
  
  # Check that the entirely missing columns match those not in the original
  checksum_align <- setdiff(not_in_new_Z,
    sort(names(which(colSums(recons_Z != 0) == 0))))
  if (length(checksum_align) > 0) {
    stop("Alignment Error")
  }
  
  Z <- recons_Z
  rm(recons_Z); gc()

  ####
  
  total_obs <- rownames(newdata)
  obs_in_both <- intersect(rownames(X), rownames(Z))

  if (type == 'terms'){
    if (samples != 0){stop('"terms" only enabled for samples=0.')}
    # Calculate the linear predictor separately for each random effect
    # (and fixed effects) and report a matrix of those predictions.
    
    X <- X[match(obs_in_both, rownames(X)), , drop = F]
    Z <- Z[match(obs_in_both, rownames(Z)), , drop = F]
    lp_FE <- as.vector(X %*% object$beta$mean)
    vi_alpha_mean <- object$alpha$mean
    lp_terms <- lapply(object$internal_parameters$cyclical_pos, FUN=function(i){
      as.vector(Z[,i,drop=F] %*% vi_alpha_mean[i,drop=F])
    })
    lp_terms <- do.call('cbind', lp_terms)
    colnames(lp_terms) <- names(object$internal_parameters$names_of_RE)
    lp_terms <- cbind('FE' = lp_FE, lp_terms)
    lp_terms <- lp_terms[match(total_obs, obs_in_both), , drop = F]
    gc()
    return(lp_terms)
    
  }else{
    XZ <- cbind(
      X[match(obs_in_both, rownames(X)), , drop = F],
      Z[match(obs_in_both, rownames(Z)), , drop = F]
    )
  }
  gc()
  
  factorization_method <- object$control$factorization_method
  if (is.matrix(samples)) {
    if (ncol(samples) != ncol(XZ)) {
      stop("Samples must be {m, ncol(Z) + ncol(X)}")
    }
    samples <- t(samples)
    only.lp <- FALSE
  } else {
    if (samples == 0) {
      only.lp <- TRUE
    } else {
      only.lp <- FALSE
    }
    if (factorization_method %in% c("strong", "partial")) {
      vi_alpha_mean <- object$alpha$mean
      vi_alpha_decomp <- object$alpha$decomp_var

      p.Z <- nrow(vi_alpha_mean)

      vi_beta_mean <- object$beta$mean
      vi_beta_decomp <- object$beta$decomp_var

      p.X <- nrow(vi_beta_mean)

      if (!only.lp) {
        sim_init_alpha <- matrix(rnorm(samples * p.Z), ncol = samples)
        sim_init_alpha <- t(vi_alpha_decomp) %*% sim_init_alpha
        sim_init_alpha <- sim_init_alpha + kronecker(vi_alpha_mean, t(matrix(1, samples)))

        sim_init_beta <- matrix(rnorm(samples * p.X), ncol = samples)
        sim_init_beta <- t(vi_beta_decomp) %*% sim_init_beta
        sim_init_beta <- sim_init_beta + kronecker(vi_beta_mean, t(matrix(1, samples)))
      } else {
        sim_init_alpha <- vi_alpha_mean
        sim_init_beta <- vi_beta_mean
      }
    } else if (factorization_method == "weak") {
      vi_alpha_mean <- object$alpha$mean
      p.Z <- nrow(vi_alpha_mean)

      vi_beta_mean <- object$beta$mean
      p.X <- nrow(vi_beta_mean)

      if (!only.lp) {
        vi_joint_decomp <- object$joint$decomp_var
        sim_init_joint <- matrix(rnorm(samples * (p.X + p.Z)), ncol = samples)
        sim_init_joint <- t(vi_joint_decomp) %*% sim_init_joint

        sim_init_beta <- sim_init_joint[1:p.X, , drop = F]
        sim_init_alpha <- sim_init_joint[-1:-p.X, , drop = F]

        rm(sim_init_joint)

        sim_init_alpha <- sim_init_alpha + kronecker(vi_alpha_mean, t(matrix(1, samples)))
        sim_init_beta <- sim_init_beta + kronecker(vi_beta_mean, t(matrix(1, samples)))
      } else {
        sim_init_alpha <- vi_alpha_mean
        sim_init_beta <- vi_beta_mean
      }
    } else {
      stop("")
    }

    samples <- rbind(sim_init_beta, sim_init_alpha)
    rm(sim_init_beta, sim_init_alpha)
  }

  if (samples_only) {
    return(t(samples))
  }

  lp <- XZ %*% samples
  if (summary) {
    if (!only.lp) {
      lp <- t(apply(lp, MARGIN = 1, FUN = function(i) {
        c(mean(i), var(i))
      }))
      lp <- data.frame(mean = lp[, 1], var = lp[, 2])
      lp <- lp[match(total_obs, obs_in_both), ]
      rownames(lp) <- NULL
    } else {
      
      if (ncol(lp) != 1){
        lp <- as.vector(lp)
      }else{
        lp <- as.vector(t(apply(lp, MARGIN = 1, FUN = function(i) {
          mean(i)
        })))
      }
      lp <- lp[match(total_obs, obs_in_both)]
      rownames(lp) <- NULL
    }
    return(lp)
  } else {
    lp <- lp[match(total_obs, obs_in_both), , drop = F]
    rownames(lp) <- NULL
    return(t(lp))
  }
}


#' @inheritParams MAVB
#' @inheritParams vglmer_predict
#' @rdname vglmer_predict
#' @export
predict_MAVB <- function(object, newdata, samples = 0, samples_only = FALSE,
                         var_px = Inf, summary = TRUE, allow_missing_levels = FALSE) {
  pxSamples <- MAVB(object = object, samples = samples, var_px = var_px)
  lp <- predict.vglmer(object,
                       newdata = newdata, samples = pxSamples, samples_only = samples_only,
                       summary = summary, allow_missing_levels = allow_missing_levels
  )
  return(lp)
}
