#' Perform MAVB after fitting vglmer
#'
#' Given a model from vglmer, perform marginally augmented variational Bayes
#' (MAVB) to improve approximation quality; Goplerud (2022) provides details. At
#' present, it is only enabled for binomial models.
#'
#' This should only be used if the parameters of the model are of interest; to
#' use MAVB when generating predictions, use \link{predict_MAVB}.
#'
#' @param object Model fit using vglmer
#' @param samples Samples to draw from MAVB distribution.
#' @param var_px Default (Inf); variance of working prior. Higher is more
#'   diffuse and thus likely better.
#' @param verbose Show progress of MAVB.
#' @import CholWishart
#' @importFrom mvtnorm rmvnorm
#' 
#' @references 
#' Goplerud, Max. 2022. "Fast and Accurate Estimation of Non-Nested Binomial
#' Hierarchical Models Using Variational Inference." Bayesian Analysis. 17(2):
#' 623-650.
#' @export
MAVB <- function(object, samples, verbose = FALSE, var_px = Inf) {
  
  if (!inherits(object, "vglmer")) {
    stop("Must provide object from vglmer")
  }
  
  if (object$family != "binomial") {
    stop("MAVB only implemented for binomial at present.")
  }
  
  M_prime <- object$MAVB_parameters$M_prime
  M_prime_one <- object$MAVB_parameters$M_prime_one
  M_mu_to_beta <- object$MAVB_parameters$M_mu_to_beta
  B_j <- object$MAVB_parameters$B_j
  if (!isDiagonal(B_j)){
    stop('MAVB not set up for non-diagonal mean expansion; do all REs have a corresponding FE?')
  }else{
    if (!isTRUE(all.equal(B_j@x, rep(1, nrow(B_j))))){
      stop('B_j is diagonal but not identity matrix; do all REs have a corresponding FE?')
    }
  }
  d_j <- object$MAVB_parameters$d_j
  g_j <- object$MAVB_parameters$g_j
  outer_alpha_RE_positions <- object$MAVB_parameters$outer_alpha_RE_positions
  factorization_method <- object$control$factorization_method

  if (factorization_method == "weak") {
    decomp_joint <- object$joint
    joint_mean <- rbind(object$beta$mean, object$alpha$mean)
    p.XZ <- ncol(decomp_joint)
    p.X <- nrow(object$beta$mean)
  } else {
    decomp_varA <- object$alpha$decomp_var
    decomp_varB <- object$beta$decomp_var
    p.XZ <- ncol(decomp_varA) + ncol(decomp_varB)
    p.X <- nrow(object$beta$mean)
    p.Z <- nrow(object$alpha$mean)
  }

  MAVB_sims <- matrix(NA, nrow = samples, ncol = p.XZ)
  regen_store <- MAVB_diff <- matrix(NA, nrow = samples, ncol = sum(d_j))

  n_MAVB <- samples

  alpha_mean <- object$alpha$mean
  beta_mean <- object$beta$mean
  sigma_df <- object$sigma$df
  sigma_cov <- object$sigma$cov

  all_sigma <- mapply(sigma_df, sigma_cov, SIMPLIFY = FALSE, FUN = function(i, j) {
    rInvWishart(n = n_MAVB, df = i, Sigma = j)
  })

  for (it in 1:n_MAVB) {
    if (it %% 1000 == 0 & verbose) {
      message(".", appendLF = F)
    }
    # Sim from VARIATIONAL approximation to posterior.
    if (factorization_method == "weak") {
      sim_joint <- joint_mean + t(decomp_joint) %*% rnorm(p.XZ)
      sim_beta <- sim_joint[1:p.X, , drop = F]
      sim_a <- sim_joint[-1:-p.X, , drop = F]
    } else {
      sim_a <- alpha_mean + t(decomp_varA) %*% rnorm(p.Z)
      sim_beta <- beta_mean + t(decomp_varB) %*% rnorm(p.X)
    }
    sim_sigma <- lapply(all_sigma, FUN = function(i) {
      i[, , it]
    })
    # Working Prior
    if (var_px == Inf) {
      sim_px <- rep(0, sum(d_j))
    } else {
      sim_px <- rnorm(sum(d_j), 0, sd = sqrt(var_px))
    }
    # Transform t^{-1}_a(z) = w
    sim_atilde <- sim_a + M_prime_one %*% sim_px
    sim_btilde <- sim_beta - t(M_mu_to_beta) %*% sim_px
    # Draw sim_px AGAIN from its full conditional
    if (var_px == Inf) {
      # Use the LIMITING transition.
      var_redux <- as.matrix(bdiag(mapply(sim_sigma, g_j, SIMPLIFY = FALSE, FUN = function(S, g) {
        S / g
      })))
      # Get the MEAN not the SUM
      mean_redux <- t(M_prime) %*% sim_atilde
    } else {
      var_redux <- solve(diag(x = 1 / var_px, ncol = sum(d_j), nrow = sum(d_j)) + as.matrix(bdiag(mapply(sim_sigma, g_j, SIMPLIFY = FALSE, FUN = function(S, g) {
        solve(S) * g
      }))))
      # Use the SUM
      mean_redux <- var_redux %*% solve(bdiag(sim_sigma)) %*% t(M_prime_one) %*% sim_atilde
    }
    regen_px <- t(rmvnorm(1, mean_redux, var_redux))
    regen_store[it, ] <- regen_px
    MAVB_diff[it, ] <- regen_px - sim_px
    final_atilde <- sim_atilde - M_prime_one %*% regen_px
    final_btilde <- sim_btilde + t(M_mu_to_beta) %*% regen_px

    MAVB_sims[it, ] <- c(as.vector(final_btilde), as.vector(final_atilde))
  }
  colnames(MAVB_sims) <- c(rownames(object$beta$mean), 
                           rownames(object$alpha$mean))
  
  return(MAVB_sims)
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


#' @import lme4
get_RE_groups <- function(formula, data) {
  
  if (inherits(formula, 'formula')){
    bars <- findbars(formula)
  }else{
    bars <- formula
  }
  if (is.null(bars)){# Usually if only splines used, then NA.
    return(list(factor = NA, design = NA))
  }  
  
  barnames <- utils::getFromNamespace('barnames', 'lme4')
  names(bars) <- barnames(bars)
  
  fr <- data
  blist <- lapply(bars, simple_blist, fr, drop.unused.levels = F, reorder.vars = FALSE)
  blist <- lapply(blist, FUN=function(i){i[c('ff', 'mm')]})

  ff <- lapply(blist, FUN=function(i){i$ff})
  ff <- lapply(ff, FUN=function(i){match(i, levels(i))})
  mm <- lapply(blist, FUN=function(i){i$mm})
  return(list(factor = ff, design = mm))
}

#' @import lme4 
#' @importFrom utils getFromNamespace
simple_blist <- function(x, frloc, drop.unused.levels = TRUE, reorder.vars = FALSE) {
  frloc <- factorize(x, frloc)
  makeFac <- utils::getFromNamespace('makeFac', 'lme4')
  if (is.null(ff <- tryCatch(eval(substitute(makeFac(fac),
                                             list(fac = x[[3]])), frloc), error = function(e) NULL)))
    stop("couldn't evaluate grouping factor ", deparse(x[[3]]),
         " within model frame:", " try adding grouping factor to data ",
         "frame explicitly if possible", call. = FALSE)
  if (all(is.na(ff)))
    stop("Invalid grouping factor specification, ", deparse(x[[3]]),
         call. = FALSE)
  if (drop.unused.levels)
    ff <- factor(ff, exclude = NA)
  nl <- length(levels(ff))
  mm <- model.matrix(eval(substitute(~foo, list(foo = x[[2]]))),
                     frloc)
  if (reorder.vars) {
    
    colSort <- utils::getFromNamespace("colSort", "lme4")
    mm <- mm[colSort(colnames(mm)), ]
  }
  list(ff = ff, nl = nl, mm = mm, cnms = colnames(mm))
}


#' Variance of Rows or Columns of Matrices
#'
#' Base R implementation for variance. Analogue of rowMeans.
#' @name var_mat
#' @keywords internal
#' @param matrix Matrix of numeric inputs.
rowVar <- function(matrix) {
  apply(matrix, MARGIN = 1, var)
}

#' @importFrom stats var
#' @rdname var_mat
colVar <- function(matrix) {
  apply(matrix, MARGIN = 2, var)
}

#' Get samples from GLMER
#'
#' Order samples from glmer to match names from vglmer.
#'
#' @param glmer object fitted using glmer
#' @param samples number of samples to draw
#' @param ordering order of output
#' @keywords internal
#' @importFrom stats rnorm
custom_glmer_samples <- function(glmer, samples, ordering) {
  fmt_glmer <- format_glmer(glmer)

  glmer_samples <- mapply(fmt_glmer$mean, fmt_glmer$var, FUN = function(m, v) {
    rnorm(samples, mean = m, sd = sqrt(v))
  })
  colnames(glmer_samples) <- fmt_glmer$name

  glmer_samples <- glmer_samples[, match(ordering, colnames(glmer_samples))]
  return(glmer_samples)
}

#' Draw samples from the (non-MAVB) posterior
#' 
#' @param object Model fit using vglmer.
#' @param samples Number of samples to draw from the variational distribution.
#' @param verbose Print progress of drawing the samples.
#' @export
posterior_samples.vglmer <- function (object, samples, verbose = FALSE) 
{
  if (!inherits(object, "vglmer")) {
    stop("Must provide object from vglmer")
  }
  
  M_prime <- object$MAVB_parameters$M_prime
  M_prime_one <- object$MAVB_parameters$M_prime_one
  M_mu_to_beta <- object$MAVB_parameters$M_mu_to_beta
  d_j <- object$MAVB_parameters$d_j
  g_j <- object$MAVB_parameters$g_j
  outer_alpha_RE_positions <- object$MAVB_parameters$outer_alpha_RE_positions
  factorization_method <- object$control$factorization_method
  if (factorization_method == "weak") {
    decomp_joint <- object$joint
    joint_mean <- rbind(object$beta$mean, object$alpha$mean)
    p.XZ <- ncol(decomp_joint)
    p.X <- nrow(object$beta$mean)
  }
  else {
    decomp_varA <- object$alpha$decomp_var
    decomp_varB <- object$beta$decomp_var
    p.XZ <- ncol(decomp_varA) + ncol(decomp_varB)
    p.X <- nrow(object$beta$mean)
    p.Z <- nrow(object$alpha$mean)
  }
  post_sims <- matrix(NA, nrow = samples, ncol = p.XZ)
  alpha_mean <- object$alpha$mean
  beta_mean <- object$beta$mean
  sigma_df <- object$sigma$df
  sigma_cov <- object$sigma$cov
  for (it in 1:samples) {
    if (it%%1000 == 0 & verbose) {
      message(".", appendLF = F)
    }
    if (factorization_method == "weak") {
      sim_joint <- joint_mean + t(decomp_joint) %*% rnorm(p.XZ)
      sim_beta <- sim_joint[1:p.X, , drop = F]
      sim_a <- sim_joint[-1:-p.X, , drop = F]
    }
    else {
      sim_a <- alpha_mean + t(decomp_varA) %*% rnorm(p.Z)
      sim_beta <- beta_mean + t(decomp_varB) %*% rnorm(p.X)
    }
    post_sims[it, ] <- c(as.vector(sim_beta), as.vector(sim_a))
  }
  
  colnames(post_sims) <- c(rownames(object$beta$mean), 
    rownames(object$alpha$mean))
  
  return(post_sims)
  
}
