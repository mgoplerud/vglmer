
#' Generic Functions after Running vglmer
#'
#' \code{vglmer} uses many standard methods from \code{lm} and \code{lme4} with
#' limited changes. These provide summaries of the estimated variational
#' distributions.
#'
#' @details The accompanying functions are briefly described below. 
#' 
#' \code{coef} and \code{vcov} return the mean and variance of the fixed effects
#' (\eqn{\beta}). \code{fixef} returns the mean of the fixed effects.
#'
#' \code{ranef} extracts the random effects (\eqn{\alpha}) in a similar,
#' although slightly different format, to \code{lme4}. It includes the estimated
#' posterior mean and variance in a list of data.frames with one entry per
#' random effect \eqn{j}.
#' 
#' \code{fitted} extracts the estimated expected \emph{linear predictor}, i.e.
#' \eqn{E_{q(\theta)}[x_i^T \beta + z_i^T \alpha]} at convergence.
#' 
#' \code{summary} reports the estimates for all fixed effects as in \code{lm} as
#' well as some summaries of the random effects (if \code{display_re=TRUE}).
#' 
#' \code{format_vglmer} collects the mean and variance of the fixed and random
#' effects into a single data.frame. This is useful for examining all of the
#' posterior estimates simultaneously. \code{format_glmer} converts an object
#' estimated with \code{[g]lmer} into a comparable format.
#'
#' \code{ELBO} extracts the ELBO from the estimated model. \code{type} can be
#' set equal to \code{"trajectory"} to get the estimated ELBO at each iteration
#' and assess convergence.
#'   
#' \code{sigma} extracts the square root of the posterior mode of
#' \eqn{q(\sigma^2)} if a linear model is used.
#' 
#' \code{formula} extracts the formula associated with the \code{vglmer} object.
#' By default, it returns the formula provided. The fixed and random effects
#' portions can be extracted separately using the \code{form} argument.
#' 
#' @name vglmer-class
#' @param object Model fit using vglmer
#'
#' @return The functions here return a variety of different objects depending on
#'   the specific function. "Details" describes the behavior of each one. Their
#'   output is similar to the typical behavior for the corresponding generic
#'   functions.

#' @rdname vglmer-class
#' @export
fixef.vglmer <- function(object, ...) {
  out <- object$beta$mean
  rn <- rownames(out)
  out <- as.vector(out)
  names(out) <- rn
  return(out)
}

# Load fixef, ranef, sigma from lme4
#' @export
lme4::fixef
#' @export
lme4::ranef

#' @importFrom stats sigma

#' @rdname vglmer-class
#' @export
sigma.vglmer <- function(object, ...){
  #{\displaystyle \frac{\sqrt{2}}{2} \left(\frac{(2m-1)\Omega}{m}\right)^{1/2}}
  
  if (object$family != 'linear'){
    stop('sigma from vglmer is only defined for linear models')
  }
  if (length(list(...)) > 0){
    stop('... not used for sigma.vglmer')
  }
  naive_sigma <- with(object$sigmasq, sqrt(b/(a+1)))
  return(naive_sigma)
}

#' @rdname vglmer-class
#' @export
ranef.vglmer <- function(object, ...) {
  if (length(list(...)) > 0) {
    stop("... not used for ranef.vglmer")
  }

  d_j <- object$internal_parameters$d_j
  g_j <- object$internal_parameters$g_j
  J <- length(d_j)

  vi_alpha_mean <- as.vector(object$alpha$mean)
  vi_alpha_var <- as.vector(object$alpha$dia.var)

  re_pos <- rep(1:J, d_j * g_j)

  vi_id <- gsub(rownames(object$alpha$mean), pattern = "^.* @ .* @ ", replacement = "")
  vi_id <- split(vi_id, re_pos)
  vi_alpha_mean <- split(vi_alpha_mean, re_pos)
  vi_alpha_var <- split(vi_alpha_var, re_pos)

  vi_parsed <- mapply(d_j, g_j, vi_alpha_mean, vi_alpha_var, vi_id, object$internal_parameters$names_of_RE,
    SIMPLIFY = F,
    FUN = function(d, g, mean_j, var_j, id_j, name_j) {
      mat_id <- matrix(id_j, byrow = TRUE, nrow = g, ncol = d)
      mat_mean <- matrix(mean_j, byrow = TRUE, nrow = g, ncol = d)
      mat_var <- matrix(var_j, byrow = TRUE, nrow = g, ncol = d)
      colnames(mat_mean) <- colnames(mat_var) <- name_j
      id <- mat_id[, 1]
      mat_mean <- data.frame(id, mat_mean, check.names = FALSE, stringsAsFactors = F)
      mat_var <- data.frame(id, mat_var, check.names = FALSE, stringsAsFactors = F)
      attributes(mat_mean)$"variance" <- mat_var
      return(mat_mean)
    }
  )
  return(vi_parsed)
}

#' @rdname vglmer-class
#' @method coef vglmer
#' @export
coef.vglmer <- function(object, ...) {
  if (length(list(...)) > 0) {
    stop("... not used for coef.vglmer")
  }
  out <- as.vector(object$beta$mean)
  names(out) <- rownames(object$beta$mean)
  return(out)
}
#' @rdname vglmer-class
#' @export
vcov.vglmer <- function(object, ...) {
  if (length(list(...)) > 0) {
    stop("... not used for vcov.vglmer")
  }
  return(as.matrix(object$beta$var))
}

#' @rdname vglmer-class
#' @method fitted vglmer
#' @export
fitted.vglmer <- function(object, ...){
  if (length(list(...)) > 0) {
    stop("... not used for vcov.vglmer")
  }
  return(object$internal_parameters$lp)
}

#' @rdname vglmer-class
#' @param x Model fit using \code{vglmer}.
#' @param ... Not used; included to maintain compatibility with existing
#'   methods.
#' @method print vglmer
#' @export
print.vglmer <- function(x, ...) {
  if (length(list(...)) > 0) {
    "print.vglmer does not use ..."
  }
  N_obs <- x$internal_parameters$N
  missing_obs <- x$internal_parameters$missing_obs
  it_used <- x$internal_parameters$it_used
  it_max <- x$internal_parameters$it_max
  final_param_change <- round(max(x$internal_parameters$parameter.change), 6)
  final_ELBO_change <- round(tail(diff(x$ELBO_trajectory$ELBO), 1), 8)
  converged <- it_max != it_used
  p.X <- nrow(x$beta$mean)
  p.Z <- nrow(x$alpha$mean)
  J <- length(x$sigma$cov)

  cat(paste0("Formula: J = ", J, ", |Z| = ", p.Z, ", |X| = ", p.X, "\n\n"))
  cat(paste(format(formula(x, form = 'original')), collapse = "\n\n"))
  cat("\n\n")
  if (missing_obs > 0) {
    missing_info <- paste0("after ", missing_obs, " deleted because of missing data and")
  } else {
    missing_info <- " and"
  }
  cat(paste0("Model fit with ", N_obs, " observations", missing_info))
  if (converged) {
    cat(paste0(" converged after ", it_used, " iterations."))
  } else {
    cat(paste0(" *failed* to converge after ", it_max, " iterations."))
  }
  cat("\n\n")
  cat(paste0("ELBO: ", round(x$ELBO[1], 2), "\n\n"))
  cat(paste0("Factorization Method: ", x$control$factorization_method, "\n"))
  cat(paste0("Parameter Expansion: ", x$control$parameter_expansion, "\n\n"))
  cat(paste0("Largest Parameter Change at Convergence: ", formatC(final_param_change, format = "e", digits = 2), "\n"))
  cat(paste0("ELBO Change at Convergence: ", formatC(final_ELBO_change, format = "e", digits = 2), "\n"))


  invisible(list(paramater = final_param_change, ELBO = final_ELBO_change))
}

#' @rdname vglmer-class
#' @param display_re Default (\code{TRUE}) prints a summary of the
#'   random effects alongside the fixed effects.
#' @importFrom lmtest coeftest
#' @method summary vglmer
#' @export
summary.vglmer <- function(object, display_re = TRUE, ...) {
  sum_obj <- coeftest(x = object)

  sum_sigma <- mapply(object$sigma$cov, object$sigma$df, SIMPLIFY = FALSE, FUN = function(a, b) {
    fmt_IW_mean(a, b)
  })
  sum_sigma <- mapply(sum_sigma, object$internal_parameters$names_of_RE, SIMPLIFY = FALSE, FUN = function(i, j) {
    rownames(i) <- colnames(i) <- j
    return(i)
  })
  re_names <- names(object$internal_parameters$names_of_RE)
  cat(paste0("Output from vglmer using '", object$control$factorization_method, "' factorization.\n"))
  cat("\nSummary of Fixed Effects\n")
  print(sum_obj)
  cat("\n")
  if (display_re) {
    cat("Summary of Random Effects: Mean of Sigma_j (Variance)")
    for (v in seq_len(length(re_names))) {
      cat("\n")
      cat(re_names[v])
      cat("\n")
      print(sum_sigma[[v]], quote = FALSE)
    }
    cat("\n")
  }
  if (object$family == "negbin") {
    r_output <- object$r
    # fmt_r <- function(x){formatC(x, format = 'e', digits = 2)}
    fmt_r <- function(x) {
      round(x, digits = 2)
    }
    r_ci <- exp(r_output$mu + sqrt(2 * r_output$sigma) * erfinv(c(0.05, 0.95)))
    r_ci <- paste0("[", paste(fmt_r(r_ci), collapse = ", "), "]")
    r_mean <- fmt_r(exp(r_output$mu + r_output$sigma / 2))

    cat("Summary of Auxiliary Parameters:\n")
    cat("Dispersion Parameter r:\n")
    if (object$r$method == "VI") {
      cat(paste0("Mean (90% Interval): ", r_mean, " ", r_ci))
    } else {
      cat(paste0("Mean: ", r_mean))
    }
    cat("\n")
    cat("\n")
  }
  invisible()
}

#' @rdname vglmer-class
#' @param form Describes the type of formula to report:
#'   \code{"original"} returns the user input, \code{"fe"} returns the fixed
#'   effects only, \code{"re"} returns the random effects only.
#' @export
formula.vglmer <- function(x, form = "original", ...) {
  
  if (form == 'original'){
    x$formula$formula
  }else if (form == 'fe'){
    x$formula$fe
  }else if (form == 're'){
    x$formula$re
  }else{stop('form must be "original", "fe", or "re".')}
  
}

#' @importFrom stats qnorm
erfinv <- function(x) {
  qnorm((1 + x) / 2) / sqrt(2)
}

# Internal function to tidy-up
# inverse Wishart to extract mean
fmt_IW_mean <- function(Phi, nu, digits = 2) {
  mean <- solve(as.matrix(Phi)) / (nu - nrow(Phi) - 1)
  if (nu - nrow(Phi) - 1 < 0) {
    return(matrix(NA, nrow = nrow(Phi), ncol = ncol(Phi)))
  } else {
    return(formatC(mean, format = "e", digits = 2))
  }
}

#' @rdname vglmer-class
#' @export
format_vglmer <- function(object) {
  beta.output <- data.frame(name = rownames(object$beta$mean), mean = as.vector(object$beta$mean), var = diag(object$beta$var), stringsAsFactors = F)
  alpha.output <- data.frame(name = rownames(object$alpha$mean), mean = as.vector(object$alpha$mean), var = as.vector(object$alpha$dia.var), stringsAsFactors = F)
  output <- rbind(beta.output, alpha.output)
  return(output)
}

#' @rdname vglmer-class
#' @importFrom stats vcov
#' @export
format_glmer <- function(object) {
  
  output <- do.call('rbind', mapply(ranef(object), names(ranef(object)), SIMPLIFY = FALSE, FUN = function(i,j) {
    obj <- data.frame(
      var = as.vector(apply(attributes(i)$postVar, MARGIN = 3, FUN = function(i) {
        diag(i)
      })),
      mean = as.vector(t(as.matrix(i))),
      name = paste0(rep(colnames(i), nrow(i)), " @ ", rep(rownames(i), each = ncol(i))), stringsAsFactors = F
    )
    obj[[".re"]] <- j
    return(obj)
  }))
  output$name <- paste0(output[[".re"]], ' @ ', output[["name"]])
  output_fe <- data.frame(mean = fixef(object), var = diag(stats::vcov(object)))
  output_fe$name <- rownames(output_fe)
  output_fe[[".re"]] <- NA
  output <- rbind(output, output_fe)
  output <- output[, (names(output) != ".re")]
  
  rownames(output) <- NULL
  
  return(output)
}

#' @rdname vglmer-class
#' @param object Model fit using \code{vglmer}.
#' @param type Default (\code{"final"}) gives the ELBO at convergence.
#'   \code{"trajectory"} gives the ELBO estimated at each iteration. This is
#'   used to assess model convergence.
#' @export
ELBO <- function(object, type = c('final', 'trajectory')){
  
  type <- match.arg(type)

  if (type == 'final'){
    object$ELBO$ELBO
  }else{
    object$ELBO_trajectory$ELBO
  }
  
}
