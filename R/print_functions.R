
#' Generic Functions after Running vglmer
#'
#' Allows the use of standard methods from lm or lmer to sumarize the posterior
#' output.
#'
#' @return \code{coef} and \code{vcov} return the mean and variance of the fixed
#' effects. \code{fixef} is a synonym for \code{coef}. \code{ranef} extracts the
#' random effects in a similar, although slightly different, format to
#' \code{lme4}.
#'
#' \code{format_vglmer} collects the mean and variance of the fixed and random
#' effects into a data.frame
#'
#' @name vglmer-class
#' @param object Model fit using vglmer
#'

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
#' @export
ELBO <- function(object) UseMethod("ELBO")

ELBO.vglmer <- function(object){
  return(object$ELBO$ELBO)
}

#' @rdname vglmer-class
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
#' @param x Model fit using vglmer
#' @param ... Not used.
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
  final_param_change <- round(max(x$parameter.change), 6)
  final_ELBO_change <- round(tail(diff(x$ELBO_trajectory$ELBO), 1), 8)
  converged <- it_max != it_used
  p.X <- nrow(x$beta$mean)
  p.Z <- nrow(x$alpha$mean)
  J <- length(x$sigma$cov)

  cat(paste0("Formula: J = ", J, ", |Z| = ", p.Z, ", |X| = ", p.X, "\n\n"))
  cat(paste(format(x$formula), collapse = "\n\n"))
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
#' @param display_re Print summary of random effects. Default is TRUE
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
  cat(paste0("Output from vglmer using ", object$factorization_method, " factorization.\n"))
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
#' @export
formula.vglmer <- function(x, ...) {
  x$formula
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
  output <- bind_rows(beta.output, alpha.output)
  return(output)
}
