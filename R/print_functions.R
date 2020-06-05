
#' Summarize output from vglmer
#' 
#' Allows the use of standard methods from lm or lmer to sumarize the posterior
#' output. The format for ranef differs from lmer and is shown below. 
#' 
#' \code{coef} and \code{vcov} return the mean and variance of the fixed effects. \code{fixef} is a
#' synonym for \code{coef}.
#' 
#' @name summary_vglmer
#' @param object Model fit using vglmer

#' @rdname summary_vglmer
#' @export
fixef.vglmer <- function(object, ...){
  return(as.vector(object$beta$mean))
}

#' @rdname summary_vglmer
ranef.vglmer <- function(object, ...){
  
  d_j <- object$internal_parameters$d_j
  g_j <- object$internal_parameters$g_j
  J <- length(d_j)
  
  vi_alpha_mean <- as.vector(object$alpha$mean)
  vi_alpha_var <- as.vector(object$alpha$dia.var)

  re_pos <- rep(1:J, d_j * g_j)
  
  vi_id <- gsub(rownames(object$alpha$mean), pattern='^.* @ .* @ ', replacement = '')
  vi_id <- split(vi_id, re_pos)
  vi_alpha_mean <- split(vi_alpha_mean, re_pos)
  vi_alpha_var <- split(vi_alpha_var, re_pos)

  vi_parsed <- mapply(d_j, g_j, vi_alpha_mean, vi_alpha_var, vi_id, object$internal_parameters$names_of_RE, SIMPLIFY = F, 
     FUN=function(d,g, mean_j, var_j, id_j, name_j){
        mat_id <- matrix(id_j, byrow = TRUE, nrow = g, ncol = d)
        mat_mean <- matrix(mean_j, byrow = TRUE, nrow = g, ncol = d)
        mat_var <- matrix(var_j, byrow = TRUE, nrow = g, ncol = d)
        colnames(mat_mean) <- colnames(mat_var) <- name_j
        id <- mat_id[,1]
        mat_mean <- data.frame(id, mat_mean, check.names = FALSE, stringsAsFactors = F)
        mat_var <- data.frame(id, mat_var, check.names = FALSE, stringsAsFactors = F)
        attributes(mat_mean)$'variance' <- mat_var
        return(mat_mean)
  })
  return(vi_parsed)
}

#' @rdname summary_vglmer
#' @export
coef.vglmer <- function(object, ...){
  return(as.vector(object$beta$mean))
}
#' @rdname summary_vglmer
#' @export
vcov.vglmer <- function(object, ...){
  return(as.matrix(object$beta$var))
}

#' @rdname summary_vglmer 
#' @param x Model fit using vglmer
#' @param ... Not used.
#' @export
print.vglmer <- function(x, ...){
  if (length(list(...)) > 0){'print.vglmer does not use ...'}
  it_used <- x$internal_parameters$it_used
  it_max <- x$internal_parameters$it_max
  final_param_change <- round(max(x$parameter.change), 6)
  final_ELBO_change <- round(tail(diff(x$ELBO_trajectory$ELBO), 1), 8)
  converged <- it_max != it_used
  p.X <- nrow(x$beta$mean)
  p.Z <- nrow(x$alpha$mean)
  J <- length(x$sigma$cov)
  
  cat(paste0('Formula: J = ', J,', |Z| = ', p.Z, ', |X| = ', p.X,'\n\n'))
  cat(paste(format(x$formula), collapse = '\n'))
  cat('\n\n')


  if (converged){
    cat(paste0('Model converged after ', it_used, ' iterations.'))
  }else{
    cat(paste0('Model *failed* to converge after ', it_max, ' iterations.'))
  } 
  cat('\n\n')
  cat(paste0('ELBO: ', round(x$ELBO[1],2), '\n\n'))
  cat(paste0('Factorization Method: ', x$factorization_method, '\n'))
  cat(paste0('Parameter Expansion: ', x$parameter_expansion, '\n\n'))
  cat(paste0('Largest Parameter Change at Convergence: ', format(final_param_change, nsmall = 2), '\n'))
  cat(paste0('ELBO Change at Convergence: ', format(final_ELBO_change, nsmall = 2), '\n'))
  
  invisible(list(paramater = final_param_change, ELBO = final_ELBO_change))
  
}

#' @importFrom lmtest coeftest
#' @export
summary.vglmer <- function(object, display_re = TRUE, ...){
  sum_obj <- coeftest(x = object)
  
  sum_sigma <- mapply(object$sigma$cov, object$sigma$df, SIMPLIFY = FALSE, FUN=function(a,b){fmt_IW_mean(a,b)})
  sum_sigma <- mapply(sum_sigma, object$internal.parameters$names_of_RE, SIMPLIFY = FALSE, FUN=function(i,j){
    rownames(i) <- colnames(i) <- j
    return(i)
  })
  re_names <- names(object$internal.parameters$names_of_RE)
  cat(paste0('Output from vglmer using ', object$factorization_method, ' factorization.\n'))
  cat('\nSummary of Fixed Effects\n')
  print(sum_obj)
  cat('\n')
  if (display_re){
    cat('Summary of Random Effects: Mean of Sigma_j (Variance)')
    for (v in seq_len(length(re_names))){
      cat('\n')
      cat(re_names[v])
      cat('\n')
      print(sum_sigma[[v]])
    }
  }
  invisible()
}

# Internal function to tidy-up
# inverse Wishart to extract mean
fmt_IW_mean <- function(Phi, nu, digits = 2){
  mean <- solve(Phi)/(nu - nrow(Phi) - 1)
  if (nu - nrow(Phi) - 1 < 0){
    return(matrix(NA, nrow = nrow(Phi), ncol = ncol(Phi)))
  }else{
    return(round(mean, digits))
  }
}
