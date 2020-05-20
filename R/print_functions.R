
#' Summarize output from vglmer
#' 
#' Allows the use of standard methods from lm or lmer to sumarize the posterior
#' output. The format for ranef differs from lmer and is shown below. 
#' 
#' \code{coef} and \code{vcov} return the mean and variance of the fixed effects. \code{fixef} is a
#' synonym for \code{coef}.
#' 
#' @name summary_vglmer
#' @param object Model fitted from vglmer

#' @rdname summary_vglmer
#' @export
fixef.vglmer <- function(object, ...){
  return(as.vector(object$beta$mean))
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

# Print status of vglmer fit
#' @export
print.vglmer <- function(x, ...){
  cat('Formula:\n\n')
  cat(paste(format(x$formula), collapse = '\n'))
  cat('\n\n')

  cat(paste0('ELBO: ', round(x$ELBO[1],2), '\n\n'))
  cat(paste0('Factorization Method: ', x$factorization_method, '\n'))
  cat(paste0('Parameter Expansion: ', x$parameter_expansion, '\n\n'))
  cat(paste0('Largest Parameter Change at Convergence: ', format(round(max(x$parameter.change), 6), nsmall = 2), '\n'))
  cat(paste0('ELBO Change at Convergence: ', '', '\n'))
  
  warning('Set up maxit, did it convergence, ELBO change, and and PX')
  invisible()
  
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
