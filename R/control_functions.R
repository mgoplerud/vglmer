#' Control for vglmer estimation
#'
#' This function controls various estimation options for \code{vglmer}.
#'
#' @param iterations Default of 1000; this sets the maximum number of iterations
#'   used in estimation.
#' @param factorization_method Factorization assumption for the variational
#'   approximation. Default of \code{"strong"}, i.e. a fully factorized model.
#'   Described in detail in Goplerud (2022). \code{"strong"}, \code{"partial"},
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
#' 
#' @return This function returns a named list with class \code{vglmer_control}.
#'   It is passed to \code{vglmer} in the argument \code{control}. This argument
#'   only accepts objects created using \code{vglmer_control}.
#' 
#' @references 
#' 
#' Goplerud, Max. 2022. "Fast and Accurate Estimation of Non-Nested Binomial
#' Hierarchical Models Using Variational Inference." \emph{Bayesian Analysis}.
#' 17(2): 623-650.
#'
#' Goplerud, Max. 2024. "Re-Evaluating Machine Learning for MRP Given the
#' Comparable Performance of (Deep) Hierarchical Models." \emph{American
#' Political Science Review}. 118(1): 529-536.
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
                           px_method = 'dynamic', px_numerical_it = 10,
                           hw_inner = 10,
                           init = "EM_FE") {
  
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
  
  if (factorization_method != "strong" & parameter_expansion != "mean"){
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
