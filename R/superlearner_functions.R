#' SuperLearner with (Variational) Hierarchical Models
#' 
#' Functions to integrate \code{vglmer} (or \code{glmer}) into
#' \code{SuperLearner}. Most of the arguments are standard for
#' \code{SuperLearner} functions. See "Details" for how \code{formula} should be
#' used.
#'
#' @param Y From SuperLearner: The outcome in the training data set.
#' @param X From SuperLearner: The predictor variables in the training data.
#' @param newX From SuperLearner: The prediction variables in validation data.
#' @param formula The formula used for estimation, e.g. for \code{vglmer} or
#'   \code{[g]lmer}.
#' @param family Family as character: Only "gaussian" and "binomial" accepted.
#' @param id From SuperLearner: Passed from SuperLearner.
#' @param obsWeights From SuperLearner: Weights for each observation. Not
#'   permitted for \code{SL.vglmer}.
#' @param control Adjust the additional options for \code{vglmer} and
#'   \code{[g]lmer}. See the corresponding estimation functions for details.
#' @param object Used in \code{predict} for \code{SL.glmer} and
#'   \code{SL.vglmer}. A model estimated using either \code{SL.vglmer} or
#'   \code{SL.glmer}.
#' @param ... Not used; added for compatibility with existing methods.
#' @param learner Character vector of model from SuperLearner. See "Details"
#'   for how this is used.
#' @param env Environment to look for model. See "Details" for how this is used.
#' @name sl_vglmer
#' 
#' @details This documentation describes two types of function. They are
#'   described below.
#'   
#'   \bold{Formula with SuperLearner}: The \code{vglmer} package provides a way
#'   to estimate models that require or use formula with \code{SuperLearner}.
#'   This allows for a design to be passed that contains variables that are
#'   \emph{not} used in estimation. This can be used as follows (see
#'   "Examples"). One calls the function \code{add_formula_SL} around the quoted
#'   name of a \code{SuperLearner} model, e.g. \code{add_formula_SL(learner =
#'   "SL.knn")}. This creates a new model and predict function with the suffix
#'   \code{"_f"}. This \bold{requires} a formula to be provided.
#'   
#'   With this in hand, \code{"SL.knn_f"} can be passed to SuperLearner with the
#'   accompanying formula argument and thus one can compare models with
#'   different formula or design on the same ensemble.
#'   
#'   \bold{Estimating Hierarchical Models in SuperLearner}: Two methods for
#'   estimating hierarchical models are provided. One allows for the use of
#'   variational hierarchical models using \code{SL.vglmer}. Non-variational
#'   hierarchical models estimated via \code{lme4} can be used in an ensemble
#'   with \code{SL.glmer}. The accompanying predict methods are also provided.
#'   
#' @examples
#' 
#' if (requireNamespace('SuperLearner', quietly = TRUE)){
#' require(SuperLearner)
#' sim_data <- data.frame(
#'   x = rnorm(100),
#'   g = sample(letters, 100, replace = TRUE)
#' )
#' sim_data$y <- rbinom(nrow(sim_data), 
#'   1, plogis(runif(26)[match(sim_data$g, letters)]))
#' sim_data$g <- factor(sim_data$g)
#' sl_vglmer <- function(...){SL.vglmer(..., formula = y ~ x + (1 | g))}
#' SL.glm <- SuperLearner::SL.glm
#' add_formula_SL('SL.glm')
#' sl_glm_form <- function(...){SL.glm_f(..., formula = ~ x)}
#  
#' \donttest{
#'    SuperLearner::SuperLearner(
#'      Y = sim_data$y, family = 'binomial',
#'      X = sim_data[, c('x', 'g')],
#'      cvControl = list(V = 2),
#'      SL.library = c('sl_vglmer', 'sl_glm_form')
#'    )
#' }
#' }
#' @export
SL.vglmer <- function(Y, X, newX, formula, family, id, obsWeights, control = vglmer_control()) {
  if(!requireNamespace('vglmer', quietly = FALSE)) {stop("SL.vglmer requires the vglmer package, but it isn't available")} 
  
  if (is.character(formula)){
    formula <- as.formula(formula)
  }
  # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
  getResponseFromFormula = function(formula) {
    if (attr(terms(as.formula(formula))    , which = 'response'))
      all.vars(formula)[1]
    else
      NULL
  }
  rformula <- getResponseFromFormula(formula)
  if (!is.null(rformula)){
    if (rformula %in% names(X)){
      warning(paste0('Outcome "', rformula, '" seems to be in "X". This is likely ill-advised'))
    }
  }
  if ('...Y' %in% names(X)){
    stop('SL.vglmer cannot accept a column in "X" called "...Y". Please rename.')
  }
  if (!all(obsWeights == 1)){
    warning('SL.vglmer does not use weights')
  }
  if (family$family == 'binomial'){
    family <- 'binomial'
  }else if (family$family == 'gaussian'){
    family <- 'linear'
  }else{stop('Family must be binomial or Gaussian for SL.vglmer.')}
  X[['...Y']] <- Y
  formula <- update.formula(formula, '`...Y` ~ .')
  
  fit.vglmer <- vglmer::vglmer(formula, data = X, family = family, control = control)
  pred <- predict(fit.vglmer, newdata = newX, allow_missing_levels = TRUE)
  if (family == 'binomial'){
    pred <- plogis(pred)
  }else if (family != 'linear'){
    stop('SuperLearner not set up for non-linear, non-binomial families.')
  }
  fit <- list(object = fit.vglmer)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.vglmer")
  return(out)
}

#' @rdname sl_vglmer
#' @param newdata Model for predicting on test data
#' @param allow_missing_levels Allowing missing levels not in training data.
#' @export
predict.SL.vglmer <- function(object, newdata, allow_missing_levels = TRUE, ...){
  if(!requireNamespace('vglmer', quietly = FALSE)) {stop("SL.vglmer requires the vglmer package, but it isn't available")} 
  
  pred <- predict(object$object, newdata = newdata, allow_missing_levels = allow_missing_levels)
  if (object$object$family == 'binomial'){
    pred <- plogis(pred)
  }else if (object$object$family != 'linear'){
    stop('SuperLearner not set up for non-linear, non-binomial families.')
  }
  return(pred)
}

#' @importFrom stats predict
#' @rdname sl_vglmer
#' @export
SL.glmer <- function(Y, X, newX, formula, family, id, obsWeights, control = NULL) {
  if(!requireNamespace('lme4', quietly = FALSE)) {stop("SL.glmer requires the lme4 package, but it isn't available")} 
  
  if (is.character(formula)){
    formula <- as.formula(formula)
  }
  # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
  getResponseFromFormula = function(formula) {
    if (attr(terms(as.formula(formula))    , which = 'response'))
      all.vars(formula)[1]
    else
      NULL
  }
  rformula <- getResponseFromFormula(formula)
  if (!is.null(rformula)){
    if (rformula %in% names(X)){
      warning(paste0('Outcome "', rformula, '" seems to be in "X". This is likely ill-advised'))
    }
  }
  if ('...Y' %in% names(X)){
    stop('SL.glmer cannot accept a column in "X" called "...Y". Please rename.')
  }
  X[['...Y']] <- Y
  formula <- update.formula(formula, '`...Y` ~ .')
  environment(formula) <- environment()
  
  if (family$family == 'gaussian'){
    if (is.null(control)){
      control <- lmerControl()
    }
    fit.glmer <- lme4::lmer(formula, data = X, weights = obsWeights, control = control)
  }else{
    if (is.null(control)){
      control <- glmerControl()
    }
    fit.glmer <- lme4::glmer(formula, data = X, weights = obsWeights, family = family, control = control)
  }
  pred <- stats::predict(fit.glmer, newdata = newX, allow.new.levels = TRUE, type = 'response')
  fit <- list(object = fit.glmer)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.glmer")
  return(out)
}

#' @rdname sl_vglmer
#' @param allow.new.levels From \code{lme4}, allow levels in prediction that are
#'   not in the training data. Default (\code{TRUE}) for \code{SuperLearner}.
#' @export
predict.SL.glmer <- function(object, newdata, allow.new.levels = TRUE, ...){
  if(!requireNamespace('lme4', quietly = FALSE)) {stop("SL.glmer requires the lme4 package, but it isn't available")} 
  
  pred <- predict(object$object, newdata = newdata, allow.new.levels = allow.new.levels, type = 'response')
  return(pred)
}

#' @rdname sl_vglmer
#' @export
add_formula_SL <- function(learner, env = parent.frame()){

  base_learner <- get(learner, envir = env)
  base_learner_predict <- get(paste0('predict.', learner), envir = env)
  # Add an argument for "formula"
  f_formals <- c(alist(formula = ), formals(base_learner, envir = env))
  f_formals_predict <- c(formals(base_learner_predict, envir = env))
  # Use model.matrix formula *first*
  
  # Placeholder to pass CRAN checks
  object <- newdata <- X <- newX <- NULL
  
  f_learner <- function(formula, ...){
    args <- mget(ls())
    args$X <- model.frame(as.formula(formula), X)
    args$newX <- model.frame(as.formula(formula), newX)
    args$formula <- NULL
    out <- do.call("base_learner", args)
    out$fit$SL_formula <- formula
    class(out$fit) <- 'base_learner_f'
    return(out)
  }
  f_learner <- deparse(f_learner)
  f_learner <- eval(parse(text = paste(gsub(f_learner, pattern='base_learner', replacement = learner), collapse = '\n')))
  formals(f_learner) <- f_formals
  
  f_learner_predict <- function(...){
    args <- mget(ls())
    args$newdata <- model.frame(as.formula(object$SL_formula), newdata)
    args$formula <- NULL
    out <- do.call("predict.base_learner", args)
    return(out)
  }
  f_learner_predict <- deparse(f_learner_predict)
  f_learner_predict <- eval(parse(text = paste(gsub(f_learner_predict, pattern='base_learner', replacement = learner), collapse = '\n')))
  formals(f_learner_predict) <- f_formals_predict
  
  assign(x = paste0(learner, '_f'), value = f_learner, envir = env)
  assign(x = paste0('predict.', learner, '_f'), value = f_learner_predict, envir = env)
  return(paste0(learner, '_f'))
}
