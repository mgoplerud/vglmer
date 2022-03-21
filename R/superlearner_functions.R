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
#' @export
SL.glmer <- function(Y, X, newX, formula, family, id, obsWeights, control = glmerControl()) {
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
  
  fit.glmer <- lme4::glmer(formula, data = X, weights = obsWeights, family = family, control = control)
  pred <- stats::predict(fit.glmer, newdata = newX, allow.new.levels = TRUE, type = 'response')
  fit <- list(object = fit.glmer)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.glmer")
  return(out)
}

#' @export
predict.SL.glmer <- function(object, newdata, allow_missing_levels = TRUE, ...){
  if(!requireNamespace('lme4', quietly = FALSE)) {stop("SL.glmer requires the lme4 package, but it isn't available")} 
  
  pred <- predict(object$object, newdata = newdata, allow.new.levels = TRUE, type = 'response')
  return(pred)
}

add_formula_SL <- function(learner, env = parent.frame()){

  base_learner <- get(learner, envir = env)
  base_learner_predict <- get(paste0('predict.', learner), envir = env)
  # Add an argument for "formula"
  f_formals <- c(alist(formula = ), formals(base_learner, envir = env))
  f_formals_predict <- c(formals(base_learner_predict, envir = env))
  # Use model.matrix formula *first*
  
  # Placeholder to pass CRAN
  newdata <- X <- newX <- NULL
  
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
