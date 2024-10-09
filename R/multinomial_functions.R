
multiclass_formula <- function(formula, choice_names, exempt = NULL){
  
  expand_fmla <- mgcv::interpret.gam(formula, extra.special = NULL)
  expand_terms <- setdiff(expand_fmla$fake.names, exempt)
  add_terms <- c()
  for (i in expand_terms){

    if (grepl(i, pattern='v_s')){
      if (grepl(i, pattern='by ?=')){
        stop('multiclass with pooled not set up for v_s(x, by = g) yet...')
      }
      out <- sapply(choice_names, FUN=function(j){
        gsub(i, pattern='\\)$', replacement=paste0(', by = ', j, ')'))
      })
      add_terms <- c(add_terms, out)
    }else if (grepl(i, pattern='\\|')){
      out <- sapply(choice_names, FUN=function(j){
        paste0('(' ,
               gsub(i, pattern='^(.*\\|+)(.*)$', perl = T, replacement = paste0('\\1 ', paste0(j, ' :\\2'))),
               ')')
      })  
      add_terms <- c(add_terms, paste0('(', i, ')'), out)
    }else{
      out <- sapply(choice_names, FUN=function(j){
        paste0('(0 + ', i, ' | ', j, ')')
      })
      add_terms <- c(add_terms, i, out)
    }
    
  }
  add_terms <- c(exempt, paste0('(1 | ', choice_names, ')'), add_terms)
  formula <- update(formula, paste0('. ~ ', paste0(add_terms, collapse = ' + ')))
  return(formula)
}
  
vglmer_multiclass <- function(
    formula, data,
    control = vglmer_control(),
    pooled = TRUE,
    family){
  
  interpret.gam <- mgcv:::interpret.gam
  
  parse_formula <- vglmer_interpret.gam0(subbars(formula),
                                         extra.special = c('v_s', 'v_fe'))
  
  if (any(!sapply(parse_formula$smooth.spec, inherits, what = 'vglmer_special'))){
    stop('gam specials are not permitted; use v_s(...) or v_fe(...) and see documentation.')
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
  
  mf <- model.frame(parse_formula$fake.formula, data,
                      drop.unused.levels = TRUE)
  mf_response <- model.response(mf)
  choices <- unique(mf_response)
  message(paste0('Beginning Multiclass Estimation with ',
                 length(choices), ' Categories'))
  
  family <- match.arg(family, c('binomial', 'poisson', 'multinomial'))
  
  exempt_terms <- NULL
  
  if (family == 'multinomial'){
    fit_family <- 'poisson'
    message('Adding observation FE as family="multinomial" is chosen')
    if (pooled){
      formula <- update(formula, '. ~ v_fe(id__) + .')
    }else{
      stop('family="multinomial" requires pooled=TRUE; for separate Poisson regressions, use "poisson"')
    }
    exempt_terms <- c(exempt_terms, 'v_fe(id__)')
  }else{
    fit_family <- family
  }
  
  if (pooled){
    message('Using Single, Pooled, Model')  
    if ('id__' %in% colnames(mf)){
      stop("'id__' cannot be in formula; this is reserved for observation identifier.")
    }
    mf[['id__']] <- 1:nrow(mf)
    aug_data <- do.call('rbind', lapply(choices, FUN=function(i){
      copy_mf <- mf
      copy_mf$response <- mf_response
      copy_mf$choice <- i
      return(copy_mf)
    }))
    aug_data$pseudo_outcome <- as.numeric(aug_data$response == aug_data$choice)
    
    formula <- multiclass_formula(formula = formula, exempt = exempt_terms, choice_names = 'choice')
    message('Augmenting formula for pooled model; formula given to vglmer shown below')
    fmla <- update(formula, 'pseudo_outcome ~ .')
    message(deparse(fmla))
    
    fit_pooled <- vglmer(
      formula = fmla, data = aug_data, family = fit_family,
           control = control)
    out <- list(fit = fit_pooled)
  }else{
    message('Using Separate Models')
    fit_separate <- lapply(choices, FUN=function(ell){
      message(paste0('Fitting ', ell))
      data$pseudo_outcome <- as.numeric(mf_response == ell)
      fmla <- update(formula, 'pseudo_outcome ~ .')
      fit_ell <- vglmer(formula = fmla, data = data, family = fit_family,
             control = control)
      return(fit_ell)
    })
    out <- list(fit = fit_separate)
  }
  
  out$choices <- choices
  out$family <- family
  out$pooled <- pooled
  class(out) <- c('vglmer_multiclass', 'vglmer')
  return(out)
}

#' Multiclass Prediction
#' @export
predict.vglmer_multiclass <- function(
    object, newdata,
    ova_method = c('softmax', 'calibrated_softmax'),
    allow_missing_levels = FALSE){
  
  ova_method <- match.arg(ova_method, several.ok = TRUE)
  
  if (object$pooled){
    
    newdata[['id__']] <- 1:nrow(newdata)
    aug_newdata <- do.call('rbind', lapply(object$choices, FUN=function(i){
      copy_mf <- newdata
      copy_mf$choice <- i
      return(copy_mf)
    }))
    pred_matrix <- predict(object$fit,
      newdata = aug_newdata, 
      allow_missing_levels = allow_missing_levels)
    pred_matrix <- t(do.call('cbind', split(pred_matrix, aug_newdata$id__)))
    
  }else{
    pred_matrix <- sapply(object$fit, predict, 
                          newdata = newdata, 
                          allow_missing_levels = allow_missing_levels)
  }
  
  if (object$family %in% c('poisson', 'multinomial')){
    pred_out <- FactorHet:::softmax_matrix(pred_matrix)
  }else if (object$family == 'binomial'){
    
    pred_out <- lapply(ova_method, FUN=function(m){
      if (m == 'softmax'){
        out <- FactorHet:::softmax_matrix(plogis(pred_matrix, log = TRUE))
      }else if (m == 'calibrated_softmax'){
        out <- FactorHet:::softmax_matrix(log(-plogis(-pred_matrix, log = TRUE)))
      }else{stop('...')}
      return(out)
    })
    names(pred_out) <- ova_method
    if (length(ova_method) == 1){
      pred_out <- pred_out[[1]]
    }
  }else{
    stop('predict.vglmer_multiclass not set up for this family')
  }
  return(pred_out)
}
