
#' #' Estimate with vglmer for categorical outcomes
#' #' @importFrom mgcv interpret.gam
#' #' @rdname vglmer-multiclass
#' #' @export
#' expand_formula_clogit <- function(formula, choice_names, exempt = NULL){
#'   
#'   expand_fmla <- interpret.gam(formula, extra.special = NULL)
#'   expand_terms <- setdiff(expand_fmla$fake.names, exempt)
#'   add_terms <- c()
#'   for (i in expand_terms){
#' 
#'     if (grepl(i, pattern='v_s')){
#'       if (grepl(i, pattern='by ?=')){
#'         stop('multiclass with pooled not set up for v_s(x, by = g) yet...')
#'       }
#'       out <- sapply(choice_names, FUN=function(j){
#'         gsub(i, pattern='\\)$', replacement=paste0(', by = ', j, ')'))
#'       })
#'       add_terms <- c(add_terms, out)
#'     }else if (grepl(i, pattern='\\|')){
#'       out <- sapply(choice_names, FUN=function(j){
#'         paste0('(' ,
#'                gsub(i, pattern='^(.*\\|+)(.*)$', perl = T, replacement = paste0('\\1 ', paste0(j, ' :\\2'))),
#'                ')')
#'       })  
#'       add_terms <- c(add_terms, paste0('(', i, ')'), out)
#'     }else{
#'       out <- sapply(choice_names, FUN=function(j){
#'         paste0('(0 + ', i, ' | ', j, ')')
#'       })
#'       add_terms <- c(add_terms, i, out)
#'     }
#'     
#'   }
#'   add_terms <- c(exempt, paste0('(1 | ', choice_names, ')'), add_terms)
#'   formula <- update(formula, paste0('. ~ ', paste0(add_terms, collapse = ' + ')))
#'   return(formula)
#' }
#' 
#' #' @rdname vglmer-multiclass
#' #' @export
#' expand_data_clogit <- function(data, response, id_name){
#'   
#'   if ('choice' %in% names(data)){
#'     stop('"choice" cannot be column in data for expand_data_clogit')
#'   }
#'   if (id_name %in% names(data)){
#'     stop('"id_name" cannot be column in "data"')
#'   }
#' 
#'   data[[id_name]] <- 1:nrow(data)
#' 
#'   # Two Options: 
#'   # (a) Provide either a character vector, selects and expands levels
#'   # (b) Provide a matrix of choices
#'   
#'   if (is.matrix(response)){
#'     choices <- response
#'     choice_names <- colnames(choices)
#'     vec_choices <- apply(choices, MARGIN = 1, FUN = paste, collapse = ' & ')
#'     u_choices <- sort(unique(vec_choices))
#'     v_response <- NULL
#'   }else{
#'    
#'     mf_response <- data[,response]
#'     for (v in colnames(mf_response)){
#'       if (class(data[[v]]) != 'factor'){
#'         stop('All outcomes in cbind(x,y) ~ ... or x ~ ... must be factors.')
#'       }
#'       mf_response[,v] <- levels(data[[v]])[as.integer(mf_response[,v])]
#'     }
#'     
#'     if (is.matrix(mf_response) | is.data.frame(mf_response)){
#'       
#'       choices <- unique(na.omit(mf_response))
#'       rownames(choices) <- NULL
#'       L_bar <- prod(apply(mf_response, MARGIN = 2, FUN=function(i){length(unique(i))}))
#'       if (L_bar != nrow(choices)){
#'         warning('Number of observed combinations does not equal the maximal possible number', immediate. = TRUE)
#'       }
#'     }else{
#'       if (is.factor(mf_response)){
#'         choices <- levels(mf_response)
#'       }else{
#'         choices <- sort(unique(mf_response))
#'       }
#'       mf_response <- matrix(mf_response)
#'       choices <- matrix(choices)
#'       colnames(choices) <- 'choice'
#'     }
#'     if (any(grepl(as.vector(choices), pattern='&'))){
#'       stop('responses must not contain character "&"; this is used to identify unique combinations')
#'     }
#'     choice_names <- colnames(choices)
#'     vec_choices <- apply(choices, MARGIN = 1, FUN = paste, collapse = ' & ')
#'     u_choices <- sort(unique(vec_choices))
#' 
#'     if (length(choice_names) > 1){
#'       v_response <- apply(mf_response, MARGIN = 1, FUN = paste, collapse = ' & ')
#'     }else{
#'       v_response <- mf_response
#'     }
#'     
#'   }
#'   
#'   aug_data <- do.call('rbind', lapply(u_choices, FUN=function(i){
#'     copy_mf <- data
#'     copy_mf$response <- v_response
#'     copy_mf$choice <- i
#'     if (length(choice_names) > 1){
#'       for (j in choice_names){
#'         copy_mf[[j]] <- mf_response[,j]
#'       }
#'     }
#'     return(copy_mf)
#'   }))
#'   
#'   if (!is.null(v_response)){
#'     aug_data$pseudo_outcome <- as.numeric(aug_data$response == aug_data$choice)
#'   }
#'   attributes(aug_data)$id_name <- id_name
#'   attributes(aug_data)$choices <- choices
#'   attributes(aug_data)$unique_choices <- u_choices
#'   return(aug_data)
#' }
#' 
#' #' @rdname vglmer-multiclass
#' #' @export
#' vglmer_clogit <- function(
#'     formula,
#'     data, family, pooled = TRUE,
#'     id_name, choices, verify_formula = TRUE,
#'     control = vglmer_control()){
#'   
#'   family <- match.arg(family, c('poisson', 'binomial', 'multinomial'))
#'   
#'   if (family == 'multinomial' & pooled == FALSE){
#'     stop('family="multinomial" requires pooled=TRUE; 
#'       for separate Poisson regressions, use "poisson"')
#'   }
#'   
#'   if (missing(id_name)){
#'     id_name <- attr(data, 'id_name')
#'     if (is.null(id_name)){
#'       stop('id_name must be provided if data not created through "expand_data_clogit"')
#'     }
#'   }
#' 
#'   if (missing(choices)){
#'     choices <- attr(data, 'choices')
#'     if (is.null(choices)){
#'       stop('choices must be provided if data not created through "expand_data_clogit"')
#'     }
#'   }
#'   
#'   if (!pooled){
#'     if (verify_formula){
#'       parse_fmla <- vglmer_interpret.gam0(formula, extra.special = c('v_s', 'v_fe'))
#'       if (any(colnames(choices) %in% parse_fmla$pred.names)){
#'         stop('formula appears to contain variable in "choice"; this model cannot run with pooled=FALSE.
#'              If you are *sure* this is not an error, set verify_formula=FALSE')
#'       }
#'       warn_m <- 'If pooled=FALSE, the formula given to vglmer_clogit *must not* depend 
#'       on "choice" or have choice-constant covariates. 
#'       Carefully ensure this is the case or else the model will not converge. 
#'       verify_formula=TRUE does some attempts to address this but it should be checked carefully.'
#'       warning(warn_m, immediate. = TRUE)
#'     }
#'   }
#' 
#'   out <- vglmer_multiclass(
#'     formula = formula, data = data, 
#'     family = family, control = control,
#'     pooled = pooled,
#'     expand_formula = FALSE,
#'     clogit = TRUE,
#'     id_name = id_name, choices = choices)
#'   
#'   return(out)
#' }
#' 
#' #' @rdname vglmer-multiclass
#' #' @export
#' vglmer_OvR <- function(formula, data, family,
#'   control = vglmer_control(),
#'   pooled = TRUE){
#' 
#'   if (family == 'multinomial' & pooled == FALSE){
#'     stop('family="multinomial" requires pooled=TRUE; 
#'       for separate Poisson regressions, use "poisson"')
#'   }
#'   
#'   out <- vglmer_multiclass(
#'     formula = formula, data = data, 
#'     family = family, control = control,
#'     pooled = pooled)
#'   
#'   return(out)
#' }
#' 
#' vglmer_multiclass <- function(
#'     formula, data, family,
#'     control = vglmer_control(),
#'     pooled = TRUE, 
#'     expand_formula = TRUE,
#'     clogit = FALSE, verify_rank = FALSE,
#'     id_name, choices){
#'   
#'   parse_formula <- vglmer_interpret.gam0(
#'     subbars(formula), extra.special = c('v_s', 'v_fe'))
#'   
#'   if (any(!sapply(parse_formula$smooth.spec, inherits, what = 'vglmer_special'))){
#'     stop('gam specials are not permitted; use v_s(...) or v_fe(...) and see documentation.')
#'   }
#'   
#'   if (control$verify_columns){
#'     if (!all(parse_formula$pred.names %in% colnames(data))){
#'       missing_columns <- setdiff(parse_formula$pred.names, colnames(data))
#'       stop(
#'         paste0('The following columns are missing from "data". Can override with vglmer_control (not usually desirable): ', 
#'                paste(missing_columns, collapse =', '))
#'       )
#'     }
#'   }
#'   
#'   mf <- model.frame(parse_formula$fake.formula, data,
#'                       drop.unused.levels = TRUE)
#'   mf_response <- model.response(mf)
#'   for (v in colnames(mf_response)){
#'     if (class(data[[v]]) != 'factor'){
#'       stop('All outcomes in cbind(x,y) ~ ... or x ~ ... must be factors.')
#'     }
#'     mf_response[,v] <- levels(data[[v]])[as.integer(mf_response[,v])]
#'   }
#'   
#'   if (!clogit){
#'     if ('choice' %in% colnames(mf) | 'choice' %in% colnames(mf_response)){
#'       stop('"choice" cannot be a variable in the formula; it is used to denote the cartesian product of all options')
#'     }
#'     
#'     if (is.matrix(mf_response)){
#'       
#'       choices <- unique(na.omit(mf_response))
#'       rownames(choices) <- NULL
#'       L_bar <- prod(apply(mf_response, MARGIN = 2, FUN=function(i){length(unique(i))}))
#'       if (L_bar != nrow(choices)){
#'         warning('Number of observed combinations does not equal the maximal possible number')
#'       }
#'     }else{
#'       if (is.factor(mf_response)){
#'         choices <- levels(mf_response)
#'       }else{
#'         choices <- sort(unique(mf_response))
#'       }
#'       mf_response <- matrix(mf_response)
#'       choices <- matrix(choices)
#'       colnames(choices) <- 'choice'
#'     }
#'     
#'     if (any(grepl(as.vector(choices), pattern='&'))){
#'       stop('responses must not contain character "&"; this is used to identify unique combinations')
#'     }
#'     choice_names <- colnames(choices)
#'     vec_choices <- apply(choices, MARGIN = 1, FUN = paste, collapse = ' & ')
#'     u_choices <- sort(unique(vec_choices))
#'     
#'   }else{
#'     choice_names <- colnames(choices)
#'     vec_choices <- apply(choices, MARGIN = 1, FUN = paste, collapse = ' & ')
#'     u_choices <- sort(unique(vec_choices))
#'   }
#'   
#' 
#'   if (!clogit){
#'     message(paste0('Beginning Multiclass Estimation with ',
#'                    nrow(choices), ' Categories'))
#'   }else{
#'     message(paste0('Beginning Multiclass Estimation with ',
#'                    length(unique(data[,'choice'])), ' Categories'))
#'   }
#'   
#'   family <- match.arg(family, c('binomial', 'poisson', 'multinomial'))
#'   
#'   exempt_terms <- NULL
#'   
#'   if (family == 'multinomial'){
#'     fit_family <- 'poisson'
#'     message('Adding observation FE as family="multinomial" is chosen')
#'     if (pooled){
#'       if (clogit){
#'         formula <- update(formula, paste0('. ~ v_fe(', id_name, ') + .'))      
#'       }else{
#'         fe_term <- 'v_fe(id__)'
#'         formula <- update(formula, '. ~ v_fe(id__) + .')
#'         exempt_terms <- c(exempt_terms, fe_term)
#'       }
#'     }else{
#'       stop('family="multinomial" requires pooled=TRUE; 
#'       for separate Poisson regressions, use "poisson"')
#'     }
#'   }else{
#'     fit_family <- family
#'   }
#' 
#'   if (pooled){
#'     message('Using Single, Pooled, Model')  
#' 
#'     if (!clogit){
#' 
#'       if ('id__' %in% colnames(mf)){
#'         stop("'id__' cannot be in formula; this is reserved for observation identifier.")
#'       }
#'       
#'       # Expand the data to have N x bar{L} observations, i.e.
#'       # one observation for each observation-choice combination
#'       
#'       data <- expand_data_clogit(data = mf,
#'           response = mf_response,
#'           id_name = 'id__')  
#'       vec_response <- apply(mf_response, MARGIN = 1, FUN = paste, collapse = ' & ')
#'       data$pseudo_outcome <- as.numeric(vec_response == data$choice)
#'       checksum <- split(data$pseudo_outcome, data$id__)
#'       stopifnot(all(lengths(checksum) == length(u_choices)))
#'       stopifnot(all(sapply(checksum, sum) == 1))
#'       
#'       if (expand_formula){
#'         formula <- expand_formula_clogit(formula = formula, exempt = exempt_terms, choice_names = unique(c(choice_names, 'choice')))
#'       }
#'       
#'       message('Augmenting formula for pooled model; formula given to vglmer shown below')
#'       fmla <- update(formula, 'pseudo_outcome ~ .')
#'       message(deparse(fmla))
#'       
#'       fit_pooled <- vglmer(
#'         formula = fmla, data = data, family = fit_family,
#'         control = control)
#'       out <- list(fit = fit_pooled)
#'       
#'     }else{
#' 
#'       fit_pooled <- vglmer(
#'         formula = formula, data = data, family = fit_family,
#'         control = control)
#'       out <- list(fit = fit_pooled)
#'       
#'     }
#'   }else{
#'     
#'     message('Using Separate Models')
#' 
#'     if (clogit){
#' 
#'       fit_separate <- lapply(u_choices, FUN=function(ell){
#'         message(paste0('Fitting ', ell))
#'         fit_ell <- vglmer(
#'           formula = formula,
#'           data = data[data$choice == ell,,drop=FALSE], family = fit_family,
#'           control = control)
#'         return(fit_ell)
#'       })
#'       
#'     }else{
#'       
#'       vec_response <- apply(mf_response, MARGIN = 1, FUN = paste, collapse = ' & ')
#'       fit_separate <- lapply(u_choices, FUN=function(ell){
#'         message(paste0('Fitting ', ell))
#'         data$pseudo_outcome <- as.numeric(vec_response == ell)
#'         fmla <- update(formula, 'pseudo_outcome ~ .')
#'         fit_ell <- vglmer(formula = fmla, data = data, family = fit_family,
#'                           control = control)
#'         return(fit_ell)
#'       })
#'     }
#'     out <- list(fit = fit_separate)
#'   }
#'   
#'   if (!missing(id_name)){
#'     out$id_name <- id_name
#'   }
#'   out$clogit <- clogit
#'   out$choices <- choices
#'   out$choice_names <- choice_names
#'   out$unique_choices <- u_choices
#'   out$family <- family
#'   out$pooled <- pooled
#'   class(out) <- c('vglmer_multiclass', 'vglmer')
#'   return(out)
#' }
#' 
#' #' Multiclass Prediction
#' #' @export
#' predict.vglmer_multiclass <- function(
#'     object, newdata, id_name,
#'     ova_method = c('softmax', 'calibrated_softmax'),
#'     allow_missing_levels = FALSE, ...){
#'   
#'   if (length(list(...)) > 0) {
#'     stop("... not used for predict.vglmer_multiclass")
#'   }
#'   if (nrow(newdata) == 0){stop('newdata must not have zero rows.')}
#'   ova_method <- match.arg(ova_method, several.ok = FALSE)
#' 
#'   if (!object$clogit){
#'     
#'     id_name <- 'id__'
#'     newdata[[id_name]] <- 1:nrow(newdata)
#'     aug_newdata <- do.call('rbind', lapply(object$unique_choices, FUN=function(i){
#'       copy_mf <- newdata
#'       copy_mf$choice <- i
#'       if (length(object$choice_names) > 1){
#'         for (j in object$choice_names){
#'           copy_mf[[j]] <- object$choices[which(object$unique_choices == i),j]
#'         }
#'       }
#'       return(copy_mf)
#'     }))
#'     
#'   }else{
#'     
#'     if (missing(id_name)){
#'       stop('id_name must be provided for clogit prediction and data must be in "long" form.')
#'     }
#'     
#'     aug_newdata <- newdata
#'   }
#'   
#'   checksum <- split(aug_newdata$choice, aug_newdata[[id_name]])
#'   stopifnot(all(lengths(checksum) == nrow(object$choices)))
#' 
#'   if (object$pooled){
#'     
#'     pred_matrix <- predict(object$fit, 
#'       newdata = aug_newdata, skip_fe = TRUE,
#'       allow_missing_levels = allow_missing_levels)
#'     pred_matrix <- t(do.call('cbind', 
#'       split(pred_matrix, aug_newdata[[id_name]])))
#'     
#'   }else{
#'     
#'     if (object$clogit){
#'       
#'       pred_matrix <- mapply(object$fit, object$unique_choices, FUN=function(f, ch){
#'         predict(f, newdata = newdata[newdata$choice == ch,],
#'                 skip_fe = TRUE, allow_missing_levels = allow_missing_levels)
#'       })
#'       
#'     }else{
#'       pred_matrix <- sapply(
#'         object$fit, predict, 
#'         newdata = newdata, skip_fe = TRUE,
#'         allow_missing_levels = allow_missing_levels)
#'     }
#'   }
#'   
#'   if (object$family %in% c('poisson', 'multinomial')){
#'     pred_out <- FactorHet:::softmax_matrix(pred_matrix)
#'     colnames(pred_out) <- object$unique_choices
#'   }else if (object$family == 'binomial'){
#'     
#'     pred_out <- lapply(ova_method, FUN=function(m){
#'       if (m == 'softmax'){
#'         out <- FactorHet:::softmax_matrix(plogis(pred_matrix, log = TRUE))
#'       }else if (m == 'calibrated_softmax'){
#'         out <- FactorHet:::softmax_matrix(log(-plogis(-pred_matrix, log = TRUE)))
#'       }else{stop('...')}
#'       colnames(out) <- object$unique_choices
#'       return(out)
#'     })
#'     names(pred_out) <- ova_method
#'     if (length(ova_method) == 1){
#'       pred_out <- pred_out[[1]]
#'     }
#'   }else{
#'     stop('predict.vglmer_multiclass not set up for this family')
#'   }
#'   
#'   if (object$clogit){
#'     uid <- unique(aug_newdata[[id_name]])
#'     fmt_out <- split(
#'       pred_out,
#'       uid
#'     )
#'     fmt_out <- do.call('rbind', 
#'       mapply(fmt_out, names(fmt_out), SIMPLIFY = FALSE, FUN=function(i, id){
#'         data.frame(
#'           id = id,
#'           choice = colnames(pred_out),
#'           value = i,
#'           stringsAsFactors = FALSE
#'         )
#'     }))
#'     if (is.numeric(uid)){
#'       fmt_out$id <- as.numeric(fmt_out$id)
#'     }
#'     rownames(fmt_out) <- NULL
#'     colnames(fmt_out)[colnames(fmt_out) == 'id'] <- id_name
#'     return(fmt_out)
#'     
#'   }else{
#'     return(pred_out)
#'   }
#' }
