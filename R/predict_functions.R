#' Predict after vglmer
#' 
#' Get linear predictor for new observations after using vglmer.
#' 
#' @name vglmer_predict
#' @param object Object from vglmer.
#' @param newdata Data to get predictions on.
#' @param samples How many samples to draw? 0, default, gets the expected value.
#'   Two methods: 
#'   \itemize{
#'   \item Number - draw that many samples.
#'   \item Matrix - use previous samples for prediction. Currently used for MAVB
#'   but will be depreciated in updated version.
#'   }
#' @param samples_only Return only samples *not* linear predictor.
#' @param summary Return summary of linear predictor, not full posterior.
#' @param allow_missing_levels Allow prediction for random effects not in object.
#'   As is standard, give an estimate of "0" for that effect.
#' @param ... Not used; included to maintain compatability
#' 
#' @examples
#' 
#' sim_data <- data.frame(x = rnorm(100), 
#' y = rbinom(100, 1, 0.5), 
#' g =sample(letters, 100, replace = TRUE))
#' 
#' # Run with defaults
#' est_vglmer <- vglmer(y ~ x + (x|g), data = sim_data, family = "binomial")
#' 
#' # Simple prediction 
#' predict(est_vglmer, newdata = sim_data)
#' 
#' \dontrun{
#' #Fails!
#' predict(est_vglmer, newdata = data.frame(g = "AB", x = 0))
#' }
#' predict(est_vglmer, newdata = data.frame(g = "AB", x = 0), 
#' allow_missing_levels = TRUE)
#' @importFrom stats delete.response terms
#' @export
predict.vglmer <- function(object, newdata, 
                           samples = 0, samples_only = FALSE, 
                           summary = TRUE, allow_missing_levels = FALSE, ...){
  if (length(list(...)) > 0){
    stop('... not used for predict.vglmer')
  }
  
  rownames(newdata) <- as.character(1:nrow(newdata))
  
  fmla <- object$formula
  #Extract X (FE design matrix)
  X <- model.matrix(delete.response(terms(nobars(fmla))), data = newdata)
  
  orig_X_names <- rownames(object$beta$mean)
  if (!identical(colnames(X), orig_X_names)){
    print(all.equal(colnames(X), orig_X_names))
    stop('Misaligned Fixed Effects')
  }
  
  #Extract the Z (Random Effect) design matrix.
  mk_Z <- model.frame(delete.response(terms(subbars(fmla))), data = newdata)
  rownames_Z <- rownames(mk_Z)
  mk_Z <- mkReTrms(findbars(fmla), mk_Z, reorder.terms = FALSE, reorder.vars = FALSE)
  Z <- t(mk_Z$Zt)

  #RE names and names of variables included for each.
  names_of_RE <- mk_Z$cnms
  
  number_of_RE <- length(mk_Z$Gp) - 1
  #The position that demarcates each random effect.
  #That is, breaks_for_RE[2] means at that position + 1 does RE2 start.
  breaks_for_RE <- c(0, cumsum(diff(mk_Z$Gp)))
  #Dimensionality of \alpha_{j,g}, i.e. 1 if random intercept
  #2 if random intercept + random slope
  d_j <- lengths(names_of_RE)
  #Number of GROUPs for each random effect.
  g_j <- diff(mk_Z$Gp)/d_j
  
  #Empty vector to build the formatted names for each random effect.
  fmt_names_Z <- c()
  init_Z_names <- colnames(Z)
  for (v in 1:number_of_RE){
    name_of_effects_v <- names_of_RE[[v]]
    
    mod_name <- rep(name_of_effects_v, g_j[v])
    
    levels_of_re <- init_Z_names[(1+breaks_for_RE[v]):breaks_for_RE[v+1]]
    
    fmt_names_Z <- c(fmt_names_Z, paste0(names(names_of_RE)[v], ' @ ', mod_name, ' @ ', levels_of_re))
    
  }
  colnames(Z) <- fmt_names_Z
  #####
  ###Confirm Alignment of the Z
  #####
  orig_Z_names <- rownames(object$alpha$mean)
  
  not_in_original_Z <- setdiff(fmt_names_Z, orig_Z_names)
  not_in_new_Z <- setdiff(orig_Z_names, fmt_names_Z)
  
  if (length(not_in_original_Z) > 0){
    if (!allow_missing_levels){stop('New levels not allowed unless allow_missing_levels = TRUE')}
  }
  
  in_both <- intersect(fmt_names_Z, orig_Z_names)
  recons_Z <- drop0(sparseMatrix(i = 1, j = 1, x = 0, dims = c(nrow(Z), length(orig_Z_names))))
  colnames(recons_Z) <- orig_Z_names
  rownames(recons_Z) <- rownames_Z
  
  recons_Z[,match(in_both, orig_Z_names)] <- Z[,match(in_both, fmt_names_Z)] 
  
  #Check that the entirely missing columns match those not in the original
  checksum_align <- setdiff(not_in_new_Z, sort(names(which(colSums(recons_Z != 0) == 0))))
  if (length(checksum_align) > 0){
    stop('Alignment Error')
  }
  
  Z <- recons_Z
  rm(recons_Z)
  ####
  total_obs <- rownames(newdata)
  obs_in_both <- intersect(rownames(X), rownames(Z))
  
  XZ <- cbind(X[match(obs_in_both, rownames(X)),,drop=F],
              Z[match(obs_in_both, rownames(Z)),,drop=F])
  
  factorization_method <- object$control$factorization_method
  if (is.matrix(samples)){
    if (ncol(samples) != ncol(XZ)){
      stop('Samples must be {m, ncol(Z) + ncol(X)}')
    }
    samples <- t(samples)
    only.lp <- FALSE
    print('Using Provided Samples')
  }else{
    if (samples == 0){
      only.lp <- TRUE
    }else{
      only.lp <- FALSE
    }
    if (factorization_method %in% c('strong', 'partial')){
      
      vi_alpha_mean <- object$alpha$mean
      vi_alpha_decomp <- object$alpha$decomp_var
      
      p.Z <- nrow(vi_alpha_mean)
      
      vi_beta_mean <- object$beta$mean
      vi_beta_decomp <- object$beta$decomp_var
      
      p.X <- nrow(vi_beta_mean)
      
      if (!only.lp){
        sim_init_alpha <- matrix(rnorm(samples * p.Z), ncol = samples)
        sim_init_alpha <- t(vi_alpha_decomp) %*% sim_init_alpha
        sim_init_alpha <- sim_init_alpha + kronecker(vi_alpha_mean, t(matrix(1, samples)))
        
        sim_init_beta <- matrix(rnorm(samples * p.X), ncol = samples)
        sim_init_beta <- t(vi_beta_decomp) %*% sim_init_beta
        sim_init_beta <- sim_init_beta + kronecker(vi_beta_mean, t(matrix(1,samples)))
      }else{
        sim_init_alpha<- vi_alpha_mean
        sim_init_beta <- vi_beta_mean
      }
      
    }else if (factorization_method == 'weak'){
      
      vi_alpha_mean <- object$alpha$mean
      p.Z <- nrow(vi_alpha_mean)
      
      vi_beta_mean <- object$beta$mean
      p.X <- nrow(vi_beta_mean)
      
      if (!only.lp){
        vi_joint_decomp <- object$joint
        sim_init_joint <- matrix(rnorm(samples * (p.X + p.Z)), ncol = samples)
        sim_init_joint <- t(vi_joint_decomp) %*% sim_init_joint
        
        sim_init_beta <- sim_init_joint[1:p.X, , drop = F]
        sim_init_alpha <- sim_init_joint[-1:-p.X, ,drop = F]
        
        rm(sim_init_joint)
        
        sim_init_alpha <- sim_init_alpha + kronecker(vi_alpha_mean, t(matrix(1, samples)))
        sim_init_beta <- sim_init_beta + kronecker(vi_beta_mean, t(matrix(1,samples)))
      }else{
        sim_init_alpha <- vi_alpha_mean
        sim_init_beta <- vi_beta_mean
      }
      
    }else{stop('')}
    
    samples <- rbind(sim_init_beta, sim_init_alpha)
    rm(sim_init_beta, sim_init_alpha)
  }
  
  if (samples_only){
    return(samples)
  }
  
  lp <- XZ %*% samples
  if (summary){
    if (!only.lp){
      lp <- t(apply(lp, MARGIN = 1, FUN=function(i){c(mean(i), var(i))}))
      lp <- data.frame(mean = lp[,1], var = lp[,2])  
      lp <- lp[match(total_obs, obs_in_both),]
    }else{
      lp <- as.vector(t(apply(lp, MARGIN = 1, FUN=function(i){mean(i)})))
      lp <- lp[match(total_obs, obs_in_both)]
    }
    return(lp)
  }else{
    lp <- lp[match(total_obs, obs_in_both)]
    return(lp)
  }
}

