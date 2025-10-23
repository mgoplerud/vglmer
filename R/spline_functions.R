
#' Code from Wand and Ormerod (2008)
#' Found here here: 10.1111/j.1467-842X.2008.00507.x
#' @param a lower boundary
#' @param b upper boundary
#' @param intKnots internal knots
#' @keywords internal
formOmega <- function(a,b,intKnots){
  allKnots <- c(rep(a,4),intKnots,rep(b,4))
  K <- length(intKnots) ; L <- 3 * (K+8)
  xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+
               rep(allKnots,each=3)[-c(1,2,L)])/2
  wts <- rep(diff(allKnots),each=3) * rep(c(1,4,1)/6,K+7)
  Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                    outer.ok=TRUE)$design
  Omega <- drop0(t(Bdd) %*% Diagonal(x = wts) %*% Bdd)
  return(Omega)
}

#' Create splines for use in vglmer
#' 
#' This function estimates splines in \code{vglmer}, similar to \code{s(...)} in
#' \code{mgcv} albeit with many fewer options than \code{mgcv}. It allows for
#' truncated (linear) splines or O'Sullivan splines. Please see \link{vglmer}
#' for more discussion and examples.
#' 
#' @param ... Variable name, e.g. \code{v_s(x)}
#' @param type Default (\code{"tpf"}) uses truncated linear splines for the
#'   basis. The other option (\code{"o"}) uses O'Sullivan splines (Wand and
#'   Ormerod 2008).
#' @param knots Default (\code{NULL}) uses \eqn{K=min(N/4,35)} knots evenly
#'   spaced at quantiles of the covariate \code{x}. A single number specifies a
#'   specific number of knots; a vector can set custom locations for knots.
#' @param by A categorical or factor covariate to interact the spline with; for
#'   example, \code{v_s(x, by = g)}.
#' @param by_re Default (\code{TRUE}) regularizes the interactions between the
#'   categorical factor and the covariate. See "Details" in \link{vglmer} for
#'   more discussion.
#' @param force_vector Force that argument to \code{knots} is treated as vector.
#'   This is usually not needed unless \code{knots} is a single integer that
#'   should be treated as a single knot (vs. the number of knots).
#' @param outer_okay Default (\code{FALSE}) does not permit values in \code{x}
#'   to exceed the outer knots.
#' @importFrom splines bs
#' 
#' @return This function returns a list of class of \code{vglmer_spline} that is
#'   passed to unexported functions. It contains the arguments noted above where
#'   \code{...} is parsed into an argument called \code{term}.
#'   
#' @references 
#' Wand, Matt P. and Ormerod, John T. 2008. "On Semiparametric Regression with
#' O'Sullivan Penalized Splines". \emph{Australian & New Zealand Journal of
#' Statistics}. 50(2): 179-198.
#' 
#' Wood, Simon N. 2017. \emph{Generalized Additive Models: An Introduction with
#' R}. Chapman and Hall/CRC.
#' @export
v_s <- function(..., type = 'tpf', knots = NULL, by = NA,
                by_re = TRUE, force_vector = FALSE,
                outer_okay = FALSE){
  if (!(type %in% c('tpf', 'o', 'gKRLS', 'randwalk'))){stop('non tpf not set up yet...')}
  # Using mgcv's syntax for "s" to make it work with "interpret.gam"
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  if (type != 'gKRLS'){
    if (d > 1){stop('Unlike mgcv, only provide a single variable')}
  }
  by.var <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  if (by.var == "."){
    stop("by=. not allowed")
  }
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  if (term[1] == "."){
    stop("s(.) not supported.")
  }
  
  term[1] <- attr(terms(reformulate(term[1])), "term.labels")

  label <- paste0("v_s(", term[1], ")")
  
  ret <- list(term = term, outer_okay = outer_okay, force_vector = force_vector,
              by = by.var, type = type, knots = knots, mi = FALSE,
              by_re = by_re)
  class(ret) <- 'vglmer_spline'
  
  return(ret)
}

# v_mi <- function(..., rank, xt = NULL){
# 
#   vars <- as.list(substitute(list(...)))[-1]
#   if (length(vars) != 2){stop('must provide two grouping factors')}
#   
#   term <- sapply(vars, FUN=function(i){deparse(i, backtick = TRUE, width.cutoff = 500)})
#   term <- sapply(term, FUN=function(i){attr(terms(reformulate(i)), "term.labels")})
#   
#   label <- paste0("v_mi(", paste(term, collapse=', '), ")")
# 
#   ret <- list(
#     mi = TRUE,
#     hier = FALSE,
#     term = term,
#     rank = rank,
#     by = "NA",
#     xt = xt
#   )
#   class(ret) <- 'vglmer_multiplicative'
#   return(ret)
# }

mi <- function(..., rank, xt = NULL){
  
  vars <- as.list(substitute(list(...)))[-1]
  if (length(vars) != 2){stop('must provide two grouping factors')}
  
  term <- lapply(vars, FUN=function(i){deparse(i, backtick = TRUE, width.cutoff = 500)})
  is_formula <- sapply(term, FUN=function(i){
    i <- tryCatch(as.formula(i), error = function(e){NULL})
    if (inherits(i, 'formula')){
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  
  term <- mapply(term, is_formula, SIMPLIFY = FALSE,
                 FUN=function(i, is_f){
                   if (is_f){
                     all.vars(formula(i))
                   }else{
                     attr(terms(reformulate(i)), "term.labels")
                   }
                  })
  if (length(intersect(term[[1]], term[[2]])) > 0){
    stop('mi does not allow overlapping terms between "u" and "v".')
  }
  label <- paste0("mi(", paste(sapply(term, FUN=function(i){i[1]}), collapse=', '), ")")
  ret <- list(
    mi = TRUE,
    hier = TRUE,
    hier_fmla = term,
    is_formula = is_formula,
    raw_formula = vars,
    term = unlist(term),
    rank = rank,
    by = "NA",
    xt = xt
  )
  class(ret) <- 'vglmer_multiplicative'
  return(ret)
}

#' @importFrom mgcv smooth.construct
#' @importFrom splines spline.des
vglmer_build_spline <- function(object, data, add_unpenalized = FALSE){

  if (inherits(object, 'vglmer_spline')){
    spline_type <- 'vglmer'
    x <- data[object$term]
    by <- data[[object$by]]
    knots <- object$knots
    Boundary.knots <- object$Boundary.knots
    type <- object$type 
    force_vector <- object$force_vector
    outer_okay <- object$outer_okay
    override_warn <- object$override_warn
    if (is.null(override_warn)){
      override_warn <- FALSE
    }
    by_re <- object$by_re
    xt <- object$xt
    # x, knots = NULL, Boundary.knots = NULL, 
    # by, type, override_warn = FALSE, 
    # outer_okay = FALSE, by_re = NULL, force_vector = FALSE
    
    if (type %in% c('gKRLS', 'randwalk')){
      
      if (type == 'gKRLS'){
        
        if (!is.null(knots)){
          object_mgcv <- knots
        }else{
          object_mgcv <- list(
            term = colnames(x),
            xt = xt,
            p.order = NA,
            bs.dim = -1,
            fixed = FALSE,
            by = 'NA'
          )
          class(object_mgcv) <- 'gKRLS.smooth.spec'
          object_mgcv <- smooth.construct(object_mgcv, data = x)
        }
        
      }else{
        
        if (!is.null(knots)){
          object_mgcv <- knots
        }else{
          if (!is.null(xt)){stop('xt must be null for randwalk at the moment...')}
          object_mgcv <- list(
            term = colnames(x),
            xt = xt,
            p.order = NA,
            bs.dim = -1,
            fixed = FALSE,
            by = 'NA'
          )
          class(object_mgcv) <- 'randwalk.smooth.spec'
          object_mgcv <- smooth.construct(object_mgcv, data = x, knots = NULL)
        }
      }
      
      x <- Predict.matrix(object_mgcv, data = x)
      # if (outer_okay){
      #   object_mgcv$internal_override <- TRUE
      #   x <- Predict.matrix(object_mgcv, data = x)
      # }else{
      #   object_mgcv$internal_override <- FALSE
      #   x <- Predict.matrix(object_mgcv, data = x)
      # }
      # object_mgcv$internal_override <- NULL
      colnames(x) <- paste0('base @ ', 1:ncol(x))
      object_mgcv$spline_type <- spline_type
      
      spline_attr <- list(knots=object_mgcv)
      
      if (!is.null(by)){
        base_x <- x
        u_by <- sort(unique(by))
        
        if (!outer_okay){
          x_by <- sparseMatrix(i = 1:length(by), j = match(by, u_by), x = 1)
        }else{
          match_j <- match(by, u_by)
          match_i <- 1:length(by)
          
          match_i <- match_i[!is.na(match_j)]
          match_j <- match_j[!is.na(match_j)]
          x_by <- sparseMatrix(i = match_i, j= match_j, x = 1, dims = c(length(by), length(u_by)))
        }
        
        
        names_x <- as.vector(outer(1:ncol(x), u_by, FUN=function(x,y){paste(y,x, sep = ' @ ')}))
        x <- t(KhatriRao(t(x_by), t(x)))
        colnames(x) <- names_x  
        
        colnames(base_x) <- paste0('base @ ', 1:ncol(base_x))
        
        out <- list(x = x, attr = spline_attr)
        class(out) <- c('spline_sparse')
        
        base_out <- list(x = base_x, attr = spline_attr)
        class(base_out) <- c('spline_sparse')
        
        return(
          list(base_out, out)
        )
      }else{
        out <- list(x = x, attr = spline_attr)
        class(out) <- c('spline_sparse')
        return(list(out))
      }
    }else{
      x <- x[object$term]
      if (ncol(x) > 1){stop('x should be only column for vglmer-based splines')}
      x <- x[[1]]
      if (is.null(knots)){
        ux <- length(unique(x))
        if (ux < 4){stop('Cannot fit spline with fewer than 4 unique values.')}
        # Use the knot heuristic in Ruppert by default.
        # Keeps the size of the problem feasible.
        numIntKnots <- floor(c(min(ux/4, 35)))
        
        intKnots <- quantile(unique(x),
                             seq(0,1,length=(numIntKnots+2)
                             )[-c(1,(numIntKnots+2))])
        names(intKnots) <- NULL
      }else if (length(knots) == 1 & !force_vector){
        
        if (knots < 1){
          stop('If an integer, at least one knot must be provided. force_vector=TRUE may be useful here.')
        }
        if (as.integer(knots) != knots){
          warning('knots appears to be not be an integer. Using "as.integer"')
          knots <- as.integer(knots)
          message(paste0('knots argument turned into ', knots, ' by coercion.'))
        }
        
        numIntKnots <- knots
        
        intKnots <- quantile(unique(x),seq(0,1,length=
                                             (numIntKnots+2))[-c(1,(numIntKnots+2))])
        names(intKnots) <- NULL
      }else{
        # Sort user provided knots
        knots <- sort(knots)
        
        # Is any knot big above the maximum in the data?
        cond_1 <- any(knots >= max(x, na.rm=T))
        # Is any knot below the minimum in the data?
        cond_2 <- any(knots <= min(x, na.rm=T))
        # If so, issue warning
        if (!cond_1 | !cond_2){
          if (!override_warn){
            warning('observed data is outside of the self-provided knots')
          }
        }
        intKnots <- knots
      }
      
      if (is.null(Boundary.knots)){
        Boundary.knots <- range(x, na.rm=T) 
      }else{
        stopifnot(length(Boundary.knots) == 2)
      }
      
      if (type == 'tpf'){
        aug_knots <- c(Boundary.knots[1], intKnots, Boundary.knots[2])
        
        x <- outer(x, aug_knots[-c(1,length(aug_knots))], '-')
        x <- drop0(x * (x > 0))
        spline_attr <- list(D = Diagonal(n = ncol(x)), Boundary.knots = Boundary.knots,
                            knots = intKnots)
        
      }else if (type == 'o'){
        
        # Form Omega from Wand and Ormerod (2008)
        D <- formOmega(a = Boundary.knots[1], b = Boundary.knots[2], intKnots = intKnots)
        # eigen decompose
        eD <- eigen(D)
        # transform spline design
        if (override_warn){
          wrapper_bs <- function(x){suppressWarnings(x)}
        }else{
          wrapper_bs <- function(x){x}
        }
        x <- wrapper_bs(splines::bs(x = x, knots = intKnots, 
                                    degree = 3, intercept = TRUE,
                                    Boundary.knots = Boundary.knots))
        
        x <- x %*% eD$vectors[,seq_len(ncol(D)-2)] %*% 
          Diagonal(x = 1/sqrt(eD$values[seq_len(ncol(D) - 2)]))
        
        spline_attr <- list(D = Diagonal(n = ncol(x)), 
                            Boundary.knots = Boundary.knots,
                            knots = intKnots, eigen_D = eD)
        
      }else{
        stop('splines only set up for tpf and o')
      }
      
      spline_attr$spline_type <- spline_type
      spline_attr$by_re <- by_re
      
      if (!is.null(by)){
        
        base_x <- x
        u_by <- sort(unique(by))
        
        if (!outer_okay){
          x_by <- sparseMatrix(i = 1:length(by), j = match(by, u_by), x = 1)
        }else{
          match_j <- match(by, u_by)
          match_i <- 1:length(by)
          
          match_i <- match_i[!is.na(match_j)]
          match_j <- match_j[!is.na(match_j)]
          x_by <- sparseMatrix(i = match_i, j= match_j, x = 1, dims = c(length(by), length(u_by)))
        }
        
        
        names_x <- as.vector(outer(1:ncol(x), u_by, FUN=function(x,y){paste(y,x, sep = ' @ ')}))
        x <- t(KhatriRao(t(x_by), t(x)))
        colnames(x) <- names_x  
        
        colnames(base_x) <- paste0('base @ ', 1:ncol(base_x))
        
        out <- list(x = x, attr = spline_attr)
        class(out) <- c('spline_sparse')
        
        base_out <- list(x = base_x, attr = spline_attr)
        class(base_out) <- c('spline_sparse')
        
        return(
          list(base_out, out)
        )
      }else{
        colnames(x) <- paste0('base @ ', 1:ncol(x))
        out <- list(x = x, attr = spline_attr)
        class(out) <- c('spline_sparse')
        return(list(out))
      }
      
    }
  }else{
    
    spline_type <- 'mgcv'
    # Use "mgcv" function to construct spline
    parse_spline <- smooth.construct(object = object, data = data[object$term], knots = NULL)    
    by <- parse_spline$by
    if (by == "NA"){
      by <- NULL
    }else{
      by <- data[[parse_spline$by]]
    }
    
    if (length(parse_spline$S) > 1){
      stop('Only one "S" allowed...')
    }else if (length(parse_spline$S) == 0){
      x <- as.matrix(parse_spline$X[,])
      eD <- list(values = rep(0, ncol(x)), vectors = Diagonal(n=ncol(x)))
    }else{
      D <- parse_spline$S[[1]]
      nat_param <- mgcv:::nat.param(X = parse_spline$X, S = parse_spline$S[[1]], type = 1)
      parse_spline$nat_param <- nat_param
      D <- Diagonal(x=c(nat_param$D, rep(0, ncol(D) - nat_param$rank)))
      x <- drop0(nat_param$X)
      eD <- eigen(D)
      eD$vectors <- drop0(eD$vectors)
    }
    zero_ev <- which(eD$values == 0 | eD$values < max(eD$values) * sqrt(.Machine$double.eps))
    nonzero_ev <- setdiff(1:ncol(x), zero_ev)
    # Demean and whiten and thus DROP intercept...
    unpen_x <- x %*% eD$vectors[,zero_ev,drop=FALSE]
    if (ncol(unpen_x) > 0){
      mean_unpen_x <- apply(unpen_x, MARGIN = 2, mean)
      unpen_x <- sweep(unpen_x, MARGIN = 2, STATS = mean_unpen_x, FUN = '-')
      cov_unpen_x <- cov(as.matrix(unpen_x))
      eigen_cov <- eigen(cov_unpen_x)
      root <- ifelse(eigen_cov$values < max(eigen_cov$values) * sqrt(.Machine$double.eps), 0, sqrt(1/eigen_cov$values))
      whiten_unpen_x <- eigen_cov$vectors %*% Diagonal(x=root)
      unpen_x <- unpen_x %*% whiten_unpen_x
      zero_unpen <- sqrt(colSums(unpen_x^2))
      zero_unpen <- (is.na(zero_unpen) | (zero_unpen < 1e-6))
      unpen_x <- unpen_x[,!zero_unpen,drop=F]
    }else{
      mean_unpen_x <- whiten_unpen_x <- zero_unpen <- NULL
    }

    x <- x %*% eD$vectors[,nonzero_ev] %*% 
      Diagonal(x = 1/sqrt(eD$values[nonzero_ev]))
    
    spline_attr <- list(D = Diagonal(n = ncol(x)), unpen_x = unpen_x,
                        std_unpen = list(mean = mean_unpen_x, whiten = whiten_unpen_x, zero = zero_unpen),
                        mgcv_object = parse_spline, eigen_D = eD,
                        spline_type = spline_type)
    
    if (ncol(unpen_x) == 0){
      spline_attr$unpen_x <- NULL
      spline_attr$std_unpen <- NULL
    }
    
    by_re <- TRUE
    outer_okay <- FALSE
    spline_attr$by_re <- by_re
    
    if (!is.null(by)){
      
      base_x <- x
      u_by <- sort(unique(by))
      
      if (!outer_okay){
        x_by <- sparseMatrix(i = 1:length(by), j = match(by, u_by), x = 1)
      }else{
        match_j <- match(by, u_by)
        match_i <- 1:length(by)
        
        match_i <- match_i[!is.na(match_j)]
        match_j <- match_j[!is.na(match_j)]
        x_by <- sparseMatrix(i = match_i, j= match_j, x = 1, dims = c(length(by), length(u_by)))
      }
      
      
      names_x <- as.vector(outer(1:ncol(x), u_by, FUN=function(x,y){paste(y,x, sep = ' @ ')}))
      x <- t(KhatriRao(t(x_by), t(x)))
      colnames(x) <- names_x  
      if (add_unpenalized){
        browser()
      }
      
      colnames(base_x) <- paste0('base @ ', 1:ncol(base_x))
      
      out <- list(x = x, attr = spline_attr)
      class(out) <- c('spline_sparse')
      
      base_out <- list(x = base_x, attr = spline_attr)
      class(base_out) <- c('spline_sparse')
      return(list(base_out, out))
    }else{
      if (ncol(x) > 0){
        colnames(x) <- paste0('base @ ', 1:ncol(x))
      }
      out <- list(x = x, attr = spline_attr)
      class(out) <- c('spline_sparse')
      return(list(out))
    }
  }

}

vglmer_build_mi <- function(object, data){

  knots <- object$knots
  x <- data[object$term]
  rank <- object$rank
  hier <- object$hier
  raw_formula <- object$raw_formula
  is_formula <- object$is_formula
  hier_fmla <- object$hier_fmla
  
  if (hier){
    if (length(hier_fmla) != 2){
      stop('must have two factors...')
    }
    nested_hier <- lapply(hier_fmla, FUN=function(i){unique(x[,i,drop=F])})
  }else{
    nested_hier <- NULL
    hier_fmla <- as.list(colnames(x))
    if (ncol(x) != 2){stop('must have two factors...')}
  }
  
  out <- mapply(raw_formula, is_formula, SIMPLIFY = FALSE, FUN=function(i,f){
    if (f){
      
      interpret_i <- mgcv:::interpret.gam0(i, extra.special = 'v_s')

      if (length(interpret_i$smooth.spec) > 0 & !is.null(knots)){
        parse_i <- sapply(interpret_i$smooth.spec, FUN=function(i){paste(i$term, collapse=',')})
        parse_i <- mapply(interpret_i$smooth.spec, parse_i, 
          SIMPLIFY = FALSE, FUN=function(i,n_i){
            if (knots$predict_attr[[n_i]]$spline_type == 'mgcv'){
              attr_i <- knots$predict_attr[[n_i]]
              x <- mgcv::Predict.matrix(attr_i$mgcv_object, data = data)
              if (!is.null(attr_i$mgcv_object$nat_param)){
                x <- x %*% attr_i$mgcv_object$nat_param$P
              }
              eD <- attr_i$eigen_D
              zero_ev <- which(eD$values == 0 | eD$values < max(eD$values) * sqrt(.Machine$double.eps))
              nonzero_ev <- setdiff(1:ncol(x), zero_ev)
              if (!is.null(attr_i$std_unpen)){
                unpen_x <- x %*% eD$vectors[,zero_ev,drop=FALSE]
                unpen_x <- sweep(unpen_x, MARGIN = 2, STATS = attr_i$std_unpen$mean, FUN = '-')
                unpen_x <- unpen_x %*% attr_i$std_unpen$whiten
                unpen_x <- unpen_x[,!attr_i$std$zero,drop=FALSE]
              }else{
                unpen_x <- NULL
              }
              x <- x %*% eD$vectors[,nonzero_ev] %*% 
                Diagonal(x = 1/sqrt(eD$values[nonzero_ev]))

            }else{
              
              object <- i
              for (v in names(knots$predict_attr)){
                object[[v]] <- knots$predict_attr[[v]]
              }
              x <- vglmer_build_spline(
                object = object, data = data)
              if (length(x) != 1){stop('...')}
              x <- x[[1]]
              x <- x$x
              unpen_x <- NULL
            }
            
            return(list(list(x=x, attr = list(unpen_x = unpen_x))))
        })
      }else{
        parse_i <- lapply(interpret_i$smooth.spec, 
                          vglmer_build_spline, 
                          data = data,
                          add_unpenalized = TRUE)
      }
      names(parse_i) <- sapply(interpret_i$smooth.spec, `[[`, 'label')
      parse_attr_i <- setNames(lapply(parse_i, FUN=function(i){i[[1]]$attr}),
        sapply(interpret_i$smooth.spec, FUN=function(i){paste0(i$term, collapse=',')}))
      parse_unpen_i <- lapply(parse_i, FUN=function(i){
        unpen_i <- lapply(i, FUN=function(j){
          if (!is.null(j$attr$unpen_x)){
            x <- j$attr$unpen_x
            colnames(x) <- paste('base @ ', 1:ncol(x))
            return(list(x=x))
          }else{return(NULL)}
        })
        if (all(sapply(unpen_i, FUN=function(j){is.null(j)}))){
          unpen_i <- NULL
        }
        return(unpen_i)
      })
      has_unpen <- !sapply(parse_unpen_i, is.null)
      parse_unpen_i <- parse_unpen_i[has_unpen]
      if (length(parse_unpen_i) > 0){
        names(parse_unpen_i) <- sapply(interpret_i$smooth.spec, `[[`, 'label')[has_unpen]
        names(parse_unpen_i) <- paste0('unpen_', names(parse_unpen_i))
        parse_i <- c(parse_i, parse_unpen_i)
      }
      any_zero <- which(sapply(parse_i, FUN=function(i){ncol(i[[1]]$x)}) == 0)
      if (length(any_zero) > 0){
        parse_i <- parse_i[-any_zero]
      }
      
      x <- model.frame(formula = interpret_i$pf, data = data)
      if (ncol(x) > 0){
        if (is.null(knots)){
          factor_x <- lapply(x, factor)
          levels_x <- lapply(factor_x, levels)
        }else{
          factor_x <- lapply(names(x), FUN=function(i){
            factor(x[[i]], levels = knots$levels[[i]])
          })
          names(factor_x) <- names(x)
          levels_x <- lapply(factor_x, levels)
        }
        M_matrix <- lapply(factor_x, FUN=function(i){
          t(fac2sparse(i, to = 'd', drop.unused.levels = FALSE))
        })      
        M_attr <- as.list(rep(NA, length(factor_x)))
        parse_param <- mapply(M_matrix, M_attr, SIMPLIFY = FALSE, FUN=function(i,j){
          list(list(x = i, attr = j))
        })
        parse_i <- c(parse_i, parse_param)
      }
      return(list(parse = parse_i, attr = parse_attr_i))
    }else{
      if (is.null(knots)){
        browser()
      }
      x <- as.list(data[i])
      if (is.null(knots)){
        factor_x <- lapply(x, factor)
        levels_x <- lapply(factor_x, levels)
      }else{
        factor_x <- lapply(names(x), FUN=function(i){
          factor(x[[i]], levels = knots$levels[[i]])
        })
        levels_x <- lapply(factor_x, levels)
      }
      M_matrix <- lapply(factor_x, FUN=function(i){t(fac2sparse(i, to = 'd'))})      
      M_id <- mapply(factor_x, levels_x, SIMPLIFY = FALSE,
        FUN=function(i,j){match(i,j)})
      browser()
      return(M_id)
    }
  })
  
  out_attr <- lapply(out, FUN=function(i){i$attr})
  out_attr <- do.call('c', out_attr)
  
  out <- lapply(out, FUN=function(i){i$parse})
  out_names <- lapply(out, names)
  out <- do.call('c', out)
  
  if (any(lengths(out) > 1)){
    browser()
  }else{
    out <- lapply(out, FUN=function(i){i[[1]]})
  }
  coef_storage <- lapply(out, FUN=function(i){
    matrix(data = NA, nrow = ncol(i$x), ncol = rank, dimnames = list(colnames(i$x), 1:rank))
  })
  M_matrix <- lapply(out, `[[`, "x")
  levels_x <- lapply(M_matrix, colnames)
  
  hier_fmla <- mapply(out_names, object$hier_fmla, SIMPLIFY = FALSE, FUN=function(i,j){
    new_i <- setdiff(i,j)
    return(c(base::intersect(i,j), new_i))
  })
  
  relevel <- do.call('c',hier_fmla)
  coef_storage <- coef_storage[relevel]
  levels_x <- levels_x[relevel]
  M_matrix <- M_matrix[relevel]

  special_attr <- list(storage = coef_storage, 
                       spline_attr = lapply(out, `[[`, "attr"),
                       #id = M_id, 
                       predict_attr = out_attr,
                       hier_grouping = hier_fmla, hier = hier,
                       levels = levels_x,
                       nested_hier = nested_hier)
  out <- list(x = M_matrix, attr = special_attr)
  class(out) <- c('mi_sparse')
  return(list(out))
}

print.spline_sparse <- function(x){
  print(x$x)
}
image.spline_sparse <- function(x){image(x$x)}

#' Interpret a vglmer formula for splines
#' @description A modified version of interpret.gam0 from mgcv. Used when mgcv's
#'   interpret.gam fails; usually when some environment object is passed to v_s.
#' @param gf A vglmer formula
#' @param textra Unused internal argument
#' @param extra.special Allow extra special terms to be passed
#' @importFrom stats reformulate terms.formula as.formula formula update.formula
#'   quantile
#' @keywords internal
fallback_interpret.gam0 <- function(gf, textra = NULL, extra.special = NULL){
  
  p.env <- environment(gf)
  
  tf <- terms.formula(gf, specials = c("s", "te", 
                                       "ti", "t2", extra.special))
  terms <- attr(tf, "term.labels")
  nt <- length(terms)
  if (attr(tf, "response") > 0) {
    response <- as.character(attr(tf, "variables")[2])
  }
  else {
    response <- NULL
  }
  sp <- attr(tf, "specials")$s
  tp <- attr(tf, "specials")$te
  tip <- attr(tf, "specials")$ti
  t2p <- attr(tf, "specials")$t2
  if (is.null(extra.special)){
    zp <- NULL  
  }else{
    # Get this to work with multiple extra.special types
    zp <- unlist(attr(tf, "specials")[extra.special])
  }
  off <- attr(tf, "offset")
  vtab <- attr(tf, "factors")
  if (length(sp) > 0) 
    for (i in 1:length(sp)) {
      ind <- (1:nt)[as.logical(vtab[sp[i], ])]
      sp[i] <- ind
    }
  if (length(tp) > 0) 
    for (i in 1:length(tp)) {
      ind <- (1:nt)[as.logical(vtab[tp[i], ])]
      tp[i] <- ind
    }
  if (length(tip) > 0) 
    for (i in 1:length(tip)) {
      ind <- (1:nt)[as.logical(vtab[tip[i], ])]
      tip[i] <- ind
    }
  if (length(t2p) > 0) 
    for (i in 1:length(t2p)) {
      ind <- (1:nt)[as.logical(vtab[t2p[i], ])]
      t2p[i] <- ind
    }
  if (length(zp) > 0) 
    for (i in 1:length(zp)) {
      ind <- (1:nt)[as.logical(vtab[zp[i], ])]
      zp[i] <- ind
    }
  k <- kt <- kti <- kt2 <- ks <- kz <- kp <- 1
  len.sp <- length(sp)
  len.tp <- length(tp)
  len.tip <- length(tip)
  len.t2p <- length(t2p)
  len.zp <- length(zp)
  ns <- len.sp + len.tp + len.tip + len.t2p + len.zp
  pav <- av <- rep("", 0)
  smooth.spec <- list()
  
  ###################
  # Modified from "mgcv"
  ####################
  
  mgcvns <- loadNamespace("vglmer")
  
  if (nt) 
    for (i in 1:nt) {
      if (k <= ns && ((ks <= len.sp && sp[ks] == i) || 
                      (kt <= len.tp && tp[kt] == i) || (kz <= len.zp && 
                                                        zp[kz] == i) || (kti <= len.tip && tip[kti] == 
                                                                         i) || (kt2 <= len.t2p && t2p[kt2] == i))) {
        
        ################
        # Modified from "mgcv::"
        #################
        st <- try(eval(parse(text = paste("vglmer::", 
                                          terms[i], sep = "")), envir = p.env), 
                  silent = TRUE)
        if (inherits(st, "try-error")) {
          st <- eval(parse(text = terms[i]), enclos = p.env, 
                     envir = mgcvns)
        }
        if (!is.null(textra)) {
          pos <- regexpr("(", st$lab, fixed = TRUE)[1]
          st$label <- paste(substr(st$label, start = 1, 
                                   stop = pos - 1), textra, substr(st$label, 
                                                                   start = pos, stop = nchar(st$label)), sep = "")
        }
        smooth.spec[[k]] <- st
        if (ks <= len.sp && sp[ks] == i) 
          ks <- ks + 1
        else if (kt <= len.tp && tp[kt] == i) 
          kt <- kt + 1
        else if (kti <= len.tip && tip[kti] == i) 
          kti <- kti + 1
        else if (kt2 <= len.t2p && t2p[kt2] == i) 
          kt2 <- kt2 + 1
        else kz <- kz + 1
        k <- k + 1
      }
      else {
        av[kp] <- terms[i]
        kp <- kp + 1
      }
    }
  if (!is.null(off)) {
    av[kp] <- as.character(attr(tf, "variables")[1 + 
                                                   off])
    kp <- kp + 1
  }
  pf <- paste(response, "~", paste(av, collapse = " + "))
  if (attr(tf, "intercept") == 0) {
    pf <- paste(pf, "-1", sep = "")
    if (kp > 1) 
      pfok <- 1
    else pfok <- 0
  }
  else {
    pfok <- 1
    if (kp == 1) {
      pf <- paste(pf, "1")
    }
  }
  fake.formula <- pf
  if (length(smooth.spec) > 0) 
    for (i in 1:length(smooth.spec)) {
      nt <- length(smooth.spec[[i]]$term)
      ff1 <- paste(smooth.spec[[i]]$term[1:nt], collapse = "+")
      fake.formula <- paste(fake.formula, "+", ff1)
      if (smooth.spec[[i]]$by != "NA") {
        fake.formula <- paste(fake.formula, "+", 
                              smooth.spec[[i]]$by)
        av <- c(av, smooth.spec[[i]]$term, smooth.spec[[i]]$by)
      }
      else av <- c(av, smooth.spec[[i]]$term)
    }
  fake.formula <- as.formula(fake.formula, p.env)
  if (length(av)) {
    pred.formula <- as.formula(paste("~", paste(av, 
                                                collapse = "+")))
    pav <- all.vars(pred.formula)
    pred.formula <- stats::reformulate(pav)
  }
  else pred.formula <- ~1
  ret <- list(pf = as.formula(pf, p.env), pfok = pfok, smooth.spec = smooth.spec, 
              fake.formula = fake.formula, response = response, fake.names = av, 
              pred.names = pav, pred.formula = pred.formula)
  class(ret) <- "split.gam.formula"
  ret
}
