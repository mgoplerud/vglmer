# Functions for creating custom "smooth.construct" to implement
# in mgcv but also vglmer

gen.trend.D <- function(order, l){
  bcoef <- choose(n = order + 1, k = 0:(order+1))
  bcoef <- bcoef * (-1)^(0:(order +1))
  if (length(bcoef) > l){
    stop('Lower Polynomial Order: order must be >= 2  + l')
  }
  n.rows <- l - order - 1
  t.D <- do.call('rbind', lapply(1:n.rows, FUN=function(i){
    cbind(i, i - 1 + 1:length(bcoef), bcoef)
  }))
  t.D <- sparseMatrix(i = t.D[,1], j = t.D[,2], x = t.D[,3])
  if (ncol(t.D) != l){
    stop('Unusual Error')
  }
  return(t.D)
}

#' Constructor for random walk smooth
#' @import mgcv
#' @importFrom Matrix sparseMatrix
#' @keywords internal
#' @details See the function \link{v_s} for details.
#' @param object a smooth object; see documentation for other methods in
#'   \code{mgcv}.
#' @param data a data.frame; see documentation for other methods in
#'   \code{mgcv}.
#' @param knots not used
#' @export
smooth.construct.randwalk.smooth.spec <- function(object, data, knots){
  
  if (length(knots) != 0){stop('"knots" not used for random walk.')}
  if (is.null(object$xt)) {
    object$xt <- NULL
  }
  
  if (is.null(object$xt)){
    smooth_order <- 0
  }else{
    smooth_order <- object$xt$order
    if (is.null(smooth_order)){
      smooth_order <- 0
    }
  }
  
  x <- data[[object$term]]
  if (is.character(x) || is.factor(x)){
    if (!is.ordered(x)){
      stop('randwalk requires an *ordered* factor to be provided.')
    }
  }else{
    stop('randwalk only set up for character or factor; you may consider using "by" for evolution of a continuous covariate over time')
  }
  
  x <- data[[object$term]]
  lx <- levels(x)
  
  X <- sparseMatrix(i = 1:length(x), j = match(x, lx), x = 1, dims = c(length(x), length(lx)))
  
  # Get the random walk penalty matrix
  S <- crossprod(gen.trend.D(order = smooth_order, l = length(lx)))
  eS <- eigen(S)
  zero_ev <- which(eS$values < sqrt(.Machine$double.eps))
  if (length(zero_ev) != 1 || isTRUE(zero_ev != ncol(S))){
    stop('Invalid eigendecomposition...')
  }
  X <- X %*% eS$vectors
  last_col <- X[,length(lx)]
  stopifnot(var(last_col) < .Machine$double.eps)
  X <- X[,-length(lx)]
  X <- X %*% Diagonal(x=1/sqrt(eS$values[-length(lx)]))
  
  S <- Diagonal(n = ncol(X))
  X <- as.matrix(X)
  S <- as.matrix(S)
  
  object$X <- X
  object$S <- list(S)
  # Required elements for prediction
  object$transf_mat <- eS$vectors[,-ncol(eS$vectors)] %*% Diagonal(x=1/sqrt(eS$values[-length(lx)]))
  object$levels <- lx
  # Required elements for "gam"
  object$rank <- ncol(S)
  object$null.space.dim <- 0
  object$df <- ncol(S)
  object$te.ok <- 0
  object$plot.me <- FALSE
  object$C <- matrix(nrow = 0, ncol = ncol(S))
  # If rescale_penalty is NOT true, then set "no.rescale" to TRUE, i.e. to not rescale.
  # otherwise, leave null.
  if (!is.null(object$xt$rescale_penalty) && !object$xt$rescale_penalty) {
    object$no.rescale <- TRUE
  }
  class(object) <- 'randwalk.smooth'
  object
}

#' Predict Methods for random walk smooth
#' @keywords internal
#' @param object a smooth object; see documentation for other methods in
#'   \code{mgcv}.
#' @param data a data.frame; see documentation for other methods in
#'   \code{mgcv}.
#' @export
Predict.matrix.randwalk.smooth <- function(object, data){

  x <- data[[object$term]]
  lx <- object$levels
  if (!all(x %in% c(NA, lx))){
    stop('Prediction for random walk cannot work if levels in (new) data are not in original factor levels.')
    # if (object$internal_override){
    #   xj <- match(x, lx)
    #   id_valid <- which(!is.na(xj))
    #   X <- sparseMatrix(i = (1:length(x))[id_valid], j = xj[id_valid], x = 1, dims = c(length(x), length(lx)))
    # }else{
    #   stop('Prediction for random walk cannot work if levels not in original factor levels.')
    # }
  }else{
    xj <- match(x, lx)
    is_not_na <- which(!is.na(x))
    X <- sparseMatrix(i = (1:length(x))[is_not_na], j = xj[is_not_na], x = 1, dims = c(length(x), length(lx)))
  }
  X <- as.matrix(X %*% object$transf_mat)
  return(X)
}
