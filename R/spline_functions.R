
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
                by_re = TRUE,
                outer_okay = FALSE){
  if (!(type %in% c('tpf', 'o'))){stop('non tpf not set up yet...')}
  # Using mgcv's syntax for "s" to make it work with "interpret.gam"
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  if (d > 1){stop('Unlike mgcv, only provide a single variable')}
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
  
  ret <- list(term = term, outer_okay = outer_okay,
              by = by.var, type = type, knots = knots,
              by_re = by_re)
  class(ret) <- 'vglmer_spline'
  
  return(ret)
}

#' @importFrom splines spline.des
vglmer_build_spline <- function(x, knots = NULL, Boundary.knots = NULL, 
  by, type, override_warn = FALSE, outer_okay = FALSE, by_re = NULL){

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
  }else if (length(knots) == 1){
    numIntKnots <- knots
    
    intKnots <- quantile(unique(x),seq(0,1,length=
        (numIntKnots+2))[-c(1,(numIntKnots+2))])
    names(intKnots) <- NULL
  }else{
    # Sort user provided knots
    knots <- sort(knots)
    
    if (any(knots > max(x, na.rm=T)) | any(knots < min(x, na.rm=T))){
      if (!override_warn){
        warning('self-provided knots are outside of the observed data.')
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

  }else{stop('splines only set up for tpf and o')}
  
  spline_attr$by_re <- by_re
  
  if (!is.null(by)){
    
    base_x <- x
    u_by <- sort(unique(by))
    x_by <- sparseMatrix(i = 1:length(by), j = match(by, u_by), x = 1)
    
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
  zp <- if (is.null(extra.special)) 
    NULL
  else attr(tf, "specials")[[extra.special]]
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
