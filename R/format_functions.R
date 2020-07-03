#' Format Existing Objects
#' 
#' Takes a regression output from glmer, stan, or vglmer and extracts the fixed
#' and random effects.
#'
#' @name format_obj
#' 
#' @param object Object from glmer, vglmer, or stan to the respective "format_"
#'   function
#' @importFrom dplyr bind_rows mutate
#' @importFrom stats vcov
#' @keywords internal
#' @export
format_glmer <- function(object){
  output <- bind_rows(lapply(ranef(object), FUN=function(i){
    obj <- data.frame(var = as.vector(apply(attributes(i)$postVar, MARGIN = 3, FUN=function(i){diag(i)})),
                      mean = as.vector(t(as.matrix(i))),
                      name = paste0(rep(colnames(i), nrow(i)), ' @ ', rep(rownames(i), each = ncol(i))), stringsAsFactors = F)
    return(obj)
  }), .id = '.re')
  output$name <- with(output, paste0(.re, ' @ ', name))
  output_fe <- data.frame(mean = fixef(object), var = diag(vcov(object)))
  output_fe$name <- rownames(output_fe)
  
  output <- bind_rows(output, output_fe)
  output <- output[, (names(output) != '.re')]
  
  rownames(output) <- NULL
  
  return(output)
}


#' Format Stan
#' @param useSigma Return variance component parameters from STAN? Default "FALSE".
#' @rdname format_obj
#' @importFrom stats var
#' @export
format_stan <- function(object, useSigma = FALSE){
  post_stan <- as.matrix(object)
  if (!useSigma){
    post_stan <- post_stan[,!grepl(colnames(post_stan), pattern='^Sigma')]
  }
  parse_stan_names <- strsplit(x = colnames(post_stan), split = '^b\\[| |\\]', perl = T)
  
  fmt_stan_names <- sapply(parse_stan_names, FUN=function(i){
    if (length(i) == 1){
      return(i)
    }else{
      i_one <- unlist(strsplit(i[3], split=':'))
      return(paste(i_one[1], i[2], i_one[2], sep=' @ '))
    }
  })
  colnames(post_stan) <- fmt_stan_names
  output <- data.frame(var = apply(post_stan, MARGIN = 2, var), mean = colMeans(post_stan))
  output$name <- rownames(output)
  rownames(output) <- NULL
  return(output)
}
