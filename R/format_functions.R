#' Format Existing Objects
#'
#' Takes a regression output from glmer, stan, or vglmer and extracts the fixed
#' and random effects.
#'
#' @name format_obj
#'
#' @param object Object from glmer or lmer to the respective "format_"
#'   function
#' @importFrom stats vcov
#' @keywords internal
#' @export
format_glmer <- function(object) {
  
  output <- do.call('rbind', mapply(ranef(object), names(ranef(object)), SIMPLIFY = FALSE, FUN = function(i,j) {
    obj <- data.frame(
      var = as.vector(apply(attributes(i)$postVar, MARGIN = 3, FUN = function(i) {
        diag(i)
      })),
      mean = as.vector(t(as.matrix(i))),
      name = paste0(rep(colnames(i), nrow(i)), " @ ", rep(rownames(i), each = ncol(i))), stringsAsFactors = F
    )
    obj[[".re"]] <- j
    return(obj)
  }))
  output$name <- paste0(output[[".re"]], ' @ ', output[["name"]])
  output_fe <- data.frame(mean = fixef(object), var = diag(vcov(object)))
  output_fe$name <- rownames(output_fe)
  output_fe[[".re"]] <- NA
  output <- rbind(output, output_fe)
  output <- output[, (names(output) != ".re")]

  rownames(output) <- NULL

  return(output)
}
