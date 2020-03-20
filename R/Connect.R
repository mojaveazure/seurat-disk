#' @include zzz.R
#' @include h5Seurat.R
#'
NULL

#' Connect to a single-cell HDF5 dataset
#'
#' @param filename Name of on-disk file
#' @param type Type of single-cell dataset to connect as; choose from:
#' \itemize{
#'  \item h5seurat
#' }
#' Leave as \code{NULL} to guess type from file extension
#' @param mode Mode to connect to data as; choose from:
#' \describe{
#'  \item{r}{Open existing dataset in read-only mode}
#'  \item{r+}{Open existing dataset in read/write mode}
#' }
#'
#' @return An object of class \code{type}, opened in mode \code{mode}
#'
#' @importFrom tools file_ext
#'
#' @export
#'
Connect <- function(filename, type = NULL, mode = c('r', 'r+')) {
  type <- tolower(x =  type %||% file_ext(x = filename))
  type <- match.arg(arg = type, choices = c('h5seurat'))
  mode <- match.arg(arg = mode)
  if (!file.exists(filename)) {
    stop("Cannot find ", type, " file ", filename, call. = FALSE)
  }
  return(switch(
    EXPR = type,
    'h5seurat' = h5Seurat$new(filename = filename, mode = mode)
  ))
}
