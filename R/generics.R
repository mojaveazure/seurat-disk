#' Save an object as an H5Seurat file
#'
#' @param x An object
#' @param filename Name of file to save \code{x} to
#' @param overwrite Overwrite \code{filename} if present?
#' @param verbose Show progress updates
#' @param ... Extra parameters passed to and from other methods
#'
#' @return ...
#'
#' @export
#'
as.h5Seurat <- function(
  x,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  UseMethod(generic = 'as.h5Seurat', object = x)
}

#' Save an object as a loom file
#'
#' @inheritParams as.h5Seurat
#'
#' @return ...
#'
#' @export
#'
as.loom <- function(
  x,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  UseMethod(generic = 'as.loom', object = x)
}
