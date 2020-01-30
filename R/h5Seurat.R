#' @include scdisk.R
#'
NULL

#' A class for connections to h5Seurat files
#'
#' @docType class
#' @name h5Seurat-class
#' @rdname h5Seurat-class
#' @aliases h5Seurat
#' @format An \code{\link[R6]{R6Class}} object
#' @seealso \code{\link[hdf5r]{H5File}}
#'
#' @importFrom R6 R6Class
#' @importFrom hdf5r H5File
#'
#' @export
#'
h5Seurat <- R6Class(
  classname = 'h5Seurat',
  inherit = scdisk,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    # Methods
  ),
  private = list(
    # Methods
    validate = function(...) {
      message("Validating h5Seurat file")
    }
  )
)

#' @rdname as.h5Seurat
#' @method as.h5Seurat Seurat
#' @export
#'
as.h5Seurat.Seurat <- function(
  x,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  .NotYetImplemented()
}

#' Save an object as an h5seurat file
#'
#' @inheritParams as.h5Seurat
#'
#' @return Invisibly returns h5seurat file name
#'
#' @export
#'
SaveH5Seurat <- function(
  object,
  filename = paste0(Project(object = object), '.h5seurat'),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  h5seurat <- as.h5Seurat(
    x = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose
  )
  h5seurat$close_all()
  return(invisible(x = NULL))
}
