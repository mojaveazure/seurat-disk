#' @include scdisk.R
#'
NULL

loom.validate <- list(
  '0.1.0' = function(lfile, verbose = TRUE) {
    .NotYetImplemented()
  },
  '3.0.0' = function(lfile, verbose = TRUE) {
    .NotYetImplemented()
  }
)

#' A class for connections to loom files
#'
#' @docType class
#' @name loom-class
#' @rdname loom-class
#' @aliases loom
#' @format An \code{\link[R6]{R6Class}} object
#' @seealso \code{\link[hdf5r]{H5File}}
#'
#' @importFrom R6 R6Class
#' @importFrom hdf5r H5File
#'
#' @export
#'
loom <- R6Class(
  classname = 'loom',
  inherit = scdisk,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    # Methods
  ),
  private = list(
    # Methods
    validate = function(validate = TRUE) {
      invisible(x = NULL)
    }
  )
)

#' @rdname as.loom
#' @method as.loom Seurat
#' @export
#'
as.loom.Seurat <- function(x, filename, overwrite = FALSE, verbose = TRUE, ...) {
  .NotYetImplemented()
}

#' Save an object as a loom file
#'
#' @inheritParams as.h5Seurat
#'
#' @return invisibly returns the loom file name
#'
#' @export
#'
SaveLoom <- function(
  object,
  filename = 'object.loom',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  lfile <- as.loom(
    x = object,
    filename = filename,
    overwrite = overwrite,
    verobse = verbose,
    ...
  )
  lfile$close_all()
  return(invisible(x = lfile$filename))
}
