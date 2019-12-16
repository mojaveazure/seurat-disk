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
  inherit = H5File,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    # Methods
    initialize = function(
      filename = NULL,
      mode = c('a', 'r', 'r+', 'w', 'w-'),
      ...
    ) {
      super$initialize(filename = filename, mode = mode,)
    },
    finalizer = function() {
      self$close_all(close_self = TRUE)
    }
  )
)
