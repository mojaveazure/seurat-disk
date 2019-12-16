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
  inherit = H5File,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    # Methods
    initialize = function(
      filename,
      mode = c('a', 'r', 'r+','w', 'w-'),
      ...
    ) {
      super$initialize(filename = filename, mode = mode, ...)
    },
    finalizer = function() {
      self$close_all(close_self = TRUE)
    }
  )
)
